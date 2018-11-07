#!/usr/bin/env python

"""APIs to query Genes in pathway databases"""

from __future__ import print_function
from sys import argv
from sys import exit
import os
from StringIO import StringIO
from pprint import PrettyPrinter
from json import loads
from csv import reader
import sqlite3
from pycurl import Curl
from pycurl import HTTP_CODE
import configparser
from wikipathways_api_client import WikipathwaysApiClient

def wikipathways(gene, exclude_reactome=True):
    """WikiPathways"""
    wp_client = WikipathwaysApiClient()
    kwargs = {'query': gene,
              'organism': 'http://identifiers.org/taxonomy/9606'
             }

    pathways = {pathway['identifier']: pathway['name']
                for pathway in wp_client.find_pathways_by_text(**kwargs)}

    if not exclude_reactome:
        return pathways

    wp_base_url = 'https://www.wikipathways.org/index.php'
    client = Curl()
    query_buffer = StringIO()
    unique_pathways = {}
    for pathway_id in pathways:
        client.setopt(client.URL, '%s/Pathway:%s' % (wp_base_url, pathway_id))
        client.setopt(client.WRITEDATA, query_buffer)
        client.perform()

        wp_response = query_buffer.getvalue()
        query_buffer.truncate(0)
        if 'http://www.reactome.org/PathwayBrowser/#DIAGRAM=' not in wp_response:
            unique_pathways[pathway_id] = pathways[pathway_id]

    return unique_pathways


def reactome(gene):
    """Reactome"""
    api_url = 'https://reactome.org/ContentService'
    query_buffer = StringIO()
    client = Curl()
    client.setopt(client.URL,
                  '%s/search/query?query=%s&cluster=true' % (api_url, gene))
    client.setopt(client.WRITEDATA, query_buffer)
    client.perform()

    if client.getinfo(HTTP_CODE) != 200:
        return {}

    gene_response = loads(query_buffer.getvalue())
    query_buffer.truncate(0)

    gene_ids = []
    for result in gene_response['results']:
        if result['typeName'] == 'Protein':
            for entry in result['entries']:
                if 'Homo sapiens' in entry['species']:
                    gene_ids.append(entry['stId'])

    pathways = {}
    for gene_id in gene_ids:
        for url_variant in ('', 'diagram/'):
            client.setopt(client.URL,
                          '%s/data/pathways/low/%sentity/%s/allForms?species=9606'
                          % (api_url, url_variant, gene_id))
            client.setopt(client.WRITEDATA, query_buffer)
            client.perform()

            if client.getinfo(HTTP_CODE) != 200:
                continue

            pathway_response = loads(query_buffer.getvalue())
            query_buffer.truncate(0)
            for pathway in pathway_response:
                pathways[pathway['stId'].encode('utf-8')] = pathway[
                    'displayName'].encode('utf-8')

    return pathways

def reactome_in_wikipathways(gene):
    """Identify pathways present in WikiPathways and Reactome"""
    wp_base_url = 'https://www.wikipathways.org/index.php'
    query_buffer = StringIO()
    client = Curl()

    pathways = wikipathways(gene, False)
    common_pathways = {}
    for pathway_id in pathways:
        client.setopt(client.URL, '%s/Pathway:%s' % (wp_base_url, pathway_id))
        client.setopt(client.WRITEDATA, query_buffer)
        client.perform()
        wp_response = query_buffer.getvalue()
        query_buffer.truncate(0)

        if 'http://www.reactome.org/PathwayBrowser/#DIAGRAM=' not in wp_response:
            continue

        position = wp_response.find('http://www.reactome.org/PathwayBrowser/#DIAGRAM=')
        position += 48
        reactome_id = 'R-HSA-'
        while wp_response[position].isdigit():
            reactome_id += wp_response[position]
            position += 1

        try:
            common_pathways[reactome_id][pathway_id] = pathways[pathway_id]
        except KeyError:
            common_pathways[reactome_id] = {}
            common_pathways[reactome_id][pathway_id] = pathways[pathway_id]

    return common_pathways


def kegg(gene):
    """KEGG"""
    api_url = 'http://rest.kegg.jp'
    query_buffer = StringIO()
    client = Curl()

    client.setopt(client.URL, '%s/find/genes/%s' % (api_url, gene))
    client.setopt(client.WRITEDATA, query_buffer)
    client.perform()
    gene_entries = [entry for entry in
                    query_buffer.getvalue().split('\n') if entry]

    if not gene_entries:
        return {}

    gene_ids = []
    for entry in gene_entries:
        entry = entry.split('\t')
        if  entry[0].startswith('hsa:') and gene in entry[1].split(
                ';')[0].split(','):
            gene_ids.append(entry[0])

    if not gene_ids:
        return {}

    pathways = {}
    for gene_id in gene_ids:
        query_buffer.truncate(0)
        client.setopt(client.URL, '%s/link/pathway/%s' % (api_url, gene_id))
        client.setopt(client.WRITEDATA, query_buffer)
        client.perform()
        pathway_ids = [entry.split('\t')[1] for entry
                       in query_buffer.getvalue().split('\n') if entry]

        if not pathway_ids:
            continue

        for pathway_id in pathway_ids:
            query_buffer.truncate(0)
            client.setopt(client.URL, '%s/get/%s' % (api_url, pathway_id))
            client.setopt(client.WRITEDATA, query_buffer)
            client.perform()
            pathway_response = query_buffer.getvalue().split('\n')
            pathway_name = pathway_response[1][5:-23].strip()
            pathways[pathway_id] = pathway_name

    query_buffer.truncate(0)
    return pathways

def build_consensus_db():
    """Build, respectively update local ConsensusDB database"""
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__),  '../docs/codes3d.conf'))
    tab_dir = os.path.join(os.path.dirname(__file__), config.get('Defaults', 'CONSENSUS_DB_TAB_FP'))
    db_dir = os.path.join(os.path.dirname(__file__), config.get('Defaults', 'CONSENSUS_DB_FP'))

    try:
        pathway_file = open(tab_dir, 'r')
    except IOError:
        print('%s could not be updated.' % db_dir)
        exit()

    db_connect = sqlite3.connect(db_dir)
    db_cursor = db_connect.cursor()
    db_cursor.execute("""DROP TABLE IF EXISTS pathways""")
    db_cursor.execute("""CREATE TABLE pathways (database text, gene text, pathway_id text, pathway_label text)""")

    tsv_reader = reader(pathway_file, delimiter='\t')
    next(tsv_reader, None)

    for line in tsv_reader:
        pathway_label = unicode(line[0], 'utf-8')
        pathway_id = unicode(line[1], 'utf-8')
        database = unicode(line[2], 'utf-8')
        genes = line[3].split(',')

        for gene in genes:
            gene = unicode(gene, 'utf-8')
            entry = (database, gene, pathway_id, pathway_label)
            db_cursor.execute("""INSERT INTO pathways
                                 VALUES (?, ?, ?, ?)""", entry)

    db_connect.commit()
    db_connect.close()
    print('%s successfully built.' % db_dir)

def consensus_db(gene):
    """ConsensusDB"""
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__),
        '../docs/codes3d.conf'))
    db_dir = os.path.join(os.path.dirname(__file__), config.get('Defaults',
        'CONSENSUS_DB_FP'))

    db_connect = sqlite3.connect(db_dir)
    db_cursor = db_connect.cursor()

    pathways = {}
    query = (gene,)
    for entry in db_cursor.execute("""SELECT database, pathway_id, pathway_label
                                      FROM pathways
                                      WHERE gene=?""", query):
        entry = [subentry.encode('utf-8') for subentry in entry]
        try:
            pathways[entry[0]][entry[1]] = entry[2]
        except KeyError:
            pathways[entry[0]] = {}
            pathways[entry[0]][entry[1]] = entry[2]

    db_connect.close()
    return pathways

def test(genes):
    """Test functionality"""
    pprinter = PrettyPrinter()

    build_consensus_db()
    for gene in genes:
        print("\nConsensusPathDB: %s" % gene)
        pprinter.pprint(consensus_db(gene))
        print("\nKEGG: %s" % gene)
        pprinter.pprint(kegg(gene))
        print("\nWikiPathways: %s" % gene)
        pprinter.pprint(wikipathways(gene, False))
        print("\nReactome: %s" % gene)
        pprinter.pprint(reactome(gene))
        print("\nReactome in WikiPathways: %s" % gene)
        pprinter.pprint(reactome_in_wikipathways(gene))

def main():
    """Main"""
    with open(argv[1]) as gene_file:
        gene_identifiers = gene_file.read().splitlines()
    test(gene_identifiers)

if __name__ == '__main__':
    main()
