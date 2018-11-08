#!/usr/bin/env python

"""APIs to query Genes in pathway databases"""

from __future__ import print_function
from sys import argv
import os
from StringIO import StringIO
from pprint import PrettyPrinter
from json import loads
import sqlite3
from pycurl import Curl
from pycurl import HTTP_CODE
import configparser
from requests import post
from wikipathways_api_client import WikipathwaysApiClient

EXCLUDE_FROM_CPDB = set(['KEGG', 'Reactome', 'Wikipathways'])
EXCLUDE_FROM_PC = set(['KEGG', 'Reactome', 'Wikipathways'])

def wikipathways(gene, exclude_reactome=True):
    """Wikipathways"""
    wp_client = WikipathwaysApiClient()
    kwargs = {'query': gene,
              'organism': 'http://identifiers.org/taxonomy/9606'
             }

    pathways = {('Wikipathways', pathway['identifier']): pathway['name']
                for pathway in wp_client.find_pathways_by_text(**kwargs)}

    if not exclude_reactome:
        return pathways

    base_url = 'https://www.wikipathways.org/index.php'
    client = Curl()
    query_buffer = StringIO()
    unique_pathways = {}
    for pathway_id in pathways:
        client.setopt(client.URL, '%s/Pathway:%s' % (base_url, pathway_id[1]))
        client.setopt(client.WRITEDATA, query_buffer)
        client.perform()

        response = query_buffer.getvalue()
        query_buffer.truncate(0)
        if 'http://www.reactome.org/PathwayBrowser/#DIAGRAM=' not in response:
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
                pathways[('Reactome', pathway['stId'].encode('utf-8'))] = pathway[
                    'displayName'].encode('utf-8')

    return pathways

def reactome_in_wikipathways(gene):
    """Identify pathways present in WikiPathways and Reactome"""
    base_url = 'https://www.wikipathways.org/index.php'
    query_buffer = StringIO()
    client = Curl()

    pathways = wikipathways(gene, False)
    common_pathways = {}
    for pathway_id in pathways:
        client.setopt(client.URL, '%s/Pathway:%s' % (base_url, pathway_id[1]))
        client.setopt(client.WRITEDATA, query_buffer)
        client.perform()
        response = query_buffer.getvalue()
        query_buffer.truncate(0)

        if 'http://www.reactome.org/PathwayBrowser/#DIAGRAM=' not in response:
            continue

        position = 48 + response.find('http://www.reactome.org/PathwayBrowser/#DIAGRAM=')
        reactome_id = 'R-HSA-'
        while response[position].isdigit():
            reactome_id += response[position]
            position += 1

        common_pathways[('Reactome', reactome_id)] = pathways[pathway_id]

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
    query_buffer.truncate(0)

    if not gene_entries:
        return {}

    gene_ids = []
    for entry in gene_entries:
        entry = entry.split('\t')
        if  entry[0].startswith('hsa:') and gene in entry[1].split(';')[0].split(','):
            gene_ids.append(entry[0])

    if not gene_ids:
        return {}

    pathways = {}
    for gene_id in gene_ids:
        client.setopt(client.URL, '%s/link/pathway/%s' % (api_url, gene_id))
        client.setopt(client.WRITEDATA, query_buffer)
        client.perform()
        pathway_ids = [entry.split('\t')[1] for entry
                       in query_buffer.getvalue().split('\n') if entry]
        query_buffer.truncate(0)

        if not pathway_ids:
            continue

        for pathway_id in pathway_ids:
            client.setopt(client.URL, '%s/get/%s' % (api_url, pathway_id))
            client.setopt(client.WRITEDATA, query_buffer)
            client.perform()
            pathway_response = query_buffer.getvalue().split('\n')
            pathway_name = pathway_response[1][5:-23].strip()
            pathways[('KEGG', pathway_id)] = pathway_name

    return pathways

def build_consensuspathdb():
    """Build inverted local ConsensusPathDB database"""
    url = 'http://cpdb.molgen.mpg.de/CPDB/getPathwayGenes'
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), '../docs/codes3d.conf'))
    db_dir = os.path.join(os.path.dirname(__file__), config.get('Defaults', 'CONSENSUS_DB_FP'))

    file_request = post(url, data={'idtype': 'hgnc-symbol'})
    if file_request.status_code != 200:
        return False

    db_connect = sqlite3.connect(db_dir)
    db_cursor = db_connect.cursor()
    db_cursor.execute("""DROP TABLE IF EXISTS pathways""")
    db_cursor.execute("""CREATE TABLE pathways
                         (database text, gene text, pathway_id text, pathway_label text)""")

    pathway_data = file_request.content.split('\n')
    del pathway_data[0]

    for line in pathway_data:
        line = line.split('\t')
        database = unicode(line[2], 'utf-8')
        if database in EXCLUDE_FROM_CPDB:
            continue
        pathway_id = unicode(line[1], 'utf-8')
        pathway_label = unicode(line[0], 'utf-8')
        genes = line[3].split(',')

        for gene in genes:
            gene = unicode(gene, 'utf-8')
            entry = (database, gene, pathway_id, pathway_label)
            db_cursor.execute("""INSERT INTO pathways
                                 VALUES (?, ?, ?, ?)""", entry)

    db_connect.commit()
    db_connect.close()

    return True

def remove_consensuspathdb():
    """Remove inverted local ConsensusPathDB database"""
    config = configparser.ConfigParser()
    config.read(os.path.join(os.path.dirname(__file__), '../docs/codes3d.conf'))
    db_dir = os.path.join(os.path.dirname(__file__), config.get('Defaults', 'CONSENSUS_DB_FP'))
    try:
        os.remove(db_dir)
    except OSError:
        return False
    return True

def consensuspathdb(gene):
    """ConsensusPathDB"""
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
        pathways[(entry[0], entry[1])] = entry[2]

    db_connect.close()

    return pathways


def get_pathways(gene):
    """All"""
    pathways = {}
    databases = (kegg, reactome, wikipathways, consensuspathdb, pathwaycommons)

    for database in databases:
        pathways.update(database(gene))

    return pathways

def pathwaycommons(gene):
    """Pathway Commons"""
    pass


def main():
    """Test functionality"""
    with open(argv[1]) as gene_file:
        genes = gene_file.read().splitlines()

    pprinter = PrettyPrinter()
    build_consensuspathdb()

    for gene in genes:
        print('\n%s' % gene)
        pprinter.pprint(pathwaycommons(gene))

    remove_consensuspathdb()

if __name__ == '__main__':
    main()
