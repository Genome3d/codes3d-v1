#!/usr/bin/env python

"""APIs to query Genes in pathway databases"""

from __future__ import print_function
from sys import argv
import os
from StringIO import StringIO
from json import loads
import sqlite3
from pycurl import Curl
from pycurl import HTTP_CODE
import argparse
import configparser
import csv
import progressbar
from wikipathways_api_client import WikipathwaysApiClient

def wikipathways(gene, exclude_reactome=True):
    """Wikipathways"""
    wp_client = WikipathwaysApiClient()
    kwargs = {'query': gene,
              'organism': 'http://identifiers.org/taxonomy/9606'
             }

    pathways = set((gene, 'Wikipathways', pathway['identifier'], pathway['name'])
                for pathway in wp_client.find_pathways_by_text(**kwargs))

    if not exclude_reactome:
        return sorted(list(pathways), key=lambda pathway: pathway[2])

    base_url = 'https://www.wikipathways.org/index.php'
    client = Curl()
    query_buffer = StringIO()
    unique_pathways = set()
    for pathway in pathways:
        client.setopt(client.URL, '%s/Pathway:%s' % (base_url, pathway[2]))
        client.setopt(client.WRITEDATA, query_buffer)
        client.perform()

        if client.getinfo(HTTP_CODE) != 200:
            continue

        response = query_buffer.getvalue()
        query_buffer.truncate(0)
        if 'http://www.reactome.org/PathwayBrowser/#DIAGRAM=' not in response:
            unique_pathways.add(pathway)


    return sorted(list(unique_pathways), key=lambda pathway: pathway[2])

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
        return []

    gene_response = loads(query_buffer.getvalue())
    query_buffer.truncate(0)

    gene_ids = []
    for result in gene_response['results']:
        if result['typeName'] == 'Protein':
            for entry in result['entries']:
                if 'Homo sapiens' in entry['species']:
                    gene_ids.append(entry['stId'])

    pathways = set()
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
                pathways.add((gene, 'Reactome', pathway['stId'].encode('utf-8'),
                                pathway['displayName'].encode('utf-8')))

    return sorted(list(pathways), key=lambda pathway: pathway[2])

def kegg(gene):
    """KEGG"""
    api_url = "http://rest.kegg.jp"
    query_buffer = StringIO()
    client = Curl()

    client.setopt(client.URL, "%s/find/genes/%s" % (api_url, gene))
    client.setopt(client.WRITEDATA, query_buffer)
    client.perform()

    if client.getinfo(HTTP_CODE) != 200:
       return []

    gene_entries = [entry for entry in
                    query_buffer.getvalue().split('\n') if entry]
    query_buffer.truncate(0)

    if not gene_entries:
        return []

    gene_ids = []
    for entry in gene_entries:
        entry = entry.split('\t')
        if  entry[0].startswith('hsa:') and gene in entry[1].split(';')[0].split(','):
            gene_ids.append(entry[0])

    if not gene_ids:
        return []

    pathways = set()
    for gene_id in gene_ids:
        client.setopt(client.URL, '%s/link/pathway/%s' % (api_url, gene_id))
        client.setopt(client.WRITEDATA, query_buffer)
        client.perform()

        if client.getinfo(HTTP_CODE) != 200:
            return []

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
            pathways.add((gene, 'KEGG', pathway_id[5:], pathway_name))

    return sorted(list(pathways), key=lambda pathway: pathway[2])


def build_pathway_db():
    with open(GENE_ID_FP) as gene_file:
        genes = gene_file.read().splitlines()
    bar = progressbar.ProgressBar(max_value=len(genes))
    bar.update(0)
    genes = sorted(dict.fromkeys(genes).keys(),
                   key=lambda gene: gene.capitalize())
    tab_file = open(PATHWAY_DB_TAB_FP, 'w')
    tab_writer = csv.writer(tab_file, delimiter='\t')
    pathway_db = sqlite3.connect(PATHWAY_DB_FP)
    pathway_db_cursor = pathway_db.cursor()
    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS pathways
        """)
    pathway_db_cursor.execute("""
        CREATE TABLE pathways (gene text, database text, pathway_id text, pathway text,
        PRIMARY KEY (gene, database, pathway_id))
        """)

    header = ['HGNC Symbol', 'Database', 'Pathway ID', 'Pathway Designation', 'Reference URL']
    tab_writer.writerow(header)
    databases = ((kegg, KEGG_BASE_URL),
                 (reactome, REACTOME_BASE_URL),
                 (wikipathways, WIKIPATHWAYS_BASE_URL))
    for i, gene in enumerate(genes):
        for database, url in databases:
            for pathway in database(gene):
                pathway_db_cursor.execute("""
                    INSERT INTO pathways
                    VALUES (?, ?, ?, ?)
                    """, pathway)
                doc_line = list(pathway)
                doc_line.append('%s%s' % (url, doc_line[2]))
                tab_writer.writerow(doc_line)
        bar.update(i)

    pathway_db.commit()
    pathway_db.close()
    tab_file.close()
    return None

if __name__ == '__main__':
    KEGG_BASE_URL = 'https://www.genome.jp/dbget-bin/www_bget?'
    REACTOME_BASE_URL = 'https://reactome.org/PathwayBrowser/#DIAGRAM='
    WIKIPATHWAYS_BASE_URL = 'https://www.wikipathways.org/index.php/Pathway:'
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-c", "--config", default=os.path.join(os.path.dirname(__file__), "../docs/codes3d.conf"), help="The configuration file to be used.")
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)

    GENE_ID_FP = os.path.join(os.path.dirname(__file__), config.get("Defaults","GENE_ID_FP"))
    PATHWAY_DB_FP = os.path.join(os.path.dirname(__file__), config.get("Defaults","PATHWAY_DB_FP"))
    PATHWAY_DB_TAB_FP = os.path.join(os.path.dirname(__file__), config.get("Defaults","PATHWAY_DB_TAB_FP"))
    build_pathway_db()
