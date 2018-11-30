#!/usr/bin/env python

"""APIs to query HGNC-Gene Symbols in KEGG, Reactome and WikiPathways"""

from __future__ import print_function
import os
from time import sleep
from StringIO import StringIO
from json import loads as json
import logging
import sqlite3
import argparse
import unicodecsv as csv
from pycurl import Curl
from pycurl import error as PycurlError
from pycurl import HTTP_CODE
import configparser
import progressbar
import xml.etree.ElementTree as ET

KEGG_API_URL = "http://rest.kegg.jp"
KEGG_INFO_URL = "http://rest.kegg.jp/info/pathway"
REACTOME_API_URL = "https://www.reactome.org/ContentService"
REACTOME_INFO_URL = "https://reactome.org/ContentService/data/database/version"
REACTOME_ORIGIN_SUFFIX = "reactome.org/PathwayBrowser/#DIAGRAM="
WIKIPATHWAYS_API_URL = "https://webservice.wikipathways.org"
WIKIPATHWAYS_ENTRY_URL = "https://www.wikipathways.org/index.php/Pathway:"
WIKIPATHWAYS_NAMESPACES = {"ns1": "http://www.wso2.org/php/xsd",
                           "ns2": "http://www.wikipathways.org/webservice"}

def kegg(gene):
    """Query Kyoto Encyclopedia of Genes and Genomes for HGNC-Symbol"""
    query_buffer = StringIO()
    client = Curl()

    client.setopt(client.URL, "%s/find/genes/%s" % (KEGG_API_URL, gene))
    client.setopt(client.WRITEDATA, query_buffer)

    try:
        client.perform()
    except PycurlError:
        return []

    if client.getinfo(HTTP_CODE) != 200:
        return []

    gene_response = query_buffer.getvalue().split("\n")
    query_buffer.truncate(0)
    gene_entries = [entry for entry in gene_response if entry]

    if not gene_entries:
        return []

    gene_ids = set()
    for entry in gene_entries:
        entry = entry.split("\t")
        if  (entry[0].startswith("hsa:") and
             gene in entry[1].split(";")[0].split(",")):
            gene_ids.add(entry[0])

    if not gene_ids:
        return []

    pathways = set()
    for gene_id in gene_ids:
        client.setopt(client.URL, "%s/link/pathway/%s"
                      % (KEGG_API_URL, gene_id))
        client.setopt(client.WRITEDATA, query_buffer)

        try:
            client.perform()
        except PycurlError:
            continue

        if client.getinfo(HTTP_CODE) != 200:
            continue

        pathway_response = query_buffer.getvalue().split("\n")
        query_buffer.truncate(0)
        pathway_ids = [entry.split("\t")[1] for entry
                       in pathway_response if entry]

        if not pathway_ids:
            continue

        for pathway_id in pathway_ids:
            client.setopt(client.URL, "%s/get/%s" % (KEGG_API_URL, pathway_id))
            client.setopt(client.WRITEDATA, query_buffer)

            try:
                client.perform()
            except PycurlError:
                continue

            if client.getinfo(HTTP_CODE) != 200:
                continue

            pathway_response = query_buffer.getvalue().split("\n")
            query_buffer.truncate(0)

            for entry in pathway_response:
                if entry.startswith("NAME"):
                    pathway_name = entry.replace("NAME", "")
                    break
                pathway_name = "NA"

            pathway_id = pathway_id.replace("path:", "")
            pathway_name = pathway_name.split("-")[0]

            pathways.add((unicode(gene, "utf-8"),
                          unicode("KEGG", "utf-8"),
                          unicode(pathway_id, "utf-8"),
                          unicode("NA"),
                          unicode(pathway_name, "utf-8")))

    return sorted(list(pathways), key=lambda pathway: pathway[2])


def reactome(gene):
    """Query Reactome for HGNC-Symbol"""
    query_buffer = StringIO()
    client = Curl()
    client.setopt(client.URL,
                  "%s/search/query?query=%s&cluster=true"
                  % (REACTOME_API_URL, gene))
    client.setopt(client.WRITEDATA, query_buffer)

    try:
        client.perform()
    except PycurlError:
        return []

    if client.getinfo(HTTP_CODE) != 200:
        return []

    try:
        gene_response = json(query_buffer.getvalue())
    except ValueError:
        query_buffer.truncate(0)
        return []
    query_buffer.truncate(0)

    gene_ids = set()
    for result in gene_response["results"]:
        if result["typeName"] == "Protein":
            for entry in result["entries"]:
                if "Homo sapiens" in entry["species"]:
                    gene_ids.add(entry["stId"])

    pathways = set()
    for gene_id in gene_ids:
        for url_variant in ("", "diagram/"):
            client.setopt(
                client.URL,
                "%s/data/pathways/low/%sentity/%s/allForms?species=9606"
                % (REACTOME_API_URL, url_variant, gene_id))
            client.setopt(client.WRITEDATA, query_buffer)

            try:
                client.perform()
            except PycurlError:
                continue

            if client.getinfo(HTTP_CODE) != 200:
                continue

            try:
                pathway_response = json(query_buffer.getvalue())
            except ValueError:
                return []

            query_buffer.truncate(0)

            for pathway in pathway_response:
                pathways.add((
                    unicode(gene, "utf-8"),
                    unicode("Reactome", "utf-8"),
                    pathway["stId"],
                    pathway["stIdVersion"].split(".")[1],
                    pathway["displayName"]))

    return sorted(list(pathways), key=lambda pathway: pathway[2])

def wikipathways(gene, exclude_reactome=True):
    """Query WikiPathways for HGNC-Symbol"""
    query_buffer = StringIO()
    client = Curl()
    client.setopt(client.URL,
                  "%s/findPathwaysByText?query=%s&species=Homo_sapiens"
                  % (WIKIPATHWAYS_API_URL, gene))
    client.setopt(client.WRITEDATA, query_buffer)

    try:
        client.perform()
    except PycurlError:
        return []

    if client.getinfo(HTTP_CODE) != 200:
        return []

    try:
        response = ET.fromstring(query_buffer.getvalue())
    except ET.ParseError:
        return []

    query_buffer.truncate(0)

    pathways = []
    for result in response.findall("ns1:result", WIKIPATHWAYS_NAMESPACES):
        pathway = {}
        pathway["id"] = result.findtext("ns2:id",
            namespaces=WIKIPATHWAYS_NAMESPACES)
        pathway["revision"] = result.findtext("ns2:revision",
            namespaces=WIKIPATHWAYS_NAMESPACES)
        pathway["name"] = result.findtext("ns2:name",
            namespaces=WIKIPATHWAYS_NAMESPACES)
        pathways.append(pathway)

    client.setopt(client.URL,
                  "%s/findInteractions?query=%s"
                  % (WIKIPATHWAYS_API_URL, gene))
    client.setopt(client.WRITEDATA, query_buffer)

    try:
        client.perform()
    except PycurlError:
        return []

    if client.getinfo(HTTP_CODE) != 200:
        return []

    try:
        response = ET.fromstring(query_buffer.getvalue())
    except ET.ParseError:
        return []

    query_buffer.truncate(0)

    interactions = {pathway["id"]: {"left": set(),
                                    "right": set()}
                    for pathway in pathways}

    for result in response.findall("ns1:result", WIKIPATHWAYS_NAMESPACES):
        species = result.findtext("ns2:species",
                                  namespaces=WIKIPATHWAYS_NAMESPACES)
        if species == "Homo sapiens":
            pathway_id = result.findtext("ns2:id",
                                         namespaces=WIKIPATHWAYS_NAMESPACES)
            for field in result.findall("ns2:fields", WIKIPATHWAYS_NAMESPACES):
                name = field.findtext("ns2:name",
                                      namespaces=WIKIPATHWAYS_NAMESPACES)
                if name in ("source", "indexerId"):
                    continue
                values = set()
                for value in field.findall("ns2:values",
                                           namespaces=WIKIPATHWAYS_NAMESPACES):
                    try:
                       values.add(value.text.upper())
                    except AttributeError:
                        continue
                if gene.upper() in values:
                    continue
                interactions[pathway_id][name] |= values

    for pathway_id in interactions:
        interactions[pathway_id]["left"] =\
                sorted(list(interactions[pathway_id]["left"]))
        interactions[pathway_id]["right"] =\
                sorted(list(interactions[pathway_id]["right"]))

    pathways = set((unicode(gene, "utf-8"),
                    unicode("WikiPathways", "utf-8"),
                    unicode(pathway["id"], "utf-8"),
                    unicode(pathway["revision"], "utf-8"),
                    unicode(pathway["name"], "utf-8"),
                    unicode(" ".join(interactions[pathway["id"]]["left"])),
                    unicode(" ".join(interactions[pathway["id"]]["right"])))
                   for pathway in pathways)

    if not exclude_reactome:
        return sorted(list(pathways), key=lambda pathway: pathway[2])

    client = Curl()
    query_buffer = StringIO()

    unique_pathways = set()
    for pathway in pathways:
        client.setopt(client.URL, WIKIPATHWAYS_ENTRY_URL + pathway[2])
        client.setopt(client.WRITEDATA, query_buffer)

        try:
            client.perform()
        except PycurlError:
            continue

        if client.getinfo(HTTP_CODE) != 200:
            continue

        response = query_buffer.getvalue()
        query_buffer.truncate(0)
        if REACTOME_ORIGIN_SUFFIX not in response:
            unique_pathways.add(pathway)

    return sorted(list(unique_pathways), key=lambda pathway: pathway[2])

def log_kegg_release():
    """Write queried KEGG release information to log file"""
    query_buffer = StringIO()
    client = Curl()

    client.setopt(client.URL, KEGG_INFO_URL)
    client.setopt(client.WRITEDATA, query_buffer)

    try:
        client.perform()
    except PycurlError:
        logging.warning("\tKEGG release not retreived.")
        return

    if client.getinfo(HTTP_CODE) != 200:
        logging.warning("\tKEGG release not retreived.")
        return

    release_response = query_buffer.getvalue().split("\n")
    query_buffer.truncate(0)
    release = release_response[1].replace("path", "").strip()
    logging.info("KEGG version: %s", release)

def log_reactome_release():
    """Write queried Reactome release information to log file"""
    query_buffer = StringIO()
    client = Curl()

    client.setopt(client.URL, REACTOME_INFO_URL)
    client.setopt(client.WRITEDATA, query_buffer)

    try:
        client.perform()
    except PycurlError:
        logging.warning("\tReactome release not retreived.")
        return

    if client.getinfo(HTTP_CODE) != 200:
        logging.warning("\tReactome release not retreived.")
        return

    release = query_buffer.getvalue()
    query_buffer.truncate(0)
    logging.info("Reactome version: %s", release)

def log_wikipathways_release():
    """Write WikiPathways release information to log file"""
    logging.info("WikiPathways version: NA")

def build_pathway_db():
    """Build the pathway database"""
    logging.info("Starting build: %s", PATHWAY_DB_TMP_FP)
    logging.info("HGNC source: %s", GENE_ID_FP)
    log_kegg_release()
    log_reactome_release()
    log_wikipathways_release()

    with open(GENE_ID_FP) as gene_file:
        genes = gene_file.read().splitlines()

    gene_num_bar = progressbar.ProgressBar(max_value=len(genes))

    genes = sorted(dict.fromkeys(genes).keys(),
                   key=lambda gene: gene.capitalize())

    tab_file = open(PATHWAY_DB_TAB_FP, "w")
    tab_writer = csv.writer(tab_file, delimiter="\t")

    pathway_db = sqlite3.connect(PATHWAY_DB_TMP_FP)
    pathway_db_cursor = pathway_db.cursor()
    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS genes
        """)
    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS pathways
        """)
    pathway_db_cursor.execute("""
        CREATE TABLE genes (
            database TEXT,
            pathway TEXT,
            gene TEXT,
            PRIMARY KEY (
                database,
                pathway,
                gene)
            ON CONFLICT IGNORE)
        """)
    pathway_db_cursor.execute("""
        CREATE TABLE pathways (
            database TEXT,
            pathway TEXT,
            pathway_name TEXT,
            PRIMARY KEY (
                database,
                pathway,
                pathway_name)
            ON CONFLICT IGNORE)
        """)
    for i, gene in enumerate(genes):
        for database in (kegg, reactome, wikipathways):
            for pathway in database(gene):
                pathway_db_cursor.execute("""
                    INSERT INTO genes
                    VALUES (?, ?, ?)
                    """, (pathway[1], pathway[2], pathway[0]))

                pathway_db_cursor.execute("""
                    INSERT INTO pathways
                    VALUES (?, ?, ?)
                    """, (pathway[1], pathway[2], pathway[4]))

                tab_writer.writerow(pathway)
        gene_num_bar.update(i)

    pathway_db.commit()
    pathway_db.close()
    tab_file.close()

    os.rename(PATHWAY_DB_TMP_FP, PATHWAY_DB_FP)
    logging.info("Build complete: %s", PATHWAY_DB_FP)


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description="")
    PARSER.add_argument("-c", "--config",
                        default=os.path.join(os.path.dirname(__file__),
                                             "../docs/codes3d.conf"),
                        help="The configuration file to be used.")

    ARGS = PARSER.parse_args()

    CONFIG = configparser.ConfigParser()
    CONFIG.read(ARGS.config)

    GENE_ID_FP = os.path.join(os.path.dirname(__file__),
                              CONFIG.get("Defaults", "GENE_ID_FP"))
    PATHWAY_DB_FP = os.path.join(os.path.dirname(__file__),
                                 CONFIG.get("Defaults", "PATHWAY_DB_FP"))
    PATHWAY_DB_TMP_FP = os.path.join(os.path.dirname(__file__),
                                 CONFIG.get("Defaults", "PATHWAY_DB_TMP_FP"))
    PATHWAY_DB_TAB_FP = os.path.join(os.path.dirname(__file__),
                                     CONFIG.get("Defaults",
                                                "PATHWAY_DB_TAB_FP"))
    PATHWAY_DB_LOG_FP = os.path.join(os.path.dirname(__file__),
                                     CONFIG.get("Defaults",
                                                "PATHWAY_DB_LOG_FP"))

    logging.getLogger("requests").setLevel(logging.WARNING)
    logging.basicConfig(
        filename=PATHWAY_DB_LOG_FP,
        level=logging.INFO,
        format="%(asctime)s\t%(levelname)s\t%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")
    build_pathway_db()
