#!/usr/bin/env python

"""Clients to query HGNC-Gene-IDs in KEGG, Reactome and WikiPathways"""

from __future__ import print_function
import os
import logging
import sqlite3
import argparse
import xml.etree.ElementTree
import time
import unicodecsv
import configparser
import progressbar
import requests

def query(url):
    """Send request to specified URL"""
    while True:
        try:
            response = requests.get(url)
            break
        except requests.exceptions.HTTPError:
            logging.critical("Unsucessful status code %s: %s",
                             response.status_code, url)
            return None
        except requests.exceptions.ConnectionError:
            logging.critical("Connection failure: %s", url)
            time.sleep(60)
            continue
        except requests.exceptions.RequestException:
            logging.critical("General Error: %s", url)
            return None

    if response.headers["content-type"].split(";")[0] == "application/json":
        try:
            content = response.json()
        except ValueError:
            logging.critical("Invalid JSON content: %s", url)
            content = None
        return content

    elif response.headers["content-type"].split(";")[0] == "application/xml":
        try:
            content = xml.etree.ElementTree.fromstring(response.text)
        except xml.etree.ElementTree.ParseError:
            logging.critical("Invalid XML content: %s", url)
            content = None
        return content

    return response.text

def kegg(gene):
    """Query Kyoto Encyclopedia of Genes and Genomes for HGNC-Gene-ID"""
    kegg_api_url = "http://rest.kegg.jp"

    gene_response = query("{}/find/genes/{}".format(kegg_api_url, gene))

    if not gene_response:
        return []

    gene_entries = gene_response.split("\n")[:-1]

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
        pathway_response = query("{}/link/pathway/{}".format(
            kegg_api_url, gene_id))

        if not pathway_response:
            continue

        pathway_entries = pathway_response.split("\n")[:-1]
        pathway_ids = [entry.split("\t")[1] for entry in pathway_entries]

        for pathway_id in pathway_ids:
            pathway_response = query("{}/get/{}".format(
                kegg_api_url, pathway_id))

            if not pathway_response:
                continue

            pathway_entries = pathway_response.split("\n")
            pathway_id = pathway_id.replace("path:", "")

            for entry in pathway_entries:
                if entry.startswith("NAME"):
                    pathway_name =\
                            entry.replace("NAME", "").split("-")[0].strip()
                    break
                pathway_name = "NA"

            pathways.add((unicode(gene, "utf-8"),
                          unicode("KEGG", "utf-8"),
                          pathway_id,
                          unicode("NA"),
                          pathway_name))

    return sorted(list(pathways), key=lambda pathway: pathway[2])


def reactome(gene):
    """Query Reactome for HGNC-Gene-ID"""
    reactome_api_url = "https://www.reactome.org/ContentService"

    response = query(
        "{}/search/query?query={}&species=Homo%%20sapiens&cluster=true".format(
            reactome_api_url, gene))

    if not response:
        return []

    pathway_ids = set()

    for result in response["results"]:
        if result["typeName"] == "Pathway":
            for entry in result["entries"]:
                pathway_ids.add(entry["stId"])

    pathways = dict.fromkeys(pathway_ids)

    for pathway_id in pathway_ids:
        response = query("{}/data/query/{}".format(
            reactome_api_url, pathway_id))

        if not response:
            del pathways[pathway_id]
            continue

        pathways[pathway_id] = {}
        pathways[pathway_id]["name"] = response["displayName"]
        pathways[pathway_id]["version"] = response["stIdVersion"].split(".")[1]

    reactions_to_pathways = dict.fromkeys(pathways, set())
    reaction_ids = set()
    for pathway_id in reactions_to_pathways:
        response = query("{}/data/pathway/{}/containedEvents/stId".format(
            reactome_api_url, pathway_id))

        if not response:
            continue

        response_values = set(response[1:-1].split(", "))

        reactions_to_pathways[pathway_id] = response_values
        reaction_ids |= response_values

    reactions = dict.fromkeys(reaction_ids)
    for reaction_id in reactions:
        reactions[reaction_id] = {}
        reactions[reaction_id]["input"] = set()
        reactions[reaction_id]["output"] = set()

        response = query("{}/data/query/{}".format(
            reactome_api_url, reaction_id))

        if not response:
            continue

        for metabolite in response["input"]:
            if (isinstance(metabolite, dict) and
                    metabolite["className"] == "Protein"):
                reactions[reaction_id]["input"].add(
                    metabolite["displayName"].split(" [")[0].split("(")[0])

        for metabolite in response["output"]:
            if (isinstance(metabolite, dict) and
                    metabolite["className"] == "Protein"):
                reactions[reaction_id]["output"].add(
                    metabolite["displayName"].split(" [")[0].split("(")[0])

    for pathway_id in pathways:
        pathways[pathway_id]["input"] = set()
        pathways[pathway_id]["output"] = set()
        for reaction_id in reactions_to_pathways[pathway_id]:
            if (gene in reactions[reaction_id]["input"] and
                    gene not in reactions[reaction_id]["output"]):
                pathways[pathway_id]["output"] |=\
                        reactions[reaction_id]["output"]
            if (gene in reactions[reaction_id]["output"] and
                    gene not in reactions[reaction_id]["input"]):
                pathways[pathway_id]["input"] |=\
                        reactions[reaction_id]["input"]
        pathways[pathway_id]["input"] =\
                sorted(list(pathways[pathway_id]["input"]))
        pathways[pathway_id]["output"] =\
                sorted(list(pathways[pathway_id]["output"]))

    pathways = set((unicode(gene, "utf-8"),
                    unicode("Reactome", "utf-8"),
                    pathway_id,
                    pathways[pathway_id]["version"],
                    pathways[pathway_id]["name"],
                    tuple(pathways[pathway_id]["input"]),
                    tuple(pathways[pathway_id]["output"]))
                   for pathway_id in pathways
                   if (pathways[pathway_id]["input"] or
                       pathways[pathway_id]["output"]))

    return sorted(list(pathways), key=lambda pathway: pathway[2])

def wikipathways(gene, exclude_reactome=True):
    """Query WikiPathways for HGNC-Gene-ID"""
    reactome_reference_suffix = "reactome.org/PathwayBrowser/#DIAGRAM="
    wikipathways_api_url = "https://webservice.wikipathways.org"
    wikipathways_entry_url = "https://www.wikipathways.org/index.php/Pathway:"
    xml_ns = {"ns1": "http://www.wso2.org/php/xsd",
              "ns2": "http://www.wikipathways.org/webservice"}

    response = query(
        "{}/findPathwaysByText?query={}&species=Homo_sapiens".format(
            wikipathways_api_url, gene))

    if response is None:
        return []

    pathways = []
    for result in response.findall("ns1:result", xml_ns):
        pathway = {}
        pathway["id"] = result.findtext("ns2:id", namespaces=xml_ns)
        pathway["revision"] = result.findtext("ns2:revision",
                                              namespaces=xml_ns)
        pathway["name"] = result.findtext("ns2:name", namespaces=xml_ns)
        pathways.append(pathway)

    response = query("{}/findInteractions?query={}".format(
        wikipathways_api_url, gene))

    if response is None:
        return []

    interactions = {pathway["id"]: {"left": set(), "right": set()}
                    for pathway in pathways}

    for result in response.findall("ns1:result", xml_ns):
        species = result.findtext("ns2:species", namespaces=xml_ns)
        if species == "Homo sapiens":
            pathway_id = result.findtext("ns2:id", namespaces=xml_ns)
            for field in result.findall("ns2:fields", xml_ns):
                name = field.findtext("ns2:name", namespaces=xml_ns)
                if name in ("source", "indexerId"):
                    continue
                values = set()
                for value in field.findall("ns2:values", namespaces=xml_ns):
                    if value.text and len(value.text.split(" ")) == 1:
                        values.add(value.text.upper())
                if gene.upper() in values:
                    continue
                interactions[pathway_id][name] |= values

    for pathway_id in interactions:
        interactions[pathway_id]["left"] =\
                sorted([unicode(metabolite, "utf-8")
                        for metabolite in interactions[pathway_id]["left"]])
        interactions[pathway_id]["right"] =\
                sorted([unicode(metabolite, "utf-8")
                        for metabolite in interactions[pathway_id]["right"]])

    pathway_results = set((unicode(gene, "utf-8"),
                           unicode("WikiPathways", "utf-8"),
                           unicode(pathway["id"], "utf-8"),
                           unicode(pathway["revision"], "utf-8"),
                           unicode(pathway["name"], "utf-8"),
                           tuple(interactions[pathway["id"]]["left"]),
                           tuple(interactions[pathway["id"]]["right"]))
                          for pathway in pathways
                          if (interactions[pathway["id"]]["left"] or
                              interactions[pathway["id"]]["right"]))


    if not exclude_reactome:
        return sorted(list(pathway_results), key=lambda pathway: pathway[2])

    unique_pathways = set()
    for pathway in pathway_results:
        response = query("{}{}".format(wikipathways_entry_url, pathway[2]))

        if reactome_reference_suffix not in response:
            unique_pathways.add(pathway)

    return sorted(list(unique_pathways), key=lambda pathway: pathway[2])

def log_kegg_release():
    """Write queried KEGG release information to log file"""
    kegg_info_url = "http://rest.kegg.jp/info/pathway"

    response = query(kegg_info_url)

    release_info = response.split("\n")
    release = release_info[1].replace("path", "").strip()
    logging.info("KEGG version: %s", release)

def log_reactome_release():
    """Write queried Reactome release information to log file"""
    reactome_info_url =\
        "https://reactome.org/ContentService/data/database/version"

    release = query(reactome_info_url)
    logging.info("Reactome version: %s", release)

def build_pathway_db():
    """Build pathway database"""
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-c", "--config",
                        default=os.path.join(os.path.dirname(__file__),
                                             "../docs/codes3d.conf"),
                        help="The configuration file to be used.")

    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config)

    gene_id_fp = os.path.join(os.path.dirname(__file__),
                              config.get("Defaults", "gene_id_fp"))
    pathway_db_fp = os.path.join(os.path.dirname(__file__),
                                 config.get("Defaults", "pathway_db_fp"))
    pathway_db_tmp_fp = os.path.join(os.path.dirname(__file__),
                                     config.get("Defaults",
                                                "pathway_db_tmp_fp"))
    pathway_db_tab_fp = os.path.join(os.path.dirname(__file__),
                                     config.get("Defaults",
                                                "pathway_db_tab_fp"))
    pathway_db_log_fp = os.path.join(os.path.dirname(__file__),
                                     config.get("Defaults",
                                                "pathway_db_log_fp"))

    logging.getLogger("requests").setLevel(logging.WARNING)
    logging.basicConfig(
        filename=pathway_db_log_fp,
        level=logging.INFO,
        format="%(asctime)s\t%(levelname)s\t%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")

    logging.info("Build started.")
    logging.info("File name during build: %s", pathway_db_tmp_fp)
    logging.info("HGNC source: %s", gene_id_fp)
    log_kegg_release()
    log_reactome_release()

    with open(gene_id_fp) as gene_file:
        genes = gene_file.read().splitlines()

    gene_num_bar = progressbar.ProgressBar(max_value=len(genes))
    gene_num_bar.update(0)

    genes = sorted(dict.fromkeys(genes).keys(),
                   key=lambda gene: gene.capitalize())

    tab_file = open(pathway_db_tab_fp, "w")
    tab_writer = unicodecsv.writer(tab_file, delimiter="\t")

    pathway_db = sqlite3.connect(pathway_db_tmp_fp)
    pathway_db_cursor = pathway_db.cursor()
    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS genes
        """)
    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS names
        """)
    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS upstream
        """)
    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS downstream
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
        CREATE TABLE names (
            database TEXT,
            pathway TEXT,
            pathway_name TEXT,
            PRIMARY KEY (
                database,
                pathway)
            ON CONFLICT IGNORE)
        """)
    pathway_db_cursor.execute("""
        CREATE TABLE upstream (
            database TEXT,
            pathway TEXT,
            gene TEXT,
            upstream TEXT,
            PRIMARY KEY (
                database,
                pathway,
                gene,
                upstream)
            ON CONFLICT IGNORE)
        """)
    pathway_db_cursor.execute("""
        CREATE TABLE downstream (
            database TEXT,
            pathway TEXT,
            gene TEXT,
            downstream TEXT,
            PRIMARY KEY (
                database,
                pathway,
                gene,
                downstream)
            ON CONFLICT IGNORE)
        """)
    for i, gene in enumerate(genes):
        for database in (reactome, wikipathways):
            for pathway in database(gene):
                pathway_db_cursor.execute("""
                    INSERT INTO genes
                    VALUES (?, ?, ?)
                    """, (pathway[1], pathway[2], pathway[0]))

                pathway_db_cursor.execute("""
                    INSERT INTO names
                    VALUES (?, ?, ?)
                    """, (pathway[1], pathway[2], pathway[4]))

                for upstream in pathway[5]:
                    pathway_db_cursor.execute("""
                        INSERT INTO upstream
                        VALUES (?, ?, ?, ?)
                        """, (pathway[1], pathway[2], pathway[0], upstream))

                for downstream in pathway[6]:
                    pathway_db_cursor.execute("""
                        INSERT INTO downstream
                        VALUES (?, ?, ?, ?)
                        """, (pathway[1], pathway[2], pathway[0], downstream))

                upstream_genes = ", ".join(pathway[5])
                downstream_genes = ", ".join(pathway[6])
                tab_writer.writerow(pathway[:5] + (upstream_genes,
                                                   downstream_genes))

        gene_num_bar.update(i)

    pathway_db.commit()
    pathway_db.close()
    tab_file.close()

    os.rename(pathway_db_tmp_fp, pathway_db_fp)
    logging.info("Renamed %s to %s.", pathway_db_tmp_fp, pathway_db_fp)
    logging.info("Build completed.")


if __name__ == "__main__":
    build_pathway_db()

