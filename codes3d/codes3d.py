#!usr/bin/env python

from __future__ import print_function, division
import itertools
import argparse
import ast
import bisect
import configparser
import csv
import json
import logging
import multiprocessing
import math
import operator
import os
import sys
import pandas
import progressbar
import psutil
import pybedtools
import re
import requests
import shutil
import sqlite3
import time
import lxml.etree
import lxml.html
from operator import itemgetter
import Bio
from Bio import SeqIO
from Bio import Restriction
from Bio.Restriction import RestrictionBatch
from Bio.Restriction import Restriction_Dictionary
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib import style
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch
import matplotlib.gridspec as gridspec
import rpy2.robjects as R
import scipy.stats as st
import scipy.cluster as cl
import base64
import re
import numpy as np
import functools
import collections

def parse_parameters(restriction_enzymes, include_cell_line, exclude_cell_line):
    """Validate user parameters -r, -n and -x.

    Args:
        restriction_enzymes: space-delimited list of restriction enzymes from
            user. Limits program to Hic libraries restricted by specified enzyme.
        include_cell_line: space-delimited list of cell_lines from -n.
        exclude_cell_line: space-delimited list of cell_lines from -x
    Returns:
        res_enzymes: A list of restriction enzymes of which HiC libraries are
            included.
        include_cells: A list of validated HiC cell lines to be querried.
        exclude_cells: A list of validated HiC cell lines to be excluded.
    """
    print('Parsing HiC library restriction enzymes...')
    res_enzymes = []
    include_cells = []
    exclude_cells = []
    if restriction_enzymes == None:
        res_enzymes = HIC_RESTRICTION_ENZYMES
    else:
        for enzyme in restriction_enzymes:
            # Handle user input case issues
            enzymes_upper = [e.upper() for e in HIC_RESTRICTION_ENZYMES]
            if not enzyme.upper() in enzymes_upper:
                print('Warning: NO HiC libraries restricted with %s.' % enzyme)
            else:
                res_enzymes.append(HIC_RESTRICTION_ENZYMES[
                    enzymes_upper.index(enzyme.upper())])
    if not res_enzymes:
        sys.exit("Program terminating: No HiC libraries are prepared with " +
                 "your restriction enzyme(s).")
    cell_lines = []
    for enzyme in res_enzymes:
        cell_lines += os.listdir(os.path.join(lib_dir, enzyme, 'hic_data'))
    cell_lines_upper = [c.upper() for c in cell_lines]
    to_help = '1.\t' + cell_lines[0]
    if len(cell_lines) > 1:
        for i, cell in enumerate(cell_lines, start=1):
            to_help += '\n' + str(i+1) + '.\t' + cell
    if include_cell_line == None:
        include_cells = None
    else:
        for cell in include_cell_line:
            if cell.upper() in cell_lines_upper:
                include_cells.append(
                    cell_lines[cell_lines_upper.index(cell.upper())])
        if include_cells != None and len(include_cells) == 0:
            sys.exit("We cannot find one or more of the cell lines you " +
                     "passed to the -n option. Following is the list of " +
                     "cell lines we have: \n" + to_help)
    if exclude_cell_line == None:
        exclude_cells = None
    else:
        for cell in exclude_cell_line:
            if cell.upper() in cell_lines_upper:
                exclude_cells.append(
                    cell_lines[cell_lines_upper.index(cell.upper())])
        if exclude_cells != None and len(exclude_cells) == 0:
            sys.exit("We cannot find one or more of the cell lines you " +
                     "passed to the -x option. Following is the list of " +
                     "cell lines we have: \n" + to_help)

    return res_enzymes, include_cells, exclude_cells


def process_inputs(inputs, snp_database_fp, lib_dir,
                   restriction_enzymes, output_dir,
                   suppress_intermediate_files=False):
    """Finds restriction fragments in which SNPs lie.

    Args:
        inputs: File(s) (or stdin) containing SNP rsIDs or genomic positions
          in the format chr<x>:<locus>
        snp_database_fp: ../../lib/snp_index_dbSNP_b151.db
        # To point to fragment databases of different hic restrictions.
        lib_dir: ../lib
        output_dir: User-specified directory for results. Defaults to inputs directory.
        suppress_intermediate_files: if 'False', snps.txt file is written to output_dir

    Returns:
        A dict named snps with the ff structure:
            {'rs9462794':{
                          'chr': '6',
                          'locus': 12445846,
                          'fragments':[
                                       {'frag': 30968, 'enzyme': 'MboI'},
                                       {'frag': 3664, 'enzyme': 'HindIII'}
                                    ]
                        },
             'rs12198798':{...},
             'rs6909834':{...}
             }

        If suppress_intermediate_files=False, a snps.txt file is written
          with the ff columns:
            1. SNP rsID
            2. SNP chromosome
            3. SNP position
            4. Fragment ID
            5. Fragment restriction enzyme
    """
    print("Processing SNP input...")
    snp_db = sqlite3.connect(snp_database_fp)
    snp_db.text_factory = str
    snp_index = snp_db.cursor()
    snps = {}
    for inp in inputs:
        if os.path.isfile(inp):
            infile = open(inp, 'r')
            for line in infile:
                id = line.strip().split(' ')[0]
                snp = None
                if id.startswith('rs'):
                    snp_index.execute("SELECT * FROM snps WHERE rsID=?", (id,))
                else:
                    if not id.strip():
                        continue
                    chr = id[id.find("chr") + 3:id.find(':')]
                    locus = int(id[id.find(':') + 1:])
                    snp_index.execute(
                        "SELECT * FROM snps WHERE chr=? and locus=?", [chr, locus])
                snp = snp_index.fetchone()
                if snp is None:
                    print("Warning: %s does not exist in SNP database." % id)
                else:
                    # Query each fragment_bed_fp for SNP fragment.
                    for enzyme in restriction_enzymes:
                        fragment_database_fp = os.path.join(
                            lib_dir, enzyme, 'dna.fragments.db')
                        fragment_index_db = sqlite3.connect(
                            fragment_database_fp)
                        fragment_index_db.text_factory = str
                        fragment_index = fragment_index_db.cursor()
                        fragment_index.execute("SELECT fragment FROM fragments " +
                                               "WHERE chr=? AND start<=? AND end>=?",
                                               [snp[1], snp[2], snp[2]])
                        snp_fragment_result = fragment_index.fetchone()
                        if snp_fragment_result is None:
                            print("Warning: error retrieving SNP fragment for SNP " +
                                  snp[0])
                        else:
                            if not snp[0] in snps:
                                snps[snp[0]] = {"chr": '',
                                                "locus": "", "fragments": []}
                                snps[snp[0]]["chr"] = snp[1]
                                snps[snp[0]]["locus"] = snp[2]
                                snps[snp[0]]["fragments"].append({"frag": snp_fragment_result[0],
                                                                  "enzyme": enzyme})
                            else:
                                snps[snp[0]]["fragments"].append({"frag": snp_fragment_result[0],
                                                                  "enzyme": enzyme})
            infile.close()
        # TODO: input is a snp id, not a txt file
        else:
            snp = None
            if inp.startswith("rs"):
                snp_index.execute("SELECT * FROM snps WHERE rsID=?", (inp,))
            else:
                chr = inp[inp.find("chr") + 3:inp.find(':')]
                locus = int(inp[inp.find(':') + 1:])
                snp_index.execute("SELECT * FROM snps WHERE chr=? and locus=?",
                                  [chr, locus])
                snp = snp_index.fetchone()
                if snp is None:
                    print(
                        "Warning: %s does not exist in SNP database." %
                        inp)
                else:
                    # Query each fragment_bed_fp for SNP fragment.
                    for enzyme in restriction_enzymes:
                        fragment_database_fp = os.path.join(
                            lib_dir, enzyme, 'dna.fragments.db')
                        fragment_index_db = sqlite3.connect(
                            fragment_database_fp)
                        fragment_index_db.text_factory = str
                        fragment_index = fragment_index_db.cursor()
                        fragment_index.execute("SELECT fragment FROM fragments " +
                                               "WHERE chr=? AND start<=? AND end>=?",
                                               [snp[1], snp[2], snp[2]])
                        snp_fragment_result = fragment_index.fetchone()
                        if snp_fragment_result is None:
                            print("Warning: error retrieving SNP fragment for SNP " +
                                  snp[0])
                        else:
                            if not snp[0] in snps:
                                snps[snp[0]] = {"chr": '',
                                                "locus": "", "fragments": []}
                                snps[snp[0]]["chr"] = snp[1]
                                snps[snp[0]]["locus"] = snp[2]
                                snps[snp[0]]["fragments"].append({"frag": snp_fragment_result[0],
                                                                  "enzyme": enzyme})
                            else:
                                snps[snp[0]]["fragments"].append({"frag": snp_fragment_result[0],
                                                                  "enzyme": enzyme})
    if not suppress_intermediate_files:
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        snpfile = open(output_dir + "/snps.txt", 'w')
        swriter = csv.writer(snpfile, delimiter='\t')
        for snp, props in snps.items():
            for frag in props["fragments"]:
                swriter.writerow((snp, props["chr"], props["locus"],
                                  frag["frag"], frag["enzyme"]))
        snpfile.close()

    return snps


def find_interactions(
        snps,
        lib_dir,
        include,
        exclude,
        output_dir,
        num_processes,
        suppress_intermediate_files=False):
    """Finds fragments that interact with SNP fragments.

    Args:
        snps: The dictionary of SNP fragments returned from process_inputs
        hic_data_dir: ../../lib/hic_data
        include: A list of cell lines to query for interactions
        exclude: A list of cell lines to exclude from interaction query
        output_dir: User-specified directory for results. Defaults to inputs directory.
        suppress_intermediate_files: if 'False', snps.txt file is written to output_dir

    Returns:
        A dict named 'interactions' containing fragments interacting with SNPs
            in each cell line e.g.:
            {'MboI':{
                'rs9462794':{  # SNP rsID
                    'GM12878':{  # Cell line
                        ('6', '30099'): {
                            },
                        # Fragment chr and ID
                        ('6', '31188'): [replicate1, replicate2]
                        }
                    'HeLa':{
                        ('12', '19151'): [replicate1]
                        }
                    }
                'rs12198798':{..}
                'rs6909834':{...}
                }
            'HindIII':{...}
            }

        If suppress_intermediate_files=False, a snp-gene_interaction.txt file is written
          with the ff columns:
            1. SNP rsID
            2. Cell line
            3. Fragment chromosome
            4. Fragment ID
            5. Number of interactions
            6. Fragment restriction enzyme
    """
    print("Finding interactions...")
    # Look for all interactions involving SNP fragments in the HiC databases
    if include:
        include = set(include)
    if exclude:
        exclude = set(exclude)

    interactions = {}
    for enzyme in restriction_enzymes:
        interactions[enzyme] = {}
        for snp in snps:
            interactions[enzyme][snp] = {}

    for enzyme in restriction_enzymes:
        print("\tSearching HiC libraries restricted with " + enzyme)
        hic_data_dir = os.path.join(lib_dir, enzyme, 'hic_data')
        manager = multiprocessing.Manager()
        snp_interactions = manager.dict()
        current_process = psutil.Process()
        num_processes = 8
        pool = multiprocessing.Pool(processes=min(num_processes,
                                                  len(current_process.cpu_affinity())))
        arglist = [(cell_line, snp_interactions, snps, enzyme, lib_dir, include, exclude)
                   for cell_line in os.listdir(hic_data_dir)]
        pool.map(get_cell_interactions, arglist)
        pool.close()
        pool.join()
        for snp_cell in snp_interactions.keys():
            snp = snp_cell[0]
            cell_line = snp_cell[1]

            if cell_line not in interactions[enzyme][snp]:
                interactions[enzyme][snp][cell_line] = {}
            interactions[enzyme][snp][cell_line] = snp_interactions[snp_cell]
    if not suppress_intermediate_files:
        intfile = open(output_dir + "/snp-gene_interactions.txt", 'w')
        iwriter = csv.writer(intfile, delimiter='\t')
        for enzyme in interactions:
            for snp in interactions[enzyme]:
                for cell_line in interactions[enzyme][snp]:
                    for interaction in interactions[enzyme][snp][cell_line]:
                        iwriter.writerow(
                            (snp, cell_line, interaction[0], interaction[1],
                             len(interactions[enzyme][snp][cell_line][interaction]), enzyme))
        intfile.close()
    return interactions


def get_cell_interactions((cell_line, snp_interactions, snps, enzyme, lib_dir, include, exclude)):
    if (include and cell_line not in include) or (
            exclude and cell_line in exclude):
        return  # Enforce the -n or -x options
    hic_data_dir = os.path.join(lib_dir, enzyme, 'hic_data')
    cell_interactions = {}
    cell_interactions = {}
    print("\t\tSearching cell line " + cell_line)
    if os.path.isdir(os.path.join(hic_data_dir, cell_line)):
        # A set of unique interactions for SNP
        for snp in snps:
            cell_interactions[snp] = {}
            print("\t\t\tFinding interactions for " + snp)
            for fragment in snps[snp]['fragments']:
                if not fragment['enzyme'] == enzyme:
                    continue
                for replicate in os.listdir(hic_data_dir + '/' + cell_line):
                    if replicate.endswith(".db"):
                        rep_db = sqlite3.connect(
                            hic_data_dir + '/' + cell_line + '/' + replicate)
                        rep_db.text_factory = str
                        rep_ints = rep_db.cursor()
                        print("\t\t\t\tSearching replicate " + replicate)
                        # Search database for fragments interacting with SNP
                        # fragment
                        from_db = rep_ints.execute(
                            "SELECT chr2," +
                            "fragment2 FROM interactions WHERE chr1=?" +
                            "AND fragment1=?",
                            [
                                snps[snp]["chr"],
                                fragment["frag"]])
                        if from_db:
                            for interaction in from_db:
                                if interaction not in cell_interactions[snp]:
                                    cell_interactions[snp][interaction] = set(
                                        [])
                                cell_interactions[snp][interaction].add(
                                    replicate)
    for snp in cell_interactions:
        snp_interactions[snp, cell_line] = cell_interactions[snp]


def find_snp_genes(
        (snp, enzyme, snp_interactions, genes_dict, fragment_database_fp,
         gene_bed_fp, gene_dict_fp, output_dir)):
    print('\t\t{}'.format(snp))
    fragment_database_fp = os.path.join(lib_dir, enzyme, 'dna.fragments.db')
    fragment_index_db = sqlite3.connect(fragment_database_fp)
    fragment_index_db.text_factory = str
    fragment_index = fragment_index_db.cursor()
    gene_dict_db = sqlite3.connect(gene_dict_fp)
    gene_dict_db.text_factory = str
    gene_index = gene_dict_db.cursor()

    # Generate BED file of all fragments interacting with SNP-containing
    # fragment
    snp_genes = {}
    for cell_line in snp_interactions:
        snpgenes_exist = False
        temp_filepath = os.path.join(output_dir, "temp_snp_bed_"
                                     + snp + "_" + enzyme + ".bed")
        temp_snp_bed = open(temp_filepath, 'w')
        twriter = csv.writer(temp_snp_bed, delimiter='\t')
        num_reps = len([rep for rep in os.listdir(
            os.path.join(lib_dir, enzyme, 'hic_data', cell_line))
            if rep.endswith('.db')])
        for interaction in snp_interactions[cell_line]:
            fragment_index.execute(
                "SELECT start, end FROM fragments WHERE chr=? and fragment=?",
                [interaction[0], interaction[1]])
            fragment_pos = fragment_index.fetchall()
            if fragment_pos is None:
                print(
                    "\tWarning: error retrieving fragment %s on chromosome %s"
                    % (interaction[1], interaction[0]))
                continue
            for f in fragment_pos:
                twriter.writerow(("chr" + interaction[0], f[0], f[1],
                                  snp_interactions[cell_line][interaction]))
            if not snpgenes_exist:
                snpgenes_exist = True
        temp_snp_bed.close()
        if snpgenes_exist:
            int_bed = pybedtools.BedTool(temp_filepath)
            # Return a list of genes with which SNP is interacting
            # and the number of HiC contacts for each cell line.
            gene_bed = int_bed.intersect(hs_gene_bed, loj=True)
            for feat in gene_bed:
                gene_name = feat[7]
                if gene_name == '.' or feat[4] == '.' or \
                   feat[5] == '-1' or feat[6] == '-1':
                    # '.' indicates a NULL overlap
                    continue
                if not gene_name in snp_genes:
                    snp_genes[gene_name] = {}
                    gene_index.execute("SELECT gene_id FROM genes WHERE " +
                                       "gene_name = ?", (gene_name,))
                    gene_id = gene_index.fetchone()[0]
                    snp_genes[gene_name]['gene_id'] = gene_id
                    snp_genes[gene_name]['cell_lines'] = {}
                if not cell_line in snp_genes[gene_name]['cell_lines']:
                    snp_genes[gene_name]['cell_lines'][cell_line] = {}
                    snp_genes[gene_name]['cell_lines'][cell_line]['interactions'] = 0
                    snp_genes[gene_name]['cell_lines'][cell_line]['rep_present'] = []
                snp_genes[gene_name]['cell_lines'][cell_line]['interactions'] += 1
                snp_genes[gene_name]['cell_lines'][cell_line]['replicates'] = num_reps
                rep_present = feat[3].replace('}', '')
                rep_present = rep_present.replace('{', '')
                rep_present = rep_present.replace("'", "")
                rep_present = rep_present.split(',')
                rep_present = [e.strip() for e in rep_present]
                snp_genes[gene_name]['cell_lines'][cell_line]['rep_present'] += rep_present
        if os.path.isfile(temp_filepath):
            os.remove(temp_filepath)
    genes_dict[snp] = snp_genes
    fragment_index.close()
    fragment_index_db.close()
    gene_index.close()
    gene_dict_db.close()


def dict_merge(dct, merge_dct):
    """ Recursive dict merge. Inspired by :meth:``dict.update()``, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``merge_dct`` is merged into
    ``dct``.
    :param dct: dict onto which the merge is executed
    :param merge_dct: dct merged into dct
    :return: dct
    """
    for k, v in merge_dct.items():
        if (k in dct and isinstance(dct[k], dict)
                and isinstance(merge_dct[k], collections.Mapping)):
            dict_merge(dct[k], merge_dct[k])
        else:
            dct[k] = merge_dct[k]
    return dct


def find_genes(
        interactions,
        fragment_database_fp,
        gene_bed_fp,
        gene_dict_fp,
        output_dir,
        suppress_intermediate_files=False):
    """Identifies genes in fragments that interact with SNP fragments.

    Args:
        interactions: The dictionary fragements that interact with SNP fragments
            returned from find_interactions
        fragment_database_fp: ../../lib/Homo_sapiens.GRCh37.75.dna.fragments.db
        gene_bed_fp: ../../lib/gene_reference.bed
        gene_dict_fp: The database containing GENCODE ids for genes
        output_dir: User-specified directory for results. Defaults to inputs directory.
        suppress_intermediate_files: if 'False', snps.txt file is written to output_dir

    Returns:
        A dict named 'genes' containing genes in fragments that interact with SNP
                fragments e.g.
            {'rs9462794':{  # SNP rsID
                'PHACTR1':{  # Gene
                    'gene_id': 'ENSG00000112137.12',
                    'cell_lines':
                        'GM12878_Rao2014': 
                             {'interactions': 2, 
                              'replicates': 23, 
                              'rep_present': 
                                  ['GSM1551552_HIC003_merged_nodups.db', 
                                   'GSM1551553_HIC004_merged_nodups.db']},
                        'KBM7_Rao2014': 
                             {'interactions': 1, 
                              'replicates': 5, 
                              'rep_present': 
                                  ['GSM1551624_HIC075_merged_nodups.db']}
                     }
                'EDN1':{...}
                }
             'rs12198798':{..}
             'rs6909834':{...}
             }

        If suppress_intermediate_files=False, SNP-gene pairs with HiC contacts
          in only one libary replicate in only one cell line are written to
          genes_to_remove.txt.The others interactions are written to genes.txt.
          Each file has the ff columns:
            1. SNP rsID
            2. Gene name
            3. Gene ID
            4. Cell line
            5. HiC contact counts
            6. Replicates in which contact is found
            7. Number of cell line replicates
    """
    print("Identifying interactions with genes...")
    global hs_gene_bed
    hs_gene_bed = pybedtools.BedTool(gene_bed_fp)
    genes = {}
    enzyme_count = 0
    for enzyme in interactions:
        print("\tin libraries restricted with " + enzyme)
        enzyme_count += 1
        manager = multiprocessing.Manager()
        genes_dict = manager.dict()
        current_process = psutil.Process()
        num_processes = 8
        pool = multiprocessing.Pool(processes=min(num_processes,
                                                  len(current_process.cpu_affinity())))
        arglist = [(snp, enzyme, interactions[enzyme][snp], genes_dict,
                    fragment_database_fp, gene_bed_fp, gene_dict_fp, output_dir)
                   for snp in interactions[enzyme]]
        pool.map(find_snp_genes, arglist)
        pool.close()
        pool.join()
        if enzyme_count < 2:
            genes.update(genes_dict)
            genes.update()
        else:
            genes = dict_merge(genes, genes_dict)

    temp_files = [os.path.join(output_dir, temp_file)
                  for temp_file in os.listdir(output_dir)
                  if temp_file.startswith('temp_snp_bed')]
    for temp_file in temp_files:
        os.remove(temp_file)
    snps_to_remove = {}
    for enzyme in interactions:
        snps_to_remove[enzyme] = []
        for snp in interactions[enzyme]:
            if not snp in genes:
                print("\tNo SNP-gene spatial interactions detected for %s, \
                      removing from analysis" % (snp,))
                snps_to_remove[enzyme].append(snp)
    for enzyme in snps_to_remove:  # Update snps and interactions mappings
        for snp in snps_to_remove[enzyme]:
            for i, frag in enumerate(snps[snp]['fragments']):
                if frag['enzyme'] == enzyme:
                    snps[snp]['fragments'].remove(snps[snp]['fragments'][i])
            del interactions[enzyme][snp]
    genes_to_remove = []
    del_genefile = open(os.path.join(output_dir, 'genes_removed.txt'), 'w')
    dwriter = csv.writer(del_genefile, delimiter='\t')
    for snp in genes:
        for gene in genes[snp]:
            num_cell_line = len(genes[snp][gene]['cell_lines'])
            for cell_line in genes[snp][gene]['cell_lines']:
                rep_present = len(
                    set(genes[snp][gene]['cell_lines'][cell_line]['rep_present']))
                interactions = genes[snp][gene]['cell_lines'][cell_line]['interactions']
                replicates = genes[snp][gene]['cell_lines'][cell_line]['replicates']
                if interactions/replicates <= 1 and rep_present < 2 and\
                        num_cell_line < 2:
                    genes_to_remove.append((snp, gene))
                    dwriter.writerow((snp, gene, genes[snp][gene]['gene_id'],
                                      cell_line, interactions,
                                      rep_present, replicates))
    del_genefile.close()

    for pair in genes_to_remove:
        del genes[pair[0]][pair[1]]
    if not suppress_intermediate_files:
        genefile = open(output_dir + "/genes.txt", 'w')
        gwriter = csv.writer(genefile, delimiter='\t')
        for snp in genes:
            for gene in genes[snp]:
                num_cell_line = len(genes[snp][gene])
                for cell_line in genes[snp][gene]['cell_lines']:
                    rep_present = len(
                        set(genes[snp][gene]['cell_lines'][cell_line]['rep_present']))
                    interactions = genes[snp][gene]['cell_lines'][cell_line]['interactions']
                    replicates = genes[snp][gene]['cell_lines'][cell_line]['replicates']
                    gwriter.writerow((snp, gene, genes[snp][gene]['gene_id'],
                                      cell_line,
                                      genes[snp][gene]['cell_lines'][cell_line]['interactions'],
                                      genes[snp][gene]['cell_lines'][cell_line]['rep_present'],
                                      genes[snp][gene]['cell_lines'][cell_line]['replicates']))
    return genes


def find_eqtls(
        snps,
        genes,
        eqtl_data_dir,
        gene_database_fp,
        fdr_threshold,
        query_databases,
        num_processes,
        output_dir,
        gene_dict_fp,
        snp_dict_fp,
        suppress_intermediate_files=False):
    """Identifies genes in fragments that interact with SNP fragments.

    Args:
        snps: The dictionary of SNP fragments returned from process_inputs
        genes: The dictionary of genes that interact with SNPs from find_genes
        eqtl_data_dir: ../../eQTLs  #For local eQTL query
        gene_database_fp: ../../gene_reference.db
        fdr_threshold: Significance threshold for the FDR. Default = 0.05
        query_databases: a string specifying whether to query only local eQTL databases
           or only online databases or both.
        num_processes: Number of processors to be used when multiprocessing.
        output_dir: User-specified directory for results. Defaults to inputs directory.
        suppress_intermediate_files: if 'False', eqtls.txt file is written to output_dir

    Returns:
        num_tests: Number of eQTL interactions with adj_p_values >= FDR threshold
        p_values

        If suppress_intermediate_files=False, a eqtls.txt file is written
          with the ff columns:
            1. SNP rsID
            2. Gene name
            3. Tissue of eQTL interaction
            4. eQTL pvalue
            5. eQTL effect size
    """
    print("Identifying eQTLs of interest...")
    eqtls = {}  # A mapping of SNP-gene eQTL relationships
    p_values = []  # A list of all p-values for computing FDR
    num_tests = 0  # Total number of tests done
    global to_online_query  # List of SNP-gene-tissue to query
    to_online_query = []
    try:
        os.remove(os.path.join(output_dir, 'eqtls.txt'))
    except OSError:
        pass
    if query_databases == 'both' or query_databases == 'local':
        num_tests, to_online_query = query_local_databases(
            eqtl_data_dir,
            genes,
            gene_dict_fp,
            snp_dict_fp,
            p_values,
            num_processes,
            output_dir)
        if query_databases == 'local':
            print('Local tests done: {}'.format(
                num_tests))
        else:
            print('Local tests done: {}\t  Online tests to do: {}'.format(
                num_tests, len(to_online_query)))
    if query_databases == 'online':
        num_tests, to_online_query = prep_GTEx_queries(genes)
    if query_databases != 'local':
        num_processes = 1  # Only one process queries GTEx
        num_tests += query_GTEx_service(snps,
                                        genes,
                                        to_online_query,
                                        p_values,
                                        num_processes,
                                        output_dir)
    return num_tests, p_values


def query_local_databases(
        eqtl_data_dir, genes, gene_dict_fp, snp_dict_fp, p_values, num_processes, output_dir):
    print("\tQuerying local databases.")
    manager = multiprocessing.Manager()
    local_eqtls = manager.list()
    not_local = manager.list()
    num_tests = 0
    current_process = psutil.Process()
    num_processes = 16
    pool = multiprocessing.Pool(processes=min(num_processes,
                                              len(current_process.cpu_affinity())))
    pool.map(query_local_tissue,
             [(tissue_db, eqtl_data_dir, local_eqtls, not_local, genes, gene_dict_fp, snp_dict_fp)
              for tissue_db in sorted(os.listdir(eqtl_data_dir))])
    pool.close()
    pool.join()
    num_tests = len(local_eqtls)
    eqtlfile = open(os.path.join(output_dir, 'eqtls.txt'), 'a')
    ewriter = csv.writer(eqtlfile, delimiter='\t')
    ewriter.writerows(local_eqtls)
    eqtlfile.close()
    p_values += [eqtl[3] for eqtl in local_eqtls]
    return num_tests, not_local


def query_local_tissue((db, eqtl_data_dir, local_eqtls, not_local, genes, gene_dict_fp, snp_dict_fp)):
    tissue = db[:db.rfind('.')]
    print("\t\tQuerying " + tissue)
    gene_dict_db = sqlite3.connect(gene_dict_fp)
    gene_dict_db.text_factory = str
    gene_dict = gene_dict_db.cursor()
    snp_dict_db = sqlite3.connect(snp_dict_fp)
    snp_dict_db.text_factory = str
    snp_dict = snp_dict_db.cursor()
    eqtl_index_db = sqlite3.connect(os.path.join(eqtl_data_dir, db))
    eqtl_index_db.text_factory = str
    eqtl_index = eqtl_index_db.cursor()
    for snp in genes.keys():
        variant_id = ''
        snp_dict.execute("SELECT variant_id FROM lookup WHERE rsID = ?;",
                         (snp,))
        snp_fetch = snp_dict.fetchone()
        if snp_fetch:
            variant_id = snp_fetch[0]
        for gene in genes[snp]:
            gene_id = genes[snp][gene]['gene_id']
            if variant_id == '':
                # SNP is not in snp_dict
                to_online_query.append((snp, gene_id, tissue))
                continue
            eqtl_index.execute(
                "SELECT pval_nominal, slope FROM associations " +
                "WHERE variant_id=? AND gene_id=?",
                (variant_id, gene_id))
            eqtl_fetch = eqtl_index.fetchall()
            if not eqtl_fetch:
                # No SNP-gene eQTL association in tissue
                not_local.append((snp, gene_id, tissue))
                continue
            else:
                for eqtl in eqtl_fetch:
                    local_eqtls.append((snp, gene, tissue, eqtl[0], eqtl[1]))
    eqtl_index.close()
    eqtl_index_db.close()
    gene_dict.close()
    gene_dict_db.close()
    snp_dict.close()
    snp_dict_db.close()


def prep_GTEx_queries(genes):
    # TODO: Keep an eye to update tissues with every GTEx release
    tissues = ["Adipose_Subcutaneous",
               "Adipose_Visceral_Omentum",
               "Adrenal_Gland",
               "Artery_Aorta",
               "Artery_Coronary",
               "Artery_Tibial",
               "Brain_Amygdala",
               "Brain_Anterior_cingulate_cortex_BA24",
               "Brain_Caudate_basal_ganglia",
               "Brain_Cerebellar_Hemisphere",
               "Brain_Cerebellum",
               "Brain_Cortex",
               "Brain_Frontal_Cortex_BA9",
               "Brain_Hippocampus",
               "Brain_Hypothalamus",
               "Brain_Nucleus_accumbens_basal_ganglia",
               "Brain_Putamen_basal_ganglia",
               "Brain_Spinal_cord_cervical_c-1",
               "Brain_Substantia_nigra",
               "Breast_Mammary_Tissue",
               "Cells_EBV-transformed_lymphocytes",
               "Cells_Transformed_fibroblasts",
               "Colon_Sigmoid",
               "Colon_Transverse",
               "Esophagus_Gastroesophageal_Junction",
               "Esophagus_Mucosa",
               "Esophagus_Muscularis",
               "Heart_Atrial_Appendage",
               "Heart_Left_Ventricle",
               "Liver",
               "Lung",
               "Minor_Salivary_Gland",
               "Muscle_Skeletal",
               "Nerve_Tibial",
               "Ovary",
               "Pancreas",
               "Pituitary",
               "Prostate",
               "Skin_Not_Sun_Exposed_Suprapubic",
               "Skin_Sun_Exposed_Lower_leg",
               "Small_Intestine_Terminal_Ileum",
               "Spleen",
               "Stomach",
               "Testis",
               "Thyroid",
               "Uterus",
               "Vagina",
               "Whole_Blood"]
    to_online_query = []
    num = 0
    for snp in genes.keys():
        for gene in genes[snp]:
            for tissue in tissues:
                to_online_query.append(
                    (snp, genes[snp][gene]['gene_id'], tissue))
    return(num, to_online_query)


def query_GTEx_service(
        snps,
        genes,
        to_query,
        p_values,
        num_processes,
        output_dir):
    """Queries GTEx for eQTL association between SNP and gene.

    Args:
        snps: The dictionary of SNP fragments returned from process_inputs
        genes: The dictionary of genes that interact with SNPs from find_genes
        eqtls: A mapping of eQTl relationships between SNPs to genes, which
            further  maps to a list of tissues in which these eQTLs occur
        p_values: A list of all p-values for use computing FDR.
        num_processes: Number of processors to be used when multiprocessing.
        output_dir: User-specified directory for results. Defaults to inputs
            directory.

    Returns:
        num_tests: The number of successful GTEx responses
        eqtls and p_values are updated.
        If GTEx query fails for SNP-gene-tissue request, the failed requests
            are written to failed_GTEx_requests.txt thus:
            [
                {
                    'variantId': 'rs9462794',
                    'gencodeId': 'ENSG00000137872.1',
                    'tissueSiteDetailId': 'Artery_Aorta'
                }
            ]
    """

    if num_processes > 2:  # Ensure GTEx is not 'flooded' with requests
        num_processes = 2
    manager = multiprocessing.Manager()
    num_tests = 0
    reqLists = [[]]
    for assoc in to_query:
        if len(reqLists[-1]) < 960:  # >1000 requests is buggy
            reqLists[-1].append({"variantId": assoc[0],
                                 "gencodeId": assoc[1],
                                 "tissueSiteDetailId": assoc[2]})
        else:
            reqLists.append([{"variantId": assoc[0],
                              "gencodeId": assoc[1],
                              "tissueSiteDetailId": assoc[2]}])
    print("\tQuerying GTEx online service...")
    gtexResponses = manager.list()
    procPool = multiprocessing.Pool(processes=num_processes)
    for i, reqList in enumerate(reqLists, start=1):
        procPool.apply_async(
            send_GTEx_query, (i, len(reqLists), reqList, gtexResponses))
        time.sleep(10)
    procPool.close()
    procPool.join()
    print("\t\tGTEx responses received: " + str(len(gtexResponses)))
    print('Writing eQTL file...')
    eqtlfile = open(output_dir + '/eqtls.txt', 'a')
    ewriter = csv.writer(eqtlfile, delimiter='\t')
    failed_requests = []
    for response in gtexResponses:
        try:
            results = response[1].json()["result"]
            for result in results:
                if (str(result["geneSymbol"]) == "gene not found" or
                    not result["snpId"] in genes or
                        result["pValue"] == "N/A"):
                    continue
                num_tests += 1
                snp = result["snpId"]
                gene = result["geneSymbol"]
                p = float(result["pValue"])
                effect_size = float(result["nes"])
                p_values.append(p)
                ewriter.writerow((snp,
                                  gene,
                                  result["tissueSiteDetailId"],
                                  p,
                                  effect_size))
        except Exception as e:
            print("\t\tWARNING: bad response (%s)" % response[1])
            failed_requests.append(response[0])
    eqtlfile.close()
    print('\n')
    if failed_requests:
        # TODO: Handle the failed requests procedure
        with open(output_dir + "/failed_GTEx_requests.txt", 'w') \
                as failed_requests_file:
            failed_requests_file.write(str(failed_requests) + '\n')
    return num_tests


def send_GTEx_query(num, num_reqLists, reqList, gtexResponses):
    """Posts and receives requests from GTEx

    Args:
        num: An integer indicating the request being processed.
        num_reqLists: The total number of requests to be submitted
        reqList: A list of < 950 SNP, gene and tissue JSON data
        gtexResponses: A manager.list() object to handle requests.

    Returns:
        gtexResponses:
           [
              (ReqList:[]
               {'beta': '0.0832273122097661',
                'gencodeId': 'ENSG00000137872.11',
                'geneSymbol': 'SEMA6D',
                'pvalue': '0.37233196039320215',
                'pvalueThreshold': None,
                'se': '0.0929115299234',
                'variantId': 'rs9462794',
                'tissueId': 'Prostate',:
                'tstat': '0.8957694731579682',
                'variantId': '6_12445847_A_T_b37'})
           ]
    """
    s = requests.Session()
    s.verify = GTEX_CERT
    try:
        while True:
            print("\t\tSending request %s of %s" % (num, num_reqLists))
            gtex_url = "https://gtexportal.org/rest/v1/association/dyneqtl"
            res = s.post(gtex_url, json=reqList)
            if res.status_code == 200:
                gtexResponses.append((reqList, res))
                time.sleep(1)
                return
            elif res.status_code == 500 or res.status_code == 400:
                print("\t\tThere was an error processing request %s. \
                      Writing to failed request log and continuing." % num)
                gtexResponses.append((reqList, "Processing error"))
                time.sleep(2)
                return
            else:
                print("\t\tRequest %s received response with status %s. "
                      "Trying again in ten seconds." % (num, res.status_code))
                time.sleep(10)
    except requests.exceptions.ConnectionError:
        try:
            print("\t\tWarning: Request %s experienced a connection error. \
                  Retrying in one minute." % num)
            time.sleep(60)
            while True:
                print("\t\tSending request %s of %s" % (num, num_reqLists))
                res = s.post(
                    "https://gtexportal.org/rest/v1/association/dyneqtl",
                    json=reqList)
                if res.status_code == 200:
                    gtexResponses.append((reqList, res))
                    time.sleep(30)
                    return
                elif res.status_code == 500:
                    print("\t\tThere was an error processing request %s. \
                          Writing to failed request log and continuing." % num)
                    gtexResponses.append((reqList, "Processing error"))
                    time.sleep(60)
                    return
                else:
                    print("\t\tRequest %s received response with status: %s. \
                          Writing to failes request log and continuing."
                          % (num, res.status_code))
                    gtexResponses.append((reqList, res.status_code))
                    time.sleep(60)
                    return
        except requests.exceptions.ConnectionError:
            print("\t\tRetry failed. Continuing, but results will be incomplete.")
            gtexResponses.append((reqList, "Connection failure"))
            time.sleep(10)
            return
    s.close()


def get_gene_expression_information(eqtls, expression_table_fp, output_dir):
    print("Getting gene expression information...")
    gene_df = pandas.read_table(expression_table_fp, index_col='Symbol')
    gene_exp = pandas.DataFrame(data=None, columns=gene_df.columns)
    for snp in eqtls:
        for gene in eqtls[snp]:
            if gene not in gene_exp.index:
                try:
                    gene_exp = gene_exp.append(gene_df.ix[gene])
                except KeyError:
                    print("Warning: No gene expression information for %s"
                          % gene)
    gene_exp.to_csv(
        path_or_buf=output_dir +
        "/gene_expression_table.txt",
        sep='\t')


def get_expression(gene_tissue_tuple):
    try:
        expression = list(GENE_DF.at[gene_tissue_tuple[0],
                                     gene_tissue_tuple[1]])
        return (expression, max(expression))
    except TypeError:
        expression = list([GENE_DF.at[gene_tissue_tuple[0],
                                      gene_tissue_tuple[1]]])
        return (expression, expression[0])
    except KeyError:
        # No expression information for gene in tissue.
        return ([], 'NA')


def get_expression_extremes(gene_exp):
    for tissue in gene_exp.keys():
        gene_exp[tissue] = [x for x in gene_exp[tissue] if x > 0.0]
        if not gene_exp[tissue]:
            del gene_exp[tissue]
    if not gene_exp:
        return ('NA', 'NA', 'NA', 'NA')
    max_expression = max(gene_exp.items(), key=lambda exp: max(exp[1]))
    min_expression = min(gene_exp.items(), key=lambda exp: min(exp[1]))
    return (max_expression[0], max(max_expression[1]), min_expression[0],
            min(min_expression[1]))


def calc_hic_contacts(snp_gene_dict):
    """Calculates score of HiC contacts between SNP and gene.

    Args:
        snp_gene_dict: The snp-gene specific entry of the dictionary
        from find_genes

    Returns:
        cell_lines: string representing a list of cell line identifiers
        scores: string representing a list of the means of contacts in each
        cell line
        hic_score: sum of the averages of contacts per cell line.
    """
    hic_score = 0
    cell_lines = sorted(snp_gene_dict.keys())
    scores = []
    for cell_line in cell_lines:
        score = snp_gene_dict[cell_line]['interactions'] / float(
            snp_gene_dict[cell_line]['replicates'])
        scores.append('{:.2f}'.format(score))
        hic_score += score
    return (', '.join(cell_lines), ', '.join(scores),
            "{:.4f}".format(hic_score))


def produce_summary(
        p_values, snps, genes, gene_database_fp,
        expression_table_fp, fdr_threshold, do_not_produce_summary, output_dir,
        buffer_size, num_processes):
    """Write final results of eQTL-eGene associations

    Args:
        snps: The dictionary of SNP fragments returned from process_inputs
        genes: The dictionary of genes that interact with SNPs from find_genes
        gene_database_fp: ../../gene_reference.db
        expression_table_fp: ../../GTEx
        batch_ran: if 'False', do not produce summary, stop when eqtls.txt
        file is written to output_dir
        output_dir:

    Returns:
        summary.txt: A file with the ff structure.
            1. SNP
            2. SNP chromosome
            3. SNP locus
            4. Gene name
            5. Gene GENCODE ID
            6. Gene chromosome
            7. Gene start
            8. Gene end
            9. Tissue
            10. p-value
            11. Adj p-value
            12. Effect size
            13. Cell_Lines
            14. Cell_Lines_HiC_scores
            15. HiC_Contact_Score
            16. cis_SNP-gene_interaction
            17. SNP-gene_Distance
            18. Expression_Level_In_eQTL_Tissue
            19. Max_Expressed_Tissue
            20. Maximum_Expression_Level
            21. Min_Expressed_Tissue
            22. Min_Expression_Level

        sig_eqtls.text: A file with eQTL associations with
            adj_p_values <= FDR threshold. Same structure as summary.txt
    """
    if do_not_produce_summary:
        print("Without producing summary files... Done")
        return

    num_sig = {}
    # Number of eQTLs deemed significant under the given threshold
    gene_index_db = sqlite3.connect(gene_database_fp)
    gene_index_db.text_factory = str
    p_values_map = {}
    adjusted_p_values = compute_adj_pvalues(p_values)
    for i in range(len(p_values)):
        p_values_map[p_values[i]] = adjusted_p_values[i]

    to_file = []
    global gene_exp
    gene_exp = {}
    snps_genes = set()
    genes_tissues = set()
    eqtlfile = open(os.path.join(output_dir, 'eqtls.txt'), 'r',
                    buffering=buffer_size)
    ereader = csv.reader(eqtlfile, delimiter = '\t')

    for i, row in enumerate(ereader):
        snp = row[0]
        gene = row[1]
        tissue = row[2]
        p_val = row[3]
        effect_size = row[4]
        qvalue = p_values_map[float(p_val)]
        snp_info = {
            "chr": snps[snp]["chr"],
            "locus": snps[snp]["locus"],
            "fragments": snps[snp]["fragments"]}
        gene_chr = "NA"
        gene_start = 100000000000000000
        gene_end = -1
        snps_genes.add((snp, gene))
        genes_tissues.add((gene, tissue))
        gene_stat = gene_index_db.execute(
            'SELECT chr, start, end FROM genes WHERE symbol=?',
            (gene,)).fetchall()

        # Consider "canonical" to be the longest record where
        # multiple records are present
        for stat in gene_stat:
            gene_start = min(gene_start, int(stat[1]))
            gene_end = max(gene_end, int(stat[2]))
            gene_chr = stat[0]
        # eQTL is cis if the SNP is within 1Mbp of the gene on the same
        # chromosome
        cis = gene_chr == snp_info["chr"] and (
            (snp_info["locus"] > gene_start - 1000000) and
            (snp_info["locus"] < gene_end + 1000000))
        cell_lines = ''
        try:
            cell_lines = list(genes[snp][gene])
        except KeyError:
            cell_lines = "NA"
        interaction = True
        if cis:
            interaction = True
        else:
            interaction = False
        distance_from_snp = 0
        if(gene_chr == "NA"):
            distance_from_snp = "NA"  # If gene location information is missing
        elif(not snp_info["chr"] == gene_chr):
            distance_from_snp = "NA"  # Not applicable to trans interactions
        elif(snp_info["locus"] < gene_start):
            distance_from_snp = gene_start - snp_info["locus"]
        elif(snp_info["locus"] > gene_end):
            distance_from_snp = snp_info["locus"] - gene_end
        to_file.append([snp,
                        snp_info["chr"],
                        snp_info["locus"],
                        gene,
                        genes[snp][gene]['gene_id'],
                        gene_chr,
                        gene_start,
                        gene_end,
                        tissue,
                        p_val,
                        qvalue,
                        effect_size,
                        interaction,
                        distance_from_snp])
        if gene not in gene_exp:
            gene_exp[gene] = {}
    eqtlfile.close()

    global GENE_DF
    GENE_DF = pandas.read_table(expression_table_fp, index_col="Description",
                                engine='c', compression=None, memory_map=True)
    all_snps = genes.keys()
    all_genes = gene_exp.keys()
    all_tissues = list(GENE_DF)
    genes_tissues = list(genes_tissues)
    snps_genes = list(snps_genes)
    current_process = psutil.Process()
    pool = multiprocessing.Pool(processes=min(num_processes,
                                              len(current_process.cpu_affinity())))
    print("Collecting gene expression rates...")
    expression = pool.map(get_expression, genes_tissues)
    del GENE_DF
    for i in range(len(genes_tissues)):
        gene_exp[genes_tissues[i][0]][genes_tissues[i][1]] = expression[i][0]
    print("Determining gene expression extremes...")
    extremes = pool.map(get_expression_extremes,
                        [gene_exp[gene] for gene in all_genes])

    for i in range(len(all_genes)):
        gene_exp[all_genes[i]]['max_tissue'] = extremes[i][0]
        gene_exp[all_genes[i]]['max_rate'] = extremes[i][1]
        gene_exp[all_genes[i]]['min_tissue'] = extremes[i][2]
        gene_exp[all_genes[i]]['min_rate'] = extremes[i][3]

    for i in range(len(genes_tissues)):
        gene_exp[genes_tissues[i][0]][genes_tissues[i][1]] = expression[i][1]

    print("Computing HiC data...")
    hic_data = pool.map(calc_hic_contacts,
                        [genes[snp][gene]['cell_lines'] for snp, gene in snps_genes])
    global hic_dict
    hic_dict = {}
    for i in range(len(snps_genes)):
        hic_dict[snps_genes[i]] = hic_data[i]
    print("Writing to summary files...")
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    summary = open(output_dir + "/summary.txt", 'w', buffering=buffer_size)
    summ_writer = csv.writer(summary, delimiter='\t')
    summ_header = ['SNP',
                   'SNP_Chromosome',
                   'SNP_Locus',
                   'Gene_Name',
                   'GENCODE_Id',
                   'Gene_Chromosome',
                   'Gene_Start',
                   'Gene_End',
                   'Tissue',
                   'p-value',
                   'Adj_p-value',
                   'Effect_Size',
                   'Cell_Lines',
                   'Cell_Lines_HiC_scores',
                   'HiC_Contact_Score',
                   'cis_SNP-gene_interaction',
                   'SNP-gene_Distance',
                   'Expression_Level_In_eQTL_Tissue',
                   'Max_Expressed_Tissue',
                   'Maximum_Expression_Level',
                   'Min_Expressed_Tissue',
                   'Min_Expression_Level']
    summ_writer.writerow(summ_header)
    to_file = [insert_hic_info(line, line[0], line[3], line[8])
               for line in to_file]
    summ_writer.writerows(to_file)
    summary.close()
    sig_file = open(os.path.join(output_dir, 'significant_eqtls.txt'), 'w',
                    buffering=buffer_size)
    sigwriter = csv.writer(sig_file, delimiter='\t')
    sigwriter.writerow(summ_header)
    sig_eqtls = [line for line in to_file if line[10] < fdr_threshold]
    sigwriter.writerows(sig_eqtls)
    sig_file.close()

    return gene_exp.keys()


def insert_hic_info(line, snp, gene, tissue):
    try:
        line.insert(12, hic_dict[(snp, gene)][2])
        line.insert(12, hic_dict[(snp, gene)][1])
        line.insert(12, hic_dict[(snp, gene)][0])
        line.extend((
            round(gene_exp[gene][tissue], 2),
            gene_exp[gene]['max_tissue'],
            round(gene_exp[gene]['max_rate'], 2),
            gene_exp[gene]['min_tissue'],
            round(gene_exp[gene]['min_rate'], 2)))
        return line
    except TypeError:
        # Gene expression in tissue is likely NA.
        line.extend((
            gene_exp[gene][tissue],
            gene_exp[gene]['max_tissue'],
            gene_exp[gene]['max_rate'],
            gene_exp[gene]['min_tissue'],
            gene_exp[gene]['min_rate']))
        return line
    except KeyError:
        # Hi-C score for SNP-Gene is likely missing.
        line.insert(12, 'NA')
        line.insert(12, 'NA')
        line.insert(12, 'NA')
        line.extend((
            gene_exp[gene][tissue],
            gene_exp[gene]['max_tissue'],
            gene_exp[gene]['max_rate'],
            gene_exp[gene]['min_tissue'],
            gene_exp[gene]['min_rate']))
        return line


def compute_adj_pvalues(p_values):
    """ A Benjamini-Hochberg adjustment of p values of SNP-gene eQTL
           interactions from GTEx.

    Args:
        p_values: List of p values of all eQTL associations

    Returns:
        adj_pvalues: A corresponding list of adjusted p values to p_values.
    """
    p_values.sort()
    return R.r['p.adjust'](p_values, method='BH')


def produce_overview(genes, eqtls, num_sig, output_dir):
    # TODO: Make compatible without 'eqtls'
    """Generates overview graphs and table

    """
    print("\nProducing overview...")
    stat_table = open(output_dir + "/overview.txt", 'w')
    stat_table.write(
        "SNP\tChromosome\tLocus\tTotal_SNP-gene_Pairs\tTotal_eQTLs\n")
    for snp in eqtls:
        stat_table.write(snp +
                         '\t' +
                         eqtls[snp]["snp_info"]["chr"] +
                         '\t' +
                         str(eqtls[snp]["snp_info"]["locus"]) +
                         '\t' +
                         str(len(genes[snp])) +
                         '\t' +
                         str(len(eqtls[snp]) -
                             1) +
                         '\n')
    stat_table.close()
    print("\tProducing graphs")
    style.use("ggplot")
    if not os.path.isdir(output_dir + "/plots"):
        os.mkdir(output_dir + "/plots")
    int_colours = "rgb"
    eqtl_colours = "myc"
    int_colours = itertools.cycle(int_colours)
    eqtl_colours = itertools.cycle(eqtl_colours)
    snps_by_chr = {}
    for snp in eqtls:
        chrm = eqtls[snp]["snp_info"]["chr"]
        if not chrm in snps_by_chr:
            snps_by_chr[chrm] = []
        snps_by_chr[chrm].append((eqtls[snp]["snp_info"]["locus"], snp))

    int_colour_list = []
    eqtl_colour_list = []
    num_snpgenes = []
    num_eqtls = []
    rsIDs = []

    chrs = []
    for c in snps_by_chr:
        chrs.append(c)
    # So that the chromosomes are in a logical order on the graph
    chrs.sort(key=natural_keys)
    chr_locs = []
    chr_ticks = []
    chrm_pos = 0
    last_count = 0
    count = 0

    for chrm in chrs:
        snp_list = snps_by_chr[chrm]
        snp_list.sort()  # Sort by locus
        int_colour = next(int_colours)
        eqtl_colour = next(eqtl_colours)
        chr_locs.append(chrm_pos)
        chr_ticks.append(chrm)
        for snp in snp_list:
            num_snpgenes.append(len(genes[snp[1]]))
            num_eqtls.append(num_sig[snp[1]] * -1)
            rsIDs.append(snp[1])
            int_colour_list.append(int_colour)
            eqtl_colour_list.append(eqtl_colour)
            count += 1
        chrm_pos = chrm_pos + count
        count = 0

    plt.clf()
    plt.bar(range(len(num_snpgenes)), num_snpgenes, color=int_colour_list)
    plt.bar(range(len(num_eqtls)), num_eqtls, color=eqtl_colour_list)
    axes = plt.gca()
    # So that eQTL values won't appear negative on plot
    axes.yaxis.set_major_formatter(FuncFormatter(abs_value_ticks))
    plt.vlines(chr_locs, axes.get_ylim()[0], axes.get_ylim()[1], colors="k")
    axes.set_xticks(range(len(rsIDs)))
    axes.set_xticklabels(rsIDs, rotation='vertical')
    ax2 = axes.twiny()
    ax2.set_xlim(axes.get_xlim())
    ax2.set_xticks(chr_locs)
    ax2.set_xticklabels(chr_ticks)
    plt.tight_layout()
    plt.savefig(
        output_dir +
        "/plots/snpgene_and_eqtl_overview.png",
        dpi=300,
        format="png")
    plt.clf()


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [atoi(c) for c in re.split('(\d+)', text)]


def abs_value_ticks(x, pos):
    return abs(x)

def parse_snps_files(snps_files):
    snps = {}
    for snp_file in snps_files:
        with open(snp_file, 'r') as snpfile:
            for line in snpfile:
                snp = line.strip().split('\t')
                snps[snp[0]] = {"chr": snp[1],
                                "locus": int(snp[2]), "frag": snp[3]}
    return snps


def parse_interactions_files(interactions_files):
    interactions = {}
    for interactions_file in interactions_files:
        with open(interactions_file, 'r') as intfile:
            for line in intfile:
                interaction = line.strip().split('\t')
                snp = interaction[0]
                cell_line = interaction[1]
                if not snp in interactions:
                    interactions[snp] = {}
                if not cell_line in interactions[snp]:
                    interactions[snp][cell_line] = set([])
                interactions[snp][cell_line].add(
                    (interaction[2], int(interaction[3])))
    return interactions


def parse_genes_files(genes_files):
    genes = {}

    for genefile in genes_files.split(' '):
        genefile = open(genefile, 'rb')
        greader = csv.reader(genefile, delimiter='\t')
        for row in greader:
            gene = row[1]
            gene_id = row[2]
            snp = row[0]
            cell_line = row[3]
            interactions = int(row[4])
            interactions_replicates = row[5].split(', ')
            interactions_replicates = [
                w.replace('[', '') for w in interactions_replicates]
            interactions_replicates = [
                w.replace(']', '') for w in interactions_replicates]
            interactions_replicates = [
                w.replace("'", "") for w in interactions_replicates]
            replicates_count = int(row[6])
            try:
                genes[snp][gene]['cell_lines'][cell_line] = {
                    'interactions': interactions,
                    'replicates': replicates_count,
                    'rep_present': interactions_replicates}
            except KeyError:
                if not snp in genes:
                    genes[snp] = {}
                if not gene in genes[snp]:
                    genes[snp][gene] = {}
                    genes[snp][gene]['gene_id'] = gene_id
                    genes[snp][gene]['cell_lines'] = {}
                genes[snp][gene]['cell_lines'][cell_line] = {
                    'interactions': interactions,
                    'replicates': replicates_count,
                    'rep_present': interactions_replicates}
    return genes

### Function to remove duplicate lines from txt files
def remove_dups(inputfile, outputfile, buffer_size):
    lines = open(inputfile, 'r', buffering=buffer_size).readlines()
    lines_set = set(lines)
    out = open(outputfile, 'w')
    for line in lines_set:
        out.write(line)
    out.close()


def parse_eqtls_files(
        eqtls_files, snp_database_fp, gene_database_fp,
        restriction_enzymes, lib_fp, output_dir, buffer_size, fdr_threshold=0.05):
    print("Merging files...")
    eqtls = {}
    snps = {}
    genes = {}
    p_values = []  # A sorted list of all p-values for use computing FDR
    num_tests = 0  # Total number of tests done
    num_sig = {}  # Number of eQTLs deemed significant under the given threshold
    snp_dp = sqlite3.connect(snp_database_fp)
    snp_dp.text_factory = str
    snp_index = snp_dp.cursor()
    genes_files = []
    snps_files = []
    joined_eqtls_file = open(os.path.join(output_dir, 'tmp_eqtls.txt'), 'w')
    for eqtls_file in eqtls_files:
        file_path = eqtls_file[:len(eqtls_file)-9]
        # Check if genes.txt and snps.txt exist for each eqtls.txt
        if os.path.isfile(os.path.join(file_path, 'genes.txt')):
            genes_files.append(os.path.join(file_path, 'genes.txt'))
        else:
            print('\t\tOops!: We can\'t find  \'genes.txt\' file in the same' +
                  ' directory as \'eqtls.txt\'. \n\t\tSome info will be missing.')
        if os.path.isfile(os.path.join(file_path, 'snps.txt')):
            snps_files.append(os.path.join(file_path, 'snps.txt'))
        else:
            print('\t\tOops!: We can\'t find  \'snps.txt\' file in the same' +
                  ' directory as \'eqtls.txt\'. \n\t\tSome info will be missing.')
        with open(eqtls_file, 'r') as efile:
            shutil.copyfileobj(efile, joined_eqtls_file)

    joined_eqtls_file.close()

    # Removing duplicates from merged eqtls
    uniq_eqtls = os.path.join(output_dir, 'eqtls.txt')
    remove_dups(os.path.join(output_dir, 'tmp_eqtls.txt'), uniq_eqtls, buffer_size)
    os.remove(os.path.join(output_dir, 'tmp_eqtls.txt'))

    # Merging genes files
    joined_genes_file = open(os.path.join(output_dir, 'tmp_genes.txt'), 'w')
    for genes_file in genes_files:
        with open(genes_file, 'rb') as gfile:
            shutil.copyfileobj(gfile, joined_genes_file)
    joined_genes_file.close()

    # Removing duplicates from merged genes
    uniq_genes = os.path.join(output_dir, 'genes.txt')
    remove_dups(os.path.join(output_dir, 'tmp_genes.txt'), uniq_genes, buffer_size)
    os.remove(os.path.join(output_dir, 'tmp_genes.txt'))
    print("Parsing genes...")
    genes = parse_genes_files(os.path.join(output_dir, 'genes.txt'))

    # Merging snps files
    joined_snps_file = open(os.path.join(output_dir, 'tmp_snps.txt'), 'w')
    for snps_file in snps_files:
        with open(snps_file, 'r') as sfile:
            shutil.copyfileobj(sfile, joined_snps_file)
    joined_snps_file.close()

    # Removing duplicates from merged snps
    uniq_snps = os.path.join(output_dir, 'snps.txt')
    remove_dups(os.path.join(output_dir, 'tmp_snps.txt'), uniq_snps, buffer_size)
    os.remove(os.path.join(output_dir, 'tmp_snps.txt'))

    print("Parsing eqtls...")
    eqtlfile = open(os.path.join(output_dir, 'eqtls.txt'), 'r')
    ereader = csv.reader(eqtlfile, delimiter='\t')
    for i, row in enumerate(ereader, 1):
        snp = row[0]
        gene = row[1]
        tissue = row[2]
        effect_size = float(row[4])
        bisect.insort(p_values, float(row[3]))
        if not snp in snps:
            if os.path.isfile(os.path.join(output_dir, 'snps.txt')):
                snp_file = open(os.path.join(output_dir, 'snps.txt'), 'r')
                sreader = csv.reader(snp_file, delimiter='\t')
                for row in sreader:
                    if row[0] == snp:
                        try:
                            snps[snp]['fragments'].append({
                                'frag': row[3], 'enzyme': row[4]})
                        except KeyError:
                            snps[snp] = {}
                            snps[snp]['chr'] = row[1]
                            snps[snp]['locus'] = int(row[2])
                            snps[snp]['fragments'] = []
                            snps[snp]['fragments'].append({
                                'frag': row[3], 'enzyme': row[4]})
            else:  # In case there is no snps.txt file
                snp_index.execute("SELECT * FROM snps WHERE rsID=?", (snp,))
                from_db = snp_index.fetchone()
                snps[snp] = {'chr': from_db[1],
                             'locus': from_db[2],
                             'fragments': []}  # Handle fragments from multiple REs
                enzymes = restriction_enzymes.split(',')
                enzymes = [e.strip() for e in enzymes]
                for enzyme in enzymes:
                    conn = sqlite3.connect(os.path.join(
                        lib_fp, enzyme, 'dna.fragments.db'))
                    conn.text_factory = str
                    cur = conn.cursor()
                    cur.execute(
                        "SELECT fragment FROM fragments WHERE chr=? " +
                        "AND start<=? AND end>=?", (snps[snp]['chr'],
                                                    snps[snp]['locus'],
                                                    snps[snp]['locus']))
                    snps[snp]['fragments'].append({
                        'frag': cur.fetchone()[0],
                        'enzyme': enzyme})
                    cur.close()
                    conn.close()
    eqtlfile.close()
    snp_index.close()
    snp_dp.close()

    return (p_values, snps, genes)


def build_snp_index(
        snp_dir,
        output_fp,
        config,
        id_col=4,
        chr_col=1,
        locus_col=2,
        do_not_tidy_up=False):
    if not output_fp:
        output_fp = config["SNP_DATABASE_FP"]

    print("Building SNP index...")
    if not os.path.isdir(snp_dir):
        print("Error: argument to build SNP index must be a directory.")
        return

    if os.path.isfile(output_fp):
        upsert = input(
            "WARNING: Upserting input to existing SNP database (%s). Continue? [y/N] " %
            output_fp)
        if not upsert.lower() == 'y':
            print("Did not write to existing SNP database.")
            return

    if not os.path.isdir(os.path.dirname(output_fp)):
        os.makedirs(os.path.dirname(output_fp))

    snp_index_db = sqlite3.connect(output_fp)
    snp_index = snp_index_db.cursor()
    snp_index.execute(
        "CREATE TABLE IF NOT EXISTS snps (rsID text unique, chr text, locus integer)")
    snp_index.execute("CREATE INDEX IF NOT EXISTS id ON snps (rsID,chr,locus)")
    res = ""
    id_col = int(id_col) - 1
    chr_col = int(chr_col) - 1
    locus_col = int(locus_col) - 1
    for bed_fp in os.listdir(snp_dir):
        print("\tProcessing " + bed_fp)
        bed = open(snp_dir + '/' + bed_fp, 'r')
        for line in bed:
            if not line.startswith("chr"):
                continue
            snp = line.strip().split('\t')
            try:
                snp_index.execute("INSERT INTO snps VALUES(?,?,?)", [
                                  snp[id_col], snp[chr_col][snp[chr_col].find("chr") + 3:], int(snp[locus_col])])
            except sqlite3.IntegrityError:
                if res == "":
                    res = input(
                        "Warning: database already contains an entry for SNP with ID %s. Overwrite?\n1: Yes\n2: No (default)\n3: Yes To All\n4: No To All\nChoice: " %
                        snp[id_col])
                if res.strip() == "1":
                    print("Overwriting SNP %s" % snp[id_col])
                    snp_index.execute(
                        "DELETE FROM snps WHERE rsID=?", [
                            snp[id_col]])
                    snp_index.execute("INSERT INTO snps VALUES(?,?,?)", [
                                      snp[id_col], snp[chr_col][snp[chr_col].find("chr") + 3:], int(snp[locus_col])])
                    res = ""
                elif res.strip() == "3":
                    print("Overwriting SNP %s" % snp[id_col])
                    snp_index.execute(
                        "DELETE FROM snps WHERE rsID=?", [
                            snp[id_col]])
                    snp_index.execute("INSERT INTO snps VALUES(?,?,?)", [
                                      snp[id_col], snp[chr_col][snp[chr_col].find("chr") + 3:], int(snp[locus_col])])
                elif res.strip() == "4":
                    print("Skipping input SNP %s" % snp[id_col])
                    pass
                else:
                    print("Skipping input SNP %s" % snp[id_col])
                    res = ""
        bed.close()
    print("\tWriting SNP index to file...")
    snp_index_db.commit()
    print("Done building SNP index.")
    if not do_not_tidy_up:
        print("Tidying up...")
        shutil.rmtree(snp_dir)


def build_hic_index(
        input_hic_fp,
        output_fp=None,
        chr1_col=3,
        chr2_col=7,
        frag1_col=5,
        frag2_col=9,
        mapq1_col=10,
        mapq2_col=11,
        mapq_cutoff=30,
        do_not_tidy_up=False):
    if not output_fp:
        if not input_hic_fp.rfind('.') == -1:
            output_fp = input_hic_fp[:input_hic_fp.rfind('.')] + ".db"
        else:
            output_fp = input_hic_fp + ".db"

    if os.path.isfile(output_fp):
        overwrite = input(
            "WARNING: Overwriting existing HiC database (%s). Continue? [y/N] " %
            output_fp)
        if not upsert.lower() == 'y':
            print("Did not overwrite existing HiC database.")
            return
        os.remove(output_fp)

    # Do line count for progress meter
    print("Determining table size...")
    hic_table = open(input_hic_fp, 'r')
    lines = 0
    for i in hic_table:
        lines += 1
    hic_table.close()
    lines = lines // 100 * 100  # Get an approximation
    do_linecount = not lines == 0

    int_db = sqlite3.connect(output_fp)
    interactions = int_db.cursor()
    interactions.execute(
        "CREATE TABLE IF NOT EXISTS interactions (chr1 text, fragment1 text, chr2 text, fragment2 text)")
    interactions.execute(
        "CREATE INDEX IF NOT EXISTS i_1 ON interactions (chr1, fragment1)")
    chr1_col = int(chr1_col) - 1
    chr2_col = int(chr2_col) - 1
    frag1_col = int(frag1_col) - 1
    frag2_col = int(frag2_col) - 1
    mapq1_col = int(mapq1_col) - 1
    mapq2_col = int(mapq2_col) - 1
    with open(input_hic_fp, 'r') as rao_table:
        print("Indexing HiC interaction table...")
        for i, line in enumerate(rao_table):
            if do_linecount:
                if i % (lines / 100) == 0:
                    print("\tProcessed %d%%..." %
                          ((float(i) / float(lines)) * 100))
            interaction = line.strip().split(' ')
            if int(
                    interaction[mapq1_col]) >= mapq_cutoff and int(
                    interaction[mapq2_col]) >= mapq_cutoff:
                chr1 = interaction[chr1_col]
                frag1 = interaction[frag1_col]
                chr2 = interaction[chr2_col]
                frag2 = interaction[frag2_col]
                interactions.execute(
                    "INSERT INTO interactions VALUES(?,?,?,?)", [
                        chr1, frag1, chr2, frag2])
                interactions.execute(
                    "INSERT INTO interactions VALUES(?,?,?,?)", [
                        chr2, frag2, chr1, frag1])
    int_db.commit()
    interactions.close()
    print("Done indexing HiC interaction table.")
    if not do_not_tidy_up:
        print("Tidying up...")
        os.remove(input_hic_fp)


def digest_genome(
        genome,
        restriction_enzyme,
        output_fp,
        output_db,
        do_not_index=False,
        linear=False):
    if not output_fp:
        if not genome.rfind('.') == -1:
            output_fp = "%s.%s.fragments.bed" % (
                genome[:genome.rfind('.')], restriction_enzyme)
        else:
            output_fp = "%s.%s.fragments.bed" % (genome, restriction_enzyme)
    if not os.path.isdir(os.path.dirname(output_fp)):
        os.makedirs(os.path.dirname(output_fp))
    if os.path.isfile(output_fp):
        overwrite = input(
            "WARNING: Overwriting existing fragment BED %s. Continue? [y/N] " %
            output_fp)
        if not overwrite.lower() == 'y':
            print("Did not overwrite existing fragment BED.")
            return
        os.remove(output_fp)

    print("Digesting")
    if "fasta" in genome or "fa" in genome:
        genome = Bio.SeqIO.parse(open(genome, "rU"), format='fasta')
    else:
        genome = Bio.SeqIO.parse(open(genome, "rU"), format='genbank')

    with open(output_fp, 'w') as bedfile:
        output = csv.writer(bedfile, delimiter="\t")
        for chromosome in genome:
            print(chromosome.id, chromosome.name)
            # Digest the sequence data and return the cut points
            fragment_num = 0
            enzyme = RestrictionBatch([restriction_enzyme])
            for enzyme, cutpoints in enzyme.search(
                    chromosome.seq, linear=linear).items():
                cutpoints.insert(0, 0)
                # Covert cut points to fragments
                for index, point in enumerate(cutpoints):
                    # Adjust for start and end of sequence and the offset for the cutpoint
                    # ATTENTION: This will only work for MspI I will have to spend some time to make it compatible with
                    #		   any restriction enzyme such as ones with blunt ends
                    #startpoint = 0 if cutpoints[index] - 1 < 0 else cutpoints[index] - 1
                    startpoint = 0 if cutpoints[index] - \
                        1 < 0 else cutpoints[index]
                    endpoint = len(chromosome.seq) if index + \
                        1 >= len(cutpoints) else cutpoints[index + 1] - 1

                    accession = chromosome.id
                    version = ''
                    if "." in chromosome.id:
                        accession, version = chromosome.id.split(".")
                    if not accession.startswith("chr"):
                        accession = "chr" + accession
                    output.writerow(
                        [accession, startpoint, endpoint, fragment_num])
                    fragment_num += 1
    if not do_not_index:
        build_fragment_index(output_fp, output_db)


def build_fragment_index(fragment_fp, output_db):
    if not output_db:
        if not fragment_fp.rfind('.') == -1:
            output_db = fragment_fp[:fragment_fp.rfind('.')] + ".db"
        else:
            output_db = fragment_fp + ".db"
    if os.path.isfile(output_db):
        overwrite = input(
            "WARNING: Overwriting existing fragment database %s. Continue? [y/N] " %
            output_db)
        if not overwrite.lower() == 'y':
            print("Did not overwrite existing fragment database.")
            return
        os.remove(output_db)

    if not os.path.isdir(os.path.dirname(output_db)):
        os.path.makedirs(os.path.dirname(output_db))

    fragment_index_db = sqlite3.connect(output_db)
    fragment_index = fragment_index_db.cursor()
    fragment_index.execute(
        "CREATE TABLE fragments (chr text, start integer, end integer, fragment integer)")
    fragment_index.execute("CREATE INDEX f_index ON fragments (chr,fragment)")

    with open(fragment_fp, 'r') as fragments_bed:
        for line in fragments_bed:
            fragment = line.strip().split('\t')
            fragment_index.execute("INSERT INTO fragments VALUES (?,?,?,?)", [
                                   fragment[0][fragment[0].find("chr") + 3:], int(fragment[1]), int(fragment[2]), fragment[3]])
    fragment_index_db.commit()


def build_gene_index(
        gene_files,
        output_bed,
        output_db,
        config,
        symbol_col=27,
        chr_col=29,
        start_col=30,
        end_col=31,
        p_thresh_col=None,
        no_header=False,
        do_not_tidy_up=False):
    genes = {}

    # Edited to cater for when only the .gtf and .conf arguments are provided.
    if symbol_col:
        symbol_col -= 1
    if chr_col:
        chr_col -= 1
    if start_col:
        start_col -= 1
    if end_col:
        end_col -= 1
    if p_thresh_col:
        p_thresh_col -= 1

    if not output_bed:
        output_bed = os.path.join(
            config.get(
                "Defaults",
                "LIB_DIR"),
            "gene_reference.bed")

    if not output_db:  # corrected output_bed to output_db
        output_db = os.path.join(
            config.get(
                "Defaults",
                "LIB_DIR"),
            "gene_reference.db")

    append_bed = False
    overwrite_bed = True
    if os.path.isfile(output_bed):
        upsert = input(
            "WARNING: Appending input to existing BED file (%s). Continue? [y/N] " %
            output_bed)
        if not upsert.lower() == 'y':
            print("Did not append to existing gene database.")
            overwrite_bed = False
        else:
            append_bed = True

    if (overwrite_bed or append_bed) and not os.path.isdir(
            os.path.dirname(output_bed)):
        os.makedirs(os.path.dirname(output_bed))

    upsert_db = True

    if os.path.isfile(output_db):
        upsert = input(
            "WARNING: Upserting input to existing gene database (%s). Continue? [y/N] " %
            output_db)
        if not upsert.lower() == 'y':
            print("Did not write to existing gene database.")
            upsert_db = False

    if upsert_db and not os.path.isdir(os.path.dirname(output_db)):
        os.makedirs(os.path.dirname(output_db))

    if (not (overwrite_bed or append_bed)) and not upsert_db:
        print("No action performed; exiting.")
        return

    if upsert_db:
        gene_index_db = sqlite3.connect(output_db)
        gene_index = gene_index_db.cursor()
        if p_thresh_col:
            gene_index.execute(
                "CREATE TABLE IF NOT EXISTS genes (symbol text, chr text, start integer, end integer, double p_thresh)")
        else:
            gene_index.execute(
                "CREATE TABLE IF NOT EXISTS genes (symbol text, chr text, start integer, end integer)")
        gene_index.execute(
            "CREATE INDEX IF NOT EXISTS g_index ON genes (symbol)")

    for gene_file in gene_files:
        # Do line count for progress meter
        lines = 0
        with open(gene_file, 'r') as genefile:
            print("Determining table size...")
            for i in genefile:
                lines += 1
            lines = lines // 100 * 100  # Get an approximation
            do_linecount = not lines == 0

        with open(gene_file, 'r') as genefile:
            # Determine if input file is a GTEx-supplied gene reference
            is_gtex_file = gene_file.endswith(".gtf")
            # If so, process headers accordingly
            if is_gtex_file:
                line = genefile.readline()
                while line.startswith("##"):
                    line = genefile.readline()
            # Otherwise, process headers according to script options
            else:
                if not no_header:
                    genefile.readline()
                line = genefile.readline()
            i = 0
            # For each line
            while line:
                if do_linecount:
                    if i % (lines / 100) == 0:
                        print("\tProcessed %d%%..." %
                              ((float(i) / float(lines)) * 100))
                # If the line is a GTEx file, extract information accordingly
                if is_gtex_file:
                    entry = line.strip().split('\t')
                    if entry[2] == "gene":
                        gene_stats = entry[8].split(';')
                        gene_symbol = gene_stats[4].strip().split(' ')[
                            1].strip('"')
                        gene_chr = entry[0]
                        gene_start = int(entry[3])
                        gene_end = int(entry[4])
                    else:
                        i += 1
                        line = genefile.readline()
                        continue  # Skip if the entry is not for a canonical gene
                # Otherwise, extract by column
                else:
                    gene = line.strip().split('\t')
                    gene_symbol = gene[symbol_col]
                    if gene[chr_col].startswith("chr"):
                        gene_chr = gene[chr_col][3:]
                    else:
                        gene_chr = gene[chr_col]
                    gene_start = int(gene[start_col])
                    gene_end = int(gene[end_col])
                    if p_thresh_col:
                        try:
                            gene_p_thresh = float(gene[p_thresh_col])
                        except ValueError:
                            gene_p_thresh = gene[p_thresh_col]
                # Enter into index, regardless of input file type
                if not gene_symbol in genes:
                    genes[gene_symbol] = {
                        "chr": gene_chr, "start": gene_start, "end": gene_end}
                    if p_thresh_col:
                        genes[gene_symbol]["p_thresh"] = gene_p_thresh
                else:
                    curr_length = genes[gene_symbol]["end"] - \
                        genes[gene_symbol]["start"]
                    if curr_length < abs(gene_end - gene_start):
                        genes[gene_symbol]["start"] = gene_start
                        genes[gene_symbol]["end"] = gene_end
                line = genefile.readline()
                i += 1

    bed_out = None
    if overwrite_bed:
        bed_out = open(output_bed, 'w')
    elif append_bed:
        bed_out = open(output_bed, 'a')
    for gene in genes:
        if bed_out:
            bed_out.write(
                'chr%s\t%s\t%s\t%s\n' %
                (genes[gene]["chr"],
                 genes[gene]["start"],
                    genes[gene]["end"],
                    gene))
        if upsert_db:
            if p_thresh_col:
                gene_index.execute(
                    "INSERT INTO genes VALUES (?,?,?,?,?)", [
                        gene, genes[gene]["chr"], genes[gene]["start"], genes[gene]["end"], genes[gene]["p_thresh"]])
            else:
                gene_index.execute(
                    "INSERT INTO genes VALUES (?,?,?,?)", [
                        gene, genes[gene]["chr"], genes[gene]["start"], genes[gene]["end"]])
    if bed_out:
        bed_out.close()
    if upsert_db:
        gene_index_db.commit()
        gene_index.close()

    if not do_not_tidy_up:
        print("Tidying up...")
        for gene_file in gene_files:
            os.remove(gene_file)


def build_eqtl_index(
        table_fp,
        output_fp=None,
        snp_col=23,
        gene_symbol_col=27,
        gene_chr_col=29,
        gene_start_col=30,
        gene_stop_col=31,
        p_val_col=6,
        effect_size_col=3,
        do_not_tidy_up=False):
    if not output_fp:
        if not table_fp.rfind('.') == -1:
            output_fp = table_fp[:table_fp.rfind('.')] + ".db"
        else:
            output_fp = table_fp + ".db"

    if not os.path.isdir(os.path.dirname(output_fp)):
        os.makedirs(os.path.dirname(output_fp))

    if os.path.isfile(output_fp):
        upsert = input(
            "WARNING: Upserting input to existing eQTL database %s. Continue? [y/N] " %
            output_fp)
        if not upsert.lower() == 'y':
            print("Did not write to existing eQTL database.")
            return

    snp_col -= 1
    gene_symbol_col -= 1
    gene_chr_col -= 1
    gene_start_col -= 1
    gene_stop_col -= 1
    p_val_col -= 1
    effect_size_col -= 1
    table_index_db = sqlite3.connect(output_fp)
    table_index = table_index_db.cursor()
    table_index.execute(
        "CREATE TABLE IF NOT EXISTS eqtls (rsID text, gene_name text, gene_chr text, gene_start integer, gene_end integer, pvalue real, effect_size real)")
    table_index.execute("CREATE INDEX IF NOT EXISTS id ON eqtls (rsID)")
    # Do line count for progress meter
    do_linecount = True
    print("Determining table size...")
    with open(table_fp, 'r') as eqtl_table:
        lines = 0
        for i in eqtl_table:
            lines += 1
    lines = lines // 100 * 100  # Get an approximation
    do_linecount = not lines == 0

    with open(table_fp, 'r') as eqtl_table:
        for i, line in enumerate(eqtl_table):
            if do_linecount:
                if i % (lines / 100) == 0:
                    print("\tProcessed %d%%..." %
                          ((float(i) / float(lines)) * 100))
            if i == 0:
                continue
            eqtl = line.strip().split('\t')
            table_index.execute(
                "INSERT INTO eqtls VALUES (?,?,?,?,?,?,?)",
                [
                    eqtl[snp_col],
                    eqtl[gene_symbol_col],
                    eqtl[gene_chr_col],
                    eqtl[gene_start_col],
                    eqtl[gene_stop_col],
                    eqtl[p_val_col],
                    eqtl[effect_size_col]])
    table_index_db.commit()
    print("Done indexing eQTL table.")
    if not do_not_tidy_up:
        print("Tidying up...")
        os.remove(table_fp)


def retry_request(url, times=10, pause=60):
    """"""
    for i in range(times):
        time.sleep(pause)
        try:
            retried_response = requests.get(url)
        except requests.exceptions.ConnectionError:
            logging.error("Connection failure: %s", url)
            continue
        except requests.exceptions.RequestException:
            logging.error("Error: %s", url)
            continue
        try:
            retried_response.raise_for_status()
            return retried_response
        except requests.exceptions.HTTPError:
            logging.error("Unsucessful status code: %s: %s",
                    retried_response.status_code, url)
            continue
    else:
        return None

def request(url):
    """"""
    while True:
        try:
            response = requests.get(url)
            break
        except requests.exceptions.ConnectionError:
            logging.error("Connection failure: %s", url)
            time.sleep(60)
            continue
        except requests.exceptions.RequestException:
            logging.error("Error: %s", url)
            return None

    try:
        response.raise_for_status()
    except requests.exceptions.HTTPError:
        if response.status_code in (502,):
            response = retry_request(url)
            if response is None:
                logging.error("Unsucessful status code: %s: %s",
                               response.status_code, url)
                return None

        else:
            logging.error("Unsucessful status code: %s: %s",
                          response.status_code, url)
            return None

    content_type = response.headers["Content-Type"].split("; ")[0]

    if content_type == "application/json":
        try:
            return response.json()
        except ValueError:
            logging.error("Invalid JSON content: %s", url)
            return None

    elif content_type == "application/xml":
        try:
            content = lxml.etree.fromstring(response.text)
        except lxml.etree.LxmlSyntaxError:
            for i in range(10):
                retried_response = retry_request(url, times=1)
                if retried_response is not None:
                    try:
                        content = lxml.etree.fromstring(retried_response.text)
                        break
                    except lxml.etree.LxmlSyntaxError:
                        continue
            else:
                logging.error("Invalid XML content: %s", url)
                return None

        reg = "^https://webservice.wikipathways.org/" +\
              "getPathwayAs\?fileType=gpml&pwId=WP[0-9]+&revision=0$"

        if re.match(reg, url):
            try:
                return lxml.etree.fromstring(base64.b64decode(
                            content.find("ns1:data", content.nsmap).text))

            except TypeError, lxml.etree.LxmlSyntaxError:
                logging.error("Invalid GPML content: %s", url)
                return None
        else:
            return content

    elif content_type == "text/xml":
        reg = "^http://rest.kegg.jp/get/hsa[0-9]+/kgml$"
        kgml = re.match(reg, url)

        try:
            content = lxml.etree.fromstring(response.text)
        except lxml.etree.XMLSyntaxError:
            if kgml:
                logging.error("Invalid KGML content: %s", url)
            else:
                logging.error("Invalid XML content: %s", url)
            content = None

        if content is not None and kgml:
            try:
                header = content.xpath("/comment()[1]")[0].text
            except IndexError:
                header = None
            return header, content
        elif content is not None:
            return content
        elif kgml:
            return None, None
        else:
            return None

    elif content_type == "text/html":
        try:
            return lxml.html.fromstring(response.text)
        except ValueError:
            logging.error("Invalid HTML content: %s", url)
            return None

    elif content_type == "text/plain":
        return response.text

    else:
        logging.error("Unsupported content type %s: %s",
                url, response.headers["Content-Type"])
        return None

def kegg(gene):
    """"""
    kegg_api_url = "http://rest.kegg.jp"
    response = request("{}/find/genes/{}".format(kegg_api_url, gene))

    if response is None:
        return []

    gene_entries = response.split("\n")
    gene_id = None

    for entry in [entry.split("\t") for entry in gene_entries if entry]:
        if (entry[0].split(":")[0] == "hsa" and
               gene in entry[1].split(";")[0].split(",")):
            gene_id = entry[0]
            break

    if gene_id is None:
        return []

    pathways = set()
    response = request("{}/link/pathway/{}".format(kegg_api_url, gene_id))

    if response is None:
        return []

    pathway_entries = response.split("\n")
    pathway_ids = [entry.split("\t")[1].split(":")[1]
                   for entry in pathway_entries if entry]

    pathways = {}
    for pathway_id in pathway_ids:
        pathways[pathway_id] = {}
        response = request("{}/get/{}".format(kegg_api_url, pathway_id))

        if response is None:
            pathways[pathway_id]["name"] = "NA"
            continue

        pathway_entries = response.split("\n")

        for entry in pathway_entries:
            if entry.startswith("NAME"):
                pathway_name = entry.replace("NAME", "").replace(
                    "- Homo sapiens (human)", "").strip()
                pathways[pathway_id]["name"] = pathway_name
                break

    for pathway_id in pathway_ids:
        response = request("{}/get/{}/kgml".format(kegg_api_url, pathway_id))

        if response is None:
            pathways[pathway_id]["version"] = "NA"
            pathways[pathway_id]["upstream"] = set()
            pathways[pathway_id]["downstream"] = set()

        else:
            header_comment, response = response

            entry_ids = set()
            for entry in response.findall("entry"):
                if (entry.get("type") == "gene" and
                        gene_id in entry.get("name").split(" ")):
                    entry_ids.add(entry.get("id"))

            upstream_entry_ids = set()
            downstream_entry_ids = set()
            for relation in response.findall("relation"):
                if relation.get("type") == "PPrel":
                    if relation.get("entry1") in entry_ids:
                        downstream_entry_ids.add(relation.get("entry2"))

                    if relation.get("entry2") in entry_ids:
                        upstream_entry_ids.add(relation.get("entry1"))

            upstream_gene_ids = set()
            downstream_gene_ids = set()
            for entry in response.findall("entry"):
                if (entry.get("id") in upstream_entry_ids and
                        entry.get("type") == "gene"):
                    upstream_gene_ids |= set(entry.get("name").split(" "))
                if (entry.get("id") in downstream_entry_ids and
                        entry.get("type") == "gene"):
                    downstream_gene_ids |= set(entry.get("name").split(" "))

            upstream_gene_ids.discard("undefined")
            downstream_gene_ids.discard("undefined")

            if upstream_gene_ids or downstream_gene_ids:
                query_genes = list(upstream_gene_ids | downstream_gene_ids)

                # Limitation of URL length avoiding denial of request
                gene_subsets = [query_genes[i:i+100]
                                for i in xrange(0, len(query_genes), 100)]

                gene_name = {}
                for gene_subset in gene_subsets:
                    genes_str = "+".join(gene_subset)
                    response = request("{}/list/{}".format(kegg_api_url,
                                    genes_str), "text/plain")

                    if response is None:
                        continue

                    for entry in [entry.split("\t")
                                for entry in response[:-1].split("\n")]:
                        if ";" in entry[1]:
                            gene_name[entry[0]] =\
                                    entry[1].split(";")[0].split(",")[0]
                        else:
                            upstream_gene_ids.discard(entry[0])
                            downstream_gene_ids.discard(entry[0])

                upstream_genes = set(gene_name[gene_id]
                                    for gene_id in upstream_gene_ids)
                downstream_genes = set(gene_name[gene_id]
                                    for gene_id in downstream_gene_ids)

                upstream_genes.discard(gene)
                downstream_genes.discard(gene)

            else:
                upstream_genes, downstream_genes = set(), set()

            pathways[pathway_id]["upstream"] = upstream_genes
            pathways[pathway_id]["downstream"] = downstream_genes

            months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                      "Aug", "Sep", "Oct", "Nov", "Dec"]

            if (header_comment and
                header_comment[:16] == " Creation date: " and
                len(header_comment.split(" ")) > 5 and
                header_comment.split(" ")[3] in months):
                header_comment = header_comment.split(" ")
                month = str(months.index(header_comment[3]) + 1)
                day = header_comment[4].replace(",", "")
                year = header_comment[5]
                pathway_version = "-".join([year, month, day])
                pathways[pathway_id]["version"] = pathway_version
            else:
                pathways[pathway_id]["version"] = "NA"

    pathways = [[unicode(gene, "utf-8"),
                 unicode(pathway_id),
                 unicode(pathways[pathway_id]["version"]),
                 unicode(pathways[pathway_id]["name"]),
                 tuple(sorted(list(pathways[pathway_id]["upstream"]))),
                 tuple(sorted(list(pathways[pathway_id]["downstream"])))]
                 for pathway_id in pathway_ids
                 if pathways[pathway_id]["upstream"] or
                    pathways[pathway_id]["downstream"]]

    return sorted(pathways, key=lambda pathway: pathway[1])

class GeneName():
    """"""
    def __init__(self, pathway_db_gene_name_fp):
        self.integer = re.compile("^[0-9]+$")

        self.synonym = {}
        with open(pathway_db_gene_name_fp) as gene_syn_file:
            gene_syn_reader = csv.reader(gene_syn_file, delimiter="\t")
            next(gene_syn_reader)

            for row in gene_syn_reader:
                if row[1] != "-":
                    for syn in row[1].split("|"):
                        self.synonym[syn.lower()] = row[0]
                self.synonym[row[2].lower()] = row[0]

    def correct(self, gene_name):
        gene_name = gene_name.replace("\n"," ").replace("  ", " ")
        gene_name = gene_name.strip("\"'()?*")
        gene_name = gene_name.split("(")[0].strip()

        if self.integer.match(gene_name):
            return None

        return self.synonym.get(gene_name.lower(), gene_name)

def reactome(gene):
    """"""
    reactome_api_url = "https://www.reactome.org/ContentService"

    response = request(
        "{}/data/mapping/UniProt/{}/pathways?species=9606".format(
            reactome_api_url, gene))

    if response is None:
        return []

    pathways = {}
    for entry in response:
        pathway_id = entry["stId"]
        pathways[pathway_id] = {}
        pathways[pathway_id]["name"] = entry["displayName"]
        pathways[pathway_id]["version"] = entry["stIdVersion"].split(".")[1]


    response = request(
            "{}/data/mapping/UniProt/{}/reactions?species=9606".format(
                reactome_api_url, gene))

    gene_reaction_ids = set()
    if response is not None:
        for entry in response:
            gene_reaction_ids.add(entry["stId"])

    reactions_to_pathways = {}
    for pathway_id in pathways:
        response = request("{}/data/pathway/{}/containedEvents/stId".format(
            reactome_api_url, pathway_id))

        if response is None:
            reactions_to_pathways[pathway_id] = set()
            continue

        response_values = set(response.strip("[]").split(", "))
        pathway_reaction_ids = response_values & gene_reaction_ids
        reactions_to_pathways[pathway_id] = pathway_reaction_ids

    reaction_ids = set(itertools.chain.from_iterable(
                                reactions_to_pathways.values()))
    reactions = {}
    for reaction_id in reaction_ids:
        reactions[reaction_id] = {}
        reactions[reaction_id]["input"] = set()
        reactions[reaction_id]["output"] = set()

        response = request("{}/data/query/{}".format(
            reactome_api_url, reaction_id))

        if response is None:
            continue

        if "input" in response:
            for metabolite in response["input"]:
                if (isinstance(metabolite, dict) and
                        metabolite["className"] == "Protein"):
                    name = metabolite["displayName"]
                    name = name.split(" [")[0].split("(")[0]
                    reactions[reaction_id]["input"].add(
                            gene_name.correct(name))

        if "output" in response:
            for metabolite in response["output"]:
                if (isinstance(metabolite, dict) and
                        metabolite["className"] == "Protein"):
                    name = metabolite["displayName"]
                    name = name.split(" [")[0].split("(")[0]
                    reactions[reaction_id]["output"].add(
                            gene_name.correct(name))

    for pathway_id in pathways:
        pathways[pathway_id]["input"] = set()
        pathways[pathway_id]["output"] = set()
        for reaction_id in reactions_to_pathways[pathway_id]:
            if gene in reactions[reaction_id]["input"]:
                pathways[pathway_id]["output"] |=\
                        reactions[reaction_id]["output"]
            if gene in reactions[reaction_id]["output"]:
                pathways[pathway_id]["input"] |=\
                        reactions[reaction_id]["input"]

        pathways[pathway_id]["input"].discard(gene)
        pathways[pathway_id]["output"].discard(gene)


    pathways = [[unicode(gene, "utf-8"),
                 pathway_id,
                 pathways[pathway_id]["version"],
                 pathways[pathway_id]["name"],
                 tuple(sorted(list(pathways[pathway_id]["input"]))),
                 tuple(sorted(list(pathways[pathway_id]["output"])))]
                 for pathway_id in pathways
                 if pathways[pathway_id]["input"] or
                    pathways[pathway_id]["output"]]

    return sorted(pathways, key=lambda pathway: pathway[1])

def wikipathways(gene,
                 exclude_reactome=False,
                 exclude_homology_mappings=True):
    """"""
    wikipathways_api_url = "https://webservice.wikipathways.org"
    discard_tags = {"Curation:Hypothetical",
                    "Curation:NoInteractions",
                    "Curation:Tutorial",
                    "Curation:UnderConstruction"}


    if exclude_reactome:
        discard_tags.add("Curation:Reactome_Approved")

    response = request(
        "{}/findPathwaysByText?query={}&species=Homo_sapiens".format(
            wikipathways_api_url, gene))

    if response is None:
        return []

    pathways = []
    for result in response.findall("ns1:result", response.nsmap):
        pathway = {}
        pathway["id"] = result.find("ns2:id", response.nsmap).text
        pathway["revision"] = result.find("ns2:revision", response.nsmap).text
        pathway["name"] = result.find("ns2:name", response.nsmap).text

        pathway["input"] = {}
        pathway["output"] = {}

        pathways.append(pathway)

    for pathway in pathways[:]:
        response = request("{}/getCurationTags?pwId={}".format(
                           wikipathways_api_url, pathway["id"]))

        if response is None:
            continue

        for tag in response.findall("ns1:tags", response.nsmap):
            if tag.find("ns2:name", response.nsmap).text in discard_tags:
                pathways.remove(pathway)
                break

    for pathway in pathways[:]:
        response = request(
                "{}/getPathwayAs?fileType=gpml&pwId={}&revision=0".format(
                wikipathways_api_url, pathway["id"]))

        if response is None:
            continue

        nsmap = {}
        nsmap["gpml"] = response.nsmap[None]

        homology_mapped_pathway = False
        if exclude_homology_mappings:
            for comment in response.findall("gpml:Comment", nsmap):
                if comment.get("Source") == "HomologyMapper":
                    homology = True
                    break
        if homology_mapped_pathway:
            pathways.remove(pathway)
            continue

        group_to_graph, graph_to_group = {}, {}
        for group in response.findall("gpml:Group", nsmap):
            group_id = group.get("GroupId")
            graph_id = group.get("GraphId")
            if graph_id:
                group_to_graph[group_id] = graph_id
                graph_to_group[graph_id] = group_id

        graph_ids = set()
        corrected_name = {}

        for data_node in response.findall("gpml:DataNode", nsmap):
            text_label = data_node.attrib["TextLabel"]

            if isinstance(text_label, str):
                corrected_name[text_label] = gene_name.correct(text_label)

                if corrected_name[text_label] == gene:
                    graph_id = data_node.get("GraphId")
                    if graph_id:
                        graph_ids.add(graph_id)

                    group_ref = data_node.get("GroupRef")
                    if group_ref in group_to_graph:
                        graph_ids.add(group_to_graph[group_ref])

        input_ids = {"non-group": {}, "group": {}}
        output_ids = {"non-group": {}, "group": {}}

        for interaction in response.findall("gpml:Interaction", nsmap):
            points = interaction.find("gpml:Graphics", nsmap).findall(
                    "gpml:Point", nsmap)

            input_ref = points[0].get("GraphRef")
            output_ref = points[-1].get("GraphRef")
            regulation_type = points[-1].get("ArrowHead")

            if input_ref in graph_ids:
                if output_ref in graph_to_group:
                    output_ref = graph_to_group[output_ref]
                    output_ids["group"][output_ref] = regulation_type
                else:
                    output_ids["non-group"][output_ref] = regulation_type

            if output_ref in graph_ids:
                if input_ref in graph_to_group:
                    input_ref = graph_to_group[input_ref]
                    input_ids["group"][input_ref] = regulation_type
                else:
                    input_ids["non-group"][input_ref] = regulation_type

        if (input_ids["non-group"] or
            input_ids["group"] or
            output_ids["non-group"] or
            output_ids["group"]):

            for data_node in response.findall("gpml:DataNode", nsmap):
                if data_node.get("Type") == "GeneProduct":
                    xref = data_node.find("gpml:Xref", nsmap)
                    if not (xref.get("Database") and xref.get("ID")):
                        continue

                    graph_id = data_node.get("GraphId")
                    group_ref = data_node.get("GroupRef")

                    if graph_id or group_ref:
                        text_label = data_node.attrib["TextLabel"]
                        name = corrected_name.get(text_label)

                        if not name or name == gene:
                            continue

                        name = unicode(name, "utf-8")

                        if graph_id in input_ids["non-group"]:
                            if input_ids["non-group"][graph_id] == "Arrow":
                                pathway["input"][name] = [True, False]
                            elif input_ids["non-group"][graph_id] == "TBar":
                                pathway["input"][name] = [False, True]
                            else:
                                pathway["input"][name] = [False, False]

                        if group_ref in input_ids["group"]:
                            if input_ids["group"][group_ref] == "Arrow":
                                pathway["input"][name] = [True, False]
                            elif input_ids["group"][group_ref] == "TBar":
                                pathway["input"][name] = [False, True]
                            else:
                                pathway["input"][name] = [False, False]

                        if graph_id in output_ids["non-group"]:
                            if output_ids["non-group"][graph_id] == "Arrow":
                                pathway["output"][name] = [True, False]
                            elif output_ids["non-group"][graph_id] == "TBar":
                                pathway["output"][name] = [False, True]
                            else:
                                pathway["output"][name] = [False, False]

                        if group_ref in output_ids["group"]:
                            if output_ids["group"][group_ref] == "Arrow":
                                pathway["output"][name] = [True, False]
                            elif output_ids["group"][group_ref] == "TBar":
                                pathway["output"][name] = [False, True]
                            else:
                                pathway["output"][name] = [False, False]

    pathways = [[unicode(pathway["id"], "utf-8"),
                 unicode(pathway["revision"], "utf-8"),
                 unicode(pathway["name"], "utf-8"),
                 pathway["input"],
                 pathway["output"]]
                 for pathway in pathways]

    return sorted(pathways, key=lambda pathway: pathway[0])

def log_kegg_release():
    """"""
    kegg_info_url = "http://rest.kegg.jp/info/pathway"
    response = request(kegg_info_url)

    if response is None:
        logging.info("KEGG version: NA")
    else:
        response = response.split("\n")
        release = response[1].replace("path", "").strip()
        release = release.split("/")[0].split(" ")[1]
        release = release.strip("+")
        logging.info("KEGG version: %s", release)

    return None

def log_reactome_release():
    """"""
    reactome_info_url =\
        "https://reactome.org/ContentService/data/database/version"

    release = request(reactome_info_url)
    if not release:
        logging.info("Reactome version: NA")
    else:
        logging.info("Reactome version: %s", release)

    return None

def get_wikipathways_release():
    """"""
    wikipathways_info_url = "http://data.wikipathways.org/"
    response = request(wikipathways_info_url)

    if response is None:
        return None
    else:
        table = response.find("body").find("div").find("table")
        for tr in table.find("tbody").findall("tr"):
            if tr.find("td").find("a").text == "current":
                break
            release = tr.find("td").find("a").text

        release = "-".join([release[0:4], release[4:6], release[6:8]])
        return release

def get_exp(args):
    """"""
    gene, tissue = args

    try:
       exp = GENE_DF.at[gene, tissue]
    except KeyError:
        return None

    if isinstance(exp, np.float64):
        return exp
    elif isinstance(exp, np.ndarray):
        return exp.max()
    else:
        return None

def normalize_distribution(exp):
    """"""
    mean = sum([exp[tissue] for tissue in exp])/len(exp)
    sdv = math.sqrt(sum([(exp[tissue] - mean)**2
                    for tissue in exp])/(len(exp) - 1))

    if sdv == 0.0:
        exp_z = {tissue: 0.0 for tissue in exp}
    else:
        exp_z = {tissue: (exp[tissue] - mean) / sdv for tissue in exp}

    exp_p = {tissue: 1 - st.norm.cdf(exp_z[tissue])
             if exp_z[tissue] > 0.0 else st.norm.cdf(exp_z[tissue])
             for tissue in exp}

    exp = {tissue: (exp[tissue], exp_z[tissue], exp_p[tissue])
           for tissue in exp}

    return exp


def build_expression_tables(
        num_processes,
        gtex_gene_exp_fp,
        hpm_map_fp,
        hpm_gene_exp_fp,
        hpm_protein_exp_fp,
        hpm,
        db_fp,
        gene_exp_tsv_fp,
        protein_exp_tsv_fp,
        genes):
    """"""
    gene_tsv_file = open(gene_exp_tsv_fp, "w", buffering=1)
    gene_tsv_writer = csv.writer(gene_tsv_file, delimiter="\t",
                                 lineterminator="\n")
    if hpm:
        protein_tsv_file = open(protein_exp_tsv_fp, "w", buffering=1)
        protein_tsv_writer = csv.writer(protein_tsv_file, delimiter="\t",
                                 lineterminator="\n")
    tsv_header = ["Gene",
                  "Tissue",
                  "Expression",
                  "z-Score",
                  "p-Value"]
    gene_tsv_writer.writerow(tsv_header)

    if hpm:
        protein_tsv_writer.writerow(tsv_header)

    pathway_db = sqlite3.connect(db_fp)
    pathway_db.text_factory = str
    pathway_db_cursor = pathway_db.cursor()

    upstream_genes = set([upstream_gene for (upstream_gene,)
                          in pathway_db_cursor.execute("""
                              SELECT
                                  upstream_gene
                              FROM
                                  UpstreamGenes""")])
    genes = set(genes) | upstream_genes

    downstream_genes = set([downstream_gene for (downstream_gene,)
                            in pathway_db_cursor.execute("""
                              SELECT
                                  downstream_gene
                              FROM
                                  DownstreamGenes""")])
    genes |= downstream_genes

    print("Number of genes including up- and downstream genes: {}".format(
          len(genes)), end=2*"\n")


    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS
            GeneExpression
        """)
    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS
            ProteinExpression
        """)

    pathway_db_cursor.execute("""
        CREATE TABLE
            GeneExpression (
                gene TEXT,
                tissue TEXT,
                gene_exp REAL,
                gene_exp_z REAL,
                gene_exp_p REAL,
                PRIMARY KEY (
                    gene,
                    tissue)
                ON CONFLICT IGNORE)
            """)

    if hpm:
        pathway_db_cursor.execute("""
            CREATE TABLE
                ProteinExpression (
                    gene TEXT,
                    tissue TEXT,
                    protein_exp REAL,
                    protein_exp_z REAL,
                    protein_exp_p REAL,
                    PRIMARY KEY (
                        gene,
                        tissue)
                    ON CONFLICT IGNORE)
                """)

    if genes:
        global GENE_DF
        GENE_DF = pandas.read_table(gtex_gene_exp_fp,
                                    index_col="Description", engine='c',
                                    compression=None, memory_map=True)

        current_process = psutil.Process()
        num_processes = min(num_processes, len(current_process.cpu_affinity()))
        pool = multiprocessing.Pool(processes=num_processes)

        gtex_tissues = list(GENE_DF)

        exp = pool.map(get_exp,
                [(gene, tissue) for gene in genes for tissue in gtex_tissues])

        gene_exp = {gene: {tissue: exp.pop(0)
                           for tissue in gtex_tissues}
                    for gene in genes}

        del GENE_DF
        pool.close()

        for gene in gene_exp.keys():
            for tissue in gene_exp[gene].keys():
                if gene_exp[gene][tissue] is None:
                    del gene_exp[gene][tissue]
            if not gene_exp[gene]:
                del gene_exp[gene]

        if hpm:
            accessions = {}
            with open(hpm_map_fp, "r") as gene_map_file:
                gene_map_reader = csv.reader(gene_map_file)
                next(gene_map_reader)

                for line in gene_map_reader:
                    for gene in line[3].split(";"):
                        if gene in genes:

                            if gene not in accessions:
                                accessions[gene] = set()

                            accessions[gene] |= set([accession.split("|")[1]
                                            for accession in line[2].split(";")])

            accession_exp = {}
            with open(hpm_protein_exp_fp, "r") as protein_exp_file:
                protein_exp_reader = csv.reader(protein_exp_file)
                header = next(protein_exp_reader)

                hpm_tissues = {i: header[i].replace("Adult ", "").replace(" ", "_")
                            for i in range(len(header))
                            if header[i][:6] == "Adult "}

                for line in protein_exp_reader:
                    accession = line[0]
                    accession_exp[accession] = {}
                    for i in hpm_tissues:
                        tissue = hpm_tissues[i]
                        accession_exp[accession][tissue] = float(line[i])

            protein_exp = {}
            for gene in accessions:
                protein_exp[gene] = {tissue: sum([accession_exp[accession][tissue]
                                    for accession in accessions[gene]])
                                    for tissue in accession_exp[accession]}

            with open(hpm_gene_exp_fp, "r") as gene_exp_file:
                gene_exp_reader = csv.reader(gene_exp_file)
                header = next(gene_exp_reader)

                hpm_tissues = {i: header[i].replace("Adult ", "").replace(" ", "_")
                            for i in range(len(header))
                            if header[i][:6] == "Adult "}

                for line in gene_exp_reader:
                    gene = line[0]
                    if gene in genes:
                        if gene not in protein_exp:
                            protein_exp[gene] = {}

                        for i in hpm_tissues:
                            tissue = hpm_tissues[i]
                            if float(line[i]) == 0.0:
                                protein_exp[gene][tissue] = 0.0

            exp = merge_gtex_hpm(gene_exp, protein_exp)

        else:
            exp = {"gene": gene_exp}

        gene_coverage = " ({:.2f}%)".format((len(exp["gene"])/len(genes))*100)

        if hpm:
            protein_coverage = " ({:.2f}%)".format((len(exp["protein"])/\
                                                        len(genes))*100)
    else:
        gene_coverage = ""
        exp = {"gene":{}}

        if hpm:
            protein_coverage = ""
            exp["protein"] = {}

    if exp["gene"]:
        print("Mapping gene expression information...")
        gene_num = progressbar.ProgressBar(max_value=len(exp["gene"]))
        gene_num.update(0)

        for i, gene in enumerate(sorted(exp["gene"]), 1):
            exp["gene"][gene] = normalize_distribution(exp["gene"][gene])
            for tissue in sorted(exp["gene"][gene]):
                pathway_db_cursor.execute("""
                    INSERT INTO
                        GeneExpression
                    VALUES
                        (?, ?, ?, ?, ?)
                    """,
                    (unicode(gene, "utf-8"),
                    unicode(tissue, "utf-8"),
                    exp["gene"][gene][tissue][0],
                    exp["gene"][gene][tissue][1],
                    exp["gene"][gene][tissue][2]))

                gene_tsv_writer.writerow([
                    gene,
                    tissue,
                    exp["gene"][gene][tissue][0],
                    exp["gene"][gene][tissue][1],
                    exp["gene"][gene][tissue][2]])

            pathway_db.commit()
            gene_num.update(i)

        gene_num.finish()
    gene_tsv_file.close()

    print("Number of genes covered by gene expression data: {}{}".format(
          len(exp["gene"]), gene_coverage), end=2*"\n")

    if hpm and exp["protein"]:
        print("Mapping protein expression information...")
        gene_num = progressbar.ProgressBar(max_value=len(exp["protein"]))
        gene_num.update(0)

        for i, gene in enumerate(sorted(exp["protein"]), 1):
            exp["protein"][gene] = normalize_distribution(exp["protein"][gene])
            for tissue in sorted(exp["protein"][gene]):
                pathway_db_cursor.execute("""
                    INSERT INTO
                        ProteinExpression
                    VALUES
                        (?, ?, ?, ?, ?)
                    """,
                    (unicode(gene, "utf-8"),
                    unicode(tissue, "utf-8"),
                    exp["protein"][gene][tissue][0],
                    exp["protein"][gene][tissue][1],
                    exp["protein"][gene][tissue][2]))

                protein_tsv_writer.writerow([
                    gene,
                    tissue,
                    exp["protein"][gene][tissue][0],
                    exp["protein"][gene][tissue][1],
                    exp["protein"][gene][tissue][2]])

            pathway_db.commit()
            gene_num.update(i)

        gene_num.finish()

    if hpm:
        protein_tsv_file.close()

    if hpm:
        print("Number of genes covered by protein expression data: {}{}".format(
            len(exp["protein"]), protein_coverage), end=2*"\n")

    pathway_db.close()

    return None

def parse_input_genes(
        buffer_size,
        gene_synonym_map_fp,
        input_gene_map_fp,
        input_genes,
        input_files,
        summary_files):

    """"""
    global gene_name
    gene_name = GeneName(gene_synonym_map_fp)

    if input_genes is None:
        input_genes = []

    if input_files:
        for input_file in input_files:
            with open(input_file, "r",
                      buffering=buffer_size) as input_gene_file:
                input_genes.extend([
                    gene.strip() for gene in input_gene_file])

    if summary_files:
        for summary_file in summary_files:
            with open(summary_file, "r", buffering=buffer_size) as summary:
                summary_reader = csv.reader(summary, delimiter="\t")
                next(summary_reader)
                input_genes.extend([line[3] for line in summary_reader])

    input_genes = [input_genes[i]
                   for i in range(len(input_genes))
                   if input_genes[i] and input_genes[i] not in input_genes[:i]]

    print("Number of unique input genes: {}".format(len(input_genes)))

    with open(input_gene_map_fp, "w", buffering=buffer_size) as gene_map:
        gene_map_writer = csv.writer(gene_map, delimiter="\t",
                                     lineterminator="\n")
        gene_map_header = [
            "Input_Gene_Name",
            "Remapped_Gene_Name"]
        gene_map_writer.writerow(gene_map_header)

        genes = []
        for i in range(len(input_genes)):
            corrected_name = gene_name.correct(input_genes[i])

            if corrected_name:

                if input_genes[i] != corrected_name:
                    gene_map_writer.writerow([
                        input_genes[i],
                        corrected_name])

                if corrected_name not in genes:
                    genes.append(corrected_name)

    if input_genes:
        coverage = " ({:.2f}%)".format((len(genes) / len(input_genes))*100)
    else:
        coverage = ""
    print("Number of unique genes excluding synonyms: {}{}".format(
          len(genes), coverage), end=2*"\n")

    return sorted(genes)


def build_pathway_tables(
        db_fp,
        log_fp,
        pathway_tsv_fp,
        genes):
    """"""
    logging.getLogger("requests").setLevel(logging.CRITICAL)
    logging.basicConfig(
        filename=log_fp,
        filemode="w",
        level=logging.ERROR,
        format="%(levelname)s\t%(asctime)s\t%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")

    wikipathways_release = get_wikipathways_release()
    if wikipathways_release:
        print("WikiPathways release: {}".format(wikipathways_release))
    else:
        print("WikiPathways release unavailable.")

    tsv_file = open(pathway_tsv_fp, "w", buffering=1)
    tsv_writer = csv.writer(tsv_file, delimiter="\t", lineterminator="\n")
    tsv_header = ["Gene_A",
                  "Pathway",
                  "Pathway_Version",
                  "Pathway_Name",
                  "Gene_B",
                  "Upstream",
                  "Downstream",
                  "Upregulating",
                  "Downregulating",
                  "Upregulated",
                  "Downregulated"]
    tsv_writer.writerow(tsv_header)


    pathway_db = sqlite3.connect(db_fp)
    pathway_db_cursor = pathway_db.cursor()

    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS
            Pathways
        """)
    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS
            PathwayNames
        """)
    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS
            UpstreamGenes
        """)
    pathway_db_cursor.execute("""
        DROP TABLE IF EXISTS
            DownstreamGenes
        """)

    pathway_db_cursor.execute("""
        CREATE TABLE
            Pathways (
                gene TEXT,
                pathway TEXT,
                PRIMARY KEY (
                    gene,
                    pathway)
                ON CONFLICT IGNORE)
        """)
    pathway_db_cursor.execute("""
        CREATE TABLE
            PathwayNames (
                pathway TEXT,
                name TEXT,
                PRIMARY KEY (
                    pathway)
                ON CONFLICT IGNORE)
        """)
    pathway_db_cursor.execute("""
        CREATE TABLE
            UpstreamGenes (
                gene TEXT,
                pathway TEXT,
                upstream_gene TEXT,
                upregulating INTEGER,
                downregulating INTEGER,
                PRIMARY KEY (
                    gene,
                    pathway,
                    upstream_gene)
                ON CONFLICT IGNORE)
        """)
    pathway_db_cursor.execute("""
        CREATE TABLE
            DownstreamGenes (
                gene TEXT,
                pathway TEXT,
                downstream_gene TEXT,
                upregulated INTEGER,
                downregulated INTEGER,
                PRIMARY KEY (
                    gene,
                    pathway,
                    downstream_gene)
                ON CONFLICT IGNORE)
        """)

    num_pathways_available = 0
    all_up_down_stream_genes = set()

    if genes:
        print("Mapping pathway information...")
        gene_num = progressbar.ProgressBar(max_value=len(genes))
        gene_num.update(0)

        for i, gene in enumerate(genes, 1):
            unicode_gene = unicode(gene, "utf-8")
            pathways = wikipathways(gene)
            if pathways:
                num_pathways_available += 1
            for pathway in pathways:
                pathway_db_cursor.execute("""
                    INSERT INTO
                        Pathways
                    VALUES
                        (?, ?)
                    """, (unicode_gene, pathway[0]))

                pathway_db_cursor.execute("""
                    INSERT INTO
                        PathwayNames
                    VALUES
                        (?, ?)
                    """, (pathway[0], pathway[2]))

                entries = [
                    (unicode_gene, pathway[0], upstream_gene,
                    int(pathway[3][upstream_gene][0]),
                    int(pathway[3][upstream_gene][1]))
                    for upstream_gene in pathway[3]]

                for entry in entries:
                    pathway_db_cursor.execute("""
                        INSERT INTO
                            UpstreamGenes
                        VALUES
                            (?, ?, ?, ?, ?)
                        """, entry)

                entries = [
                    (unicode_gene, pathway[0], downstream_gene,
                    int(pathway[4][downstream_gene][0]),
                    int(pathway[4][downstream_gene][1]))
                    for downstream_gene in pathway[4]]

                for entry in entries:
                    pathway_db_cursor.execute("""
                        INSERT INTO
                            DownstreamGenes
                        VALUES
                            (?, ?, ?, ?, ?)
                        """, entry)

                up_down_stream_genes = set(pathway[3]) | set(pathway[4])

                tsv_row = [entry.encode("utf-8") for entry in pathway[:3]]

                if not up_down_stream_genes:
                    tsv_rows = [[
                        gene,
                        tsv_row[0],
                        tsv_row[1],
                        tsv_row[2]]]
                    tsv_rows[0].extend(["NA" for j in range(7)])

                else:
                    tsv_rows = [
                        [gene,
                        tsv_row[0],
                        tsv_row[1],
                        tsv_row[2],
                        up_down_stream_gene.encode("utf-8"),
                        up_down_stream_gene in pathway[3],
                        up_down_stream_gene in pathway[4],
                        (up_down_stream_gene in pathway[3] and
                            pathway[3][up_down_stream_gene][0]),
                        (up_down_stream_gene in pathway[3] and
                            pathway[3][up_down_stream_gene][1]),
                        (up_down_stream_gene in pathway[4] and
                            pathway[4][up_down_stream_gene][0]),
                        (up_down_stream_gene in pathway[4] and
                            pathway[4][up_down_stream_gene][1])]
                        for up_down_stream_gene
                            in sorted(up_down_stream_genes)]

                for tsv_row in tsv_rows:
                    tsv_writer.writerow(tsv_row)

                all_up_down_stream_genes |= up_down_stream_genes

            pathway_db.commit()
            gene_num.update(i)

        gene_num.finish()
    logging.shutdown()
    tsv_file.close()

    if genes:
        coverage = " ({:.2f}%)".format(
                (num_pathways_available / len(genes))*100)
    else:
        coverage = ""
    print("Number of genes covered by pathway data: {}{}".format(
        num_pathways_available, coverage))

    pathway_db.close()

    return None

def merge_gtex_hpm(gene_exp, protein_exp):
    """"""
    tissue_map = {
        "Adipose_Subcutaneous": None,
        "Adipose_Visceral_Omentum": None,
        "Adrenal_Gland": "Adrenal",
        "Artery_Aorta": None,
        "Artery_Coronary": None,
        "Artery_Tibial": None,
        "Bladder": "Urinary_Bladder",
        "Brain_Amygdala": None,
        "Brain_Anterior_cingulate_cortex_BA24": None,
        "Brain_Caudate_basal_ganglia": None,
        "Brain_Cerebellar_Hemisphere": None,
        "Brain_Cerebellum": None,
        "Brain_Cortex": None,
        "Brain_Frontal_Cortex_BA9": "Frontal_Cortex",
        "Brain_Hippocampus": None,
        "Brain_Hypothalamus": None,
        "Brain_Nucleus_accumbens_basal_ganglia": None,
        "Brain_Putamen_basal_ganglia": None,
        "Brain_Spinal_cord_cervical_c-1": "Spinal_Cord",
        "Brain_Substantia_nigra": None,
        "Breast_Mammary_Tissue": None,
        "Cells_EBV-transformed_lymphocytes": None,
        "Cells_Transformed_fibroblasts": None,
        "Cervix_Ectocervix": None,
        "Cervix_Endocervix": None,
        "Colon_Sigmoid": "Colon",
        "Colon_Transverse": "Colon",
        "Esophagus_Gastroesophageal_Junction": "Esophagus",
        "Esophagus_Mucosa": "Esophagus",
        "Esophagus_Muscularis": "Esophagus",
        "Fallopian_Tube": None,
        "Heart_Atrial_Appendage": "Heart",
        "Heart_Left_Ventricle": "Heart",
        "Kidney_Cortex": "Kidney",
        "Liver": "Liver",
        "Lung": "Lung",
        "Minor_Salivary_Gland": None,
        "Muscle_Skeletal": None,
        "Nerve_Tibial": None,
        "Ovary": "Ovary",
        "Pancreas": "Pancreas",
        "Pituitary": None,
        "Prostate": "Prostate",
        "Skin_Not_Sun_Exposed_Suprapubic": None,
        "Skin_Sun_Exposed_Lower_leg": None,
        "Small_Intestine_Terminal_Ileum": None,
        "Spleen": None,
        "Stomach": None,
        "Testis": "Testis",
        "Thyroid": None,
        "Uterus": None,
        "Vagina": None,
        "Whole_Blood": None
        }

    """
    Unmapped HPM tissues:
        Gallbladder
        Rectum
        Retina
    """

    exp = {"gene": {},
           "protein": {}}

    for gene in gene_exp:
        for tissue in tissue_map:
            if tissue_map[tissue] is not None:
                if tissue in gene_exp[gene]:
                    if gene not in exp["gene"]:
                        exp["gene"][gene] = {}
                    exp["gene"][gene][tissue] = gene_exp[gene][tissue]

    for gene in protein_exp:
        for tissue in tissue_map:
            if tissue_map[tissue] is not None:
                if tissue_map[tissue] in protein_exp[gene]:
                    if gene not in exp["protein"]:
                        exp["protein"][gene] = {}
                    exp["protein"][gene][tissue] =\
                            protein_exp[gene][tissue_map[tissue]]

    return exp

def produce_pathway_summary(
        buffer_size,
        db_fp,
        summary_fp,
        hpm,
        input_genes):

    pathway_db = sqlite3.connect(db_fp)
    pathway_db.text_factory = str
    pathway_db_cursor = pathway_db.cursor()

    up_down_stream_genes = set()
    pathways = {}
    genes = set(input_genes)
    for gene in input_genes:
        pathways[gene] = {}
        for (pathway, upstream_gene, upregulating,
             downregulating) in pathway_db_cursor.execute("""
                    SELECT
                        pathway,
                        upstream_gene,
                        upregulating,
                        downregulating
                    FROM
                        UpstreamGenes
                    WHERE
                        gene = ?
                    """, (gene,)):
            if pathway not in pathways[gene]:
                pathways[gene][pathway] = {
                    "upstream": {},
                    "downstream": {}}
            pathways[gene][pathway]["upstream"][upstream_gene] = [
                    bool(upregulating),
                    bool(downregulating)]
            genes.add(upstream_gene)

        for (pathway, downstream_gene, upregulated,
             downregulated) in pathway_db_cursor.execute("""
                    SELECT
                        pathway,
                        downstream_gene,
                        upregulated,
                        downregulated
                    FROM
                        DownstreamGenes
                    WHERE
                        gene = ?
                    """, (gene,)):
            if pathway not in pathways[gene]:
                pathways[gene][pathway] = {
                    "upstream": {},
                    "downstream": {}}
            pathways[gene][pathway]["downstream"][downstream_gene] = [
                    bool(upregulated),
                    bool(downregulated)]
            genes.add(downstream_gene)

    tissues = set()
    exp = {}
    for gene in genes:
        exp[gene] = {}
        for (tissue, gene_exp, gene_exp_z, gene_exp_p
            ) in pathway_db_cursor.execute("""
                SELECT
                    tissue,
                    gene_exp,
                    gene_exp_z,
                    gene_exp_p
                FROM
                    GeneExpression
                WHERE
                   gene = ?
                """, (gene,)):
            tissues.add(tissue)
            exp[gene][tissue] = [(gene_exp, gene_exp_z, gene_exp_p)]

    if hpm:
        for gene in genes:
            for (tissue, protein_exp, protein_exp_z, protein_exp_p
                ) in pathway_db_cursor.execute("""
                    SELECT
                        tissue,
                        protein_exp,
                        protein_exp_z,
                        protein_exp_p
                    FROM
                        ProteinExpression
                    WHERE
                        gene = ?
                    """, (gene,)):
                exp[gene][tissue].append(
                        (protein_exp, protein_exp_z, protein_exp_p))


    pathway_db.close()

    tissues = sorted(tissues)

    summary_file = open(summary_fp, "w", buffering=buffer_size)
    summary_writer = csv.writer(summary_file, delimiter="\t",
                                lineterminator="\n")
    summary_header = [
        "Gene_A",
        "Pathway",
        "Gene_B",
        "Upstream",
        "Downstream",
        "Upregulating",
        "Downregulating",
        "Upregulated",
        "Downregulated",
        "Tissue",
        "Gene_A_Gene_Expression",
        "Gene_A_Gene_Expression_z-Score",
        "Gene_A_Gene_Expression_p-Value",
        "Gene_B_Gene_Expression",
        "Gene_B_Gene_Expression_z-Score",
        "Gene_B_Gene_Expression_p-Value"]

    if hpm:
        summary_header.insert(13, "Gene_A_Protein_Expression")
        summary_header.insert(14, "Gene_A_Protein_Expression_z-Score")
        summary_header.insert(15, "Gene_A_Protein_Expression_p-Value")
        summary_header.extend([
            "Gene_B_Protein_Expression",
            "Gene_B_Protein_Expression_z-Score",
            "Gene_B_Protein_Expression_p-Value"])

    summary_writer.writerow(summary_header)

    num_rows = 0
    for gene in input_genes:
        if not pathways[gene]:
            for tissue in tissues:
                num_rows += 1
        for pathway in pathways[gene]:
            up_down_stream_genes = set(itertools.chain.from_iterable(
                pathways[gene][pathway].values()))
            for up_down_stream_gene in up_down_stream_genes:
                    for tissue in tissues:
                            num_rows += 1
    if num_rows:
        print("Writing summary file...")
        row_num = progressbar.ProgressBar(max_value=num_rows)
        row_num.update(0)

        current_row = 0
        for gene in input_genes:
            if not pathways[gene]:
                for tissue in tissues:
                    to_file = [
                        gene]
                    to_file.extend(["NA" for i in range(12)])
                    to_file.append(tissue)

                    to_file.extend([
                        exp[gene][tissue][0][i]
                        if exp[gene] and exp[gene][tissue][0][i] else "NA"
                            for i in range(3)])

                    if hpm:
                        to_file.extend([
                            exp[gene][tissue][1][i]
                            if exp[gene] and exp[gene][tissue][1][i] else "NA"
                                for i in range(3)])

                    to_file.extend(["NA" for i in range(6)])

                    summary_writer.writerow(to_file)

                    current_row += 1
                    row_num.update(current_row)

            for pathway in sorted(pathways[gene],
                    key=lambda pathway: int(
                        "".join([char for char in pathway
                            if char.isdigit()]))):
                up_down_stream_genes = set(itertools.chain.from_iterable(
                    pathways[gene][pathway].values()))

                for up_down_stream_gene in sorted(up_down_stream_genes):
                    for tissue in tissues:
                        to_file = [
                            gene,
                            pathway,
                            up_down_stream_gene,
                            up_down_stream_gene in pathways[gene][pathway][
                                "upstream"],
                            up_down_stream_gene in pathways[gene][pathway][
                                "downstream"],
                            up_down_stream_gene in pathways[gene][pathway][
                                "upstream"] and
                            pathways[gene][pathway]["upstream"][
                                up_down_stream_gene][0],
                            up_down_stream_gene in pathways[gene][pathway][
                                "upstream"] and
                            pathways[gene][pathway]["upstream"][
                                up_down_stream_gene][1],
                            up_down_stream_gene in pathways[gene][pathway][
                                "downstream"] and
                            pathways[gene][pathway]["downstream"][
                                up_down_stream_gene][0],
                            up_down_stream_gene in pathways[gene][pathway][
                                "downstream"] and
                            pathways[gene][pathway]["downstream"][
                                up_down_stream_gene][1],
                            tissue]

                        to_file.extend([
                            exp[gene][tissue][0][i]
                            if exp[gene] and exp[gene][tissue][0][i] else "NA"
                                for i in range(3)])

                        if hpm:
                            to_file.extend([
                                exp[gene][tissue][1][i]
                                if exp[gene] and exp[gene][tissue][1][i] else "NA"
                                    for i in range(3)])

                        to_file.extend([
                            exp[up_down_stream_gene][tissue][0][i]
                            if (exp[up_down_stream_gene] and
                            exp[up_down_stream_gene][tissue][0][i]) else "NA"
                                for i in range(3)])

                        if hpm:
                            to_file.extend([
                                exp[up_down_stream_gene][tissue][1][i]
                                if (exp[up_down_stream_gene] and
                                exp[up_down_stream_gene][tissue][1][i]) else "NA"
                                    for i in range(3)])

                        summary_writer.writerow(to_file)

                        current_row += 1
                        row_num.update(current_row)

        row_num.finish()
        print()
    summary_file.close()

    return None

def cluster_genes(tissues, genes, expr, position):
    if len(genes) < 3:
        return genes

    no_expression = []
    no_data = []
    j = 0
    for i in range(len(genes)):
        if not any(expr[genes[i-j]].values()):
            if None in expr[genes[i-j]].values():
                no_data.append(genes.pop(i-j))
            else:
                no_expression.append(genes.pop(i-j))
            j += 1

    if position:
        no_data.extend(no_expression)
        exclude = no_data
    else:
        no_expression.extend(no_data)
        exclude = no_expression

    if len(genes) > 2:
        distr = {gene: np.array([expr[gene][tissue] for tissue in tissues])
                for gene in genes}
        dist_mat = np.nan_to_num(np.array([[abs(st.spearmanr(
                                distr[genes[i]], distr[genes[j]])[0])
                        if j < i else np.nan
                        for j in range(len(genes))]
                        for i in range(len(genes))]))
        link_mat = cl.hierarchy.linkage(dist_mat, "complete",
                                        optimal_ordering=True)

        genes = [genes[i] for i in cl.hierarchy.leaves_list(link_mat)]

    if position:
        exclude.extend(genes)
        return exclude
    else:
        genes.extend(exclude)
        return genes

def plot_heatmaps(
        db_fp,
        input_genes,
        plot_dir_z,
        plot_dir_p,
        hpm,
        p_value,
        numbers):
    """"""
    font = {"size": 8}

    matplotlib.rc("font", **font)

    pathway_db = sqlite3.connect(db_fp)
    pathway_db.text_factory = str
    pathway_db_cursor = pathway_db.cursor()

    genes = {}
    all_genes = set(input_genes)
    for gene in input_genes:
        genes[gene] = {}
        for (upstream_gene, upregulating,
             downregulating) in pathway_db_cursor.execute("""
                    SELECT
                        upstream_gene,
                        upregulating,
                        downregulating
                    FROM
                        UpstreamGenes
                    WHERE
                        gene = ?
                    """, (gene,)):
            if upstream_gene not in genes[gene]:
                genes[gene][upstream_gene] = [
                        True,
                        False,
                        bool(upregulating),
                        bool(downregulating),
                        False,
                        False]
            else:
                genes[gene][upstream_gene][0] = True
                genes[gene][upstream_gene][2] = bool(upregulating)
                genes[gene][upstream_gene][3] = bool(downregulating)
            all_genes.add(upstream_gene)

        for (downstream_gene, upregulated,
             downregulated) in pathway_db_cursor.execute("""
                    SELECT
                        downstream_gene,
                        upregulated,
                        downregulated
                    FROM
                        DownstreamGenes
                    WHERE
                        gene = ?
                    """, (gene,)):
            if downstream_gene not in genes[gene]:
                genes[gene][downstream_gene] = [
                        False,
                        True,
                        False,
                        False,
                        bool(upregulated),
                        bool(downregulated)]
            else:
                genes[gene][downstream_gene][1] = True
                genes[gene][downstream_gene][4] = bool(upregulated)
                genes[gene][downstream_gene][5] = bool(downregulated)
            all_genes.add(downstream_gene)

    expr = {}
    for gene in all_genes:
        expr[gene] = {}
        for (tissue, gene_exp_z, gene_exp_p)in pathway_db_cursor.execute("""
                SELECT
                    tissue,
                    gene_exp_z,
                    gene_exp_p
                FROM
                    GeneExpression
                WHERE
                    gene = ?
                """, (gene,)):

            expr[gene][tissue] = [(gene_exp_z, gene_exp_p)]

    if hpm:
        for gene in all_genes:
            for (tissue, protein_exp_z, protein_exp_p
                ) in pathway_db_cursor.execute("""
                    SELECT
                        tissue,
                        protein_exp_z,
                        protein_exp_p
                    FROM
                        ProteinExpression
                    WHERE
                        gene = ?
                    """, (gene,)):

                expr[gene][tissue].append((protein_exp_z, protein_exp_p))

    tissues = sorted([tissue for (tissue,) in pathway_db_cursor.execute("""
                SELECT DISTINCT
                    tissue
                FROM
                    GeneExpression
                """)])

    pathway_db.close()

    if genes:
        print("Plotting gene expression heatmaps...")
        gene_num = progressbar.ProgressBar(max_value=len(genes))
        gene_num.update(0)

        for i, gene_a in enumerate(genes, 1):
            interacting_genes = sorted([gene_b for gene_b in genes[gene_a]])

            upstream_genes = [
                gene_b for gene_b in interacting_genes
                if genes[gene_a][gene_b][0]]


            if upstream_genes:
                upstream_expr_gene = {gene: {tissue: expr[gene][tissue][0][0]
                                            if expr[gene] else None
                                            for tissue in tissues}
                                        for gene in upstream_genes}
                if hpm:
                    upstream_expr_protein = {gene: {tissue:
                                                expr[gene][tissue][1][0]
                                                if expr[gene] else None
                                                for tissue in tissues}
                                            for gene in upstream_genes}

                upstream_genes_gene = cluster_genes(tissues, upstream_genes,
                                            upstream_expr_gene,
                                            True)
                if hpm:
                    upstream_genes_protein = cluster_genes(tissues, upstream_genes,
                                                upstream_expr_protein,
                                                True)

                upstream_reg_type_data = np.array(
                    [[float(genes[gene_a][gene_b][2])
                    if genes[gene_a][gene_b][2] != genes[gene_a][gene_b][3]
                    else 0.5]
                    for gene_b in upstream_genes_gene])

            downstream_genes = [
                gene_b for gene_b in interacting_genes
                if genes[gene_a][gene_b][1]]

            if downstream_genes:
                downstream_expr_gene = {gene:
                                            {tissue: expr[gene][tissue][0][0]
                                            if expr[gene] else None
                                            for tissue in tissues}
                                        for gene in downstream_genes}
                if hpm:
                    downstream_expr_protein = {gene:
                                                {tissue: expr[gene][tissue][0][1]
                                                if expr[gene] else None
                                                for tissue in tissues}
                                            for gene in downstream_genes}

                downstream_genes_gene = cluster_genes(tissues,
                                                      downstream_genes,
                                                      downstream_expr_gene,
                                                      False)

                if hpm:
                    downstream_genes_protein = cluster_genes(tissues,
                                                            downstream_genes,
                                                        downstream_expr_protein,
                                                        False)

                downstream_reg_type_data = np.array(
                    [[float(genes[gene_a][gene_b][4])
                    if genes[gene_a][gene_b][4] != genes[gene_a][gene_b][5]
                    else 0.5]
                    for gene_b in downstream_genes])


            for (j, rmin, rmax, cmap, plot_dir, label) in (
                (0, st.norm.ppf(0.5*p_value), -st.norm.ppf(0.5*p_value),
                 "RdBu_r", plot_dir_z, "z-Score"),
                (1, p_value, 0.5, "Reds_r", plot_dir_p, "p-Value")):

                gene_reg_type_data = np.array([[0.5]])

                gene_gene_data = np.array(
                    [[expr[gene_a][tissue][0][j]
                      if (expr[gene_a] and
                          expr[gene_a][tissue][0][j] is not None)
                      else np.nan
                    for tissue in tissues]])

                if hpm:
                    gene_protein_data = np.array(
                        [[expr[gene_a][tissue][1][j]
                        if (expr[gene_a] and
                            expr[gene_a][tissue][1][j] is not None)
                        else np.nan
                        for tissue in tissues]])

                scale = 0.5
                if hpm:
                    fig = plt.figure(figsize=(scale*(len(tissues)+2),
                        scale*(2+2*(len(upstream_genes)+len(downstream_genes)))))
                else:
                    fig = plt.figure(figsize=(scale*(len(tissues)+2),
                        scale*(1+len(upstream_genes)+len(downstream_genes))))

                if hpm:
                    gs = gridspec.GridSpec(3, 1,
                                        height_ratios=[
                                            (1+len(upstream_genes)+
                                                len(downstream_genes)),
                                            (1+len(upstream_genes)+
                                                len(downstream_genes)),
                                            1],
                       hspace=1/(1+len(upstream_genes) +
                                len(downstream_genes)))

                else:
                    gs = gridspec.GridSpec(2, 1,
                                        height_ratios=[
                                            (1+len(upstream_genes)+
                                                len(downstream_genes)),
                                            1],
                       hspace=1/(1+len(upstream_genes) +
                                len(downstream_genes)))


                gs1 = gridspec.GridSpecFromSubplotSpec(3, 2,
                        subplot_spec=gs[0],
                        width_ratios=[1, len(tissues)],
                        height_ratios=[len(upstream_genes), 1,
                            len(downstream_genes)],
                       wspace= 0.05,
                       hspace=1/(2+2*(len(upstream_genes) +
                                len(downstream_genes))))

                if hpm:
                    gs2 = gridspec.GridSpecFromSubplotSpec(3, 2,
                            subplot_spec=gs[1],
                            width_ratios=[1, len(tissues)],
                            height_ratios=[len(upstream_genes), 1,
                                len(downstream_genes)],
                        wspace= 0.05,
                        hspace=1/(2+2*(len(upstream_genes) +
                                    len(downstream_genes))))
                if hpm:
                    gs3 = gridspec.GridSpecFromSubplotSpec(1, 2,
                            subplot_spec=gs[2],
                            width_ratios=[1+len(tissues)-10, 10],
                        wspace= 0.05,
                        hspace=1/(2+2*(len(upstream_genes) +
                                    len(downstream_genes))))
                else:
                    gs3 = gridspec.GridSpecFromSubplotSpec(1, 2,
                            subplot_spec=gs[1],
                            width_ratios=[1+len(tissues)-10, 10],
                        wspace= 0.05,
                        hspace=1/(2+2*(len(upstream_genes) +
                                    len(downstream_genes))))

                cbar_ax = plt.subplot(gs3[1])
                leg_ax = plt.subplot(gs3[0])

                custom_patches = [
                 Patch(facecolor="black", edgecolor="black",
                       label="Upregulatory Protein Interaction"),
                 Patch(facecolor="grey", edgecolor="black",
                       label="Not Applicable or Contradictory Data"),
                 Patch(facecolor="white", edgecolor="black",
                       label="Downregulatory Protein Interaction")]

                leg_ax.legend(handles=custom_patches, loc="lower left",
                              frameon=False, mode="expand")
                leg_ax.axis("off")

                gene_gene_reg = plt.subplot(gs1[2])
                gene_gene = plt.subplot(gs1[3])

                if hpm:
                    gene_protein_reg = plt.subplot(gs2[2])
                    gene_protein = plt.subplot(gs2[3])

                expr_plots = [gene_gene]

                if hpm:
                    expr_plots.append(gene_protein)

                reg_plots = [gene_gene_reg]

                if hpm:
                    reg_plots.append(gene_protein_reg)

                plot_data = {
                        gene_gene_reg: gene_reg_type_data,
                        gene_gene: gene_gene_data}

                if hpm:
                    plot_data[gene_protein_reg] = gene_reg_type_data
                    plot_data[gene_protein] = gene_protein_data

                plot_labels = {
                        gene_gene_reg: [gene_a],
                        gene_gene: []}

                if hpm:
                    plot_labels[gene_protein_reg] = [gene_a]
                    plot_labels[gene_protein] = []

                if upstream_genes:
                    upstream_gene_gene_data = np.array(
                        [[expr[gene_b][tissue][0][j]
                          if (expr[gene_b] and
                              expr[gene_b][tissue][0][j] is not None)
                          else np.nan
                        for tissue in tissues]
                        for gene_b in upstream_genes_gene])

                    if hpm:
                        upstream_gene_protein_data = np.array(
                            [[expr[gene_b][tissue][1][j]
                            if (expr[gene_b] and
                                expr[gene_b][tissue][1][j] is not None)
                            else np.nan
                            for tissue in tissues]
                            for gene_b in upstream_genes_protein])

                    upstream_gene_gene_reg = plt.subplot(gs1[0])
                    upstream_gene_gene = plt.subplot(gs1[1])

                    if hpm:
                        upstream_gene_protein_reg = plt.subplot(gs2[0])
                        upstream_gene_protein = plt.subplot(gs2[1])

                    expr_plots.append(upstream_gene_gene)
                    if hpm:
                        expr_plots.append(upstream_gene_protein)

                    reg_plots.append(upstream_gene_gene_reg)
                    if hpm:
                        reg_plots.append(upstream_gene_protein_reg)

                    plot_data[upstream_gene_gene] = upstream_gene_gene_data
                    plot_data[upstream_gene_gene_reg] = upstream_reg_type_data
                    if hpm:
                        plot_data[upstream_gene_protein] =\
                                upstream_gene_protein_data
                        plot_data[upstream_gene_protein_reg] =\
                                upstream_reg_type_data

                    plot_labels[upstream_gene_gene_reg] = upstream_genes_gene
                    if hpm:
                        plot_labels[upstream_gene_protein_reg] =\
                                upstream_genes_protein

                if downstream_genes:
                    downstream_gene_gene_data = np.array(
                        [[expr[gene_b][tissue][0][j]
                          if (expr[gene_b] and
                              expr[gene_b][tissue][0][j] is not None)
                          else np.nan
                        for tissue in tissues]
                        for gene_b in downstream_genes_gene])

                    if hpm:
                        downstream_gene_protein_data = np.array(
                            [[expr[gene_b][tissue][1][j]
                            if (expr[gene_b] and
                                expr[gene_b][tissue][1][j] is not None)
                            else np.nan
                            for tissue in tissues]
                            for gene_b in downstream_genes_protein])

                    downstream_gene_gene_reg = plt.subplot(gs1[4])
                    downstream_gene_gene = plt.subplot(gs1[5])
                    if hpm:
                        downstream_gene_protein_reg = plt.subplot(gs2[4])
                        downstream_gene_protein = plt.subplot(gs2[5])

                    expr_plots.append(downstream_gene_gene)
                    if hpm:
                        expr_plots.append(downstream_gene_protein)

                    reg_plots.append(downstream_gene_gene_reg)
                    if hpm:
                        reg_plots.append(downstream_gene_protein_reg)

                    plot_data[downstream_gene_gene] =\
                            downstream_gene_gene_data
                    plot_data[downstream_gene_gene_reg] =\
                            downstream_reg_type_data
                    if hpm:
                        plot_data[downstream_gene_protein] =\
                                downstream_gene_protein_data
                        plot_data[downstream_gene_protein_reg] =\
                                downstream_reg_type_data

                    plot_labels[downstream_gene_gene_reg] =\
                            downstream_genes_gene
                    if hpm:
                        plot_labels[downstream_gene_protein_reg] =\
                                downstream_genes_protein

                for plot in expr_plots:
                    plot_data[plot] = np.ma.masked_invalid(plot_data[plot])
                    im = plot.imshow(plot_data[plot], aspect="auto",
                                     vmin=rmin, vmax=rmax, cmap=cmap)

                    plot.set_xticks([])
                    plot.set_xticklabels([])
                    plot.set_yticks([])
                    plot.set_yticklabels([])

                    if numbers:
                        for j in range(np.ma.size(plot_data[plot], 0)):
                            for k in range(len(tissues)):
                                text = plot.text(k, j,
                                    "{:.2f}".format(plot_data[plot][j, k]),
                                    ha="center", va="center", color="black")

                cbar = plt.colorbar(im, ax=expr_plots, cax=cbar_ax,
                        orientation="horizontal")
                cbar.set_label(
                        "Standardized Expression Rate ({})".format(label))

                for plot in reg_plots:
                    im = plot.imshow(plot_data[plot], aspect="auto", vmin=0.0,
                                    vmax=1.0, cmap="binary")

                    plot.set_xticks([])
                    plot.set_xticklabels([])
                    plot.set_yticks(np.arange(len(plot_labels[plot])))
                    plot.tick_params(length=0.0)
                    plot.set_yticklabels(plot_labels[plot])

                if upstream_genes:
                    upstream_gene_gene.set_xticks(np.arange(len(tissues)))
                    upstream_gene_gene.tick_params(length=0, top=True,
                            bottom=False, labeltop=True, labelbottom=False)
                    upstream_gene_gene.set_xticklabels([
                        tissue.replace("_", " ") for tissue in tissues])
                    plt.setp(upstream_gene_gene.get_xticklabels(),
                            rotation=-90, ha="right", rotation_mode="anchor")

                else:
                    gene_gene.set_xticks(np.arange(len(tissues)))
                    gene_gene.tick_params(length=0)
                    gene_gene.set_xticklabels([
                        tissue.replace("_", " ") for tissue in tissues])
                    gene_gene.tick_params(length=0, top=True,
                            bottom=False, labeltop=True, labelbottom=False)
                    plt.setp(gene_gene.get_xticklabels(),
                            rotation=-90, ha="right", rotation_mode="anchor")

                plt.savefig(os.path.join(plot_dir, "{}.png".format(gene_a)),
                            format="png", dpi=300, bbox_inches="tight")

                plt.close()
            gene_num.update(i)
        gene_num.finish()
    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--inputs", nargs='+', required=True,
                        help="The the dbSNP IDs or loci of SNPs of interest in the format " +
                        "\"chr<x>:<locus>\"")
    parser.add_argument("-c", "--config",
                        default=os.path.join(os.path.dirname(__file__),
                                             "../docs/codes3d.conf"),
                        help="The configuration file to be used for this " +
                        "hiCquery run (default: conf.py)")
    parser.add_argument("-n", "--include_cell_lines", nargs='+',
                        help="Space-separated list of cell lines to include " +
                        "(others will be ignored). NOTE: Mutually exclusive " +
                        "with EXCLUDE_CELL_LINES.")
    parser.add_argument("-x", "--exclude_cell_lines", nargs='+',
                        help="Space-separated list of cell lines to exclude " +
                        "(others will be included). NOTE: Mutually exclusive " +
                        "with INCLUDE_CELL_LINES.")
    parser.add_argument("-o", "--output_dir", default="codes3d_output",
                        help="The directory in which to output results " +
                        "(\"hiCquery_output\" by default).")
    parser.add_argument("-q", "--query_databases", type=str,
                        default="both",
                        help="[local] for only local databases. " +
                        "This will only include cis-eQTLs if using downloadable " +
                        "GTEx dataset. " +
                        "[online] for only online GTEx queries. " +
                        "[both] for both local and online queries " +
                        "(default: both).")
    parser.add_argument("-d", "--do_not_produce_summary", action="store_true", default=False,
                        help="Do not produce summary files, stop process after " +
                        "querying GTEx (default: False).")
    parser.add_argument("-s", "--suppress_intermediate_files", action="store_true",
                        default=False,
                        help="Do not produce intermediate " +
                        "files. These can be used to run the pipeline from " +
                        "an intermediate stage in the event of interruption " +
                        "(default: False).")
    parser.add_argument("-p", "--num_processes", type=int, default=1,
                        help="Desired number of processes for multithreading " +
                        "(default: 1).")
    parser.add_argument("-f", "--fdr_threshold", type=float, default=0.05,
                        help="The FDR threshold to consider an eQTL " +
                        "statistically significant (default: 0.05).")
    parser.add_argument("-r", "--restriction_enzymes", nargs='+',
                        help="Space-separated list of  " +
                        "restriction enzymes used in HiC data.")
    parser.add_argument("-b","--buffer_size",type=int,default=1048576,
                        help="Buffer size applied to file I/O during " +\
                        "compilation (default: 1048576).")
    parser.add_argument("-k", "--num_processes_summary", type=int,
                        default=(psutil.cpu_count() // 2),
                        help="The number of processes for compilation of " +\
    help="The number of processes for compilation of " +

                        "the results (default: %s)." %
                        str((psutil.cpu_count() // 2)))
    parser.add_argument("--pathway_config",
                        default=os.path.join(
                        os.path.dirname(__file__),
                        "../docs/codes3d_pathway.conf"),
                        help="Configuration file to be used.")
    parser.add_argument("--hpm", action="store_true", default=False,
                        help="Include protein expression data")
    parser.add_argument("--numbers", action="store_true", default=False,
                        help="Label the heatmap entries with the " +\
                        "represented value.")
    parser.add_argument("--p_value", type=float, default=0.05,
                        help="p-value defining the range of the heatmap " +\
                        "gradient")
    parser.add_argument("--summary_files", nargs="+",
                        help="List of summary files")

    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    
    lib_dir = os.path.join(os.path.dirname(__file__),config.get("Defaults",
                                                                "LIB_DIR"))
    lib_dir = os.path.join(os.path.dirname(__file__),
                           config.get("Defaults", "LIB_DIR"))
    snp_database_fp = os.path.join(os.path.dirname(__file__),
                                   config.get("Defaults", "SNP_DATABASE_FP"))
    hic_data_dir = os.path.join(os.path.dirname(__file__),
                                config.get("Defaults", "HIC_DATA_DIR"))
    fragment_bed_fp = os.path.join(os.path.dirname(__file__),
                                   config.get("Defaults", "FRAGMENT_BED_FP"))
    fragment_database_fp = os.path.join(os.path.dirname(__file__),
                                        config.get("Defaults",
                                                   "FRAGMENT_DATABASE_FP"))
    gene_bed_fp = os.path.join(os.path.dirname(__file__),
                               config.get("Defaults", "GENE_BED_FP"))
    gene_database_fp = os.path.join(os.path.dirname(__file__),
                                    config.get("Defaults", "GENE_DATABASE_FP"))
    eqtl_data_dir = os.path.join(os.path.dirname(__file__),
                                 config.get("Defaults", "EQTL_DATA_DIR"))
    expression_table_fp = os.path.join(os.path.dirname(__file__),
                                       config.get("Defaults",
                                                  "EXPRESSION_TABLE_FP"))
    gene_dict_fp = os.path.join(os.path.dirname(__file__),
                                config.get("Defaults", "GENE_DICT_FP"))
    snp_dict_fp = os.path.join(os.path.dirname(__file__),
                               config.get("Defaults", "SNP_DICT_FP"))
    GTEX_CERT = os.path.join(os.path.dirname(__file__),
                             config.get("Defaults", "GTEX_CERT"))
    HIC_RESTRICTION_ENZYMES = [e.strip() for e in \
                                config.get("Defaults",
                                    "HIC_RESTRICTION_ENZYMES").split(',')]
    restriction_enzymes, include_cell_lines, exclude_cell_lines =\
        parse_parameters(args.restriction_enzymes,
                         args.include_cell_lines,
                         args.exclude_cell_lines)

    pathway_config = configparser.ConfigParser()
    pathway_config.read(args.pathway_config)

    gene_synonym_map_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("LIB", "gene_synonym_map_fp"))

    gtex_gene_exp_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("LIB", "gtex_gene_exp_fp"))

    hpm_gene_exp_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("LIB", "hpm_gene_exp_fp"))

    hpm_map_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("LIB", "hpm_map_fp"))

    hpm_protein_exp_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("LIB", "hpm_protein_exp_fp"))

    input_gene_map_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("OUT", "input_gene_map_fp"))

    log_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("OUT", "log_fp"))

    db_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("LIB", "db_fp"))

    pathway_tsv_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("OUT", "pathway_tsv_fp"))

    gene_exp_tsv_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("OUT", "gene_exp_tsv_fp"))

    protein_exp_tsv_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("OUT", "protein_exp_tsv_fp"))

    summary_fp = os.path.join(os.path.dirname(__file__),
        pathway_config.get("OUT", "summary_fp"))

    plot_dir = os.path.join(os.path.dirname(__file__),
        pathway_config.get("OUT", "plot_dir"))

    plot_dir_z = os.path.join(os.path.dirname(__file__),
        pathway_config.get("OUT", "plot_dir_z"))

    plot_dir_p = os.path.join(os.path.dirname(__file__),
        pathway_config.get("OUT", "plot_dir_p"))

    out_dir = os.path.join(os.path.dirname(__file__),
        pathway_config.get("OUT", "out_dir"))

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    heatmap_dir = os.path.join(args.output_dir,
            os.path.basename(os.path.normpath(plot_dir)))

    if not os.path.isdir(heatmap_dir):
        os.mkdir(heatmap_dir)

    plot_dir_z = os.path.join(heatmap_dir,
                            os.path.basename(os.path.normpath(plot_dir_z)))
    plot_dir_p = os.path.join(heatmap_dir,
                            os.path.basename(os.path.normpath(plot_dir_p)))

    for sub_dir in (plot_dir_z, plot_dir_p):
        if not os.path.isdir(sub_dir):
            os.mkdir(sub_dir)
        else:
            for heatmap in os.listdir(sub_dir):
                os.remove(os.path.join(sub_dir, heatmap))

    snps = process_inputs(
        args.inputs, snp_database_fp, lib_dir,
        restriction_enzymes, args.output_dir,
        args.suppress_intermediate_files)
    interactions = find_interactions(
        snps, lib_dir, include_cell_lines,
        exclude_cell_lines, args.output_dir, args.num_processes_summary,
        args.suppress_intermediate_files)
    genes = find_genes(interactions, lib_dir, gene_bed_fp, gene_dict_fp,
                       args.output_dir, args.suppress_intermediate_files)
    num_sig, p_values = find_eqtls(
        snps, genes, eqtl_data_dir, gene_database_fp,
        args.fdr_threshold, args.query_databases,
        args.num_processes_summary, args.output_dir, gene_dict_fp, snp_dict_fp,
        suppress_intermediate_files=args.suppress_intermediate_files)
    gene_ids = produce_summary(
         p_values, snps, genes, gene_database_fp, expression_table_fp,
         args.fdr_threshold, args.do_not_produce_summary, args.output_dir,
         args.buffer_size, args.num_processes_summary)

    gene_name = GeneName(gene_synonym_map_fp)
    build_pathway_tables(db_fp, log_fp, pathway_tsv_fp,
                         gene_ids)
    build_expression_tables(args.num_processes_summary,
            gtex_gene_exp_fp, hpm_map_fp, hpm_gene_exp_fp,
            hpm_protein_exp_fp, args.hpm, db_fp, gene_exp_tsv_fp,
            protein_exp_tsv_fp, gene_ids)
    produce_pathway_summary(args.buffer_size, db_fp, summary_fp,
                            args.hpm, gene_ids)
    plot_heatmaps(db_fp, gene_ids, plot_dir_z, plot_dir_p, args.hpm,
                  args.p_value, args.numbers)
