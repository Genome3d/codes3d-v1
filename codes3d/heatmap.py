#!/usr/bin/env python

import argparse
import codes3d
import configparser
import os
import psutil

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-b", "--buffer_size", type=int, default=1048576,
        help="Buffer size for file I/O (default: 1MB (1048576)).")
    parser.add_argument("-c", "--config",
        default=os.path.join(
            os.path.dirname(__file__), "../docs/codes3d_pathway.conf"),
        help="Configuration file to be used.")
    parser.add_argument("-f", "--files", nargs="+",
        help="List of files holding input genes")
    parser.add_argument("-g", "--genes", nargs="+",
        help="List of input genes")
    parser.add_argument("-k", "--num_processes", type=int,
        default=(psutil.cpu_count() // 2),
        help="The maximum number of processes (default: {}).".format(
        str((psutil.cpu_count() // 2))))

    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config)

    num_processes = args.num_processes

    gene_synonym_map_fp = os.path.join(os.path.dirname(__file__),
        config.get("LIB", "gene_synonym_map_fp"))

    gtex_gene_exp_fp = os.path.join(os.path.dirname(__file__),
        config.get("LIB", "gtex_gene_exp_fp"))

    hpm_gene_exp_fp = os.path.join(os.path.dirname(__file__),
        config.get("LIB", "hpm_gene_exp_fp"))

    hpm_map_fp = os.path.join(os.path.dirname(__file__),
        config.get("LIB", "hpm_map_fp"))

    hpm_protein_exp_fp = os.path.join(os.path.dirname(__file__),
        config.get("LIB", "hpm_protein_exp_fp"))

    input_gene_map_fp = os.path.join(os.path.dirname(__file__),
        config.get("OUT", "input_gene_map_fp"))

    log_fp = os.path.join(os.path.dirname(__file__),
        config.get("OUT", "log_fp"))

    db_fp = os.path.join(os.path.dirname(__file__),
        config.get("OUT", "db_fp"))

    pathway_tsv_fp = os.path.join(os.path.dirname(__file__),
        config.get("OUT", "pathway_tsv_fp"))

    gene_exp_tsv_fp = os.path.join(os.path.dirname(__file__),
        config.get("OUT", "gene_exp_tsv_fp"))

    protein_exp_tsv_fp = os.path.join(os.path.dirname(__file__),
        config.get("OUT", "protein_exp_tsv_fp"))

    summary_fp = os.path.join(os.path.dirname(__file__),
        config.get("OUT", "summary_fp"))

    plot_dir = os.path.join(os.path.dirname(__file__),
        config.get("OUT", "plot_dir"))

    out_dir = os.path.join(os.path.dirname(__file__),
        config.get("OUT", "out_dir"))

    for directory in (out_dir, plot_dir):
        if not os.path.isdir(directory):
            os.mkdir(directory)

    input_genes = codes3d.parse_input_genes(args.buffer_size,
                gene_synonym_map_fp, input_gene_map_fp, args.genes, args.files)
    genes = codes3d.build_pathway_tables(db_fp, log_fp, pathway_tsv_fp,
            input_genes)
    codes3d.build_expression_tables(args.num_processes, gtex_gene_exp_fp,
        hpm_map_fp, hpm_gene_exp_fp, hpm_protein_exp_fp, db_fp,
        gene_exp_tsv_fp, protein_exp_tsv_fp, genes)
    codes3d.produce_pathway_summary(args.buffer_size, db_fp, summary_fp,
        input_genes)
    codes3d.plot_heatmaps(db_fp, input_genes, plot_dir)



