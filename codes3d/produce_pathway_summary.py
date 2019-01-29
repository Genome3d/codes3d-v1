#!/usr/bin/env python

import argparse
import codes3d
import configparser
import os
import sys
import psutil

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
            "-s","--summary_file",required=True,
            help="summary file produced by produce_summary()")
    parser.add_argument(
            "-c","--config",
        default=os.path.join(os.path.dirname(__file__),
                             "../docs/codes3d.conf"),
            help="The configuration file to be used in this "+\
            "instance (default: conf.py)")
    parser.add_argument(
            "-b","--buffer_size",type=int,default=1048576,
            help="Buffer size applied to file input during compilation "+\
            " (default: 1048576).")
    parser.add_argument("-e", "--significant_expression", type=float,
            default=0.05,
            help="p-value of significant expression variation "+\
            "(default: 0.05).")

    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    pathway_db_fp = os.path.join(os.path.dirname(__file__),
                                 config.get("Defaults", "PATHWAY_DB_FP"))
    pathway_db_gene_name_fp = os.path.join(os.path.dirname(__file__),
                                           config.get("Defaults",
                                               "PATHWAY_DB_GENE_NAME_FP"))
    pathway_sum_fp = os.path.join(os.path.dirname(__file__),
                                  config.get("Defaults",
                                             "PATHWAY_SUM_FP"))
    pathway_sum_gene_map_fp = os.path.join(os.path.dirname(__file__),
                                           config.get("Defaults",
                                               "PATHWAY_SUM_GENE_MAP_FP"))

    output_dir = os.path.dirname(pathway_sum_fp)
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    genes = codes3d.parse_summary_file(args.summary_file,
        args.buffer_size)

    codes3d.produce_pathway_summary(genes, pathway_db_fp,
        pathway_db_gene_name_fp, pathway_sum_fp, pathway_sum_gene_map_fp,
        args.buffer_size, args.significant_expression)

