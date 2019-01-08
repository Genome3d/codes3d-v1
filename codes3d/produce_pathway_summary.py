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
            "-o","--output_dir",default="codes3d_pathway_summary",
            help="The directory in which to output results "+\
            "(\"codes3d_pathway_summary\" by default).")
    parser.add_argument(
            "-c","--config",
        default=os.path.join(os.path.dirname(__file__),
                             "../docs/codes3d.conf"),
            help="The configuration file to be used in this "+\
            "instance (default: conf.py)")
    parser.add_argument(
            "-b","--buffer_size_in",type=int,default=1048576,
            help="Buffer size applied to file input during compilation "+\
            " (default: 1048576).")
    parser.add_argument(
            "-d","--buffer_size_out",type=int,default=1048576,
            help="Buffer size applied to file output during compilation "+\
            " (default: 1048576).")
    parser.add_argument(
            "-r", "--num_processes_summary", type=int,
            default=min(psutil.cpu_count(), 32),
            help="The number of processes for compilation of the results " +\
            "(default: %s)." % str(min(psutil.cpu_count(), 32)))
    parser.add_argument("-e", "--significant_expression", type=float,
            default=0.05,
            help="p-value of significant expression variation "+\
            "(default: 0.05).")

    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    pathway_db_fp = os.path.join(os.path.dirname(__file__),
                                 config.get("Defaults", "PATHWAY_DB_FP"))
    if not os.path.isdir(args.output_dir):
	    os.mkdir(args.output_dir)

    genes, snps = codes3d.parse_summary_file(args.summary_file,
                      args.buffer_size_in)

    codes3d.produce_pathway_summary(genes, snps, pathway_db_fp, args.output_dir,
        args.buffer_size_out, args.num_processes_summary,
        args.significant_expression)

