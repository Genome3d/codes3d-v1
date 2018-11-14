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
            "-t", "--num_processes_summary", type=int,
            default=min(psutil.cpu_count(), 32),
            help="The number of processes for compilation of the results " +\
            "(default: %s)." % str(min(psutil.cpu_count(), 32)))
    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    pathway_db_fp = os.path.join(os.path.dirname(__file__),
                                 config.get("Defaults", "PATHWAY_DB_FP"))
    if not os.path.isdir(args.output_dir):
	    print('\tCreating output directory..')
	    os.mkdir(args.output_dir)
    gene_exp = codes3d.parse_summary_file(args.summary_file,
        args.buffer_size_in)
    codes3d.produce_pathway_summary(gene_exp, pathway_db_fp, args.output_dir,
            args.buffer_size_out, args.num_processes_summary)

