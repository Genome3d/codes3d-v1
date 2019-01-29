#!/usr/bin/env python

import argparse
import codes3d
import configparser
import logging
import os
import psutil

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-c", "--config",
                        default=os.path.join(os.path.dirname(__file__),
                                             "../docs/codes3d.conf"),
                        help="The configuration file to be used.")
    parser.add_argument("-k", "--num_processes_summary", type=int,
                        default=(psutil.cpu_count() // 2),
                        help="The maximum number of processes for data frame access." +\
                        "(default: %s)." % str((psutil.cpu_count() // 2)))
    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config)

    num_processes = args.num_processes_summary

    pathway_db_fp = os.path.join(os.path.dirname(__file__),
                                 config.get("Defaults", "pathway_db_fp"))
    pathway_db_tmp_fp = os.path.join(os.path.dirname(__file__),
                                     config.get("Defaults",
                                                "pathway_db_tmp_fp"))
    pathway_db_tsv_pw_fp = os.path.join(os.path.dirname(__file__),
                                        config.get("Defaults",
                                                   "pathway_db_tsv_pw_fp"))
    pathway_db_tsv_gene_exp_fp = os.path.join(os.path.dirname(__file__),
                                         config.get("Defaults",
                                                 "pathway_db_tsv_gene_exp_fp"))
    pathway_db_tsv_prot_exp_fp = os.path.join(os.path.dirname(__file__),
                                         config.get("Defaults",
                                                 "pathway_db_tsv_prot_exp_fp"))
    pathway_db_log_fp = os.path.join(os.path.dirname(__file__),
                                     config.get("Defaults",
                                                "pathway_db_log_fp"))
    expression_table_fp = os.path.join(os.path.dirname(__file__),
                                          config.get("Defaults",
                                                     "expression_table_fp"))
    pathway_db_gene_map_fp = os.path.join(os.path.dirname(__file__),
                                          config.get("Defaults",
                                                     "pathway_db_gene_map_fp"))
    pathway_db_exp_prot_fp = os.path.join(os.path.dirname(__file__),
                                          config.get("Defaults",
                                                     "pathway_db_exp_prot_fp"))
    pathway_db_gene_sym_fp = os.path.join(os.path.dirname(__file__),
                                    config.get("Defaults",
                                               "pathway_db_gene_sym_fp"))
    pathway_db_gene_syn_fp = os.path.join(os.path.dirname(__file__),
                                    config.get("Defaults",
                                               "pathway_db_gene_syn_fp"))
    pathway_db_gene_name_fp = os.path.join(os.path.dirname(__file__),
                                    config.get("Defaults",
                                               "pathway_db_gene_name_fp"))
    codes3d.build_pathway_db(pathway_db_fp, pathway_db_tmp_fp,
        pathway_db_tsv_gene_exp_fp, pathway_db_tsv_prot_exp_fp,
        pathway_db_tsv_pw_fp, pathway_db_log_fp,
        expression_table_fp, pathway_db_gene_map_fp, pathway_db_exp_prot_fp,
        pathway_db_gene_sym_fp, pathway_db_gene_syn_fp,
        pathway_db_gene_name_fp, num_processes)

