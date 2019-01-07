#!/usr/bin/env python

import argparse
import codes3d
import configparser
import logging
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-c", "--config",
                        default=os.path.join(os.path.dirname(__file__),
                                             "../docs/codes3d.conf"),
                        help="The configuration file to be used.")

    args = parser.parse_args()

    config = configparser.ConfigParser()
    config.read(args.config)

    hgnc_gene_sym_fp = os.path.join(os.path.dirname(__file__),
                                    config.get("Defaults",
                                               "pathway_db_gene_sym_fp"))
    pathway_db_fp = os.path.join(os.path.dirname(__file__),
                                 config.get("Defaults", "pathway_db_fp"))
    pathway_db_tmp_fp = os.path.join(os.path.dirname(__file__),
                                     config.get("Defaults",
                                                "pathway_db_tmp_fp"))
    pathway_db_tsv_exp_fp = os.path.join(os.path.dirname(__file__),
                                         config.get("Defaults",
                                                    "pathway_db_tsv_exp_fp"))
    pathway_db_tsv_pw_fp = os.path.join(os.path.dirname(__file__),
                                        config.get("Defaults",
                                                   "pathway_db_tsv_pw_fp"))
    pathway_db_log_fp = os.path.join(os.path.dirname(__file__),
                                     config.get("Defaults",
                                                "pathway_db_log_fp"))
    pathway_db_gene_map_fp = os.path.join(os.path.dirname(__file__),
                                          config.get("Defaults",
                                                     "pathway_db_gene_map_fp"))
    pathway_db_exp_gene_fp = os.path.join(os.path.dirname(__file__),
                                          config.get("Defaults",
                                                     "pathway_db_exp_gene_fp"))
    pathway_db_exp_pept_fp = os.path.join(os.path.dirname(__file__),
                                          config.get("Defaults",
                                                     "pathway_db_exp_pept_fp"))
    pathway_db_exp_prot_fp = os.path.join(os.path.dirname(__file__),
                                          config.get("Defaults",
                                                     "pathway_db_exp_prot_fp"))

    codes3d.build_pathway_db(pathway_db_fp, pathway_db_tmp_fp,
                             pathway_db_tsv_exp_fp, pathway_db_tsv_pw_fp,
                             pathway_db_log_fp, pathway_db_gene_map_fp,
                             pathway_db_exp_gene_fp, pathway_db_exp_pept_fp,
                             pathway_db_exp_prot_fp, hgnc_gene_sym_fp)

