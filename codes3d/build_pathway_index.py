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

    hgnc_gene_id_fp = os.path.join(os.path.dirname(__file__),
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

    codes3d.build_pathway_db(hgnc_gene_id_fp, pathway_db_fp, pathway_db_tmp_fp,
                             pathway_db_tab_fp, pathway_db_log_fp)


