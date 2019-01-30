#!/usr/bin/env python

import argparse
import codes3d
import configparser
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
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

    args = parser.parse_args()
    config = configparser.ConfigParser()
    config.read(args.config)
    pathway_sum_fp = os.path.join(os.path.dirname(__file__),
                                  config.get("Defaults",
                                             "PATHWAY_SUM_FP"))
    pathway_plot_fp = os.path.join(os.path.dirname(__file__),
                                   config.get("Defaults",
                                              "PATHWAY_PLOT_FP"))
    if not os.path.isdir(pathway_plot_fp):
        os.mkdir(pathway_plot_fp)

    codes3d.plot_pathway_summary(pathway_sum_fp, args.buffer_size,
                                 pathway_plot_fp)

