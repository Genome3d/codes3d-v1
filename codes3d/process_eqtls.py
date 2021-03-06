#!/usr/bin/env python

import argparse
import codes3d
import configparser
import os

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-s","--snps_file",nargs='+',required=True,help="The list of SNPs generated by process_inputs.py")
	parser.add_argument("-g","--genes_file",nargs='+',required=True,help="The list of genes generated by find_genes.py")
	parser.add_argument("-e","--eqtls_file",nargs='+',help="Space-separated list of one or more eQTL files generated by find_eqtls.py")
	parser.add_argument("-o","--output_dir",default="codes3d_output",help="The directory in which to output results (\"codes3d_output\" by default).")
	parser.add_argument("-c","--config",default="docs/conf.py",help="The configuration file to be used in this instance (default: conf.py)")
	args = parser.parse_args()
	config = configparser.ConfigParser()
	config.read(args.config)
	gene_database_fp = config.get("Defaults","GENE_DATABASE_FP")
	if not os.path.isdir(args.output_dir):
		os.makedirs(args.output_dir)
	
	print("Parsing SNPs...")
	snps = codes3d.parse_snps_files(args.snps_file)
	print("Parsing genes...")
	genes = codes3d.parse_genes_files(args.genes_file)
	print("Parsing eQTLs...")
	eqtls,num_sig = codes3d.parse_eqtls_files(args.eqtls_file)
	print("Re-processing eQTLs...")
	codes3d.process_eqtls(snps,genes,eqtls,gene_database_fp)
	print("Writing new eQTL file...")
	eqtlfile = open(args.output_dir + '/eqtls_corrected.txt','w')
	for snp in eqtls.keys():
		for gene in eqtls[snp].keys():
			if not gene == "snp_info":
				for tissue in eqtls[snp][gene]["tissues"].keys():
					eqtlfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (snp,eqtls[snp]["snp_info"]["chr"],eqtls[snp]["snp_info"]["locus"],gene,eqtls[snp][gene]["gene_chr"],eqtls[snp][gene]["gene_start"],eqtls[snp][gene]["gene_end"],eqtls[snp][gene]["cell_lines"],eqtls[snp][gene]["cis?"],eqtls[snp][gene]["p_thresh"],tissue,eqtls[snp][gene]["tissues"][tissue]["pvalue"]))
	