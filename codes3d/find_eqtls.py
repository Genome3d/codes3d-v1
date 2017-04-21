#!/usr/bin/env python

import argparse,codes3d,configparser,os

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-s","--snps_file",required=True,help="The list of SNPs generated by process_inputs.py")
	parser.add_argument("-g","--genes_file",required=True,help="The list of genes generated by find_genes.py")
	parser.add_argument("-c","--config",default=os.path.join(os.path.dirname(__file__),"../docs/codes3d.conf"),help="The configuration file to be used in this instance (default: conf.py)")
	parser.add_argument("-l","--local_databases_only",action="store_true",default=False,help="Consider only local databases. Will only include cis-eQTLs if using downloadable GTEx dataset.")
	parser.add_argument("-p","--num_processes",type=int,default=1,help="Desired number of processes for multiprocessing (default: 1).")
	parser.add_argument("-f","--fdr_threshold",type=float,default=0.05,help="The FDR threshold to consider an eQTL statistically significant (default: 0.05).")
	parser.add_argument("-o","--output_dir",default="hiCquery_output",help="The directory in which to output results (\"hiCquery_output\" by default).")
	args = parser.parse_args()
	config = configparser.ConfigParser();
	config.read(args.config)
	eqtl_data_dir = config.get("Defaults","EQTL_DATA_DIR")
	gene_database_fp = config.get("Defaults","GENE_DATABASE_FP")
	if not os.path.isdir(args.output_dir):
		os.makedirs(args.output_dir)
	snps = codes3d.parse_snps_files([args.snps_file])
	genes = codes3d.parse_genes_files([args.genes_file])
	eqtls,num_sig = codes3d.find_eqtls(snps,genes,eqtl_data_dir,gene_database_fp,args.fdr_threshold,args.local_databases_only,args.num_processes,args.output_dir)