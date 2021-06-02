#!/usr/bin/env python

import vcf as pyvcf
import pandas as pd
import numpy as mp
import logging
import coloredlogs
import pandas as pd
import sys, getopt

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', fmt="%(asctime)s EDM.py [%(process)d] %(levelname)s %(message)s")

# parse command line arguments and run function
def main(argv):
    infile=None
    outfile=None
    try:
        opts, args = getopt.getopt(argv,"hi:o:l:",["help","ifile=","ofile=", "log="])
    except getopt.GetoptError:
        print 'convert_vep_vcf_to_table.py -i <in.maf> -o <out.tsv> (--log DEBUG)'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            help='''  Usage: TODO
            '''
            print 'convert_vep_vcf_to_table.py -i <in.maf> -o <out.tsv> (--log DEBUG)'
            print help
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-o", "--ofile"):
            outfile = arg
        elif opt in ("-l", "--log"):
            loglevel = arg
            logger.setLevel(loglevel)
    if not infile or not outfile:
        logger.critical("Need to specify input and output files. Stop. Get help with -h")
        return 1
    # run script
    outfile = open(outfile, 'w')
    status = convert_vep_vcf_to_table(infile, outfile)
    outfile.close()
    return status

def convert_vep_vcf_to_table(infile, outfile):
    vcffile = pyvcf.VCFReader(filename=infile)
    samples = vcffile.samples
    vcffile_fp = open(infile,'r')
    count = 0
    headers = ["gene", "gene_type", "var_type", "chrom", "pos", "ref", "alt"]
    headers.extend(samples)
    outfile.write("\t".join(headers) + "\n")
    logger.info('starting')
    for line in vcffile_fp:
        if "#" in line: continue
        if (count % 10000 == 0) & (count != 0):
            logger.info('processed %s', str(count))
        l = line.split('\t',10)
        chrom = l[0]
        pos = l[1]
        allele_ref = l[3]
        allele_alt = l[4]
        info_line = l[7].split(",")
        for i in info_line:
            split_effects = i.split("|",8)
            gene_name = split_effects[3]
            variant_effect = split_effects[1]
            gene_type = split_effects[7]
            if "protein_coding" != gene_type:
                continue
            if gene_name == "":
                continue
            count += 1
            recorddict = dict()
            recorddict["gene"] = gene_name
            recorddict["gene_type"] = gene_type
            recorddict["var_type"] = variant_effect
            recorddict["chrom"] = chrom
            recorddict["pos"] = pos
            recorddict["ref"] = allele_ref
            recorddict["alt"] = allele_alt
            l = line.split('\t') # split the whole line
            gts = l[9:]
            gts_proc = [eval(g.split(":")[0].replace("/","+").replace('.','-2')) for g in gts]
            gts_proc2 = [i if i >= 0 else None for i in gts_proc]
            record_gt = {s:g for s,g in zip(samples,gts_proc2)}
            recorddict.update(record_gt)
            outfile.write("\t".join([str(recorddict[i]) for i in headers]) + "\n")
            break

if __name__ == "__main__":
   main(sys.argv[1:])