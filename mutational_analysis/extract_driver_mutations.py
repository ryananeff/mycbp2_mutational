#!/usr/bin/env python

# Ryan Neff
# Zhang Lab
# June 7th, 2016

# Notes:
#
# run by typing in "python extract_driver_mutations.py -i in.maf -o out.tsv"
#
# this program uses Ensembl Release 75 data (hg19) - do a Ctrl-F for "EnsemblRelease" to change it to hg38 version!!
#
# if there are missing packages, install them with 'pip install <packagename>'.
#	* You may want to try the Anaconda Python distribution for running this, available online.
#
# you will also likely need to download the pyensembl data - do this by running 'pyensembl install --release 75' into the terminal.

# imports 
import varcode, pyensembl
from pyensembl import EnsemblRelease
import pysam
import numpy as np
import scipy as sp
import pandas as pd
from scipy import stats
import itertools
import os.path, sys
import json, pickle
import sys, getopt
import coloredlogs, logging
import math

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', fmt="%(asctime)s EDM.py [%(process)d] %(levelname)s %(message)s")

# parse command line arguments and run function
def main(argv):
   infile=None
   outfile=None
   ref = './data/human_g1k_v37.fa'
   possible = './data/kmer_muts_all_genes_new_2.pickle'
   try:
      opts, args = getopt.getopt(argv,"hi:o:rpl",["ifile=","ofile=", "ref=", "possible=", "log="])
   except getopt.GetoptError:
      print 'extract_driver_mutations.py -i <in.maf> -o <out.tsv> (-r <ref.fa>) (-p <possible_mutations.pickle>)'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         help='''  Usage: Calls somatic driver mutations (SNVs only at the moment) from MAF files using a dN/dS method.
    inputs:
        input.maf
            variant file for input
        ref.fa
            reference genome in FASTA format (will auto-index if index not present)
        possible_mutations.pickle
            File containing all possible mutations for every gene in the reference. 
            Generated with parallelize_gene_effects.py script on the Minerva cluster.
        covariates_file.tsv (optional, not used at the time)
            file with useful gene level covariates to calculate the gene_specific mutation rate

    outputs:
        gene_liklihoods.tsv
            table of likelihoods with following columns:
                gene_name, n_syn, n_mis, n_non, n_splice, w_mis, w_non, w_splice, q_mis, q_non, q_splice, q_all
    '''

         print 'extract_driver_mutations.py -i <in.maf> -o <out.tsv> (-r <ref.fa>) (-p <possible_mutations.pickle>)'
         print help
         sys.exit()
      elif opt in ("-i", "--ifile"):
         infile = arg
      elif opt in ("-o", "--ofile"):
         outfile = arg
      elif opt in ("-r", "--ref"):
         ref = arg
      elif opt in ("-p", "--possible"):
         possible = arg
      elif opt in ("-l", "--log"):
	 loglevel = arg
	 logger.setLevel(loglevel)
   if not infile or not outfile:
       logger.critical("Need to specify input and output files. Stop. Get help with -h")
       return 1
   # run script
   mut_table = driver_gene_caller(infile, ref, possible)
   st = pd.DataFrame(mut_table).T
   st.index = st.index.map(lambda x: lookup(x))
   st.to_csv(outfile, sep="\t", na_rep="NULL")

#the DNA mutation types to consider. Tweak this based on what you want to focus on. These types are from the varcode python package (docs online)

muts_keys = ['FivePrimeUTR',
 'ExonicSpliceSite',
 'StopLoss',
 'SpliceDonor',
 'ThreePrimeUTR',
 'PrematureStop',
 'Substitution',
 'Intronic',
 'StartLoss',
 'AlternateStartCodon',
 'SpliceAcceptor',
 'IntronicSpliceSite']

sig_muts_keys = ['FivePrimeUTR',
 'ExonicSpliceSite',
 'StopLoss',
 'SpliceDonor',
 'ThreePrimeUTR',
 'PrematureStop',
 'Substitution',
 'StartLoss',
 'AlternateStartCodon',
 'SpliceAcceptor',
 'IntronicSpliceSite']

# main program. Note that "covariates" is not used currently.
def driver_gene_caller(infile, 
                       ref = './data/human_g1k_v37.fa', 
                       possible_mutations = './data/kmer_muts_all_genes_new_2.pickle',
                       covariates = None):
    logger.info("Starting...")
    sys.stdout.flush()
    
    # open the reference genome
    ref = pysam.FastaFile(ref)
    
    #Step 0 - generate the 192 mutation pairs
    all_3mers, all_3mer_combinations = get_possible_3mers()
    
    # Step 1 - calculate global rate of 3mer mutations in genome
    if os.path.isfile(possible_mutations):
        logger.info("Opening file containing global rate of 3mer mutations...")
        with open(possible_mutations, "rb") as fp: 
            kmer_muts_updated = pickle.load(fp)
    else:
        logger.error("Please provide a possible mutations file!")
        raise    
    
    logger.info("Getting 3mers around each mutation in input and classifying them...")
    
    # read input file
    varfile = pd.read_csv(infile, sep="\t", comment="#", usecols=['Hugo_Symbol', 
                               'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2'])
    
    # get strand mutations are on
    gene_strand = dict()
    ensembldb = pd.read_csv("./data/Homo_sapiens.GRCh37.75.gtf.expanded.csv")
    for n,g in ensembldb[(ensembldb["source"]=="protein_coding")].groupby("gene_name"):
        gene_strand[n] = str(g.iloc[0].strand)

    # calculate the context of each mutation and recalculate the mutation effect to match varcode prediction  
    kmer_muts_vars_updated = dict()
    es = EnsemblRelease(75) #using the hg19 ensembl release
    geneCount = 0
    for gene1,variants in varfile.groupby(varfile["Hugo_Symbol"]):
        outdict = dict(zip(all_3mer_combinations, [dict() for n in range(0,len(all_3mer_combinations))]))
        if len(variants["Chromosome"].unique()) > 1:
            logger.warning("gene %s spans more than one chromosome", gene1)
            continue
        chrom = str(variants["Chromosome"].values[0])
        start = min(variants["Start_Position"])
        stop = max(variants["End_Position"])
        start -= 2
        stop += 1
        try:
            gene = None
            gene_obj = es.genes_at_locus(chrom, start, end=stop)
            genes = [(str(g.name),g) for g in gene_obj if g.biotype == "protein_coding"]
            if genes == []:
                logger.debug("gene " + gene1 + " not protein coding, skipping")
                continue
            for g,obj in genes:
                if gene1 == g:
                    gene = gene1
                    gene_obj = obj
            if gene == None:
                gene = genes[0][0]
                gene_obj = genes[0][1]
        except IndexError:
            logger.warning("skipping gene '" + gene1 + "', not in index")
            continue
        direction = gene_strand[gene]
        if gene != gene1:
            logger.debug("renaming gene -> orig: " + gene1 + ", new: " + gene)
        if gene not in kmer_muts_updated: 
            if gene1 in kmer_muts_updated:
                gene = gene1
            else:
                logger.warning("can't find gene '"  + gene1 + "' or '" + gene + "'")
                continue
        if chrom not in ref.references: 
            logger.critical("input chrom does not match reference")
            continue
        r = ref.fetch(chrom, start=start, end=stop)    
        for ix,var in variants.iterrows():
            geneCount += 1
            if var["End_Position"] - var["Start_Position"] != 0: continue # only look at SNPs
            base = var["Start_Position"]-1
            altbase = var["Tumor_Seq_Allele2"]
            orig = ''
            # correct for strand
            orig = r[base-start-1:base-start+2]
            if direction == "-":
                orig = orig[::-1]
            # ignore variants with Ns around variant
            if "N" in orig: continue
            try:
                alt = orig[0]+altbase+orig[2]
            except IndexError:
                logger.critical("Something went wrong with getting the context for a variant.")
                logger.critical("Chrom: %s, gene: %s, start: %s, context: %s, altbase: %s", str(chrom), gene1, str(var["Start_Position"]), orig, altbase)
                raise Exception
            if orig == alt: continue
            if var["Reference_Allele"] not in ["A","C","T","G"]: continue
            if var["Reference_Allele"] != orig[1]:
               # raise Exception("ERROR: Ref allele in file does not match reference. Stop.")
            	logger.critical("ref allele in file does not match reference")
                continue
            if altbase not in ["A","C","T",'G']: continue # check for weirdness in alternate base
            try:
                peffect = varcode.Variant(chrom, base+1, orig[1],altbase,ensembl=es).effects() # classify variant effect
                peffect = str(peffect.top_priority_effect()).split("(")[0] # take top transcript and get top effect
                if peffect in outdict[(orig,alt)]:
                    outdict[(orig,alt)][peffect]+=1
                else:
                    outdict[(orig,alt)][peffect]=1
            except: 
                logger.warning("could not look up variant at pos (" + chrom + "," + str(base) + ")")
        if gene in kmer_muts_vars_updated:
            logger.info("collision detected at gene %s - merging dicts", gene)
            indict = kmer_muts_vars_updated[gene]
            for key in outdict.iterkeys():
                if key in indict:
                    for key2 in outdict[key]:
                        if key2 in indict[key]:
                            indict[key][key2] += outdict[key][key2]
                        else:
                            indict[key][key2] = outdict[key][key2]
                else:
                    indict[key] = outdict[key]
        else:
            kmer_muts_vars_updated[gene] = outdict
    logger.info("Got context for " + str(len(kmer_muts_vars_updated)) + " variants.")
    logger.debug("Input variants = " + str(len(varfile)))
    logger.debug("Genes analyzed = " + str(geneCount))
    if len(kmer_muts_vars_updated) == 0:
    	logger.critical("No variants analyzed.")
	raise
    
    logger.info("Calculating the mutation rate for each mutation type...")
    new_mut_rate = dict() 
    # maf file loading
    for gene,pairs in kmer_muts_vars_updated.iteritems(): #this was generated by putting all of the vars together
        for pair,variants in pairs.iteritems():
            if pair not in new_mut_rate:
                new_mut_rate[pair] = variants.copy()
            else:
                for vartype,value in variants.iteritems():
                    if vartype not in new_mut_rate[pair]:
                        new_mut_rate[pair][vartype] = value
                    else:
                        new_mut_rate[pair][vartype] += value
    
    # possible_variants file loading
    for gene,pairs in kmer_muts_updated.iteritems(): #we loaded this (all possible muts)
        for pair,variants in pairs.iteritems():
            for vartype,value in variants.iteritems():
		if pair not in new_mut_rate: continue #TODO - fix this!!
                if vartype not in new_mut_rate[pair]:
                    new_mut_rate[pair][vartype] = [0, value]
                else:
                    if type(new_mut_rate[pair][vartype]) == type(1):
                        new_mut_rate[pair][vartype] = [new_mut_rate[pair][vartype],value]
                    else:
                        new_mut_rate[pair][vartype][1] += value
    
    logger.info("Calculating the gene-level mutation rates...")
    gene_mut_rate = dict()
    muts_per_gene = kmer_muts_vars_updated
    mut_rates = new_mut_rate.copy()
    r = mut_rates # relative rate globally
    S = kmer_muts_updated #= pickle.load( open("kmer_muts_all_genes.pickle", 'rb') ) # max number possible for gene
    new_S = dict()
    for k in S:
        if S[k] != None:
            new_S[k] = S[k]
    for k in new_S:
        ms = muts_per_gene[k] if k in muts_per_gene else None
        if ms==None: continue
        s_gene = new_S[k]
        sites = 0
        #use silent mutations for the mutation rate prediction
        for pair in s_gene:
            if ("Silent" in s_gene[pair]) & ("Silent" in r[pair]):
                sites += s_gene[pair]["Silent"]*r[pair]["Silent"][0]/float(r[pair]["Silent"][1])
        obs_sites = 0
        for a,w in muts_per_gene[k].iteritems(): #loop over each pair in muts_per_gene[gene]
            if "Silent" in w:
                obs_sites += w["Silent"]
        if sites == 0:
            gene_mut_rate[k] = sites
        else:
            gene_mut_rate[k] = float(obs_sites)/sites
    
    logger.info("Calculating the expected vs observed mutation rates for each mutation type...")
    labels, x, y = [],dict(),dict()
    for mtype in muts_keys:
        x[mtype] = []
        y[mtype+"_obs"] = []
        for k in new_S:
            ms = muts_per_gene[k] if k in muts_per_gene else None
            if ms==None: continue
            s_gene = new_S[k]
            sites = 0
            # predict the amount of mutation for each mut_type
            for pair in s_gene:
                if (mtype in s_gene[pair]) & (mtype in r[pair]):
                    sites += s_gene[pair][mtype]*r[pair][mtype][0]/float(r[pair][mtype][1])*gene_mut_rate[k]
            labels.append(k)
            x[mtype].append(sites)
            obs_sites = 0
            for a,w in kmer_muts_vars_updated[k].iteritems(): #loop over each pair in muts_per_gene[gene]
                if mtype in w:
                    obs_sites += w[mtype]
            y[mtype+"_obs"].append(obs_sites)
    
    #combine everything into a dataframe
    data_out = pd.concat([pd.DataFrame(labels), pd.DataFrame(x), pd.DataFrame(y), ], join="inner", axis=1)
    data_out.set_index(data_out[0], inplace=True)
    del data_out[0]
    data_out["rate"] = [gene_mut_rate[i] for i in data_out.index]
    
    # step 5 - calculate significance
    logger.info("Calculating significance...")
    def find_sites(gene_name, muttype, S=S):
        subs = 0
        for pair,vals in S[gene_name].iteritems():
            if muttype in vals:
                subs += vals[muttype]
        return subs

    logger.warning("TODO:: Setting n_participants to 808 for ADNI dataset.")
    n_participants = 808 #len(set(varfile["Tumor_Sample_Barcode"]))
    threshold = 0.05/20500 # this is the BH p-value base threshold
    logger.info("Sig theshold: " + str(threshold))
    sig_sites = []
    sig_table = dict()
    for gene_name, row in data_out.iterrows():
        sig_table[gene_name] = dict()
	global_sig_value = 1
        for muttype in muts_keys:
            x = row[muttype+"_obs"]
            pred_sites = row[muttype]
            n_sites = find_sites(gene_name, muttype)
            if n_sites == 0:
                continue
            if pred_sites == 0:
                pred_sites = np.average(data_out[data_out[muttype] != 0][muttype])
            sig = sp.stats.binom_test(x, n=n_sites*n_participants, p=pred_sites/(n_sites*n_participants), 
                                alternative="greater") # using a negative binomial test
                                # x = (observed mutations in gene)/(all possible for gene)
                                # n = (number of possible mutations in gene)
                                # p = (predicted number of mutations)/(all possible for gene)
            sig_table[gene_name][muttype+"_p"] = sig
            sig_table[gene_name][muttype+"_obs"] = x
            sig_table[gene_name][muttype+"_n"] = n_sites*n_participants
            sig_table[gene_name][muttype+"_exp"] = pred_sites
            if (sig <= threshold) & (muttype in sig_muts_keys): # if significant
                logger.debug("Significant:\t%s\t%s\t%s", gene_name, muttype, str(sig))
                sig_sites.append([gene_name, muttype, sig])
	for muttype in sig_muts_keys:
	    if muttype+"_p" in sig_table[gene_name]:
	    	sig = sig_table[gene_name][muttype+"_p"]
            	if (sig < 0) | (sig > 1) | (math.isnan(sig)):
		     sig = 1
	    	global_sig_value *= sig
        sig_table[gene_name]["@ALL_p"] = global_sig_value
        if global_sig_value <= threshold: #if significant
            logger.debug("Significant:\t%s\t%s\t%s",gene_name,"ALL",str(global_sig_value))
            sig_sites.append([gene_name, "ALL", global_sig_value])
    logger.info("Completed.")
    return sig_table

def get_possible_3mers():
    '''purpose: generate the 192 mutation pairs'''
    bases = ['A', 'C', 'T', 'G']
    all_3mers = []
    all_3mer_combinations = []
    for combination in itertools.product(bases, repeat=3):
        all_3mers.append(''.join(map(str, combination)))
    for combination in itertools.product(all_3mers, all_3mers):
        if combination[0] == combination[1]: continue
        if combination[0][0] != combination[1][0]: continue
        if combination[0][2] != combination[1][2]: continue
        all_3mer_combinations.append(combination)
    return all_3mers, all_3mer_combinations 

# gene name converter
alias_to_gene = pd.read_csv("./data/alias_to_gene.txt", sep =" ", header=None)
def forward_translate(x):
    if x in alias_to_gene[0]:
        return alias_to_gene[alias_to_gene[0]==x][1]
    else:
        return None
def reverse_translate(x):
    if x in alias_to_gene[1]:
        return alias_to_gene[alias_to_gene[1]==x][0]
    else:
        return None
def lookup(x):
    if type(x) != type("string"):
        return x
    rev = reverse_translate(x)
    if rev:
        return x
    fwd = forward_translate(x)
    if fwd:
        return str(fwd)
    else:
        return x

if __name__ == "__main__":
   main(sys.argv[1:])
