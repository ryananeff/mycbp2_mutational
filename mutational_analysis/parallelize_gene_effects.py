from tqdm import tqdm
import varcode
import itertools
import pyensembl
from pyensembl import EnsemblRelease
import pysam
import os.path, json, sys, pickle
import pandas as pd

def merge_intervals(intervals):
	"""
	A simple algorithm can be used:
	1. Sort the intervals in increasing order
	2. Push the first interval on the stack
	3. Iterate through intervals and for each one compare current interval
	with the top of the stack and:
	A. If current interval does not overlap, push on to stack
	B. If current interval does overlap, merge both intervals in to one
	and push on to stack
	4. At the end return stack
	"""
	sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
	merged = []

	for higher in sorted_by_lower_bound:
		if not merged:
			merged.append(higher)
		else:
			lower = merged[-1]
			# test for intersection between lower and higher:
			# we know via sorting that lower[0] <= higher[0]
			if higher[0] <= lower[1]:
				upper_bound = max(lower[1], higher[1])
				merged[-1] = (lower[0], upper_bound)  # replace by merged interval
			else:
				merged.append(higher)
	return merged

gene_list_file = open(sys.argv[1], "r")
outfile = open(sys.argv[2],"w")

ensembldb = pd.read_csv("./data/Homo_sapiens.GRCh37.75.gtf.expanded.csv")
ref = pysam.FastaFile('./data/human_g1k_v37.fa')
gene_list = []
for line in gene_list_file:
	gene_list.append(line.strip())

es = EnsemblRelease(75)
print "starting..."
#pbar = ProgressBar(len(ensembldb[ensembldb["source"]=="protein_coding"].groupby("gene_name")))
for n,g in ensembldb[(ensembldb["source"]=="protein_coding")&(ensembldb["gene_name"].isin(gene_list))].groupby("gene_name"):
	print "Processing gene " + n
	effects = []
	eff = ''
	a = g[g['feature']=='exon'][['feature', 'seqname', 'start', 'end', 'strand']].sort(columns=["start"])
	if len(a.values) == 0:
		print "no exons"
		continue
	if str(a.values[0][1]) not in ref.references:
		gene_set = set(a["seqname"].astype(str).values)
		for name in gene_set:
			if name not in ref.references:
				a = a[a["seqname"] != name]
		if len(a.values) == 0:
			print "not in ref", gene_set
			continue
	intervals = a[["start","end"]].values
	merged = merge_intervals(intervals)
        pbar=tqdm(total=sum([y-x for x,y in merged]))
	for start, stop in merged:
		gene_seq = ref.fetch(str(a.values[0][1]),start=start,end=stop)
		for base in range(start,stop):
			pbar.update()
			try:
				ref_allele = gene_seq[base-start]
			except:
				"ref error"
				continue
			if ref_allele=="N":
				"base error"
				continue
			for alt in ["A","C","T","G"]:
				if alt!=ref_allele:
					eff = varcode.Variant(a.values[0][1],base+1,ref_allele,alt,ensembl=es).effects()
					effstr = str(type(eff.top_priority_effect()))[24:-2]
					effects.append([n, a.values[0][1], base,ref_allele, alt, effstr])
					del eff
					if effstr not in ["Silent","Substitution","StopLoss","PrematureStop",
								 "SpliceDonor","SpliceAcceptor","ExonicSpliceSite","IntronicSpliceSite","ExonLoss"]:
						break # skip these
	pbar.close()
	for line in effects:
		outfile.write("\t".join([str(l) for l in line]) + "\n")
print "done"
gene_list_file.close()
outfile.close()
