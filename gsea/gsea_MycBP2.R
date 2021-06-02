#!/usr/bin/env Rscript
## module load R/3.4.3
## using: ~/.Rlib

set.seed(12345)

library(GOtest) ##from minghui
library(msigdb) ##from minghui
#library(tidyr)
#library(dplyr)

options(stringsAsFactors=FALSE)

GOsets = c('C5.BP','C5.CC', 'C5.MF')
gosets_genes = msigdb.genesets(sets=GOsets, type='symbols', species='human', return.data.frame=T)

universe = universe=curated.genesets(c('HGNC_universe'))$Gene

####

setwd("/sc/orga/projects/zhangb03a/neffr01/MycBP2_data/processed/DEGs/")

print("changed to directory, loading data...")
diff_exp_cisplatin = read.table("mycbp2_si6cis_vs_controlcis.031519.deg.tsv",header=T)
diff_exp_cisplatin["Gene"] = rownames(diff_exp_cisplatin)
diff_exp_cisplatin[,"Phenotype"] = "si6cis_vs_controlcis" 
diff_exp_cisplatin  = diff_exp_cisplatin[,c("Gene","Phenotype","t")]
diff_exp_cisplatin = diff_exp_cisplatin[diff_exp_cisplatin[,3]!=0,]
colnames(diff_exp_cisplatin) = c("Gene","Phenotype","Z")

diff_exp_normal = read.table("mycbp2_si6_vs_control.031519.deg.tsv",header=T)
diff_exp_normal["Gene"] = rownames(diff_exp_normal)
diff_exp_normal[,"Phenotype"] = "si6_vs_control" 
diff_exp_normal  = diff_exp_normal[,c("Gene","Phenotype","t")]
diff_exp_normal = diff_exp_normal[diff_exp_normal[,3]!=0,]
colnames(diff_exp_normal) = c("Gene","Phenotype","Z")

diff_exp_control = read.table("mycbp2_controlcis_vs_control.DEG.031519.deg.tsv",header=T)
diff_exp_control["Gene"] = rownames(diff_exp_control)
diff_exp_control[,"Phenotype"] = "controlcis_vs_control" 
diff_exp_control = diff_exp_control[,c("Gene","Phenotype","t")]
diff_exp_control = diff_exp_control[diff_exp_control[,3]!=0,]
colnames(diff_exp_control) = c("Gene","Phenotype","Z")

diff_exp_si6 = read.table("mycbp2_si6cis_vs_si6.DEG.031519.deg.tsv",header=T)
diff_exp_si6["Gene"] = rownames(diff_exp_si6)
diff_exp_si6[,"Phenotype"] = "si6cis_vs_si6" 
diff_exp_si6 = diff_exp_si6[,c("Gene","Phenotype","t")]
diff_exp_si6 = diff_exp_si6[diff_exp_si6[,3]!=0,]
colnames(diff_exp_si6) = c("Gene","Phenotype","Z")

print("loaded data.")

name = "MycBP2_MSigDB"
perms = 15

print("testing cisplatin...")
result_weight_cisplatin = GOtest(x=diff_exp_cisplatin, go=gosets_genes, 
                                 name.x='Cisplatin', name.go='MSigDB', method='GSEA', permutations=perms)
print("writing cisplatin results...")
write.table(result_weight_cisplatin,"si6cis_vs_controlcis_gsea_enrichment.tsv")
print("testing normal treatment...")
result_weight_normal = GOtest(x=diff_exp_normal, go=gosets_genes, 
                              name.x='Normal', name.go='MSigDB', method='GSEA', permutations=perms)
print("writing normal results...")
write.table(result_weight_normal,"si6_vs_control_gsea_enrichment.tsv")
print("testing control cells...")
result_weight_control = GOtest(x=diff_exp_control, go=gosets_genes,
                               name.x='Control', name.go='MSigDB', method='GSEA', permutations=perms)
print("writing control results...")
write.table(result_weight_control,"controlcis_vs_control_gsea_enrichment.tsv")
print("testing si6 cells...")
result_weight_si6 = GOtest(x=diff_exp_si6, go=gosets_genes, 
                           name.x='si6', name.go='MSigDB', method='GSEA', permutations=perms)
print("writing si6 results...")
write.table(result_weight_si6,"si6cis_vs_si6_gsea_enrichment.tsv")

print("all results generated. saving image...")
save.image("GSEA_results.Rsave")
print("imaged saved. generating plots...")

pdf("result_weight_cisplatin_gsea_top30.pdf")
plotGseaEnrTable(GseaTable=result_weight_cisplatin[1:30,], x=diff_exp_cisplatin, go=gosets_genes)
dev.off()

pdf("result_weight_normal_gsea_top30.pdf")
plotGseaEnrTable(GseaTable=result_weight_normal[1:30,], x=diff_exp_normal, go=gosets_genes)
dev.off()

pdf("result_weight_control_gsea_top30.pdf")
plotGseaEnrTable(GseaTable=result_weight_control[1:30,], x=diff_exp_control, go=gosets_genes)
dev.off()

pdf("result_weight_si6_gsea_top30.pdf")
plotGseaEnrTable(GseaTable=result_weight_si6[1:30,], x=diff_exp_si6, go=gosets_genes)
dev.off()

print("all done!")