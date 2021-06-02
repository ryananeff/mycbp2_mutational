library(limma)
library(edgeR)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

#### INPUTS ###

# Read the data into R
seqdata <- read.delim("/sc/orga/projects/zhangb03a/neffr01/MycBP2_data/processed/gene_expression/MycBP2.gene_exp.primary.tsv", 
			stringsAsFactors = FALSE,sep="\t")

rownames(seqdata) = make.names(seqdata[,2],unique=T)
seqdata[,1:2] = NULL

# Read the sample information into R
sampleinfo <- read.delim("/sc/orga/projects/zhangb03a/neffr01/MycBP2_data/processed/gene_expression/mycbp2_meta.corrected.txt")

#create DGElist object
countdata <- seqdata
rownames(sampleinfo) = sampleinfo$SampleID
group = sampleinfo[colnames(seqdata),]$group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
y <- DGEList(counts=countdata,group = group)

#filter low CPM genes
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]

# normalize
y <- calcNormFactors(y,method="TMM")
y <- estimateDisp(y,design = design)
print(y$samples)

cpm <- cpm(y)
lcpm <- cpm(y, log=TRUE)

col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(y, labels=group, col=col.group)
title(main="A. Sample groups")

contr.matrix <- makeContrasts(
  si6cis_vs_controlcis = si6_Cis-control_Cis,
  si6_vs_control = si6-control,
  si6cis_vs_si6 = si6_Cis-si6,
  controlcis_vs_control = control_Cis-control,
  levels = colnames(design))
contr.matrix

v <- voom(y, design, plot=TRUE)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)

summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)

head(tfit$genes$SYMBOL[de.common], n=20)

vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))

si6cis_vs_controlcis <- topTreat(tfit, coef=1, n=Inf)
si6_vs_control <- topTreat(tfit, coef=2, n=Inf)
si6cis_vs_si6 <- topTreat(tfit, coef=3, n=Inf)
controlcis_vs_control <- topTreat(tfit, coef=4, n=Inf)

head(si6cis_vs_controlcis)
head(si6_vs_control)
head(si6cis_vs_si6)
head(controlcis_vs_control)

plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))
plotMD(tfit, column=2, status=dt[,2], main=colnames(tfit)[2], xlim=c(-8,13))

library(gplots)
si6cis_vs_controlcis.topgenes <- rownames(si6cis_vs_controlcis)[1:100]
i <- match(sort(si6cis_vs_controlcis.topgenes),rownames(lcpm))
mycol <- colorpanel(1000,"blue","white","red")

heatmap.2(lcpm[i,], scale="row",
          labRow=rownames(lcpm)[i], labCol=group,
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

####################

setwd("/sc/orga/projects/zhangb03a/neffr01/MycBP2_data/processed/")

write.table(si6cis_vs_controlcis, 
            file=paste0("DEGs/mycbp2_si6cis_vs_controlcis",".031519.deg.tsv"),
            sep="\t",quote=F,row.names = T)

write.table(si6_vs_control, 
            file=paste0("DEGs/mycbp2_si6_vs_control",".031519.deg.tsv"),
            sep="\t",quote=F,row.names = T)

write.table(si6cis_vs_si6, 
            file=paste0("DEGs/mycbp2_si6cis_vs_si6",".031519.deg.tsv"),
            sep="\t",quote=F,row.names = T)

write.table(controlcis_vs_control, 
            file=paste0("DEGs/mycbp2_controlcis_vs_control",".031519.deg.tsv"),
            sep="\t",quote=F,row.names = T)
