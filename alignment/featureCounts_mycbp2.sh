#!/bin/bash
#BSUB -J MYCBP2.$1.RNA.featureCounts
#BSUB -W 48:00
#BSUB -n 12
#BSUB -q premium
#BSUB -P acc_PBG
#BSUB -o /sc/orga/projects/zhangb03a/neffr01/MycBP2_data/processed/logs/featureCounts/$1.log.out
#BSUB -e /sc/orga/projects/zhangb03a/neffr01/MycBP2_data/processed/logs/featureCounts/$1.log.err
#BSUB -u ryan.neff@icahn.mssm.edu
#BSUB -R rusage[mem=4012]
#BSUB -R span[hosts=1]
#BSUB -L /bin/bash
#BSUB -w "done(MYCBP2.$1.RNA.star)"

cd /sc/orga/projects/zhangb03a/neffr01/MycBP2_data/processed/featureCounts/
mkdir $1
cd $1
aligndir="/sc/orga/projects/zhangb03a/neffr01/MycBP2_data/processed/star/$1"
countdir="/sc/orga/projects/zhangb03a/neffr01/MycBP2_data/processed/featureCounts/$1"

module load subread

featureCounts -T 10 \
        -t exon \
        -g transcript_id \
        -p \
        -a /sc/orga/projects/PBG/REFERENCES/hg19/ensembl/Homo_sapiens.GRCh37.75.processed.gtf \
        -o $countdir/$1.exon.transcriptID.txt \
        $aligndir/bams/$1.accepted_hits.sort.coord.bam

featureCounts -T 10 \
        -t exon \
        -g gene_id \
        -p \
        -a /sc/orga/projects/PBG/REFERENCES/hg19/ensembl/Homo_sapiens.GRCh37.75.processed.gtf \
        -o $countdir/$1.exon.geneID.txt \
        $aligndir/bams/$1.accepted_hits.sort.coord.bam

featureCounts -T 10 \
        -t exon \
        -f -O \
        -p \
        -a /sc/orga/projects/PBG/REFERENCES/hg19/ensembl/Homo_sapiens.GRCh37.75.processed.gtf \
        -o $countdir/$1.exon.txt \
        $aligndir/bams/$1.accepted_hits.sort.coord.bam

featureCounts -T 10 \
        -t exon \
        -g gene_id \
        --primary -O \
        -p \
        -a /sc/orga/projects/PBG/REFERENCES/hg19/ensembl/Homo_sapiens.GRCh37.75.processed.gtf \
        -o $countdir/$1.primary.txt \
        $aligndir/bams/$1.accepted_hits.sort.coord.bam
