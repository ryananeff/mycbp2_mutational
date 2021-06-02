#!/bin/bash

#changed ensembl release to 75

announce_program() {
    # $1 == name of program
    echo "=============================="
    echo "$1 RUNNING"
    echo "=============================="
    echo
}

test_exit_status() {
    # $1 == name of program
    # $2 == program exit status (a.k.a. $?)
    if [ $2 -ne 0 ]
    then
        echo "=============================="
        echo "$1 FAILED with exit status $2" | tee /dev/stderr
        echo "=============================="
        echo
        exit 42
    else
        echo "=============================="
        echo "$1 SUCCEEDED with exit status $2"
        echo "=============================="
        echo
    fi
}
module purge all
module load star/2.4.0g1 #yes, this is older, but the reference files are also old
module load samtools
module load java
module load picard
cd /sc/orga/projects/zhangb03a/neffr01/MycBP2_data/processed/star/
mkdir $1
cd $1
export aligndir="/sc/orga/projects/zhangb03a/neffr01/MycBP2_data/processed/star/$1"

#    | samtools view -bS - 

# Run star to generate accepted_hits.bam file
announce_program star
STAR \
    --chimSegmentMin 15  \
    --chimJunctionOverhangMin 15 \
    --outSAMstrandField intronMotif \
    --genomeDir /sc/orga/projects/PBG/REFERENCES/hg19/star/2.4.0d/Homo_sapiens.GRCh37.ensembl75.overhang75bp \
    --sjdbGTFfile /sc/orga/projects/PBG/REFERENCES/hg19/ensembl/Homo_sapiens.GRCh37.75.processed.gtf \
    --runThreadN 10 \
    --outReadsUnmapped Fastx \
    --outStd SAM \
    --outSAMmode Full \
    --outFileNamePrefix $aligndir/accepted_hits \
    --readFilesCommand zcat \
    --readFilesIn $2,$3 \
     \
    1> $aligndir/accepted_hits.sam
#test_exit_status star $?

# Move tophat output to ``bams`` directory
mkdir -p $aligndir/bams
cd $aligndir/bams
mv ../accepted_hits.sam $1.accepted_hits.sam
cd -

# Sort and index the bam by coordinate, using karyotypic order
announce_program AddOrReplaceReadGroups.jar
java -jar /sc/orga/projects/PBG/scripts/picard/AddOrReplaceReadGroups.jar \
    I=$aligndir/bams/$1.accepted_hits.sam \
    O=$aligndir/bams/$1.accepted_hits.sort.coord.bam \
    SO=coordinate \
    ID=$1 \
    LB=$1 \
    PL=ILLUMINA \
    PU=HiSeq2500 \
    SM=$1 \
    CN=MSSM \
    DS=rnaseq
test_exit_status AddOrReplaceReadGroups.jar $?
announce_program BuildBamIndex.jar
java -jar $PICARD_HOME/BuildBamIndex.jar \
    I=$aligndir/bams/$1.accepted_hits.sort.coord.bam
test_exit_status BuildBamIndex.jar $?
announce_program RemoveOldBAM
    rm $aligndir/bams/$1.accepted_hits.sam
test_exit_status RemoveOldBAM $?
cp $aligndir/bams/$1.accepted_hits.sort.coord.bai $aligndir/bams/$1.accepted_hits.sort.coord.bam.bai

# Run other programs that use the sorted bam files (optional)