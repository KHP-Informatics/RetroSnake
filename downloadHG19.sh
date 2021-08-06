
#!/bin/bash

#Usage: bash downloadHG19.sh PATH_TO_DIR_TO_PLACE_REFERENCE_GENOME

HG19_DIR=$1
cd $HG19_DIR

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
tar -zxvf chromFa.tar.gz
for i in chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrY.fa chrX.fa chrM.fa; do cat $i >> hg19.fa ; rm $i ; done
rm chr*
samtools faidx hg19.fa

#in config.yaml update the path to the downloaded HG19.fa
