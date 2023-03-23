#!/bin/bash

while getopts v:t:f: flag
do
    case "${flag}" in
        v) vcf_in=${OPTARG};;
        t) threads=${OPTARG};;
        f) filter_opts=${OPTARG};;
    esac
done

# vcf_in=montipora_SNPs.vcf
# threads=30
# filter_opts="--minQ 30"

#grab the header
head -n 1000 $vcf_in | grep "^#" > header

#grab the non header lines
grep -v "^#" $vcf_in > variants

#split for each thread
split -n l/$threads --additional-suffix=.vcf variants
rm -f variants

#reattach headers to split vcfs
parallel -j $threads \
'cat header {/.}.vcf > header_{/.}.vcf' ::: \
` find ./ -name "x*.vcf" `
rm -f x*.vcf

echo "$filter_opts $(date) ##############" >> .vcftools.log
#run filtering step on vcf
parallel -j $threads \
vcftools --vcf header_{/.}.vcf $filter_opts --recode --recode-INFO-all --out {/.}_temp 2> .vcfparallel.log ::: \
` find ./ -name "header*.vcf" | cut -d "_" -f2 | sort `

#restrip headers
parallel -j $threads \
'grep -v "^#" {/.}_temp.recode.vcf > {/.}_post_fil.vcf' ::: \
` find ./ -name "x*_temp*.vcf" | cut -d "_" -f1 | sort ` 
rm -f x*temp*vcf header*vcf

cat *post_fil.vcf | cat header - > filtered_out.vcf
rm header x*vcf