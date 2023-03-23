#!/bin/bash

while getopts v:t:f:p: flag
do
    case "${flag}" in
        v) vcf_in=${OPTARG};;
        t) threads=${OPTARG};;
        f) filter_opts=${OPTARG};;
        p) vcflib_PATH=${OPTARG};;
    esac
done

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

#run filtering step on vcf
parallel -j $threads \
"$vcflib_PATH/vcffilter -f $filter_opts header_{.}.vcf > {.}_temp.recode.vcf" ::: \
` find ./ -name "header*.vcf" | cut -d "_" -f2 | sort `  

#restrip headers
parallel -j $threads \
'grep -v "^#" {/.}_temp.recode.vcf > {/.}_post_fil.vcf' ::: \
` find ./ -name "x*_temp*.vcf" | cut -d "_" -f1 | sort ` 
rm -f x*temp*vcf header*vcf

cat *post_fil.vcf | cat header - > filtered_out.vcf
rm header x*vcf