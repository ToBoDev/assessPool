#sets up new filtering file and outputs initial data
vcf_init <- function(vcf, wd, output){
  setwd(working_dir)
  
  #set up filtering log file:
  system(paste0("if [ ! -d ",paste0(working_dir,"/",project_name)," ]; then mkdir ",paste0(working_dir,"/",project_name),
                "; else echo \"\nWarning: Project directory already exists.\"; fi"))
  
  system(paste0("if [ ! -d ",paste0(working_dir,"/",project_name,"/logs")," ]; then mkdir ",paste0(working_dir,"/",project_name, "/logs"),
                "; else echo \"\nWarning: Project directory already exists.\"; fi"))
  
  system(paste0("if [ ! -d ",paste0(working_dir,"/",project_name,"/output")," ]; then mkdir ",paste0(working_dir,"/",project_name,"/output"),
                "; fi"))
                
  write.log(paste(project_name, Sys.time()), paste(working_dir,project_name,"logs/filter.log",sep="/"))
  
  m <- paste(snp_count(vcf_file, working_dir), "SNPs before filtering.")
  message(m); write.log(paste("\n",m,sep=""), paste(working_dir,project_name,"logs/filter.log",sep="/"))
  
  # if (output){
  #    m <- system(paste(paste(working_dir,"scripts/filtering/ErrorCount.sh",sep="/"), paste(working_dir, vcf_file, sep="/"), sep=" "), intern=TRUE)
  #    for (i in m){
  #      message(i)# ; write.log(pste("\n",i,sep=""), paste0(paste0(working_dir,project_name),"logs/filter.log",sep="/"))
  #    }
  #    write.log(m, paste(paste(working_dir,project_name, sep="/"),"logs/filter.log",sep="/"))
  
  #copies vcf file
  filtered_filename <- paste0("filtered_",vcf_file)
  #sink("/dev/null"); file.copy(from=vcf_file, to=filtered_filename, overwrite = TRUE); sink()
  
  fdf <- data.frame(Filter=as.character("Prefilter"), SNPs=as.numeric(snp_count(vcf_file, working_dir)))
  return(fdf)
  }
#}

#counts SNPs in given vcf file
snp_count <- function(vcf, wd) {
  return(system(paste("grep -v '#'", paste(wd, as.character(vcf), sep="/"), "| wc -l", sep=" "), intern=T))
}

#calculates differences in # SNPs between two vcf files
snp_diff <- function(vcf1, wd, vcf2) {
  return(paste(
    snp_count(vcf2, wd),
    " SNPs kept. ",
    as.numeric(snp_count(vcf1, wd)) - as.numeric(snp_count(vcf2, wd)),
    " (",
    round(100 - (as.numeric(snp_count(vcf2, wd)) / as.numeric(snp_count(vcf1, wd))) * 100, 2),
    "%) filtered. Cumulative filtered: ",
    as.numeric(snp_count(vcf_file, wd)) - as.numeric(snp_count(vcf2, wd)),
    " (",
    round(100 - (as.numeric(snp_count(vcf2, wd)) / as.numeric(snp_count(vcf_file, wd))) * 100, 2),
    "%)"
  , sep=""))
}

# max.missing <- 0.75
# min.allele.count <- 2
# min.read.quality <- 20 
# min.depth <- 3
# min.mean.depth <- 5
# min.mapping.quality <- c(39,39)
# mapping.ratio <- c(0.75,1.25)
# read.balance <- c(0,0)
# quality.depth.ratio <- 0.25
# min.number.pools <- 2
# threads <- 10
# working_dir <- "/10tb_leviathan/evan/multispp/porites/variantcalling/vcf_test_filtering"
# project_name <- "filter_test"
# vcf_file <- "porites_TotalRawSNPs.vcf"
# vcflib_PATH <- "/10tb_leviathan/evan/assesspool_dev/scripts/vcflib/bin/vcffilter"

loss_per_filter <- function(max.missing,
                            min.allele.count,
                            min.mean.depth,
                            max.mean.depth,
                            hwe.cutoff,
                            min.mapping.quality,
                            mapping.ratio,
                            read.balance,
                            quality.depth.ratio,
                            mispaired.reads,
                            min.number.pools,
                            min.depth.all.pools,
                            max.allele.length,
                            min.qual,
                            variant.type,
                            alt.obs,
                            threads,
                            working_dir,
                            project_name,
                            vcf_file,
                            vcflib_PATH)
      {
commands <- paste0("#!/bin/bash
vcf_in=\"",paste(working_dir,vcf_file,sep="/"), "\"
NUMPROC=",threads,"
vcflib_PATH=\"", paste(vcflib_PATH), "\"
#vcflib_PATH=\"/10tb_leviathan/evan/assesspool_dev/scripts/vcflib/bin/vcffilter\" 
proj=\"", paste(project_name), "\"
wd=\"", paste(working_dir), "\"          
if [ ! -d $wd ]; then mkdir $wd; fi\ 

echo \"filter,arg,SNPs\" >  ${wd}/${proj}/${proj}_filtering_results.csv
grep -v '#' $vcf_in | wc -l | awk '{print \"0_prefilter,NA,\"$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv \n",

####VCFTOOLS set up#########
if(!is.na(max.missing)){
  paste0(
    "VT1=\"--max-missing ", max.missing, "\" 
     echo \"vcftools --vcf $vcf_in $VT1 --recode --recode-INFO-all --stdout | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VT1}\\\" '{print \\\"1_max.missing,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" > ${wd}/${proj}/${proj}_filtering_send.sh \n")
  },

if(!is.na(min.allele.count)){
  paste0(
    "VT2=\"--mac ", min.allele.count, "\" 
     echo \"vcftools --vcf $vcf_in $VT2 --recode --recode-INFO-all --stdout | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VT2}\\\" '{print \\\"2_minor.allele.count,\\\" FLAG \\\",\\\" \\$1}' >>${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},

if(!is.na(min.mean.depth)){
  paste0(
    "VT3=\"--min-meanDP ", min.mean.depth, "\" 
     echo \"vcftools --vcf $vcf_in $VT3 --recode --recode-INFO-all --stdout | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VT3}\\\" '{print \\\"3_min.mean.depth,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},

if(!is.na(max.mean.depth)){
  paste0(
    "VT4=\"--max-meanDP ", max.mean.depth, "\" 
     echo \"vcftools --vcf $vcf_in $VT4 --recode --recode-INFO-all --stdout | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VT4}\\\" '{print \\\"4_max.mean.depth,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},

if(!is.na(hwe.cutoff)){
  paste0(
    "VT5=\"--hwe ", hwe.cutoff, "\" 
     echo \"vcftools --vcf $vcf_in $VT5 --recode --recode-INFO-all --stdout | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VT5}\\\" '{print \\\"5_hwe,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},

#####VCFLIB set up##########
if(!is.na(min.mapping.quality)){
  paste0(
    "VL1=(\"MQM > ", min.mapping.quality[1], " & MQMR > ", min.mapping.quality[2], "\") 
    echo \"${vcflib_PATH}/vcffilter -s -f \\\"${VL1[@]}\\\" $vcf_in | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VL1[@]}\\\" '{print \\\"6_min.mapping.quality,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},

if(!is.na(mapping.ratio)){
  paste0(
    "VL2=(\"MQM / MQMR > ", mapping.ratio[1], " & MQM / MQMR < ", mapping.ratio[2], "\") 
    echo \"${vcflib_PATH}/vcffilter -f \\\"${VL2[@]}\\\" $vcf_in | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VL2[@]}\\\" '{print \\\"7_mapping.ratio,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},
   
if(!is.na(read.balance)){
  paste0( 
    "VL3=(\"RPR > ", read.balance[1], " & RPL > ", read.balance[2], "\") 
    echo \"${vcflib_PATH}/vcffilter -f \\\"${VL3[@]}\\\" $vcf_in | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VL3[@]}\\\" '{print \\\"8_read.balance,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},

if(!is.na(quality.depth.ratio)){
  paste0(
    "VL4=(\"QUAL / DP > ", quality.depth.ratio, "\") 
    echo \"${vcflib_PATH}/vcffilter -f \\\"${VL4[@]}\\\" $vcf_in | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VL4[@]}\\\" '{print \\\"9_quality.depth.ratio,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},

if(!is.na(mispaired.reads)){
  paste0(
    "VL5=(\"PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25\")
    echo \"${vcflib_PATH}/vcffilter -f \\\"${VL5[@]}\\\" $vcf_in | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VL5[@]}\\\" '{print \\\"10_mispaired.reads,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},

if(!is.na(min.number.pools)){
  paste0(
    "VL6=(\"NS > ", min.number.pools-1, "\") 
    echo \"${vcflib_PATH}/vcffilter -f \\\"${VL6[@]}\\\" $vcf_in | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VL6[@]}\\\" '{print \\\"11_min.number.pools,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},

if(!is.na(min.depth.all.pools)){
  paste0(
    "VL7=(\"DP > ", min.depth.all.pools, "\") 
    echo \"${vcflib_PATH}/vcffilter -f \\\"${VL7[@]}\\\" $vcf_in | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VL7[@]}\\\" '{print \\\"12_min.depth.all.pools,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},

if(!is.na(max.allele.length)){
  paste0(
    "VL8=(\"LEN < ", max.allele.length+1, "\") 
    echo \"${vcflib_PATH}/vcffilter -f \\\"${VL8[@]}\\\" $vcf_in | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VL8[@]}\\\" '{print \\\"13_max.allele.length,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
  
},

if(!is.na(min.qual)){
  paste0(
    "VL9=(\"QUAL > ", min.qual, "\") 
    echo \"${vcflib_PATH}/vcffilter -f \\\"${VL9[@]}\\\" $vcf_in | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VL9[@]}\\\" '{print \\\"14_min.QUAL,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
  
},

if(!is.na(variant.type)){
  paste0(
    "VL10=(\"TYPE = ", variant.type, "\")
    echo \"${vcflib_PATH}/vcffilter -f \\\"${VL10[@]}\\\" $vcf_in | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VL10[@]}\\\" '{print \\\"15_variant.type,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
  
},
if(!is.na(alt.obs)){
  paste0(
    "VL11=(\"AO > ", alt.obs, "\")
echo \"${vcflib_PATH}/vcffilter -f \\\"${VL11[@]}\\\" $vcf_in | grep -v '#' - | wc -l | awk -v FLAG=\\\"${VL11[@]}\\\" '{print \\\"16_alt.obs,\\\" FLAG \\\",\\\" \\$1}' >> ${wd}/${proj}/${proj}_filtering_results.csv\" >> ${wd}/${proj}/${proj}_filtering_send.sh \n")
},

"parallel --jobs $NUMPROC :::: ${wd}/${proj}/${proj}_filtering_send.sh")[1]

write(commands, file = paste0(working_dir,"/",project_name,"/",project_name,"_filter_wrapper.sh"))
system(paste("bash ", paste0(working_dir,"/",project_name,"/",project_name,"_filter_wrapper.sh")))
}

execute_filtering <- function(max.missing,
                              min.allele.count,
                              min.mean.depth,
                              max.mean.depth,
                              hwe.cutoff,
                              min.mapping.quality,
                              mapping.ratio,
                              read.balance,
                              quality.depth.ratio,
                              mispaired.reads,
                              min.number.pools,
                              min.depth.all.pools,
                              max.allele.length,
                              min.qual,
                              variant.type,
                              alt.obs,
                              working_dir,
                              project_name,
                              vcf_file,
                              vcflib_PATH)
      {
commands <- paste0("#!/bin/bash
vcf_in=\"",paste(working_dir, vcf_file, sep="/"), "\"
vcf_out=\"",paste0(working_dir,"/filtered_",vcf_file), "\"

vcflib_PATH=\"", paste(vcflib_PATH), "\"
#vcflib_PATH=\"/10tb_leviathan/evan/assesspool_dev/scripts/vcflib/bin/vcffilter\"
proj=\"", paste(project_name), "\" 
wd=\"", paste(working_dir), "\"
if [ ! -d $wd ]; then mkdir $wd; fi \n",
                   
####VCFTOOLS set up#########
if(!is.na(max.missing)){
  paste0("VT1=\"--max-missing ", max.missing, "\" \n")},

if(!is.na(min.allele.count)){
  paste0("VT2=\"--mac ", min.allele.count, "\" \n")},

if(!is.na(min.mean.depth)){
  paste0("VT3=\"--min-meanDP ", min.mean.depth, "\" \n")},

if(!is.na(max.mean.depth)){
  paste0("VT4=\"--max-meanDP ", max.mean.depth, "\" \n")},

if(!is.na(hwe.cutoff)){
  paste0("VT5=\"--hwe ", hwe.cutoff, "\" \n")},

#####VCFLIB set up##########
if(!is.na(min.mapping.quality)){
  paste0("VL1=(\"MQM > ", min.mapping.quality[1], " & MQMR > ", min.mapping.quality[2], "\") \n")},

if(!is.na(mapping.ratio)){
  paste0("VL2=(\"MQM / MQMR > ", mapping.ratio[1], " & MQM / MQMR < ", mapping.ratio[2], "\") \n")},

if(!is.na(read.balance)){
  paste0("VL3=(\"RPR > ", read.balance[1], " & RPL > ", read.balance[2], "\") \n")},

if(!is.na(quality.depth.ratio)){
  paste0("VL4=(\"QUAL / DP > ", quality.depth.ratio, "\") \n")},

if(!is.na(mispaired.reads)){
  paste0("VL5=(\"PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 \") \n")},

if(!is.na(min.number.pools)){
  paste0("VL6=(\"NS > ", min.number.pools-1, "\") \n")},

if(!is.na(min.depth.all.pools)){
  paste0("VL7=(\"DP > ", min.depth.all.pools, "\") \n")},

if(!is.na(max.allele.length)){
  paste0("VL8=(\"LEN < ", max.allele.length+1, "\") \n")},

if(!is.na(min.qual)){
  paste0("VL9=(\"QUAL > ", min.qual, "\") \n")},

if(!is.na(variant.type)){
  paste0("VL10=(\"TYPE = ", variant.type, "\") \n")},

if(!is.na(alt.obs)){
  paste0("VL11=(\"AO > ", alt.obs, "\") \n")},                   

"${vcflib_PATH}/vcffilter -f \"",
      if(!is.na(min.mapping.quality)){paste0("${VL1} & ")},
      if(!is.na(mapping.ratio)){paste0("${VL2} & ")},
      if(!is.na(read.balance)){paste0("${VL3} & ")},
      if(!is.na(quality.depth.ratio)){paste0("${VL4} & ")},
      if(!is.na(mispaired.reads)){paste0("${VL5} & ")},
      if(!is.na(min.number.pools)){paste0("${VL6} & ")},
      if(!is.na(min.depth.all.pools)){paste0("${VL7} & ")},
      if(!is.na(max.allele.length)){paste0("${VL8} & ")},
      if(!is.na(min.qual)){paste0("${VL9} & ")},
      if(!is.na(variant.type)){paste0("${VL10} & ")},
      if(!is.na(alt.obs)){paste0("${VL11}")},
       "\" $vcf_in | \\
vcftools --vcf - ",
      if(!is.na(max.missing)){paste=("$VT1 ")},
      if(!is.na(min.allele.count)){paste=("$VT2 ")},
      if(!is.na(min.mean.depth)){paste=("$VT3 ")},
      if(!is.na(max.mean.depth)){paste=("$VT4 ")},
      if(!is.na(hwe.cutoff)){paste=("$VT5 ")},
"--recode --recode-INFO-all --stdout > ${vcf_out}

echo \"SNPs passing all filters:\"
grep -v '#' $vcf_out | wc -l")

write(commands, file = paste(working_dir,project_name,"run_filter_wrapper.sh",sep="/"))
system(paste("bash ", paste(working_dir,project_name,"run_filter_wrapper.sh",sep="/")))

#Save # SNPs to filtering log
numSNPs <- system(paste("grep -v '#'", paste0("filtered_",vcf_file), "| wc -l", sep=" "), intern=T)
tmp_fdf <- data.frame(filter=as.character(paste("all_filters",sep="_")),arg=as.character(paste("NA")), SNPs=as.numeric(numSNPs))
post_filter_df <<- dplyr::full_join(filter_df, tmp_fdf)
return(fdf)        
}




### Thin vcf to avoid linked loci
filter_thinning <- function(working_dir, project_name, vcf_file, x, filter_df){
  setwd(working_dir)
  #system(paste("sh", paste0(working_dir,"/scripts/filtering/vcf_par.sh"), "-v", paste0("filtered_",vcf_file), "-t", threads, "-f \"--thin",x, "\""))
  system(paste("vcftools --vcf", paste0("filtered_",vcf_file), "--thin", x, "--recode --recode-INFO-all --stdout > ", paste0("final_filter_",vcf_file)))
  
  m  <- paste("After thinning to min ",x,"bp between loci:", sep="")
  message(""); message(m); write.log(paste("\n\n\n\n",m,sep=""), paste(working_dir,project_name,"logs/filter.log",sep="/"))
  
  m <- snp_diff(paste0("filtered_",vcf_file), working_dir, paste0("final_filter_",vcf_file))
  message(m); write.log(paste(m,"\n\n",sep=""), paste(working_dir,project_name,"logs/filter.log",sep="/"))
  
  #Save # SNPs to filtering log
  numSNPs <- system(paste("grep -v '#'", paste0("final_filter_",vcf_file), "| wc -l", sep=" "), intern=T)
  tmp_fdf <- data.frame(filter=as.character(paste("thinned",x,sep="_")),arg=as.character(paste("--thin", x)), SNPs=as.numeric(numSNPs))
  fdf <- dplyr::full_join(filter_df, tmp_fdf)
  return(fdf)
}


#run ErrorCount.sh from ddocent to estimate potential genotyping errors due to low read depth
run_ddocentErrorCount <- function(working_dir, vcf){
  message("#########\nPrefiltering:\n#########")
  m <- system(paste(paste(working_dir,"/scripts/filtering/ErrorCount.sh", sep=""), paste0(vcf), sep=" "), intern=TRUE)
  for (i in m){
    message(i)
  }
  write.log(m, paste(working_dir,project_name,"logs/filter.log",sep="/"))
  
  message("#########\nAfter filtering:\n#########")
  m <- system(paste(paste(working_dir,"/scripts/filtering/ErrorCount.sh", sep=""), paste0("final_filter_",vcf), sep=" "), intern=TRUE)
  for (i in m){
    message(i)
  }
  write.log(m, paste(working_dir,project_name,"logs/filter.log",sep="/"))
}

#snp_histogram
snp_hist <- function(depth_table) {
  ggplotly(read.table(depth_table) %>% dplyr::rename(.,depth=V1) %>%
             ggplot() + geom_histogram(aes(depth/max.missing), binwidth=5) + 
             xlab("Mean Depth") + ylab("# SNPs") + theme_bw() + ggtitle("Choose cutoff"))
}

count_SNP_loss <- function(max.missing,
                           min.allele.count,
                           min.mean.depth,
                           max.mean.depth,
                           hwe.cutoff,
                           min.mapping.quality,
                           mapping.ratio,
                           read.balance,
                           quality.depth.ratio,
                           mispaired.reads,
                           min.number.pools,
                           min.depth.all.pools,
                           max.allele.length,
                           min.qual,
                           variant.type,
                           alt.obs,
                           threads,
                           working_dir,
                           project_name,
                           vcf_file,
                           vcflib_PATH)
{
  filter_df <<- vcf_init(vcf_file, working_dir, show.filter.output)
  system(paste0('echo "vcf initialization done \\m/"'))

  loss_per_filter(max.missing=max.missing,
                  min.allele.count=min.allele.count,
                  min.mean.depth=min.mean.depth,
                  max.mean.depth=max.mean.depth,
                  hwe.cutoff=hwe.cutoff,
                  min.mapping.quality=min.mapping.quality,
                  mapping.ratio=mapping.ratio,
                  read.balance=read.balance,
                  quality.depth.ratio=quality.depth.ratio,
                  mispaired.reads=mispaired.reads,
                  min.number.pools=min.number.pools,
                  min.depth.all.pools=min.depth.all.pools,
                  max.allele.length=max.allele.length,
                  min.qual=min.qual,
                  variant.type=variant.type,
                  alt.obs=alt.obs,
                  threads=threads,
                  working_dir=working_dir,
                  project_name=project_name,
                  vcf_file=vcf_file,
                  vcflib_PATH=vcflib_PATH)
}


run_selected_filters <- function(max.missing,
                                 min.allele.count,
                                 min.mean.depth,
                                 max.mean.depth,
                                 hwe.cutoff,
                                 min.mapping.quality,
                                 mapping.ratio,
                                 read.balance,
                                 quality.depth.ratio,
                                 mispaired.reads,
                                 min.number.pools,
                                 min.depth.all.pools,
                                 max.allele.length,
                                 min.qual,
                                 variant.type,
                                 alt.obs,
                                 working_dir,
                                 project_name,
                                 vcf_file,
                                 vcflib_PATH)
{
  
   if(filter_SNPS==T){    
    execute_filtering(max.missing=max.missing,
                      min.allele.count=min.allele.count,
                      min.mean.depth=min.mean.depth,
                      max.mean.depth=max.mean.depth,
                      hwe.cutoff=hwe.cutoff,
                      min.mapping.quality=min.mapping.quality,
                      mapping.ratio=mapping.ratio,
                      read.balance=read.balance,
                      quality.depth.ratio=quality.depth.ratio,
                      mispaired.reads=mispaired.reads,
                      min.number.pools=min.number.pools,
                      min.depth.all.pools=min.depth.all.pools,
                      max.allele.length=max.allele.length,
                      min.qual=min.qual,
                      variant.type=variant.type,
                      alt.obs=alt.obs,
                      working_dir=working_dir,
                      project_name=project_name,
                      vcf_file=vcf_file,
                      vcflib_PATH=vcflib_PATH)
}
}
