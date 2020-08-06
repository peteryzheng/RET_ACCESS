library(data.table)
library(tidyr)
library(stringr)
library(dplyr)
source('/ifs/work/bergerm1/zhengy1/RET_all/Code/genotype_source_code.R')
source('/ifs/work/bergerm1/zhengy1/RET_all/Code/table_to_maf.R')

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_121519.csv')
# results.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_',format(Sys.time(),'%m%d%y'))
results.dir <- "/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_021320/"
CH.calls = fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/signedout_CH.txt')


# criteria <- 'permissive'
criteria <- 'stringent'

if(criteria == 'permissive'){
  hotspot.support <- 1
  non.hotspot.support <- 5
}else{
  hotspot.support <- 3
  non.hotspot.support <- 5
}

dir.create(paste0(results.dir,'/results_',criteria))

# DMP stuff ---------------------------------------------------------------
DMP.key <- fread('/ifs/dmprequest/12-245/key.txt')

# Pooled normals ----------------------------------------------------------
duplexsupport <- function(x) {
  print(x)
  length(which(x >= 2))
}
pooled.normal.mafs <-
  fread(paste0(results.dir,'/pooled/all_all_unique_fillout_maf2maf.maf')) %>% mutate(Tumor_Sample_Barcode = paste0(Tumor_Sample_Barcode,'___pooled')) %>%
  select(Hugo_Symbol,Tumor_Sample_Barcode,Chromosome,Start_Position,End_Position,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2,t_alt_count_fragment) %>% 
  group_by(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Reference_Allele,Tumor_Seq_Allele2) %>% 
  summarise(duplex_support_num = length(which(t_alt_count_fragment >= 2))) %>% filter(duplex_support_num >= 2) %>% data.table()

# for each patient produce the correct results ----------------------------
x <- 'C-0EU9LX'
x <- unique(master.ref$`CMO patient ID`)[3]
all.sample.maf <- do.call(rbind,lapply(unique(master.ref$`CMO patient ID`),function(x){
  print(x)
  # Inputs and sanity checks ------------------------------------------------
  fillouts.filenames <- list.files(paste0(results.dir,'/',x,'/fillout'),full.names = T)  
  fillouts.dt <- do.call(rbind,lapply(fillouts.filenames,function(y){
    maf.file <- fread(y) 
    print(unique(maf.file$Tumor_Sample_Barcode))
    # fragment counts replacing actual allele counts
    if(!grepl('normal_DMP_fillout_maf2maf|Tumor_fillout_maf2maf',y)){
      maf.file <- maf.file[,c('t_alt_count','t_ref_count','t_depth') := list(t_alt_count_fragment,t_ref_count_fragment,t_total_count_fragment)][,!c('t_alt_count_fragment','t_ref_count_fragment','t_total_count_fragment'),with = F] %>%
        # sample name cleanup
        mutate(Tumor_Sample_Barcode = gsub('_cl.*.FX-|_S.._..._cl.*.FX-|_cl.*.FX|_S.._..._cl.*.FX',
                                           '___',Tumor_Sample_Barcode)) %>% 
        mutate(Tumor_Sample_Barcode = ifelse(grepl('-N[0-9][0-9][0-9]-d',Tumor_Sample_Barcode),paste0(Tumor_Sample_Barcode,'unfilterednormal'),Tumor_Sample_Barcode)) %>% data.table()
    }else{
      maf.file <- maf.file %>% merge(DMP.key[,.(V1,V2)],by.x = 'Tumor_Sample_Barcode',by.y = 'V2',all.x = T) %>%
        mutate(Tumor_Sample_Barcode = gsub(' ','',ifelse(grepl('-T',V1),paste0(V1,'___Tumor'),paste(V1,'___normal_DMP')))) %>% 
        select(-c(V1,t_alt_count_fragment,t_ref_count_fragment,t_total_count_fragment)) %>% data.table()
    }
    return(maf.file)
  }))
  # # Testing (commented out because testing finished)  -----------------------
  #   if(length(unique(table(fillouts.dt$Tumor_Sample_Barcode))) == 1){
  #     print(paste0('All mafs for patient ',x,' are the same!'))
  #     return(TRUE)
  #   }else{
  #     print(paste0('All mafs for patient ',x,' are not the same!'))
  #     return(FALSE)
  #   }
  # merging and melting -----------------------------------------------------
  sample.sheet <- fread(paste0(results.dir,'/',x,'/',x,'_sample_sheet.tsv')) %>% rowwise() %>%
    mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                                ifelse(Sample_Type == 'normal','unfilterednormal',
                                       ifelse(Sample_Type == 'plasma_simplex','simplex',Sample_Type)))) %>%
    mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
  hotspot.maf <- fread(paste0(results.dir,'/',x,'/',x,'_all_unique_calls_hotspots.maf')) %>% rowwise() %>%
    transmute(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,
              # HGVSp_Short,
              Reference_Allele,Tumor_Seq_Allele2,Hotspot = ifelse(hotspot_whitelist,'Hotspot',NA)) %>% data.table()
  dmp.maf <- fread(paste0(results.dir,'/',x,'/',x,'_impact_calls.maf')) %>%
    select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,
           # HGVSp_Short,
           Reference_Allele,Tumor_Seq_Allele2) %>% mutate(DMP = 'Signed out') %>% unique() %>% data.table()
  
  fillouts.dt <-
    fillouts.dt %>% mutate(t_var_freq = paste0(t_alt_count,'/',t_depth,'(',round(t_alt_count/t_depth,4),')')) %>% 
    select(Hugo_Symbol,Tumor_Sample_Barcode,Chromosome,Start_Position,End_Position,Variant_Classification,HGVSp_Short,Reference_Allele,
           Tumor_Seq_Allele2,t_var_freq,ExAC_AF) %>% spread(Tumor_Sample_Barcode,t_var_freq) %>%
    # hotspot information 
    merge(hotspot.maf,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification',
                             # 'HGVSp_Short',
                             'Reference_Allele','Tumor_Seq_Allele2'),all.x = T) %>%
    # Identifying signed out calls
    merge(dmp.maf,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2'),all.x = T) %>%
    # pooled normal for systemic artifacts
    merge(pooled.normal.mafs,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','Reference_Allele','Tumor_Seq_Allele2'),all.x = T) %>% 
    filter(is.na(duplex_support_num) | !is.na(DMP)) %>% data.table()
  # Interesting cases where DMP signed out calls are artifacets
  if(any(!is.na(fillouts.dt$DMP) & !is.na(fillouts.dt$duplex_support_num))){
    print(paste0('Look at ',x,' for DMP signed out plasma artifacts...'))
  }
  # germline filtering for matched and unmatched ----------------------------
  plasma.samples <- sample.sheet[Sample_Type %in% c('duplex')]$column.names
  normal.samples <- sample.sheet[Sample_Type %in% c('unfilterednormal','normal_DMP')]$column.names
  fillouts.dt[,c(
    paste0(plasma.samples,'.called')
    # paste0(gsub('duplex','simplex',plasma.samples),'.called')
    
  ) := 'Not Called']
  # preliminary calling
  # tmp.col.name <- plasma.samples[1]
  lapply(plasma.samples,function(tmp.col.name){
    # genotyping (signed out stuff)
    fillouts.dt[(as.numeric(gsub("/.*.$",'',get(tmp.col.name))) >= 1 | as.numeric(gsub("/.*.$",'',get(paste0(gsub('duplex','simplex',tmp.col.name))))) > 1) & DMP == 'Signed out',
                eval(paste0(tmp.col.name,'.called')) := 'Genotyped']
    # c(eval(paste0(tmp.col.name,'.called')),eval(paste0(gsub('duplex','simplex',tmp.col.name),'.called'))) := list('Called','Called')]
    # hotspot reads
    if(criteria == 'stringent'){
      fillouts.dt[as.numeric(gsub("/.*.$",'',get(tmp.col.name))) >= hotspot.support & Hotspot == 'Hotspot',
                  eval(paste0(tmp.col.name,'.called')) := 'Called']
    }else{
      fillouts.dt[(as.numeric(gsub("/.*.$",'',get(tmp.col.name))) > hotspot.support | 
                     (as.numeric(gsub("/.*.$",'',get(tmp.col.name))) >= hotspot.support &
                      as.numeric(gsub("/.*.$",'',get(gsub('duplex','simplex',tmp.col.name)))) >= hotspot.support)) &
                    Hotspot == 'Hotspot',
                  eval(paste0(tmp.col.name,'.called')) := 'Called']
    }
    # c(eval(paste0(tmp.col.name,'.called')),eval(paste0(gsub('duplex','simplex',tmp.col.name),'.called'))) := list('Called','Called')]
    # non hotspot reads
    fillouts.dt[as.numeric(gsub("/.*.$",'',get(tmp.col.name))) >= non.hotspot.support & is.na(Hotspot),
                eval(paste0(tmp.col.name,'.called')) := 'Called']
    # c(eval(paste0(tmp.col.name,'.called')),eval(paste0(gsub('duplex','simplex',tmp.col.name),'.called'))) := list('Called','Called')]
    print(table(fillouts.dt[,get(paste0(tmp.col.name,'.called'))]))
  })
  if(all(!c('unfilterednormal','normal_DMP') %in% sample.sheet$Sample_Type)){
    # tmp.col.name <- plasma.samples[1]
    lapply(plasma.samples,function(tmp.col.name){
      # germline filtering
      fillouts.dt[(as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name),"\\(.*.\\)"))) >= 0.3 | 
                     ExAC_AF >= 0.0001) & Hugo_Symbol != 'RET',
                  eval(paste0(tmp.col.name,'.called')) := 'Not Called']
      # not covered
      fillouts.dt[get(tmp.col.name)  == '0/0(NaN)',eval(paste0(tmp.col.name,'.called')) := 'Not Covered']
    })
  }else{
    lapply(plasma.samples,function(tmp.col.name){
      lapply(normal.samples,function(tmp.col.name.normal){
        # duplex tvar/nvar > 5
        fillouts.dt[((as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name),"\\(.*.\\)")))/as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name.normal),"\\(.*.\\)"))) < 2) |
                       # if duplex have no reads, use simplex tvar
                       (as.numeric(gsub("\\(|\\)",'',str_extract(get(gsub('duplex','simplex',tmp.col.name)),"\\(.*.\\)")))/as.numeric(gsub("\\(|\\)",'',str_extract(get(tmp.col.name.normal),"\\(.*.\\)"))) < 2 &
                          as.numeric(gsub("/.*.$",'',get(tmp.col.name)))  == 0)) &
                      (Hugo_Symbol != 'RET' | DMP != 'Signed out'),
                    eval(paste0(tmp.col.name,'.called')) := 'Not Called']
        # annotating not covered calls
        fillouts.dt[get(tmp.col.name)  == '0/0(NaN)',eval(paste0(tmp.col.name,'.called')) := 'Not Covered']
      })
    })
  }
  
  # final processing --------------------------------------------------------
  # Save only the useful column
  fillouts.dt <-  fillouts.dt[DMP == 'Signed out' | fillouts.dt[,apply(.SD,1,function(x){any(x == 'Called')})]]  
  # combining duplex and simplex counts
  lapply(plasma.samples,function(tmp.col.name){
    # hotspot reads
    fillouts.dt[,eval(gsub('duplex','total',tmp.col.name)) := paste0(
      as.numeric(gsub("/.*.$",'',get(tmp.col.name)))+as.numeric(gsub("/.*.$",'',get(gsub('duplex','simplex',tmp.col.name)))),'/',
      as.numeric(gsub("^.*./|\\(.*.$",'',get(tmp.col.name)))+as.numeric(gsub("^.*./|\\(.*.$",'',get(gsub('duplex','simplex',tmp.col.name)))),'(',
      round((as.numeric(gsub("/.*.$",'',get(tmp.col.name)))+as.numeric(gsub("/.*.$",'',get(gsub('duplex','simplex',tmp.col.name)))))/
              (as.numeric(gsub("^.*./|\\(.*.$",'',get(tmp.col.name)))+as.numeric(gsub("^.*./|\\(.*.$",'',get(gsub('duplex','simplex',tmp.col.name))))),4),')'
    )]
    fillouts.dt[,c(eval(gsub('duplex','simplex',tmp.col.name)),eval(tmp.col.name)):= list(NULL,NULL)]
  })
  fillouts.dt <- fillouts.dt[,order(colnames(fillouts.dt)),with = F] %>%
    # filter for artifacts
    mutate(call_confidence = ifelse(
      (Hugo_Symbol == 'TERT' & Hotspot != 'Hotspot') | (Hugo_Symbol == 'ERBB2' & grepl('[A-Z]90[0-9][A-Z]',HGVSp_Short)) | 
        (Hugo_Symbol == 'APC' & grepl('142',HGVSp_Short)) | (Hugo_Symbol == 'NF1' & grepl('[A-Z]106[0-9][A-Z]',HGVSp_Short)), 'Low',''
    )) %>% filter(call_confidence != 'Low') %>%
    merge(CH.calls[,.(Hugo_Symbol = Gene,Chromosome = Chrom,Start_Position = Start,Reference_Allele = Ref,Tumor_Seq_Allele2 = Alt,HGVSp_Short = AAchange,Variant_Classification = VariantClass,CH = 'Yes')],
          by = c('Hugo_Symbol','Chromosome','Start_Position','Variant_Classification','HGVSp_Short','Reference_Allele','Tumor_Seq_Allele2'),
          all.x = T) %>%
    select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2,
           ExAC_AF,Hotspot,DMP,CH,duplex_support_num,sort(everything())) %>% unique()
  # Adding CH flags
  
  
  write.csv(fillouts.dt,paste0(results.dir,'/results_',criteria,'/ret_',
                               gsub('_|-','',str_extract(master.ref[`CMO patient ID` == x]$`Investigator sample ID.plasma`[1],'[0-9]+-|[0-9]+_')),
                               '_SNV_table.csv'),row.names = F)
  if(nrow(fillouts.dt) > 0){
    return_table = table_to_maf(fillouts.dt,sample.sheet)
    return_table$Patient_ID = paste0('ret_',gsub('_|-','',str_extract(master.ref[`CMO patient ID` == x]$`Investigator sample ID.plasma`[1],'[0-9]+-|[0-9]+_')))
    return(return_table)
  }
}))

write.table(all.sample.maf,paste0(results.dir,'/results_',criteria,'/all_sample_maf.maf'),quote = F,sep = '\t',row.names = F)

system(paste0(
  ' /opt/common/CentOS_7-dev/bin/perl /opt/common/CentOS_6-dev/vcf2maf/v1.6.17/maf2maf.pl ',
  ' --input-maf ',results.dir,'/results_',criteria,'/all_sample_maf.maf',
  ' --ncbi-build GRCh37 --tum-vad-col t_alt_count --tum-depth-col t_total_count --custom-enst /opt/common/CentOS_6-dev/vcf2maf/v1.6.17/data/isoform_overrides_at_mskcc ',
  ' --tmp-dir ',tempfile(pattern = "file", tmpdir = paste0(results.dir,'/tmp'), fileext = ""),
  ' --vep-path /opt/common/CentOS_6-dev/vep/v95 --vep-forks 4 --nrm-rad-col n_ref_count --buffer-size 5000 ',
  ' --filter-vcf /opt/common/CentOS_6-dev/vep/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --nrm-vad-col n_alt_count --max-filter-ac 10 ',
  ' --vep-data /opt/common/CentOS_6-dev/vep/cache/ --tum-rad-col t_ref_count --nrm-depth-col n_depth ',
  ' --ref-fasta /ifs/depot/pi/resources/genomes/GRCh37/fasta/b37.fasta ',
  ' --output-maf ',results.dir,'/results_',criteria,'/all_sample_maf2maf.maf',
  ' --retain-cols Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,',
  'Tumor_Sample_UUID,Matched_Norm_Sample_UUID,Caller,Hotspot,DMP,duplex_support_num,call_confidence,Patient_ID --species homo_sapiens  '
))

system(paste0(
  'python /ifs/work/bergerm1/zhengy1/oncokb-annotator/MafAnnotator.py -i ',
  results.dir,'/results_',criteria,'/all_sample_maf2maf.maf','  -o ',
  results.dir,'/results_',criteria,'/all_sample_maf2maf_oncokb.maf',' -b ',
  'fbf69d61-bebb-44cd-9214-d6e867d4bde5'
))  

# part of debugging section in the main lapply function
# if(all(unlist(all.fillout.dim))){
#   print('All dimension of  fillout mafs for each patient looks correct')
# }
