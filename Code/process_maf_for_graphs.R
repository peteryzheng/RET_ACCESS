library(data.table)
library(tidyverse)
source('/ifs/work/bergerm1/zhengy1/RET_all/Code/collapse_AF.R')

# convert naming to timepoint, get rid of uncovered impact and access calls
process_maf_for_graph = function(tmp.maf){
  print('convert naming to timepoint, get rid of uncovered impact and access calls')
  # tmp.maf = ret.054.calls
  # tumor sample
  tumor.sample = structure(gsub('-','',str_extract(unique(tmp.maf$Tumor_Sample_Barcode[grep('IM[0-9]$',tmp.maf$Tumor_Sample_Barcode)]),'-T..-')),
                           names = as.character(unique(tmp.maf$Tumor_Sample_Barcode[grep('IM[0-9]$',tmp.maf$Tumor_Sample_Barcode)])))
  print(tumor.sample)
  # rest of the samples are plasma
  plasma.sample = setdiff(tmp.maf$Tumor_Sample_Barcode,names(tumor.sample))
  # filter for plasma sample only
  tmp.maf = tmp.maf[Tumor_Sample_Barcode %in% plasma.sample]
  # change samples into timepoint information
  plasma.sample = structure(case_when(
    # some of the DA-ret sample need to be renamed
    grepl('-T0._',plasma.sample) ~ gsub('-|_','',gsub('T','L0',str_extract(plasma.sample,'-T0._'))),
    # otherwise extract L00 something
    TRUE ~ gsub('-','',str_extract(plasma.sample,'-L...-'))
  ),names = plasma.sample)
  print(plasma.sample)
  sample.name.conversion = c(tumor.sample,plasma.sample)
  print(sample.name.conversion)
  # get all not covered calls
  not.covered.df = unique(tmp.maf[call_confidence == 'Not Covered',.N,.(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,
                                                                     HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2)])[N > length(plasma.sample)/2]
  only.covered.tmp.maf = anti_join(tmp.maf,not.covered.df,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification',
                                                                 'HGVSp_Short','Reference_Allele','Tumor_Seq_Allele2')) %>% data.table()
  # only the converted timepoint names
  only.covered.tmp.maf$Tumor_Sample_Barcode = sample.name.conversion[as.character(only.covered.tmp.maf$Tumor_Sample_Barcode)]
  if(any(grepl('__',only.covered.tmp.maf$Hugo_Symbol))){
    fusion.only.covered.tmp.maf = data.table(only.covered.tmp.maf)[grepl('__',Hugo_Symbol)]
    # process original hugo_symbol column (sort two genes by name)
    fusion.only.covered.tmp.maf$Hugo_Symbol = unlist(lapply(fusion.only.covered.tmp.maf$Hugo_Symbol,function(x){paste0(sort(str_split(x,'__')[[1]]),collapse = '-')}))
    fusion.only.covered.tmp.maf$Chromosome = unlist(lapply(fusion.only.covered.tmp.maf$Chromosome,function(x){paste0(sort(str_split(x,'__')[[1]]),collapse = '-')}))
    # collapsing AF for rows of the same events (i.e. reciprocal rearrangement) while perserving the sample level seaparation in AF
    fusion.only.covered.tmp.maf = fusion.only.covered.tmp.maf[,.(Start_Position = Start_Position[1],End_Position = End_Position[1],HGVSp_Short = HGVSp_Short[1],
                                                                 Reference_Allele = Reference_Allele[1], Tumor_Seq_Allele2 = Tumor_Seq_Allele2[1],
                                                                 ExAC_AF = ExAC_AF[1], Hotspot = Hotspot[1], DMP = DMP[1], duplex_support_num = duplex_support_num[1],
                                                                 call_confidence = ifelse(any(call_confidence == 'Called'),'Called','Not Called'),
                                                                 call_info = paste0(call_info,collapse = ' | '),CH = 'No',
                                                                 t_alt_count = sum(t_alt_count,na.rm = T),t_total_count = sum(t_total_count,na.rm = T)),
                                                              .(Hugo_Symbol,Chromosome,Variant_Classification,Tumor_Sample_Barcode)]
    only.covered.tmp.maf = only.covered.tmp.maf[-grep('__',Hugo_Symbol)]
    only.covered.tmp.maf = rbind(only.covered.tmp.maf,fusion.only.covered.tmp.maf)
  }
  only.covered.tmp.maf$t_alt_count = ifelse(is.na(only.covered.tmp.maf$t_alt_count),0,only.covered.tmp.maf$t_alt_count)
  only.covered.tmp.maf$t_total_count = ifelse(is.na(only.covered.tmp.maf$t_total_count),0,only.covered.tmp.maf$t_total_count)
  return(only.covered.tmp.maf)
}