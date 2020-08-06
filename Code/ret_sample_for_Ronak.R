library(data.table)
library(stringr)
library(tidyverse)

source('/ifs/work/bergerm1/zhengy1/RET_all/Code/collapse_AF.R')

results.dir <- "/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_120919/"
condition.dir = 'results_combined_stringent'
results.dir = paste(results.dir,condition.dir,sep = '/')
# ret.fusion.pt = fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/RET_fusion_patients.txt')
ret.alt.pt = fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/RET_access_baseline.txt')

x = "/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_091819//results_combined_permissive/ret_004_table.csv"
x = "/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_120919/results_combined_stringent/ret_006_table.csv"
do.call(rbind,lapply(list.files(results.dir,pattern = 'table.csv',full.names = T),function(x){
  print(x)
  tmp.data = fread(x)[grepl('__RET|RET__',Hugo_Symbol)]
  fusion_expected = FALSE
  if(nrow(tmp.data) > 0){
    return(tmp.data %>% select(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,contains('___total')) %>%
             melt.data.table(id.vars = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification'),
                             variable.name = 'Tumor_Sample_Barcode',value.name = 'VAF',na.rm = F) %>%
             mutate(Patient_ID = gsub('ret_','EDD_ret_pt',str_extract(x,'ret_...')),Hugo_Symbol = gsub('__','-',Hugo_Symbol),
                    Tumor_Sample_Barcode = gsub('___total','',Tumor_Sample_Barcode),VAF = gsub('.*.\\(|\\)','',VAF)) %>%
             separate(Chromosome,into = c('Chr1','Chr2'),sep = '__') %>% data.table() %>%
             setnames(c('Start_Position','End_Position'),c('Pos1','Pos2')))
  }
})) %>% data.table() -> ret.detection.table

ret.detection.table = merge(ret.detection.table,ret.alt.pt[Alt_Type == 'Fusion' & grepl('-',genetic_alterations),.(Study_ID,DMP_ID,genetic_alterations,indicated = 'Y')],
      by.x = c('Hugo_Symbol','Patient_ID'),by.y = c('genetic_alterations','Study_ID'),all = T)

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_030920.csv')[,.(Tumor_Sample_Barcode = Tumor_Sample_Barcode.plasma,Investigator_sample_ID=`Investigator sample ID.plasma`,BAM_path.plasma,Patient_ID = EDD_Patient_ID)]

merge(ret.detection.table,master.ref,by = c('Tumor_Sample_Barcode','Patient_ID'),all = T) %>% 
  filter(!is.na(Tumor_Sample_Barcode)) %>% rowwise() %>%
  mutate(request_id = gsub('/.*.','',
                           gsub('.*.ACCESS-Projects/','',
                                system(paste0('readlink -m ',BAM_path.plasma),intern = T)
                           ))) -> final.table

write.csv(final.table,'/ifs/work/bergerm1/zhengy1/RET_all/For_Ronak/SV_table.csv',row.names = F)