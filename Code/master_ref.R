library(data.table)
library(tidyr)
library(dplyr)

dmp.mapping.data <- do.call(rbind,lapply(list.files('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_mapping/',full.names = T),fread)) %>% rowwise %>%
  mutate(BAM_path = ifelse(length(list.files('/ifs/work/bergerm1/zhengy1/RET_all/all_bams/duplex_bams/',pattern = `CMO sample ID`)) == 0,
                           list.files('/ifs/work/bergerm1/zhengy1/RET_all/all_bams/duplex_bams/',pattern = paste0(`Investigator sample ID`,'.*.bam'),full.names = T),
                           list.files('/ifs/work/bergerm1/zhengy1/RET_all/all_bams/duplex_bams/',pattern = paste0(`CMO sample ID`,'.*.bam'),full.names = T))) %>% rowwise %>%
  mutate(Tumor_Sample_Barcode = ifelse(grepl('DA-ret',`Investigator sample ID`),`Investigator sample ID`,`CMO sample ID`)) %>%
  mutate(Sample_Type = ifelse(grepl('L00.',`CMO sample ID`),'Plasma','Normal')) %>%
  # unite(united_column,BAM_path,`CMO sample ID`, `Investigator sample ID`,Tumor_Sample_Barcode,sep = '___') %>%
  # spread(Sample_Type,united_column) %>% separate(Plasma,paste0(c('BAM_path','CMO sample ID', 'Investigator sample ID','Tumor_Sample_Barcode'),'.plasma'),sep = '__')
  data.table()

dmp.mapping.data <-  merge(dmp.mapping.data[Sample_Type == 'Plasma'],dmp.mapping.data[Sample_Type == 'Normal'],
                           by = c('CMO patient ID','MRN','Sex'),all = T, suffixes = c('.plasma','.normal')) %>% 
  mutate(Paired = ifelse(!is.na(Sample_Type.plasma) & !is.na(Sample_Type.normal),'Paired','Unpaired')) %>% 
  mutate(DMP_ID = MRN,BAM_path.normal = gsub('duplex_bams','unfiltered_bams',gsub('-duplex','',BAM_path.normal)),
         EDD_Patient_ID = paste0('EDD_ret_pt',gsub('EDD_ret_pt|EDD-ret-pt|DA-ret-','',
                                                 str_extract(master.ref$`Investigator sample ID.plasma`,
                                                             'EDD\\-ret\\-pt...|EDD_ret_pt...|DA-ret-...')))) %>%
  select(-MRN) %>% data.table()



write.csv(dmp.mapping.data,
          paste0('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_',format(Sys.time(),'%m%d%y'),'.csv'),
          row.names  = F)

write.table(dmp.mapping.data[,.(tumor_id = Tumor_Sample_Barcode.plasma,normal_id = Tumor_Sample_Barcode.normal)],
            paste0('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/mutation_calling_pairs_',format(Sys.time(),'%m%d%y'),'.tsv'),
            row.names  = F,sep = '\t',quote = F)
a