library(readxl)
library(data.table)
library(dplyr)

manifest.7th <- read_xlsx('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/7th_cohort/5500_GD_manifest.xlsx') %>%
  # total manifest, so have to filter and pick out EDD ret samples. 
  filter(STUDY_ID == 'EDD_ret') %>%
  transmute(`CMO sample ID` = gsub('_IGO.*.$','',CMO_SAMPLE_ID),Sex = SEX,
            `CMO patient ID` = CMO_PATIENT_ID, `Investigator sample ID` = INVESTIGATOR_SAMPLE_ID,
            `Patient ID` = gsub('-','_',gsub('_bc.*.$|-cf.*.$|_cf.*.$','',INVESTIGATOR_SAMPLE_ID))) %>% data.table()

DMP.mapping.3rd <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_RET_ID/DMP_RET_ID_processed.txt') %>%
  transmute(MRN = gsub('No Impact',NA,MRN),`Patient ID`) %>% unique() %>% data.table()

manifest.7th <- merge(manifest.7th,DMP.mapping.3rd,by = 'Patient ID',all.x = T) %>% select(MRN,everything()) %>% select(-`Patient ID`)

write.csv(manifest.7th,'/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_mapping/DMP_mapping_7.csv',row.names = F)
