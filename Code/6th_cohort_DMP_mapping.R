library(readxl)
library(data.table)
library(dplyr)

manifest.6th <- read_xlsx('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/6th_cohort/Proj_05500_GB_sample_manifest18July2019.xlsx') %>%
  # total manifest, so have to filter and pick out EDD ret samples. 
  filter(STUDY_ID == 'EDD-ret') %>%
  transmute(`CMO sample ID` = gsub('_IGO.*.$','',CMO_SAMPLE_ID),Sex = SEX,
            `CMO patient ID` = CMO_PATIENT_ID, `Investigator sample ID` = INVESTIGATOR_SAMPLE_ID,
            `Patient ID` = gsub('-','_',gsub('_bc.*.$|-cf.*.$|_cf.*.$','',INVESTIGATOR_SAMPLE_ID))) %>% data.table()

DMP.mapping.3rd <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_RET_ID/DMP_RET_ID_processed.txt') %>%
  transmute(MRN = gsub('No Impact',NA,MRN),`Patient ID`) %>% unique() %>% data.table()

manifest.6th <- merge(manifest.6th,DMP.mapping.3rd,by = 'Patient ID',all.x = T) %>% select(MRN,everything()) %>% select(-`Patient ID`)

write.csv(manifest.6th,'/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_mapping/DMP_mapping_6.csv',row.names = F)
