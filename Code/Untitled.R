library(readxl)
library(data.table)
library(dplyr)

manifest.5th <- read_xlsx('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/5th_cohort/5500_FU_manifest.xlsx') %>%
  # total manifest, so have to filter and pick out EDD ret samples. 
  filter(STUDY_ID == 'EDD-ret') %>%
  transmute(`CMO sample ID` = gsub('_IGO.*.$','',CMO_SAMPLE_ID),Sex = SEX,
            `CMO patient ID` = CMO_PATIENT_ID, `Investigator sample ID` = INVESTIGATOR_SAMPLE_ID,
            `Patient ID` = gsub('-','_',gsub('_bc.*.$|-cf.*.$|_cf.*.$','',INVESTIGATOR_SAMPLE_ID))) %>% data.table()

DMP.mapping.3rd <- read_xlsx('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/3rd_cohort/17-256deidentified-master.xlsx') %>%
  transmute(MRN = gsub('No Impact',NA,`DMP Patient ID (P-7digit)`),`Patient ID`) %>% unique() %>% data.table()

manifest.5th <- merge(manifest.5th,DMP.mapping.3rd,by = 'Patient ID',all.x = T) %>% select(MRN,everything()) %>% select(-`Patient ID`)

write.csv(manifest.5th,'/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_mapping/DMP_mapping_5.csv',row.names = F)