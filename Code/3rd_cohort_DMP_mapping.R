library(readxl)
library(data.table)
library(dplyr)

manifest_3rd <- read_xlsx('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/3rd_cohort/Proj_05500_FP_sample_manifest.xlsx') %>% 
  transmute(`CMO sample ID` = gsub('_IGO.*.$','',CMO_SAMPLE_ID),Sex = SEX,
            `CMO patient ID` = CMO_PATIENT_ID, `Investigator sample ID` = INVESTIGATOR_SAMPLE_ID,
            `Patient ID` = gsub('-','_',gsub('_bc.*.$|-cf.*.$','',INVESTIGATOR_SAMPLE_ID))) %>% data.table()

DMP.mapping.3rd <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_RET_ID/DMP_RET_ID_processed.txt') %>%
  transmute(MRN = gsub('No Impact',NA,MRN),`Patient ID`) %>% unique() %>% data.table()


manifest_3rd <- merge(manifest_3rd,DMP.mapping.3rd,by = 'Patient ID',all.x = T) %>% select(MRN,everything()) %>% select(-`Patient ID`)

write.csv(manifest_3rd,'/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_mapping/DMP_mapping_3.csv',row.names = F)