library(readxl)
library(data.table)
library(dplyr)

manifest.4th.a <- read_xlsx('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/4th_cohort/5500_FS_sample_manifest.xlsx') %>%
  # total manifest, so have to filter and pick out EDD ret samples. 
  filter(STUDY_ID == 'EDD-ret') %>% data.table()

manifest.4th.b <- read_xlsx('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/4th_cohort/Proj_05500_FP_sample_manifes_1lane.xlsx') %>%data.table()


manifest_4th <- rbind(manifest.4th.a,manifest.4th.b) %>% 
  transmute(`CMO sample ID` = gsub('_IGO.*.$','',CMO_SAMPLE_ID),Sex = SEX,
            `CMO patient ID` = CMO_PATIENT_ID, `Investigator sample ID` = INVESTIGATOR_SAMPLE_ID,
            `Patient ID` = gsub('-','_',gsub('_bc.*.$|-cf.*.$|_cf.*.$','',INVESTIGATOR_SAMPLE_ID))) %>% data.table()

DMP.mapping.3rd <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_RET_ID/DMP_RET_ID_processed.txt') %>%
  transmute(MRN = gsub('No Impact',NA,MRN),`Patient ID`) %>% unique() %>% data.table()


manifest_4th <- merge(manifest_4th,DMP.mapping.3rd,by = 'Patient ID',all.x = T) %>% select(MRN,everything()) %>% select(-`Patient ID`)

write.csv(manifest_4th,'/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_mapping/DMP_mapping_4.csv',row.names = F)