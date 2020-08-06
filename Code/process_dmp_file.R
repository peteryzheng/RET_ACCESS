library(data.table)

DMP.mapping.file <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_RET_ID/DMP_RET_ID_080219.txt') %>%
  transmute(MRN = ifelse(grepl('^P-[0-9]+',`DMP Patient ID (P-7digit)`),str_extract(`DMP Patient ID (P-7digit)`,'^P-[0-9]+'),NA),
            `Patient ID` = gsub('_cf[0-9]+$','',`Unique sample ID recommended: PI initials-cohort-patient number-pl/bc-T/N number (pl-T for cfDNAs, bc-N for buffy coats)`)) %>% 
  group_by(`Patient ID`) %>% summarise(MRN = paste0(unique(MRN[!is.na(MRN)]),collapse = '')) %>% mutate(MRN = ifelse(MRN == '',NA,MRN)) %>% data.table()

write.csv(DMP.mapping.file,'/ifs/work/bergerm1/zhengy1/RET_all/Original_files/DMP_RET_ID/DMP_RET_ID_processed.txt',row.names = F)
