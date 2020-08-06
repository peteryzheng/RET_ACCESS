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
x = "/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_091819//results_combined_permissive/ret_006_table.csv"
do.call(rbind,lapply(list.files(results.dir,pattern = 'table.csv',full.names = T),function(x){
  print(x)
  tmp.data = fread(x)[grepl('__RET|RET__|^RET$',Hugo_Symbol)]
  fusion_expected = FALSE
  if(nrow(tmp.data) > 0){
    # # filtered for "called" rows
    # row.nums.called = unlist(lapply(1:nrow(tmp.data),function(y){
    #   # any RET event called at any timepoint is preserved for output
    #   if(any(tmp.data[y,(colnames(tmp.data)[grepl('duplex.called',colnames(tmp.data))]),with = F] == 'Called')) return(y)
    # }))
    # # filter for all called events
    # tmp.data = tmp.data[row.nums.called,]
    # 
    if(nrow(tmp.data) > 0){
      tmp.data = tmp.data %>% rowwise() %>% mutate(Hugo_Symbol = ifelse(Hugo_Symbol == 'RET',gsub('p.','',HGVSp_Short),Hugo_Symbol)) %>% data.table()
      return(data.frame(Patient_ID = str_extract(x,'ret_...'),
                        genetic_alterations = tmp.data$Hugo_Symbol,
                        AF = unite(tmp.data[,(colnames(tmp.data)[grepl('___total',colnames(tmp.data))]),with = F],col = 'AFs',sep = '|'),
                        Status = unite(tmp.data[,(colnames(tmp.data)[grepl('___duplex.called',colnames(tmp.data))]),with = F],col = 'Status',sep = '|')
      ))
    }
  }
})) %>% data.table() -> ret.detection.table

# separate the tables -- SNVs AFs dont need to be 'collapsed'
fusion.detection.table = ret.detection.table[grepl('__',genetic_alterations)]
snv.detection.table = ret.detection.table[!grepl('__',genetic_alterations)]

# process original hugo_symbol column (sort two genes by name)
fusion.detection.table$genetic_alterations = unlist(lapply(fusion.detection.table$genetic_alterations,function(x){paste0(sort(str_split(x,'__')[[1]]),collapse = '-')}))
# collapsing AF for rows of the same events (i.e. reciprocal rearrangement) while perserving the sample level seaparation in AF
fusion.detection.table = fusion.detection.table[,.(collapsed_AF = collapse_AF(AFs)),.(Patient_ID,genetic_alterations)]
# fusion call status also needs to be re-coded (0 == not called, !0 == Called)
fusion.detection.table$Status = unlist(lapply(fusion.detection.table$collapsed_AF,function(x){
  paste0(case_when(strsplit(x,'\\|')[[1]] == '0' ~ 'Not Called',TRUE ~ 'Called'),collapse = '|')
}))

# snv tables AFs need to be extracted -- only AFs in the 
snv.detection.table$collapsed_AF =  gsub('\\d+\\/\\d+|\\(|\\)','',snv.detection.table$AFs)

ret.detection.table = rbind(fusion.detection.table,snv.detection.table[,!c('AFs'),with = F])

ret.detection.table[,.SD %>% separate('collapsed_AF',sep = '\\|',paste0('T',1:(length(str_extract_all(.SD$collapsed_AF,'\\|')[[1]])+1),'_AF')) %>%
                      # list of list is the only way to return the per event melted table
                      # logic is to group by patient event, and melt the table by timepoint and then later dcast the table back so we have timepoint per column
                      melt(id.vars = c('Patient_ID','genetic_alterations','Status'),variable.name = 'time_point',value.name = 'AFs') %>% data.table() %>% list() %>% list(),
                    .SDcol = c('Patient_ID','genetic_alterations','Status','collapsed_AF'),by = .(Patient_ID,genetic_alterations,Status,collapsed_AF)] -> ret.detection.table
ret.detection.table = do.call(rbind,ret.detection.table$V1) %>% dcast(Patient_ID+genetic_alterations+Status ~ time_point, value.var = 'AFs')

# ret.detection.table %>% group_by(Patient_ID,genetic_alterations,Status,collapsed_AF) %>% mutate(timepoints = length(str_extract_all(collapsed_AF,'\\|')[[1]])) %>%
#   separate(collapsed_AF,paste0('T',1:(timepoints+1),'_AF'),sep = '\\|') #%>%
#   melt(id.vars = c('Patient_ID','genetic_alterations','Status'),variable.name = 'time_point',value.name = 'AFs')


# change patient ID
ret.detection.table$Patient_ID = gsub('ret_','EDD_ret_pt',ret.detection.table$Patient_ID)

# clean up -- merging, renaming
merge(ret.alt.pt[,c('Study_ID','DMP_ID','Alt_Type','genetic_alterations')],ret.detection.table,
      by.x = c('Study_ID','genetic_alterations'), by.y = c('Patient_ID','genetic_alterations'),all = T) %>% rowwise() %>%
  setnames(c('T1_AF','Status'),c('Baseline_AF','Notes')) %>%
  mutate(Notes = ifelse(is.na(Notes),'',Notes),DMP_ID = ifelse(is.na(DMP_ID),'',DMP_ID),Alt_Type = ifelse(is.na(Alt_Type),'',Alt_Type),
         Baseline_Alt = case_when(
           is.na(Notes) ~ 'No',grepl('^Called',Notes) ~ 'Yes',TRUE ~ 'No'
         )) %>% 
  select('Study_ID','DMP_ID','Alt_Type','genetic_alterations','Baseline_Alt','Baseline_AF',ends_with('AF'),'Notes') %>% data.table() -> ret.meta
# ret.meta[is.na(ret.meta)] = 0

ret.meta.indicated.ret = ret.meta[Alt_Type != '']
ret.meta.other.ret = ret.meta[Alt_Type == '']

# detection rate of all known targets
table(ret.meta.indicated.ret$Baseline_Alt)
# detection rate of all known fusion targets
table(ret.meta.indicated.ret[grepl('-',genetic_alterations)]$Baseline_Alt)
# 53.2%
# detection rate of all known SNV targets
table(ret.meta.indicated.ret[Alt_Type == 'Point']$Baseline_Alt)
# 68.4%

# manual review coded in 
ret.meta.indicated.ret[Study_ID %in% c('EDD_ret_pt061','EDD_ret_pt069','EDD_ret_pt024','EDD_ret_pt037',
                                       'EDD_ret_pt063','EDD_ret_pt038','EDD_ret_pt021','EDD_ret_pt027',
                                       'EDD_ret_pt028','EDD_ret_pt058','EDD_ret_pt065','EDD_ret_pt044',
                                       'EDD_ret_pt033')]$Baseline_Alt = 'Yes_review'

write.csv(ret.meta.indicated.ret,paste0("/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/RET_Detection/RET_indicated_detected_",format(Sys.time(),'%m%d%y'),'.csv'),quote = F,row.names = F)
write.csv(ret.meta.other.ret,paste0("/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/RET_Detection/RET_other_detected_",format(Sys.time(),'%m%d%y'),'.csv'),quote = F,row.names = F)

