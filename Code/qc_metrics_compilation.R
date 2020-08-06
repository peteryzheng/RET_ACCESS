library(data.table)
library(tidyverse)

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/ret_study_id.csv') %>%
  filter(grepl('^C-',cmo_sample_id) & !grepl('^P-',study_id) & study_id != "#N/A") %>%
  mutate(cmo_sample_id = ifelse(grepl('DA-ret',inv_sample_id),inv_sample_id,cmo_sample_id)) %>% data.table()

# project.directories = c(
#   "/ifs/work/bergerm1/ACCESS-Projects/5500-FF/5500-FF_EDD_ret-1.1.1",
#   "/ifs/work/bergerm1/ACCESS-Projects/5500-FJ/5500-FJ_EDD_ret-1.1.1",
#   "/ifs/work/bergerm1/ACCESS-Projects/5500-FP/*",
#   "/ifs/work/bergerm1/ACCESS-Projects/5500-FS/5500-FS_EDD_ret-1.1.1",
#   "/ifs/work/bergerm1/ACCESS-Projects/5500-FU/5500-FU-EDD-ret-1.1.1",
#   "/work/bergerm1/MSK-ACCESS/ACCESS-Projects/05500_GB/EDD-ret-1.1.14",
#   "/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/05500_GD/EDD_ret-1.1.14"
# )

project.directories = c(
  "/ifs/work/bergerm1/ACCESS-Projects/5500-FF/5500-FF_EDD_ret-1.1.1",
  "/ifs/work/bergerm1/ACCESS-Projects/5500-FJ/5500-FJ_EDD_ret-1.1.1",
  "/juno/work/access/production/runs/Project_05500_FP/bam_qc/5500_FP_EDD_RET_lane_1-1.1.1",
  "/juno/work/access/production/runs/Project_05500_FP/bam_qc/5500_FP_EDD_RET_lane_3-1.1.1",
  "/juno/work/access/production/runs/Project_05500_FP/bam_qc/5500_FP_EDD_RET_lane_4-1.1.1",
  "/juno/work/access/production/runs/Project_05500_FP/bam_qc/5500_FP_EDD_RET_lane_5-1.1.1",
  "/juno/work/access/production/runs/Project_05500_FP/bam_qc/5500_FP_EDD_RET_lane_6-1.1.1",
  "/juno/work/access/production/runs/Project_05500_FP/bam_qc/5500_FP_EDD_RET_lane_7-1.1.1_fixed",
  "/juno/work/access/production/runs/Project_05500_FP/bam_qc/5500_FP_EDD_RET_lane_8-1.1.1",
  "/juno/work/access/production/runs/Project_05500_FP/bam_qc/5500-FQ_ret-1.1.1",
  "/ifs/work/bergerm1/ACCESS-Projects/5500-FS/5500-FS_EDD_ret-1.1.1",
  "/ifs/work/bergerm1/ACCESS-Projects/5500-FU/5500-FU-EDD-ret-1.1.1",
  "/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/05500_GB/EDD-ret-1.1.14",
  "/juno/work/bergerm1/MSK-ACCESS/ACCESS-Projects/05500_GD/EDD_ret-1.1.14"
)



# lapply(unique(paste0(gsub('/+[_0-9A-Za-z\\.\\-]+$','',project.directories),'/metadata')),function(x){list.files(x,'txt',full.names = T)})

title.paths = c(
  "/ifs/work/bergerm1/ACCESS-Projects/5500-FF/metadata/FF_RET_title_file.txt",
  "/ifs/work/bergerm1/ACCESS-Projects/5500-FJ/metadata/EDD-ret_FJ_title_file.txt",
  "/juno/work/access/production/runs/Project_05500_FP/metadata/EDD_RET_lane_1_title_file.txt",
  "/juno/work/access/production/runs/Project_05500_FP/metadata/EDD_RET_lane_3_title_file.txt",
  "/juno/work/access/production/runs/Project_05500_FP/metadata/EDD_RET_lane_4_title_file.txt",
  "/juno/work/access/production/runs/Project_05500_FP/metadata/EDD_RET_lane_5_title_file.txt",
  "/juno/work/access/production/runs/Project_05500_FP/metadata/EDD_RET_lane_6_title_file.txt",
  "/ifs/work/bergerm1/zhengy1/RET_all/lane7/metadata/EDD-ret_lane-7_title_file.txt",
  "/juno/work/access/production/runs/Project_05500_FP/metadata/EDD_RET_lane_8_title_file.txt",
  "/juno/work/access/production/runs/Project_05500_FP/metadata/FQ_ret_title_file.txt",
  "/ifs/work/bergerm1/ACCESS-Projects/5500-FS/metadata/EDD-ret_FS_title_file.txt",
  "/ifs/work/bergerm1/ACCESS-Projects/5500-FU/metadata/EDD-ret_FU_title_file.txt",
  "/ifs/work/bergerm1/zhengy1/RET_all/Original_files/ProJ_05500_GB_sample_manifest.xlsx",
  "/ifs/work/bergerm1/zhengy1/RET_all/Original_files/Proj_05500_GD_sample_manifest_03092019.xlsx" 
)


# title files -------------------------------------------------------------

title.file = do.call(rbind,lapply(title.paths,function(x){
  if(grepl('xlsx',x)){
    data.table(read_xlsx(x))[grepl('^EDD',INVESTIGATOR_SAMPLE_ID),.(CMO_SAMPLE_ID = gsub('_IGO.*.','',CMO_SAMPLE_ID),
                                CMO_PATIENT_ID,SAMPLE_TYPE,`LIBRARY_INPUT[ng]`,`LIBRARY_YIELD[ng]`,`CAPTURE_INPUT[ng]`,
                            project_id = gsub('.*.ACCESS-Projects/','',gsub('/+[_0-9A-Za-z\\.\\-]+$','',x)))]
  }else{
    fread(x)[,.(CMO_SAMPLE_ID,CMO_PATIENT_ID,SAMPLE_TYPE,`LIBRARY_INPUT[ng]`,`LIBRARY_YIELD[ng]`,`CAPTURE_INPUT[ng]`,
                project_id = gsub('.*.ACCESS-Projects/','',gsub('/+[_0-9A-Za-z\\.\\-]+$','',x)))]
  }
}))
dim(title.file)

# coverage tables ---------------------------------------------------------

agg.tables = do.call(rbind,lapply(list.files(paste0(project.directories,'/QC_Results/aggregate_tables/'),'coverage_agg.txt',full.names = T),function(x){
  fread(x) %>% filter(method %in% c('Duplex','Simplex','TotalCoverage','All Unique')) %>% mutate(pool = gsub(' Targets','',pool),method = case_when(
    method %in% c('Duplex','Simplex') ~ paste0(method,'Coverage'),
    method == 'All Unique' ~ 'UnfilteredCoverage',TRUE ~ method
  )) %>% unite('column_name',method,pool) %>% data.table %>% dcast.data.table(CMO_SAMPLE_ID ~ column_name,value.var = 'average_coverage')
}))
dim(agg.tables)

# family size -------------------------------------------------------------
# family type x pool a/b -- duplex, simplex, singletons, sub_simplex
pool_a.family.size.tables = do.call(rbind,lapply(list.files(paste0(project.directories,'/QC_Results'),'family-types-A.txt',full.names = T),function(x){
  dcast.data.table(fread(x),CMO_SAMPLE_ID ~ Type,value.var = 'Count')
}))
pool_b.family.size.tables = do.call(rbind,lapply(list.files(paste0(project.directories,'/QC_Results'),'family-types-B.txt',full.names = T),function(x){
  dcast.data.table(fread(x),CMO_SAMPLE_ID ~ Type,value.var = 'Count')
}))

family.size.tables = merge(pool_a.family.size.tables,pool_b.family.size.tables,by = 'CMO_SAMPLE_ID',all = T,suffixes = c('_size_pool_a','_pool_b'))
dim(family.size.tables)

# family size median and medium
family.stats.table = do.call(rbind,lapply(list.files(paste0(project.directories,'/QC_Results'),'family-sizes.txt',full.names = T),function(x){
  fread(x)[,.(FamilySize_medium = sum(.SD$Frequency*.SD$FamilySize)/sum(.SD$Frequency),
              FamilySize_median = as.numeric(median(rep(.SD$FamilySize,.SD$Frequency)))
              ),.(CMO_SAMPLE_ID,FamilyType)] %>% melt.data.table(id.vars = c('CMO_SAMPLE_ID','FamilyType'),variable.name = 'stats',value.name = 'value') %>%
    unite('column_name',FamilyType,stats) %>% dcast.data.table(CMO_SAMPLE_ID ~ column_name, value.var = 'value') %>% data.table()
}))
dim(family.stats.table)


# noise profile -----------------------------------------------------------
# by duplex/standard x pool a
# AltCount -- # of variant genotype + site combination
# AltPercent -- % of variant genotype + site combination out of all genotype + site combination
# ContributingSites -- # of noise contributing sites 
duplex_pool_a.noise.tables = do.call(rbind,lapply(list.files(paste0(project.directories,'/QC_Results/waltz_duplex_pool_a'),'noise.txt',full.names = T),function(x){
  fread(x)[Method == 'Total',.(CMO_SAMPLE_ID = gsub('_cl_aln.*.','',CMO_SAMPLE_ID),AltCount,AltPercent,ContributingSites)]
}))
standard_pool_a.noise.tables = do.call(rbind,lapply(list.files(paste0(project.directories,'/QC_Results/waltz_standard_pool_a'),'noise.txt',full.names = T),function(x){
  fread(x)[Method == 'Total',.(CMO_SAMPLE_ID = gsub('_cl_aln.*.','',CMO_SAMPLE_ID),AltCount,AltPercent,ContributingSites)]
}))
noise.tables = merge(duplex_pool_a.noise.tables,standard_pool_a.noise.tables,by = 'CMO_SAMPLE_ID',all = T,suffixes = c('_duplex_a','_standard_a')) 
dim(noise.tables)

# alignment ---------------------------------------------------------------
# by duplex/standard x pool a/b
# duplication rate, total/unique on target rate
duplex_pool_a.alignment.tables = do.call(rbind,lapply(list.files(paste0(project.directories,'/QC_Results/waltz_duplex_pool_a'),'read-counts.txt',full.names = T),function(x){
  fread(x)[,.(CMO_SAMPLE_ID,DuplicateFraction,TotalOnTargetFraction,UniqueOnTargetFraction)]
}))
duplex_pool_b.alignment.tables = do.call(rbind,lapply(list.files(paste0(project.directories,'/QC_Results/waltz_duplex_pool_b'),'read-counts.txt',full.names = T),function(x){
  fread(x)[,.(CMO_SAMPLE_ID,DuplicateFraction,TotalOnTargetFraction,UniqueOnTargetFraction)]
}))
standard_pool_a.alignment.tables = do.call(rbind,lapply(list.files(paste0(project.directories,'/QC_Results/waltz_standard_pool_a'),'read-counts.txt',full.names = T),function(x){
  fread(x)[,.(CMO_SAMPLE_ID,DuplicateFraction,TotalOnTargetFraction,UniqueOnTargetFraction)]
}))
standard_pool_b.alignment.tables = do.call(rbind,lapply(list.files(paste0(project.directories,'/QC_Results/waltz_standard_pool_b'),'read-counts.txt',full.names = T),function(x){
  fread(x)[,.(CMO_SAMPLE_ID,DuplicateFraction,TotalOnTargetFraction,UniqueOnTargetFraction)]
}))
alignment.tables = merge(duplex_pool_a.alignment.tables,duplex_pool_b.alignment.tables,by = 'CMO_SAMPLE_ID',all = T,suffixes = c('_duplex_a','_duplex_b')) %>%
  merge(standard_pool_a.alignment.tables,by = 'CMO_SAMPLE_ID',all = T,suffixes = c('','_standard_a')) %>%
  merge(standard_pool_b.alignment.tables,by = 'CMO_SAMPLE_ID',all = T,suffixes = c('','_standard_b')) %>% data.table()
dim(alignment.tables)


# total QC stats ----------------------------------------------------------
total.qc = merge(master.ref[,.(CMO_SAMPLE_ID = cmo_sample_id,study_id,inv_sample_id)],title.file[,!c('project_id'),with = F],all = T,by = 'CMO_SAMPLE_ID') %>%
  merge(agg.tables,all.x = T,by = 'CMO_SAMPLE_ID') %>%
  merge(family.size.tables,all.x = T,by = 'CMO_SAMPLE_ID') %>%
  merge(family.stats.table,all.x = T,by = 'CMO_SAMPLE_ID') %>%
  merge(noise.tables,all.x = T,by = 'CMO_SAMPLE_ID') %>%
  merge(alignment.tables,all.x = T,by = 'CMO_SAMPLE_ID')

write.csv(total.qc,'/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/QC_stats.csv',row.names = F)


