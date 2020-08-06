library(data.table)

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_021019.csv')
gene.list <- readLines('/ifs/work/bergerm1/zhengy1/access_resources/MSK-ACCESS-v1_0-probe-A_genelist.txt')

delly.outputs <- list.files('/ifs/work/bergerm1/SV_Testing/juber/iCallSV/RET/gene-filtered','final',full.names = T)
delly.calls <- do.call(rbind,lapply(delly.outputs,function(x){fread(x)})) %>% 
  merge(master.ref[,.(Tumor_Sample_Barcode.plasma,DMP_ID)],by.x = 'TumorId',by.y = 'Tumor_Sample_Barcode.plasma',all.x = T) %>% data.table()

dmp.SV <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/all_fusions.txt') %>% filter(EventInfo != 'Antisense Fusion') %>%
  mutate(DMP_ID = gsub('-T0.-IM.','',DMP_SAMPLE_ID)) %>% filter(Gene1 %in% gene.list | Gene2 %in% gene.list) %>%
  dplyr::select(-c(DMP_SAMPLE_ID,PairedReadCount,SplitReadCount)) %>% unique() %>%
  data.table()

diff.threshold <- 5
merge(delly.calls,dmp.SV,by = 'DMP_ID',all.y = T,allow.cartesian=TRUE,suffixes = c('.delly','')) %>%
  # threshold for calling two events the same
  mutate(close.to.chr1pos = ifelse((Chr1.delly == Chr1 & abs(Pos1.delly - Pos1) <= diff.threshold) |
                                     (Chr2.delly == Chr1 & abs(Pos2.delly - Pos1) <= diff.threshold),TRUE,FALSE),
         close.to.chr2pos = ifelse((Chr1.delly == Chr2 & abs(Pos1.delly - Pos2) <= diff.threshold) |
                                     Chr2.delly == Chr2 & abs(Pos2.delly - Pos2) <= diff.threshold,TRUE,FALSE)) %>% 
  # filter for identical events
  mutate(identical.event = close.to.chr1pos&close.to.chr2pos) %>% filter(identical.event) %>% rowwise() %>%
  group_by(DMP_ID,EventType,Chr1,Pos1,Chr2,Pos2,Gene1,Gene2,Gene1desc,Gene2desc,EventInfo,Breakpoint_Type) %>% 
  summarise(called.in.plasma = length(unique(TumorId))) %>%
  merge(dmp.SV[,.(DMP_ID,EventType,Chr1,Pos1,Chr2,Pos2,Gene1,Gene2,Gene1desc,Gene2desc,EventInfo,Breakpoint_Type)],all.y = T) %>%
  mutate(called.in.plasma = ifelse(is.na(called.in.plasma),0,called.in.plasma)) %>% merge(plasma.samples.of.patients,by = 'DMP_ID',all.x = T) %>% data.table() -> plasma.called

# table(plasma.called$called.in.plasma)
print(paste0(
  'Detection rate of a event present in a tumor sample at any timepoint in all possible corresponding plasma sample -- ',
  sum(plasma.called$called.in.plasma)/sum(plasma.called$plasma.samples),' (very likely undercalling detection rate)'
))
print(paste0(
  'Detection rate of an event present in a tumor sample at any timepoint in at least 1 corresponding plasma sample -- ',
  length(which(plasma.called$called.in.plasma != 0))/nrow(plasma.called)
))
# lenient.calling.df <- plasma.called[,.(called.num = max(called.in.plasma)),by = .(DMP_ID,plasma.samples)]
# sum(lenient.calling.df$called.num)/sum(lenient.calling.df$plasma.samples)
# nrow(lenient.calling.df[called.num != 0])/nrow(lenient.calling.df)

write.csv(plasma.called,paste0('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/called_in_plasma_SV_delly_',format(Sys.time(),'%m%d%y'),'.csv'),row.names = F)

# look at the ret fusion calls reported in dmp ----------------------------
merge(plasma.SV[!is.na(DMP_ID)],dmp.SV,by = 'DMP_ID',all.y = T,allow.cartesian=TRUE) %>%
  # threshold for calling two events the same
  mutate(close.to.chr1pos = ifelse((Chr1.delly == Chr1 & abs(Pos1.delly - Pos1) <= diff.threshold) |
                                     (Chr2.delly == Chr1 & abs(Pos2.delly - Pos1) <= diff.threshold),TRUE,FALSE),
         close.to.chr2pos = ifelse((Chr1.delly == Chr2 & abs(Pos1.delly - Pos2) <= diff.threshold) |
                                     Chr2.delly == Chr2 & abs(Pos2.delly - Pos2) <= diff.threshold,TRUE,FALSE)) %>% 
  # filter for identical events
  mutate(identical.event = close.to.chr1pos&close.to.chr2pos) %>% filter(identical.event) %>% rowwise() %>%
  # concatenate call and sample ID
  mutate(ID = paste0(unlist(str_extract_all(ID,"([0-9]{2,5})")),collapse = ':')) %>% unite(call_ID,Tumor_Sample_Barcode,ID,remove = F) %>% data.table() -> vcf.plasma.called
unique(vcf.plasma.called[,.(DMP_ID,Tumor_Sample_Barcode,call_ID,Tumor_Support,Pos1,Pos2,Gene1,Gene2,Gene1desc,Gene2desc)]) %>% rowwise() %>% 
  mutate(FAF = sum(rbindlist(lapply(strsplit(unlist(strsplit(Tumor_Support,':')),','),function(x){data.frame(alt = as.numeric(x[2]),ref = as.numeric(x[1]))}))$alt)/
           sum(rbindlist(lapply(strsplit(unlist(strsplit(Tumor_Support,':')),','),function(x){data.frame(alt = as.numeric(x[2]),ref = as.numeric(x[1]))}))$ref)) %>%
  group_by(DMP_ID,Tumor_Sample_Barcode,Pos1,Pos2,Gene1,Gene2,Gene1desc,Gene2desc) %>% summarise(FAF = sum(FAF)/n()) %>% data.table() -> called.fusion.AF
write.csv(called.fusion.AF,paste0('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/plasma_SV_FAF_',format(Sys.time(),'%m%d%y'),'.csv'),row.names = F)

# all possible ret mutations ----------------------------------------------
merge(plasma.SV[!is.na(DMP_ID)],dmp.SV,by = 'DMP_ID',all = T,allow.cartesian=TRUE) %>%
  mutate(close.to.chr1pos = ifelse(Chr == Chr1 & abs(POS - Pos1) <= diff.threshold,TRUE,FALSE),
         close.to.chr2pos = ifelse(Chr == Chr2 & abs(POS - Pos2) <= diff.threshold,TRUE,FALSE)) %>% 
  mutate(identical.event = ifelse(is.na(close.to.chr1pos|close.to.chr2pos),FALSE,close.to.chr1pos|close.to.chr2pos)) %>% 
  select(-c(EventType,Chr1,Pos1,Chr2,Pos2,Gene1,Gene2,Gene1desc,Gene2desc,EventInfo,Breakpoint_Type,close.to.chr1pos,close.to.chr2pos)) %>%
  group_by(!!!syms(colnames(plasma.SV))) %>% summarise(detected = any(identical.event)) %>% filter(!detected) %>% ungroup() %>% rowwise() %>%
  mutate(ID = paste0(unlist(str_extract_all(ID,"([0-9]{2,5})")),collapse = ':')) %>% unite(call_ID,Tumor_Sample_Barcode,ID) %>% filter(!call_ID %in% vcf.plasma.called$call_ID) %>%
  filter(Chr == 10 & POS < 43625799 & POS > 43572475) %>% data.table() -> vcf.plasma.not.called
plasma.SV %>% rowwise() %>% mutate(ID = paste0(unlist(str_extract_all(ID,"([0-9]{2,5})")),collapse = ':')) %>% unite(call_ID,Tumor_Sample_Barcode,ID)  %>%
  filter(call_ID %in% vcf.plasma.not.called$call_ID) %>% data.table() -> plasma.not.called

