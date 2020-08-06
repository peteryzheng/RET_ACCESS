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

diff.threshold <- 1000
merge(delly.calls,dmp.SV,by = 'DMP_ID',all.y = T,allow.cartesian=TRUE,suffixes = c('.delly','')) %>%
  # threshold for calling two events the same
  mutate(close.to.chr1pos = ifelse((Chr1.delly == Chr1 & abs(Pos1.delly - Pos1) <= diff.threshold) |
                                     (Chr2.delly == Chr1 & abs(Pos2.delly - Pos1) <= diff.threshold),TRUE,FALSE),
         close.to.chr2pos = ifelse((Chr1.delly == Chr2 & abs(Pos1.delly - Pos2) <= diff.threshold) |
                                     Chr2.delly == Chr2 & abs(Pos2.delly - Pos2) <= diff.threshold,TRUE,FALSE)) %>% 
  # filter for identical events
  mutate(identical.event = close.to.chr1pos&close.to.chr2pos) %>% filter(identical.event) %>% rowwise() %>%
  group_by(DMP_ID,EventType,Chr1,Pos1,Chr2,Pos2,Gene1,Gene2,Gene1desc,Gene2desc,EventInfo,Breakpoint_Type) %>% data.table() -> delly.plasma.called


manta.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/manta_032719/')
gene.list <- readLines('/ifs/work/bergerm1/zhengy1/access_resources/MSK-ACCESS-v1_0-probe-A_genelist.txt')

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_020519.csv') %>%
  mutate(BAM_path.plasma = gsub('__aln.*.duplex','',gsub('duplex_bams//|duplex_bams/','standard_bams/',BAM_path.plasma)),
         BAM_path.normal = gsub('__aln.*.duplex','',gsub('duplex_bams//|duplex_bams/','standard_bams/',BAM_path.normal))) %>% data.table()
master.ref[is.na(BAM_path.normal)]$BAM_path.normal <- '/ifs/work/bergerm1/Innovation/projects/zhengy1/MSK-Norm-F12-pl-T01_IGO_05500_FE_1_S1_001_cl_aln_srt_MD_IR_FX_BR.bam'
master.ref %>% group_by(DMP_ID) %>% summarise(plasma.samples = length(unique(`CMO sample ID.plasma`))) %>% filter(!is.na(DMP_ID)) %>% data.table()-> plasma.samples.of.patients

# plasma data and sorting -------------------------------------------------
vcf.filenames <- list.files(paste0(manta.dir,'vcf_inv_corrected_dir'),full.names = T)
vcf.filenames <- list.files(paste0(manta.dir,'vcf_dir'),full.names = T)
# vcf.filenames.to.read <- vcf.filenames[-which(file.info(vcf.filenames)[['size']] == 0)]
do.call(rbind,lapply(vcf.filenames,function(x){
  tmp.vcf <- fread(x,skip = '#CHROM')
  print(dim(tmp.vcf))
  if(nrow(tmp.vcf) > 0){
    tmp.vcf %>% mutate(Tumor_Sample_Barcode = colnames(tmp.vcf)[ncol(tmp.vcf)]) %>% data.table() %>%
      dplyr::rename(Tumor_Support := !!colnames(tmp.vcf)[ncol(tmp.vcf)],Normal_Support := !!colnames(tmp.vcf)[ncol(tmp.vcf)-1],Chr = `#CHROM`) %>%
      filter(Chr != 'GL000220.1') %>% data.table()
  }else{NULL}
})) %>% merge(master.ref[,.(Tumor_Sample_Barcode = Tumor_Sample_Barcode.plasma,DMP_ID)],all = T) %>% filter(!is.na(Chr)) %>% data.table() -> plasma.SV



# detection rate of fusion ------------------------------------------------
diff.threshold <- 5
merge(plasma.SV[!is.na(DMP_ID)],dmp.SV,by = 'DMP_ID',all.y = T,allow.cartesian=TRUE) %>%
  # threshold for calling two events the same
  mutate(close.to.chr1pos = ifelse(Chr == Chr1 & abs(POS - Pos1) <= diff.threshold,TRUE,FALSE),
         close.to.chr2pos = ifelse(Chr == Chr2 & abs(POS - Pos2) <= diff.threshold,TRUE,FALSE)) %>% 
  # filter for identical events
  # mutate(ID = paste0(unlist(str_extract_all(ID,"([0-9]{2,5})")),collapse = ':')) %>% unite(call_ID,Tumor_Sample_Barcode,ID) %>%
  group_by(DMP_ID,EventType,Chr1,Pos1,Chr2,Pos2,Gene1,Gene2,Gene1desc,Gene2desc,EventInfo,Breakpoint_Type,Tumor_Sample_Barcode) %>% 
  summarise(close.to.chr1pos = any(close.to.chr1pos),close.to.chr2pos = any(close.to.chr2pos)) %>%
  mutate(identical.event = close.to.chr1pos|close.to.chr2pos) %>% filter(identical.event) %>% data.table() -> manta.plasma.called


merge(delly.plasma.called[,.(TumorId,Pos1,Pos2,Chr1,Chr2,Gene1,Gene2,Gene1desc,Gene2desc)],
      manta.plasma.called[,.(Tumor_Sample_Barcode,Pos1,Pos2,Chr1,Chr2,Gene1,Gene2,Gene1desc,Gene2desc)],
      by = c('Pos1','Pos2','Chr1','Chr2','Gene1','Gene2','Gene1desc','Gene2desc'),suffixes = c('.delly','.manta'),all = T) %>% rowwise() %>%
  mutate(TumorId = ifelse(Tumor_Sample_Barcode == 'DA-ret-004-pl-T01_IGO_05500_FF_18' & TumorId == 'DA-ret-004-pl-T02_IGO_05500_FF_19',NA,TumorId)) %>%
  filter(ifelse(any(is.na(TumorId)|any(is.na(Tumor_Sample_Barcode))),T,TumorId == Tumor_Sample_Barcode)) %>% data.table() -> summary.table

library(VennDiagram)
common.items <- nrow(summary.table[TumorId == Tumor_Sample_Barcode])
delly.item <- nrow(summary.table[is.na(Tumor_Sample_Barcode)])
manta.item <- nrow(summary.table[is.na(TumorId)])
venn.diagram(list(c(1:(delly.item+common.items)),c((1+delly.item):(delly.item+common.items+manta.item))),
             category.names = c('Delly','Manta'),fill = c('purple','blue'),
             filename = '/ifs/work/bergerm1/zhengy1/RET_all/For_comparison/Venn/VennDiagram.png',
             output = TRUE ,imagetype="png" ,height = 1000 , width = 1000 , resolution = 300)
