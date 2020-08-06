library(data.table)

manta.dir <- '/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/manta_041519_2/'
annotate.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/annotate_',format(Sys.Date(),'%m%d%y'))
dir.create(annotate.dir)

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_020519.csv') 

# # excute ------------------------------------------------------------------
# 
# apply(master.ref,1,function(x){
#   tumor_path <- x[which(colnames(master.ref) == 'BAM_path.plasma')]
#   tumor.sample.name <- gsub('^.*.duplex_bams//|_IGO.*.$|_cl.*.$','',tumor_path)
#   dir.create(paste0(annotate.dir,'/',tumor.sample.name))
#   system(paste0(
#     'bsub -cwd ',annotate.dir,'/',tumor.sample.name,' -oo %J.o -eo %J.e -We 00:30',
#     '/ifs/work/bergerm1/zhengy1/RET_all/Code/manta_suite/iAnnotateSV.sh ',
#     manta.dir,tumor.sample.name,'/results/variants/somaticSV.vcf.gz ',
#     tumor.sample.name,' ',annotate.dir,'/',tumor.sample.name,' ',
#     '/ifs/work/bergerm1/zhengy1/RET_all/IannotateSV_testing/data'
#   ))
# })
# 



# Filters? ----------------------------------------------------------------

gene.list <- c('ALK','BRAF','ERBB2','EGFR','FGFR1','FGFR2','FGFR3','KIT','MET','NTRK1','NTRK2','NTRK3','PDGFRA','RET','ROS1')
annotate.filenames <- list.files(annotate.dir,'Annotated_Evidence.txt',full.names = T,recursive = T)
annotate.data <- do.call(rbind,lapply(annotate.filenames,fread)) %>% 
  filter(gene1 %in% gene.list | gene2 %in% gene.list) %>% data.table()
