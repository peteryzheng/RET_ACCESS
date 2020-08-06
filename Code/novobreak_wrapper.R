library(data.table)
library(dplyr)
library(tidyr)

nb.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/nb_',format(Sys.time(),'%m%d%y'))
dir.create(nb.dir)

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_020519.csv') %>%
  mutate(BAM_path.plasma = gsub('__aln.*.duplex','',gsub('duplex_bams//|duplex_bams/','standard_bams/',BAM_path.plasma)),
         BAM_path.normal = gsub('__aln.*.duplex','',gsub('duplex_bams//|duplex_bams/','standard_bams/',BAM_path.normal))) %>% data.table()
master.ref[is.na(BAM_path.normal)]$BAM_path.normal <- '/ifs/work/bergerm1/Innovation/projects/zhengy1/MSK-Norm-F12-pl-T01_IGO_05500_FE_1_S1_001_cl_aln_srt_MD_IR_FX_BR.bam'
# # Tumor bams --------------------------------------------------------------
# tumor.bam.filenames <- list.files('/ifs/work/bergerm1/ACCESS-Projects/5500-FJ/5500-FJ_EDD-ret-0.0.47/standard_bams/',
#                                   pattern = '.bam$',full.names = T)
# normal.bam.filenames <- '/ifs/work/bergerm1/Innovation/projects/zhengy1/MSK-Norm-F12-pl-T01_IGO_05500_FE_1_S1_001_cl_aln_srt_MD_IR_FX_BR.bam'
# 
# 
# excute ------------------------------------------------------------------

apply(master.ref,1,function(x){
  tumor_path <- x[which(colnames(master.ref) == 'BAM_path.plasma')]
  normal_path <- x[which(colnames(master.ref) == 'BAM_path.normal')]
  dir.create(paste0(nb.dir,'/',gsub('^.*.standard_bams/|_IGO.*.$|_cl.*.$','',tumor_path)))
  print(paste0(
    'bsub -n 8 -cwd ',nb.dir,'/',gsub('^.*.standard_bams/|_IGO.*.$|_cl.*.$','',tumor_path),' -o %J.o -e %J.e -R "rusage[mem=40]" ',
    'bash /ifs/work/bergerm1/SV_Testing/SV_callers/novoBreak/nb_distribution/run_novoBreak.sh ',
    '/ifs/work/bergerm1/SV_Testing/SV_callers/novoBreak/nb_distribution/ ',
    '/ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta ',
    tumor_path,' ',normal_path,' 8 ',nb.dir,'/',gsub('^.*.standard_bams/|_IGO.*.$|_cl.*.$','',tumor_path)
  ))
})

