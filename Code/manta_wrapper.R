library(data.table)
library(dplyr)
library(tidyr)

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_091719.csv') %>%
  mutate(BAM_path.plasma = gsub('__aln.*.duplex','_filter_sc',gsub('duplex_bams//|duplex_bams/','standard_bams_filtered/',BAM_path.plasma)),
         BAM_path.normal = gsub('__aln.*.duplex','_filter_sc',gsub('duplex_bams//|duplex_bams/','standard_bams_filtered/',BAM_path.normal))) %>% data.table()
  master.ref[is.na(BAM_path.normal)]$BAM_path.normal <- '/ifs/work/bergerm1/zhengy1/RET_all/Code/test_filter_sc/MSK-Norm-F12-pl-T01_IGO_05500_FE_1_S1_001_cl_aln_srt_MD_IR_FX_BR_filter_sc.bam'
manta.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/manta_',format(Sys.time(),'%m%d%y'))
dir.create(manta.dir)

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
  tumor_name <- gsub('^.*standard_bams_filtered/|_IGO.*.$|_cl.*.$','',tumor_path)
  dir.create(paste0(manta.dir,'/',tumor_name))
  system(paste0(
    '/ifs/work/bergerm1/zhengy1/manta/manta-1.5.0.centos6_x86_64/bin/configManta.py --normalBam ',
    normal_path ,' --tumorBam ',tumor_path,
    ' --runDir ',manta.dir,'/',tumor_name,' --exome ',
    '--referenceFasta /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta'
  ))
  system(paste0(
    'bsub -cwd ',manta.dir,'/',tumor_name,' -J ',tumor_name,
    ' -oo ',tumor_name,'.o ',' -eo ',tumor_name,'.e ',
    ' -W 4:00 -M 16 /opt/common/CentOS_7-dev/python/python-2.7.10/bin/python ',manta.dir,'/',
    tumor_name,'/runWorkflow.py -m local -j 8'
  ))
})

