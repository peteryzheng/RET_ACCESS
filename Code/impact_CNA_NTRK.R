library(data.table)
library(dplyr)
library(tidyr)

cna.dir <- paste0('/ifs/work/bergerm1/wonh/MSK-ACCESS/5500-FU/cna_',format(Sys.time(),'%m%d%y'))
# cna.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/cna_',format(Sys.time(),'%m%d%y'),'_no_ret')
dir.create(cna.dir)


# manifests for tumor and normals -----------------------------------------
dir.create(paste0(cna.dir,'/manifests/'))
master.ref <- data.table(BAM_path.plasma = list.files('/ifs/work/bergerm1/ACCESS-Projects/5500-FU/EDD-ntrk_FU-1.1.1/unfiltered_bams/',pattern = '-L0.*.bam$',full.names = T),Sex = 'Female')
write.table(master.ref[,.(BAM_path.plasma = gsub('-duplex','',gsub('duplex_bams','unfiltered_bams',BAM_path.plasma)),Sex)],
            paste0(cna.dir,'/manifests/tumor_manifest.txt'),sep = '\t',quote = F,row.names = F,col.names = F)

# high depth normal (not perfect because no spike in ret probes) ----------
write.table(data.frame(bam_paths = list.files('/ifs/work/bergerm1/brannona/ACCESS_M1.8/ACCESSv1-VAL-20190004/unfiltered',pattern = 'DONOR[0-9]+-T.*.bam',full.names = T),Sex = 'Male'),
            paste0(cna.dir,'/manifests/normal_manifest.txt'),sep = '\t',quote = F,row.names = F,col.names = F)

# # curated low/no tumor content plasma bams --------------------------------
# write.table(master.ref[c(5:8),.(BAM_path.plasma = gsub('-duplex','',gsub('duplex_bams','unfiltered_bams',BAM_path.plasma)),Sex)],
#             paste0(cna.dir,'/manifests/normal_manifest.txt'),sep = '\t',quote = F,row.names = F,col.names = F)
# 
# running pipeline --------------------------------------------------------

system(paste0(
  'bsub -sla Berger -q sol -cwd ',cna.dir,' -J ','cna_',format(Sys.time(),'%m%d%y'),' -o %J.o -e %J.e',
  ' -We 24:00 -R "rusage[mem=8]" -M 8 -n 1 ',
  ' /home/ptashkir/.conda/envs/py27/bin/python /home/ptashkir/CNV_ACCESS/cBX_pipeline/scripts/cfdna_scna.py',
  ' -t ',cna.dir,'/manifests/tumor_manifest.txt',
  ' -n ',cna.dir,'/manifests/normal_manifest.txt',
  ' -tr 25 -b /ifs/work/bergerm1/zhengy1/RET_all/Code/ACCESS_CNV/data/ACCESS_targets_coverage.bed',
  ' -g /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta',
  ' -r /opt/common/CentOS_6/R/R-3.2.0/bin/R -q sol -o ',cna.dir,
  ' -bsub /common/lsf/9.1/linux2.6-glibc2.3-x86_64/bin/bsub -id EDD_ret',
  ' -l /home/ptashkir/CNV_ACCESS/cBX_pipeline/scripts/loessnormalize_nomapq_cfdna.R',
  ' -cn /home/ptashkir/CNV_ACCESS/cBX_pipeline/scripts/copynumber_tm.batchdiff_cfdna.R',
  ' -ta /ifs/work/bergerm1/zhengy1/RET_all/Code/ACCESS_CNV/data/ACCESS_targets_coverage.txt'
))

# bed.file <- fread('/ifs/work/bergerm1/zhengy1/RET/Original_file/MSK-ACCESS-v1_0-probe-A.sorted.RET.bed')
# write.table(bed.file[,!c('V4'),with = F],
#             '/ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0-probe-A.sorted.RET.bed',
#             sep = '\t',quote = F,row.names = F,col.names = F)
# system('bedtools nuc -fi /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta -bed /ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0-probe-A.sorted.RET.bed > /ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0-probe-A.sorted.RET.txt')


# # using non ret bed file --------------------------------------------------
# 
# system(paste0(
#   'bsub -sla Berger -q sol -cwd ',cna.dir,' -J ','cna_',format(Sys.time(),'%m%d%y'),'_no_ret -o %J.o -e %J.e',
#   ' -We 24:00 -R "rusage[mem=8]" -M 8 -n 1 ',
#   ' /home/ptashkir/.conda/envs/py27/bin/python /home/ptashkir/CNV_ACCESS/cBX_pipeline/scripts/cfdna_scna.py',
#   ' -t ',cna.dir,'/manifests/tumor_manifest.txt',
#   ' -n ',cna.dir,'/manifests/normal_manifest.txt',
#   ' -tr 25 -b /ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0.sorted.woRET.bed',
#   ' -g /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta',
#   ' -r /opt/common/CentOS_6/R/R-3.2.0/bin/R -q sol -o ',cna.dir,
#   ' -bsub /common/lsf/9.1/linux2.6-glibc2.3-x86_64/bin/bsub -id EDD_ret',
#   ' -l /home/ptashkir/CNV_ACCESS/cBX_pipeline/scripts/loessnormalize_nomapq_cfdna.R',
#   ' -cn /home/ptashkir/CNV_ACCESS/cBX_pipeline/scripts/copynumber_tm.batchdiff_cfdna.R',
#   ' -ta /ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0.sorted.woRET.merged.txt'
# ))
# 
