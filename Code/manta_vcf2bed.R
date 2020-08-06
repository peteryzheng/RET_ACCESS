library(data.table)
library(dplyr)

manta.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/manta_050919/')

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_020519.csv') %>%
  mutate(BAM_path.plasma = gsub('__aln.*.duplex','',gsub('duplex_bams//|duplex_bams/','standard_bams/',BAM_path.plasma))) %>% data.table()

# Tumor bams --------------------------------------------------------------
# tumor.bam.filenames <- list.files('/ifs/work/bergerm1/ACCESS-Projects/5500-FJ/5500-FJ_EDD-ret-0.0.47/standard_bams',
#                                   pattern = '.bam$',full.names = T)

dir.create(paste0(manta.dir,'/vcf_dir'))
dir.create(paste0(manta.dir,'/vcf_inv_corrected_dir'))
dir.create(paste0(manta.dir,'/vcf_ct_inv'))
apply(master.ref,1,function(x){
  tumor_path <- x[which(colnames(master.ref) == 'BAM_path.plasma')]
  normal_path <- x[which(colnames(master.ref) == 'BAM_path.normal')]
  vcf.gz.filenames <- list.files(paste0(manta.dir,'/',gsub('^.*.standard_bams/|_IGO.*.$|_cl.*.$','',tumor_path),
                                        '/results/variants'),pattern = 'somaticSV.vcf.gz$',full.names = T)
  system(paste0('zcat ',vcf.gz.filenames,' > ',manta.dir,'/vcf_dir','/',gsub(paste0(c(manta.dir,'/results/variants','/','.gz'),collapse = '|'),'',vcf.gz.filenames)))
  # system(paste0(
  #   'python  ~/bin/svtools/vcfToBedpe -i ',manta.dir,'/vcf_dir','/',gsub(paste0(c(manta.dir,'/results/variants','/','.gz'),collapse = '|'),'',vcf.gz.filenames),
  #   ' -o ',manta.dir,'/bedpe_dir/',gsub('^.*.manta_....../|/results.*.$','',vcf.gz.filenames),'_vcf2maf.bedpe'
  #   # ' --ref-fasta /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta ',
  #   # ' --tmp-dir ',manta.dir,'/vcf2maf_tmp_dir/'
  # ))
  system(paste0(
    '/opt/common/CentOS_6-dev/manta/1.5.0/libexec/convertInversion.py /opt/common/CentOS_6-dev/bin/current/samtools ',
    ' /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta ',
    manta.dir,'/vcf_dir','/',gsub(paste0(c(manta.dir,'/results/variants','/','.gz'),collapse = '|'),'',vcf.gz.filenames),' > ',
    manta.dir,'/vcf_inv_corrected_dir','/',gsub(paste0(c(manta.dir,'/results/variants','/','.gz'),collapse = '|'),'',vcf.gz.filenames)
  ))
  # system(paste0(
  #   'Rscript /ifs/work/bergerm1/zhengy1/RET_all/Code/vcf_ct_edits.R -v ',
  #   manta.dir,'/vcf_inv_corrected_dir','/',gsub(paste0(c(manta.dir,'/results/variants','/','.gz'),collapse = '|'),'',vcf.gz.filenames),' -o ',
  #   manta.dir,'/vcf_ct_inv','/',gsub(paste0(c(manta.dir,'/results/variants','/','.gz'),collapse = '|'),'',vcf.gz.filenames)
  # ))
})



