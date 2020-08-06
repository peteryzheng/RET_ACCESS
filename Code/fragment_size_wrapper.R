library(data.table)
library(tidyverse)

fragment.ref = fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/RET_fragment.txt')
master.ref = fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_030920.csv')
fragment.dir = paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/fragment_',format(Sys.time(),'%m%d%y'))
dir.create(fragment.dir)

fragment.ref = merge(fragment.ref,master.ref[,.(Tumor_Sample_Barcode.plasma,BAM_path.plasma)],
                     by.x = 'Sample',by.y = 'Tumor_Sample_Barcode.plasma',all.x = T)

# fragment.ref[Sample == 'C-N440DJ-L001-d']
# samtools view -bh /ifs/work/bergerm1/zhengy1/RET_all/all_bams/duplex_bams//C-N440DJ-L001-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam '10:43609933-43609947' > /ifs/work/bergerm1/zhengy1/RET_all/test_fragment_size/C-N440DJ-L001-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex_subset.bam
# samtools index /ifs/work/bergerm1/zhengy1/RET_all/test_fragment_size/C-N440DJ-L001-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex_subset.bam
# samtools view /ifs/work/bergerm1/zhengy1/RET_all/test_fragment_size/C-N440DJ-L001-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex_subset.bam | awk '{FS=OFS="\t"}{print $9}' > /ifs/work/bergerm1/zhengy1/RET_all/test_fragment_size/C-N440DJ-L001-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex_subset_tlen.txt
# samtools view /ifs/work/bergerm1/zhengy1/RET_all/test_fragment_size/C-N440DJ-L001-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex_subset.bam | awk '{FS=OFS="\t"}{print $6}' > /ifs/work/bergerm1/zhengy1/RET_all/test_fragment_size/C-N440DJ-L001-d_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex_subset_mstat.txt

# bam_str = args_list[0]
# sample_name = args_list[1]
# output_dir = args_list[2]
# chrom = args_list[3]
# start = args_list[4]
# end = args_list[5]
# ref_allele = args_list[6]
# alt_allele = args_list[7]
x = unlist(fragment.ref[1,])
apply(fragment.ref, 1, function(x){
  sample.name = x[which(colnames(fragment.ref) == 'Sample')]
  Gene = x[which(colnames(fragment.ref) == 'Gene')]
  Mutation = x[which(colnames(fragment.ref) == 'Mutation')]
  dir.create(paste0(fragment.dir,'/',sample.name,'__',Gene,'__',Mutation))
  system(paste0(
    'python /ifs/work/bergerm1/zhengy1/RET_all/Code/fragment_size.py ',
    x[which(colnames(fragment.ref) == 'BAM_path.plasma')],' ',
    sample.name,' ',
    fragment.dir,'/',sample.name,'__',Gene,'__',Mutation,' ',
    x[which(colnames(fragment.ref) == 'Chr')],' ',
    x[which(colnames(fragment.ref) == 'Start')],' ',
    x[which(colnames(fragment.ref) == 'End')],' ',
    x[which(colnames(fragment.ref) == 'Ref')],' ',
    x[which(colnames(fragment.ref) == 'Alt')],' '
  ))
  # mt.frag = readLines(paste0(fragment.dir,'/',sample.name,'/',sample.name,'_mt_frag.txt'))
  # wt.frag = readLines(paste0(fragment.dir,'/',sample.name,'/',sample.name,'_wt_frag.txt'))
})


