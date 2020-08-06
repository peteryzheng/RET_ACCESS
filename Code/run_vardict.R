library(data.table)

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_080219.csv')
default.normal.bam <- '/ifs/work/bergerm1/zhengy1/RET_all/all_bams/duplex_bams//DA-ret-004-pl-T01_IGO_05500_FF_18_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam'
other.default.normal.bam <- '/ifs/work/bergerm1/zhengy1/RET_all/all_bams/duplex_bams//DA-ret-028-pl-T02_IGO_05500_FF_25_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam'

vardict.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/vardict_output/vardict_run_',format(Sys.time(),'%m%d%y'))
dir.create(vardict.dir)
system(paste0('cp /ifs/work/bergerm1/zhengy1/RET_all/Code/runVardict.sh ',vardict.dir))
setwd(paste0(vardict.dir))
system(paste0('echo $PWD'))

# running vardict (paired and unpaired)
# x <- as.character(master.ref[1,])
apply(master.ref,1,function(x){
  plasma.bam.path <- x[which(colnames(master.ref) == 'BAM_path.plasma')]
  # for using high-depth donor normal without spike in OR buffy coat when available
  # if(x[which(colnames(master.ref) == 'Paired')] == 'Paired'){
  #   normal.bam.path <- x[which(colnames(master.ref) == 'BAM_path.normal')]
  # }else{
  #   normal.bam.path <- default.normal.bam
  # }
  # for using spiked in low purity plasma sample
  normal.bam.path <- default.normal.bam
  if(plasma.bam.path == normal.bam.path){
    normal.bam.path <- other.default.normal.bam
  }
  system(paste0(
    vardict.dir,'/runVardict.sh ',plasma.bam.path,' ',normal.bam.path
  ))
})



