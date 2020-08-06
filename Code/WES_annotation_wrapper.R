library(data.table)
library(dplyr)
library(tidyr)

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_091719.csv') 
manta.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/WES_manta_',format(Sys.time(),'%m%d%y'))

# excute ------------------------------------------------------------------
manta.results <- list.files('/ifs/work/bergerm1/zhengy1/RET_all/WES_RET/somatic/structural_variants/manta/',pattern = 'gz$',full.names = T)
dir.create(manta.dir)
lapply(manta.results,function(x){
  sample_id <- gsub('.*./manta//|__.*.','',x)
  dir.create(paste0(manta.dir,'/',sample_id))
  system(paste0(
    'bsub -cwd ',manta.dir,'/',sample_id,' -eo %J.o -oo %J.o -W 00:59 bash /ifs/work/bergerm1/zhengy1/RET_all/Code/ACCESS_SV/scripts/iAnnotateSV.sh ',
    x,' ',sample_id,' ',manta.dir,'/',sample_id,' /opt/common/CentOS_6-dev/manta/1.5.0 /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta ',
    ' /home/johnsoni/virtualenvs/pipeline_1.2.4/bin/python'
  ))
})



# # clean up step -----------------------------------------------------------
# out.log.filenames <- unlist(lapply(list.files(paste0(manta.dir,'/results'),full.names = T),function(x){list.files(x,pattern = '\\.o$',full.names = T)}))
# out.log <- unlist(lapply(out.log.filenames,function(x){
#   file.content <- readLines(x)
#   job.status <- gsub('^.*.cluster <juno. ','',file.content[which(grepl('Subject: ',file.content))])
# }))
# if(all(out.log == 'Done') & length(out.log) == nrow(master.ref)){
#   print('Everything ran correctly')
# }
# 
# manta.annot.results <- unlist(lapply(list.files(paste0(manta.dir,'/results'),full.names = T),function(x){list.files(x,pattern = 'Annotated_Evidence-annotated',full.names = T)}))
# dir.create(paste0(manta.dir,'/final_results/'))
# lapply(manta.annot.results, function(x){
#   system(paste0('cp ',x,' ',manta.dir,'/final_results/'))
# })
