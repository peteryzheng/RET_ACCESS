library(data.table) 
library(tidyverse)
library(stringr)

master.ref = fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_121519.csv')
SV.pt = paste0('EDD_ret_',c('pt014','pt023','pt032','pt041','pt047','pt048','pt058','pt060'))

wes.manta.results = list.files('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/WES_manta_011520/final_results/',full.names = T)

lapply(wes.manta.results[grepl(paste0(unique(master.ref[EDD_Patient_ID %in% SV.pt]$`CMO patient ID`),collapse = '|'),wes.manta.results)],function(x){
  fread(x)
})

