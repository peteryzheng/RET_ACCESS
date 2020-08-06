library(data.table)
library(tidyverse)

master.ref = fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_121519.csv')
baseline.impact = fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/RET_DMP_tissue.txt')

write.table(merge(master.ref[grepl('L001',`CMO sample ID.plasma`),.(`CMO sample ID.plasma`,Study_ID = EDD_Patient_ID)],
                  baseline.impact,by = 'Study_ID'),
            paste0('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/baseline_pairs_',format(Sys.time(),'%m%d%y'),'.tsv'),
            sep = '\t',quote = F,row.names = F)
