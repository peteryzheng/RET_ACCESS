library(data.table)
library(stringr)
library(tidyverse)
source('/ifs/work/bergerm1/zhengy1/GDD/MolecularDiagnosis/feature_categories.R')

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_091719.csv')
data.CNA = fread('/ifs/work/bergerm1/zhengy1/dmp/mskimpact/data_CNA.txt')
data.seg = fread('/ifs/work/bergerm1/zhengy1/dmp/mskimpact/mskimpact_data_cna_hg19.seg')[grep(paste0(unique(master.ref$DMP_ID),collapse = '|'),ID)]
broad.cn = data.table(broad_cn(seg = data.seg)) %>% column_to_rownames('SAMPLE_ID') %>% t() %>% data.frame() %>% rownames_to_column('chrom_arm')

# focal change
focal.ret.impact.cna = data.CNA[,colnames(data.CNA)[which(grepl(paste(master.ref$DMP_ID[which(!is.na(master.ref$DMP_ID))],collapse = '|'),colnames(data.CNA)))],with = FALSE]
focal.flat.ret.impact = names(which(apply(focal.ret.impact.cna,2,function(x){ length(which(x == 0)) }) == nrow(focal.ret.impact.cna)))

# broad level change
broad.flat.ret.impact = gsub('\\.','-',names(which(apply(broad.cn,2,function(x){ length(which(x == 0)) }) == nrow(broad.cn))))

# intersecting samples with neither broad nor focal changes
flat.ret.impact = intersect(focal.flat.ret.impact,broad.flat.ret.impact)


# manual review of the results combining SNV, SV, and CNA
candidate.df = master.ref[DMP_ID %in% gsub('-T0.-IM.','',flat.ret.impact),.(`CMO sample ID.plasma`,Sex,`Investigator sample ID.plasma`)] %>%
  mutate(notes = case_when(
    grepl('02WK6K',`CMO sample ID.plasma`) ~ 'No SNV, RET SV not detected, CNV not conclusive',
    grepl('2XM3KC',`CMO sample ID.plasma`) ~ 'No SNV, RET SV not detected, CNV suspicious of chr1 arm level loss',
    grepl('30N2P4-L001',`CMO sample ID.plasma`) ~ 'RET SNV 1.2%, no SV, CNV flat',
    grepl('30N2P4-L002',`CMO sample ID.plasma`) ~ 'RET SNV 0.2%, no SV, a few probes might be elevated on chr2',
    grepl('7AVWDK',`CMO sample ID.plasma`) ~ 'A few SNV genotyped, RET SV not detected, CNV flat',
    grepl('DF0LWR',`CMO sample ID.plasma`) ~ 'No SNV, RET SV not detected, CNV flat',
    grepl('EW1Y3E',`CMO sample ID.plasma`) ~ 'Suspicious SNV, RET SV not detected, CNV flat',
    grepl('HD5CTA',`CMO sample ID.plasma`) ~ 'SNV below genotype threshold, RET SV not detected, CNV flat',
    grepl('LT4FT2',`CMO sample ID.plasma`) ~ 'SNV below genotype threshold, RET SV not detected, CNV flat',
    grepl('TF7VUA',`CMO sample ID.plasma`) ~ 'RET SNV 17.7%, no SV, CNV slightly suspicious',
    grepl('YP5R0K',`CMO sample ID.plasma`) ~ 'No SNV, RET SV not detected, CNV flact',
    grepl('YW82CY',`CMO sample ID.plasma`) ~ 'No SNV, no SV, CNV flat',
  ))

fake.normals = c(
  'C-30N2P4-L001-d','C-30N2P4-L002-d','C-7AVWDK-L001-d','C-EW1Y3E-L001-d', # Male
  'C-DF0LWR-L001-d','C-HD5CTA-L001-d','C-LT4FT2-L001-d','C-YP5R0K-L001-d','C-YW82CY-L001-d' # Female
)
