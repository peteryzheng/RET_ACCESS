library(data.table)
library(tidyverse)

baseline.sample.pairs = fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/baseline_pairs_121519.tsv')

# pull data for both samples

# DMP.maf <- fread('/ifs/work/bergerm1/zhengy1/dmp/mskimpact/data_mutations_extended.txt') %>% 
#   filter(Mutation_Status == 'SOMATIC' & Tumor_Sample_Barcode %in% baseline.sample.pairs$DMP_ID) %>% data.table()
maf.dir = '/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_011520/results_combined_permissive/'
access.maf = fread(paste0(gsub('_combined','',maf.dir),'all_sample_maf2maf_oncokb.maf')) %>%
  filter(Tumor_Sample_Barcode %in% c(baseline.sample.pairs$`CMO sample ID.plasma`,baseline.sample.pairs$DMP_ID)) %>% 
  data.table()

x = baseline.sample.pairs$Study_ID[1]
x = 'EDD_ret_pt013'
total.detection.comp = do.call(rbind,lapply(baseline.sample.pairs$Study_ID,function(x){
  print(x)
  # get all possible mutation called in either or both samples
  all.unique.calls = access.maf[Tumor_Sample_Barcode %in% c(baseline.sample.pairs[Study_ID == x]$`CMO sample ID`,
                                                            baseline.sample.pairs[Study_ID == x]$DMP_ID)] %>% 
    transmute(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Reference_Allele,
              Tumor_Seq_Allele2,HGVSp_Short,oncogenic,call_confidence) %>%
    filter(call_confidence %in% c('Called','Genotyped')) %>% select(-call_confidence) %>% unique() %>% data.table()
  reviewed.calls = fread(paste0(maf.dir,'/',gsub('EDD_|pt','',x),'_table.csv'))[call_confidence == 'High',
                                                                                .(Hugo_Symbol,Chromosome = as.character(Chromosome),
                                                                                  Start_Position,End_Position,Variant_Classification,
                                                                                  Reference_Allele,Tumor_Seq_Allele2,HGVSp_Short)]
  all.unique.calls = merge(all.unique.calls,reviewed.calls,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification',
                                                                  'Reference_Allele','Tumor_Seq_Allele2','HGVSp_Short'))
  # get those all possible calls' genotypes in access sample
  access.tmp.maf = access.maf[Tumor_Sample_Barcode == baseline.sample.pairs[Study_ID == x]$`CMO sample ID`] %>% 
    transmute(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Reference_Allele,
              Tumor_Seq_Allele2,HGVSp_Short,t_alt_count,t_ref_count,t_depth,Tumor_Sample_Barcode,oncogenic,call_confidence) %>%
    merge(all.unique.calls,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification',
                                  'Reference_Allele','Tumor_Seq_Allele2','HGVSp_Short','oncogenic'))%>% data.table()
  # get those all possible calls' genotypes in DMP sample
  dmp.tmp.maf = access.maf[Tumor_Sample_Barcode == baseline.sample.pairs[Study_ID == x]$DMP_ID] %>% 
    transmute(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Reference_Allele,
              Tumor_Seq_Allele2,HGVSp_Short,t_alt_count,t_ref_count,t_depth,Tumor_Sample_Barcode,oncogenic,call_confidence) %>%
    merge(all.unique.calls,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification',
                                  'Reference_Allele','Tumor_Seq_Allele2','HGVSp_Short','oncogenic'))%>% data.table()
  detection.comp = merge(dmp.tmp.maf,access.tmp.maf,
                         by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification',
                                'Reference_Allele','Tumor_Seq_Allele2','HGVSp_Short','oncogenic'),
                         all = T,suffixes = c('.tumor','.plasma')) %>% rowwise() %>% 
    mutate(EDD_Patient_ID = x,detection_status = case_when(
      call_confidence.tumor != 'Not Called' & call_confidence.plasma != 'Not Called' ~ 'Both',
      call_confidence.tumor != 'Not Called' ~ 'Tumor',
      call_confidence.plasma != 'Not Called' ~ 'Plasma',
      TRUE ~ 'neither'
    )) %>% data.table()
}))


# graphing section --------------------------------------------------------
# artifact filtering for plasma
artifacts = total.detection.comp[detection_status == 'Plasma',.N,.(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Reference_Allele,Tumor_Seq_Allele2,oncogenic)]
total.detection.comp = merge(total.detection.comp,artifacts,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification',
                                                                   'Reference_Allele','Tumor_Seq_Allele2','oncogenic'), all = T) %>% 
  filter(N == 1 | is.na(N)) %>% data.table()

# overview of detection status
table(total.detection.comp$detection_status)

# graphing snippet
detection.graph.df = total.detection.comp[,.(oncogenic_mutation = length(which(oncogenic %in% c('Oncogenic','Likely Oncogenic'))),
                        non_oncogenic_mutation = length(which(!oncogenic %in% c('Oncogenic','Likely Oncogenic')))),detection_status] %>%
  melt(id.vars = 'detection_status',variable.name = 'Oncogenicity',value.name = 'Count')

pdf('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/RET_Detection/baseline_concordance.pdf',width = 7,height = 7)
ggplot(detection.graph.df) +
  geom_bar(aes(x = detection_status,y = Count,fill = Oncogenicity),stat = 'identity') + 
  scale_fill_brewer(palette = "Pastel2") + theme_classic()
dev.off()

detection.graph.df = total.detection.comp[,.(Count = .N),.(detection_status,Variant_Classification)] 
pdf('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/RET_Detection/baseline_concordance_class.pdf',width = 7,height = 7)
ggplot(detection.graph.df) +
  geom_bar(aes(x = detection_status,y = Count,fill = Variant_Classification),stat = 'identity') + 
  scale_fill_brewer(palette = "Pastel2") + theme_classic()
dev.off()


total.detection.comp[,.(tumor_only = nrow(.SD[detection_status == 'Tumor']),
                        plasma_only = nrow(.SD[detection_status == 'Plasma']),
                        both = nrow(.SD[detection_status == 'Both'])),EDD_Patient_ID]

colourCount = nrow(unique(total.detection.comp[,.(EDD_Patient_ID)]))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
pdf('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/RET_Detection/baseline_VAFs.pdf',width = 10,height = 7)
ggplot(total.detection.comp[detection_status == 'Both']) + 
  geom_point(aes(x = as.numeric(t_alt_count.tumor)/as.numeric(t_depth.tumor), 
                 y = as.numeric(t_alt_count.plasma)/as.numeric(t_depth.plasma),
                 color = EDD_Patient_ID,shape = oncogenic)) +
  scale_color_manual(values = getPalette(colourCount),name = 'Alteration') + theme_classic() + 
  xlab('Tumor_VAF') + ylab('Plasma_VAF')
dev.off()

pdf('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/RET_Detection/baseline_VAFs_all.pdf',width = 10,height = 7)
ggplot(total.detection.comp) + 
  geom_point(aes(x = as.numeric(t_alt_count.tumor)/as.numeric(t_depth.tumor), 
                 y = as.numeric(t_alt_count.plasma)/as.numeric(t_depth.plasma),
                 color = EDD_Patient_ID,shape = oncogenic)) +
  scale_color_manual(values = getPalette(colourCount),name = 'Alteration') + theme_classic() + 
  xlab('Tumor_VAF') + ylab('Plasma_VAF')
dev.off()


# analysis of each category -----------------------------------------------

table(total.detection.comp[detection_status == 'Plasma']$Hugo_Symbol)




