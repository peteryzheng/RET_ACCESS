library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
source('/ifs/work/bergerm1/zhengy1/RET_all/Code/table_to_maf.R')
# convert naming to timepoint, get rid of uncovered impact and access calls
source('/ifs/work/bergerm1/zhengy1/RET_all/Code/process_maf_for_graphs.R')

# three patients have on target acquired resistance mechanisms
output_dir = '/ifs/work/bergerm1/zhengy1/RET_all/resistance_cases/'

# for plotting consistency
status_id = c('Called' = 19, 'Not Called' =  4, 'Signed out' = 15,
              'Not Signed out' = 13, 'Not Covered' = 8, 'Genotyped' = 17)


# results directory
results.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_120919/')
# criteria <- 'stringent'
criteria <- 'permissive'

combine.name = paste0('results_combined_',criteria)


# 006 ---------------------------------------------------------------------
ret.006.table = fread(list.files(paste0(results.dir,'/',combine.name),pattern = '006',full.names = T))
ret.006.sample.sheet = fread(list.files(paste0(results.dir,'/C-8F9XD2/'),pattern = 'sample_sheet',full.names = T)) %>% rowwise() %>%
  mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                              ifelse(Sample_Type == 'normal','duplexnormal',
                                     ifelse(Sample_Type == 'plasma_simplex','simplex',Sample_Type)))) %>%
  mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
ret.006.calls = table_to_maf(ret.006.table,ret.006.sample.sheet)
ret.006.calls = data.table(process_maf_for_graph(ret.006.calls))

# likely oncgenic
## MLH1 splice site, NF1 splice site, NF1 nonsense (2), ALK,
## ATM (consistently high AF, signedout as well, not in normal), NFE2L2 (only permissive criteria)
## KRAS
# MLH and MSH2 missense unknown oncogenicity
processed.ret.006.calls = ret.006.calls[grepl('L',Tumor_Sample_Barcode) & 
                                          # RET indel is not real
                                          HGVSp_Short != 'p.V804Rfs*51' &
                                          HGVSp_Short != 'p.P640L' &
                                          (Hugo_Symbol %in% c('KRAS','MLH1','MSH2','RET','ALK','ATM') | 
                                             (Hugo_Symbol == 'NF1' & Variant_Classification %in% c('Splice_Site','Nonsense_Mutation'))
                                          )]
colourCount = nrow(unique(processed.ret.006.calls[,.(Hugo_Symbol,HGVSp_Short)]))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))

pdf(paste0(output_dir,'/SNV_plot_006_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.006.calls) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 006',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()

pdf(paste0(output_dir,'/SNV_plot_006_condensed_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.006.calls[grepl('RET|KRAS|NRAS',Hugo_Symbol)]) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 006',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()


# 047 --------------------------------------------------------------------
ret.047.table = fread(list.files(paste0(results.dir,'/',combine.name),pattern = '047',full.names = T))
ret.047.sample.sheet = fread(list.files(paste0(results.dir,'/C-898XKN/'),pattern = 'sample_sheet',full.names = T)) %>% rowwise() %>%
  mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                              ifelse(Sample_Type == 'normal','duplexnormal',
                                     ifelse(Sample_Type == 'plasma_simplex','simplex',Sample_Type)))) %>%
  mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
ret.047.calls = table_to_maf(ret.047.table,ret.047.sample.sheet)
ret.047.calls = data.table(process_maf_for_graph(ret.047.calls))

# likely oncogenic
## RB1, TP53, PIK3CA
# others are not oncogenic
processed.ret.047.calls = ret.047.calls[grepl('L',Tumor_Sample_Barcode) & 
                                          # this ret fusion is not real
                                          Hugo_Symbol != 'RET__ZNF487']
colourCount = nrow(unique(processed.ret.047.calls[,.(Hugo_Symbol,HGVSp_Short)]))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))

pdf(paste0(output_dir,'/SNV_plot_047_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 5)
print(ggplot(processed.ret.047.calls) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 047',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()

pdf(paste0(output_dir,'/SNV_plot_047_condensed_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 5)
print(ggplot(processed.ret.047.calls[grepl('RET|KRAS|NRAS',Hugo_Symbol)]) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 047',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()


# 054 --------------------------------------------------------------------
ret.054.table = fread(list.files(paste0(results.dir,'/',combine.name),pattern = '054',full.names = T))
ret.054.sample.sheet = fread(list.files(paste0(results.dir,'/C-6WY7PC/'),pattern = 'sample_sheet',full.names = T)) %>% rowwise() %>%
  mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                              ifelse(Sample_Type == 'normal','duplexnormal',
                                     ifelse(Sample_Type == 'plasma_simplex','simplex',Sample_Type)))) %>%
  mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
ret.054.calls = table_to_maf(ret.054.table,ret.054.sample.sheet)
ret.054.calls = data.table(process_maf_for_graph(ret.054.calls))

# likely oncogenic
## KRAS, MAP2K1, PIK3CA, RB1, SMARCA4, TP53
# others are not oncogenic
processed.ret.054.calls = ret.054.calls[grepl('L',Tumor_Sample_Barcode) & 
                                          Hugo_Symbol %in% c('KRAS', 'MAP2K1', 'PIK3CA', 'RB1', 'SMARCA4', 'TP53','RET',
                                                             'KIF5B-RET','CTNNBL1-ERBB2')]
colourCount = nrow(unique(processed.ret.054.calls[,.(Hugo_Symbol,HGVSp_Short)]))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))

pdf(paste0(output_dir,'/SNV_plot_054_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 5)
print(ggplot(processed.ret.054.calls) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 054',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()

pdf(paste0(output_dir,'/SNV_plot_054_condensed_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 5)
print(ggplot(processed.ret.054.calls[grepl('RET|KRAS|NRAS',Hugo_Symbol)]) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 054',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()





































