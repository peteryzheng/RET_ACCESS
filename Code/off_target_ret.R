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


# 019 ---------------------------------------------------------------------
ret.019.table = fread(list.files(paste0(results.dir,'/',combine.name),pattern = '019',full.names = T))
ret.019.sample.sheet = fread(list.files(paste0(results.dir,'/C-F5H7YA/'),pattern = 'sample_sheet',full.names = T)) %>% rowwise() %>%
  mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                              ifelse(Sample_Type == 'normal','duplexnormal',
                                     ifelse(Sample_Type == 'plasma_simple  x','simplex',Sample_Type)))) %>%
  mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
ret.019.calls = table_to_maf(ret.019.table,ret.019.sample.sheet)
ret.019.calls = data.table(process_maf_for_graph(ret.019.calls))

# likely oncgenic
# in 'Hugo_Symbol %in% ...' clause
# KRAS p.A59del not oncogenic
processed.ret.019.calls = ret.019.calls[grepl('L',Tumor_Sample_Barcode) & 
                                          HGVSp_Short != 'p.A59del' &
                                          (Hugo_Symbol %in% c('ARID1A','FGFR1','KEAP1','KRAS','PIK3CA','PTEN','SF3B1','TGFBR1'))]
colourCount = nrow(unique(processed.ret.019.calls[,.(Hugo_Symbol,HGVSp_Short)]))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))

pdf(paste0(output_dir,'/SNV_plot_019_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.019.calls) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 019',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()

pdf(paste0(output_dir,'/SNV_plot_019_condensed_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.019.calls[grepl('RET|KRAS|NRAS',Hugo_Symbol)]) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 019',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()


# 070 ---------------------------------------------------------------------
ret.070.table = fread(list.files(paste0(results.dir,'/',combine.name),pattern = '070',full.names = T))
ret.070.sample.sheet = fread(list.files(paste0(results.dir,'/C-DU7TVF/'),pattern = 'sample_sheet',full.names = T)) %>% rowwise() %>%
  mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                              ifelse(Sample_Type == 'normal','duplexnormal',
                                     ifelse(Sample_Type == 'plasma_simple  x','simplex',Sample_Type)))) %>%
  mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
ret.070.calls = table_to_maf(ret.070.table,ret.070.sample.sheet)
ret.070.calls = data.table(process_maf_for_graph(ret.070.calls))

# likely oncgenic
# in 'Hugo_Symbol %in% ...' clause
# KRAS p.A59del not oncogenic
processed.ret.070.calls = ret.070.calls[grepl('L',Tumor_Sample_Barcode) & 
                                          HGVSp_Short != 'p.A59del' &
                                          (Hugo_Symbol %in% c('ARID2','KEAP1','KRAS','STK11','TP53') | 
                                             grepl('RET',Hugo_Symbol))]
colourCount = nrow(unique(processed.ret.070.calls[,.(Hugo_Symbol,HGVSp_Short)]))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))

pdf(paste0(output_dir,'/SNV_plot_070_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.070.calls) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 070',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()

pdf(paste0(output_dir,'/SNV_plot_070_condensed_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.070.calls[grepl('RET|KRAS|NRAS',Hugo_Symbol)]) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 070',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()



# 058 ---------------------------------------------------------------------
ret.058.table = fread(list.files(paste0(results.dir,'/',combine.name),pattern = '058',full.names = T))
ret.058.sample.sheet = fread(list.files(paste0(results.dir,'/C-2FM5RX/'),pattern = 'sample_sheet',full.names = T)) %>% rowwise() %>%
  mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                              ifelse(Sample_Type == 'normal','duplexnormal',
                                     ifelse(Sample_Type == 'plasma_simple  x','simplex',Sample_Type)))) %>%
  mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
ret.058.calls = table_to_maf(ret.058.table,ret.058.sample.sheet)
ret.058.calls = data.table(process_maf_for_graph(ret.058.calls))

# likely oncgenic
# in 'Hugo_Symbol %in% ...' clause
# KRAS p.A59del not oncogenic
processed.ret.058.calls = ret.058.calls[grepl('L',Tumor_Sample_Barcode) & 
                                          (Hugo_Symbol %in% c('KRAS','PIK3CA','TP53') | 
                                             grepl('RET|DOT1L',Hugo_Symbol))]
colourCount = nrow(unique(processed.ret.058.calls[,.(Hugo_Symbol,HGVSp_Short)]))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))

pdf(paste0(output_dir,'/SNV_plot_058_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.058.calls) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 058',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()

pdf(paste0(output_dir,'/SNV_plot_058_condensed_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.058.calls[grepl('RET|KRAS|NRAS',Hugo_Symbol)]) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 058',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()


# 025 ---------------------------------------------------------------------
ret.025.table = fread(list.files(paste0(results.dir,'/',combine.name),pattern = '025',full.names = T))
ret.025.sample.sheet = fread(list.files(paste0(results.dir,'/C-4YT24X/'),pattern = 'sample_sheet',full.names = T)) %>% rowwise() %>%
  mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                              ifelse(Sample_Type == 'normal','duplexnormal',
                                     ifelse(Sample_Type == 'plasma_simple  x','simplex',Sample_Type)))) %>%
  mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
ret.025.calls = table_to_maf(ret.025.table,ret.025.sample.sheet)
ret.025.calls = data.table(process_maf_for_graph(ret.025.calls))

# likely oncgenic
# in 'Hugo_Symbol %in% ...' clause
# PMS2 mutations are artifact splice site mtuations
processed.ret.025.calls = ret.025.calls[grepl('L',Tumor_Sample_Barcode) &
                                          (Hugo_Symbol %in% c('BRAF','KDM5C','ZRSR2'))]
colourCount = nrow(unique(processed.ret.025.calls[,.(Hugo_Symbol,HGVSp_Short)]))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))

pdf(paste0(output_dir,'/SNV_plot_025_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.025.calls) +
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) +
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) +
        labs(title='RET patient 025',x='Time Point', y='VAF') +
        scale_shape_manual(values=status_id,name = 'Call Status') +
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()

# pdf(paste0(output_dir,'/SNV_plot_025_condensed_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
# print(ggplot(processed.ret.025.calls[grepl('RET|KRAS|NRAS',Hugo_Symbol)]) +
#         geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
#                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) +
#         geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
#                        color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) +
#         labs(title='RET patient 025',x='Time Point', y='VAF') +
#         scale_shape_manual(values=status_id,name = 'Call Status') +
#         scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
#         theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
# dev.off()





# 041 ---------------------------------------------------------------------
ret.041.table = fread(list.files(paste0(results.dir,'/',combine.name),pattern = '041',full.names = T))
ret.041.sample.sheet = fread(list.files(paste0(results.dir,'/C-J8DMAP/'),pattern = 'sample_sheet',full.names = T)) %>% rowwise() %>%
  mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                              ifelse(Sample_Type == 'normal','duplexnormal',
                                     ifelse(Sample_Type == 'plasma_simple  x','simplex',Sample_Type)))) %>%
  mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
ret.041.calls = table_to_maf(ret.041.table,ret.041.sample.sheet)
ret.041.calls = data.table(process_maf_for_graph(ret.041.calls))

# likely oncgenic
# in 'Hugo_Symbol %in% ...' clause
# AR not oncogenic
processed.ret.041.calls = ret.041.calls[grepl('L',Tumor_Sample_Barcode) & 
                                          (Hugo_Symbol %in% c('PTEN','TP53') | 
                                             grepl('RET|-',Hugo_Symbol))]
colourCount = nrow(unique(processed.ret.041.calls[,.(Hugo_Symbol,HGVSp_Short)]))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))

pdf(paste0(output_dir,'/SNV_plot_041_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.041.calls) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 041',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()

pdf(paste0(output_dir,'/SNV_plot_041_condensed_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.041.calls[grepl('RET|KRAS|NRAS',Hugo_Symbol)]) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 041',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()


# 028 ---------------------------------------------------------------------
ret.028.table = fread(list.files(paste0(results.dir,'/',combine.name),pattern = '028',full.names = T))
ret.028.sample.sheet = fread(list.files(paste0(results.dir,'/C-F6H84U/'),pattern = 'sample_sheet',full.names = T)) %>% rowwise() %>%
  mutate(Sample_Type = ifelse(Sample_Type == 'plasma','duplex',
                              ifelse(Sample_Type == 'normal','duplexnormal',
                                     ifelse(Sample_Type == 'plasma_simple  x','simplex',Sample_Type)))) %>%
  mutate(column.names = paste0(Sample_Barcode,'___',Sample_Type)) %>% data.table()
ret.028.calls = table_to_maf(ret.028.table,ret.028.sample.sheet)
ret.028.calls = data.table(process_maf_for_graph(ret.028.calls))

# likely oncgenic
# in 'Hugo_Symbol %in% ...' clause
# KRAS p.A59del not oncogenic
processed.ret.028.calls = ret.028.calls[grepl('L',Tumor_Sample_Barcode) & 
                                          (Hugo_Symbol %in% c('NRAS') | 
                                             grepl('RET',Hugo_Symbol))]
colourCount = nrow(unique(processed.ret.028.calls[,.(Hugo_Symbol,HGVSp_Short)]))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))

pdf(paste0(output_dir,'/SNV_plot_028_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.028.calls) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 028',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()

pdf(paste0(output_dir,'/SNV_plot_028_condensed_',format(Sys.time(),'%m%d%y'),'.pdf'),width = 10,height = 7)
print(ggplot(processed.ret.028.calls[grepl('RET|KRAS|NRAS',Hugo_Symbol)]) + 
        geom_line(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                      color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
        geom_point(aes(x = Tumor_Sample_Barcode, y = ifelse(t_total_count == 0, 0, as.numeric(t_alt_count/t_total_count)),
                       color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 1.5) + 
        labs(title='RET patient 028',x='Time Point', y='VAF') + 
        scale_shape_manual(values=status_id,name = 'Call Status') + 
        scale_color_manual(values = getPalette(colourCount),name = 'Alteration') +
        theme_minimal() + theme(panel.grid.major = element_blank()) + scale_y_log10() )
dev.off()

