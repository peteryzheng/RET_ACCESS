library(data.table)
library(tidyverse)
library(ggplot2)
library(randomcoloR)

results.dir <- "/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_120919/results_stringent"
all.sample.maf = fread(paste0(results.dir,'/all_sample_maf2maf_oncokb.maf'))
output.dir = paste0('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/doublet_ret/graphs_',format(Sys.time(),'%m%d%y'))
dir.create(output.dir)
# ret_006,ret_035,ret_038,ret_044,
# ret_047,ret_054(fusion and mutation)
# ret_048 fusion and mutation but 1 timepoint

# doublet mutations
x = 'ret_035'
lapply(c('ret_006','ret_035','ret_038','ret_044','ret_054'),function(x){
  print(x)
  tmp.maf = all.sample.maf[Patient_ID == x & 
                             # oncogenic %in% c('Oncogenic','Likely Oncogenic') &
                             Hugo_Symbol == 'RET' &
                             # after review -- these two genes have the artifacts in the four patients
                             !Hugo_Symbol %in% c('PMS2','GATA3')]  
  # tumor sample 
  tumor.sample = structure(gsub('-','',str_extract(unique(tmp.maf$Tumor_Sample_Barcode[grep('IM[0-9]$',tmp.maf$Tumor_Sample_Barcode)]),'-T..-')),
                           names = unique(tmp.maf$Tumor_Sample_Barcode[grep('IM[0-9]$',tmp.maf$Tumor_Sample_Barcode)]))
  # rest of the samples are plasma
  plasma.sample = setdiff(tmp.maf$Tumor_Sample_Barcode,names(tumor.sample))
  # filter for plasma sample only
  tmp.maf = tmp.maf[Tumor_Sample_Barcode %in% plasma.sample]
  # change samples into timepoint information
  plasma.sample = structure(case_when(
    # some of the DA-ret sample need to be renamed
    grepl('-T0._',plasma.sample) ~ gsub('-|_','',gsub('T','L0',str_extract(plasma.sample,'-T0._'))),
    # otherwise extract L00 something
    TRUE ~ gsub('-','',str_extract(plasma.sample,'-L...-'))
  ),names = plasma.sample)
  sample.name.conversion = c(tumor.sample,plasma.sample)
  # get all not covered calls
  not.covered.df = unique(tmp.maf[call_confidence == 'Not Covered',.(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,
                                                                     HGVSp_Short,Reference_Allele,Tumor_Seq_Allele2)])
  # anti join for 'in tmp maf but not in not covered maf' --> only covered by access 
  only.covered.tmp.maf = anti_join(tmp.maf,not.covered.df,by = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification',
                                          'HGVSp_Short','Reference_Allele','Tumor_Seq_Allele2'))
  # only the converted timepoint names
  only.covered.tmp.maf$Tumor_Sample_Barcode = sample.name.conversion[only.covered.tmp.maf$Tumor_Sample_Barcode]
  # plotting
  status_id = c('Called' = 19, 'Not Called' =  4, 'Signed out' = 15,
                'Not Signed out' = 13, 'Not Covered' = 8, 'Genotyped' = 17)
  pdf(paste0(output.dir,'/SNV_plot_',x,'.pdf'))
  print(ggplot(only.covered.tmp.maf) + 
          geom_line(aes(x = Tumor_Sample_Barcode, y = as.numeric(t_alt_count/t_depth),
                        color = paste0(Hugo_Symbol,' ',HGVSp_Short),group = paste0(Hugo_Symbol,'_',HGVSp_Short))) + 
          geom_point(aes(x = Tumor_Sample_Barcode, y = as.numeric(t_alt_count/t_depth),
                         color = paste0(Hugo_Symbol,' ',HGVSp_Short),shape = call_confidence),size = 3) + 
          labs(title=unique(only.covered.tmp.maf$Patient_ID),x='Time Point', y='VAF') + 
          scale_shape_manual(values=status_id,name = 'Call Status') + 
          scale_color_brewer(palette = "Set1",name = 'Alteration') + 
          theme_minimal() + theme(panel.grid.major = element_blank()))
  dev.off()
})
