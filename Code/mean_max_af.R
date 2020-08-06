library(data.table)
library(ggplot2)
library(plotly)
library(stringr)

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_020519.csv')
results.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_021319/')
manta.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/manta_021919/')

fusionAF <- fread('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/plasma_SV_FAF_022619.csv')

fusion.might.be.detected.plasma <- master.ref[DMP_ID %in% fread('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/called_in_plasma_SV022619.csv')$DMP_ID]$Tumor_Sample_Barcode.plasma

x <- list.files(paste0(results.dir,'results/'),full.names = T)[4]
mean.af <- do.call(rbind,lapply(list.files(paste0(results.dir,'results_stringent/'),full.names = T),function(x){
  tmp.table <- fread(x)
  print(tmp.table)
  plasma.samples <- gsub('___duplex.called','',colnames(tmp.table)[grep('___duplex.called',colnames(tmp.table))])
  print(plasma.samples)
  mean.af.df <- data.table(Tumor_Sample_Barcode = plasma.samples,
                          mean.AF = unlist(lapply(plasma.samples,function(y){
                            mean(as.numeric(gsub('\\(|\\)','',str_extract(unlist(tmp.table[(Hugo_Symbol != 'APC' & HGVSp_Short != 'p.S142T') &
                                                                                            (Hugo_Symbol != 'GATA3' & !grepl('p..35',HGVSp_Short) & 
                                                                                               DMP == 'Signed out')
                                                                                          ,
                                                                                          eval(paste0(y,'___total')),with = F]),"\\(.*.\\)"))),na.rm = T)
                          })),
                          max.AF = unlist(lapply(plasma.samples,function(y){
                            max(as.numeric(gsub('\\(|\\)','',str_extract(unlist(tmp.table[get(paste0(y,'___duplex.called')) %in% c('Called','Genotyped') & 
                                                                                            (Hugo_Symbol != 'APC' & HGVSp_Short != 'p.S142T') &
                                                                                            (Hugo_Symbol != 'GATA3' & !grepl('p..35',HGVSp_Short))
                                                                                          ,
                                                                                          eval(paste0(y,'___total')),with = F]),"\\(.*.\\)"))),na.rm = T)
                          })))
  mean.af.df[is.na(mean.af.df)]  <- 0 
  return(mean.af.df)
}))

x <- list.files(paste0(results.dir,'results_stringent/'),full.names = T)[2]
x <- '/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_021319/results_stringent//ret_056_SNV_table.csv'
max.af <- do.call(rbind,lapply(list.files(paste0(results.dir,'results_stringent/'),full.names = T),function(x){
  print(x)
  tmp.table <- fread(x)
  plasma.samples <- gsub('___duplex.called','',colnames(tmp.table)[grep('___duplex.called',colnames(tmp.table))])
  # max.af.df <- data.table(Tumor_Sample_Barcode = plasma.samples,
  #                         max.AF = unlist(lapply(plasma.samples,function(x){
  #                           max(as.numeric(gsub('\\(|\\)','',str_extract(unlist(tmp.table[get(paste0(x,'___duplex.called')) %in% c('Called','Genotyped') & 
  #                                                                                           (Hugo_Symbol != 'APC' & HGVSp_Short != 'p.S142T') &
  #                                                                                           (Hugo_Symbol != 'GATA3' & !grepl('p..35',HGVSp_Short))
  #                                                                                         ,
  #                                                                                         eval(paste0(x,'___total')),with = F]),"\\(.*.\\)"))),na.rm = T)
  #                         })))
  max.af.df <- do.call(rbind,lapply(plasma.samples,function(y){
    tmp.tmp.table <- tmp.table[get(paste0(y,'___duplex.called')) %in% c('Called','Genotyped') & 
                                 (Hugo_Symbol != 'APC' & HGVSp_Short != 'p.S142T') &
                                 (Hugo_Symbol != 'GATA3' & !grepl('p..35',HGVSp_Short))]
    af.of.sample <- as.numeric(gsub('\\(|\\)','',str_extract(unlist(tmp.tmp.table[,eval(paste0(y,'___total')),with = F]),"\\(.*.\\)")))
    print(which(af.of.sample == max(af.of.sample,na.rm = T))[1])
    if(!is.na(which(af.of.sample == max(af.of.sample,na.rm = T))[1])){
      return(tmp.tmp.table[which(af.of.sample == max(af.of.sample,na.rm = T))[1],1:12] %>% mutate(VAF = max(af.of.sample,na.rm = T),Sample = y))
    }else{
      empty.table <- data.frame(matrix(nrow = 1,ncol = 13))
      colnames(empty.table) <- c(colnames(tmp.tmp.table)[1:12],'VAF') 
      empty.table <- empty.table %>% mutate(Sample = y)
      return(empty.table)
    }
  }))
  max.af.df[max.af.df == -Inf]  <- 0 
  return(max.af.df)
}))

write.csv(max.af,'/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/max_vaf.csv',row.names = F)
