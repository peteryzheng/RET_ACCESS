library(data.table)
library(ggplot2)
library(plotly)
library(stringr)
Sys.setenv("plotly_username"="zhengy1")
Sys.setenv("plotly_api_key"="gOuyJbekjN5xbzLKXARu")

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_020519.csv')
results.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_021319/')
manta.dir <- paste0('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/manta_021919/')

fusionAF <- fread('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/plasma_SV_FAF_022619.csv')

fusion.might.be.detected.plasma <- master.ref[DMP_ID %in% fread('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/called_in_plasma_SV022619.csv')$DMP_ID]$Tumor_Sample_Barcode.plasma

x <- list.files(paste0(results.dir,'results/'),full.names = T)[4]
max.af <- do.call(rbind,lapply(list.files(paste0(results.dir,'results_stringent/'),full.names = T),function(x){
  tmp.table <- fread(x)
  plasma.samples <- gsub('___duplex.called','',colnames(tmp.table)[grep('___duplex.called',colnames(tmp.table))])
  max.af.df <- data.table(Tumor_Sample_Barcode = plasma.samples,
                          max.AF = unlist(lapply(plasma.samples,function(x){
                            max(as.numeric(gsub('\\(|\\)','',str_extract(unlist(tmp.table[get(paste0(x,'___duplex.called')) %in% c('Called','Genotyped') & 
                                                                                            (Hugo_Symbol != 'APC' & HGVSp_Short != 'p.S142T') &
                                                                                            (Hugo_Symbol != 'GATA3' & !grepl('p..35',HGVSp_Short))
                                                                                          ,
                                                                                          eval(paste0(x,'___total')),with = F]),"\\(.*.\\)"))),na.rm = T)
                          })))
  max.af.df[max.af.df == -Inf]  <- 0 
  return(max.af.df)
}))

combine.af <- merge(max.af,fusionAF[,.(Tumor_Sample_Barcode,FAF)],all = T) %>%
  mutate(fusion.detected = ifelse(is.na(FAF),'No','Yes')) %>% filter(Tumor_Sample_Barcode %in% fusion.might.be.detected.plasma) %>% data.table()
combine.af[is.na(combine.af)] <- 0
pdf(paste0('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/AF_cor',format(Sys.time(),'%m%d%y'),'.pdf'),width = 8,height = 8)
ggplot(combine.af,aes(x = max.AF,y = FAF,color = fusion.detected)) + geom_point() +
  xlim(0,1) + ylim(0,1)
dev.off()
pdf(paste0('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/AF_cor_zoomed',format(Sys.time(),'%m%d%y'),'.pdf'),width = 8,height = 8)
ggplot(combine.af,aes(x = max.AF,y = FAF,color = fusion.detected)) + geom_point() +
  xlim(0,0.01) + ylim(0,0.05)
dev.off()

cor.test(combine.af$max.AF,combine.af$FAF,method = 'spearman')
cor.test(combine.af[fusion.detected == 'Yes']$max.AF,combine.af[fusion.detected == 'Yes']$FAF,method = 'spearman')

plot_ly(combine.af, x = ~max.AF, y = ~FAF, color = ~fusion.detected, type = 'scatter',
        mode = 'text', name = ~gsub('DA-ret-','',Tumor_Sample_Barcode)) %>% add_markers() %>%
  layout(title = 'AF Comparison',
         xaxis = list(title = 'max.AF',
                      zeroline = TRUE,
                      range = c(-0.1, 1)),
         yaxis = list(title = 'FAF',
                      range = c(-0.1,1)))
# paste0(merge(combine.af[FAF == 0],master.ref[,.(Tumor_Sample_Barcode = Tumor_Sample_Barcode.plasma,
#                                                 EDD.ID = str_extract(`Investigator sample ID.plasma`,'EDD-ret-pt...|EDD_ret_pt...|DA-ret-...'),
#                                                 DMP_ID)],
#              by = 'Tumor_Sample_Barcode',all.x = T) %>% arrange(max.AF) %>% unite(Tumor_Sample_Barcode,Tumor_Sample_Barcode,EDD.ID,DMP_ID,sep = ' -- ') %>%
#          select(Tumor_Sample_Barcode) %>% unlist() %>% as.vector(),collapse = '   ')
