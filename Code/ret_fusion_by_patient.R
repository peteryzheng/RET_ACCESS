library(data.table)
library(stringr)

x = c("32-117-190|53-0-959","NA|16-0-1035","63-0-954|NA" )
x = c('3-5-1300')
x = c('NA|NA|NA|0-56-4183','NA|NA|NA|3-4-76')
collapse_AF = function(x){
  print(paste0(x,collapse = "','"))
  # number of samples
  samples.num = attr(gregexpr('\\|',x[1])[[1]],"match.length")
  print(samples.num)
  # when not found, return -1
  if(samples.num != -1){
    paste0(round(apply(separate(data.frame(AF = x),'AF',paste0('timpoint_',c(1:(length(samples.num)+1))),sep = '\\|'),2,function(one.tp.afs){
      mean(unlist(lapply(one.tp.afs,function(tmp.af){
        if(tmp.af == 'NA'){return(0)}
        else{
          read.counts = as.numeric(str_split(tmp.af,'-')[[1]])
          return((read.counts[1]+read.counts[2])/read.counts[3])
        }
      })))
    }),digits = 3),collapse = '|')
  }else{
    print(x)
    as.character(round(mean(unlist(lapply(x,function(tmp.af){
      if(tmp.af == 'NA'){return(0)}
      else{
        read.counts = as.numeric(str_split(tmp.af,'-')[[1]])
        return((read.counts[1]+read.counts[2])/read.counts[3])
      }
    }))),digits = 3))
  }
}

results.dir <- "/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_091819/"
condition.dir = 'results_combined_permissive'
results.dir = paste(results.dir,condition.dir,sep = '/')
ret.fusion.pt = fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/RET_fusion_patients.txt')


x = "/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_091819//results_combined_permissive/ret_004_table.csv"
do.call(rbind,lapply(list.files(results.dir,full.names = T),function(x){
  print(x)
  tmp.data = fread(x)[grepl('__RET|RET__',Hugo_Symbol)]
  fusion_expected = FALSE
  if(nrow(tmp.data) > 0){
    row.nums.called = unlist(lapply(1:nrow(tmp.data),function(y){
      if(any(tmp.data[y,(colnames(tmp.data)[grepl('duplex.called',colnames(tmp.data))]),with = F] == 'Called')) return(y)
    }))
    tmp.data = tmp.data[row.nums.called,]
  }
  if(any(grepl('RET',tmp.data$Hugo_Symbol))){
    return(data.frame(Patient_ID = str_extract(x,'ret_...'),
                      genetic_alterations = tmp.data$Hugo_Symbol,AF = unite(tmp.data[,(colnames(tmp.data)[grepl('___total',colnames(tmp.data))]),with = F],col = 'AFs',sep = '|')
                      # Fusion_Event = unique(tmp.data[,.(Fusion_Event = paste(Hugo_Symbol,Chromosome,Start_Position,End_Position,sep = '---'))]$Fusion_Event)
    ))
  }
})) %>% data.table() -> fusion.detection.table

# process original hugo_symbol column (sort two genes by name)
fusion.detection.table$genetic_alterations = unlist(lapply(fusion.detection.table$genetic_alterations,function(x){paste0(sort(str_split(x,'__')[[1]]),collapse = '-')}))
# collapsing AF for rows of the same events (i.e. reciprocal rearrangement) while perserving the sample level seaparation in AF
fusion.detection.table = fusion.detection.table[,.(collapsed_AF = collapse_AF(AFs)),.(Patient_ID,genetic_alterations)]
# change patient ID
fusion.detection.table$Patient_ID = gsub('ret_','EDD_ret_pt',fusion.detection.table$Patient_ID)

# clean up -- merging, renaming
merge(ret.fusion.pt,fusion.detection.table,by.x = c('Study_ID','genetic_alterations'), by.y = c('Patient_ID','genetic_alterations'),all = T) %>%
  transmute(Study_ID, genetic_alterations, DMP_ID, ACCESS_detected = ifelse(is.na(collapsed_AF),'NO','YES'),ACCESS_AF = collapsed_AF) %>% data.table() -> fusion.detection.table
# detection rate 
table(fusion.detection.table[!is.na(DMP_ID)]$ACCESS_detected)

write.csv(fusion.detection.table,paste0("/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/Fusion_Detection/Fusions_detected_",format(Sys.time(),'%m%d%y'),'.csv'),quote = F,row.names = F)

