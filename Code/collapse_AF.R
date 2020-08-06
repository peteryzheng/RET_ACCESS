library(data.table)
library(stringr)
library(tidyverse)

# x = c("32-117-190(0.7842)|53-0-959(0.0553)","NA|16-0-1035(0.0155)","63-0-954(0.066)|NA" )
# x = c('3-5-1300')
# x = c('NA|NA|NA|0-56-4183','NA|NA|NA|3-4-76')
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
          read.counts = as.numeric(str_split(gsub('\\(.*','',tmp.af),'-')[[1]])
          return((read.counts[1]+read.counts[2])/read.counts[3])
        }
      })))
    }),digits = 3),collapse = '|')
  }else{
    print(x)
    as.character(round(mean(unlist(lapply(x,function(tmp.af){
      if(tmp.af == 'NA'){return(0)}
      else{
        read.counts = as.numeric(str_split(gsub('\\(.*','',tmp.af),'-')[[1]])
        return((read.counts[1]+read.counts[2])/read.counts[3])
      }
    }))),digits = 3))
  }
}