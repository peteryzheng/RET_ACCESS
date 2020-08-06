library(data.table)
library(tidyverse)

table_to_maf = function(tmp.table,sample.table){
  # tmp.table = fillouts.dt
  # sample.table = sample.sheet
  # tmp.table = ret.006.table
  # sample.table = ret.006.sample.sheet
  # extract information for plasma and tumor
  tmp.table = data.table(tmp.table)
  lapply(sample.table[Sample_Type %in% c('duplex')]$Sample_Barcode,function(y){
    sample.call.status.colname = paste0(y,'___duplex.called')
    sample.af.colname = paste0(y,'___total')
    tmp.table[,eval(y) := paste0(get(sample.call.status.colname),' | ',get(sample.af.colname))]
  })
  lapply(sample.table[Sample_Type %in% c('Tumor')]$column.names,function(y){
    tmp.table[,eval(gsub('___.*.','',y)) := paste0(case_when(
      !is.na(get('DMP')) & get(paste0(sample.table[Sample_Type %in% c('duplex')]$Sample_Barcode[1],'___duplex.called')) != 'Not Covered' ~ 'Called',
      !is.na(get('DMP')) & get(paste0(sample.table[Sample_Type %in% c('duplex')]$Sample_Barcode[1],'___duplex.called')) == 'Not Covered' ~ 'Called (but not covered in ACCESS)',
      is.na(get('DMP')) & as.numeric(gsub('/.*','',get(y))) > 3 ~ 'Genotyped',
      TRUE ~ 'Not Called'
    ),' | ',get(y))]
  })
  processed.tmp.table = tmp.table[,!grep('___',colnames(tmp.table)),with = F] %>%
    # melting data frame by tumor samples
    melt(id.vars = c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification','HGVSp_Short',
                     'Reference_Allele','Tumor_Seq_Allele2','ExAC_AF','Hotspot','DMP','duplex_support_num','call_confidence','CH'),
         variable.name = "Tumor_Sample_Barcode", value.name = "call_info") %>%
    mutate(call_confidence =  gsub(' \\| ','',str_extract(call_info,'.*.\\| ')),call_info = gsub('.*.\\| ','',call_info)) %>% rowwise() %>%
    mutate(t_alt_count = ifelse(grepl('-[0-9]+-',call_info),
                                # SV parsing
                                sum(as.numeric(str_split(call_info,'-|\\(')[[1]][1:2])),
                                # SNV parsing
                                as.numeric(gsub(' |\\/.*.','',call_info))),
           t_total_count = ifelse(grepl('-[0-9]+-',call_info),
                                  # SV parsing
                                  as.numeric(str_split(call_info,'-|\\(')[[1]][3]),
                                  # SNV parsing
                                  as.numeric(gsub('.*.\\/|\\(.*.','',call_info)))) %>% data.table()
  return(processed.tmp.table)
}