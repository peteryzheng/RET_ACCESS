graph_AFs = function(table_path,covered,normals,simplified,output_dir){
  table_path = '/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_120919/results_combined/'
  print(x)
  tmp.data = fread(x)[grepl('__RET|RET__|^RET$',Hugo_Symbol)]
  fusion_expected = FALSE
  if(nrow(tmp.data) > 0){
    # filtered for "called" rows
    row.nums.called = unlist(lapply(1:nrow(tmp.data),function(y){
      # any RET event called at any timepoint is preserved for output
      if(any(tmp.data[y,(colnames(tmp.data)[grepl('duplex.called',colnames(tmp.data))]),with = F] == 'Called')) return(y)
    }))
    # filter for all called events
    tmp.data = tmp.data[row.nums.called,]
    
    if(nrow(tmp.data) > 0){
      tmp.data = tmp.data %>% rowwise() %>% mutate(Hugo_Symbol = ifelse(Hugo_Symbol == 'RET',gsub('p.','',HGVSp_Short),Hugo_Symbol)) %>% data.table()
      return(data.frame(Patient_ID = str_extract(x,'ret_...'),
                        genetic_alterations = tmp.data$Hugo_Symbol,
                        AF = unite(tmp.data[,(colnames(tmp.data)[grepl('___total',colnames(tmp.data))]),with = F],col = 'AFs',sep = '|'),
                        Status = unite(tmp.data[,(colnames(tmp.data)[grepl('___duplex.called',colnames(tmp.data))]),with = F],col = 'Status',sep = '|')
      ))
    }
}