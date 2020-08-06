library(data.table)
library(stringr)

master.ref = fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_121519.csv')
results.by.pt = list.files('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_011520/results_combined_permissive/',full.names = T)
output.dir = paste0('/ifs/work/bergerm1/zhengy1/RET_all/SV_AF/SV_AF_',format(Sys.time(),'%m%d%y'))
dir.create(output.dir)

x = results.by.pt[4]
lapply(results.by.pt,function(x){
  tmp.table = fread(x)
  if(any(grepl('RET__|__RET',tmp.table$Hugo_Symbol))){
    pt.num = str_extract(x,'ret_...')
    dir.create(paste0(output.dir,'/',pt.num))
    # for every ret event
    y = grep('RET__|__RET',tmp.table$Hugo_Symbol)[1]
    lapply(grep('RET__|__RET',tmp.table$Hugo_Symbol),function(y){
      dir.create(paste0(output.dir,'/',pt.num,'/',tmp.table[y,1],'-',tmp.table[y,3],'-',tmp.table[y,4]))
      # for every ret positive sample
      detected.samples = colnames(tmp.table)[grep('Called',tmp.table[y,])]
      detected.duplex.samples = gsub('___.*.','',detected.samples[grep('___duplex',detected.samples)])
      print(detected.duplex.samples)
      lapply(detected.duplex.samples,function(z){
        dir.create(paste0(output.dir,'/',pt.num,'/',tmp.table[y,1],'-',tmp.table[y,3],'-',tmp.table[y,4],'/',z))
        write.table(data.frame(
          chr = c(gsub('__.*.','',tmp.table[y,2]),gsub('.*.__','',tmp.table[y,2])),
          start = as.numeric(c(tmp.table[y,3],tmp.table[y,4])),end = as.numeric(c(tmp.table[y,3]+1,tmp.table[y,4]+1))
        ),paste0(output.dir,'/',pt.num,'/',tmp.table[y,1],'-',tmp.table[y,3],'-',tmp.table[y,4],'/',z,'/input.bed'),
        sep = '\t',quote = F,row.names = F,col.names = F)
        system(paste0(
          'bsub -We 00:59 -cwd ',output.dir,'/',pt.num,'/',tmp.table[y,1],'-',tmp.table[y,3],'-',tmp.table[y,4],'/',z,
          ' -eo filtering_results.e -oo filtering_results.o ',
          ' bash /ifs/work/bergerm1/zhengy1/RET_all/Code/fragment_bam_filter.sh ',
          gsub('__aln.*.duplex','_filter_sc',gsub('duplex_bams//|duplex_bams/','standard_bams_filtered/',
                                                  master.ref[Tumor_Sample_Barcode.plasma == z]$BAM_path.plasma)),
          ' ',output.dir,'/',pt.num,'/',tmp.table[y,1],'-',tmp.table[y,3],'-',tmp.table[y,4],'/',z
        ))
      })
    })
  }
})

