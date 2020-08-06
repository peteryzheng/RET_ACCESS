library(data.table)
library(readxl)
library(dplyr)
library(tidyr)

old.bed.files <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Code/ACCESS_CNV/data/ACCESS_targets_coverage.bed')
new.bed.files <- read_xlsx('/ifs/work/bergerm1/zhengy1/RET/Original_file/RET probes list.xlsx') %>%
  mutate(tiling.name = paste0('panel_A:RET_',Chromosome,':',Start)) %>%
  select(Chromosome,Start,Stop,tiling.name) %>% data.table() 

colnames(old.bed.files) <- colnames(new.bed.files)

total.bed.files <- rbind(new.bed.files,old.bed.files) %>% arrange(Chromosome,Start) %>% data.table()

write.table(total.bed.files,
            '/ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0.sorted.RET.bed',
            sep = '\t',quote = F,row.names = F,col.names = F)

write.table(total.bed.files[!grepl('RET',tiling.name),],
            '/ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0.sorted.woRET.bed',
            sep = '\t',quote = F,row.names = F,col.names = F)


# gc section, normally don't need it --------------------------------------
system(paste0(
  'bedtools nuc -fi /ifs/depot/resources/dmp/data/pubdata/hg-fasta/VERSIONS/hg19/Homo_sapiens_assembly19.fasta ',
  ' -bed /ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0.sorted.RET.bed > ',
  ' /ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0.sorted.RET.txt'
)) 

ret.bed.gc.file <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0.sorted.RET.txt')  %>%
  setnames(c('#1_usercol','2_usercol','3_usercol','4_usercol','6_pct_gc'),
           c('Chrom','Start','End','Target','GC_150bp')) %>%
  select(Chrom,Start,End,Target,GC_150bp) %>% filter(grepl('RET',Target)) %>% data.table()

# example input for impact cna pipeline without RET information
bed.gc.model.file <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Code/ACCESS_CNV/data/ACCESS_targets_coverage.txt')

# get the extra probes no existing in the standard access pools
bed.gc.file <- merge(ret.bed.gc.file,bed.gc.model.file,by = c('Chrom','Start','End','Target','GC_150bp'),all.x = T) %>%
  filter(is.na(GeneExon)) %>%
  mutate(GeneExon = 'RET:exon:spiked-in',Gene = 'RET',Cyt = '10q11.21',Interval = paste0(Chrom,':',Start,'-',End))

bed.gc.ret.file <- rbind(bed.gc.model.file,bed.gc.file)
# all ret probes are pool A
bed.gc.ret.file[grepl('RET_intron',Target)]$Target <- gsub('panel_B','panel_A',bed.gc.ret.file[grepl('RET_intron',Target)]$Target)

write.table(bed.gc.ret.file,
            '/ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0.sorted.RET.merged.txt',
            sep = '\t',quote = F,row.names = F)

write.table(bed.gc.ret.file[!grepl('RET',Target),],
            '/ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0.sorted.woRET.merged.txt',
            sep = '\t',quote = F,row.names = F)



# bed file for mutation calling -------------------------------------------
old.bed.files <- fread('/ifs/work/bergerm1/Innovation/Resources/MSK-ACCESS-v1_0/MSK-ACCESS-v1_0_panelA_canonicaltargets.bed')
new.bed.files <- read_xlsx('/ifs/work/bergerm1/zhengy1/RET/Original_file/RET probes list.xlsx') %>%
  mutate(tiling.name = paste0('panel_A:RET_',Chromosome,':',Start),strand = '+') %>%
  select(Chromosome,Start,Stop,strand,tiling.name) %>% data.table() 

colnames(old.bed.files) <- colnames(new.bed.files)

total.bed.files <- rbind(new.bed.files,old.bed.files) %>% arrange(Chromosome,Start) %>% data.table()

write.table(total.bed.files,
            '/ifs/work/bergerm1/zhengy1/RET_all/Original_files/MSK-ACCESS-v1_0.canonicaltargets.RET.bed',
            sep = '\t',quote = F,row.names = F,col.names = F)


