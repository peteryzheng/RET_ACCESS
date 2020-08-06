library(data.table)

cna.dir <- '/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/cna_053019/'

probe.data <- fread(paste0(cna.dir,'/EDD_ret_copynumber_segclusp.probes.txt'))
gene.data <- fread(paste0(cna.dir,'/EDD_ret_copynumber_segclusp.genes.txt'))



# probe+normal artifact ---------------------------------------------------

# total times each probe+normal appeared
artifact.df <- merge(probe.data,probe.data[abs(fc) > 1.25,.(freq = .N),.(norm_used_auto,region)],by = c('norm_used_auto','region')) %>%
  # Get total times normal was used
  merge(unique(probe.data[,.(sample,norm_used_auto)])[,.(norm_used_freq = .N),.(norm_used_auto)],by = 'norm_used_auto') %>% data.table()
  
blacklist.probes <- unique(artifact.df[freq >= norm_used_freq* 0.75 & freq != 1,.(norm_used_auto,region,freq,norm_used_freq)])

blacklist.probes[,.(.N),.(region)]
blacklist.probes[,.(.N),.(region)][N == 6]
blacklist.probes[,.(.N),.(norm_used_auto)]

# met amp patients --------------------------------------------------------

gene.data[region == 'MET']



# pten loss 898XKN (pt 47) max AF RB1/TP53 ~45% ---------------------------

pten.898xkn <- probe.data[sample == 'C.898XKN.L001.d_mean_cvg' & chr == 7 & fc < -2]
merge(pten.898xkn,blacklist.probes,by = c('norm_used_auto','region'),all.x = T)

# other CNAs
cna.898xkn <- gene.data[sample == 'C.898XKN.L001.d_mean_cvg' & chr %in% c(1,5,9,16,17,18,20,'X')]
# Looks like CNA on 1,5,9, 16p, 17p (TP53),18.20

# 17p stat. sign. probes
probe.data[sample == 'C.898XKN.L001.d_mean_cvg' & chr %in% c(17) & p.adj < 0.05]

# focal loss at intron region of ret before exon 12 -- before kinase domain
probe.data[sample == 'C.898XKN.L001.d_mean_cvg' & chr %in% c(10) & p.adj < 0.30]


# 8F9XD2 patient 6  -------------------------------------------------------
# Driver mutation -- VTCN1, gained KRAS

# plasma timepoint1 
# max af 60% ATM, doublet RET mutation, SMARCA4 (40%)
gene.data[sample == 'DA.ret.006.pl.T01_IGO_05500_FF_20_mean_cvg' & region == 'ATM'] # allelic imbalance?
# loss @ 3,4,5q, 8,9,11,13
probe.data[sample == 'DA.ret.006.pl.T01_IGO_05500_FF_20_mean_cvg' & chr %in% c(3,4,5,8,9,11,13)] 
gene.data[sample == 'DA.ret.006.pl.T01_IGO_05500_FF_20_mean_cvg' & chr %in% c(3,4,5,8,9,11,13)] 

# plasma timepoint2 
# max af  33% ATM, ~40% doublet RET mutation, SMARCA4 23%
# loss @ 1,3,4,5q, 8,9,11,13
probe.data[sample == 'DA.ret.006.pl.T02_IGO_05500_FF_21_mean_cvg' & chr %in% c(1,3,4,5,8,9,11,13)] 
gene.data[sample == 'DA.ret.006.pl.T02_IGO_05500_FF_21_mean_cvg' & chr %in% c(1,3,4,5,8,9,11,13)] 

# 003 max AF 25% ATM, RET 804 disappeared, 918 stayed at 20%, SMARCA4 15%
# gain of 7, 14
probe.data[sample == 'C.8F9XD2.L003.d_mean_cvg']
