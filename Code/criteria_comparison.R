library(data.table)

stringent.maf = fread('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_120919//results_stringent/all_sample_maf2maf_oncokb.maf')
permissive.maf = fread('/ifs/work/bergerm1/zhengy1/RET_all/Analysis_files/run_120919//results_permissive/all_sample_maf2maf_oncokb.maf')

# stringent called by patient
stringent.called = stringent.maf[call_confidence == 'Called' & grepl('-L...-|DA-ret',Tumor_Sample_Barcode),
                                 .(samples_called_in = .N,sample_names = paste0(Tumor_Sample_Barcode,collapse = '|')),
                                 .(Patient_ID,Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,
                                   Reference_Allele,Tumor_Seq_Allele2,HGVSp_Short,oncogenic,Hotspot)]

permissive.called = permissive.maf[call_confidence == 'Called' & grepl('-L...-|DA-ret',Tumor_Sample_Barcode),
                                   .(samples_called_in = .N,reads_called_by = paste0(t_alt_count,collapse = '|'),
                                     sample_names = paste0(Tumor_Sample_Barcode,collapse = '|')),
                                   .(Patient_ID,Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,
                                     Reference_Allele,Tumor_Seq_Allele2,HGVSp_Short,oncogenic,Hotspot)]

permissive.only = anti_join(permissive.called,stringent.called,by = c('Patient_ID','Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification',
                                                                      'Reference_Allele','Tumor_Seq_Allele2','HGVSp_Short','oncogenic','Hotspot')) %>% data.table()
# should be empty
stringent.only = anti_join(stringent.called,permissive.called,by = c('Patient_ID','Hugo_Symbol','Chromosome','Start_Position','End_Position','Variant_Classification',
                                                                     'Reference_Allele','Tumor_Seq_Allele2','HGVSp_Short','oncogenic','Hotspot')) %>% data.table()


View(permissive.only[oncogenic %in% c('Oncogenic','Likely Oncogenic')])
