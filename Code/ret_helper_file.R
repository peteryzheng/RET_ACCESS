library(data.table)

dmp.maf = fread('/ifs/work/bergerm1/zhengy1/dmp/mskimpact/data_mutations_extended.txt')

# Make the helper maf file for genotyping some specific RET mutations
K666N.maf = fread('/ifs/work/bergerm1/zhengy1/RET_all/Original_files/helper_maf/ret_K666N.csv')
dmp.ret.maf = dmp.maf[Hugo_Symbol == 'RET' & HGVSp_Short %in% c('p.C634R','p.C634Y','p.K666N','p.M918T','p.V804M'),.SD[1,],HGVSp_Short]
dmp.ret.maf = dmp.ret.maf[,colnames(K666N.maf),with = FALSE]

helper.ret.maf = rbind(dmp.ret.maf,K666N.maf)

write.table(helper.ret.maf,'/ifs/work/bergerm1/zhengy1/RET_all/Original_files/helper_maf/ret_helper.maf',quote = F,row.names = F,sep = '\t')

system(paste0(
  'bsub -W 12:00  -R "rusage[mem=8]" -cwd /ifs/work/bergerm1/zhengy1/RET_all/Original_files/helper_maf/',
  ' -oo %J_fillout.o -eo %J_fillout.e',
  ' /opt/common/CentOS_7-dev/bin/perl /opt/common/CentOS_6-dev/vcf2maf/v1.6.17/maf2maf.pl ',
  ' --input-maf /ifs/work/bergerm1/zhengy1/RET_all/Original_files/helper_maf/ret_helper.maf',
  ' --ncbi-build GRCh37 --tum-vad-col t_alt_count --custom-enst /opt/common/CentOS_6-dev/vcf2maf/v1.6.17/data/isoform_overrides_at_mskcc ',
  ' --tmp-dir /ifs/work/bergerm1/zhengy1/RET_all/Original_files/helper_maf/',
  ' --vep-path /opt/common/CentOS_6-dev/vep/v95 --vep-forks 4 --nrm-rad-col n_ref_count --buffer-size 5000 ',
  ' --filter-vcf /opt/common/CentOS_6-dev/vep/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --nrm-vad-col n_alt_count --max-filter-ac 10 ',
  ' --vep-data /opt/common/CentOS_6-dev/vep/cache/ --tum-rad-col t_ref_count --nrm-depth-col n_depth ',
  ' --ref-fasta /ifs/depot/pi/resources/genomes/GRCh37/fasta/b37.fasta ',
  ' --output-maf /ifs/work/bergerm1/zhengy1/RET_all/Original_files/helper_maf/ret_helper_maf2maf.maf',
  ' --retain-cols Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,',
  'Tumor_Sample_UUID,Matched_Norm_Sample_UUID --species homo_sapiens --tum-depth-col t_depth'
))