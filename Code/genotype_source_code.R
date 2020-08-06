genotype_plasma_bams <- function(bam_path,sample_barcode,CMO_ID,maf_dir,sample_type,tmp_dir){
  dir.create(paste0(maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/'))
  tmp.maf.file <- list.files(paste0(maf_dir,'/',CMO_ID,'/'),pattern = 'temp.maf',full.names = T)
  # get frag count
  fillout.job.id <- system(paste0(
    'bsub -W 12:00  -R "rusage[mem=8]" -cwd ',maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',
    ' -oo %J_genotyping.o -eo %J_genotyping.e /home/hasanm/Innovation/software/maysun/GetBaseCountsMultiSample/GetBaseCountsMultiSample ',
    '--omaf --filter_duplicate 0 --thread 10 --maq 20 --fasta /ifs/depot/assemblies/H.sapiens/b37/b37.fasta ',
    '--maf ',tmp.maf.file,' --fragment_count 1 --output ',
    maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',sample_barcode,'___',sample_type,'_fillout.maf',
    ' --bam ',gsub('^.*.duplex_bams//|^.*.simplex_bams//|^.*.duplex_bams/|^.*.simplex_bams/|^.*.unfiltered_bams/|^.*.unfiltered_bams//|.bam$','',bam_path),
    ':',bam_path
  ),intern = T)
  fillout.job.id <- as.numeric(gsub('Job <|> is.*.$','',fillout.job.id))
  print(fillout.job.id)
  # cmo maf2maf
  # maf2mafid <- system(paste0(
  #   'bsub -W 12:00 -cwd ',maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',
  #   ' -oo %J_fillout.o -eo %J_fillout.e -w \"done(',fillout.job.id,')\"',' cmo_maf2maf --input-maf ',
  #   maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',sample_barcode,'___',sample_type,'_fillout.maf',
  #   ' --output-maf ',
  #   maf_dir,'/',CMO_ID,'/fillout/',sample_barcode,'___',sample_type,'_fillout_maf2maf.maf',
  #   ' --ref-fasta /ifs/depot/assemblies/H.sapiens/b37/b37.fasta ',
  #   ' --retain-cols Center,Verification_Status,Validation_Status,Mutation_Status,',
  #   'Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,',
  #   'Tumor_Sample_UUID,Matched_Norm_Sample_UUID,Caller,t_total_count_fragment,',
  #   't_ref_count_fragment,t_alt_count_fragment'
  # ),intern = T)
  # maf2maf
  maf2mafid <- system(paste0(
    'bsub -W 12:00  -R "rusage[mem=8]" -cwd ',maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',
    ' -oo %J_fillout.o -eo %J_fillout.e -w \"done(',fillout.job.id,')\"',
    ' /opt/common/CentOS_7-dev/bin/perl /opt/common/CentOS_6-dev/vcf2maf/v1.6.17/maf2maf.pl ',
    ' --input-maf ',maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',sample_barcode,'___',sample_type,'_fillout.maf',
    ' --ncbi-build GRCh37 --tum-vad-col t_alt_count --custom-enst /opt/common/CentOS_6-dev/vcf2maf/v1.6.17/data/isoform_overrides_at_mskcc ',
    ' --tmp-dir ',tempfile(pattern = "file", tmpdir = tmp_dir, fileext = ""),
    ' --vep-path /opt/common/CentOS_6-dev/vep/v95 --vep-forks 4 --nrm-rad-col n_ref_count --buffer-size 5000 ',
    ' --filter-vcf /opt/common/CentOS_6-dev/vep/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --nrm-vad-col n_alt_count --max-filter-ac 10 ',
    ' --vep-data /opt/common/CentOS_6-dev/vep/cache/ --tum-rad-col t_ref_count --nrm-depth-col n_depth ',
    ' --ref-fasta /ifs/depot/pi/resources/genomes/GRCh37/fasta/b37.fasta ',
    ' --output-maf ',maf_dir,'/',CMO_ID,'/fillout/',sample_barcode,'___',sample_type,'_fillout_maf2maf.maf',
    ' --retain-cols Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,',
    'Tumor_Sample_UUID,Matched_Norm_Sample_UUID,Caller,t_total_count_fragment,t_ref_count_fragment,t_alt_count_fragment --species homo_sapiens --tum-depth-col t_depth'
  ),intern = T)
  maf2mafid <- as.numeric(gsub('Job <|> is.*.$','',maf2mafid))
  return(maf2mafid)
}

genotype_impact_bams <- function(bam_path,sample_barcode,CMO_ID,maf_dir,sample_type,tmp_dir){
  dir.create(paste0(maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/'))
  fillout.job.id <- system(paste0(
    'bsub -W 12:00  -R "rusage[mem=8]" -cwd ',maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',
    ' -oo %J_genotyping.o -eo %J_genotyping.e  cmo_fillout -m ',maf_dir,'/',CMO_ID,'/',CMO_ID,'_all_unique_calls.maf',
    ' -g GRCh37 -f 1 -b ',bam_path,' -p ',maf_dir,'/',CMO_ID,'/portal/',sample_barcode,'___',sample_type,'_fillout.portal.maf',
    ' -o ',maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',sample_barcode,'___',sample_type,'_fillout.maf',
    ' -v 1.2.1'
  ),intern = T)
  fillout.job.id <- as.numeric(gsub('Job <|> is.*.$','',fillout.job.id))
  print(fillout.job.id)
  # cmo maf2maf
  # maf2mafid <- system(paste0(
  #   'bsub -W 12:00 -cwd ',maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',
  #   ' -oo %J_fillout.o -eo %J_fillout.e -w \"done(',fillout.job.id,')\"',' cmo_maf2maf --input-maf ',
  #   maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',sample_barcode,'___',sample_type,'_fillout.maf',
  #   ' --output-maf ',
  #   maf_dir,'/',CMO_ID,'/fillout/',sample_barcode,'___',sample_type,'_fillout_maf2maf.maf',
  #   ' --ref-fasta /ifs/depot/assemblies/H.sapiens/b37/b37.fasta ',
  #   ' --retain-cols Center,Verification_Status,Validation_Status,Mutation_Status,',
  #   'Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,',
  #   'Tumor_Sample_UUID,Matched_Norm_Sample_UUID,Caller'
  # ),intern = T)
  # maf2maf
  maf2mafid <- system(paste0(
    'bsub -W 12:00  -R "rusage[mem=8]"  -cwd ',maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',
    ' -oo %J_fillout.o -eo %J_fillout.e -w \"done(',fillout.job.id,')\"',
    ' /opt/common/CentOS_7-dev/bin/perl /opt/common/CentOS_6-dev/vcf2maf/v1.6.17/maf2maf.pl ',
    ' --input-maf ',maf_dir,'/',CMO_ID,'/',sample_barcode,'_',sample_type,'/',sample_barcode,'___',sample_type,'_fillout.maf',
    ' --ncbi-build GRCh37 --tum-vad-col t_alt_count --custom-enst /opt/common/CentOS_6-dev/vcf2maf/v1.6.17/data/isoform_overrides_at_mskcc ',
    ' --tmp-dir ',tempfile(pattern = "file", tmpdir = tmp_dir, fileext = ""),
    ' --vep-path /opt/common/CentOS_6-dev/vep/v95 --vep-forks 4 --nrm-rad-col n_ref_count --buffer-size 5000 ',
    ' --filter-vcf /opt/common/CentOS_6-dev/vep/cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz --nrm-vad-col n_alt_count --max-filter-ac 10 ',
    ' --vep-data /opt/common/CentOS_6-dev/vep/cache/ --tum-rad-col t_ref_count --nrm-depth-col n_depth ',
    ' --ref-fasta /ifs/depot/pi/resources/genomes/GRCh37/fasta/b37.fasta ',
    ' --output-maf ',maf_dir,'/',CMO_ID,'/fillout/',sample_barcode,'___',sample_type,'_fillout_maf2maf.maf',
    ' --retain-cols Center,Verification_Status,Validation_Status,Mutation_Status,Sequencing_Phase,Sequence_Source,Validation_Method,Score,BAM_file,Sequencer,',
    'Tumor_Sample_UUID,Matched_Norm_Sample_UUID,Caller,t_total_count_fragment,t_ref_count_fragment,t_alt_count_fragment --species homo_sapiens --tum-depth-col t_depth'
  ),intern = T)
  maf2mafid <- as.numeric(gsub('Job <|> is.*.$','',maf2mafid))
  return(maf2mafid)
}
