library(data.table)

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_111919.csv')
ret.detected.file = fread('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/RET_Detection/RET_indicated_detected_120219_review.csv')

# number of patients with multiple time points
nrow(master.ref[,.N,`CMO patient ID`][N > 1])

# number of access run
nrow(master.ref)

# how many patients
length(unique(master.ref$`CMO patient ID`))


unique.patient.data = unique(master.ref[,.(.N,study_ID = paste0(unique(gsub('EDD-ret-pt|DA-ret-','EDD_ret_pt',
                                        str_extract(`Investigator sample ID.plasma`,'EDD-ret-pt...|EDD_ret_pt...|DA-ret-...')
)),collapse = ', ')),.(`CMO patient ID`,DMP_ID)])

unique.patient.data[!study_ID %in% ret.detected.file$Study_ID]

# detection rate of all known targets
table(ret.detected.file$Baseline_Alt)
1-length(which(ret.detected.file$Baseline_Alt == 'No'))/nrow(ret.detected.file)

# detection rate of all known fusion targets
table(ret.detected.file[Alt_Type == 'Fusion']$Baseline_Alt)
1-length(which(ret.detected.file[Alt_Type == 'Fusion']$Baseline_Alt == 'No'))/nrow(ret.detected.file[Alt_Type == 'Fusion'])

# detection rate of all known SNV targets
table(ret.detected.file[Alt_Type == 'Point']$Baseline_Alt)
1-length(which(ret.detected.file[Alt_Type == 'Point']$Baseline_Alt == 'No'))/nrow(ret.detected.file[Alt_Type == 'Point'])
