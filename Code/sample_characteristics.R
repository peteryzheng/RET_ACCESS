library(data.table)
library(tidyverse)

master.ref <- fread('/ifs/work/bergerm1/zhengy1/RET_all/Sample_mapping/master_ref_121519.csv')

# number of samples sequenced
nrow(master.ref)

# number of unique patinets
length(unique(master.ref$`CMO patient ID`))

# breakdown of samples per patient
master.ref[,.N,EDD_Patient_ID]
table(master.ref[,.N,EDD_Patient_ID]$N)

# number of patients with paired DMP IMPACT
length(unique(master.ref$DMP_ID))

DMP.key <- fread('/ifs/dmprequest/12-245/key.txt')
DMP.key[grepl(paste0(unique(master.ref$DMP_ID),collapse = '|'),V1)]$V1


# DMP.key[grepl(paste0(unique(master.ref$DMP_ID),collapse = '|'),V1)]$V1[grep('XS',DMP.key[grepl(paste0(unique(master.ref$DMP_ID),collapse = '|'),V1)]$V1)]

# detection rate overview
detection.rate.data = fread('/ifs/work/bergerm1/zhengy1/RET_all/For_Ezra/RET_Detection/RET_indicated_detected_120219_review.csv')
# fusions
table(detection.rate.data[Alt_Type == 'Fusion']$genetic_alterations)
table(detection.rate.data[Alt_Type == 'Fusion']$Baseline_Alt)
detection.rate.data[Alt_Type == 'Fusion',.(Yes = length(which(Baseline_Alt == 'Yes')),
                                           No = length(which(Baseline_Alt == 'No')),
                                           Reviewed = length(which(Baseline_Alt == 'Yes_review'))),genetic_alterations]
nrow(detection.rate.data[Alt_Type == 'Point'])
table(detection.rate.data[Alt_Type == 'Point']$genetic_alterations)
table(detection.rate.data[Alt_Type == 'Point']$Baseline_Alt)
detection.rate.data[Alt_Type == 'Point',.(Yes = length(which(Baseline_Alt == 'Yes')),
                                           No = length(which(Baseline_Alt == 'No')),
                                           Reviewed = length(which(Baseline_Alt == 'Yes_review'))),genetic_alterations]
