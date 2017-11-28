data1 <- read.table("/hpc/users/nguyet26/psychgen/methods/extTADA/data/scz_data_Sep_2016_exac_noexac_oldDataBases_and_allSilentFCPK_andUK10KandFIN.txt",                    header = TRUE)

dirFile <- "~/psychgen/methods/extTADA/Re_annotate/PvalueGeneList/genesListAllPaper/"
tempData <- readLines(paste0(dirFile, "fmrp.txt"))

t1 <- pmatch(tempData, data1[, 1])
t1 <- t1[!is.na(t1)]
data2 <- data[t1,]


(sum(data1[, "case_silentCFPK_noexac_group1nca_3157_ncn_4672"])/3157)/(sum(data1[, "control_silentCFPK_noexac_group1nca_3157_ncn_4672"])/4672)
(sum(data2[, "case_silentCFPK_noexac_group1nca_3157_ncn_4672"])/3157)/(sum(data2[, "control_silentCFPK_noexac_group1nca_3157_ncn_4672"])/4672)
(sum(data2[, "case_lof_noexac_group1nca_3157_ncn_4672"])/3157)/(sum(data2[, "control_lof_noexac_group1nca_3157_ncn_4672"])/4672)
(sum(data2[, "case_damaging_noexac_group1nca_3157_ncn_4672"])/3157)/(sum(data2[, "control_damaging_noexac_group1nca_3157_ncn_4672"])/4672)


(sum(data2[, "case_lof_UK10K"])/1353)/(sum(data2[, "control_lof_UK10K"])/4769)
(sum(data2[, "case_damaging_UK10K"])/1353)/(sum(data2[, "control_damaging_UK10K"])/4769)
(sum(data2[, "case_missense_UK10K"])/1353)/(sum(data2[, "control_missense_UK10K"])/4769)

(sum(data2[, "case_lof_UK10K_FIN"])/392)/(sum(data2[, "control_lof_UK10K_FIN"])/2020)
(sum(data2[, "case_damaging_UK10K_FIN"])/392)/(sum(data2[, "control_damaging_UK10K_FIN"])/2020)
(sum(data2[, "case_missense_UK10K_FIN"])/392)/(sum(data2[, "control_missense_UK10K_FIN"])/2020)

#tempData <- head(tempData, 500)


