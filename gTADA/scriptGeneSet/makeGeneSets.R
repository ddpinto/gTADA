
data <- read.table("/hpc/users/nguyet26/psychgen/methods/extTADA/Re_annotate/DenovoData/epi_id_dd_scz_data_Sep_2016_addNewID.txt",
                   header = TRUE, as.is = TRUE, sep = ",")
dirGeneSet <- "../PvalueGeneList/genesListAllPaper/"
geneSetName <- c("fmrp.txt", "constrained.txt", "essentialGenes.txt", "cahoy.txt")

gAll <- data[, 1]
geneSetM <- matrix(0, ncol = length(geneSetName), nrow = dim(data)[1])
for (ii in 1:length(geneSetName)){
    g1 <- readLines(paste0(dirGeneSet, geneSetName[ii]))
    g11 <- pmatch(g1, gAll)
    g11 <- g11[!is.na(g11)]
    geneSetM[g11, ][, ii] <- 1
    print(length(g11)) #table(g11)
    print(table(geneSetM[, ii]))

    }

