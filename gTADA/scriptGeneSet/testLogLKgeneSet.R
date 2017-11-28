args <- commandArgs(TRUE)
geneSetName0 <- as.character(args[1])
OUTDIR <- as.character(args[2])

source("U_publishModel/extTADAaddGeneSets.R")
load("W_AUTgeneSet/AUTgeneSet.0..adjustRatioMut.1.1.nIteration1000.nThin.1.index.18_30_31_January.nGroupDN.3.nGroupCC.3.damaging.SameBeta.RData")

m <- stan_model(model_code = DNandCCextTADA)
#m <- stan_model(model_code = DNextTADA)

dirGeneSet <- "/hpc/users/nguyet26/psychgen/methods/extTADA/Re_annotate/PvalueGeneList/genesListAllPaper/"
allGeneSet <- dir(dirGeneSet, ".txt$")

geneSetName <- c(geneSetName0) #c("fmrp.txt", "constrained.txt", "essentialGenes.txt", "cahoy.txt", "psd95.txt", "mir137.txt", "synaptome.txt",
#                 "pLI09.txt", "nmdarc.txt", "chd8.hNSC+human_brain.txt", "zhu_gwas_eqtl.txt")
#geneSetName <- c(geneSetName,              allGeneSet[grep("^MGI", allGeneSet)])

gAll <- data[, 1]
geneSetM <- matrix(0, ncol = length(geneSetName), nrow = dim(data)[1])

for (ii in 1:length(geneSetName)){
    g1 <- readLines(paste0(dirGeneSet, geneSetName[ii]))
    g11 <- pmatch(g1, gAll)
    g11 <- g11[!is.na(g11)]
    if (!is.null(dim(geneSetM[g11, ])))
        geneSetM[g11, ][, ii] <- 1
    else
        geneSetM[g11, ] <- 1
    print(length(g11)) #table(g11)
    print(table(geneSetM[, ii]))

    }

mixDataKclasses <- list(K = 2, Kcc = Kcc, Kdn = Kdn,  hyperBetaMax = 100,
                        NCcc = dim(allCaseData)[2], NCdn = dim(allDenovoData)[2],
                        adjustHyperBeta = adjustHyperBeta0,
                        useDataCC = 1, useDataDN = 1,
                        NN = dim(allDenovoData)[1],
                        Ncontrol = array(allNcontrol),
                        Ntotal = array(allNcontrol + allNcase),
                        Ndn = array(c(rep(ntrio, 2))),
                        Ncase = array(allNcase),
                        dataDN = allDenovoData,
                        mutRate = allMutationData,

                        dataCCcase = allCaseData,
                        dataCCtotal = allCCData,
                        thetaH0 = array(allNcase/(allNcase + allNcontrol)),
                        upperPi0 = 0.5, lowerPi0 = 0,
                        lowerHyperGamma = 1, lowerBeta = 1, lowerGamma = 1,
                        hyperBetaDN0 = hyperBetaDN0,
                        hyperBetaCC0 = hyperBetaCC0,
                        betaPars = c(6.7771073, -1.7950864, -0.2168248),
                        Tgs = 0,
                        Ngs = dim(geneSetM)[2],
                        geneSet = data.frame(geneSetM))


ff <- optimizing(m, data = mixDataKclasses)

outLLK <- c(geneSetName0, round(ff$value, 4), round(head(ff$par, 4), 2))

write.table(t(outLLK), paste0(OUTDIR, "/LogLLK.AUT.", geneSetName0, ".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
#save.image(t(outLLK), paste0(OUTDIR, "/LogLLK.AUT.", geneSetName0, ".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
###Model with the gene set
geneSetName <- c("listMcRae2016.txt") #c("fmrp.txt", "constrained.txt", "essentialGenes.txt", "cahoy.txt", "psd95.txt", "mir137.txt", "synaptome.txt",
#                 "pLI09.txt", "nmdarc.txt", "chd8.hNSC+human_brain.txt", "zhu_gwas_eqtl.txt")
#geneSetName <- c(geneSetName,              allGeneSet[grep("^MGI", allGeneSet)])

gAll <- data[, 1]
geneSetM <- matrix(0, ncol = length(geneSetName), nrow = dim(data)[1])

for (ii in 1:length(geneSetName)){
    g1 <- readLines(paste0(dirGeneSet, geneSetName[ii]))
    g11 <- pmatch(g1, gAll)
    g11 <- g11[!is.na(g11)]
    if (!is.null(dim(geneSetM[g11, ])))
        geneSetM[g11, ][, ii] <- 1
    else
        geneSetM[g11, ] <- 1
    print(length(g11)) #table(g11)
    print(table(geneSetM[, ii]))

    }

mixDataKclasses <- list(K = 2, Kcc = Kcc, Kdn = Kdn,  hyperBetaMax = 100,
                        NCcc = dim(allCaseData)[2], NCdn = dim(allDenovoData)[2],
                        adjustHyperBeta = adjustHyperBeta0,
                        useDataCC = 1, useDataDN = 1,
                        NN = dim(allDenovoData)[1],
                        Ncontrol = array(allNcontrol),
                        Ntotal = array(allNcontrol + allNcase),
                        Ndn = array(c(rep(ntrio, 2))),
                        Ncase = array(allNcase),
                        dataDN = allDenovoData,
                        mutRate = allMutationData,

                        dataCCcase = allCaseData,
                        dataCCtotal = allCCData,
                        thetaH0 = array(allNcase/(allNcase + allNcontrol)),
                        upperPi0 = 0.5, lowerPi0 = 0,
                        lowerHyperGamma = 1, lowerBeta = 1, lowerGamma = 1,
                        hyperBetaDN0 = hyperBetaDN0,
                        hyperBetaCC0 = hyperBetaCC0,
                        betaPars = c(6.7771073, -1.7950864, -0.2168248),
                        Tgs =  dim(geneSetM)[2],
                        Ngs = dim(geneSetM)[2],
                        geneSet = data.frame(geneSetM))


ff1 <- optimizing(m, data = mixDataKclasses)

################################
###ADD for FMRP
geneSetName <- c("listMcRae2016.txt", "listLoFtolerantEXAC.txt") #geneSetName0) #c("fmrp.txt", "constrained.txt", "essentialGenes.txt", "cahoy.txt", "psd95.txt", "mir137.txt", "synaptome.txt",
#                 "pLI09.txt", "nmdarc.txt", "chd8.hNSC+human_brain.txt", "zhu_gwas_eqtl.txt")
#geneSetName <- c(geneSetName,              allGeneSet[grep("^MGI", allGeneSet)])

gAll <- data[, 1]
geneSetM <- matrix(0, ncol = length(geneSetName), nrow = dim(data)[1])

for (ii in 1:length(geneSetName)){
    g1 <- readLines(paste0(dirGeneSet, geneSetName[ii]))
    g11 <- pmatch(g1, gAll)
    g11 <- g11[!is.na(g11)]
    if (!is.null(dim(geneSetM[g11, ])))
        geneSetM[g11, ][, ii] <- 1
    else
        geneSetM[g11, ] <- 1
    print(length(g11)) #table(g11)
    print(table(geneSetM[, ii]))

    }

mixDataKclasses <- list(K = 2, Kcc = Kcc, Kdn = Kdn,  hyperBetaMax = 100,
                        NCcc = dim(allCaseData)[2], NCdn = dim(allDenovoData)[2],
                        adjustHyperBeta = adjustHyperBeta0,
                        useDataCC = 1, useDataDN = 1,
                        NN = dim(allDenovoData)[1],
                        Ncontrol = array(allNcontrol),
                        Ntotal = array(allNcontrol + allNcase),
                        Ndn = array(c(rep(ntrio, 2))),
                        Ncase = array(allNcase),
                        dataDN = allDenovoData,
                        mutRate = allMutationData,

                        dataCCcase = allCaseData,
                        dataCCtotal = allCCData,
                        thetaH0 = array(allNcase/(allNcase + allNcontrol)),
                        upperPi0 = 0.5, lowerPi0 = 0,
                        lowerHyperGamma = 1, lowerBeta = 1, lowerGamma = 1,
                        hyperBetaDN0 = hyperBetaDN0,
                        hyperBetaCC0 = hyperBetaCC0,
                        betaPars = c(6.7771073, -1.7950864, -0.2168248),
                                                               Tgs = dim(geneSetM)[2],
                                       Ngs = dim(geneSetM)[2],
                  geneSet = data.frame(geneSetM))


ff2 <- optimizing(m, data = mixDataKclasses)

geneSetName <- c("listMcRae2016.txt", "listLoFtolerantEXAC.txt", geneSetName0)
gAll <- data[, 1]
geneSetM <- matrix(0, ncol = length(geneSetName), nrow = dim(data)[1])

for (ii in 1:length(geneSetName)){
    g1 <- readLines(paste0(dirGeneSet, geneSetName[ii]))
    g11 <- pmatch(g1, gAll)
    g11 <- g11[!is.na(g11)]
    if (!is.null(dim(geneSetM[g11, ])))
        geneSetM[g11, ][, ii] <- 1
    else
        geneSetM[g11, ] <- 1
    print(length(g11)) #table(g11)
    print(table(geneSetM[, ii]))

    }

mixDataKclasses <- list(K = 2, Kcc = Kcc, Kdn = Kdn,  hyperBetaMax = 100,
                        NCcc = dim(allCaseData)[2], NCdn = dim(allDenovoData)[2],
                        adjustHyperBeta = adjustHyperBeta0,
                        useDataCC = 1, useDataDN = 1,
                        NN = dim(allDenovoData)[1],
                        Ncontrol = array(allNcontrol),
                        Ntotal = array(allNcontrol + allNcase),
                        Ndn = array(c(rep(ntrio, 2))),
                        Ncase = array(allNcase),
                        dataDN = allDenovoData,
                        mutRate = allMutationData,

                        dataCCcase = allCaseData,
                        dataCCtotal = allCCData,
                        thetaH0 = array(allNcase/(allNcase + allNcontrol)),
                        upperPi0 = 0.5, lowerPi0 = 0,
                        lowerHyperGamma = 1, lowerBeta = 1, lowerGamma = 1,
                        hyperBetaDN0 = hyperBetaDN0,
                        hyperBetaCC0 = hyperBetaCC0,
                        betaPars = c(6.7771073, -1.7950864, -0.2168248),
                                                               Tgs = dim(geneSetM)[2],
                                       Ngs = dim(geneSetM)[2],
                  geneSet = data.frame(geneSetM))


ff3 <- optimizing(m, data = mixDataKclasses)


outLLK1 <- c(geneSetName0, round(ff1$value, 4), round(head(ff1$par, 5), 2),
             round(ff1$value - (ff$value), 4),
             round(ff2$value - (ff1$value), 4),
             round(ff3$value - (ff2$value), 4)
             )


write.table(t(outLLK1), paste0(OUTDIR, "/LogLLK.AUT.", geneSetName0, ".txt.FMRP"), col.names = FALSE, row.names = FALSE, quote = FALSE)
