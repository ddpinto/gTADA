args <- commandArgs(TRUE)
nUpdate <- as.numeric(args[1])
nIteration <-   as.numeric(args[2])
inDex <- as.character(args[3])
nThin <- as.numeric(args[4])
nChain <- as.numeric(args[5])
annotationType <- as.character(args[6])
annotationType2 <- as.character(args[7])
outDir <- as.character(args[9])
lowerHyperGamma <- as.numeric(args[10])
annotationType3 <- as.character(args[11])
resultDir <- as.character(args[12])
nGroupDN <- as.numeric(args[13])
nGroupCC <- as.numeric(args[14])
nCore <-  nChain ###cores = chains
 adjustHyperBeta <- as.numeric(args[15])
misType <- as.character(args[16])
sigmaPrior <- as.numeric(args[17])
geneSetName0 <- as.character(args[18])

.libPaths("~/InstallSoftware/Rlibs/")

library("rstan")
citation("rstan")
library("coda")

data <- read.table("/hpc/users/nguyet26/psychgen/methods/extTADA/data/scz_data_Sep_2016_exac_noexac_oldDataBases_and_allSilentFCPK.txt",
                   header = TRUE)
data <- read.table("/hpc/users/nguyet26/psychgen/methods/extTADA/data/scz_data_Sep_2016_exac_noexac_oldDataBases_and_allSilentFCPK_andUK10KandFIN_fullMclust.txt",
                   header = TRUE)
data <- read.table("/hpc/users/nguyet26/psychgen/methods/extTADA/data/scz_data_Nov_2016_silent_Guipponi.txt", header = TRUE)


adjustRatioMut <- sum(sum(data[, 'dn_silent']))/(2*1024*sum(data[, 'mut_silent']))

dirGeneSet <- "/hpc/users/nguyet26/psychgen/methods/extTADA/Re_annotate/PvalueGeneList/genesListAllPaper/"
allGeneSet <- dir(dirGeneSet, ".txt$")

geneSetName <- c("fmrp.txt", "constrained.txt", "essentialGenes.txt",
                 "cahoy.txt", "psd95.txt", "mir137.txt", "synaptome.txt")
               #,                  allGeneSet[grep("^MGI", allGeneSet)])
g0 <- read.table(paste0("Test_LLK/TestlogLLKxem.AddGenes.lof_SCZ.damaging_SCZ.txt2"))
geneSetName <- as.character(g0[, 2])[-1]

geneSetName <- geneSetName0

print(geneSetName)

gAll <- data[, 1]
geneSetM <- matrix(0, ncol = length(geneSetName), nrow = dim(data)[1])
for (ii in 1:length(geneSetName)){
    g1 <- readLines(paste0(dirGeneSet, geneSetName[ii]))
    g11 <- pmatch(g1, gAll)
    g11 <- g11[!is.na(g11)]
    geneSetM[g11, ii] <- 1
    print(length(g11)) #table(g11)
    print(table(geneSetM[, ii]))

    }



t1 <- as.numeric(data$mut_lof)
t1 <- t1[t1 > 0]
data[data$mut_lof == 0, ]$mut_lof <- min(t1)/10
#data <- data[data$mut_missense >0, ]
t2 <- as.numeric(data$mut_missense)
t2 <- t2[t2 > 0]
data[data$mut_missense == 0, ]$mut_missense <- min(t2)/10
dim(data)

#data <- data[data$mut_missense >0, ]
t2 <- as.numeric(data$mut_damaging)
t2 <- t2[t2 > 0]
data[data$mut_damaging == 0, ]$mut_damaging <- min(t2)/10

##Silent
t3 <- as.numeric(data$mut_silentCFPK)
t3 <- t3[t3 > 0]
data[data$mut_silentCFPK == 0, ]$mut_silentCFPK <- min(t3)/10


dim(data)


ntrio <- 617 #Menachem
ntrio <- ntrio + 14 #Girard
ntrio <- ntrio + 105 #Gulsuner
ntrio <- ntrio + 57 #McCarthy
ntrio <- ntrio + 231 #Xu
ntrio <- ntrio + 53

ncase = 4253 #3064 #4954
ncontrol <- 5882 #4536 #6239
nTrans <- 617 ##Transmitted data are from Menachem publication



ncase1 = 3157 ##Swed
ncontrol1 = 4672 ##Swed

ncase2=1091
ncontrol2=1193

ncase3 = 1353
ncontrol3 = 4769

ncase4 = 392
ncontrol4 = 2020

#ncase <- nTrans # + ncase
#ncontrol <- nTrans # + ncontrol

ncaseA <- ncase 
ncontrolA <- ncontrol 


data <- data[data$mut_lof != 0, ]
dim(data)

#counts <- cbind(sCountCC2[, 1], sCountCC2[, 2:6])



#adjustHyperBeta0 <- 0
adjustHyperBeta0 <- adjustHyperBeta

betaPars = c(6.7771073, -1.7950864, -0.2168248)

#exp(betaPars[1]*hyperGammaMeanCC[j1j]^(betaPars[2]) + betaPars[3])
data0 <- data
fmrp <- readLines("../PvalueGeneList/genesListAllPaper/fmrp.txt")
t11 <- pmatch(fmrp, data[, 1])
t11 <- t11[!is.na(t11)]
#data <- data0[t11,]
caseData1 <- cbind(data[, 'case_damaging_ncase3157_ncontrol4672_NoExAC_Mclust']
                   +                    data[, 'case_LoF_ncase3157_ncontrol4672_NoExAC_Mclust']                   )
controlData1 <- cbind(data[, 'control_damaging_ncase3157_ncontrol4672_NoExAC_Mclust'] + 
                      data[, 'control_LoF_ncase3157_ncontrol4672_NoExAC_Mclust'])

caseData2 <- cbind(
                   data[, "case_damaging_ncase1091_ncontrol1193_NoExAC_Mclust"] + 
                   data[, "case_LoF_ncase1091_ncontrol1193_NoExAC_Mclust"])
controlData2 <- cbind(
                      data[, "control_damaging_ncase1091_ncontrol1193_NoExAC_Mclust"] + 
                      data[, "control_LoF_ncase1091_ncontrol1193_NoExAC_Mclust"])

caseData3 <- cbind(data[, 'case_damaging_UK10K_UK_NoExAC'] + data[, 'case_lof_UK10K_UK_NoExAC'])
controlData3 <- cbind(data[, 'control_damaging_UK10K_UK_NoExAC'] + data[, 'control_lof_UK10K_UK_NoExAC'])

caseData4 <- cbind(       data[, 'case_damaging_UK10K_FIN_NoExAC'] + data[, 'case_lof_UK10K_FIN_NoExAC'])

controlData4 <- cbind(
                      data[, 'control_damaging_UK10K_FIN_NoExAC'] + data[, 'control_lof_UK10K_FIN_NoExAC'])



allCaseData <- data.frame(caseData1, caseData2,  caseData3)
allControlData <-  data.frame(controlData1, controlData2, controlData3)

allNcase <- c(rep(ncase1, dim(caseData1)[2]), rep(ncase2, dim(caseData2)[2]), rep(ncase3, dim(caseData3)[2]))
allNcontrol <- c(rep(ncontrol1, dim(controlData1)[2]), rep(ncontrol2, dim(controlData2)[2]), rep(ncontrol3, dim(controlData3)[2]))

(apply(allCaseData, 2, sum)/allNcase)/(apply(allControlData, 2, sum)/allNcontrol)

misType
allDenovoData <- data[, c('dn_silentCFPK', paste0('dn_', misType), 'dn_lof')]
allMutationData <- data[, c('mut_silentCFPK', paste0('mut_', misType), 'mut_lof')]

adjustRatioMut <- 1
allMutationData <- adjustRatioMut*allMutationData

y.case.lof = allCaseData[, 1]
NPdn = 1
nGroupCC = 3 #dim(allControlData)[2]
nGroupDN = 3
NPcc = 3; NPdn = 1
Kdn = 3; Kcc = 2 #dim(allControlData)[2]

allCCData <- allCaseData + allControlData
NCcc = dim(allCaseData)[2]

hyperBetaDN0 <- rep(5, 3)
hyperBetaCC0 <- rep(5, NCcc)
nIteration1 <- nIteration
resultDir1 <- resultDir
nThin1 <- nThin
adjustHyperBeta01 <- adjustHyperBeta0
 hyperBetaDN0 <- rep(4, dim(allDenovoData)[2])
 hyperBetaCC0 <- rep(4, dim(allCaseData)[2])

 
adjustHyperBeta <- adjustHyperBeta0
#load("W_SCZgeneSet3/GuipponiNoFIN.DNandCC.adjustHyperBeta.1..adjustRatioMut.1.3.nIteration5000.nThin.5.index.11_34_17_April.nGroupDN.3.nGroupCC.3.damaging.SameBeta.RData")


mixDataKclasses <- list(K = 2, Kcc = Kcc, Kdn = Kdn,  hyperBetaMax = 100,
                        NCcc = dim(allCaseData)[2], NCdn = dim(allDenovoData)[2],
                        adjustHyperBeta = adjustHyperBeta0,
                        useDataCC = 1, useDataDN = 1,
                        NN = dim(allDenovoData)[1],
                        Ncontrol = array(allNcontrol),
                        Ntotal = array(allNcontrol + allNcase),
                        Ndn = array(c(rep(ntrio, 3))),
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
                                       Ngs = dim(geneSetM)[2],
                  geneSet = data.frame(geneSetM),
                                 Tgs = dim(geneSetM)[2],
                  sigmaPrior = sigmaPrior)


head(mixDataKclasses$dataCCcase, 2)
head(mixDataKclasses$dataCCtotal, 2)
head(mixDataKclasses$mutRate, 2)

mixDataKclasses$hyperBetaDN

apply(mixDataKclasses$dataDN, 2, sum)
apply(mixDataKclasses$mutRate, 2, sum)
apply(mixDataKclasses$dataDN, 2, sum)
apply(mixDataKclasses$dataCCcase, 2, sum)
apply(mixDataKclasses$dataCCtotal, 2, sum)
mixDataKclasses$thetaH0
mixDataKclasses$hyperBetaDN0

nCore <- nChain
 #load("W_SCZgui12/GuipponiBothNoFINcombined.DNandCC.adjustHyperBeta.1..adjustRatioMut.1.3.nIteration5000.nThin.5.index.22_46_22_November.nGroupDN.3.nGroupCC.3.damaging.SameBeta.RData"); nChain = 1; nCore  = 1; nThin = 1; nIteration = 1000
#detach("package:rstan", unload = TRUE); detach("package:rstan", unload = TRUE)
#source("U_publishModel/extTADAforMultiPopsNotDividePops.R")
source("U_publishModel/extTADAaddGeneSets.R")
#library("rstan")

if (nIteration < 1000)
    nThin = 1

testIntegratedModel <- stan(model_code = DNandCCextTADA,
                        data = mixDataKclasses, iter = nIteration,
                                                chains = nChain, cores = nCore,
                                                thin = nThin,
                                                pars = c('alpha0', 'hyperGammaMeanDN', 'hyperBetaDN',
                                                         'hyperGammaMeanCC', 'hyperBetaCC'))

inDex <- format(Sys.time(), "%H_%M_%d_%B")
save.image(paste0(resultDir, "/GuipponiNoFIN.DNandCC.adjustHyperBeta.", adjustHyperBeta, ".",
                  ".adjustRatioMut.",  round(adjustRatioMut, 2), ".",
                         nChain, ".nIteration", nIteration, ".nThin.",
        nThin, ".index.", inDex,  ".nGroupDN.", nGroupDN,
        ".nGroupCC.", nGroupCC, ".", misType,
        ".SameBeta.RData"))

############Estimate pars
source("TADA/TADA.R")
fileR <- dir("../extTADA/script", ".R")
for (ii in fileR)
    source(paste0("../extTADA/script/", ii))

fileR <- dir("../gTADA/script", ".R")
for (ii in fileR)
    source(paste0("../gTADA/script/", ii))
source("../gTADA/script/estimatePars.R")

parTemp <- c('alpha0[1]', 'alpha0[2]',
             'hyperGammaMeanDN[1]', 'hyperGammaMeanDN[2]', 'hyperGammaMeanDN[3]',
                               'hyperBetaDN[1]', 'hyperBetaDN[2]', 'hyperBetaDN[3]',
                               'hyperGammaMeanCC[1]', 'hyperGammaMeanCC[2]',
                               'hyperBetaCC[1]', 'hyperBetaCC[2]')

pars0 <- estimatePars(pars = parTemp,
                      mcmcResult = testIntegratedModel)

pars0 <- apply(pars0, 2, function(x) round(x, 3))
pars1 <- cbind(rep("SCZ", dim(pars0)[1]),
               rep(geneSetName0, dim(pars0)[1]),
               rownames(pars0),
               pars0)



#load("W_SCZgeneSetEachGeneSet/GuipponiNoFIN.DNandCC.adjustHyperBeta.1..adjustRatioMut.1.3.nIteration2000.nThin.2.index.09_38_19_April.nGroupDN.3.nGroupCC.3.damaging.SameBeta.RData")
write.table(pars1, paste0(resultDir, "/EstimatedPars.adjustBeta.", adjustHyperBeta, ".",
                         nChain, ".nI", nIteration, ".nThin.",
        nThin, ".index.", inDex, 
        ".nGroupCC.", nGroupCC, ".adjustRatio.", round(adjustRatioMut, 2), ".",
        geneSetName0,
        ".sigmaPrior.", sigmaPrior, ".", geneSetName0,
        ".txt"), row.names = FALSE, quote = FALSE, col.names= FALSE)

#################
#################

b1 <- as.data.frame(testIntegratedModel)
bMCMC <- mcmc(b1)
bHPD <-  HPDinterval(bMCMC)
medianHPD <- NULL





for (ii in 1:dim(bHPD)[1]){
    t1 <- b1[, ii]

    t2 <- t1[(t1>=bHPD[ii, 1]) & (t1<=bHPD[ii, 2])]

    medianHPD[ii] <- mean(t2)

     }

names(medianHPD) <- colnames(b1)

t3 <- c(medianHPD)

write.table(t(t3), paste0(resultDir, "/MCMC.", ntrio,
                         nChain, ".nIteration", nIteration, ".nThin.",
        nThin, ".index.", inDex, ".", annotationType, ".", annotationType2, ".nGroup.", nGroupDN,
        ".nGroupCC.", nGroupCC,  ".", sigmaPrior, ".",
        ".notChangeDNgammaPrior.txt"),
         col.names = FALSE, quote= FALSE, row.names = FALSE)

save.image(paste0(resultDir, "/GuipponiNoFIN.DNandCC.adjustHyperBeta.", adjustHyperBeta, ".",
                  ".adjustRatioMut.",  round(adjustRatioMut, 2), ".",
                         nChain, ".nIteration", nIteration, ".nThin.",
        nThin, ".index.", inDex,  ".nGroupDN.", nGroupDN,
        ".nGroupCC.", nGroupCC, ".", misType, ".", sigmaPrior, ".",
        ".SameBeta.RData"))















#save.image(paste0(resultDir, "/test_rstan_2016.Pio.nChain.",
 #                        nChain, ".nIteration", nIteration, ".nThin.",
  #      nThin, ".index.", inDex, ".", annotationType, ".", annotationType2, ".nGroupDN.", nGroupDN,
   #     ".nGroupCC.", nGroupCC,
    #    ".SameBeta.RData"))
