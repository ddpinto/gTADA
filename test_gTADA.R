library("rstan")
fileR <- dir("../script", ".R$")
for (ii in fileR)
    source(paste0("../script/", ii))

data <- read.table("../data/data_mut_DD.csv", header = TRUE, as.is = TRUE)
allDNData <- data[, paste0("dn_", c("damaging", "lof"), "_DD")]
allMutData <- data[,paste0("mut_", c("damaging", "lof"))]
head(data.frame(allMutData, allDNData))

mcmcDD <- gTADA(modelName = DNextTADA,
                  dataDN = allDNData, mutRate = allMutData,
                            Ndn = rep(4293, 2),
                nIteration = 1000,
                geneSet = NULL
                )

##Make a random gene set
rgeneSet <- data.frame(sample(0:1, dim(allDNData)[1], replace = TRUE))
mcmcDD1 <- gTADA(modelName = DNextTADA,
                  dataDN = allDNData, mutRate = allMutData,
                            Ndn = rep(4293, 2),
                nIteration = 1000,
                geneSet = rgeneSet
                )


geneNDD <- as.character(read.table("../data/NDDgenesList_DataForAllDisorders.txt", header = TRUE)[, 1])

geneSet1 <- as.numeric(data[, 1] %in% geneNDD)
table(geneSet1)


mcmcDD2 <- gTADA(modelName = DNextTADA,
                  dataDN = allDNData, mutRate = allMutData,
                            Ndn = rep(4293, 2),
                nIteration = 1000,
                geneSet = data.frame(geneSet1)
                )


listDD <- read.table("../scriptGeneSet/Test_LLK/TestlogLLKxem.AddGenes.lof_DD.damaging_DD.txt2", as.is = TRUE)

geneSet2 <- NULL

#geneSet1 <- as.numeric(data[, 1] %in% geneNDD)
for (j in 2:(dim(listDD)[1] - 1)) {
    t1 <- readLines(paste0("~/Documents/Github/TestR/HBproject/GeneSetAll/genesListAllPaper/", listDD[j, 2]))
    geneSet2 <- cbind(geneSet2, as.numeric(data[, 1] %in% t1))
}

mcmcDD3 <- gTADA(modelName = DNextTADA,
                  dataDN = allDNData, mutRate = allMutData,
                            Ndn = rep(4293, 2),
                nIteration = 1000,
                geneSet = data.frame(geneSet2)
                )


gTemp <- prcomp(geneSet2)
geneSet3 <- gTemp$x[, 1:4]
mcmcDD4 <- gTADA(modelName = DNextTADA,
                  dataDN = allDNData, mutRate = allMutData,
                            Ndn = rep(4293, 2),
                nIteration = 1000,
                geneSet = data.frame(geneSet3)
                )


#options(repr.plot.width=5, repr.plot.height=5)
stan_trace(mcmcDD)

stan_trace(mcmcDD)
stan_trace(mcmcDD2)
stan_trace(mcmcDD3)

pars0 <- estimatePars(pars = c('alpha0[1]',
                               'hyperGammaMeanDN[1]', 'hyperGammaMeanDN[2]',
                               'hyperBetaDN[1]', 'hyperBetaDN[2]'),
                     mcmcResult = mcmcDD)

pars1 <- estimatePars(pars = c('alpha0[1]', 'alpha0[2]',
                               'hyperGammaMeanDN[1]', 'hyperGammaMeanDN[2]',
                               'hyperBetaDN[1]', 'hyperBetaDN[2]'),
                     mcmcResult = mcmcDD1)

pars2 <- estimatePars(pars = c('alpha0[1]', 'alpha0[2]',
                               'hyperGammaMeanDN[1]', 'hyperGammaMeanDN[2]',
                               'hyperBetaDN[1]', 'hyperBetaDN[2]'),
                     mcmcResult = mcmcDD2)

allPars = c('alpha0[1]', 'alpha0[2]',
                               'hyperGammaMeanDN[1]', 'hyperGammaMeanDN[2]',
                               'hyperBetaDN[1]', 'hyperBetaDN[2]')
#aT <- as.data.frame(mcmcDD3)
aTN <- colnames(aT)
aTN[grep("alpha0|hyperGammaMean|hyperBetaDN", aTN)]

pars3 <- estimatePars(pars = aTN[grep("alpha0|hyperGammaMean|hyperBetaDN", aTN)],
                     mcmcResult = mcmcDD3)

aTN1 <- colnames(as.data.frame(mcmcDD4))
pars4 <- estimatePars(pars = aTN1[grep("alpha0|hyperGammaMean|hyperBetaDN", aTN1)],
                     mcmcResult = mcmcDD4)

#save.image("mcmcDDdata.RData")

ntrioDD = 4293

parsFDR <- list(gammaMeanDN = pars0[, 1][2:3],
                betaDN = pars0[, 1][4:5],
                                        #pi0 = pars0[, 1][1],
                alpha0 = pars0[1, 1],
            nfamily = rep(ntrioDD, 2))

geneName = data[, 1]
dataFDR <- calculateFDR(pars = parsFDR, geneSet = NULL,
                       dnData = allDNData, mutData = allMutData,
                       geneName = geneName)

parsFDR1 <- list(gammaMeanDN = pars1[, 1][3:4],
                betaDN = pars1[, 1][5:6],
                                        #pi0 = pars0[, 1][1],
                alpha0 = pars1[, 1][1:2],
            nfamily = rep(ntrioDD, 2))

geneName = data[, 1]
dataFDR1 <- calculateFDR(pars = parsFDR1, geneSet = rgeneSet,
                       dnData = allDNData, mutData = allMutData,
                       geneName = geneName)

###
parsFDR2 <- list(gammaMeanDN = pars2[, 1][3:4],
                betaDN = pars2[, 1][5:6],
                                        #pi0 = pars0[, 1][1],
                alpha0 = pars2[, 1][1:2],
            nfamily = rep(ntrioDD, 2))

parsFDR3 <- list(gammaMeanDN = pars3[, 1][18:19],
                betaDN = pars3[, 1][20:21],
                                        #pi0 = pars0[, 1][1],
                alpha0 = pars3[, 1][1:17],
            nfamily = rep(ntrioDD, 2))

parsFDR4 <- list(gammaMeanDN = pars4[, 1][6:7],
                betaDN = pars4[, 1][8:9],
                                        #pi0 = pars0[, 1][1],
                alpha0 = pars3[, 1][1:5],
            nfamily = rep(ntrioDD, 2))

geneName = data[, 1]
dataFDR2 <- calculateFDR(pars = parsFDR2, geneSet = geneSet1,
                       dnData = allDNData, mutData = allMutData,
                       geneName = geneName)

dataFDR3 <- calculateFDR(pars = parsFDR3, geneSet = geneSet2,
                       dnData = allDNData, mutData = allMutData,
                       geneName = geneName)

dataFDR4 <- calculateFDR(pars = parsFDR4, geneSet = geneSet3,
                       dnData = allDNData, mutData = allMutData,
                       geneName = geneName)


pThreshold = 0.8
dim(dataFDR[dataFDR$PP > pThreshold, ])
dim(dataFDR1[dataFDR1$PP > pThreshold, ])
dim(dataFDR2[dataFDR2$PP > pThreshold, ])
dim(dataFDR3[dataFDR3$PP > pThreshold, ])
dim(dataFDR4[dataFDR4$PP > pThreshold, ])

jC = 7
b11 <- merge(dataFDR[, c(1, jC)], dataFDR3[, c(1, jC)], by = "geneName")
dim(b11[(b11[, 2] > 0.8) & (b11[, 2] > 0.8),])

plot(b11[, 2], b11[, 3])
abline(a = 0, b = 1)
