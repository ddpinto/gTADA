 .libPaths("~/InstallSoftware/Rlibs/")
library("rstan")
dirFile <- "W_EPIgeneSet3/"
dirFile <- "W_EPIgeneSet/"
dirFile <- "W_EPIgeneSet_R/"
dirFile <- "W_SCZgeneSet1/"
dirFile <- "W_AUTgeneSet/"
dirFile <- "W_AUTgeneSet1/"
dirFile <- "W_EPIgeneSetN/"
dirFile <- "W_EPIgeneSetN1/"
dirFile <- "W_EPIgeneSetN3/"
dirFile <- "W_EPIgeneSetN2/"
dirFile <- "W_EPIgeneSetN2_1/"
dirFile <- "W_EPIgeneSetN2_ElasticN/"
dirFile <- "W_EPIgeneSetN2_ElasticN2/"
dirFile <- "W_EPIgeneSetN2_ElasticN1/"
dirFile <- "W_EPIgeneSetN2_ElasticNcolinear/"
dirFile <- "W_EPIgeneSetN2_BayesianPrior/"
dirFile <- "W_EPIgeneSetN2_BayesianPriorLaplace/"

dirFile <- "W_AddGeneSet1/"
dirFile <- "W_AddGeneSetKF3/"
dirFile <- "W_AddGeneSetKF5/"
dirFile <- "W_SCZgeneSet3/"
#dirFile <- "W_AUTgeneSet3/"
dirFile <- "W_AddGeneSetKF5_LaplacePrior/"
dirFile <- "W_AddGeneSetKF6_NormalPrior/"
 fileN <- dir(dirFile, "RData")
fileN

#fileN <- fileN[grep("2000|5000", fileN)]
#fileN <- fileN[grep("1000", fileN)]
#fileN <- fileN[grep("Feb_08_2017", fileN)]
#fileN <- fileN[grep("5000", fileN)]

nameAn <- "damaging" #"missense" #
fileN <- fileN[grep(nameAn, fileN)]
fileN
oTheta0 <- NULL
outLog <- NULL

llk1 <- function(pars){
        hyperGammaMean <- pars[1:nCol]
        hyperBeta <- pars[(nCol + 1):(2*nCol)]
        alpha0 <- pars[(2*nCol+1):(2*nCol + nGeneSet)]
        alphaIntercept <- tail(pars, 1)
        f0 <- f1 <- 1
        for (j in 1:nCol){
            f0 <- f0*dpois(dD[, j], 2*nFamily[j]*muRate[, j])
            f1<- f1*dnbinom(dD[, j], hyperGammaMean[j]*hyperBeta[j], hyperBeta[j]/(hyperBeta[j] + 2*nFamily[j]*muRate[, j]))
            }
        pi1 <- apply(geneSet, 1, function(x){
            #pi1 = exp(geneSet)/(1 + exp(geneSet))
            exp0 <- exp(alphaIntercept + sum(alpha0*x))
            return(exp0/(1 + exp0))})
        tllk <- sum(log(f1*pi1 + (1 - pi1)*f0))
##lasso
#        tllk <- tllk - pT*sum(abs(alpha0))
#ridge 
#        tllk <- tllk - pT*sum((alpha0^2))

        return(-tllk)
    }

#pdf("~/www/Test/xemD.pdf")

#fileN <- fileN[grep("sigmaPrior.3|sigmaPrior.2|sigmaPrior.10", fileN)]
fileN <- fileN[grep("sigmaPrior.3.2", fileN)]
fileN <- fileN[grep("sigmaPrior.3.2kF.5", fileN)]
fileN <- fileN[grep("10000", fileN)]


 pdf(paste0("~/www/Test/xemD", nameAn, ".", gsub("/", "", dirFile), ".pdf"))
NK <- length(fileN)

for (ja in 1:NK){

#    par(mfrow = c(3, 2))
   load(paste0(dirFile, fileN[ja]))
b1 <- as.data.frame(testIntegratedModel)
t1 <- grep('lambda0', colnames(b1))
#    if (length(t1) == 1)
 #       oTheta0[ja] <- median(b1[, 'lambda0[1]'])
  #  else
   #     oTheta0[ja] <- median(b1[, 'lambda0[1]'])
    #print(oTheta0[ja])
    #plot(density(b1[, 'lambda0[1]']))
#14892
    if (dim(b1)[1] > 0){
        apply(b1[1:3, 1:8], 2, quantile)

        tempLog <- c(annotationType, annotationType2,
                     sigmaPrior, round(median(b1[, 'lp__']), 3))

         nameAlpha <- colnames(b1)
    nameAlpha <- nameAlpha[grep("alpha0", nameAlpha)]
#tail(b1[, 'pi0[14892]'])
#tail(b1[, 'pi0[1]'])
#tail(b1[, 'pi0[2]'])

        geneSetNameN <- geneSetName
    balpha0 <- t(apply(b1[, nameAlpha], 2, function(x) quantile(x, c(0.025, 0.5, 0.975))))
    rownames(balpha0) <- c( "alpha0", geneSetNameN)

        if (is.element("alpha0Intercept", nameAlpha))
                rownames(balpha0) <- c( geneSetNameN, "alpha0")

        ###########Calculate logLK
        exp0 <- balpha0[1, 2]
        for (iG1 in 1:dim(geneSetM)[2]){
            exp0 <- exp0 + balpha0[iG1 + 1, 2]*geneSetM[, iG1]
                    }
        pE0 <- exp(exp0)/(1 + exp(exp0))
        b1E <- apply(b1, 2, median)
        gM <- b1E[grep("hyperGammaMeanDN", names(b1E))]
        bM <- b1E[grep("hyperBetaDN", names(b1E))]
        nFamily = mixDataKclasses$Ndn
        muRate <- mixDataKclasses$mutRate
        geneSet = geneSetM
        nCol = length(gM)
        nGeneSet = dim(geneSetM)[2]
        dD <- mixDataKclasses$dataDN

        llkE <- llk1(pars = c(gM, bM, balpha0[, 2][-1], balpha0[1,2]))

        outLog[[ja]] <- c(tempLog, llkE)

        ############################
    balpha0 <- balpha0[order(balpha0[, 2]),]
    head(balpha0)
    tail(balpha0)

    if (dim(balpha0)[1] > 55)
        balpha0 <- balpha0[((balpha0[, 1] > 0) & (balpha0[, 3] > 0)) |
                           ((balpha0[, 1] < 0) & (balpha0[, 3] < 0)),, drop = FALSE]


#plot(density(b1[, "fmrp.txt"]))
plot(c(-6, 6), c(0, dim(balpha0)[1] + 1), xlim = range(as.numeric(balpha0)),
     main = paste0(annotationType, ' & ', annotationType2, '\nsigmaPrior: ', sigmaPrior),
     xlab = '', ylab = '', col = 'white')

for (ik in 1:dim(balpha0)[1]){
    xTemp <- as.numeric(balpha0[ik, ]) #b1[, paste0('alpha0[', ik, ']')]
    lines(xTemp[c(1,3)], rep(ik, 2), lwd = 2, col = 'blue')
    text(xTemp[2], ik+0.2, paste0(rownames(balpha0)[ik]), cex = 0.6)
    text(xTemp[2], ik, 'x', col = ifelse((xTemp[1]>0) | (xTemp[3]<0), 'red', 'blue'))
    abline(v = 0, col = 'blue', cex = 2)

}
    } else {
        file.remove(paste0(dirFile, fileN[ja]))
    }
}
dev.off()

######################
outLog1 <- do.call(rbind, outLog)
outLog1 <- outLog1[order(outLog1[, 1]), ]
outLog1
colnames(outLog1) <- paste0("V", 1:dim(outLog1)[2])

outLog2 <- data.frame(outLog1)
outLog2 <- outLog2[outLog2[, 1] == "lof_AST_TADA",]

outLog2[, 3] <- as.numeric(as.character(outLog2[, 3]))
outLog2[, 5] <- round(as.numeric(as.character(outLog2[, 5])), 2)

p1 <- ggplot(outLog2, aes(V3, V5, col = as.factor(V1))) + geom_line()

pdf("~/www/Test/testSigmaPrior.pdf")
p1
dev.off()

b22 <- apply(b1, 2, median)

tail(b22)

###Test CR
xC <- read.table("result.CR10.RealData.txt", as.is = TRUE)
xC <- read.table("result.CR10.RealData.LaplacePrior.txt", as.is = TRUE)
xC <- read.table("result.CR10.RealData.normal.txt", as.is = TRUE)
#xC <- read.table("result.CR10.RealData.txt2", as.is = TRUE)
library("dplyr")
xC1 <- xC %>% group_by(V1, V2, V3) %>% summarise(mK = mean(V5), mL = mean(V4))
xC2 <- as.data.frame(xC1)


xC2[grep("damaging_CHD", xC2[, 2]),]
xC2[grep("missense_CHD", xC2[, 2]),]

pName = as.factor(paste0(xC1[, 1], xC1[, 2]))
p1 <- ggplot(xC1, aes(V3, mK, col = V2)) + geom_line()
p1 <- ggplot(xC1, aes(V3, mK)) + geom_line()
p1 <- p1 + facet_grid(. ~ V2)

pdf("~/www/Test/CVplot.pdf")
p1
dev.off()

xC2 <- data.frame(xC1)
xC2[, c(4, 5)] <- apply(xC2[, c(4, 5)], 2, function(x) round(x, 0))

pdf("~/www/Test/CVplot.pdf")
for (ii in names(table(xC2[, 2]))){
    temp1 <- xC2[xC2[, 2] == ii,]
    plot(temp1[, 3], temp1[, 5], main = ii, xlab = 'sigmaPrior', ylab = 'logLLK', type = 'l')
    t10 <-  temp1[which(temp1[, 5] == max(temp1[, 5])), 3][1]
    t11 <-  temp1[is.element(temp1[, 3], t10),][, 5][1]
    abline(v = t10, col = 'red')
    text(t10, t11 -2,
         paste0(t10, "-", t11), srt = 45, lwd = 0.7)
    
}
dev.off()



#load("W_AddGeneSetKF1/DN.adjustBeta.1.1.nI1000.nThin.1.index.10_17_Mar_03_2017.lof_CHD_allN1213.missense_CHD_allN1213.nGroupDN.2.nGroupCC.2.adjustRatio.1.0.05.sigmaPrior.4kF.8.RData")
str(testIntegratedModel@sim$samples)

annotationType
annotationType2
bA <- apply(b1, 2, median)
bAN <- names(bA)

hyperGammaMean <- bAN[grep("hyperGammaMeanDN", bAN)]

g1 <- apply(b1[, hyperGammaMean], 2, median)

source("TADA/TADA.R")
gammaDN <- bA[hyperGammaMean]
betaDN <- bA[grep("hyperBetaDN", bAN)]
pName <- bAN[grep("pi0", bAN)]
pi0A <- bA[pName]

bfDN <- matrix(0, nrow = dim(data)[1], ncol = length(gammaDN))
for (i2 in 1:length(gammaDN))
    bfDN[, i2] <- bayes.factor.denovo(x = mixDataKclasses$dataDN[, i2],
                                      N = mixDataKclasses$Ndn[i2],
                                      gamma.mean = gammaDN[i2], beta = betaDN[i2],
                                      mu = mixDataKclasses$mutRate[, i2])

bfAll <- apply(bfDN, 1, prod)
pp <- pi0A*bfAll/((1 - pi0A) + pi0A*bfAll)
zpp <- ifelse(pp > 0.5, 1, 0)
table(zpp)
dataOut <- data.frame(data[, 1],  mixDataKclasses$dataDN, pp)
tCC <- mixDataKclasses$dataCCtotal - mixDataKclasses$dataCCcase
x11 <- mixDataKclasses$Ntotal*round(mixDataKclasses$dataCCcase/mixDataKclasses$Ncase - tCC/(mixDataKclasses$Ntotal - mixDataKclasses$Ncase), 3); colnames(x11) = paste0("V", 1:dim(x11)[2])
dataOut <- cbind(dataOut, x11)

fdr0 = 0.1 #0.1 #0.05
dataOut[pp > (1 - fdr0),]
dim(dataOut[pp > (1 - fdr0),])


###Use optimizing
#m <- stan_model(model_code = DNandCCextTADA);
#m1 <- stan_model(model_code = DNextTADA);

dirGeneSet <- "/hpc/users/nguyet26/psychgen/methods/extTADA/Re_annotate/PvalueGeneList/genesListAllPaper/"
allGeneSet <- dir(dirGeneSet, ".txt$")

geneSetName <- c("fmrp.txt", "constrained.txt", "essentialGenes.txt",
                 "cahoy.txt", "psd95.txt", "mir137.txt", "synaptome.txt",
                 "pLI09.txt", "nmdarc.txt", "chd8.hNSC+human_brain.txt", "zhu_gwas_eqtl.txt")
#geneSetName <- c(geneSetName,              allGeneSet[grep("^MGI", allGeneSet)])

geneSetName <- c("constrained.txt","kirov_psych.ed.set.ARC.txt","chd8.human_brain.txt","GO:0070925.txt")
gAll <- data[, 1]
geneSetM <- matrix(0, ncol = length(geneSetName), nrow = dim(data)[1])

for (ii in 1:length(geneSetName)){
    g1 <- readLines(paste0(dirGeneSet, geneSetName[ii]))
    g11 <- pmatch(g1, gAll)
    g11 <- g11[!is.na(g11)]
    if (dim(geneSetM)[2] == 1)
        geneSetM[g11, ][ii] <- 1
    else
        geneSetM[g11, ][, ii] <- 1
    print(length(g11)) #table(g11)
    print(table(geneSetM[, ii]))

    }

Ngs = dim(geneSetM)[2]
Ngs = 0
#load("/hpc/users/nguyet26/psychgen/methods/extTADA/Re_annotate/SimulatedData/TestDNandCC0025and02_4293/TestArray.100001.gammaMeanDN.2.5.gammaMeanDenovo.15.gammaMeanDN2.3.5.gammaMeanDenovo2.5.pi0.0.05.rhocC.0.0002161572.rhoCC2.0.000123642505820651.lowerHyperGamma.1.casecontrolAndDenovo.1classes.RData")

#allDenovoData <- (mixDataKclasses$dataDN); allMutationData <- (mixDataKclasses$mutRate); rGene <- data.frame(data[, 1], sDenovoCC2[[2]][, 2])

#load("W_EPIgeneSetN2_BayesianPrior/DN.adjustBeta.1.1.nI5000.nThin.1.index.14_20_Feb_12_2017.lof_SCZ.damaging_SCZ.nGroupDN.2.nGroupCC.2.adjustRatio.1.0.05.sigmaPrior.4.951.RData")

str(sDenovoCC2)
mixDataKclasses1 <- list(K = 2, Kcc = Kcc, Kdn = Kdn,  hyperBetaMax = 100,
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
                        Tgs = 0, ##Test Gene Set
                        Ngs = dim(geneSetM)[2],
                  geneSet = data.frame(geneSetM))
mixDataKclasses1$Ndn <- mixDataKclasses$Ndn

mixDataKclasses$Tgs = dim(geneSetM)[2]
mixDataKclasses$Ngs = dim(geneSetM)[2]
mixDataKclasses$geneSet <- geneSetM

source("U_publishModel/extTADAaddGeneSets.R")
#source("U_publishModel/extTADAaddGeneSetsPenalizedAlphaUseSTAN.R")

ff <- optimizing(m, data = mixDataKclasses, iter = 10000, tol_obj = 10^-14, tol_grad = 10^-14)
#m1 <- stan_model(model_code = DNextTADA);
#ff <- optimizing(m1, data = mixDataKclasses, iter = 10000, tol_obj = 10^-18, tol_grad = 10^-18)
ff$value -ff2
ff$value -ff1

bA <- ff$par
bAa <- bA[grep("alpha", names(bA))]
bAN <- names(bA)
hyperGammaMean <- bAN[grep("hyperGammaMeanDN", bAN)]

source("TADA/TADA.R")
gammaDN <- bA[hyperGammaMean]
betaDN <- bA[grep("hyperBetaDN", bAN)]
pName <- bAN[grep("pi0", bAN)]
pi0A <- bA[pName]
mean(pi0A)


bfDN <- matrix(0, nrow = dim(data)[1], ncol = length(gammaDN))
for (i2 in 1:length(gammaDN))
    bfDN[, i2] <- bayes.factor.denovo(x = mixDataKclasses$dataDN[, i2],
                                      N = mixDataKclasses$Ndn[i2],
                                      gamma.mean = gammaDN[i2], beta = betaDN[i2],
                                      mu = mixDataKclasses$mutRate[, i2])
bfAll <- apply(bfDN, 1, prod)
pp <- pi0A*bfAll/((1 - pi0A) + pi0A*bfAll)
zpp <- ifelse(pp > 0.5, 1, 0)
table(zpp)
dataOut <- data.frame(data[, 1],  mixDataKclasses$dataDN, pp)

fdr0 = 0.01 #0.1 #0.05
dataOut[pp > (1 - fdr0),]
dim(dataOut[pp > (1 - fdr0),])

d11 <- dataOut[pp > (1 - fdr0),]
sum(d11[, 2] + d11[, 3])

cGene <- as.character(d11[, 1])
table(rGene[pmatch(cGene, rGene[, 1]),][, 2])
pGene <- readLines("ListGeneFDR001fromTADAonlyDenovo.txt")
length(pGene)
length(intersect(cGene, pGene))

setdiff(pGene, cGene)

dataOut[pmatch(setdiff(pGene, cGene), dataOut[, 1]),]
################
rAUT <- read.table("W_logLK_AUT/allGeneSet.LLK.txt")
rAUT <- read.table("W_logLK_AUT/allGeneSet.LLK.txtFMRP")
rAUT <- rAUT[order(-rAUT[, 2]),]
head(rAUT)
##############

rG <- rAUT[rAUT[, dim(rAUT)[2]] > 2,]
rG <- rG[grep("listMcRae2016.txt|denovo|FDR|1.essentila.txt", rG[, 1], ignore.case = TRUE, invert = TRUE),]
rG[, 1]

head(dataOut[order(-dataOut[, dim(dataOut)[2]]),])

for (i2 in c(14892, 16407)){
    
    
}

p1 <- NULL
y <- data$dn_lof_SCZ + data$dn_damaging_SCZ

y <- data$dn_lof + data$dn_damaging
for (i in 1:dim(geneSetM)[2]){
    z <- geneSetM[, i]
    a1 <- glm(z ~ y)
    p1[i] <- (summary(a1)$coefficients)[2, 4]
}

geneSetName[p1 < 0.01]


dirFile <- "W_EPIgeneSetN2_BayesianPrior/"
fileN <- dir(dirFile, "RData")
fileN
fileN <- fileN[grep("2000|5000", fileN)]

fileN
oTheta0 <- NULL
outLog <- NULL

#pdf("~/www/Test/xemD.pdf")
N2 <- length(fileN)
for (ja in 1:10){

load(paste0(dirFile, fileN[ja]))
b1 <- as.data.frame(testIntegratedModel)
        tempLog <- c(annotationType, annotationType2,
                     sigmaPrior, round(median(b1[, 'lp__']), 3))

    write.table(t(tempLog), paste0(resultDir, "/LogLK.adjustBeta.", adjustHyperBeta, ".",
                         nChain, ".nI", nIteration, ".nThin.",
        nThin, ".index.", inDex, ".", annotationType, ".", annotationType2, ".nGroupDN.", nGroupDN,
        ".nGroupCC.", nGroupCC, ".adjustRatio.", round(adjustRatioMut, 2), ".", beta22, ".sigmaPrior.", sigmaPrior, 
        ".text"), row.names = FALSE, quote = FALSE, col.names= FALSE)

        outLog[[ja]] <- c(tempLog, llkE)

}


for ii in $(seq 1 1 244)
do
bsub -q premium -P acc_psychgen -R "rusage[mem=10000]" -W 0:10 -oo W_EPIgeneSetN2_BayesianPrior/OUT \
  R CMD BATCH --no-save --no-restore "--args ${ii}" writeLogLK.R W_EPIgeneSetN2_BayesianPrior/RoutOUT

done

for ii in $(seq 1 1 244)
do
bsub -q premium -P acc_psychgen -R "rusage[mem=10000]" -W 0:10 -oo W_EPIgeneSetN2_BayesianPrior/OUT \
  R CMD BATCH --no-save --no-restore "--args ${ii}" writeLogLK.R W_EPIgeneSetN2_BayesianPrior/RoutOUT

done

###Use penalized likelihood
for ii in $(seq 0.00001 0.05 10)
do
index=$(date|sed 's/ //g'|sed 's/:/_/g')
bsub -q premium -P acc_psychgen -R "rusage[mem=3000]" -W 0:50 -oo Test_LLK/OUT \
  R CMD BATCH --no-save --no-restore "--args ${ii} ${index}" pLikelihood.R Test_LLK/R.${ii}.RoutOUT

done



rFile=increaseLLK1_addColinear.R
NN=146
rFile=increaseLLK1.R
for ii in $(seq 1 1 146)

do
index=$(date|sed 's/ //g'|sed 's/:/_/g')
bsub -q premium -P acc_psychgen -M 20000 -R "rusage[mem=3000]" -W 3:50 -oo Test_LLK/OUT \
  R CMD BATCH --no-save --no-restore "--args ${ii} ${index}" ./${rFile} Test_LLK/R.${ii}.RoutOUT

done


############
rFile=increaseLLK1_removeLowLK1.R
for ii in "DD" "ASD" "ID" "SCZ"
do
index=$(date|sed 's/ //g'|sed 's/:/_/g')
bsub -q premium -P acc_psychgen -M 20000 -R "rusage[mem=3000]" -W 22:50 -oo Test_LLK/OUT \
  R CMD BATCH --no-save --no-restore "--args ${ii} ${index}" ./${rFile} Test_LLK/R.${ii}.RoutOUT

done

##################33
x1 <- read.table("testLLK.txt", header = FALSE, as.is = TRUE)
x12 <- unique(x1[, 1:2])

pdf("~/www/Test/TestLLK.pdf")
for (i in 1:dim(x12)[1]){
    x1T <- x1[(x1[, 1] == x12[i, 1]) & (x1[, 2] == x12[i, 2]),]
    plot(x1T[, 3], x1T[, 4] - min(x1T[, 4]),
         main = paste0(x12[i, 1], "-", x12[i, 2]))
    abline(v = x1T[which(max(x1T[, 4]) == x1T[, 4]), 3], col = 'red')
    print(i)
}
dev.off()
    
x12 <- split(x1, cbind(x1[, 1], x1[, 2]))


fileN <- "xemLLK.DD.txt"
#fileN <- "xemLLK.DD.txt1"
#fileN <- "xemLLK.SCZ.txt"
#fileN <- "xemLLK.ID.txt"
#fileN <- "xemLLK.AST.txt"
x1 <- read.table(paste0("Test_LLK/", fileN))
x1 <- x1[order(x1[, 3], decreasing = FALSE),]
x1[, 1] <- -x1[, 1]
x1[, 2] <- -x1[, 2]
tT <- NULL
pdf(paste0("~/www/Test/", fileN, ".pdf"))
plot(x1[, 3], x1[, 2], xlab = 'Add Pars', ylab = 'LLK')
abline(h = x1[1, 1], col = 'red')
for (ii in 2:dim(x1)[1]){
    if ((x1[ii, 2] - x1[ii-1, 2]) < 0){
        abline(v = x1[ii-1, 3], col = 'blue')
        text(x1[ii-1, 3], x1[ii-1, 2], x1[ii-1, 2], col = 'red', cex = 0.4)
        tT <- c(tT, x1[ii - 1, 3])
    }}
    
        
dev.off()

print(tT)


load("../scriptXB/W_SCZthesameHyperGamma/GuipponiNoFIN.DNandCC.adjustHyperBeta.0..adjustRatioMut.1.3.nIteration10000.nThin.10.index.18_17_29_March.nGroupDN.3.nGroupCC.3.damaging.SameBeta.RData")

dirFile <- "W_AUTgeneSet1/"
fileN <- dir(dirFile, "RData")
ja = 1
load(paste0(dirFile, fileN[ja]))

#########################
fileR <- dir("../extTADA/script", ".R")
for (ii in fileR)
    source(paste0("../extTADA/script/", ii))

parTemp <- c('alpha', 'hyperGammaMeanDN[1]', 'hyperGammaMeanDN[2]',
                               'hyperGammaMeanCC[1]', 'hyperGammaMeanCC[2]',
                               'hyperBetaDN[1]', 'hyperBetaDN[2]',
                               'hyperBetaCC[1]', 'hyperBetaCC[2]')

pars0 <- estimatePars(pars = NULL,
                      mcmcResult = testIntegratedModel)

caseData <- mixDataKclasses$dataCCcase
controlData <- mixDataKclasses$dataCCtotal - mixDataKclasses$dataCCcase
dnData <- mixDataKclasses$dataDN
mutData <- mixDataKclasses$mutRate
geneName <- data[, 1]
geneSet <- mixDataKclasses$geneSet
ncase = mixDataKclasses$Ncase
ncontrol = mixDataKclasses$Ncontrol
parsPP <- list(gammaMeanDN = pars0[grep("hyperGammaMeanDN", rownames(pars0)), 1],
               betaDN = pars0[grep("hyperBetaDN", rownames(pars0)), 1],
               nfamily = rep(ntrio, dim(dnData)[2]),
               gammaMeanCC = pars0[grep("hyperGammaMeanCC", rownames(pars0)), 1],
               betaCC = pars0[grep("hyperGammaMeanCC", rownames(pars0)), 1],
               ncase = ncase, ncontrol = ncontrol,
               alphaValue = pars0[grep("alpha0", rownames(pars0)), 1])

dataPP <- calculatePP(pars = parsPP,
                      dnData = dnData, mutData = mutData,
                      caseData = caseData, controlData = controlData,
                      geneName = as.character(geneName),
                      geneSet = geneSet)

dim(dataPP[dataPP$PP > 0.95, ])
dim(dataPP[dataPP$PP > 0.9, ])
dim(dataPP[dataPP$PP > 0.7, ])

dataPP[dataPP$pi0 > 0.5, ]
dataPP[(dataPP$BF > 50) & (dataPP$pi0 < 0.09), ]
##Bad data at the 11807

bayes.factor.CC(x = list(ca = 16, cn = 103), N = list(cn = ncontrol[1], ca = ncase[1]),
                gamma.mean = parsPP$gammaMeanCC[2], beta = parsPP$betaCC[2],
                rho1 = e.rho, nu1 = e.nu, rho0 = e.rho, nu0 = e.nu)

###Re-calculate Bayes Factors
BayesFactorCC3 <- function(x.case, x.control, Nsample,
                           gamma.meanCC, betaCC, rhoCC, nuCC){
    gAll <- range(rgamma(10000, gamma.meanCC*betaCC, rate = betaCC))
    gLower = gAll[1]; gUpper = gAll[2]

                                        #    print(range(gAll))
    altCC <- apply(cbind(x.case, x.control), 1, function(y){
        x2 <- list(ca = y[1], cn = y[2])
        evidence.alt.cc3 <- function(x = x2, N = Nsample, gamma.mean = gamma.meanCC, beta = betaCC,
                                     rho1 = rhoCC, nu1 = nuCC) {
            bControl <- log(dnbinom(x$cn, size = rho1, prob = nu1/(N$cn + nu1)))
                                        #                print(bControl)
            fCase <- function(gGamma) {
                dnbinom(x$ca, size = rho1 + x$cn, prob = (N$cn + nu1)/(N$cn + nu1 + N$ca*gGamma))*dgamma(gGamma, gamma.mean*betaCC, rate = betaCC)
            }
            bCase <- log(integrate(fCase, lower = gLower, upper = gUpper, stop.on.error = FALSE)$value)
                                        #print(bCase)
            return(exp(bCase + bControl))
        }
        t1 <- evidence.alt.cc3()

        return(t1)
    })
                                        #    print(altCC)
    nullCC <- apply(cbind(x.case, x.control), 1, function(y, rho1 = rhoCC, nu1 = nuCC, N = Nsample){
        x <- list(ca = y[1], cn = y[2])
        bControl <- log(dnbinom(x$cn, size = rho1, prob = nu1/(N$cn + nu1)))
        bCase <- log(dnbinom(x$ca, rho1 + x$cn, (N$cn + nu1)/(N$cn + nu1 + N$ca)))


        t1 <- exp(bCase + bControl)

        return(t1)
    })
                                        #    print(nullCC)

    tempBF <- altCC/nullCC #ifelse((x.case == 0) & (x.control == 0), 1, altCC/nullCC)
    return(tempBF)
}


#####################
#####################

# pdf(paste0("~/www/Test/xemD", nameAn, ".", gsub("/", "", dirFile), ".pdf"))
#par(mfrow = c(3, 3))
NK <- length(fileN)

for (ja in 1:NK){

#    par(mfrow = c(3, 2))
   load(paste0(dirFile, fileN[ja]))
    pdf(paste0("~/www/Test/AAplot_", annotationType, ".pdf"))
b1 <- as.data.frame(testIntegratedModel)
t1 <- grep('lambda0', colnames(b1))
#    if (length(t1) == 1)
 #       oTheta0[ja] <- median(b1[, 'lambda0[1]'])
  #  else
   #     oTheta0[ja] <- median(b1[, 'lambda0[1]'])
    #print(oTheta0[ja])
    #plot(density(b1[, 'lambda0[1]']))
#14892
    if (dim(b1)[1] > 0){
        apply(b1[1:3, 1:8], 2, quantile)

        tempLog <- c(annotationType, annotationType2,
                     sigmaPrior, round(median(b1[, 'lp__']), 3))

         nameAlpha <- colnames(b1)
    nameAlpha <- nameAlpha[grep("alpha0", nameAlpha)]
#tail(b1[, 'pi0[14892]'])
#tail(b1[, 'pi0[1]'])
#tail(b1[, 'pi0[2]'])

        geneSetNameN <- geneSetName
    balpha0 <- t(apply(b1[, nameAlpha], 2, function(x) quantile(x, c(0.025, 0.5, 0.975))))
    rownames(balpha0) <- c( "alpha0", geneSetNameN)

        if (is.element("alpha0Intercept", nameAlpha))
                rownames(balpha0) <- c( geneSetNameN, "alpha0")

        ###########Calculate logLK
        exp0 <- balpha0[1, 2]
        for (iG1 in 1:dim(geneSetM)[2]){
            exp0 <- exp0 + balpha0[iG1 + 1, 2]*geneSetM[, iG1]
                    }
        pE0 <- exp(exp0)/(1 + exp(exp0))
        b1E <- apply(b1, 2, median)
        gM <- b1E[grep("hyperGammaMeanDN", names(b1E))]
        bM <- b1E[grep("hyperBetaDN", names(b1E))]
        nFamily = mixDataKclasses$Ndn
        muRate <- mixDataKclasses$mutRate
        geneSet = geneSetM
        nCol = length(gM)
        nGeneSet = dim(geneSetM)[2]
        dD <- mixDataKclasses$dataDN

        llkE <- llk1(pars = c(gM, bM, balpha0[, 2][-1], balpha0[1,2]))

        outLog[[ja]] <- c(tempLog, llkE)

        ############################
    balpha0 <- balpha0[order(balpha0[, 2]),]
    head(balpha0)
    tail(balpha0)

    if (dim(balpha0)[1] > 55)
        balpha0 <- balpha0[((balpha0[, 1] > 0) & (balpha0[, 3] > 0)) |
                           ((balpha0[, 1] < 0) & (balpha0[, 3] < 0)),, drop = FALSE]


#plot(density(b1[, "fmrp.txt"]))
plot(c(-6, 6), c(0, dim(balpha0)[1] + 1), xlim = range(as.numeric(balpha0)),
     main = paste0(annotationType, ' & ', annotationType2, '\nsigmaPrior: ', sigmaPrior),
     xlab = '', ylab = '', col = 'white')

for (ik in 1:dim(balpha0)[1]){
    xTemp <- as.numeric(balpha0[ik, ]) #b1[, paste0('alpha0[', ik, ']')]
    lines(xTemp[c(1,3)], rep(ik, 2), lwd = 2, col = 'blue')
    text(xTemp[2], ik+0.2, paste0(rownames(balpha0)[ik]), cex = 0.4)
    text(xTemp[2], ik, 'x', col = ifelse((xTemp[1]>0) | (xTemp[3]<0), 'red', 'blue'))
    abline(v = 0, col = 'blue', cex = 2)

}
    } else {
        file.remove(paste0(dirFile, fileN[ja]))
    }
    dev.off()
}
#dev.off()

#################################
 dirFile <- "W_SCZgeneSetEachGeneSet/"
fileN <- dir(dirFile, ".RData")

source("TADA/TADA.R")
fileR <- dir("../extTADA/script", ".R")
for (ii in fileR)
    source(paste0("../extTADA/script/", ii))

fileR <- dir("../gTADA/script", ".R")
for (ii in fileR)
    source(paste0("../gTADA/script/", ii))
source("../gTADA/script/estimatePars.R")

dirFile <- "W_EACHannotation/"
fileN <- dir(dirFile, ".RData")

NN <- length(fileN)
#NN <- 2
for (i in 1:NN){
    load(paste0(dirFile, fileN[i]))
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
    print(geneSetName0)
}

####For de novo mutations
source("TADA/TADA.R")
fileR <- dir("../extTADA/script", ".R")
for (ii in fileR)
    source(paste0("../extTADA/script/", ii))

fileR <- dir("../gTADA/script", ".R")
for (ii in fileR)
    source(paste0("../gTADA/script/", ii))
source("../gTADA/script/estimatePars.R")

dirFile <- "W_EACHannotation/"
fileN <- dir(dirFile, ".RData")

NN <- length(fileN)
#NN <- 2
for (i in 1:NN){
    load(paste0(dirFile, fileN[i]))


parTemp <- c('alpha0[1]', 'alpha0[2]', 'hyperGammaMeanDN[1]', 'hyperGammaMeanDN[2]',
#                               'hyperGammaMeanCC[1]', 'hyperGammaMeanCC[2]',
                               'hyperBetaDN[1]', 'hyperBetaDN[2]')
            #,                               'hyperBetaCC[1]', 'hyperBetaCC[2]')

pars0 <- estimatePars(pars = parTemp,
                      mcmcResult = testIntegratedModel)

pars0 <- apply(pars0, 2, function(x) round(x, 3))
    pars0 <- cbind(rownames(pars0), pars0)
pars1 <- cbind(rep(geneSetName0, dim(pars0)[1]),
               rep(annotationType, dim(pars0)[1]),
               rep(annotationType2, dim(pars0)[1]),
               pars0)

pars1 <- cbind(pars1, rep(table(geneSetM[, 1])[2], dim(pars0)[1]))
#load("W_EACHannotation/DN.adjustBeta.1.3.nI2000.nThin.1.index.15_42_Apr_17_2017.lof_SCZ.damaging_SCZ.nGroupDN.2.nGroupCC.2.adjustRatio.1.0.05.sigmaPrior.2kF.1.RData")
write.table(pars1, paste0(resultDir, "/EstimatedPars.adjustBeta.", adjustHyperBeta, ".",
                         nChain, ".nI", nIteration, ".nThin.",
        nThin, ".index.", inDex, ".", annotationType, ".", annotationType2, ".nGroupDN.", nGroupDN,
        ".nGroupCC.", nGroupCC, ".adjustRatio.", round(adjustRatioMut, 2), ".", beta22,
        ".sigmaPrior.", sigmaPrior, ".", geneSetName0, 
        ".txt"), row.names = FALSE, quote = FALSE, col.names= FALSE)
}
################
load("W_AUTgeneSetEachGeneSet/AUTgeneSet.1..adjustRatioMut.1.3.nIteration2000.nThin.2.index.23_59_20_April.nGroupDN.3.nGroupCC.3.damaging.SameBeta.RData")

