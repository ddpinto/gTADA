dirFile <- "TestGeneSetDNmcmc4_April/"
dirFile <- "TestGeneSetDNmcmc4_April2/"
dirFile <- "TestGeneSetDNmcmc4_April3/"
fileN <- dir(dirFile, "RData")


ja = 1
NK <- length(fileN) #10
dOut <- NULL
for (ja in 1:NK){
    
load(paste0(dirFile, fileN[ja]))

sGeneSet0 <- sDenovoCC2[[2]][, 2]
    names(sGeneSet0) <- 1:length(sGeneSet0)

    s1 <- as.numeric(names(sGeneSet0[sGeneSet0 > 0]))

    names(geneSetM) <- 1:length(geneSetM)
    s2 <- as.numeric(names(geneSetM[geneSetM > 0]))
    
table(geneSet)

bInter <- data.frame(sGeneSet0, geneSet)

tt <- c(length(s1), length(s2), length(intersect(s1, s2))) ##c(length(sGeneSet0[sGeneSet0 == 1]), length(geneSet[geneSet == 1]), length(sGeneSet0[geneSet == 1]))
dOut <- rbind(dOut, tt)

}

dOut1 <- cbind(dOut, dOut[, 3]/dOut[, 1])
dOut1[dOut1[, 4] > 1, ]
