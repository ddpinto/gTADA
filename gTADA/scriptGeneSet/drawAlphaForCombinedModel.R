 .libPaths("~/InstallSoftware/Rlibs/")
library("rstan")
dirFile <- "W_AddGeneSetKF6_NormalPrior/"
 fileN <- dir(dirFile, "RData")
fileN
nameAn <- "damaging" #"missense" #
fileN <- fileN[grep(nameAn, fileN)]
fileN <- fileN[grep("sigmaPrior.3.2kF.5", fileN)]
fileN <- fileN[grep("10000", fileN)]


NK <- length(fileN)
##
for (ja in 1:NK){

#    par(mfrow = c(3, 2))
   load(paste0(dirFile, fileN[ja]))
    
    diseaseName <- substr(annotationType, 5, 7)
    if (diseaseName == "AST")
        diseaseName = "AUT"
    singleGene <- readLines(paste0("GeneListNameFromSingleModel.", diseaseName, ".txt"))
#is.element(rownames(balpha0), singleGene)    
    pdf(paste0("~/www/Test/AAplot_", annotationType, ".pdf"))
b1 <- as.data.frame(testIntegratedModel)

    t1 <- grep('lambda0', colnames(b1))
    if (dim(b1)[1] > 0){
        apply(b1[1:3, 1:8], 2, quantile)

        tempLog <- c(annotationType, annotationType2,
                     sigmaPrior, round(median(b1[, 'lp__']), 3))

         nameAlpha <- colnames(b1)
    nameAlpha <- nameAlpha[grep("alpha0", nameAlpha)]
        geneSetNameN <- geneSetName
    balpha0 <- t(apply(b1[, nameAlpha], 2, function(x) quantile(x, c(0.025, 0.5, 0.975))))
    rownames(balpha0) <- gsub(".txt", "", c( "alpha0", geneSetNameN))

        if (is.element("alpha0Intercept", nameAlpha))
                rownames(balpha0) <- c( geneSetNameN, "alpha0")

        ###########Calculate logLK
        bT2 <- balpha0[is.element(rownames(balpha0), singleGene), ]
        balpha0 <- rbind(balpha0[1, ], bT2)
        rownames(balpha0)[1] <- "alpha0"
        
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
     xlab = '',
      ylab = 'Gene set name',
      col = 'white', axes = FALSE)

        Axis(side=1, labels=TRUE)


for (ik in 1:dim(balpha0)[1]){
    xTemp <- as.numeric(balpha0[ik, ]) #b1[, paste0('alpha0[', ik, ']')]
    lines(xTemp[c(1,3)], rep(ik, 2), lwd = 2, col = 'blue')
    text(xTemp[2], ik+0.2, gsub(".txt", "", paste0(rownames(balpha0)[ik])), cex = 0.4)
    text(xTemp[2], ik, 'x', col = ifelse((xTemp[1]>0) | (xTemp[3]<0), 'red', 'blue'))
    abline(v = 0, col = 'blue', cex = 2)

}
    } else {
        file.remove(paste0(dirFile, fileN[ja]))
    }
    dev.off()
}
