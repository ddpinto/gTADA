dirFile <- "W_EPIgeneSet3/"
fileN <- dir(dirFile, "RData")
fileN
#fileN <- fileN[grep("5000", fileN)]
nameAn <- "damaging" #"missense" #"damaging" #
fileN <- fileN[grep(nameAn, fileN)]
fileN

#pdf("~/www/Test/xemD.pdf")
pdf(paste0("xemD", nameAn, ".", gsub("/", "", dirFile), ".pdf"))
for (ja in 1:length(fileN)){

load(paste0(dirFile, fileN[ja]))
b1 <- as.data.frame(testIntegratedModel)

#14892
apply(b1[1:3, 1:8], 2, quantile)
    nameAlpha <- colnames(b1)
    nameAlpha <- nameAlpha[grep("alpha0", nameAlpha)]
tail(b1[, 'pi0[14892]'])
tail(b1[, 'pi0[1]'])
tail(b1[, 'pi0[2]'])

    balpha0 <- t(apply(b1[, nameAlpha], 2, function(x) quantile(x, c(0.025, 0.5, 0.975))))
    rownames(balpha0) <- c("alpha0", geneSetName)
    balpha0 <- balpha0[order(balpha0[, 2]),]
    head(balpha0)
    tail(balpha0)

    if (dim(balpha0)[1] > 15)
        balpha0 <- balpha0[(balpha0[, 1] > 0) | (balpha0[, 3] < 0),]

    
plot(c(-4, 4), c(0, dim(balpha0)[1] + 1), main = paste0(annotationType, ' & ', annotationType2),
     xlab = '', ylab = '', col = 'white')
for (ik in 1:dim(balpha0)[1]){
    xTemp <- as.numeric(balpha0[ik, ]) #b1[, paste0('alpha0[', ik, ']')]
    lines(xTemp[c(1,3)], rep(ik, 2))
    text(xTemp[2], ik+0.2, paste0(rownames(balpha0)[ik]))
    text(xTemp[2], ik, 'x', col = 'red')
    abline(v = 0, col = 'blue', cex = 2)

}
}

dev.off()

b22 <- apply(b1, 2, median)

tail(b22)

str(testIntegratedModel@sim$samples)
