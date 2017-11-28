 fileName <- "result.SingleGeneSet.DD.txt"
#fileName <- "result.SingleGeneSet.EPI.txt"
#fileName <- "result.SingleGeneSet.ID.txt"
#fileName <- "result.SingleGeneSet.CHD.txt"
d1 <- read.table(fileName, as.is = TRUE, header = FALSE)
t1 <- d1[, 1]
d1[, 1] <- d1[, 2]
d1[, 2] <- t1
d1[, 3:7] <- d1[, 4:8]
#d1[, 1]
#fileName <- "result.SingleGeneSet.SCZ.txt"; d1 <- read.table(fileName, as.is = TRUE, header = FALSE)
#fileName <- "result.SingleGeneSet.AUT.txt"; d1 <- read.table(fileName, as.is = TRUE, header = FALSE)
dName <- unlist(strsplit(fileName, ".", fixed = TRUE))[3]
d1 <- d1[order(d1[, 4]),]

aT <- range(unlist(as.list(d1[, c(4, 5, 6)])))

d1 <- d1[grep("denovo.epi", d1[, 2], invert = TRUE), ]
d1 <- d1[grep("SNPsINdel.denovo.epi", d1[, 2], invert = TRUE), ]
d1 <- d1[grep("1.essentila", d1[, 2], invert = TRUE), ]
d1 <- d1[grep("denovo.scz", d1[, 2], invert = TRUE),]
d1 <- d1[d1[, 2] != "dd.txt", ] ##The same list with McRae2016
d1 <- d1[d1[, 2] != "denovo.aut.txt", ] ##The same list with SNPsINdel.denovo.aut.txt
#d1 <- d1[d1[, 2] != "SNPsINdel.denovo.dd.txt", ]
#d1 <- d1[d1[, 2] != "SNPsINdel.denovo.aut.txt", ]
d1 <- d1[d1[, 2] != "newListFromSupTable2.txt", ]
d1 <- d1[d1[, 2] != "listMcRae2016.txt", ]
d1 <- d1[d1[, 2] != "listLoFtolerantEXAC.txt", ]
d1 <- d1[grep("FDR005", d1[, 2], invert = TRUE),]
d1 <- d1[grep("FDR01", d1[, 2], invert = TRUE),]
d1 <- d1[grep("extTADA", d1[, 2], invert = TRUE),]
d1[, 2] <- gsub(".txt", "", d1[, 2])
d1[, 2] <- gsub("listgenePanous2016", "hiPSCneurons", d1[, 2])

write.table(d1[, 2],
            paste0("GeneListNameFromSingleModel.", dName, ".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)


pdf(paste0("~/www/Test/", dName, "EnrichmentResultsForSingleGeneSet.pdf"))
plot(c(aT[1] - 0.5, aT[2] + 0.5), c(0, dim(d1)[1] + 1),
     col = 'white',
     xlab = expression(paste("Gene-set ", alpha)), ylab = 'Gene set name',
     yaxt='n', main = paste0(dName, "\n(", dim(d1)[1], " gene sets)"), axes = FALSE)
Axis(side = 1, labels = TRUE)

for (i in 1:dim(d1)[1]){
    lines(c(d1[i, 5], d1[i, 6]), c(i, i), col = 'blue')
    text(d1[i, 4], i, "*", col = 'red')
    text(d1[i, 6], i + 0.3, gsub("AST", "ASD", d1[i, 2]), cex = 0.25)
}
abline(v = 0)
dev.off() 

###########################
############################
##ASD
d1 <- read.table("result.SingleGeneSet.AUT.txt", as.is = TRUE, header = FALSE)
d1 <- d1[order(d1[, 4]),]

d1 <- d1[grep("denovo.epi", d1[, 2], invert = TRUE), ]
d1 <- d1[grep("SNPsINdel.denovo.epi", d1[, 2], invert = TRUE), ]
d1 <- d1[grep("1.essentila", d1[, 2], invert = TRUE), ]
d1 <- d1[grep("denovo.scz", d1[, 2], invert = TRUE),]
d1 <- d1[d1[, 2] != "dd.txt", ] ##The same list with McRae2016
d1 <- d1[d1[, 2] != "denovo.aut.txt", ] ##The same list with SNPsINdel.denovo.aut.txt
#d1 <- d1[d1[, 2] != "SNPsINdel.denovo.dd.txt", ]
#d1 <- d1[d1[, 2] != "SNPsINdel.denovo.aut.txt", ]
d1 <- d1[d1[, 2] != "newListFromSupTable2.txt", ]
d1 <- d1[d1[, 2] != "listMcRae2016.txt", ]
d1 <- d1[d1[, 2] != "listLoFtolerantEXAC.txt", ]
d1 <- d1[grep("FDR005", d1[, 2], invert = TRUE),]
d1 <- d1[grep("FDR01", d1[, 2], invert = TRUE),]
d1 <- d1[grep("extTADA", d1[, 2], invert = TRUE),]
d1[, 2] <- gsub(".txt", "", d1[, 2])
d1[, 2] <- gsub("listgenePanous2016", "hiPSCneurons", d1[, 2])


aT <- range(unlist(as.list(d1[, c(4, 5, 6)])))

pdf(paste0("~/www/Test/AUT.EnrichmentResultsForSingleGeneSet.pdf"))
plot(c(aT[1] - 0.5, aT[2] + 0.5), c(0, dim(d1)[1] + 1),
     col = 'white',
     xlab = '', ylab = 'Gene set name',
     yaxt='n', main = 'ASD', axes = FALSE)

Axis(side=1, labels=TRUE)

for (i in 1:dim(d1)[1]){
    lines(c(d1[i, 5], d1[i, 6]), c(i, i), col = 'blue')
    text(d1[i, 4], i, "*", col = 'red')
    text(d1[i, 6], i + 0.3, gsub("AST", "ASD", d1[i, 2]), cex = 0.35)
    abline(v = 0)
}

dev.off()

