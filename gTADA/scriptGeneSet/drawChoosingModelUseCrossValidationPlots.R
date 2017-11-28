###Test CR
xC <- read.table("result.CR10.RealData.txt", as.is = TRUE)
xC <- read.table("result.CR10.RealData.LaplacePrior.txt", as.is = TRUE)
xC <- read.table("result.CR10.RealData.normal.txt", as.is = TRUE)
#xC <- read.table("result.CR10.RealData.txt2", as.is = TRUE)
library("dplyr")
xC1 <- xC %>% group_by(V1, V2, V3) %>% summarise(mK = median(V5), mL = median(V4))
xC2 <- as.data.frame(xC1)

xC2 <- data.frame(xC1)
xC2[, c(4, 5)] <- apply(xC2[, c(4, 5)], 2, function(x) round(x, 0))

pdf("~/www/Test/CVplot.pdf")
par(mfrow = c(2, 3))
for (ii in names(table(xC2[, 2]))){
    temp1 <- xC2[xC2[, 2] == ii,]
    plot(temp1[, 3], temp1[, 5], main = ii, xlab = 'sigmaPrior', ylab = 'logLLK', type = 'l',
         xlim = c(0, 15))
    t10 <-  temp1[which(temp1[, 5] == max(temp1[, 5])), 3][1]
    t11 <-  temp1[is.element(temp1[, 3], t10),][, 5][1]
    abline(v = t10, col = 'red')
    text(t10, t11 -5,
         paste0(t10, "-", t11), srt = 45, lwd = 0.7, cex = 0.7)

    if (ii == "damaging_DD")
        text(3.2, -4756, "*", col = 'red')
}
dev.off()

#####################
###
xC <- read.table("result.CR10.RealData.txt", as.is = TRUE)
xC <- read.table("result.CR10.RealData.LaplacePrior.txt", as.is = TRUE)

pdf("~/www/Test/CVplotAllSection.pdf")

for (i in 1:10) {
    
    xC <- read.table(paste0("result.CR10.RealData.normal.", i, ".txt"), as.is = TRUE)
#xC <- read.table("result.CR10.RealData.txt2", as.is = TRUE)
    xC <- xC[grep("damaging", xC[, 2]), ]
library("dplyr")
xC1 <- xC %>% group_by(V1, V2, V3) %>% summarise(mK = median(V5), mL = median(V4))
xC2 <- as.data.frame(xC1)

xC2 <- data.frame(xC1)
xC2[, c(4, 5)] <- apply(xC2[, c(4, 5)], 2, function(x) round(x, 0))


par(mfrow = c(3, 3))
for (ii in names(table(xC2[, 2]))){
    temp1 <- xC2[xC2[, 2] == ii,]
    plot(temp1[, 3], temp1[, 5], main = ii, xlab = 'sigmaPrior', ylab = 'logLLK', type = 'l',
         xlim = c(0, 15))
    t10 <-  temp1[which(temp1[, 5] == max(temp1[, 5])), 3][1]
    t11 <-  temp1[is.element(temp1[, 3], t10),][, 5][1]
    abline(v = t10, col = 'red')
    text(t10, t11 -5,
         paste0(t10, "-", t11), srt = 45, lwd = 0.7, cex = 0.7)

    if (ii == "damaging_DD")
        text(3.2, -4756, "*", col = 'red')
}

    }
    
dev.off()

#####################

##############
xCn <- xC[grep("DD", xC[, 1]),]
xCn <- xCn[grep("damaging", xCn[, 2]),]
xCn

x1 <- split(xCn[, 5], xCn[, 3])
