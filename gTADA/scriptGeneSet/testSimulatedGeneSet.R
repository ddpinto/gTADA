x1 <- read.table("xemGeneSet.txt")

#tOut <- c(          dim(d11)[1],
#          dim(d11[d11[, 2] == 1,])[1],
 #   gammaMeanDenovo, gammaMeanDenovo2, betaDN,
  #        returnValue0$par)

#t11 <- c(dim(d11)[1],
#          dim(d11[d11[, 2] == 1,])[1],
 #         returnValue1$par)



#tOut1 <- c(tOut, t11, length(sGeneSet)/dim(data)[1], p11)
#tOut1 <- round(tOut1, 4)

x1 <- x1[x1[, 1] < 5000,]

library("ggplot2")

pG1 <- ggplot(x1, aes(V3, V7, col = as.factor(V21))) + geom_line()
pG2 <- ggplot(x1, aes(V4, V8, col = as.factor(V21))) + geom_line()

pG3 <- ggplot(x1, aes(V3, V14, col = as.factor(V21))) + geom_line()
pG4 <- ggplot(x1, aes(V4, V15, col = as.factor(V21))) + geom_line()


p1 <- ggplot(x1, aes(V1, V2/V1, col = as.factor(V21))) + geom_line()
p2 <- ggplot(x1, aes(V12, V13/V12, col = as.factor(V21))) + geom_line()


pdf("~/www/Test/TgeneSet.pdf")
print(pG1)
print(pG2)
print(pG3)
print(pG4)

print(p1)
print(p2)
dev.off()
