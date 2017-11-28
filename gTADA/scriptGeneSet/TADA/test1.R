source("TADA.R")

d1 <- read.table("~/Documents/Github/TestR/TestDisorders/OCD/DNinputFromOCD.txt", as.is = TRUE, header = TRUE)

cList <- NULL
for (ii in 1:1000){
aD <- simulator.denovo(N = 184, mu = d1$mut.cls1, mu.frac = 1, pi = 0.0170,
                       gamma.mean = 31.734, #63.8,
                       beta = 1)

cList[ii] <- sum(aD[[1]])
}
hist(cList, 100)
abline(v = sum(d1$dn.cls1), col = 'red')
abline(v = quantile(cList, c(0.025, 0.975)))


cList1 <- NULL
for (ii in 1:1000){
aD <- simulator.denovo(N = 184, mu = d1$mut.cls2, mu.frac = 1, pi = 0.0170,
                       gamma.mean = 26.70157668, #30.29,
                       beta = 1)

cList1[ii] <- sum(aD[[1]])
}

hist(cList1, 100, xlim = c(50, 200))
abline(v = sum(d1$dn.cls2), col = 'red')
abline(v = quantile(cList1, c(0.025, 0.975)))

##TOURETTE DISORDER
d1 <- read.table("/Users/hoangnguyen/Documents/Github/TestR/TestDisorders/Tourette/Tourette.Willsey2017.txt",
                 header = TRUE, as.is = TRUE)

cList <- NULL
for (ii in 1:1000){
aD <- simulator.denovo(N = 484, mu = d1$mut.lof, mu.frac = 1, pi = 0.02369,
                 gamma.mean = 49.1155, beta = 1)

cList[ii] <- sum(aD[[1]])
}

hist(cList, 100)
abline(v = sum(d1$dn.lof), col = 'red')
abline(v = quantile(cList, c(0.025, 0.975)))

cList1 <- NULL
for (ii in 1:1000){
aD <- simulator.denovo(N = 484, mu = d1$mut.mis3, mu.frac = 1, pi = 0.02369,
                 gamma.mean = 10.39, beta = 1)

cList1[ii] <- sum(aD[[1]])
}

hist(cList1, 100)
abline(v = sum(d1$dn.mis3), col = 'red')
abline(v = quantile(cList1, c(0.025, 0.975)))

######################
#####################DATA from TADAv1.1
d1 <- read.table("~/Documents/Github/TestR/TADAv1.1/TADA_demo_counts_de-novo_only.txt", as.is = TRUE, header = TRUE)

cList <- NULL
for (ii in 1:1000){
aD <- simulator.denovo(N = 4500, mu = d1$mut.cls1, mu.frac = 1, pi = 0.05,
                       gamma.mean = 20, #31.734, #63.8,
                       beta = 1)

cList[ii] <- sum(aD[[1]])
}
hist(cList, 100)
abline(v = sum(d1$dn.cls1), col = 'red')
abline(v = quantile(cList, c(0.025, 0.975)))


cList1 <- NULL
for (ii in 1:1000){
aD <- simulator.denovo(N = 4500, mu = d1$mut.cls2, mu.frac = 1, pi = 0.05,
                       gamma.mean = 4.7, #26.70157668, #30.29,
                       beta = 1)

cList1[ii] <- sum(aD[[1]])
}

hist(cList1, 100) #, xlim = c(50, 200))
abline(v = sum(d1$dn.cls2), col = 'red')
abline(v = quantile(cList1, c(0.025, 0.975)))


