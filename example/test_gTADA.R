## Load codes
dirFile <- "../script/"
fileR <- dir(dirFile, ".R$")

for (ii in fileR){
    print(ii)
    source(paste0(dirFile, ii))
}


## Get data
data <- read.table("../data/data_DD.txt", header = TRUE, as.is = TRUE)

ntrio = rep(4293, 2)
allDNData <- data[, paste0("dn_", c("damaging_", "lof_"), "DD")]
allMutData <- data[, paste0("mut_", c("damaging", "lof"))]

head(data.frame(allMutData, allDNData))

## Format data
inputData <- data.frame(Gene = data[, 1], allMutData, allDNData)

## Get gene sets
fmrp <- read.table("../data/FMRP_targets.txt")
fmrp <- data.frame(V1 = fmrp[, 1], rep(1, dim(fmrp)[1]))
fmrp <- merge(inputData, fmrp, by.x = 'Gene', by.y = 'V1', all.x = TRUE)
colnames(fmrp) <- paste0("V", 1:dim(fmrp)[2])
fmrp <- fmrp[, c('V1', 'V6')]
fmrp[, 2] <- ifelse(is.na(fmrp[, 2]), 0, 1)
write.table(fmrp, paste0("../data/fmrp_formatted.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)

fmrp <- read.table("../data/fmrp_formatted.txt", header = FALSE)

## Set parameters
nIteration = 1000
nThin = 1
nCore = 1

mcmcDD <- gTADA(modelName = DNextTADA, geneSet = data.frame(fmrp[, 2]),
                inputData = inputData, ## Input data should be formatted as above
                Ndn = array(c(ntrio)), #rep(ntrio, 1), ##Two de novo categories
                                        #                                    Ncase = array(ncase), #rep(N$ca, 1), ##Two case categories
                                        #                                   Ncontrol = array(ncontrol), #rep(N$cn, 1), ##Two control categories
                nIteration = nIteration ## Number of iterations: should be upto higher for real data
                )
