NN = 100

d1 <- NULL
for (i in 1:4){
    x1 <- sample(0:1, NN, replace = TRUE)
    d1 <- cbind(d1, x1)
}

d1pca <- prcomp(d1)

library(ggfortify)

d2 <- cbind(d1, d1[, c(1, 2)])
d2pca <- prcomp(d1)

autoplot(d2pca)
