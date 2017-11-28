simulator2 <- function(N, muAll, mu.frac, pi,
                       alpha0, geneSet,
                       gamma.mean.dn, beta.dn, gamma.mean.CC, beta.CC, rho1, nu1, rho0, nu0, tradeoff=FALSE) {
      m <- dim(muAll)[1] # number of genes
        K <- length(mu.frac) # number of mutational categories
      exp0 <- alpha0[1]
      for (i in 1:dim(geneSet)[2])
          exp0 <- exp0 + alpha0[1+i]*geneSet[, i]
      pp0 <- exp(exp0)/(1 + exp(exp0))
      z0 <- sapply(pp0, function(x)           sample(0:1, 1, prob = c(1 - x, x)))
      

      z = z0
        gamma.dn <- array(1, dim=c(m,K))
        gamma.CC <- array(1, dim=c(m,K))
        q <- array(0, dim=c(m,K))
        x <- array(0, dim=c(m,3*K))
        k <- sum(z==1)
        for (j in 1:K) {
                mu <- muAll[, j]
                    # sample de novo
                    gamma.dn[z==1, j] <- rgamma(k, gamma.mean.dn[j]*beta.dn[j], beta.dn[j])
                    col <- 3*(j-1)+1
                    x[,col] <- rpois(m, 2 * mu * mu.frac[j] * gamma.dn[,j] * N$dn)

                    # sample case-control
                    gamma.CC[z==1, j] <- rgamma(k, gamma.mean.CC[j]*beta.CC[j], beta.CC[j])
                    q[z==0, j] <- rgamma(m-k, rho0[j], nu0[j])
                    if (tradeoff==FALSE) {
                              q[z==1, j] <- rgamma(k, rho1[j], nu1[j])
                    } else {
                              q[z==1, j] <- mu[z==1] * mu.frac[j] / (delta[j] * gamma.CC[z==1, j])
                    }
                    x[,col+1] <- rpois(m, q[,j] * gamma.CC[,j] *N$ca)
                    x[,col+2] <- rpois(m, q[,j] * N$cn)

        }

        sample.info <- cbind(mu, z, gamma.dn, gamma.CC, q, x)

        return (list(sample=x, sample.info=sample.info))
}
