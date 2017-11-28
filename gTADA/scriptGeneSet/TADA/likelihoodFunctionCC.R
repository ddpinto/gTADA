
calculateProb <- function(x, N, mu){

    gamma.mean.dn <- hyperpar[1] #gamma.dn|H1 ~ Gamma(gamma.mean.dn*beta.dn, beta.n)
    beta.dn <- hyperpar[2]
    gamma.mean.CC <- hyperpar[3] #gamma.mean.CC|H1 ~ Gamma(gamma.mean.CC*beta.CC, beta.CC)
    beta.CC <- hyperpar[4]
    rho1 <- hyperpar[5] #q|H1 ~ Gamma(rho1, nu1)
    nu1 <- hyperpar[6]
    rho0 <- hyperpar[7] #q|H0 ~ Gamma(rho0, nu0)
    nu0 <- hyperpar[8]


    x <- as.numeric(x)

    k <- length(mu)

    prob.gene.CC <- rep(1, k)
    prob.gene.dn <- rep(1, k)
    gamma.mean.dn <- hyperpar[3]

###Calculate for all categories


    evidence.alt.CC <- function(x, N, gamma.mean, beta, rho1, nu1, q.lower=1e-8, q.upper=0.1, debug=FALSE) {
  integrand <- function(u) {
    q <- exp(u)
    return (dnbinom(x$ca, gamma.mean*beta, beta/(beta+N$ca*q)) * dgamma(q, rho1+x$cn, nu1+N$cn) * exp(u))
  }

  marglik1.ctrl <- dnbinom(x$cn, rho1, nu1/(nu1+N$cn))
  marglik1.case <- integrate(integrand, lower=log(q.lower), upper=log(q.upper))$value

  marglik1 <- marglik1.ctrl * marglik1.case

  #   return (exp(marglik1.ctrl.log+marglik1.case.log))
  return (list(ctrl=marglik1.ctrl, case=marglik1.case, total=marglik1))
}



}

########################################
logLik <- function(hyperpar, pi, counts, N, mu, prior.weight=0, denovo.only=FALSE,
                   q.lower=1e-8, q.upper=0.1, debug=FALSE) {
    gamma.mean.dn <- hyperpar[1, ]
  beta.dn <- hyperpar[2, ]
  gamma.mean.CC <- hyperpar[3, ]
  beta.CC <- hyperpar[4, ]
  rho1 <- hyperpar[5, ]
  nu1 <- hyperpar[6, ]
  rho0 <- hyperpar[7, ]
  nu0 <- hyperpar[8, ]



  n <- nrow(counts)
  prob <- numeric(n)
  posterior <- numeric(n)
  bf <- numeric(n)

    ###Number of annotation types
    nT <- dim(mu)[2]

    prob.M1.list <- prob.M0.list <- NULL

  for (i in 1:n)  {

      ###Extract the genes
      xRow <- as.numeric(counts[i, ])
      prob.M1.Type <- NULL#rep(1, nT)
      prob.M0.Type <- NULL#rep(1, nT)


      if (denovo.only==TRUE) {
#          prob.M1 <- evidence.alt.dn(counts[i,"dn"], N$dn, mu[i, ], gamma.mean.dn, beta.dn)

 #     prob.M0 <- evidence.null.dn(counts[i,"dn"], N$dn, mu[i, ])

          for (jD in 1:nT){
              start <- 3*(jD-1) + 1


              prob.M1.Type[jD]  <- dnbinom(xRow[start], gamma.mean.dn[jD]*beta.dn[jD],
                                      beta.dn[jD]/(beta.dn[jD] + 2*N$dn*mu[i, jD]))

              prob.M0.Type[jD] <- dpois(xRow[start], 2*N$dn*mu[i, jD])
          }
      }
      #############Calculate for de novo + case-control
      else {

          for (jD in 1:nT){
              start <- 3*(jD-1) + 1

              ##For alternative##################
              message("xRow: ", xRow, " nT: ", nT)

              integrand <- function(u) {
                  q <- exp(u)

                  return (dnbinom(xRow[start + 1], gamma.mean.CC[jD]*beta.CC[jD], beta.CC[jD]/(beta.CC[jD] + N$ca*q)) * dgamma(q, rho1[jD] + xRow[start + 2], nu1[jD] + N$cn) * exp(u))
              }

              marglik1.ctrl <- dnbinom(xRow[start + 1], rho1[jD], nu1[jD]/(nu1[jD] + N$cn))

              #marglik1.case <- integrate(integrand, lower=log(q.lower), upper=log(q.upper))$value

              marglik1.case <- quad(integrand, xa = log(q.lower), xb = log(q.upper))

              marglik1 <- marglik1.ctrl * marglik1.case
              prob.CC.alt <- marglik1

              prob.DN  <- dnbinom(xRow[start], gamma.mean.dn[jD]*beta.dn[jD],
                                      beta.dn[jD]/(beta.dn[jD] + 2*N$dn*mu[i, jD]))

              prob.M1.Type[jD] <- exp(log(sum(prob.CC.alt, prob.DN)))

              ##For null############################
                      #  marglik0.ctrl.log <- log(dnbinom(x$cn, rho0, nu0/(nu0+N$cn)))
          marglik0.ctrl.log <- log(dnbinom(xRow[start + 2], rho0[jD], nu0[jD]/(nu0[jD]  + N$cn)))
          # marglik0.case.log <- log(dnbinom(x$ca, rho0+x$cn, (nu0+N$cn)/(nu0+N$cn+N$ca)))
          marglik0.case.log <- log(dnbinom(xRow[start + 1], rho0[jD] + xRow[start + 2], (nu0[jD] + N$cn)/(nu0[jD] + N$cn + N$ca)))
          marglik0.log <- marglik0.ctrl.log + marglik0.case.log

              prob.CC.null <- exp(marglik0.log)
              prob.DN.null <- dpois(xRow[start], 2*N$dn*mu[i, jD])

              prob.M0.Type[jD] <- exp(log(sum(prob.CC.null, prob.DN.null)))





          }
      } ###End else

      if (debug) cat("prob.M1 = ", prob.M1,"\n")
      if (debug) cat("Prob.M0 = ", prob.M0, "\n")

      prob.M1 <- exp(sum(log(prob.M1.Type)))
      prob.M0 <- exp(sum(log(prob.M0.Type)))

      prob.M0.list[i] <- prob.M0
      prob.M1.list[i] <- prob.M1

      prob[i] <- pi * prob.M1 + (1 - pi) * prob.M0
    posterior[i] <- pi*prob.M1 / (prob[i])
    bf[i] <- prob.M1 / prob.M0
  }
  marg.nll <- -sum(log(prob[!is.na(prob)]))

  # add the prior
  if (prior.weight > 0) {
    u <- prior.weight
    q1.mean <- rho1/nu1
    q0.mean <- rho0/nu0
    marg.nll <- marg.nll - log(dgamma(q1.mean, q0.mean*u, u))
  }

  result <- list(prob=cbind(prob, prob.M0.list, prob.M1.list), marginal=marg.nll, posterior=posterior, bf=bf)
  return (result)
}
