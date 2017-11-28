gTADA <- function(modelName  ,
                         dataDN = NULL, mutRate = NULL, Ndn = NULL,
                  dataCCcase = NULL, dataCCcontrol = NULL, Ncase = NULL, Ncontrol = NULL,
                  geneSet = NULL,
                    nIteration = NULL, nIteration2 = NULL,
                    nThin = NULL, nThin2 = NULL, nCore = 1, nChain = 1,
                         hyperBetaDN0 = NULL,
                         hyperBetaCC0 = NULL,
                         hyper2GammaMeanDN = NULL, hyper2BetaDN = NULL, ##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                         hyper2GammaMeanCC = NULL, hyper2BetaCC = NULL,

                         upperPi0 = 0.5, lowerPi0 = 0, lowerBeta = 0, ##Lower and upper limits of pi: should be default
                         lowerHyperGamma = 1, lowerGamma = 1, #Should be default
                         betaPars = c(6.7771073, -1.7950864, -0.2168248), #Adjust beta's values: should be default
                    adjustHyperBeta = as.integer(1), ##1 if want to adjust beta, 0 if don't want to adjust beta
                  autoAdjustHyperBeta = FALSE,
                  sigmaPrior = 2)
     {

         if (is.null(nIteration))
             stop("======\nNo input for the nIteration parameter\n=====")
         if ((is.null(dataDN) | is.null(mutRate)) & (is.null(dataCCcase) | is.null(dataCCcontrol)))
             stop("Need to have input data: only DN, only CC or DN + CC")
         if (is.null(nThin))
             nThin = ifelse(nIteration > 1000, floor(nIteration/1000), 1)
         if (nCore != nChain)
             warning("nCore is different from nChain")
         if (!is.null(dataDN))
             Ngene <- dim(dataDN)[1]
         if (!is.null(dataCCcase))
             Ngene <- dim(dataCCcase)[1]
         message("\nThere are ", Ngene, " genes in this analysis")

         if (is.null(hyper2GammaMeanDN) & !is.null(dataDN))
             hyper2GammaMeanDN <- rep(1, dim(dataDN)[2])
         if (is.null(hyper2BetaDN) & !is.null(dataDN))
             hyper2BetaDN <- rep(0.025, dim(dataDN)[2])
        if (is.null(hyper2GammaMeanCC) & !is.null(dataCCcase))
             hyper2GammaMeanCC <- rep(1, dim(dataCCcase)[2])
         if (is.null(hyper2BetaCC) & !is.null(dataCCcase))
             hyper2BetaCC <- rep(0.2, dim(dataCCcase)[2])
         if ((adjustHyperBeta == 0) & is.null( hyperBetaDN0))
             hyperBetaDN0 <- rep(1, dim(dataDN)[2])
         if ((adjustHyperBeta == 0) & is.null( hyperBetaCC0) & !is.null(dataCCcase))
             hyperBetaCC0 <- rep(1, dim(dataCCcase)[2])

         NCdn <- dim(dataDN)[2]
         NCcc <- dim(dataCCcase)[2]
         if (is.null(Ncase))              Ncase = 0
         if (is.null(Ncontrol)) Ncontrol = 1
         if (is.null(Ndn)) Ndn = 0
         if (is.null(hyperBetaDN0)) hyperBetaDN0 <- rep(1, length(dataDN[1, ]))
         if (is.null(hyperBetaCC0)) hyperBetaCC0 <- rep(4, length(dataCCcase[1, ]))
         if (is.null(hyper2GammaMeanDN)) hyper2GammaMeanDN = rep(1, length(dataDN[1, ]))
         if (is.null(hyper2GammaMeanCC)) hyper2GammaMeanCC = rep(4, length(dataCCcase[1, ]))
         if (is.null(hyper2BetaDN)) hyper2BetaDN = rep(0.05, length(dataDN[1, ]))
         if (is.null(hyper2BetaCC)) hyper2BetaCC = rep(0.2, length(dataCCcase[1, ]))


         Tgs <- 0
         if (!is.null(geneSet)) {
             Tgs = dim(geneSet)[2]
             Ngs = Tgs
         } else {
             geneSet <- data.frame(sample(0:1, Ngene, replace = TRUE))
             Ngs = dim(geneSet)[2]
             message("No gene set input\n")

         }


         modelData <- list(NN = Ngene, #Gene numbers
                           K = 2, #Hypothesis numbers: should be default
                           NCdn = NCdn, #Number of de novo classes
                           NCcc = NCcc,
                           Ndn = array(Ndn), # Family numbers
                           Ncase = array(Ncase), Ncontrol = array(Ncontrol), Ntotal = array(Ncase + Ncontrol),
                           dataDN = dataDN, # Denovo data
                           mutRate = mutRate, # Mutation rates
                           dataCCcase = data.frame(dataCCcase), dataCCtotal = data.frame(dataCCcase + dataCCcontrol),
                           thetaH0 = array(Ncase/(Ncase + Ncontrol)),
                           betaPars = array(betaPars), #Adjust beta's values: should be default
                           adjustHyperBeta = adjustHyperBeta, ##1 if want to adjust beta, 0 if don't want to adjust beta
                           upperPi0 = upperPi0, lowerPi0 = lowerPi0, lowerBeta = lowerBeta, ##Lower and upper limits of pi: should be default
                          lowerHyperGamma = lowerHyperGamma, lowerGamma = lowerGamma, #Should be default
                          hyperBetaDN0 = array(hyperBetaDN0),
                          hyperBetaCC0 = array(hyperBetaCC0),
                            hyper2GammaMeanDN = array(hyper2GammaMeanDN),##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                          hyper2BetaDN = array(hyper2BetaDN),
                          hyper2GammaMeanCC = array(hyper2GammaMeanCC),
                          hyper2BetaCC = array(hyper2BetaCC),
                          Tgs = Tgs,
                          Ngs = Ngs,
                          geneSet = geneSet,
                          sigmaPrior = sigmaPrior)

         message("\n=============FIRST TIME ===============\n")
         message("\nSampling with nter = ", nIteration, " and nThin = ", nThin, "\n")
         message("\nThe model ", deparse(substitute(modelName)), " is used\n")
         mcmcModel <- stan(model_code = modelName,
                        data = modelData, ##Model data as described above
                        iter = nIteration, chains = nChain, cores = nCore, thin = nThin)

         ####################Finish the FIRST TIME
####Re-sample using new hyper betas
         if ( autoAdjustHyperBeta){
             if (!is.null(nIteration2))
                 nIteration = nIteration2

             if (is.null(nThin2))
                 nThin2 = ifelse(nIteration > 1000, floor(nIteration/1000), 1)
             nThin = nThin2

         mcmcData <- as.data.frame(mcmcModel)
         cName <- colnames(mcmcData)
         hyperGammaMeanNameCC <- cName[grep("hyperGammaMeanCC", cName)]
         hyperGammaMeanNameDN <- cName[grep("hyperGammaMeanDN", cName)]

             hyperGammaMeanDN0 <- as.numeric(apply(data.frame(mcmcData[, hyperGammaMeanNameDN]), 2, median))
             hyperGammaMeanCC0 <- as.numeric(apply(data.frame(mcmcData[, hyperGammaMeanNameCC]), 2, median))

             hyperBetaDN0 <- exp(betaPars[1]*hyperGammaMeanDN0^(betaPars[2]) + betaPars[3])
             hyperBetaCC0 <- exp(betaPars[1]*hyperGammaMeanCC0^(betaPars[2]) + betaPars[3])

             adjustHyperBeta <- 0

         modelData <- list(NN = Ngene, #Gene numbers
                           K = 2, #Hypothesis numbers: should be default
                           NCdn = NCdn, #Number of de novo classes
                           NCcc = NCcc,
                           Ndn = array(Ndn), # Family numbers
                           Ncase = array(Ncase), Ncontrol = array(Ncontrol), Ntotal = array(Ncase + Ncontrol),
                           dataDN = dataDN, # Denovo data
                           mutRate = mutRate, # Mutation rates
                           dataCCcase = data.frame(dataCCcase), dataCCtotal = data.frame(dataCCcase + dataCCcontrol),
                           thetaH0 = array(Ncase/(Ncase + Ncontrol)),
                           betaPars = array(betaPars), #Adjust beta's values: should be default
                           adjustHyperBeta = adjustHyperBeta, ##1 if want to adjust beta, 0 if don't want to adjust beta
                           upperPi0 = upperPi0, lowerPi0 = lowerPi0, lowerBeta = lowerBeta, ##Lower and upper limits of pi: should be default
                          lowerHyperGamma = lowerHyperGamma, lowerGamma = lowerGamma, #Should be default
                          hyperBetaDN0 = array(hyperBetaDN0),
                           hyperBetaCC0 = array(hyperBetaCC0),
                            hyper2GammaMeanDN = array(hyper2GammaMeanDN),##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                          hyper2BetaDN = array(hyper2BetaDN),
                          hyper2GammaMeanCC = array(hyper2GammaMeanCC),
                          hyper2BetaCC = array(hyper2BetaCC)  )

             message("\n=============SECOND TIME ===============\n")
             message("\nThis process is using beta values estimated from the previous process\n")
         message("\nSampling with nter = ", nIteration, " and nThin = ", nThin, "\n")
         message("\nThe model ", deparse(substitute(modelName)), " is used\n")
         mcmcModel <- stan(model_code = modelName,
                        data = modelData, ##Model data as described above
                        iter = nIteration, chains = nChain, cores = nCore, thin = nThin)
         }



         return(mcmcModel)
          }



###Bayes Factor For Case-control
BayesFactorCC3 <- function(x.case, x.control, Nsample,
                           gamma.meanCC, betaCC, rhoCC, nuCC){
    gAll <- range(rgamma(10000, gamma.meanCC*betaCC, rate = betaCC))
    gLower = gAll[1]; gUpper = gAll[2]

                                        #    print(range(gAll))
    altCC <- apply(cbind(x.case, x.control), 1, function(y){
        x2 <- list(ca = y[1], cn = y[2])
        evidence.alt.cc3 <- function(x = x2, N = Nsample, gamma.mean = gamma.meanCC, beta = betaCC,
                                     rho1 = rhoCC, nu1 = nuCC) {
            bControl <- log(dnbinom(x$cn, size = rho1, prob = nu1/(N$cn + nu1)))
                                        #                print(bControl)
            fCase <- function(gGamma) {
                dnbinom(x$ca, size = rho1 + x$cn, prob = (N$cn + nu1)/(N$cn + nu1 + N$ca*gGamma))*dgamma(gGamma, gamma.mean*betaCC, rate = betaCC)
            }
            bCase <- log(integrate(fCase, lower = gLower, upper = gUpper, stop.on.error = FALSE)$value)
                                        #print(bCase)
            return(exp(bCase + bControl))
        }
        t1 <- evidence.alt.cc3()

        return(t1)
    })
                                        #    print(altCC)
    nullCC <- apply(cbind(x.case, x.control), 1, function(y, rho1 = rhoCC, nu1 = nuCC, N = Nsample){
        x <- list(ca = y[1], cn = y[2])
        bControl <- log(dnbinom(x$cn, size = rho1, prob = nu1/(N$cn + nu1)))
        bCase <- log(dnbinom(x$ca, rho1 + x$cn, (N$cn + nu1)/(N$cn + nu1 + N$ca)))


        t1 <- exp(bCase + bControl)

        return(t1)
    })
                                        #    print(nullCC)

    tempBF <- altCC/nullCC #ifelse((x.case == 0) & (x.control == 0), 1, altCC/nullCC)
    return(tempBF)
}

################


library(locfit)
loc2plot <- function(x,y,cprob=0.5, xlim,gxlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,scale=c(sc1,sc2), xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	lev <- sort(fitted(fit))[floor(cprob*length(x))]
	plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
	xlim=gxlim,...)
}

# finds the mode for a bivariate density
loc2mode <- function(x,y,alpha=0.5,xlim,maxk,...)
{
	sc1 <- sqrt(var(x))
	sc2 <- sqrt(var(y))
	if(missing(maxk)) maxk <- 100
	if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2), maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),xlim=xlim,maxk=maxk,mint=100,maxit=100)
#,cut=0.8
	tt <- max(fitted(fit))
	wt <- fitted(fit) == tt
	c(x[wt][1],y[wt][1])
}

# this works for univariate data; gives a vector with
# mode (global), hpd1_low, hpd1_high, hpd2_low, hpd2_hi, etc
#The reason for multiple hpd limits is if you have multimodal data.
#prob = prob *outside* the limit; i.e for a normal distribution 0.05 is expected to give
#     the 0.025 and 0.975 quantiles.
# this won't work for weighted data, use loc1statsx instead.
# xlim is optional - use it to define the support of the density.
loc1stats <- function(x,xlim,prob=0.05,...)
{
	if(missing(xlim)){
		fit <- locfit(~x)
	}
	else {
		fit <- locfit(~x,xlim=xlim)
	}
	fx <- fitted(fit)
	x.modef <- max(fx)
	x.mode <- x[fx == x.modef]
	if(!missing(xlim)){
		if(predict(fit,xlim[1]) > x.modef){
			x.modef <- predict(fit,xlim[1])
			x.mode <- xlim[1]
		}
		if(predict(fit,xlim[2]) > x.modef){
			x.modef <- predict(fit,xlim[2])
			x.mode <- xlim[2]
		}
	}

	if(length(x.mode)>1)x.mode <- x.mode[1]
	lev <- sort(fx)[floor(prob*length(x))]
#	print("log difference from max is ")
#	print(log(x.modef)-log(lev))
	l1 <- list()
	l1[[1]] <- x.mode
	indx <- order(x)
	ii <- 2
	flip <- TRUE
	for(j in 2:length(x)){
		if(flip && fx[indx[j]] > lev){
			l1[[ii]] <- x[indx[j-1]]
			if(j==2 && !missing(xlim)){
				if(predict(fit,xlim[1]) >= lev)l1[[ii]] <- xlim[1]
			}
			flip <- FALSE
			ii <- ii+1
		}
		else if(!flip && fx[indx[j]] < lev){
			l1[[ii]] <- x[indx[j]]
			flip <- TRUE
			ii <- ii+1
		}
		if(!flip && !missing(xlim) && j == length(x)){
			l1[[ii]] <- xlim[2]
			flip <- TRUE
		}
	}
	if(!flip)stop("HPD interval not closed")
	as.numeric(l1)
}



