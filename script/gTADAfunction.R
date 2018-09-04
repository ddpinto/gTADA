gTADA <- function(modelName  , inputData,
                         Ndn = NULL,
                         Ncase = NULL, Ncontrol = NULL,
                    geneSet = NULL,

                    nIteration = NULL, nIteration2 = NULL,
                    
                    nThin = NULL, nThin2 = NULL, nCore = 1, nChain = 1,
                         hyperBetaDN0 = NULL,
                         hyperBetaCC0 = NULL,
                         hyper2GammaMeanDN = NULL, hyper2BetaDN = NULL, ##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                         hyper2GammaMeanCC = NULL, hyper2BetaCC = NULL,
                  alpha01Mean = 0, alpha01SD = 2,

                         upperPi0 = 0.5, lowerPi0 = 0, lowerBeta = 0, ##Lower and upper limits of pi: should be default
                         lowerHyperGamma = 1, lowerGamma = 1, #Should be default
                         betaPars = c(6.7771073, -1.7950864, -0.2168248), #Adjust beta's values: should be default
                    adjustHyperBeta = as.integer(1), ##1 if want to adjust beta, 0 if don't want to adjust beta
                    autoAdjustHyperBeta = FALSE,
                    drawHeatMap = FALSE, 
                    writeResult = FALSE,
                  resultDir = NULL)
     {


###MCMC process
         geneName <- data.frame(inputData[, "Gene"])
         dataDN <- data.frame(inputData[, grep("dn_", colnames(inputData))])
#         colnames(dataDN) <- paste0("dn_", 1:dim(dataDN)[2])

         
         mutRate <- data.frame(inputData[, grep("mut_", colnames(inputData))])
 #        colnames(mutRate) <- paste0("mut_", 1:dim(mutRate)[2])
         dataCCcase <- data.frame(inputData[, grep("cc_case", colnames(inputData))])
  #       colnames(dataCCcase) <- paste0("cc_case", 1:dim(dataCCcase)[2])
         dataCCcontrol <- data.frame(inputData[, grep("cc_control", colnames(inputData))])
   #      colnames(dataCCcontrol) <- paste0("cc_control", 1:dim(dataCCcontrol)[2])

         if (!is.null(resultDir))
             resultDir <- "."
         if (dim(dataCCcontrol)[2] == 0)
             dataCCcontrol = NULL
         if (dim(dataCCcase)[2] == 0)
             dataCCcase = NULL
         if (dim(dataDN)[2] == 0)
             dataDN = NULL
         if (dim(mutRate)[2] == 0)
             mutRate = NULL

                  message("print(head(dataDN)) :", head(dataDN))
         
         message("MCMC is running")
         mcmcData <- gTADAmcmc(modelName = modelName,
                               geneSet = geneSet,
                                 dataDN = dataDN,
                               mutRate = mutRate, Ndn = Ndn,
                                 dataCCcase = dataCCcase, dataCCcontrol = dataCCcontrol, Ncase = Ncase, Ncontrol = Ncontrol,
                                 nIteration = nIteration, nIteration2 = nIteration2,
                               alpha01Mean = alpha01Mean, alpha01SD = alpha01SD,
                                 nThin = nThin, nThin2 = nThin2, nCore = nCore, nChain = nChain,
                                 hyperBetaDN0 = hyperBetaDN0, hyperBetaCC0 = hyperBetaCC0,
                                 hyper2GammaMeanDN = hyper2GammaMeanDN, hyper2BetaDN = hyper2BetaDN, ##Priors for mean RRs: meanRR ~ gamma(hyper2GammaMeanDN, hyper2BetaDN)
                                 hyper2GammaMeanCC = hyper2GammaMeanCC, hyper2BetaCC = hyper2BetaCC,
                                 upperPi0 = upperPi0, lowerPi0 = lowerPi0, lowerBeta = lowerBeta, ##Lower and upper limits of pi: should be default
                                 lowerHyperGamma = lowerHyperGamma, lowerGamma = lowerGamma, #Should be default
                                 betaPars = betaPars, #Adjust beta's values: should be default
                                 adjustHyperBeta = adjustHyperBeta, ##1 if want to adjust beta, 0 if don't want to adjust beta
                    autoAdjustHyperBeta =  autoAdjustHyperBeta)



###############Estimate genetic parameters
         message("\nEstimate genetic parameters from MCMC results")
         ############Add parameters
         parNames <- c("alpha0[1]")
         if (!is.null(geneSet))
             parNames <- paste0('alpha0[', 1:(dim(geneSet)[2]+1), ']')
         if (!is.null(dataCCcontrol)){
             parNames <- c(parNames, paste0('hyperGammaMeanCC[', 1:dim(dataCCcontrol)[2], ']'))
             parNames <- c(parNames, paste0('hyperBetaCC[', 1:dim(dataCCcontrol)[2], ']'))
         }
         if (!is.null(dataDN)){
             parNames <- c(parNames, paste0('hyperGammaMeanDN[', 1:dim(dataDN)[2], ']'))
             parNames <- c(parNames, paste0('hyperBetaDN[', 1:dim(dataDN)[2], ']'))
         }
         

         message("Estimate parameters from MCMC results\n")
         
         pars0 <- estimatePars(pars = parNames,
                     mcmcResult = mcmcData)
         print(pars0)
         pars1 <- as.numeric(pars0[, 1])
         names(pars1) <- rownames(pars0)

          parsFDR <- list(alpha0 = as.numeric(pars1[grep("alpha0", names(pars1))]),
             gammaMeanDN = as.numeric(pars1[grep("hyperGammaMeanDN", names(pars1))]),
                         betaDN = as.numeric(pars1[grep("hyperBetaDN", names(pars1))]),
                         gammaMeanCC = as.numeric(pars1[grep("hyperGammaMeanCC", names(pars1))]),
                         betaCC = as.numeric(pars1[grep("hyperBetaCC", names(pars1))]),

                         nfamily = Ndn,
                         ncase = Ncase,
                         ncontrol = Ncontrol
                         )

         message("\nCalculate posterior probabilities and FDRs")
         colnames(geneName) <- "Gene"
         dataOut <- calculateFDR(pars = parsFDR, geneSet = geneSet,
                                 dnData = dataDN, mutData = mutRate,
                                 caseData = dataCCcase, controlData = dataCCcontrol,
                                 geneName = geneName)

         dataFDR <- dataOut$dataFDR
         loglk <- dataOut$loglk
         message("\nDraw heatmaps")
         outTime <-  format(Sys.time(), "%a_%b_%d_%H_%M_%S_%Y")
         if (drawHeatMap) {
             pdf(paste0(resultDir, "/heatMap", outTime, ".pdf"))
                          allHyperGamma <- rownames(pars0[grep("hyperGammaMean", rownames(pars0)), ])
             for (i1 in 1:length(allHyperGamma))
                 plotParHeatmap(pars = c("alpha0[1]", allHyperGamma[i1]), mcmcResult = mcmcData)
             dev.off()
         }
         if (writeResult){
             write.table(dataFDR, paste0(resultDir, "/Result_extTADA_PosteriorAndFDR", outTime, ".txt"),
                         row.names = FALSE, quote = FALSE)
             write.table(pars0,   paste0("Result_extTADA_estimatedPars", outTime, ".txt"), quote = FALSE)
         }
########################
         message(paste0("\nThe analysis is completed.\nIf you want to analyse steps seperately, please take a look at the example in the manual"))

         pars0 <- cbind(rownames(pars0), pars0)
         colnames(pars0)[1] <- c("Parameter")
         return(list(dataPP = dataFDR, pars = pars0, gTADAmcmc = mcmcData, loglk = loglk))


          }

