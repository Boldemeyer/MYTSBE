
#' @title Multi-year time stratified Bayesian estimator (cohort summary)

#' @description MYTSBE utilizes hierarchical Bayesian methods to estimate abundances from capture-mark-recapture (CMR) data using a temporally-stratified Lincoln-Petersen estimator. The between-year hierarchical
#' structures allows for annually recurring species characteristics to be incorporated into capture probabilities and abundance estimates.
#'
#' Summary output is by cohort year.
#'
#' Peer reviewed publication can be found at ... yadda yadda
#'
#' @param data capture-mark-recapture data frame
#' @param burnin number of initial MCMC chain iterations to be discarded
#' @param chains number of MCMC chains (>1)
#' @param iterations number of MCMC iterations per chain
#' @param thin thin rate
#' @param sel.years selected year(s) to compute abundance estimates
#' @param boot number of boot strap iterations to calculate yearly and life stage abundance estimates
#' @param species character string used for titles and descriptions of reports
#' @param model.params parameters to be returned from MCMC simulation
#' @param trap.name character string used for  titles and descriptions of reports
#' @param effort.cor expands the number of unmarked fish captured by a sampling effort
#' @param strata.length number of days in strata
#' @param smolt.parr.date "MM-DD" date to partition smolt life stage
#' @param parr.presmolt.date "MM-DD" date to partition parr life stage
#' @param rm.bad.GR remove U posterior distributions for strata with Gelman-Rubins test statistic >1.1 from yearly summary
#' @param strata.op.min minimum number of years data need to have been collected in a stratum to be included in summary
#' @param den.plot return density plots of MCMC chains
#' @param trace.plot return trace plots of MCMC chains
#' @import plyr
#' @importFrom lubridate hour
#' @importFrom R2jags jags
#' @import coda
#' @import lattice
#' @import superdiag
#' @import mcmcplots
#' @import ggmcmc
#' @importFrom data.table as.data.table
#' @export
#' @return NULL


MYTSBE_Cohort <- function(data,
                   effort.cor = FALSE,
                   sel.years = currentyear,
                   strata.op.min = 1,
                   smolt.parr.date = "07-01",
                   parr.presmolt.date = "09-01",
                   species = "",
                   trap.name = "",
                   den.plot = TRUE,
                   trace.plot = TRUE,
                   strata.length = s.length,
                   burnin = 100000,
                   chains = 3,
                   iterations = 400000,
                   thin = 100,
                   boot = 5000,
                   rm.bad.GR = FALSE,
                   model.params = c("p", "U", "etaP1", "etaU1", "sigmaU", "sigmaP")) {


  currentyear <- max(data$year)-2
  s.length <- max(data$days)
  smolt.date <- paste("2010-",smolt.parr.date, sep ="")
  parr.date <- paste("2010-",parr.presmolt.date, sep ="")

  model <- function() {
    for(i in 1:strata){
      for(j in 1:year){
        m[i,j] ~ dbin(p[i,j],n[i,j])
        logit(p[i,j]) <- etaP[i,j]
        etaP[i,j] ~ dnorm(etaP1[i],tauP)

        u[i,j] ~ dbin(p[i,j],U[i,j])
        U[i,j] <- round(exp(etaU[i,j]))
        etaU[i,j]~ dnorm(etaU1[i],tauU)

      }
      etaP1[i] ~ dnorm(np,tp)

      etaU1[i] ~ dnorm(nu,tu)
    }
    tauP ~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmaP <- 1/sqrt(tauP)
    tauU~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmaU <- 1/sqrt(tauU) #variance

    np ~ dnorm(-2,.666) #variance=1.5, precision=1/1.5->.666
    tp ~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmatp <- 1/sqrt(tp) #variance

    nu ~ dnorm(10,.25) #variance=4, precision=1/4=.25 , Increased mean to 10 due to initializing issues
    tu ~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmatu <- 1/sqrt(tu) #variance
  }

  #Input Data######################################################################################

  #U updated for effort correction factor
  data$cor.factor <-data$effort/data$days
  data$u.cor <- round(data$u*(1/(data$effort/data$days)))

  if(effort.cor == TRUE) {
    datau=matrix(round(data$u*(1/(data$effort/data$days))), nrow=length(unique(data$strata)), ncol=length(unique(data$year)))
  } else {
    datau=matrix(data$u,nrow=length(unique(data$strata)), ncol=length(unique(data$year)))
  }

  datan=matrix(data$n,nrow=length(unique(data$strata)), ncol=length(unique(data$year)))
  datam=matrix(data$m,nrow=length(unique(data$strata)), ncol=length(unique(data$year)))

  strata <- as.numeric(nrow(datan))# number of unique stratum
  year <- as.numeric(ncol(datan)) # number of unique years
  n    <- as.matrix(datan) # number of marked individuals available for recapture during the stratum
  m    <- as.matrix(datam) # number of marked individuals recaptured during the stratum
  u    <- as.matrix(datau) # number of unmarked individuals captured during the stratum

  options(max.print=100000)

  model.data <- list("strata", "year", "n", "m", "u")

  model.fit <- jags(data = model.data, inits = NULL,
                    parameters.to.save = model.params, n.chains = chains, n.iter = iterations,
                    n.burnin = burnin, n.thin = thin, model.file = model)

  #######################################################################################################################

  #Add in progress bar here!

  model.fit.mcmc <- as.mcmc(model.fit) #call model.fit from rjags
  model.fit.gg <- suppressWarnings(ggs(model.fit.mcmc)) #convert model.fit.mcmc to gg output (useable dataframe)

  options(scipen = 999) #change scientific notation

  #export data used for the model
  sink("Input Data.txt", append=FALSE, split=FALSE)
  print(data)
  sink()

  #export mcmc chains

  #need a for loop here############################################
  for(i in 1:chains){
    chain <- model.fit.gg[model.fit.gg$Chain == i,]
    write.table(chain, file = paste("MCMC Chains", i,".txt"), sep="\t")
  }

  #############
  #Diagnostics#
  #############
  gd <- gelman.diag(model.fit.mcmc, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
                    multivariate=FALSE)

  gdata <-as.data.frame(gd[1])
  gdata <- cbind(Parameter = rownames(gdata), gdata)
  colnames(gdata) <- c("Parameter", "Point.est", "Upper.CI")
  nonconv.gd <- gdata[gdata$Point.est > 1.1,]
  nonconv.gd <- droplevels(nonconv.gd)

  nc.gd<-as.character(nonconv.gd$Parameter)

  sink("GR Diagnostic Test Results.txt", append=FALSE, split=FALSE)

  writeLines(paste(Sys.Date(),"\n",
                   "\n",
                   "Number of chains =", chains,"\n",
                   "Number of iterations per chain =", iterations, "\n",
                   "Burn in =", burnin,"\n",
                   "Thin rate =", thin, "\n",
                   "\n",
                   "********** The Gelman-Rubin diagnostic: **********","\n","\n",
                   "\n",
                   "\n",
                   "   *  Flagged parameters with point est. >1.1  *  ", "\n",
                   "\n",
                   paste(nc.gd, collapse=', ' ),
                   "\n",
                   "\n"))
  print(gd)

  sink()

  ###################
  #Manual Statistics#
  ###################

  #apply descriptive statistic functions to parameter values over all chains & iterations
  means<-ddply(model.fit.gg, "Parameter", summarise,
               mode=names(which.max(table(value))),
               mean = mean(value),sd = sd(value),
               niaveSE = sd(value)/sqrt(length(value)))

  #get quantiles for each parameter
  means$mode<-as.numeric(means$mode)
  quantiles<-ddply(model.fit.gg, "Parameter", function (x) quantile(x$value, c(0.025, 0.25, 0.5, 0.75, 0.975)))

  #add year variable
  quantiles$Year<-NA
  for(j in 1:length(unique(data$year))){
    quantiles$Year<-ifelse(grepl(paste(",",j,"]$" , sep = ""), quantiles$Parameter), min(data$year-1)+j, quantiles$Year)
  }
  #add strata variable
  quantiles$strata<-NA
  for(s in 1:length(unique(data$strata))){
    quantiles$strata<-ifelse(grepl(paste(s,",", sep = ""), quantiles$Parameter), min(data$strata-1)+s, quantiles$strata)
  }

  #merge summary statistics and quantiles
  outputsummary <- merge(quantiles, means, by.x="Parameter") # Merge statistics
  is.num <- sapply(outputsummary, is.numeric) #reduce digits
  outputsummary[is.num] <- lapply(outputsummary[is.num], round, 3) #reduce digits
  outputsummary<-outputsummary[with(outputsummary, order(Year,Parameter)), ] #reorder by year and parameter

  ag.strata <- aggregate(ceiling(data$cor.factor), by=list(strata=data$strata), FUN=sum)
  excl.strata <- ag.strata[ag.strata$x < strata.op.min,]
  x.strata<-as.character(excl.strata$strata)

  outputsummary <- outputsummary[,c(7:8,1:6,9:12)]

  sink(paste("All",species,trap.name,"Results.txt"), append=FALSE, split=FALSE)
  writeLines(paste( Sys.Date(),"\n",
                    "\n",
                    species, trap.name,"\n",
                    "Number of chains =", chains,"\n",
                    "Number of iterations per chain =", iterations, "\n",
                    "Burn in =", burnin,"\n",
                    "Thin rate =", thin, "\n",
                    "* Strata with <", strata.op.min, "year(s) of operation  * ", "\n",
                    paste(x.strata, collapse=', ' ),
                    "\n",
                    "*  Flagged parameters with point est. >1.1  *  ", "\n",
                    paste(nc.gd, collapse=', ' ),
                    "\n"
  )
  )
  print(outputsummary)
  sink()

  write.csv(outputsummary, file = paste("All",species,trap.name,"Results.csv"))

  #########################################################################################
  #############                      Cohort Summary Loop                          #########
  #########################################################################################

  smolt.strata <- round(yday(as.Date(smolt.date))/strata.length) - (min(data$strata)-1)# convert to smolt cutoff date to strata
  parr.strata <- round(yday(as.Date(parr.date))/strata.length) - (min(data$strata)-1) # convert parr cutoff date to strata

  for(selectyr in sel.years){
    #remove U parameters for parr and presmolt
    pyear <- (selectyr+2) - (as.numeric(min(data$year)))

    usep <- subset(model.fit.gg,grepl("^U", Parameter)) #subset all U's
    usep.p <- subset(usep,grepl(paste(",",pyear,"]$" , sep = ""), Parameter)) #subset U's for first year
    usep.p$strata <- sub(".*?U(.*?)(,.*|$)", "\\1", usep.p$Parameter) #create strata variable
    usep.p$strata <- as.numeric(substr(usep.p$strata, 2, 5)) #clean strata variable

    parr<- subset(usep.p, strata >= (smolt.strata) & strata <= (parr.strata-1))

    parr<-subset(parr, !(parr$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold

    if(rm.bad.GR == TRUE){
      parr<-subset(parr, !(parr$Parameter %in% nc.gd))
    }

    parrdt<-as.data.table(parr)

    parrUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {parrUdist[i] <- sum(parrdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(parrdt$strata))
    }

    parrUdist<-as.data.frame(parrUdist) #change output to dataframe
    parrUdist$parrUdist<-as.numeric(parrUdist$parrUdist) #change output to numeric

    #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
    parrUsum<-adply(parrUdist, 2, summarise,
                    mode=names(which.max(table(parrUdist))),
                    mean = mean(parrUdist),sd = sd(parrUdist),
                    niaveSE = sd(parrUdist)/sqrt(length(parrUdist)))
    parrUsum$mode <- as.numeric(parrUsum$mode)

    #Get quantiles for U bootstrap distribution
    parrUquantiles<-adply(parrUdist, 2, function (x) quantile(x$parrUdist, c(0.025, 0.25, 0.5, 0.75, 0.975)))

    parrUoutputsummary <- merge(parrUquantiles, parrUsum, by.x="X1")
    parrUoutputsummary$X1 <- NULL #removed X1 variable
    parrUoutputsummary$Parameter<-paste("Parr U Cohort",selectyr) # add Parameter variable
    parrUoutputsummary$Year<- selectyr+1 #add year variable

    #########Pre smolt###############

    presmolt<- subset(usep.p, strata >= parr.strata)

    presmolt<-subset(presmolt, !(presmolt$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold

    if(rm.bad.GR == TRUE){
      presmolt<-subset(presmolt, !(presmolt$Parameter %in% nc.gd))
    }

    presmoltdt<-as.data.table(presmolt)

    presmoltUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {presmoltUdist[i] <- sum(presmoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(presmoltdt$strata))
    }

    presmoltUdist<-as.data.frame(presmoltUdist) #change output to dataframe
    presmoltUdist$presmoltUdist<-as.numeric(presmoltUdist$presmoltUdist) #change output to numeric

    #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
    presmoltUsum<-adply(presmoltUdist, 2, summarise,
                        mode=names(which.max(table(presmoltUdist))),
                        mean = mean(presmoltUdist),sd = sd(presmoltUdist),
                        niaveSE = sd(presmoltUdist)/sqrt(length(presmoltUdist)))
    presmoltUsum$mode <- as.numeric(presmoltUsum$mode)

    #Get quantiles for U bootstrap distribution
    presmoltUquantiles<-adply(presmoltUdist, 2, function (x) quantile(x$presmoltUdist, c(0.025, 0.25, 0.5, 0.75, 0.975)))

    presmoltUoutputsummary <- merge(presmoltUquantiles, presmoltUsum, by.x="X1")
    presmoltUoutputsummary$X1 <- NULL #removed X1 variable
    presmoltUoutputsummary$Parameter<-paste("Presmolt U Cohort",selectyr) # add Parameter variable
    presmoltUoutputsummary$Year<- selectyr+1 #add year variable

    ########## Smolt #############

    syear <- (selectyr+3) - (as.numeric(min(data$year)))

    usep.s <- subset(usep,grepl(paste(",",syear,"]$" , sep = ""), Parameter)) #subset U's for first year
    usep.s$strata <- sub(".*?U(.*?)(,.*|$)", "\\1", usep.s$Parameter) #create strata variable
    usep.s$strata <- as.numeric(substr(usep.s$strata, 2, 5)) #clean strata variable

    smolt<- subset(usep.s, strata < smolt.strata)

    smolt<-subset(smolt, !(smolt$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold

    if(rm.bad.GR == TRUE){
      smolt<-subset(smolt, !(smolt$Parameter %in% nc.gd))
    }

    smoltdt<-as.data.table(smolt)

    smoltUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {smoltUdist[i] <- sum(smoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(smoltdt$strata))
    }

    smoltUdist<-as.data.frame(smoltUdist) #change output to dataframe
    smoltUdist$smoltUdist<-as.numeric(smoltUdist$smoltUdist) #change output to numeric

    #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
    smoltUsum<-adply(smoltUdist, 2, summarise,
                     mode=names(which.max(table(smoltUdist))),
                     mean = mean(smoltUdist),sd = sd(smoltUdist),
                     niaveSE = sd(smoltUdist)/sqrt(length(smoltUdist)))
    smoltUsum$mode <- as.numeric(smoltUsum$mode)

    #Get quantiles for U bootstrap distribution
    smoltUquantiles<-adply(smoltUdist, 2, function (x) quantile(x$smoltUdist, c(0.025, 0.25, 0.5, 0.75, 0.975)))

    smoltUoutputsummary <- merge(smoltUquantiles, smoltUsum, by.x="X1")
    smoltUoutputsummary$X1 <- NULL #removed X1 variable
    smoltUoutputsummary$Parameter<-paste("Smolt U Cohort",selectyr) # add Parameter variable
    smoltUoutputsummary$Year<- selectyr+2 #add year variable

    #setup bootstrap to randomly draw samples from each distribution in order to obtain total U statistics

    totUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {totUdist[i] <- sum(parrdt[, value[sample.int(.N, 1, TRUE)], by = strata]) +
      sum(presmoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) +
      sum(smoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) -
      sum(unique(parrdt$strata)) -
      sum(unique(presmoltdt$strata)) -
      sum(unique(smoltdt$strata))
    }

    #Get summary statistics for U bootstrapped distribution
    totUdist<-as.data.frame(totUdist) #change output to dataframe
    totUdist$totUdist<-as.numeric(totUdist$totUdist) #change output to numeric

    #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
    totUsum<-adply(totUdist, 2, summarise,
                   mode=names(which.max(table(totUdist))),
                   mean = mean(totUdist),sd = sd(totUdist),
                   niaveSE = sd(totUdist)/sqrt(length(totUdist)))
    totUsum$mode <- as.numeric(totUsum$mode)

    #Get quantiles for U bootstrap distribution
    totUquantiles<-adply(totUdist, 2, function (x) quantile(x$totUdist, c(0.025, 0.25, 0.5, 0.75, 0.975)))

    #merge all summary statistics together
    totUoutputsummary <- merge(totUquantiles, totUsum, by.x="X1")
    totUoutputsummary$X1 <- NULL #removed X1 variable
    totUoutputsummary$Parameter<-paste("Total Cohort",selectyr) # add Parameter variable
    totUoutputsummary$Year<- NA #add year variable

    ######### Cohort strata used #########
    parr1<- subset(usep.p, strata >= (smolt.strata) & strata <= (parr.strata-1))
    presmolt1<- subset(usep.p, strata >= parr.strata)
    smolt1<- subset(usep.s, strata < smolt.strata)

    cohort.strata <- rbind(parr1,presmolt1,smolt1)

    cmeans<-ddply(cohort.strata, "Parameter", summarise,
                  mode=names(which.max(table(value))),
                  mean = mean(value),sd = sd(value),
                  niaveSE = sd(value)/sqrt(length(value)))

    #get quantiles for each parameter
    cmeans$mode<-as.numeric(cmeans$mode)
    cquantiles<-ddply(cohort.strata, "Parameter", function (x) quantile(x$value, c(0.025, 0.25, 0.5, 0.75, 0.975)))

    #add year variable
    cquantiles$Year<-NA
    for(j in 1:length(unique(data$year))){
      cquantiles$Year<-ifelse(grepl(paste(",",j,"]$" , sep = ""), cquantiles$Parameter), min(data$year-1)+j, cquantiles$Year)
    }

    cquantiles$strata<-NA
    for(s in 1:length(unique(data$strata))){
      cquantiles$strata<-ifelse(grepl(paste(s,",", sep = ""), cquantiles$Parameter), min(data$strata-1)+s, cquantiles$strata)
    }

    #merge summary statistics and quantiles
    coutputsummary <- merge(cquantiles, means, by.x="Parameter") # Merge statistics
    is.num <- sapply(coutputsummary, is.numeric) #reduce digits
    coutputsummary[is.num] <- lapply(coutputsummary[is.num], round, 3) #reduce digits
    coutputsummary<-coutputsummary[with(coutputsummary, order(Year,Parameter)), ] #reorder by year and parameter

    cohortsummary<-rbind.fill(totUoutputsummary,parrUoutputsummary,presmoltUoutputsummary, smoltUoutputsummary, coutputsummary) #merge outputs

    cohortsummary <- cohortsummary[,c(10:12,1:9)]

    # bind usep.s and usep.p to find what strata were excluded
    if(rm.bad.GR == TRUE){
      parr.usep.not.used<-subset(parr1, (parr1$Parameter %in% nc.gd))
      parr.usep.not.used<-sort(unique(parr.usep.not.used$strata)) + (min(data$strata)-1)
      parr.usep.not.used<-paste(parr.usep.not.used, collapse= ", ")
      presmolt.usep.not.used<-subset(presmolt1, (presmolt1$Parameter %in% nc.gd))
      presmolt.usep.not.used<-sort(unique(presmolt.usep.not.used$strata)) + (min(data$strata)-1)
      presmolt.usep.not.used<-paste(presmolt.usep.not.used, collapse=", ")
      smolt.usep.not.used<-subset(smolt1, (smolt1$Parameter %in% nc.gd))
      smolt.usep.not.used<-sort(unique(smolt.usep.not.used$strata)) + (min(data$strata)-1)
      smolt.usep.not.used<-paste(smolt.usep.not.used, collapse=", ")

    }  else {
      parr.usep.not.used <- "none"
      presmolt.usep.not.used <- "none"
      smolt.usep.not.used <- "none"
    }

    #read out files
    sink(paste(selectyr, species, trap.name, "Cohort Results.txt"), append=FALSE, split=FALSE)
    writeLines(paste( Sys.Date(),"\n",
                      "\n",
                      selectyr, species, trap.name,"\n",
                      "Boot strap iterations =", boot,"\n",
                      "Specified smolt date:", smolt.parr.date, "\n",
                      "Specified parr date:", parr.presmolt.date, "\n",
                      "\n",
                      "Strata EXCLUDED from results due to <", strata.op.min, "year(s) of operation ","\n",
                      paste(x.strata, collapse=', '), "\n",
                      "\n",
                      "Strata EXCLUDED from results due to GR diagnostics >1.1  ",
                      "\n",
                      "PARR (strata ",smolt.strata+(min(data$strata)-1),"-", parr.strata-1+(min(data$strata)-1),");",
                      parr.usep.not.used,
                      "\n",
                      "PRESMOLT (strata > =", parr.strata+(min(data$strata)-1),");",
                      presmolt.usep.not.used,
                      "\n",
                      "SMOLT (strata < = ", smolt.strata-1+(min(data$strata)-1) ,");",
                      smolt.usep.not.used,
                      "\n"
    )
    )
    print(cohortsummary)
    sink()

    write.csv(cohortsummary, file = paste(selectyr, species, trap.name, "Cohort Results.csv"))

    options(scipen = 0)

  }
  if(den.plot == TRUE) {
    ggmcmc(model.fit.gg, file="U Density Plots.pdf", family="U",plot="density")
  }
  if(trace.plot == TRUE) {
    ggmcmc(model.fit.gg, file="Trace Plot.pdf", plot="traceplot")
  }
}
