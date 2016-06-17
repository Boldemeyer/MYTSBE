
#' @title Multi-year time stratified Bayesian estimator (calendar year summary)

#' @description MYTSBE utilizes hierarchical Bayesian methods to estimate abundances from capture-mark-recapture (CMR) data using a temporally-stratified Lincoln-Petersen estimator. The between-year hierarchical
#' structures allows for annually recurring species characteristics to be incorporated into capture probabilities and abundance estimates.
#'
#' Summary output is by calendar year.
#'
#' Peer reviewed publication can be found at ...
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
#' @param trap.name character string used for titles and descriptions of reports
#' @param effort.cor expands the number of unmarked fish captured by a sampling effort
#' @param strata.length number of days in strata
#' @param smolt.juv.date "MM-DD" date to partition smolt life stage
#' @param rm.bad.GR remove U posterior distributions for strata with Gelman-Rubins test statistic >1.1 from yearly summary
#' @param strata.op.min minimum number of years data need to have been collected in a stratum to be included in summary
#' @param den.plot return density plots of MCMC chains
#' @param trace.plot return trace plots of MCMC chains
#' @import plyr
#' @importFrom lubridate yday
#' @importFrom R2jags jags
#' @import coda
#' @import lattice
#' @import superdiag
#' @import mcmcplots
#' @import ggmcmc
#' @importFrom data.table as.data.table
#' @export
#' @return NULL
#'
#'
MYTSBE_Calendar <- function(data,
                         effort.cor = FALSE,
                         sel.years = currentyear,
                         strata.op.min = 1,
                         smolt.juv.date = "06-01",
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
  smolt.date <- paste("2010-",smolt.juv.date, sep ="")


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

  writeLines(paste( Sys.Date(),"\n",
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
  quantiles$Strata<-NA
  for(s in 1:length(unique(data$strata))){
    quantiles$Strata<-ifelse(grepl(paste(s,",", sep = ""), quantiles$Parameter), min(data$strata-1)+s, quantiles$Strata)
  }

  #merge summary statistics and quantiles
  outputsummary <- merge(quantiles, means, by.x="Parameter") # Merge statistics
  is.num <- sapply(outputsummary, is.numeric) #reduce digits
  outputsummary[is.num] <- lapply(outputsummary[is.num], round, 3) #reduce digits
  outputsummary<-outputsummary[with(outputsummary, order(Year,Parameter)), ] #reorder by year and parameter

  ag.strata<- aggregate(ceiling(data$cor.factor), by=list(strata=data$strata), FUN=sum)
  excl.strata <-ag.strata[ag.strata$x < strata.op.min,]
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


  ###########################################################
  #                      Life Stage                         #
  ###########################################################

  smolt.strata <- ceiling(yday(as.Date(smolt.date))/strata.length) - (min(data$strata)-1) # convert to smolt cutoff date to strata

  for(selectyr in sel.years){
    syear <- (selectyr+1) - (as.numeric(min(data$year)))

    usep <- subset(model.fit.gg,grepl("^U", Parameter)) #subset all U's
    usep <- subset(usep,grepl(paste(",",syear,"]$" , sep = ""), Parameter)) #subset U's for first year
    usep$strata <- sub(".*?U(.*?)(,.*|$)", "\\1", usep$Parameter)
    usep$strata <- as.numeric(substr(usep$strata, 2, 5))

    juv<- subset(usep, strata >= (smolt.strata+1))

    juv<-subset(juv, !(juv$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold

    if(rm.bad.GR == TRUE){
      juv<-subset(juv, !(juv$Parameter %in% nc.gd))
    }

    juvdt<-as.data.table(juv)

    juvUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {juvUdist[i] <- sum(juvdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(juvdt$strata))
    }

    juvUdist<-as.data.frame(juvUdist) #change output to dataframe
    juvUdist$juvUdist<-as.numeric(juvUdist$juvUdist) #change output to numeric

    #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
    juvUsum<-adply(juvUdist, 2, summarise,
                    mode=names(which.max(table(juvUdist))),
                    mean = mean(juvUdist),sd = sd(juvUdist),
                    niaveSE = sd(juvUdist)/sqrt(length(juvUdist)))
    juvUsum$mode <- as.numeric(juvUsum$mode)

    #Get quantiles for U bootstrap distribution
    juvUquantiles<-adply(juvUdist, 2, function (x) quantile(x$juvUdist, c(0.025, 0.25, 0.5, 0.75, 0.975)))

    juvUoutputsummary <- merge(juvUquantiles, juvUsum, by.x="X1")
    juvUoutputsummary$X1 <- NULL #removed X1 variable
    juvUoutputsummary$Parameter<-paste("Juvenile Utot") # add Parameter variable
    juvUoutputsummary$Year<- selectyr #add year variable

    ################## SMOLT ##########################


    smolt<- subset(usep, strata <= smolt.strata)

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
    smoltUoutputsummary$Parameter<-paste("Smolt Utot") # add Parameter variable
    smoltUoutputsummary$Year<- selectyr #add year variable

    #setup bootstrap to randomly draw samples from each distribution in order to obtain total U statistics
    usepdt<-as.data.table(usep) # turn usep into a datable for bootstrap

    usepdt<-subset(usepdt, !(usepdt$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold
    if(rm.bad.GR == TRUE){
      usepdt<-subset(usepdt, !(usepdt$Parameter %in% nc.gd))
    }

    totUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {totUdist[i] <- sum(usepdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(usepdt$strata))
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
    totUoutputsummary$Parameter<-"Utot" # add Parameter variable
    totUoutputsummary$Year<- selectyr #add year variable

    year_int <- outputsummary[outputsummary$Year == selectyr ,]
    year_int <- year_int[complete.cases(year_int),]
    year_int <- subset(year_int,grepl("^U", Parameter)) #subset all U's

    cal.summary<-rbind.fill(totUoutputsummary, smoltUoutputsummary, juvUoutputsummary, year_int) #merge outputs

    cal.summary <- cal.summary[,c(10:12,1:9)]

    juv1<- subset(usep, strata >= (smolt.strata+1))
    smolt1<- subset(usep, strata <= smolt.strata)


    # bind usep.s and usep.p to find what strata were excluded
    if(rm.bad.GR == TRUE){
      juv.usep.not.used<-subset(juv1, (juv1$Parameter %in% nc.gd))
      juv.usep.not.used<-sort(unique(juv.usep.not.used$strata)) + (min(data$strata)-1)
      juv.usep.not.used<-paste(juv.usep.not.used, collapse= ", ")
      smolt.usep.not.used<-subset(smolt1, (smolt1$Parameter %in% nc.gd))
      smolt.usep.not.used<-sort(unique(smolt.usep.not.used$strata)) + (min(data$strata)-1)
      smolt.usep.not.used<-paste(smolt.usep.not.used, collapse=", ")
    }  else {
      juv.usep.not.used <- "none"
      smolt.usep.not.used <- "none"
    }

    #read out files
    sink(paste(selectyr, species, trap.name, "Calendar Results.txt"), append=FALSE, split=FALSE)
    writeLines(paste( Sys.Date(),"\n",
                      "\n",
                      selectyr, species, trap.name,"\n",
                      "Boot strap iterations =", boot,"\n",
                      "Specified smolt date:", smolt.juv.date, "\n",
                      "\n",
                      "Strata EXCLUDED from results due to <", strata.op.min, "year(s) of operation ","\n",
                      paste(x.strata, collapse=', '), "\n",
                      "\n",
                      "Strata EXCLUDED from results due to GR diagnostics >1.1  ", "\n",

                      "SMOLT (strata <= ", smolt.strata+(min(data$strata)-1) ,");",
                      smolt.usep.not.used,
                      "\n",
                      "JUVENILE (strata >", smolt.strata+1+(min(data$strata)-1),");",
                      juv.usep.not.used,
                      "\n",

                      "\n"
    )
    )
    print(cal.summary)
    sink()

    write.csv(cal.summary, file = paste(selectyr, species, trap.name, "Calendar Results.csv"))

    options(scipen = 0)
  }
  if(den.plot == TRUE) {
    ggmcmc(model.fit.gg, file="U Density Plots.pdf", family="U",plot="density")
    }
  if(trace.plot == TRUE) {
    ggmcmc(model.fit.gg, file="Trace Plot.pdf", plot="traceplot")
    }
}



