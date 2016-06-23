#########################################################
## mortalitymodels.R
##
## this file contains the code which assembles hazard
## objects and likelihoods together into mortalityModel
## objects.
## this file uses the classes defined in
## oldage-s4-classdefns.R to create
## mortalityModel objects
## in the first section, we create a bunch of
## mortalityHazard objects
## in the second section, we create generic likelihood
## functions (Binomial and Poisson)
## in the third section, we put these together to
## create mortalityModel objects
##
## @include mortfit-help.R
## @include mortalityhazard-consthaz.R
## @include mortalityhazard-eta.R
## @include mortalityhazard-exweibull.R
## @include mortalityhazard-gompertz.R
## @include mortalityhazard-kannisto.R
## @include mortalityhazard-logistic.R
## @include mortalityhazard-lynchbrown.R
## @include mortalityhazard-makeham.R
## @include mortalityhazard-quadratic.R
## @include mortalityhazard-weibull.R

###########################################################
##' createBinomialModel
##'
##' factory method for models that use
##' the binomial log-likelihood fn;
##' takes hazard being used as an input
##' NB: do not add 0.5 to the ages or
##' anything like that here...
##'
##' @param haz.obj the mortalityHazard object to use in
##'                creating a mortalityModel
##' @return a mortalityModel object
createBinomialModel <- function(haz.obj) {

    ## build the likelihood function, using the
    ## hazard from the hazard object
    ## note that here we want the ages to be the
    ## start of the age interval we are interested in
    this.loglik <- function(theta,
                            Dx, Nx, ages) {

      rawhaz <- haz.obj@haz.fn(theta, ages)

      if (any (is.na(rawhaz) | (! is.finite(rawhaz)))) {
        return(NA)
      }
      
      pi <- haz.obj@haz.to.prob.fn(theta, ages )

      res <- sum((Dx*log(pi) + (Nx-Dx)*log1p(-pi)))

      return(res)
    }


    ## make a function to draw data from this
    ## model, given parameter values and a number of
    ## people who start in each age group
    ## here, the ages passed in should be the
    ## start of the age interval
    this.simfn <- function(theta,
                           Nx,
                           age,
                           seed=NULL,
                           means=FALSE) {

      if(! is.null(seed)) {
        set.seed(seed)
      }

      ## start age at 0 for computational purposes...
      age.offset <- min(age)
      age <-  age - min(age) + 1

      pi <- haz.obj@haz.to.prob.fn(theta, age)

      if( means ) {
        Dx <- Nx*pi
      } else {
        Dx <- rbinom(n=age,
                     size=Nx,
                     prob=pi)
      }

      this.dat <- new("mortalityData",
                      name=paste("Simulated from Binomial - ",
                                 haz.obj@name, sep=""),
                      ## TODO -- be sure age.offset and interval width
                      ##         are both right, or make sense as
                      ##         defaults in a simulated dataset
                      age.offset=age.offset,
                      age.interval.width=1L,
                      data=data.frame(age=age,
                                      Nx=Nx,
                                      Dx=Dx))
      return(this.dat)

    }

    ## make a function to predict from this
    ## model, given parameter values and a number of
    ## people who start in each age group
    ## here, the ages passed in should be the
    ## start of the age interval    
    this.predfn <- function(theta,
                            Nx, age,
                            age.offset=79) {

      fitted.qx <- haz.obj@haz.to.prob.fn(theta,
                                          age)
      fitted.Dx <- fitted.qx*Nx

      return(new("mortalityPrediction",
                 ## TODO -- eventually, fill in a sensible name here
                 name=paste("TODO -- need to fill in names"),
                 fitted.Dx=fitted.Dx,
                 fitted.qx=fitted.qx,
                 theta=theta,
                 age=age,
                 age.offset=age.offset,
                 Nx=Nx))

    }

    ##function that gives model-specific evaluations
    ## of a fit (for example, deviance, chi-squared stats, etc)
    ## NB: this function is intended to be passed a mortalityFit object
    ## observed - a mortalityData object
    ## fittedValues - a mortalityFit object
    this.evalfit <- function(fit.obj, observed, fittedValues) {

      data <- observed@data

      ### be sure that the ages line up; if not,
      ### exit with an error
      if (! all(data$age==fittedValues@age)) {
        stop("Error evaluating fit - fitted ages and observed ages do not line up.\n")
      }

      obs.Mx <- data$Dx/(data$Nx-0.5*data$Dx)
      obs.Dx <- data$Dx
      obs.qx <- data$Dx/data$Nx

      fitted.Dx <- fittedValues@fitted.Dx
      fitted.qx <- fittedValues@fitted.qx

      this.AIC <- (-2*fit.obj@log.likelihood+2*fit.obj@model@num.param)
      this.BIC <- (-2*fit.obj@log.likelihood+2+log(nrow(data))*
                                               fit.obj@model@num.param)
      
      this.SSE.Dx <- sum( obs.Dx*((1-fitted.qx)^2) +
                     (data$Nx-obs.Dx)*(fitted.qx^2))
      ##this.SSE.Dx <- sum( obs.Dx*((1-fitted.qx)^2) +
      ##               (data$Nx-obs.Dx)*(fitted.qx^2))      

      ## this is german rodriguez's way to compute the deviance
      ## (see http://data.princeton.edu/wws509/notes/a2s4.html
      ##  and also mccullagh and nelder, 1989)
      deviance.gr <- sum( 2*
                         (data$Dx*log(data$Dx/fitted.Dx) +
                          (data$Nx-data$Dx)*log((data$Nx-data$Dx)/
                                                (data$Nx-fitted.Dx))))
      ## chi-square stats
      chisq.byage <- (data$Dx-fitted.Dx)^2/
                     (data$Nx*fitted.qx*(1-fitted.qx))
      chisq.byage.p <- 1-pchisq(chisq.byage, df=1)

      chisq <- sum(chisq.byage)
      chisq.p <- 1-pchisq(chisq, df=length(chisq.byage))

      return(new("fitSummary",
                 AIC=this.AIC,
                 BIC=this.BIC,
                 SSE.Dx=this.SSE.Dx,
                 chisq=chisq,
                 chisq.p=chisq.p,
                 chisq.byage=chisq.byage,
                 chisq.byage.p=chisq.byage.p,
                 deviance=deviance.gr))
    }

    ## if the hazard object has a function for
    ## chosing parameter starting values, then
    ## include it in our mortalityModel object;
    ## otherwise, don't
    if (! is.null(haz.obj@theta.start.fn)) {
      tifn <- haz.obj@theta.start.fn
    } else {
      tifn <- NULL
    }

    ## if the hazard object has a function for
    ## computing the binomial log-likelihood gradient, then
    ## include it in our mortalityModel object;
    ## otherwise, don't
    if (! is.null(haz.obj@binomial.grad.fn)) {
      bgfn <- haz.obj@binomial.grad.fn
    } else {
      bgfn <- NULL
    }

    
    this.obj <- new("mortalityModel",
                    name=paste("Binomial", "-", haz.obj@name),
                    loglik.fn=this.loglik,
                    binomial.grad.fn=bgfn,
                    num.param=haz.obj@num.param,
                    theta.default=haz.obj@theta.default,
                    theta.range=haz.obj@theta.range,
                    optim.default=haz.obj@optim.default,
                    theta.start.fn=tifn,
                    predict.fn=this.predfn,
                    simulate.fn=this.simfn,
                    eval.fn=this.evalfit,
                    hazard=haz.obj)

    return(this.obj)

}

###########################################################
##' factory method for models using the Poisson likelihood
##'
##' factory method for models that use
##' the binomial log-likelihood fn;
##' takes hazard being used as an input
##' NB: this function expects the ages
##'     to have 0.5 already added, if
##'     that is appropriate
##'
##' @param haz.obj the mortalityHazard object to use in
##'                creating a mortalityModel
##' @return a mortalityModel object
createPoissonModel <- function(haz.obj) {

  this.loglik <- function(theta,
                          Dx, Nx, ages) {

    ##mu <- haz.obj@haz.fn( theta,
    ##                      ages )
    ##lambda <- Nx*mu

    ## TODO -- think this through again later,
    ## to be sure this is right...
    mu <- haz.obj@haz.fn(theta,
                         ages + 0.5 )


    lambda <- (Nx-0.5*Dx)*mu

    res <- sum( Dx*log(lambda) - lambda )

    return(res)
  }

  ## also make a function to draw data from this
  ## model, given parameter values and a number of
  ## people who start in each age group
  ## if means is TRUE, then return the expected
  ## values rather than draws from the distn
  this.simfn <- function(theta,
                         Nx, age,
                         seed=NULL,
                         means=FALSE) {

    if(! is.null(seed)) {
      set.seed(seed)
    }

    ## start at 0 for computational purposes...
    age.offset <- min(age)
    age <- age - min(age)

    mu <- haz.obj@haz.fn( theta, age + 0.5 )

    ## TODO -- does this have to be exposure?
    ## think about this...
    lambda <- Nx*mu

    ## TODO -- should these Dx come from probs instead
    ##         of from exposure?
    if(means) {
      Dx <- lambda
    } else {
      Dx <- rpois(length(lambda), lambda)
    }

    this.dat <- new("mortalityData",
                    name=paste("Simulated from Poisson - ",
                               haz.obj@name, sep=""),
                    ## TODO -- be sure age.offset and interval width
                    ##         are both right...
                    age.offset=age.offset,
                    age.interval.width=1L,
                    data=data.frame(age=age,
                                    Nx=Nx,
                                    Dx=Dx))
    return(this.dat)

  }

  ## and make a function that will predict using
  ## fitted parameters
  this.predfn <- function(theta,
                          Nx, age, age.offset=79) {

    ## TODO -- think about this: should I be adding 0.5?
    ## OR should we compute the average hazard, maybe?
    ##mu <- haz.obj@haz.fn( theta,
    ##                      age )

    fitted.qx <- haz.obj@haz.to.prob.fn(theta,
                                        age)
    fitted.Dx <- fitted.qx*Nx

    return(new("mortalityPrediction",
               name="TODO - eventually fill this in!",
               fitted.Dx=fitted.Dx,
               fitted.qx=fitted.qx,
               theta=theta,
               age=age,
               age.offset=age.offset,
               Nx=Nx))
  }

  ## and a function that gives model-specific evaluations
  ## of a fit (for example, deviance, chi-squared stats, etc)
  ## NB: this function is intended to be passed a mortalityFit object
  this.evalfit <- function(fit.obj, observed, fittedValues) {

    ## NB: observed is the observed mortalityData
    ##     fitted is the mortalityPrediction object we're evaluating
    ##     against the truth
    data <- observed@data

    ### be sure that the ages line up; if not,
    ### exit with an error
    if (! all(data$age==fittedValues@age)) {
      stop("Error evaluating fit - fitted ages and observed ages do not line up.\n")
    }

    obs.Mx <- data$Dx/(data$Nx-0.5*data$Dx)
    obs.Dx <- data$Dx
    obs.qx <- data$Dx/data$Nx

    fitted.Dx <- fittedValues@fitted.Dx
    fitted.qx <- fittedValues@fitted.qx

    this.AIC <- (-2*fit.obj@log.likelihood+2*fit.obj@model@num.param)
    this.BIC <- (-2*fit.obj@log.likelihood+2+log(nrow(data))*
                                             fit.obj@model@num.param)
    this.SSE.Dx <- sum( obs.Dx*((1-fitted.qx)^2) +
                       (data$Nx-obs.Dx)*(fitted.qx^2))


    ## this way of computing the deviance is from German's notes
    ## TODO -- go through and think about this again. would it be
    ##         better to explicity compute the saturated log-likelihood,
    ##         etc?
    deviance.gr <- sum(2*
                       (data$Dx*log(data$Dx/fitted.Dx)-
                       (data$Dx-fitted.Dx)))

    ## chi-square stats
    chisq.byage <- (data$Dx-fitted.Dx)^2/
                   (fitted.Dx)
    chisq.byage.p <- 1-pchisq(chisq.byage, df=1)

    chisq <- sum(chisq.byage)
    chisq.p <- 1-pchisq(chisq, df=length(chisq.byage))

    return(new("fitSummary",
               AIC=this.AIC,
               BIC=this.BIC,
               SSE.Dx=this.SSE.Dx,
               chisq=chisq,
               chisq.p=chisq.p,
               chisq.byage=chisq.byage,
               chisq.byage.p=chisq.byage.p,
               deviance=deviance.gr))
  }

  ## if the hazard object has a function for
  ## chosing parameter starting values, then
  ## include it in our mortalityModel object;
  ## otherwise, don't
  if (! is.null(haz.obj@theta.start.fn)) {
    tifn <- haz.obj@theta.start.fn
  } else {
    tifn <- NULL
  }
  
  this.obj <- new("mortalityModel",
                  name=paste("Poisson", "-",
                             haz.obj@name),
                  loglik.fn=this.loglik,
                  num.param=haz.obj@num.param,
                  theta.default=haz.obj@theta.default,
                  theta.range=haz.obj@theta.range,
                  optim.default=haz.obj@optim.default,
                  theta.start.fn=tifn,
                  predict.fn=this.predfn,
                  simulate.fn=this.simfn,
                  eval.fn=this.evalfit,
                  hazard=haz.obj)

  return(this.obj)
}

#########################################################
## SECTION 3 -- create mortalityModel objects

poissonWeibullModel <- createPoissonModel(weibull.haz)
binomialWeibullModel <- createBinomialModel(weibull.haz)

poissonExtendedWeibullModel <- createPoissonModel(exweibull.haz)
binomialExtendedWeibullModel <- createBinomialModel(exweibull.haz)

poissonConstantHazardModel <- createPoissonModel(consthaz.haz)
binomialConstantHazardModel <- createBinomialModel(consthaz.haz)

poissonQuadraticModel <- createPoissonModel(quad.haz)
binomialQuadraticModel <- createBinomialModel(quad.haz)

##TEMP - TEST C++ versions
##poissonQuadraticModelCpp <- createPoissonModel(quad.haz.cpp)
##binomialQuadraticModelCpp <- createBinomialModel(quad.haz.cpp)
##END TEMP

poissonLynchBrownModel <- createPoissonModel(lb.haz)
binomialLynchBrownModel <- createBinomialModel(lb.haz)

poissonMakehamModel <- createPoissonModel(mak.haz)
binomialMakehamModel <- createBinomialModel(mak.haz)

poissonGompertzModel <- createPoissonModel(gomp.haz)
binomialGompertzModel <- createBinomialModel(gomp.haz)

poissonLogisticModel <- createPoissonModel(logistic.haz)
binomialLogisticModel <- createBinomialModel(logistic.haz)

poissonKannistoModel <- createPoissonModel(kannisto.haz)
binomialKannistoModel <- createBinomialModel(kannisto.haz)

poissonETAModel <- createPoissonModel(eta.haz)
binomialETAModel <- createBinomialModel(eta.haz)

##' this is a list, provided for convenience, of all of
##' the models based on the poisson likelihood that are
##' automatically created
##' @name poisson.models
##' @format list
##' @export
poisson.models <- c(poissonWeibullModel,
                    ##poissonExtendedWeibullModel,
                    poissonConstantHazardModel,
                    poissonQuadraticModel,
                    ##poissonQuadraticModelCpp,
                    poissonLynchBrownModel,
                    poissonMakehamModel,
                    poissonGompertzModel,
                    poissonLogisticModel,
                    poissonKannistoModel
                    ##poissonETAModel
                    )

##' this is a list, provided for convenience, of all of
##' the models based on the binomial likelihood that are
##' automatically created
##' @name binomial.models
##' @format list
##' @export
binomial.models <- c(binomialWeibullModel,
                     ##binomialExtendedWeibullModel,
                     binomialConstantHazardModel,
                     binomialQuadraticModel,
                     ##binomialQuadraticModelCpp,
                     binomialLynchBrownModel,
                     binomialMakehamModel,
                     binomialGompertzModel,
                     binomialLogisticModel,
                     binomialKannistoModel
                     ##binomialETAModel
                     )


##' this is a list, provided for convenience, of all of
##' the models automatically created
##' @name all.models
##' @format list
##' @export
all.models <- c(binomial.models, poisson.models)
##all.models <- c(binomial.models)
names(all.models) <- unlist(plyr::llply(all.models, function(x) { x@name }))
