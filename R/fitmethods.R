#######################################################
##
## fitmethods.R
##
## working version of the code for the oldage
## project, re-factored to use s4 classes
##
## this file uses the classes defined in
## oldage-s4-classdefns.R to actually create
## fit method objects; for example, there will
## be one object for grid search and another
## for nelder-mead
##
## @include fitmethod-optim.R
## @include fitmethod-gridSearch.R

##########################################################
##' summarize.fit
##'
##' summarize the fit of a model to data,
##' create a fitSummary object with
##' the results, and store the fit
##' summary object in the fit object
##'
##' @param fit.obj the mortalityFit object to summarize
summarize.fit <- function(fit.obj) {

  data <- fit.obj@data@data
  
  obs.Mx <- data$Dx/(data$Nx-0.5*data$Dx)
  obs.Dx <- data$Dx
  obs.qx <- data$Dx/data$Nx

  fits <- fit.obj@model@predict.fn(fit.obj@theta.hat,
                                   data$Nx,
                                   data$age)
  
  fitted.qx <- fits$fitted.qx  
  fitted.Dx <- fits$fitted.Dx

  # the number of people starting in the earliest age group is the
  # cohort size
  num.obs <- data$Nx[1]

  this.AIC <- (-2*fit.obj@log.likelihood+2*fit.obj@model@num.param)
  this.BIC <- (-2*fit.obj@log.likelihood+2+log(num.obs)*
                                               fit.obj@model@num.param)
  #this.BIC <- (-2*fit.obj@log.likelihood+2+log(nrow(data))*
  #                                             fit.obj@model@num.param)
  this.SSE.Dx <- sum( obs.Dx*((1-fitted.qx)^2) +
                     (data$Nx-obs.Dx)*(fitted.qx^2))


  ## TODO -- add RMSE?
  
  this.modelspec <- fit.obj@model@eval.fn(fit.obj,fitted.Dx,fitted.qx)
  
  this.fit.summary <- new("fitSummary",
                          AIC=this.AIC,
                          BIC=this.BIC,
                          SSE.Dx=this.SSE.Dx,
                          deviance=this.modelspec$deviance,
                          chisq=this.modelspec$chisq,
                          chisq.p=this.modelspec$chisq.p,
                          chisq.byage=this.modelspec$chisq.byage,
                          chisq.byage.p=this.modelspec$chisq.byage.p
                          ## TODO -- how to handle
                          ## optional cross validation?
                          ##CV5.Dx=this.CV5.Dx,
                          )
  
  fit.obj@fit.summary <- this.fit.summary

  return(fit.obj)
}

#########################################
##' summarize.prediction
##'
##' summarize the accuracy of predictions, by which
##' we mean predictions from a model to data where
##' we know the true value (though it doesn't
##' necessarily have tobe the data that we used
##' to fit the model)
##'
##' @param pred.obj a mortalityPrediction object
##' @param known a mortalityData object with the
##'              known values to compare the predictions
##'              to
summarize.prediction <- function(pred.obj, known) {

  data <- known@data
  
  obs.Mx <- data$Dx/(data$Nx-0.5*data$Dx)
  obs.Dx <- data$Dx
  obs.qx <- data$Dx/data$Nx

  SSE.Dx <- sum((pred.obj@fitted.Dx - obs.Dx)^2)
  RMSE.Dx <- sqrt( mean((pred.obj@fitted.Dx - obs.Dx)^2))

  return(list(SSE.Dx=SSE.Dx,
              RMSE.Dx=RMSE.Dx))
  
}



