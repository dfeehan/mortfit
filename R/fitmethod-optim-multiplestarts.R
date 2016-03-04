#########################################################
## optimMultipleFit
##
## a fitMethod object for fitting
## mortality hazards by using several different
## starting values for optim and picking the one
## with the best fit

##' @include fitmethod-optim.R

## first, a new extension of mortalityFit with the extra
## info provided by an optim fit...
setClass("mortalityFitOptimMultiple",
         contains="mortalityFit",
         ## the first entry in this list will always be the fit using
         ## the model defaults
         representation(optim.fits="list",
                        all.lls="numeric",
                        best.fit.idx="integer"))

#############################################
##' optim.multiple.fit
##'
##' fit optim using the hazard default thetas and
##' several random starts; pick the best result
##'
##' @param model.obj a mortalityModel object
##' @param data a mortalityData object
##' @param M the number of random starts (in addition to defaults)
##' @param verbose set to TRUE to print extra information while fitting (useful for debugging)
##' @param ignore.folded TODO
##' @param ... other params to pass to optimFit fit method
##' @return a mortalityFitOptimMultiple object
##' @export
optim.multiple.fit <- function(model.obj,
                               data,
                               M=5,
                               verbose=TRUE,
                               ignore.folded=FALSE,
                               random.start=TRUE,
                               ...)
{

  fits <- as.list(rep(NA, M + 1))

  ## step 1: produce a preliminary fit using optim
  fits[[1]] <- mort.fit(model.obj,
                        data,
                        method=optimFit,
                        ignore.folded=ignore.folded,
                        verbose=verbose,
                        ...)  

  if(verbose) {
    cat("Fit with default thetas optim finished; log-likelihood: ", 
        fits[[1]]@log.likelihood, "\n")
    cat("Now using random starting values...\n")
  }

  for(i in 2:(M+1)) {

      fits[[i]] <- mort.fit(model.obj,
                            data,
                            method=optimFit,
                            random.start=TRUE,
                            ignore.folded=ignore.folded,
                            verbose=verbose,
                            ...)  
      
      if(verbose) {
        cat("Fit number ", i, " with random start finished; log-likelihood: ", 
            fits[[i]]@log.likelihood, "\n")
        cat("Now using random starting values...\n")
      }

  }
  
  ## step 3: summarize the results and return

  all.lls <- laply(fits, function(x) { x@log.likelihood })

  best.fit.idx <- which(all.lls == max(all.lls))

  ## we take the min in case there are multiple fits with the biggest log-likelihood
  best.fit <- fits[[min(best.fit.idx)]]

  this.fit <- new("mortalityFitOptimMultiple",
                  name=paste(data@name, "-",
                             model.obj@name, "-",
                             "fit via multiple start optim"),
                  model=model.obj,
                  data=data,
                  method=optimMultipleFit,
                  theta.init=best.fit@theta.init,
                  ## TODO -- do we need this theta.start.fn slot?
                  theta.start.fn=NULL,
                  theta.hat=best.fit@theta.hat,
                  log.likelihood=best.fit@log.likelihood,
                  optim.fits=fits,
                  all.lls=all.lls,
                  best.fit.idx=best.fit.idx)
}
                        
##' fitMethod object for fitting via optim using default 
##' parameters and several random starts, then picking best result
##' @name optimMultipleFit
##' @format fitMethod object
##' @export
optimMultipleFit <- new("fitMethod",
                        name="optimMultipleFit",
                        fit=optim.multiple.fit)



