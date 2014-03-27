#########################################################
## combinedFit
##
## a fitMethod object for fitting
## mortality hazards that combines optim and
## gridSearch; it produces initial parameter estimates
## from the optim approach, and then refines them with
## a grid search

##' @include fitmethod-optim.R
##' @include fitmethod-gridSearch.R

## first, a new extension of mortalityFit with the extra
## info provided by an optim fit...
setClass("mortalityFitCombined",
         contains="mortalityFit",
         ### TODO -- THINK ABOUT THIS MORE...
         representation(optim.fit="mortalityFitOptim",
                        grid.fit="mortalityFitGridSearch",
                        ll.gain="numeric"))

#############################################
##' combined.fit
##'
##' use a combination of optim followed by a grid search
##' to fit a mortality model
##'
##' @param model.obj a mortalityModel object
##' @param data a mortalityData object
##' @param M TODO
##' @param delta TODO
##' @param grid.dim TODO
##' @param verbose TODO
##' @param ignore.folded TODO
##' @param ... other params to pass to optimFit and gridSearchFit's
##'            fit methods
##' @return a mortalityFitOptim object
##' @export
combined.fit <- function(model.obj,
                         data,
                         M=5,
                         delta=0.1,
                         grid.dim=15,
                         verbose=TRUE,
                         ignore.folded=FALSE,
                         ...)
{

  ### TODO -- params to incorporate
  ### (not sure why ... is not behaving as expected):
  ###   - ignore.folded
  ###   - verbose
  ###   - also need to handle errors more robustly below...
  
  ## step 1: produce a preliminary fit using optim
  prelim.res <- mort.fit(model.obj,
                         data,
                         method=optimFit,
                         ignore.folded=ignore.folded,
                         verbose=verbose,
                         ...)  

  if(verbose) {
    cat("Preliminary fit via optim finished; log-likelihood: ", prelim.res@log.likelihood, "\n")
    cat("Now refining with grid fit...\n")
  }
  
  ## step 2: refine the fit with grid search
  grid.res <- mort.fit(model.obj,
                       data,
                       method=gridSearchFit,
                       ignore.folded=ignore.folded,
                       theta.init=prelim.res@theta.hat,
                       M=M,
                       delta=delta,
                       grid.dim=grid.dim,
                       verbose=verbose)
  
  ## step 3: summarize the results and return

  ## improvement in log likelihood due to grid search
  ll.gain <- grid.res@log.likelihood - prelim.res@log.likelihood
  
  ## check to see that the results from the grid search are at
  ## least as good as the ones from the optimization routine

  this.fit <- new("mortalityFitCombined",
                  name=paste(data@name, "-",
                             model.obj@name, "-",
                             "fit via combination optim/grid search"),
                  model=model.obj,
                  data=data,
                  method=combinedFit,
                  theta.init=prelim.res@theta.init,
                  ## TODO -- do we need this theta.start.fn slot?
                  theta.start.fn=NULL,
                  theta.hat=grid.res@theta.hat,
                  log.likelihood=grid.res@log.likelihood,
                  optim.fit=prelim.res,
                  grid.fit=grid.res,
                  ll.gain=ll.gain)
}
                        
##' fitMethod object for fitting via combination of
##' optim and grid search
##' @name combinedFit
##' @format fitMethod object
##' @export
combinedFit <- new("fitMethod",
                name="combinedFit",
                fit=combined.fit)



