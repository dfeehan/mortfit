#########################################################
## optimFit - a fitMethod object for fitting
## mortality hazards using R's function optim
##

## first, a new extension of mortalityFit with the extra
## info provided by an optim fit...
setClass("mortalityFitOptim",
         contains="mortalityFit",
         representation(optim.out="list",
                        start.gradient="numeric",
                        opt.gradient="numeric",
                        opt.hessian="matrix",
                        ## these have the eigenvalues of the Hessian...
                        opt.eigenvalues="numeric"))

#############################################
##' optim.fit
##'
##' use optim to fit a mortality model
##'
##' @param model.obj a mortalityModel object
##' @param data a mortalityData object
##' @param theta.init if not NULL, starting values to use for theta
##' @param random.start if TRUE, and the model.obj has
##'                     a non-NULL theta.range, then
##'                     draw a random starting value
##'                     uniformly from the values specified
##'                     by theta.range
##' @param verbose if TRUE, print extra info
##' @param ... TODO 
##' @return a mortalityFitOptim object
##' @export
optim.fit <- function(model.obj,
                      data,
                      theta.init=NULL,
                      random.start=FALSE,
                      verbose=FALSE,
                      ...)
{

  ## get initial parameters, unless they were
  ## explicitly specified
  if (is.null(theta.init)) {

    ## use the defaults by default...
    theta.init <- model.obj@theta.default
    
    ## if no function for choosing starting parameters was specified...
    ##if (is.null(model.obj@theta.start.fn)) {

    ## ... unless we're choosing random starting values, in which case
    ## we should do so
    if (random.start) {      

       if (! is.null(model.obj@theta.range)) {


          testhaz <- rep(0, length(data@data$age))
          numdraw <- 1

          ## draw random starting values; keep doing this until the initial
          ## hazards are all positive (for some hazards, like lynch-brown, 
          ## it's pretty easy to get random starting params that produce negative hazards)
          while (! all(testhaz > 0)) {

              if(numdraw > 50) {
                  stop("Failed to get starting values after 50 tries.")
              }

              if(verbose) {
                cat("randomly drawing starting values ", numdraw, "...\n")
              }
              numdraw <- numdraw + 1

              theta.init <- laply(model.obj@theta.range,
                                  function(tr) {
                                    if (tr[1] < tr[2]) {
                                      return(runif(1, min=tr[1], max=tr[2]))
                                    } else {
                                      return(runif(1, min=tr[2], max=tr[1]))
                                    }
                                  })
              testhaz <- model.obj@hazard@haz.fn(theta.init, data@data$age)
          }
        } else {
          if(verbose) {
            cat("no parameter ranges for this model, using defaults (not random).\n")
          }
        }

    ## otherwise, if we do have a function for choosing starting parameters,
    ## then do that. (note that this function takes a mortalityData object
    ## as an argument)
    } else if (! is.null(model.obj@theta.start.fn)) {

      if(verbose) {
        cat("calling model fn for starting values...\n")
      }
      theta.init <- model.obj@theta.start.fn(data)
    }

  }

  if(verbose) {
    cat("starting values: \n")
    cat(theta.init, "\n")
  }

  ## compute the gradient of the objective function at the
  ## starting values; this will be a useful comparison point later,
  ## when we want to figure out whether or not this has converged
  ## NB: this requires the numDeriv() library

  out <- try(start.gradient <- numDeriv::grad(func=model.obj@loglik.fn,
                                              x=theta.init,
                                              Dx=data@data$Dx,
                                              Nx=data@data$Nx,
                                              ages=data@data$age,
                                              method="simple"))
  if (class(out)=="try-error") {
    if(verbose) {
      cat("Error in computing gradient at starting values!\n")
      cat("theta = ", paste(theta.init, "\n"))
    }
    #browser()
  }

  ## use the numerical gradient 
  #num_grad <- Curry(numDeriv::grad, func=model.obj@loglik.fn, method="simple")
  num_grad <- functional::Curry(numDeriv::grad, func=model.obj@loglik.fn, method="Richardson")

  out <- try(op.out <- optim(par=theta.init,
                   fn=model.obj@loglik.fn,
                   gr=num_grad,
                   Dx=data@data$Dx,
                   Nx=data@data$Nx,
                   ages=data@data$age,
                   method=model.obj@optim.default$method,
                   control=c(list(fnscale=-1,
                                  trace=verbose,
                                  maxit=10000),
                             model.obj@optim.default$control),
                    ## NB: doing the hessian using numDeriv
                   hessian=FALSE
                   ))

  if (class(out)=="try-error") {
    if(verbose) {
      cat("Error in running optim!\n")
    }
    #browser()
  }

  if (op.out$convergence != 0) {
    stop("optimFit's optim.fit: optimization did not converge!\ncalled with ", model.obj@name, " - ", data@name, "\nstarting values: ", theta.init, "\n")
  }

  ## NB: this requires the numDeriv() library
  ## NB: I've found the method 'simple' to be somewhat
  ## more robust than 'Richardson', though th latter sounds
  ## more rigorous. Since our main focus isn't on estimates
  ## of the gradient or the Hessian, I'm using the 'simple'
  ## approach for now.
  opt.gradient <- numDeriv::grad(func=model.obj@loglik.fn,
                                 x=op.out$par,
                                 Dx=data@data$Dx,
                                 Nx=data@data$Nx,
                                 ages=data@data$age,
                                 method="simple")

  opt.hessian <- numDeriv::hessian(func=model.obj@loglik.fn,
                                   x=op.out$par,
                                   Dx=data@data$Dx,
                                   Nx=data@data$Nx,
                                   ages=data@data$age)

  opt.eigenvalues <- as.numeric(NA)
  try(opt.eigenvalues <- eigen(opt.hessian)$values)

  this.fit <- new("mortalityFitOptim",
                  name=paste(data@name, "-",
                             model.obj@name, "-",
                             "fit via optim"),
                  model=model.obj,
                  data=data,
                  method=optimFit,
                  theta.init=theta.init,
                  theta.hat=op.out$par,
                  log.likelihood=op.out$value,
                  start.gradient=start.gradient,
                  optim.out=op.out,
                  opt.gradient=opt.gradient,
                  opt.hessian=opt.hessian,
                  opt.eigenvalues=opt.eigenvalues)
                             
                             
  
  ## TODO -- need to create a fit object to return...
  ##stop("not finished writing this function yet...")
  return(this.fit)

}

##' fitMethod object for fitting via optim
##' @name optimFit
##' @format fitMethod object
##' @export
optimFit <- new("fitMethod",
                name="optimFit",
                fit=optim.fit)
