#################################################################
## gridSearch - a fitMethod object and functions for
## performing grid search

## first, a new extension of mortalityFit with the extra
## info provided by a grid search fit...
setClass("mortalityFitGridSearch",
         contains="mortalityFit",
         representation(grid.search.out="list"))

#########################################
##' grid.search
##'
##' find the max (or min) of a function
##' by brute-force search through a grid
##' TODO
##'   - handle multiple maxima
##'
##' @param grid.vals either a list whose entries are
##'             vectors w/ parameter ranges, or
##'             a table of parameter values to
##'             evaluate
##' @param func the function to maximize
##' @param return.grid if TRUE, return the grid used
##' @param cpp if TRUE, use the c++ code; otherwise,
##'            use R's apply()
##' @param  ... other args, to be passed
##'             along to func
grid.search <- function(grid.vals,
                        func,
                        return.grid=FALSE,
                        #env,
                        cpp=FALSE,
                        verbose=FALSE,
                        ...)
{


    if (is.list(grid.vals)) {
      if(verbose) {
        cat("... building grid...")
      }
      grid.vals <- expand.grid(grid.vals)
      if(verbose) {
        cat("done.\n")
      }
    }

    if(verbose) {
      cat("... evaluating function on grid...")
    }

    if (cpp) {
      stop("c++ version not yet implemented...\n")
      ##if (missing(env)) {
      ##  env <- new.env()
      ##}
      ##vals <- .Call("grid_search_cpp", as.matrix(grid.vals), func, env, ...)
    } else {
      ## TODO -- the seg faults (or, one instance of them) seem to be coming
      ## from here...      
      vals <- apply(grid.vals, 1, func, ...)      
      ##browser()
      ##vals <- aaply(grid.vals,
      ##              1,
      ##              func,
      ##              ...)      
    }
    
    if(verbose) {
      cat("done.\n")
    }

    grid.vals$val <- vals

    ## handle duplicate maxima
    maxidx <- which(grid.vals$val==max(grid.vals$val,na.rm=TRUE))
    ##maxidx <- which.max(grid.vals$val)

    ## TODO -- changing this
    ##   -> also need to change non .cpp version
    ##   -> in calling function, need to detect if single or multiple maxima
    ##   -> if multiple maxima, then figure out what to do in general, and especially
    ##      when one of them is on the edge of the grid
    
    maxval <- grid.vals[maxidx[1],ncol(grid.vals)]

    maxparams <- grid.vals[maxidx,-ncol(grid.vals)]
    
    if(return.grid) {
        return(list(max=maxval,       # actual max
                    max.idx=maxidx,    # rows in grid.vals corresponding to maxima
                    grid=grid.vals))   # the whole grid
    } else {
        return(maxvals)
    }

}


#########################################
##' optimize.grid
##'
##' maximize a function using grid search
##' TODO -- not finished documenting this
##'
##' @param grid.vals a list with an entry for
##'             each param, which has a vector
##'             of all of the values the param
##'             should take on
##' @param M the number of times to run
##'        grid search (refining results
##'        each time...)
##' @param delta the size of the window to construct
##'        grid over. if delta=0.1, make a grid from
##'        (0.9*alpha, 1.1*alpha)
##' @param func the function to optimize
##' @param epsilon the minimum difference to consider significant;
##'        if the difference between the current optimum and the previous
##'        one is less than epsilon (in absolute value), then stop
##' @param verbose if TRUE, print out info about optimization
##' @param cpp if TRUE, use the C++ version... (default is FALSE)
##' @param ... other arguments, which get passed along to \link{grid.search}
##' @return TODO
##'
optimize.grid <- function(grid.vals,
                          M=3,
                          delta=.1,
                          func,
                          epsilon=1e-6,
                          verbose=FALSE,
                          cpp=FALSE,
                          ...)
{

    ## initialize some of the variables we'll use below
    next.grid.vals <- grid.vals
    this.M <- 1
    this.delta <- delta
    ## to start, pick an extremely small number for the previous maximum
    last.max <- -1e50
    cur.max <- last.max
    
    glength <- length(grid.vals[[1]])

    if(verbose) {
      cat("optimize grid called with M =", M,
          "delta =", delta, "\n")
    }
    
    ## be sure that the grid parameter entries are in
    ## order, least to greatest
    for(i in 1:length(grid.vals)) {
      if(grid.vals[[i]][1] > grid.vals[[i]][glength]) {
        if(verbose) {
          cat("Reversing grid dimension ", i,
              " so that entries range from smallest to largest.\n")
        }
        grid.vals[[i]] <- rev(grid.vals[[i]])
      }
    }

    ## this vector will hold the sequence of parameter estimates
    ests <- c()

    itnum <- 0
    maxit <- 100
      
    ## optimize on next grid
    while (this.M <= M) {
      
      if(verbose) {
        cat("running grid search round ", this.M, "..., delta =", this.delta,
            "last.max =", last.max, "\n")
      }

      itnum <- itnum + 1

      ## to avoid looping forever, be sure we haven't exceeded the (pre-specified)
      ## maximum number of possible iterations
      if (itnum > maxit) {
        browser()
        stop("reached maximum possible number of iterations without finding a solution.\n")
      }

      ## call the workhorse function, grid.search, which takes a set of grid
      ## values and the function to optimize, and figures out which grid entry
      ## actually corresponds to the maximum of the function over the grid
      gres <- grid.search(next.grid.vals,
                          func,
                          return.grid=TRUE,
                          verbose=verbose,
                          cpp=cpp,
                          ...)

      last.max <- cur.max
      cur.max <- gres$max

      ## note that if there was more than one maximum,
      ## these.params will be a data.frame
      these.params <- gres$grid[gres$max.idx,-ncol(gres$grid)]

      ## figure out whether or not these params
      ## are on the edge of the grid; this will be
      ## true if the values we've found match the
      ## grid entries at index 1 or glength
      ## in any direction

      ## go through and match the estimates from this iteration
      ## with their indices in the grid
      these.idx <- mapply(FUN=match,
                          these.params,
                          next.grid.vals)

      ## if all optima are on the edge of the grid, translate
      ## and repeat this step
      on.edge <- aaply(matrix(these.idx,byrow=FALSE,ncol=ncol(gres$grid)-1),
                       1,
                       function(each.idx) {
                         if (any(each.idx==1) || any(each.idx==glength)) {
                           return(TRUE)
                         } else {
                           return(FALSE)
                         }
                       })


      ## if all of the optima are estimated to be at the edge of
      ## a grid dimension, we'll translate and repeat this step
      if (all(on.edge)) {

        if (verbose) {
          cat("Found optimum at edge of grid. Translating...\n")
          cat("Repeating step ", this.M, " with grid:\n")          
        }

        new.grid <- next.grid.vals

        ## TODO -- eventually think about this more: arbitrarily pick
        ## the first of the maxima (since all were found on
        ## grid's edge)
        if (! is.null(dim(these.idx))) {
          these.idx <- these.idx[1,]
        }
          
        ## go through each dimension and
        ##   i) figure out if this optimum was at the grid's edge
        ##  ii) if so, translate it over in the right direction
        for( i in 1:length(these.idx)) {
            
          if (these.idx[i] == 1 || these.idx[i] == glength) {

            new.grid[[i]] <- grid.around(next.grid.vals[[i]][these.idx[i]],
                                         delta=this.delta,
                                         numsteps=glength)
              
          }
          
          if(verbose) {
            cat(range(new.grid[[i]]),"\n")
          }
            
        }

        ## update next.grid.vals to have the translated
        ## values for the dimensions whose optimum had been found
        ## to be at the grid's edge
        next.grid.vals <- new.grid

        ## re-start this iteration through the loop
        next

      } else {
        ## pick  the optima that is not on the edge that is closest
        ## to the current estimate
        ## TODO -- THINK THROUGH THIS!
        if (! is.null(dim(these.params)) && nrow(these.params) > 1) {
          if(verbose) {
            cat("picking one of multiple optima...\n")
          }
          these.params <- (these.params[!on.edge,])[1,,drop=FALSE]
        }

      }
        
      ##ests[[this.M]] <- cur.max
      ests[[this.M]] <- c(these.params, cur.max)

      if (verbose) {
        cat("    found: ", paste(cur.max), "\n")
      }
      
      ## build up refined param list
      next.grid.vals <- lapply(these.params,
                               function(x) {

                                 return(grid.around(x,
                                                    delta=this.delta,
                                                    numsteps=glength))
                                 
                               })

      this.M <- this.M + 1
      this.delta <- this.delta/2

    }

    ### TODO -- return a warning if the optimized value
    ### is at the edge of the grid in some dimension
    ### (this would suggest that it may not have found an
    ###  optimum in the steps needed...)
    return(list(M=M,
                ests=ests,
                max=ests[[length(ests)]]))

}

##' helper function which constructs a vector grid around
##' the value x; the grid is guaranteed to include x in its middle
##'
##' @param x the value to center the grid around
##' @param delta the range of the grid (eg .1 for .9x to 1.1x)
##' @param numsteps the number of entries on the grid
##' @return the grid slice
##' @keywords internal
grid.around <- function(x, delta, numsteps) {
  
  halfseq1 <- seq(from=(1-delta)*x,
                  to=x,
                  length=ceiling((numsteps+1)/2))
  halfseq2 <- seq(from=x,
                  to=(1+delta*x),
                  length=ceiling((numsteps+1)/2))
  fullseq <- c(halfseq1,
               halfseq2[-1])

  return(fullseq)
}

#########################################
##' grid.fit
##'
##' use grid search to fit a model to data
##'
##' @param model.obj the model object
##' @param data a mortalityData object
##' @param verbose if TRUE, print details of
##'                fit to the screen
##' @param M the number of rounds of grid search to use (see \link{optimize.grid})
##' @param delta the delta parameter (see \link{optimize.grid})
##' @param grid.dim the size of each dimension of the grid (see \link{optimize.grid})
##' @param theta.init if not NULL, the values around which the grid should be centered.
##'                   if NULL, defaults will be used
##' @param cpp if TRUE, use C++ version of the grid search
##'        algorithm (defaults to FALSE)
grid.fit <- function(model.obj,
                     data,
                     verbose=FALSE,
                     M = 3,
                     delta = 0.1,
                     grid.dim = 20,
                     theta.init=NULL,
                     cpp=FALSE,
                     ...)
{

  if(verbose) {
    cat("grid.fit called with M =", M, "delta =", delta,
        "and grid.dim =", grid.dim, "\n")
  }
  
  this.M <- M
  this.delta <- delta

  if (is.null(theta.init)) {
    theta0 <- model.obj@theta.default
  } else {
    theta0 <- theta.init
  }
  
  these.grid.vals <- lapply(theta0,
                            function(x) {

                              return(grid.around(x,
                                                 delta=this.delta,
                                                 numsteps=grid.dim))
                        })

  
  this.res <- optimize.grid(these.grid.vals,
                            M=this.M,
                            delta=this.delta,
                            func=model.obj@loglik.fn,
                            Dx=data@data$Dx,
                            Nx=data@data$Nx,
                            ages=data@data$age,
                            verbose=verbose,
                            cpp=cpp)

  this.fit <- new("mortalityFitGridSearch",
                  name=paste(data@name, "-",
                             model.obj@name, "-",
                             "fit via grid search"),
                  model=model.obj,
                  data=data,
                  method=gridSearchFit,
                  theta.init=theta0,
                  theta.start.fn=NULL,
                  theta.hat=as.numeric(this.res$ests[[this.M]])[1:model.obj@num.param],
                  log.likelihood=as.numeric(this.res$ests[[this.M]][(model.obj@num.param+1)]),
                  grid.search.out=this.res)

  return(this.fit)
                            
}

##' fitMethod object for fitting via grid search
##' @name gridSearchFit
##' @format fitMethod object
##' @export
gridSearchFit <- new("fitMethod",
                     name="gridSearch",
                     fit=grid.fit)
