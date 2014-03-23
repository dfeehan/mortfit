##########################################################
##
## library.R
##
## this file has a bunch of methods that make use
## of the classes defined elsewhere in the mortfit
## package
##

#########################################
##' partition.into.folds
##'
##' given a set of indexes, partition it
##' into k folds
##'
##' @param num.folds the number of folds, k
##' @param data a mortalityData object
##' @return a mortalityDataFolded object
##' @export
partition.into.folds <- function(num.folds, data) {

  ##if (! require("rmnom")) {
  ##  stop("the rmnom package is required.\n")
  ##}

  cv.data <- c()

  fold.probs <- rep(1/num.folds, num.folds)

  mdat <- data@data
  
  mdat$Sx <- mdat$Nx-mdat$Dx

  fold.count.Dx <- t(sapply(mdat$Dx,
                            ##rmnom::rmnom,
                            rmnom,
                            matrix(fold.probs,ncol=num.folds)))
  fold.count.Sx <- t(sapply(mdat$Sx,
                            ##rmnom::rmnom,
                            rmnom,
                            matrix(fold.probs,ncol=num.folds)))

  fold.data <- list()
  fold.complement.data <- list()
  
  for(f in 1:num.folds) {

    this.fold.data <- data.frame(age=mdat$age,
                                 Dx=fold.count.Dx[,f],
                                 Nx=(fold.count.Dx[,f]+fold.count.Sx[,f]))

    ## only need to keep rows for which we've retained
    ## some number of observations
    this.fold.data <- subset(this.fold.data,
                             Nx > 0 )


    ## also get the complement of each fold
    this.fold.complement.data <- fold.complement(mdat,
                                                 this.fold.data)

    this.fold.complement.data <- subset(this.fold.complement.data,
                                        Nx > 0)

    this.fold.data <- new("mortalityData",
                          name=paste(data@name, " - fold",f),
                          age.offset=data@age.offset,
                          age.interval.width=data@age.interval.width,
                          data=this.fold.data,
                          tags=c(data@tags,fold=f))
    
    this.fold.complement.data <- new("mortalityData",
                                     name=paste(data@name, " - fold complement ",f),
                                     age.offset=data@age.offset,
                                     age.interval.width=data@age.interval.width,
                                     data=this.fold.complement.data,
                                     tags=c(data@tags, fold=f))

    fold.data[[f]] <- this.fold.data
    fold.complement.data[[f]] <- this.fold.complement.data

  }

  ## create a mortalityDataFolded object
  ## to return
  toret <- new("mortalityDataFolded",
               name=data@name,
               age.offset=data@age.offset,
               age.interval.width=data@age.interval.width,
               data=data@data,
               num.folds=as.integer(num.folds),
               folds=fold.data,
               fold.complements=fold.complement.data,
               tags=data@tags)

  return(toret)

}

#########################################
##' fold.complement
##'
##' return all of the obs in data that
##' aren't in the given fold
##'
##' @param data a data.frame with the entire dataset
##'             which was split into folds;
##'             must have columns 'Nx' and 'Dx'
##' @param fold a data.frame with the data for the fold
##'             we wish to complement; must have
##'             columns 'Nx' and 'Dx'
##' @return a data.frame with the complement of
##'           fold in data
fold.complement <- function(data, fold) {

    comp <- data

    fidx <- match(fold$age, data$age)

    comp[fidx,]$Nx <- comp[fidx,]$Nx-fold$Nx
    comp[fidx,]$Dx <- comp[fidx,]$Dx-fold$Dx

    return(comp)
}




#########################################
##' mort.fit
##'
##' use a fitMethod to fit a mortalityModel
##' to a set of mortalityData.
##' TODO -- more description
##'
##' @param model.obj a mortalityModel object with
##'           the mortality model to be fit
##' @param data a mortalityData object with the
##'           data to fit to
##' @param method a fitMethod object with the
##'           fit method to use
##' @param verbose if TRUE, print out extra details
##'           while the fit proceeds
##' @param ignore.folded if TRUE, then even if the
##'            data passed in are a mortalityDataFolded object,
##'            only run fits on the entire dataset (not each fold)
##' @param ... other parameters to pass to method
##' @return a mortalityFit object
##' @export
#########################################
mort.fit <- function(model.obj,
                     data,
                     method,
                     verbose=TRUE,
                     ignore.folded=FALSE,
                     ...)
{

  ## first, figure out whether the data we were
  ## passed has folds or not; if it does, then do
  ## mort.fit.cv
  if (is(data, "mortalityDataFolded") && (! ignore.folded)) {

      ## if this is called with num.folds
      ## as something other than NULL, then
      ## we want to run mort.fit.cv and
      ## return the results of that function...
      cl <- match.call()
      cl[[1]] <- as.name("mort.fit.cv")

      return(eval( cl, parent.frame() ))
  }

  if(verbose) {
    cat("================================================\n")
    cat("mort.fit: Fitting ", model.obj@name, " using ",
        method@name, " on ", data@name, ".\n")
  }

  this.fit <- method@fit(model.obj, data, verbose=verbose,
                         ignore.folded=ignore.folded, ...)

  if(! is(this.fit, "mortalityFit")) {
    stop("mortalityFit not returned for ", data@name, " - ", model.obj@name)
  }
  
  ## also add fitted values
  these.pred <- model.obj@predict.fn(this.fit@theta.hat,
                                     data@data$Nx,
                                     data@data$age)
  
  ## and also evaluate the fit
  this.eval <- model.obj@eval.fn(this.fit,
                                 data,
                                 these.pred)

  this.fit@fitted.values <- these.pred
  this.fit@fit.summary <- this.eval

  if(verbose) {
    cat("================================================\n")
  }
  
  return(this.fit)
  
}


#########################################
##
##' mort.fit.cv
##'
##' fit one of the functions to mortality data
##' using k-fold cross-validation
##' this function basically chooses k folds,
##' runs mort.fit on each one, computes the
##' error, and then re-runs the fit on the
##' entire data; finally, it returns the result.
##' it is designed to be called automatically
##' from mort.fit when that method is passed
##' a mortalityDataFolded object.
##'
##' @param model.obj the mortalityModel object
##' @param data a mortalityDataFolded object
##' @param method a fitMethod
##' @param verbose if TRUE, print extra info
##' @param keep.folds.fits if TRUE, retain all of the
##'              results from fitting each individual fold.
##'              this is useful for debugging, but can take
##'              up a lot of memory.
##' @param ... other arguments to mort.fit
##' @return a mortalityFit object
##' @export

#########################################
mort.fit.cv <- function(model.obj,
                        data,
                        method,
                        verbose=TRUE,
                        keep.folds.fits=FALSE,
                        ...)
{

  if(verbose) {
    cat("================================================\n")
    cat("mort.fit.cv: Fitting ", model.obj@name, " using ",
        method@name, " on ", data@name, ".\n")
    cat("starting values: \n")
    cat(model.obj@theta.default, "\n")
    cat("================================================\n")
  }
  
  num.folds <- data@num.folds

  cv.err.Mx <- list(rep(NA, num.folds))
  cv.err.Dx <- list(rep(NA, num.folds))
  cv.err.qx <- list(rep(NA, num.folds))

  if(keep.folds.fits) {
    folds.fits <- list(NA)
  }

  for(this.fold in 1:num.folds) {

    ## grab the next fold
    oos.data <- data@folds[[this.fold]]

    ## this is the data we'll use to fit the
    ## model (everything but the fold)
    fit.data <- data@fold.complements[[this.fold]]

    ## fit the model to everything but the held out data
    ## NB: this would be more elegant if we took the list of arguments,
    ## etc etc (TODO later)
    this.out <- try(mort.fit(model.obj,
                             fit.data,
                             method,
                             verbose=verbose,
                             ...))

    ## computed predicted values for the held out data
    these.pred <- model.obj@predict.fn(this.out@theta.hat,
                                       oos.data@data$Nx,
                                       oos.data@data$age)

    ## and summarize the prediction (we'll use the
    ## sum of squared errors in the predicted number of deaths)
    oos.fit <- summarize.prediction(these.pred, oos.data)

    cv.err.Dx[[this.fold]] <- oos.fit$SSE.Dx

    ## and finally, keep track of the fit for each fold so
    ## that we can look at it later...
    if(keep.folds.fits) {

      folds.fits[[this.fold]] <- c(this.out,
                                   list(fold.fit=this.out,
                                        fold.eval.fit=oos.fit,
                                        fold.err.Dx=cv.err.Dx[this.idx]
                                        ))
    }

  }
  
  ## now run on the entire dataset, by forcing ignore.folded to be TRUE
  ## and calling mort.fit...
  ubercall <- match.call(expand.dots=TRUE)
  ubercall[[1]] <- as.name("mort.fit")  
  ubercall[["ignore.folded"]] <- TRUE

  uber.fit <- eval(ubercall, parent.frame())
  
  ## also add fitted values
  these.pred <- model.obj@predict.fn(uber.fit@theta.hat,
                                     data@data$Nx,
                                     data@data$age)
  
  ## and also evaluate the fit
  this.eval <- model.obj@eval.fn(uber.fit,
                                 data,
                                 these.pred)

  uber.fit@fitted.values <- these.pred
  uber.fit@fit.summary <- this.eval
  
  ## now summarize CV results and add
  ## them to the fit summary object
  uber.fit@fit.summary <- as(uber.fit@fit.summary,
                             "fitSummaryCV")
  uber.fit@fit.summary@cv.err.Dx <- cv.err.Dx
  uber.fit@fit.summary@cv.rmse.Dx <- sqrt(mean(unlist(cv.err.Dx)^2))

  return(uber.fit)

}

####################################
##' create.mortalityFits
##'
##' given a list of mortalityFit objects,
##' create a mortalityFits object to
##' represent them. first, check to be sure
##' that all of the models were fit to the same dataset
##'
##' @param fits the list of mortalityFit objects
##' @return a mortalityFits object
##' @export
create.mortalityFits <- function( fits ) {

  ## first, double check that the datasets that were fitted
  ## to are all the same
  ubertmp <- length(unique(laply(fits, function(x) {
                                          paste(x@data@name)
                                        })))
  if (ubertmp !=1) {
    stop("Cannot make fits into a mortalityFits object. The models do not all appear to have been fitted to the same dataset.\n")
  }

  listoffits <- fits
  names(listoffits) <- laply(fits,
                             function(x) { x@name })
  
  ## then, create a new mortalityFits object to hold all of
  ## the fits in one place
  res <- new("mortalityFits",
             name=paste("fits to", fits[[1]]@data@name),
             fits=listoffits,
             data=fits[[1]]@data)

  return(res)
  
}

####################################
##' convert from hazards to probs
##'
##' this function takes a hazard function
##' and an age range and returns an approximation
##' of the probability of death implied by the
##' hazard over that age range. it does this either
##' through approximations or through numerical
##' integration.
##' TODO - add more details here
##' NB: this function expects ages
##' to be the starting age (no need
##' to add 0.5 or anything)
##' TODO - document warning for integration try-error
##' TODO - tests here - not too hard to think of...
##'
##' @param haz.fn the hazard function
##' @param theta parameter values for the hazard function
##' @param ages a vector of ages with entries (z1,z2,...,zk). we'll return
##'             probabilities of survival for (z1,z2), (z2,z3), ..., (z(k-1),zk)
##' @param nozero if TRUE, replace values that would be returned as zero with very small numbers instead. (this prevents errors when taking logs.)
##' @param approx if TRUE, use analytical approximation; if FALSE,
##'               use numerical integration
##' @return a vector of probabilities
haz.to.prob <- function( haz.fn,
                         theta,
                         ages,
                         nozero=TRUE,
                         approx=FALSE )
{

   ## get the hazards for the ages
   ## and then convert to a prob with exp(-hazard)
   ## TODO -- should this use ages+0.5? or something other
   ##         than the left endpoint? look at preston et al
   if(approx) {
      this.q <- 1 - exp( -1 * haz.fn( theta, ages + 0.5) )
   } else {

      ## in order to convert hazards to probabilities over
      ## an age range (z,z+1), we will want to integrate this function,
      ## hazard(x)
      toint <- function(z) {
        haz.fn(theta,z)
      }

      ##this.p <- sapply(ages,
      this.p <-  apply(as.matrix(ages),
                       1,
                       function(z) {
                         ## check to be sure that all of the ages yield finite
                         ## values; otherwise, the integral won't converge
                         ## (note that this would have to be Inf and never -Inf, since
                         ##  hazards are always positive)
                         if (! is.finite(toint(z+1))) {
                           return(Inf)
                         }
                         
                         out <- try(integrate(toint,z,z+1,
                                              subdivisions=1000,
                                              stop.on.error=TRUE),
                                    silent=TRUE)
                         if (class(out) == "try-error") {
                           
                           tmp <- as.character(out)

                           ## check to see if this error comes up b/c the function
                           ## we're integrating isn't finite. if so, we can return Inf
                           ## (this will happen because some of the optimization algorithms
                           ##  wander into ridiculously implausible regions of parameter
                           ##  space)
                           ##if (grepl("non-finite function value", tmp)) {
                           ##  return(Inf)
                           ##}

                           msg <- paste("try-error in integrating to conver haz to prob...\n",
                                        out,"\n")
                           warning(msg)

                           ## TODO -- put this back in in a sec...
                           ##browser()

                           return(NA)                           
                         }

                         return(exp(-out$value))
                       })

      this.q <- 1 - this.p
  }

   ## replace values that are 1 or 0 with very close
   ## approximations (since they'll mess up logs)
   if (nozero) {
     tinyval <- .Machine$double.eps
     this.q[this.q<=tinyval] <- tinyval
   }

   return( this.q )
}


######################################################
##' plot likelihood profile for all parameters
##'
##' given a fitted mortalityModel, plot the profile
##' of the log-likelihood function for all pairs
##' of parameters. this can take a while if there are
##' many parameters
##'
##' @param fit a mortalityFit object
##' @param delta window around the fitted max to profile. if delta
##'              is 0.1, then plot values for (0.9*max, 1.1*max), where
##'              max is the fitted maximum parameter value
##' @param res the resolution, or size of the grid on which the
##'            likelihood is evaluated. large grids may take a long time
##' @param heatmap if TRUE, return a heatmap for each pair of parameters;
##'                otherwise, return a contour plot
##' @param show if TRUE, actually plot the list of plots (as a side effect),
##'             in addition to returning it
##' @param ... other arguments TODO
##' @return a list whose entries contain ggplot2 plots of
##'         the profiles
plot.ll.prof.all <- function(fit,
                             delta=.1,
                             res=100,
                             heatmap=FALSE,
                             show=FALSE,
                             ...) {

    num.param <- length(fit@theta.hat)

    ## get all pairs of indices of the parameter vector
    theta.pairs <- combn(1:num.param, 2)

    ## make a likelihood profile plot for each of these pairs,
    ## and return it
    pairwise.plots <- alply(theta.pairs,2,
                            function(pair) {
                              return(plot.ll.prof(fit,theta.idx=pair,
                                                  delta=delta,res=res,
                                                  heatmap=heatmap))
                            })

    if(show) {
      do.call("grid.arrange", pairwise.plots)
    }
    
    return(pairwise.plots)

}

######################################################
##' plot likelihood profile for a pair of parameters
##'
##' given a fitted mortalityModel, plot the profile
##' of the log-likelihood function for a given pair
##' of parameters. this can take a while if the likelihood
##' takes a while to evaluate and/or the resolution is
##' chosen to be large
##'
##' @param fit a mortalityFit object
##' @param theta.idx vector with the indexes of the two parameters to plot against
##'                  each other; for example, in a four-parameter model,
##'                  to examine the profile of the likelihood as
##'                  parameters one and three vary, pass in theta.idx=c(1,3)
##' @param delta window around the fitted max to profile. if delta
##'              is 0.1, then plot values for (0.9*max, 1.1*max), where
##'              max is the fitted maximum parameter value
##' @param res the resolution, or size of the grid on which the
##'            likelihood is evaluated. large grids may take a long time
##' @param heatmap if TRUE, return a heatmap for each pair of parameters;
##'                otherwise, return a contour plot
##' @return a ggplot2 plots of the profile
plot.ll.prof <- function(fit,
                         theta.idx=c(1:2),
                         delta=.1,
                         res=200,
                         heatmap=FALSE) {

    ## if asked for a profile of more
    ## than two params, send to
    ## plot.ll.prof.all
    if (length(theta.idx) != 2) {
        cl <- match.call()
        cl[[1]] <- as.name("plot.ll.prof.all")
        return(eval(cl, parent.frame()))
    }

    ## grab what we'll use for plotting the
    ## profile from the object...
    theta.hat <- fit@theta.hat
    this.dat <- fit@data@data
    this.lik <- fit@model@loglik.fn
    
    this.grid <- as.list(theta.hat)

    this.grid[[theta.idx[1]]] <- seq(from=(1-delta)*
                                          theta.hat[theta.idx[1]],
                                     to=(1+delta)*
                                        theta.hat[theta.idx[1]],
                                       length=res)
    this.grid[[theta.idx[2]]] <- seq(from=(1-delta)*
                                          theta.hat[theta.idx[2]],
                                     to=(1+delta)*
                                        theta.hat[theta.idx[2]],
                                       length=res)

    this.grid <- expand.grid(this.grid)

    colnames(this.grid) <- paste("theta",1:ncol(this.grid),sep="")

    this.grid$ll <- apply(this.grid,
                          1,
                          function(theta) {
                              this.lik(theta=theta,
                                       Dx=this.dat$Dx,
                                       Nx=this.dat$Nx,
                                       ages=this.dat$age)
                          })

    ## TODO -- also check max against grid search here?

    plotgrid <- this.grid[,theta.idx]
    plotgrid$ll <- this.grid$ll

    pgn <- colnames(plotgrid)

    ll.plot <- ggplot(this.grid, aes_string(x=pgn[1],
                                            y=pgn[2],
                                            z=pgn[3])) +
               geom_contour(bins=40,
                            aes(colour=..level..)) +
               opts(legend.position="none")

    if (heatmap) {
        ll.plot <- ll.plot + geom_tile(aes(fill=ll)) +
                   scale_fill_gradient(low="red", high="green")
    }

    ll.plot <- ll.plot + geom_point(x=theta.hat[theta.idx[1]],
                                    y=theta.hat[theta.idx[2]],
                                    color="blue",size=3)

    return(ll.plot)


}

##############################################
##'
##' lre
##'
##' compute the log10 relative error in an
##' estimated versus a true parameter vector.
##' this is useful for unit tests to be sure
##' optimization routines are well-behaved
##' see Altman et al, "Numerical Issues in
##' Statistical Computing...", pg 54
##' TODO -- come back to this; j is not defined
##' in their expression (3.4)!
##'
##' @param theta.true the true parameter vector
##' @param theta.hat the estimated parameter vector
##' @param j scaling factor
##' @return a vector of the same length as theta.true
lre <- function(theta.true, theta.hat, j=1) {
  
  relerr <- abs(j*(theta.hat-theta.true)/(theta.true*j))

  relerr[theta.true==0] <- j
  
  return(-1*log10(relerr))
  
}

##############################################
##'
##' load.kt.data
##'
##' load the period and cohort mortality rates
##' from the kannisto-thatcher database
##' based on the years, cohorts, and countries
##' given in the passed-in file
##' TODO - the load.hmd files are not
##' yet documented
##'
##' @param kt.desc.file location of .csv file which describes
##'                   which country-years/cohorts to use
##' @param kt.dir the directory which has the data files for the
##'               kannisto-thatcher database
##' @param kt.format which of the files to use (for period data)
##'                   defaults to "1"; no foreseeable need to change this
##' @return an object containing (i) a two-dimensional list of data frames,
##'           indexed by (country, sex); (ii) country.years, a data frame with the
##'           countries and years in the object
##' @export
##############################################
load.kt.data <- function(kt.desc.file,
                         kt.dir,
                         kt.format="1" )
{

    ## load the description file
    touse <- read.csv(kt.desc.file)  

    ## now go through and pull all of
    ## the relevant files from kt database
    countries <- touse$country
    num.countries <- length(countries)

    #################################################
    ## grab the cohort data
    #################################################
    ##blank.md <- new("mortalityData", name="blank")

    kt.data <- list()
    
    for(this.sex in c("m","f")) {
        for(this.country in countries) {

            this.file <- paste(this.sex, this.country, ".txt", sep="")
            this.dat <- read.table(paste(kt.dir, "/", this.file, sep=""),
                                   header=TRUE, na.strings=".", sep=",",
                                   strip.white=TRUE)

            this.dat <- melt(this.dat,
                             id.vars=c("Cohort", "Age", "Year", "Triangle"))

            this.dat <- dcast(this.dat,
                              Cohort + Age ~ variable,
                              fun.aggregate=sum,
                              na.rm=TRUE)

            colnames(this.dat) <- c("cohort", "age", "Nx", "Dx")

            these.cohorts <- seq(from=touse[touse$country==this.country,
                                            "min.cohort"],
                                 to=touse[touse$country==this.country,
                                          "max.cohort"],
                                 by=1)


            ## finally, make a dataframe and stick it in the list
            for(this.cohort in these.cohorts) {

                this.cdat <- subset(this.dat, cohort==this.cohort)

                min.age <- min(this.cdat$age)
                this.cdat$age <- this.cdat$age - min.age + 1

                this.obj <- new("mortalityData",
                                name=paste(this.country, this.sex,
                                           "cohort", this.cohort, sep="-"),
                                age.offset=as.integer(min.age-1),
                                age.interval.width=1L,
                                data=as.data.frame(this.cdat),
                                tags=list(country=paste(this.country),
                                          sex=this.sex,
                                          type="cohort",
                                          time=paste(this.cohort)))

                kt.data <- c(kt.data, list(this.obj))

              }
          }                                        
    }

    #################################################
    ## grab the period data
    #################################################
    for(this.sex in c("m","f")) {
        for(this.country in countries) {

            this.file <- paste(this.sex, this.country, "_",
                               kt.format, ".txt", sep="")
            this.dat <- read.table(paste(kt.dir, "/", this.file, sep=""),
                                   header=TRUE, na.strings=".", sep=",",
                                   strip.white=TRUE)

            this.dat <- subset(this.dat,
                               select=c("FirstYear", "Age", "Dx", "Nx"))
            colnames(this.dat) <- c("year", "age", "Dx", "Nx")

            these.years <- seq(from=touse[touse$country==this.country,
                                          "min.year"],
                               to=touse[touse$country==this.country,
                                        "max.year"],
                               by=1)

            ## finally, make a dataframe and stick it in the list
            for(this.year in these.years) {

                this.pdat <- subset(this.dat, year==this.year)

                min.age <- min(this.pdat$age)
                this.pdat$age <- this.pdat$age - min.age + 1

                this.obj <- new("mortalityData",
                                name=paste(this.country, this.sex,
                                           "period", this.year, sep="-"),
                                age.offset=as.integer(min.age-1),
                                age.interval.width=1L,
                                data=as.data.frame(this.pdat),
                                tags=list(country=paste(this.country),
                                          sex=this.sex,
                                          type="period",
                                          time=paste(this.year)))
                
                kt.data <- c(kt.data, list(this.obj))                
            }
        }
    }

    ## the range of cohorts and of periods covered
    ## (across countries)
    cohorts <- seq(from=min(touse$min.cohort),to=max(touse$max.cohort))
    num.cohorts <- length(cohorts)
    years <- seq(from=min(touse$min.year), to=max(touse$max.year))
    num.years <- length(years)
    types <- c("cohort", "period")
    
    return(list(kt.data=kt.data,
                countries=countries,
                types=types,
                num.countries=num.countries,
                sexes=c("m","f"),
                cohorts=cohorts,
                num.cohorts=num.cohorts,
                years=years,
                num.years=num.years))
}

##############################################
##'
##' grab.kt.data
##'
##' this is a helper function which grabs
##' datasets from the kannisto-thatcher
##' database using their tags
##'
##' @param data the list containing the mortalityData
##'             objects from the kannisto-thatcher
##'             database
##' @param tags the tags to search for
##' @return a list of mortalityData objects from
##'         the database that match the tags
##'         passed into this function
##' @export
grab.kt.data <- function(data,
                         tags)
{

  res <- laply(data,
               function(x) {

                 ## if this entry doesn't have any tags
                 ## recorded, don't return it
                 if (! is.null(x@tags)) {

                   ## if this entry does have tags and
                   ## all of the ones specified in the
                   ## argument passed in match, then
                   ## return it
                   if(all(x@tags[names(tags)]==tags)) {
                     return(TRUE)

                   ## otherwise, don't return it...
                   } else {
                     return(FALSE)
                   }
                 } else {
                   return(FALSE)
                 }
               })

  return(data[res])
  
}
                         

##############################################
#
# load.hmd.Mx
#
# function to load mortality rates (Mx)
# from the human mortality database
#
# arguments
#      country - the short-form of the country name
#       cohort - if TRUE, cohort data; otherwise, period
#       format - the format to use (eg "1x1")
#      hmd.dir - the directory w/ the hmd data
# returns
#      a data frame with the mortality rates;
#   adds a numeric age variable and renames
#   Age to agelab
#
##############################################
load.hmd.Mx <- function( country, format, cohort=FALSE,
                         hmd.dir = "/Volumes/LaCie/data/hmd/hmd_countries") {

  chrt <- ""
  if (cohort) {
    chrt <- "c"
  }

  filename <- paste(hmd.dir, "/", country, "/STATS/", chrt, "Mx_", format, ".txt", sep="")

  this.dat <- read.table( filename, skip=1, header=TRUE, na.strings="." )

  # TODO -- this may need to be made more general...
  #  (in particular, is oldest age the same for everyone?)
  colnames(this.dat)[colnames(this.dat)=="Age"] <- "agelab"
  this.dat$numage <- as.numeric(paste(this.dat$agelab))
  maxage <- max(this.dat$numage, na.rm=TRUE) + 1
  this.dat$numage[is.na(this.dat$numage)] <- maxage

  return(this.dat)
}

##############################################
#
# load.hmd.Dx
#
# function to load death counts (Dx)
# from the human mortality database
#
# arguments
#      country - the short-form of the country name
#       format - the format to use (eg "1x1")
#      hmd.dir - the directory w/ the hmd data
# returns
#      a data frame with the mortality rates;
#   adds a numeric age variable and renames
#   Age to agelab
#
##############################################
load.hmd.Dx <- function( country, format,
                         hmd.dir = "/Volumes/LaCie/data/hmd/hmd_countries") {

  filename <- paste(hmd.dir, "/", country, "/STATS/Deaths_", format, ".txt", sep="")

  this.dat <- read.table( filename, skip=1, header=TRUE, na.strings="." )

  # TODO -- this may need to be made more general...
  #  (in particular, is oldest age the same for everyone?)
  colnames(this.dat)[colnames(this.dat)=="Age"] <- "agelab"
  this.dat$numage <- as.numeric(paste(this.dat$agelab))
  maxage <- max(this.dat$numage, na.rm=TRUE) + 1
  this.dat$numage[is.na(this.dat$numage)] <- maxage

  return(this.dat)
}

##############################################
#
# load.hmd.Nx
#
# function to load death counts (Dx)
# from the human mortality database
#
# arguments
#      country - the short-form of the country name
#       format - the format to use (eg "" or "5";
#                default is "" for single-year of age)
#      hmd.dir - the directory w/ the hmd data
# returns
#      a data frame with the mortality rates;
#   adds a numeric age variable and renames
#   Age to agelab
#
##############################################
load.hmd.Nx <- function( country, format="",
                         hmd.dir = "/Volumes/LaCie/data/hmd/hmd_countries") {

  filename <- paste(hmd.dir, "/", country, "/STATS/Population", format, ".txt", sep="")

  this.dat <- read.table( filename, skip=1, header=TRUE, na.strings="." )

  # TODO -- this may need to be made more general...
  #  (in particular, is oldest age the same for everyone?)
  colnames(this.dat)[colnames(this.dat)=="Age"] <- "agelab"
  this.dat$numage <- as.numeric(paste(this.dat$agelab))
  maxage <- max(this.dat$numage, na.rm=TRUE) + 1
  this.dat$numage[is.na(this.dat$numage)] <- maxage

  return(this.dat)
}

##############################################
#
# load.hmd.Exp
#
# function to load Exposures
# from the human mortality database
#
# arguments
#      country - the short-form of the country name
#       format - the format to use (eg "1x1")
#      hmd.dir - the directory w/ the hmd data
# returns
#      a data frame with the mortality rates;
#   adds a numeric age variable and renames
#   Age to agelab
#
##############################################
load.hmd.Exp <- function( country, format,
                          hmd.dir = "/Volumes/LaCie/data/hmd/hmd_countries") {

  filename <- paste(hmd.dir, "/", country, "/STATS/Exposures_", format, ".txt", sep="")

  this.dat <- read.table( filename, skip=1, header=TRUE, na.strings="." )

  # TODO -- this may need to be made more general...
  #  (in particular, is oldest age the same for everyone?)
  colnames(this.dat)[colnames(this.dat)=="Age"] <- "agelab"
  this.dat$numage <- as.numeric(paste(this.dat$agelab))
  maxage <- max(this.dat$numage, na.rm=TRUE) + 1
  this.dat$numage[is.na(this.dat$numage)] <- maxage

  return(this.dat)
}

##############################################
#
# load.hmd.lt
#
# function to load life tables
# from the human mortality database
#
# arguments
#      country - the short-form of the country name
#          sex - "male", "female", or "total"
#       cohort - if TRUE, cohort data; otherwise, period
#       format - the format to use (eg "1x1")
#      as.list - return a list with each entry a year?
#                (defaults to FALSE)
#                (NOT YET IMPLEMENTED)
#      hmd.dir - the directory w/ the hmd data
# returns
#      a data frame with the mortality rates;
#   adds a numeric age variable and renames
#   Age to agelab
#
##############################################
load.hmd.lt <- function( country, sex, format, cohort=FALSE,
                         as.list=FALSE,
                         hmd.dir = "/Volumes/LaCie/data/hmd/hmd_countries") {

  chrt <- "per"
  if (cohort) {
    chrt <- "coh"
  }

  sexchar <- ""
  if (sex == "male") {
    sexchar <- "m"
  } else if (sex == "female") {
    sexchar <- "f"
  } else {
    sexchar <- "b"
  }

  filename <- paste(hmd.dir, "/", country, "/STATS/", sexchar, "lt", chrt, "_", format, ".txt", sep="")

  this.dat <- read.table( filename, skip=1, header=TRUE, na.strings="." )

  # TODO -- this may need to be made more general...
  #  (in particular, is oldest age the same for everyone?)
  colnames(this.dat)[colnames(this.dat)=="Age"] <- "agelab"
  this.dat$numage <- as.numeric(paste(this.dat$agelab))
  maxage <- max(this.dat$numage, na.rm=TRUE) + 1
  this.dat$numage[is.na(this.dat$numage)] <- maxage



  return(this.dat)
}



