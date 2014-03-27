#######################################################
##
## plotmethods.R
##
## this file has various plotting functions used
## in the mortfit library
## it also sets them up as generic methods with the
## appropriate signatures.
##
## @include mortfit-help.R
## @include classdefns.R

##########################################################
##' method to plot mortality data
##'
##' @param x the mortalityData object to plot
##' @param y (not used)
##' @param ... TODO
##' @export
plot.mortalityData <- function(x,y=NULL,...) {

    x@data$plotage <- x@data$age + x@age.offset

    toret <- ggplot(x@data) +
        geom_point(aes(x=plotage,y=(Dx/Nx),size=Nx),pch=1) +
        xlab("age") + ylab("Dx/Nx") + opts(title=x@name) + scale_area()
    return(toret)
}

##########################################################
##' method to plot mortality data
##'
##' @param x the mortalityData object to plot
##' @param y (not used)
##' @param ... TODO 
##' @export
plot.mortalityDataFolded <- function(x,y=NULL,...) {

    x@data$plotage <- x@data$age + x@age.offset

    toret <- ggplot(x@data) +
        geom_point(aes(x=plotage,y=(Dx/Nx),size=Nx),pch=1) +
        xlab("age") + ylab("Dx/Nx") + opts(title=x@name) + scale_area()
    return(toret)
}


#########################################################
##' Plot mortality hazard, represented by a
##' mortalityHazard object and possibly a
##' vector of parameters
##'
##' @param x the mortalityHazard object to plot
##' @param age the age range to plot over
##' @param age.offset the offset used (often 80)
##' @param theta the vector of parameters to use. if
##'    none is specified, use the mortalityHazard object's defaults
##' @return a ggplot2 object containing the plot
##' @export
plot.mortalityHazard <- function(x,
                                 age=1:19,
                                 age.offset=80,
                                 theta=NULL) {

  if (is.null(theta)) {
    theta <- x@theta.default
  }

  mu <- x@haz.fn(theta,
                 age)
  
  plotdat <- data.frame(age=(age+age.offset), mu=mu)
  
  toret <- ggplot(plotdat) +
           geom_line(aes(x=age,y=mu)) +
           opts(title=paste(x@name, " - (",
                            paste(theta, collapse=","),
                            ") - age offset: ", age.offset, sep=""))

  return(toret)
}


## TODO -- add plot method for mortalityDataFolded?


##############################################################
##' method to plot mortality hazard on top
##' of data
##' method to plot mortality hazard
##
##' @param  x mortalityData object
##' @param  y mortalityHazard object
##' @param theta the parameter values to use (defaults are chosen if NULL)
##' @export
plot.mortalityDataWithHazard <- function(x,y,
                                         ##age=80:99,
                                         theta=NULL) {

  ## if no paramater values specified, use
  ## this hazard's defaults
  if (is.null(theta)) {
    theta <- y@theta.default
  }

  age <- x@data$age
  age.offset <- x@age.offset
  
  mu <- y@haz.fn(theta,
                 age)

  pi <- haz.to.prob(y@haz.fn, theta, age)

  haz.dat <- data.frame(age=age + age.offset, mu=mu, pi=pi)
  
  plot.dat <- x@data
  ##mort.dat <- x@data  

  ##plot.dat <- merge(haz.dat,
  ##                  mort.dat,
  ##                  by="age",
  ##                  all=TRUE)

  ## plot central death rates halfway through each
  ## age interval...
  plot.dat$age <- plot.dat$age + age.offset + 0.5

  ## TODO -- add legend and improve color scheme. for now,
  ##   - black line is hazard
  ##   - blue circles are central death rates, Dx/(Nx-0.5*Dx)
  ##     (which sometimes don't estimate the hazard super well, eg if
  ##      hazard is very large)
  ##   - red line is true probabilities
  ##   - yellow triangles are estimated probabilities, via Dx/Nx
  ##             
  
  toret <- ggplot(plot.dat) +
           geom_line(aes(x=age,y=mu), data=haz.dat) +
           geom_point(aes(x=age,y=(Dx/(Nx-.5*Dx)),size=Nx),pch=1,color="blue") +
           ##geom_line(aes(x=(age),y=pi),color="red", data=haz.dat)+
           ##geom_point(aes(x=(age-0.5), y=(Dx/Nx)), color="yellow", pch=2)+
           labs(x="age", y="hazard /\ncentral death rate") +
           opts(title=paste(x@name, "\n", y@name, "\n",
                            "(", paste(round(theta,4),collapse=", "), ")", sep="")) +
           scale_area()

  return(toret)
    
}

##########################################################
##' method to plot a mortalityFit object and the data
##' it was fitted to
##'
##' @param  x mortalityFit object
##' @param  y (needed to match generic signature -- not used)
##' @param Dx if TRUE, plot deaths; otherwise, plot central death rates
##' @export
plot.mortalityDataWithFit <- function(x,
                                      y=NULL,
                                      Dx=FALSE) {

  if(! Dx) {

    dataplot <- plot(x@data,
                     x@model@hazard,
                     theta=x@theta.hat)

  } else {

    obsDx <- x@data@data$Dx
    obsNx <- x@data@data$Nx
    fitDx <- x@fitted.values@fitted.Dx
    ages <- x@fitted.values@age + x@fitted.values@age.offset

    dataplot <- ggplot(data=data.frame(obsDx=obsDx,fitDx=fitDx,ages=ages,obsNx=obsNx)) +
                geom_point(aes(x=ages,y=obsDx, size=obsNx),color="blue",pch=1) +
                geom_line(aes(x=ages, y=fitDx),color="red") +
                geom_point(aes(x=ages, y=fitDx),color="red",pch=3) +
                xlab("age") + ylab("number of deaths") +
                scale_area() +
                opts(title=x@name)
    
  }

  return(dataplot)
                     
}

##########################################################
##' method to plot a set of mortalityFit objects and the
##' data they were fitted to (results in a grid of plots)
##'
##'  @param x mortalityFit object
##'  @param y not used
##'  @param Dx if TRUE, plot in terms of deaths; otherwise, use central death rates
##'  @param ncol if not NULL, the number of columns in the array of plots produced
##'  @param pdffile the name of a .pdf file to save to (if NULL, plot is displayed to the screen)
##'  @param ... other args, which are passed on to the pdf device
##'  @export
plot.mortalityDataWithFits <- function(x,
                                       y=NULL,
                                       Dx=TRUE,
                                       ncol=NULL,
                                       pdffile=NULL,
                                       ...) {

  resplots <- llply(x@fits,
                    function(x) {
                      return(plot(x,Dx=Dx) + opts(legend.position="none"))
                    })
                    ##plot,
                    ##Dx=Dx)

  ## NB: see help page for grid.arrange in the gridExtra library
  ##legGrob <- ggplotGrob(resplots[[1]] + opts(keep="legend_box"))
  ##legend <- gTree(children=gList(legGrob), cl="legendGrob")
  ##widthDetails.legendGrob <- function(x) unit(2, "cm")

  ##resplots <- llply(resplots,
  ##                  function(x) { return(x + opts(legend.position="none")) })

  
  if (! is.null(pdffile)) {
    pdf(file=pdffile, ...)
  }
  ##uberplot <- do.call(grid.arrange,c(resplots,list(legend=legend,
  ##                                                 ncol=ncol)))

  uberplot <- do.call(grid.arrange, c(resplots, list(ncol=ncol)))
  
  if (! is.null(pdffile)) {
    dev.off()
  }
    
  return(uberplot)

}

############################################################
## TODO -- WRITE THIS: plot comparing mortality data and
##         predicted values
plot.mortalityDataWithPrediction <- function(x,
                                             y=NULL,
                                             theta=NULL) {
  ## TODO FINISH THIS...
  stop("not yet written.")
  
}

##setGeneric("plot",
##           signature=c("x","y"))
##setGeneric("plot",
##           function(x,y,...) {
##             standardGeneric("plot")
##             })
## for reasons I don't entirely understand, when moving
## to a package, this is necessary. perhaps it's because
## once namespaces are involved, i need to explicitly point
## out which function is the default plot? see
## http://r.789695.n4.nabble.com/Overloading-S4-methods-td3565588.html
## and also http://r.789695.n4.nabble.com/S4-plot-generic-documentation-td2543096.html

#####################################################
##' generic plot method for mortfit package
##'
##' @param x req'd as part of the generic defn of plot
##' @param y req'd as part of the generic defn of plot
##' @param ... used in some contexts
##' @seealso \code{\link{plot.mortalityDataWithFit},
##'                \link{plot.mortalityData},
##'                \link{plot.mortalityDataWithHazard},
##'                \link{plot.mortalityHazard},
##'                \link{plot.mortalityDataWithFit}}
##' @export
##' @docType methods
setGeneric("plot", 
           function(x, y, ...) {
             standardGeneric("plot")
           })

setMethod("plot",
          signature=c(x="mortalityFitOptim",
                      y="missing"),
          function(x,y,...) {
            plot.mortalityDataWithFit(x,y,...)
          })

setMethod("plot",
          signature=c(x="mortalityData",y="missing"),
          function(x,y,...) {
            plot.mortalityDataFolded(x,y,...)
          })

## FOR NOW, also use this method for folded data.
## eventually, we should have a custom method
setMethod("plot",
          signature=c(x="mortalityDataFolded",y="missing"),
          function(x,y,...) {
            plot.mortalityDataFolded(x,y,...)
          })


setMethod("plot",
          signature=c(x="mortalityData",
                      y="mortalityHazard"),
          function(x,y,...) {
            plot.mortalityDataWithHazard(x,y,...)
            })

setMethod("plot",
          signature=c(x="mortalityHazard",
                      y="missing"),
          function(x,y,...) {
            plot.mortalityHazard(x,...)
            })

setMethod("plot",
          signature=c(x="mortalityFit",
                      y="missing"),
          function(x,y,...) {
            plot.mortalityDataWithFit(x,...)
          })

setMethod("plot",
          signature=c(x="mortalityFits",
                      y="missing"),
          function(x,y,...) {
            plot.mortalityDataWithFits(x,...)
          })

