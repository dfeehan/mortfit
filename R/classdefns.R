#######################################################
##
## classdefns.R
##
## this file has most of the S4 class definitions
## used in the mortfit library
##

## see chambers (1998), pg 351

setClassUnion("listOrNULL",
              c("list", "NULL"))

setClassUnion("numericOrList",
              c("numeric", "list"))

setClassUnion("functionOrNULL",
              c("function", "NULL"))

#######################################################
##'
##' mortalityHazard
##'
##' representation of a mortality hazard function
##' \describe{
##'     \item{name}{the hazard's name}
##'     \item{num.param}{the number of parameters the hazard uses}
##'     \item{theta.default}{default values for the parameters}
##'     \item{theta.start.fn}{if not NULL, a function which takes a mortalityData object and returns starting values to be used in a numerical optimization routine}
##'     \item{theta.range}{if not NULL, a list with one entry for each parameter. Each parameter's entry has a range of plausible values for the parameter, which can be used if an optimization routine wants to choose a random starting value}
##'     \item{optim.default}{default settings to use with R's \code{\link{optim}} function}
##'     \item{haz.fn}{a function which takes two arguments, a vector of parameter values and a vector of ages, and returns the value of the hazard function with the given parameter values at each age.}
##'     \item{haz.to.prob.fn}{if not NULL a function which takes two arguments, a vector of parameter values and a vector of ages, and returns the value of \eqn{exp(-\int_{z}^{z+1} -\mu(x)dx)}, ie the result of converting the hazard to conditional probabilities of dying between ages z and z+1 for each age z.}##' 
##' }
##' @name mortalityHazard-class
##' @rdname mortalityHazard-class
##' @exportClass mortalityHazard
setClass("mortalityHazard",
         representation( name="character",
                         num.param="integer",
                         theta.default="numeric",
                         theta.start.fn="functionOrNULL",
                         theta.range="listOrNULL",
                         optim.default="list",
                         haz.fn="function",
                         haz.to.prob.fn="functionOrNULL"),
         ## TODO -- add prototype?
         validity=function(object) {
             ## TODO -- add validity checks
             return(TRUE)
         }
         )

setClassUnion("mortalityHazardOrNULL",
              c("mortalityHazard", "NULL"))

#######################################################
##'
##' mortalityData
##'
##' representation of a dataset that contains
##'
##' \describe{
##'     \item{name}{the dataset's name}
##'     \item{age.offset}{the amount to add to the age column of \code{data} in order to obtain the actual ages. For example, the ages in \code{data} might range from 1 to 20 and represent ages 80 to 99; in that case, \code{age.offset} would be 79}
##'     \item{age.interval.width}{the width, in years, of each age interval. for now, this is always assumed to be 1. future versions may accommodate different age intervals.}
##'     \item{data}{a data.frame with at least three columns:
##'           \describe{
##'               \item{age}{the age (last birthday)}
##'               \item{Nx}{the amount of exposure at the given age}
##'               \item{Dx}{the number of deaths at the given age}
##'           }
##'     }
##'     \item{tags}{a list whose key/value entries can be used to store metadata; useful for grabbing subsets of a big database. For example, this might have \code{list(country=\"belgi\", sex=\"m\")}}
##' }
##' @name mortalityData-class
##' @rdname mortalityData-class
##' @exportClass mortalityData
setClass("mortalityData",
         representation(name="character",
                        age.offset="integer",
                        age.interval.width="numeric",
                        data="data.frame",
                        tags="list"),
         ## TODO -- add prototype?
         validity=function(object) {
             ## TODO -- add validity checks
             return(TRUE)
         }
         )

#######################################################
##'
##' mortalityDataFolded
##'
##' representation of a dataset that contains
##' data on deaths and populations and that has
##' been split into K folds (for cross validation)
##' extends the \code{\link{mortalityData-class}} class
##'
##' \describe{
##'     \item{num.folds}{the number of folds the dataset has been divided into}
##'     \item{folds}{a list whose entries contain the \code{num.folds} folds}
##'     \item{fold.complements}{a list whose entries contain, for each fold, the data not in that fold}
##' }
##' @name mortalityDataFolded-class
##' @rdname mortalityDataFolded-class
##' @exportClass mortalityDataFolded
setClass("mortalityDataFolded",
         representation(num.folds="integer",
                        folds="list",
                        fold.complements="list"),
         contains="mortalityData",
         ## TODO -- add prototype?
         validity=function(object) {
             ## TODO -- add validity checks
             ## FOR EXAMPLE, should check that totals
             ##   add up across folds...
             return(TRUE)
         }
         )

#######################################################
##'
##' fitMethod
##'
##' representation of a strategy for fitting a
##' function to data. for example, one fitMethod included
##' with the \code{mortfit} package is \code{optimFit}, which
##' uses the BFGS algorithm in R's \link{stats:optim} package
##'
##' \describe{
##'     \item{name}{the name of the fit method}
##'     \item{fit}{a function which actually fits a model to data. it takes as its arguments
##'        \describe{
##'          \item{model.obj}{the mortalityModel object to fit}
##'          \item{data}{the mortalityData to fit the model to}
##'          \item{...}{other arguments, specific to the fit method}
##'        }}
##' }
##' @name fitMethod-class
##' @rdname fitMethod-class
##' @exportClass fitMethod
setClass("fitMethod",
         representation(name="character",
                        fit="function"
                        ),
         ## TODO -- add prototype?
         validity=function(object) {
             ## TODO -- add validity checks
             return(TRUE)
         }
         )

#######################################################
##'
##' fitSummary
##'
##' Summary of a the fit of a mortality model to
##' a mortality dataset. Note that many of these criteria
##' only make sense when compared to other models. 
##'
##' \describe{
##'     \item{AIC}{Akaike's Information Criterion}
##'     \item{BIC}{Bayesian Information Criterion}
##'     \item{SSE.Dx}{The sum of squared errors in the predicted number of deaths by age}
##'     \item{chisq}{The chi-square value}
##'     \item{chisq.p}{The p-value associated with the chi-square value}
##'     \item{chisq.byage}{Chi-square value for each age's prediction}
##'     \item{chisq.byage.p}{The p-values associated with the age-specific chi-square values}
##'     \item{deviance}{The deviance}
##' }
##' @name fitSummary-class
##' @rdname fitSummary-class
##' @exportClass fitSummary
setClass("fitSummary",
         representation(AIC="numeric",
                        BIC="numeric",
                        SSE.Dx="numeric",
                        chisq="numeric",
                        chisq.p="numeric",
                        chisq.byage="numeric",
                        chisq.byage.p="numeric",
                        deviance="numeric"
                        ),
         ## TODO -- add prototype?
         validity=function(object) {
             ## TODO -- add validity checks
             return(TRUE)
         }
         )

setClassUnion("fitSummaryOrNULL",
              c("fitSummary", "NULL"))

#######################################################
##'
##' fitSummaryCV
##'
##' summary of a the fit of a mortality model to
##' a mortality dataset, including cross-validation;
##' extends \code{\link{fitSummary-class}}
##'
##' \describe{
##'     \item{cv.err.Dx}{the sum of squared errors in using cross validation to predict out of sample}
##'     \item{cv.rmse.Dx}{the root mean squared error derived from \code{cv.err.Dx}}
##'     \item{fold.fits}{a list whose entries have the results of fitting the model for each fold (see \code{mort.fit.cv}}
##' }
##' @name fitSummaryCV-class
##' @rdname fitSummaryCV-class
##' @exportClass fitSummaryCV
setClass("fitSummaryCV",
         contains="fitSummary",
         representation(cv.err.Dx="list",
                        cv.rmse.Dx="numeric",
                        fold.fits="list"),
         ## TODO -- add prototype?
         validity=function(object) {
             ## TODO -- add validity checks
             return(TRUE)
         }
         )

#######################################################
##'
##' mortalityModel
##'
##' this object represents a likelihood function for
##' use with mortality data. generally, it will be
##' the result of combining a likelihood with a hazard;
##' for example, it could be the binomial model with the gompertz hazard
##' embedded in it. many of the slots of this object correspond
##' to those found in \code{\link{mortalityHazard-class}}
##'
##' \describe{
##'     \item{name}{the model's name}
##'     \item{loglik.fn}{the log-likelihood function, which takes as its arguments
##'        \describe{
##'           \item{theta}{a vector of parameter values}
##'           \item{Dx}{vector of deaths by age}
##'           \item{Nx}{vector of exposure by age}
##'           \item{ages}{vector of ages that corresponds to the entries of \code{Nx} and \code{Dx}}
##'       } }
##'    \item{num.param}{the number of parameters in the model}
##'    \item{theta.default}{default values for the parameters}
##'    \item{theta.range}{if not NULL, a range of plausible values that the parameters can assume}
##'    \item{optim.default}{a list of default parameters to pass into R's \code{optim} function}
##'    \item{predict.fn}{a function which, given values for the parameters and a dataset of exposures, produces predicted numbers of deaths by age. it takes as its arguments:
##'         \describe{
##'           \item{theta}{a vector of parameter values}
##'           \item{Nx}{vector of exposure by age}
##'           \item{age}{vector of ages that the \code{Nx} correspond to}
##'           \item{age.offset}{the offset used for \code{age}. (See \code{\link{mortalityData-class}}.)}
##'         } }
##'    \item{theta.start.fn}{if not NULL, a function which takes a mortalityData object and returns very preliminary estimates of the parameters which can then be used in an optimization routine}
##'    \item{simulate.fn}{if not NULL, a function which takes as its arguments:
##'         \describe{
##'           \item{theta}{a vector of parameter values}
##'           \item{Nx}{vector of exposure by age}
##'           \item{age}{vector of ages that the \code{Nx} correspond to}
##'           \item{seed}{if not NULL, the value to seed the pseudorandom number generator with}
##'           \item{means}{TODO}
##'         } }
##'     \item{eval.fn}{if not NULL, a function which produces an evaluation of the model's fit to a given dataset. it takes as its arguments:
##'         \describe{
##'           \item{fit.obj}{the mortalityFit object that resulted from fitting the model to a dataset}
##'           \item{observed}{the mortalityData object with the observed data}
##'           \item{fittedValues}{the mortalityPrediction object with fitted values, usually resulting from a call to \code{predict.fn}}
##'         } }
##'     \item{hazard}{if not NULL, the mortalityHazard object used in making this mortalityModel}
##' }
##' @name mortalityModel-class
##' @rdname mortalityModel-class
##' @exportClass mortalityModel
setClass("mortalityModel",
         representation(name="character",
                        loglik.fn="function",
                        num.param="integer",
                        theta.default="numeric",
                        theta.range="listOrNULL",
                        optim.default="list",
                        predict.fn="function",
                        theta.start.fn="functionOrNULL",
                        simulate.fn="functionOrNULL",
                        eval.fn="functionOrNULL",
                        hazard="mortalityHazardOrNULL"),
         ## TODO -- add prototype?
         validity=function(object) {
             ## TODO -- add validity checks
             return(TRUE)
         }
         )

#######################################################
##'
##' mortalityPrediction
##'
##' this object contains a set of fitted values,
##' generated from a model's prediction
##'
##' \describe{
##'     \item{name}{the name of this mortalityPrediction object}
##'     \item{fitted.Dx}{vector with the estimated number of deaths by age}
##'     \item{fitted.qx}{vector with the estimated probabilities of death by age}
##'     \item{theta}{vector of parameter values used for these predictions}
##'     \item{age}{vector of ages corresponding to the entries of \code{fitted.Dx} and \code{fitted.Dx}}
##'     \item{age.offset}{the offset used for ages; for example, if \code{age} has values 1 to 20 and \code{age.offset} is 79, then the real ages being studied range from 80 to 99.}
##'     \item{Nx}{the amount of exposure, for each age, used in producing these predictions}
##' }
##' @name mortalityPrediction-class
##' @rdname mortalityPrediction-class
##' @exportClass mortalityPrediction
setClass("mortalityPrediction",
         representation(name="character",
                        fitted.Dx="numeric",
                        fitted.qx="numeric",
                        theta="numeric",
                        age="numeric",
                        age.offset="numeric",
                        Nx="numeric"),
         ## TODO -- add prototype?
         validity=function(object) {
             ## TODO -- add validity checks
             return(TRUE)
         }
         )

#######################################################
##'
##' mortalityFit
##'
##' representation of a fit of a hazard function
##' (given as a mortalityHazard object)
##'  to mortality data (given as a
##'  mortalityData object)
##'
##' \describe{
##'     \item{name}{the name of this mortality fit}
##'     \item{model}{the mortalityModel object used to produce this fit}
##'     \item{data}{the mortalityData object used to produce this fit}
##'     \item{method}{the fitMethod used to produce this fit}
##'     \item{theta.init}{vector with the initial estimates of the parameters used to produce this fit}
##'     \item{theta.hat}{the final estimated values of the parameters}
##'     \item{log.likelihood}{the log-likelihood of this model fit to this data.}
##'     \item{fitted.values}{a mortalityPrediction object containing the actual predicted values for this model and dataset}
##'     \item{fit.summary}{if not NULL, a fitSummary object containing a summary of various measures of how well this model fits this dataset}
##' }
##' @name mortalityFit-class
##' @rdname mortalityFit-class
##' @exportClass mortalityFit
setClass("mortalityFit",
         representation(name="character",
                        model="mortalityModel",
                        data="mortalityData",
                        method="fitMethod",
                        theta.init="numericOrList",
                        theta.start.fn="functionOrNULL",
                        theta.hat="numeric",
                        log.likelihood="numeric",
                        fitted.values="mortalityPrediction",
                        fit.summary="fitSummaryOrNULL"),
         ## TODO -- add prototype?
         validity=function(object) {
             ## TODO -- add validity checks
             return(TRUE)
         }
         )

#######################################################
##'
##' mortalityFits
##'
##' representation of a list of several
##' fits of mortalityHazard objects to
##' the same data
##'
##' \describe{
##'     \item{name}{the name of this set of fits}
##'     \item{data}{the mortalityData object that all of these models were fit to}
##'     \item{fits}{a list whose entries contain the fits of all of the models}
##' }
##' @name mortalityFits-class
##' @rdname mortalityFits-class
##' @exportClass mortalityFits
setClass("mortalityFits",
         representation(name="character",
                        data="mortalityData",
                        fits="list"),
         ## TODO -- add prototype?
         validity=function(object) {
             ## TODO -- add validity checks
             return(TRUE)
         }
         )

########################################################
##'
##' as.mortalityData.mortalityDataFolded
##'
##' convert a mortalityDataFolded object into
##' a mortalityData one by losing the info on the folds
##'
##' @param obj the mortalityDataFolded object to convert
##'            to a mortalityData one
##' @param ... other params (not yet used)
##' @return a mortalityData object
##' @export
as.mortalityData.mortalityDataFolded <- function(obj,...) {

  return(new("mortalityData",
             name=paste("as.mortalityData:", obj@name),
             age.offset=obj@age.offset,
             age.interval.width=obj@age.interval.width,
             data=obj@data,
             tags=obj@tags
             ))  
}

########################################################
##'
##' as.data.frame.mortalityFit
##'
##' convert a mortalityFit object to a data.frame
##' this function is intended to be called from
##' as.data.frame.mortalityFits (ie, from a list of
##' data frames)
##'
##' @param obj the mortalityFit object to convert to a data.frame
##' @return a data.frame containing the info from the
##'         mortalityFit object
as.data.frame.mortalityFit <- function(obj) {

  if (! is(obj@fit.summary, "fitSummaryCV")) {
    this.cv.rmse.Dx <- NA
  } else {
    this.cv.rmse.Dx <- obj@fit.summary@cv.rmse.Dx
  }
  
  return(data.frame(model=obj@model@name,
                    data=obj@data@name,
                    method=obj@method@name,
                    num.param=obj@model@num.param,
                    log.likelihood=obj@log.likelihood,
                    AIC=obj@fit.summary@AIC,
                    BIC=obj@fit.summary@BIC,
                    SSE.Dx=obj@fit.summary@SSE.Dx,
                    chisq=obj@fit.summary@chisq,
                    chisq.p=obj@fit.summary@chisq.p,
                    deviance=obj@fit.summary@deviance,
                    CV.rmse.Dx=this.cv.rmse.Dx))
  
}

########################################################
##'
##' as.data.frame.mortalityFits
##'
##' convert a mortalityFits object, which is basically
##' a list of mortalityFit objects, to a data.frame
##'
##' @param obj the list of mortalityFit objects
##' @param order.by if not NULL, the name of a column to order
##'                 the results by
##' @param ... other arguments to pass to order
##' @return a data.frame containing the info from the
##'         mortalityFit objects listed in obj
##' @export
as.data.frame.mortalityFits <- function(obj,
                                        order.by=NULL,
                                        ...) {

  res <- ldply(obj@fits,
               as.data.frame.mortalityFit)

  if(! is.null(order.by)) {
    res <- res[order(res[,order.by],...),]
  }
  
  return(res)
  
}


########################################################
##'
##' compare.fits.mortalityFits
##'
##' given a list of mortalityFit objects, contained in a
##' mortalityFits object, produce comparisons
##' between them along a number of dimensions.
##' NB/TODO: for now, assuming that all were fit
##' to the same data. eventually, it would be nice to explicitly
##' check this
##'
##' @param mf a data frame with the fit results, as produced
##'               by as.data.frame.mortalityFits
##' @param ... other params, not currently  used
##' @return a list with entries df, which is a dataframe with
##'        the fit results, and summary
##' @export
compare.fits.mortalityFits <- function(mf,...) {

  df <- as.data.frame(mf, order.by="AIC")

  df$delta.AIC <- df$AIC - df$AIC[1]
  
  df$rank.delta.AIC <- rank(df$delta.AIC)
  df$rank.log.likelihood <- rank(-df$log.likelihood)
  df$rank.SSE.Dx <- rank(df$SSE.Dx)
  df$rank.CV.rmse.Dx <- rank(df$CV.rmse.Dx)
  df$rank.chisq <- rank(df$chisq)
  df$rank.deviance <- rank(df$deviance)

  ## TODO -- LEFT OFF HERE: developing this function...
  
  return(list(df=df,
              summary="this will be a summary"))
}

##setGeneric("compare.fits",
##           useAsDefault=compare.fits.mortalityFits)

##setMethod("compare.fits",
##          signature(mf="mortalityFits"),
##          function(mf,...) {
##            compare.fits.mortalityFits(x)
##          })

########################################################
##'
##' [[ function for mortalityFits object
setMethod("[[",
          signature(x = "mortalityFits"),
          function(x, i, j, ...) {
            return(x@fits[[i]])
          })

#########################################################
##'
##' function for as.data.frame applied to a mortalityFits
##' object
setMethod("as.data.frame",
          "mortalityFits",
          function(x,...) {
            as.data.frame.mortalityFits(x,...)
          })

##setGeneric("as.mortalityData",
##           useAsDefault=as.mortalityData.mortalityDataFolded)
##setGeneric("as.mortalityData")

