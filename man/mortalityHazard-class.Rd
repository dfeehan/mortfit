\name{mortalityHazard-class}
\alias{mortalityHazard-class}
\title{mortalityHazard}
\description{
  representation of a mortality hazard function \describe{
  \item{name}{the hazard's name} \item{num.param}{the
  number of parameters the hazard uses}
  \item{theta.default}{default values for the parameters}
  \item{theta.start.fn}{if not NULL, a function which takes
  a mortalityData object and returns starting values to be
  used in a numerical optimization routine}
  \item{theta.range}{if not NULL, a list with one entry for
  each parameter. Each parameter's entry has a range of
  plausible values for the parameter, which can be used if
  an optimization routine wants to choose a random starting
  value} \item{optim.default}{default settings to use with
  R's \code{\link{optim}} function} \item{haz.fn}{a
  function which takes two arguments, a vector of parameter
  values and a vector of ages, and returns the value of the
  hazard function with the given parameter values at each
  age.} \item{haz.to.prob.fn}{if not NULL a function which
  takes two arguments, a vector of parameter values and a
  vector of ages, and returns the value of
  \eqn{exp(-\int_{z}^{z+1} -\mu(x)dx)}, ie the result of
  converting the hazard to conditional probabilities of
  dying between ages z and z+1 for each age z.}##' }
}

