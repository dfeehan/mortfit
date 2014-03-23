\name{optim.fit}
\alias{optim.fit}
\title{optim.fit}
\usage{
  optim.fit(model.obj, data, theta.init = NULL,
    random.start = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{model.obj}{a mortalityModel object}

  \item{data}{a mortalityData object}

  \item{theta.init}{if not NULL, starting values to use for
  theta}

  \item{random.start}{if TRUE, and the model.obj has a
  non-NULL theta.range, then draw a random starting value
  uniformly from the values specified by theta.range}

  \item{verbose}{if TRUE, print extra info}
}
\value{
  a mortalityFitOptim object
}
\description{
  use optim to fit a mortality model
}

