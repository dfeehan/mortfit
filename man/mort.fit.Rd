\name{mort.fit}
\alias{mort.fit}
\title{mort.fit}
\usage{
  mort.fit(model.obj, data, method, verbose = TRUE,
    ignore.folded = FALSE, ...)
}
\arguments{
  \item{model.obj}{a mortalityModel object with the
  mortality model to be fit}

  \item{data}{a mortalityData object with the data to fit
  to}

  \item{method}{a fitMethod object with the fit method to
  use}

  \item{verbose}{if TRUE, print out extra details while the
  fit proceeds}

  \item{ignore.folded}{if TRUE, then even if the data
  passed in are a mortalityDataFolded object, only run fits
  on the entire dataset (not each fold)}

  \item{...}{other parameters to pass to method}
}
\value{
  a mortalityFit object
}
\description{
  use a fitMethod to fit a mortalityModel to a set of
  mortalityData. TODO -- more description
}

