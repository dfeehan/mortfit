\name{compare.fits.mortalityFits}
\alias{compare.fits.mortalityFits}
\title{compare.fits.mortalityFits}
\usage{
  compare.fits.mortalityFits(mf, ...)
}
\arguments{
  \item{mf}{a data frame with the fit results, as produced
  by as.data.frame.mortalityFits}

  \item{...}{other params, not currently used}
}
\value{
  a list with entries df, which is a dataframe with the fit
  results, and summary
}
\description{
  given a list of mortalityFit objects, contained in a
  mortalityFits object, produce comparisons between them
  along a number of dimensions. NB/TODO: for now, assuming
  that all were fit to the same data. eventually, it would
  be nice to explicitly check this
}

