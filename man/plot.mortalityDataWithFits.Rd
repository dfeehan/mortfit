\name{plot.mortalityDataWithFits}
\alias{plot.mortalityDataWithFits}
\title{method to plot a set of mortalityFit objects and the
data they were fitted to (results in a grid of plots)}
\usage{
  plot.mortalityDataWithFits(x, y = NULL, Dx = TRUE,
    ncol = NULL, pdffile = NULL, ...)
}
\arguments{
  \item{x}{mortalityFit object}

  \item{y}{not used}

  \item{Dx}{if TRUE, plot in terms of deaths; otherwise,
  use central death rates}

  \item{ncol}{if not NULL, the number of columns in the
  array of plots produced}

  \item{pdffile}{the name of a .pdf file to save to (if
  NULL, plot is displayed to the screen)}

  \item{...}{other args, which are passed on to the pdf
  device}
}
\description{
  method to plot a set of mortalityFit objects and the data
  they were fitted to (results in a grid of plots)
}

