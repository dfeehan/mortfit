\name{as.data.frame.mortalityFits}
\alias{as.data.frame.mortalityFits}
\title{as.data.frame.mortalityFits}
\usage{
  as.data.frame.mortalityFits(obj, order.by = NULL, ...)
}
\arguments{
  \item{obj}{the list of mortalityFit objects}

  \item{order.by}{if not NULL, the name of a column to
  order the results by}

  \item{...}{other arguments to pass to order}
}
\value{
  a data.frame containing the info from the mortalityFit
  objects listed in obj
}
\description{
  convert a mortalityFits object, which is basically a list
  of mortalityFit objects, to a data.frame
}

