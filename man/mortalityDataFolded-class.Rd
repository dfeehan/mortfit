\name{mortalityDataFolded-class}
\alias{mortalityDataFolded-class}
\title{mortalityDataFolded}
\description{
  representation of a dataset that contains data on deaths
  and populations and that has been split into K folds (for
  cross validation) extends the
  \code{\link{mortalityData-class}} class
}
\details{
  \describe{ \item{num.folds}{the number of folds the
  dataset has been divided into} \item{folds}{a list whose
  entries contain the \code{num.folds} folds}
  \item{fold.complements}{a list whose entries contain, for
  each fold, the data not in that fold} }
}

