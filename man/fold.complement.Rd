\name{fold.complement}
\alias{fold.complement}
\title{fold.complement}
\usage{
  fold.complement(data, fold)
}
\arguments{
  \item{data}{a data.frame with the entire dataset which
  was split into folds; must have columns 'Nx' and 'Dx'}

  \item{fold}{a data.frame with the data for the fold we
  wish to complement; must have columns 'Nx' and 'Dx'}
}
\value{
  a data.frame with the complement of fold in data
}
\description{
  return all of the obs in data that aren't in the given
  fold
}

