\name{lre}
\alias{lre}
\title{lre}
\usage{
  lre(theta.true, theta.hat, j = 1)
}
\arguments{
  \item{theta.true}{the true parameter vector}

  \item{theta.hat}{the estimated parameter vector}

  \item{j}{scaling factor}
}
\value{
  a vector of the same length as theta.true
}
\description{
  compute the log10 relative error in an estimated versus a
  true parameter vector. this is useful for unit tests to
  be sure optimization routines are well-behaved see Altman
  et al, "Numerical Issues in Statistical Computing...", pg
  54 TODO -- come back to this; j is not defined in their
  expression (3.4)!
}

