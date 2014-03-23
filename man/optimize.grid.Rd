\name{optimize.grid}
\alias{optimize.grid}
\title{optimize.grid}
\usage{
  optimize.grid(grid.vals, M = 3, delta = 0.1, func,
    epsilon = 1e-06, verbose = FALSE, cpp = FALSE, ...)
}
\arguments{
  \item{grid.vals}{a list with an entry for each param,
  which has a vector of all of the values the param should
  take on}

  \item{M}{the number of times to run grid search (refining
  results each time...)}

  \item{delta}{the size of the window to construct grid
  over. if delta=0.1, make a grid from (0.9*alpha,
  1.1*alpha)}

  \item{func}{the function to optimize}

  \item{epsilon}{the minimum difference to consider
  significant; if the difference between the current
  optimum and the previous one is less than epsilon (in
  absolute value), then stop}

  \item{verbose}{if TRUE, print out info about
  optimization}

  \item{cpp}{if TRUE, use the C++ version... (default is
  FALSE)}

  \item{...}{other arguments, which get passed along to
  \link{grid.search}}
}
\value{
  TODO
}
\description{
  maximize a function using grid search TODO -- not
  finished documenting this
}

