\name{grid.search}
\alias{grid.search}
\title{grid.search}
\usage{
  grid.search(grid.vals, func, return.grid = FALSE,
    cpp = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{grid.vals}{either a list whose entries are vectors
  w/ parameter ranges, or a table of parameter values to
  evaluate}

  \item{func}{the function to maximize}

  \item{return.grid}{if TRUE, return the grid used}

  \item{cpp}{if TRUE, use the c++ code; otherwise, use R's
  apply()}

  \item{...}{other args, to be passed along to func}
}
\description{
  find the max (or min) of a function by brute-force search
  through a grid TODO - handle multiple maxima
}

