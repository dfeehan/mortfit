% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/library.R
\name{plot.ll.prof}
\alias{plot.ll.prof}
\title{plot likelihood profile for a pair of parameters}
\usage{
plot.ll.prof(fit, theta.idx = c(1:2), delta = 0.1, res = 200,
  heatmap = FALSE)
}
\arguments{
\item{fit}{a mortalityFit object}

\item{theta.idx}{vector with the indexes of the two parameters to plot against
each other; for example, in a four-parameter model,
to examine the profile of the likelihood as
parameters one and three vary, pass in theta.idx=c(1,3)}

\item{delta}{window around the fitted max to profile. if delta
is 0.1, then plot values for (0.9*max, 1.1*max), where
max is the fitted maximum parameter value}

\item{res}{the resolution, or size of the grid on which the
likelihood is evaluated. large grids may take a long time}

\item{heatmap}{if TRUE, return a heatmap for each pair of parameters;
otherwise, return a contour plot}
}
\value{
a ggplot2 plots of the profile
}
\description{
given a fitted mortalityModel, plot the profile
of the log-likelihood function for a given pair
of parameters. this can take a while if the likelihood
takes a while to evaluate and/or the resolution is
chosen to be large
}

