##########################################################
##
## rmnom.R
##
## this file wraps the c++ code for the rmnom function

##' rmnom
##'
##' take a draw from a multinomial distribution
##' TODO -- example
##' 
##' @param n a vector with the number of trials for each draw
##' @param p a matrix with one row for each draw (so number of rows should be the same as the length of the vector n). Across each row, p has the probability of each category (so the rows all sum to 1).
##' @export
rmnom <- function(n, p) {
  res <- .Call("rmnom", n, p)
  return(res)
}
