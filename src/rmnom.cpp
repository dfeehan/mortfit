/************************************************
 * rmnom.cpp
 *
 * simple fn to take draws from
 * a multinomial distribution, given
 * the sample size and the probs of
 * each category
 *
 * see "Writing a package that uses Rcpp"
 * by Edelbuettel and Francois, Sept 29 2011
 *
 * dennis, dec 2011
 ************************************************/

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rmnom(SEXP n, SEXP p)
{

  using namespace Rcpp;

  NumericVector xn(n);
  NumericMatrix xp(p);

  int len_n = xn.size();
  int num_cat = xp.ncol();

  int i, j, r;
  double thisp = 0;
  double cumsum = 0;

  // note that this is filled with 0s by default
  NumericMatrix result(xp.nrow(), xp.ncol());

  // TODO -- should check that number of rows in xp
  // is equal to the length of xn...

  /* seed the random number generator with the time
   * this may need fixing if used in a cluster
   */
  srand(time(0));

  /*
   * loop through each of the draws to make
   */
  for(j=0; j < len_n; j++) {

    /*
     * draw a category for each unit in the sample
     */
    for(i=0; i < xn[j]; i++) {

      /* get a random draw from the uniform distribution
       * between 0 and 1 for this unit
       */
      thisp = (double)rand()/ (double)RAND_MAX;

      /* cumulative sum is the lower bound to check
       * against; start it at zero
       */
      cumsum = 0;

      for( r=0; r < num_cat; r++) {

        if( thisp > cumsum & thisp <= cumsum + xp(j,r)) {
	  result(j,r)++;
	  cumsum = cumsum + xp(j,r);
        } else {
	  cumsum = cumsum + xp(j,r);
	}
      }

    }
  }

  return(result);

}


