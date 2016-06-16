/************************************************
 * mortalityhazard-makeham.cpp
 *
 *
 *
 * see "Writing a package that uses Rcpp"
 * by Edelbuettel and Francois, Sept 29 2011
 *
 * dennis, dec 2011
 ************************************************/

#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mortalityhazard_makeham_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  int parlen = theta.size();

  if (parlen != 3) {

  }

  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = theta[1];
  double gamma = exp(theta[2]);

  for(int i=0; i < len; i++) {

    res[i] = gamma + alpha*exp(beta*z[i]);

    if (res[i] < 0) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

}

// [[Rcpp::export]]
NumericVector mortalityhazard_to_prob_makeham_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = theta[1];
  double gamma = exp(theta[2]);

  double temp1 = 0.0, temp2=0.0;

  /*
   * we want
   *     1 - exp(-gamma -(alpha/beta) * [exp(beta*(z+1)) - exp(beta*z)])
   * but we'll perform this as
   *     1 - [exp(gamma) * exp(-(alpha/beta)*exp(beta*(z+1))) * exp((alpha/beta)*exp(beta*z))]
   * to try to minimize possible float cancellation problems with the difference of the exps
   */

  for(int i=0; i < len; i++) {
    temp1 = -1.0 * (alpha/beta)*exp(beta*(z[i]+1));
    temp2 = (alpha/beta)*(exp(beta*z[i]));
    res[i] = 1 - (exp(-1.0*gamma) * exp(temp1) * exp(temp2));
  }

  return(res);

}


