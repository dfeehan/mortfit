/************************************************
 * mortalityhazard-gompertz.cpp
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
NumericVector mortalityhazard_gompertz_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = theta[1];

  for(int i=0; i < len; i++) {
    res[i] = alpha*exp(beta*z[i]);
  }

  return(res);

}

// [[Rcpp::export]]
NumericVector mortalityhazard_to_prob_gompertz_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = theta[1];
  double temp = 0.0;

  for(int i=0; i < len; i++) {
    temp = (alpha/beta)*(exp(beta*(z[i]+1))-exp(beta*z[i]));
    res[i] = 1 - exp(-1*temp);
  }

  return(res);

}


