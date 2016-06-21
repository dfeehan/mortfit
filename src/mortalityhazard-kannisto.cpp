/************************************************
 * mortalityhazard-kannisto.cpp
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
NumericVector mortalityhazard_kannisto_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = exp(theta[1]);
  //double beta = theta[1];

  double num = 0.0;
  double denom = 0.0;

  for(int i=0; i < len; i++) {
    num = alpha*exp(beta*z[i]);
    denom = 1 + alpha*exp(beta*z[i]);
    res[i] = num/denom;
    if (res[i] < 0) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

}

// [[Rcpp::export]]
NumericVector mortalityhazard_to_prob_kannisto_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = exp(theta[1]);
  //double beta = theta[1];

  double temp = 0.0;

  for(int i=0; i < len; i++) {
    temp = (alpha * exp(beta * z[i]) + 1.0) / (alpha*exp(beta*(z[i]+1.0)) + 1.0);
    res[i] = 1.0 - pow(temp, 1.0/beta);
  }

  return(res);

}


