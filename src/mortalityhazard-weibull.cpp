/************************************************
 * mortalityhazard-weibull.cpp
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
NumericVector mortalityhazard_weibull_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = theta[0];
  double beta = theta[1];

  for(int i=0; i < len; i++) {

    res[i] = alpha*pow(z[i], beta);

    if (res[i] < 0) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

}

// [[Rcpp::export]]
NumericVector mortalityhazard_to_prob_weibull_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = theta[0];
  double beta = theta[1];

  double temp = 0.0;

  for(int i=0; i < len; i++) {
    temp = (alpha*(pow((z[i]+1),(1+beta))-pow(z[i],(1+beta))))/(1+beta);
    res[i] = 1 - exp(-1*temp);
  }

  return(res);

}


