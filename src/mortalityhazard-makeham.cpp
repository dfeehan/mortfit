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

  double temp = 0.0;

  for(int i=0; i < len; i++) {
    temp = (alpha/beta)*(exp(beta*(z[i]+1))-exp(beta*z[i])) + gamma;
    res[i] = 1 - exp(-1*temp);
  }

  return(res);

}


