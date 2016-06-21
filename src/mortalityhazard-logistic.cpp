/************************************************
 * mortalityhazard-logistic.cpp
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
NumericVector mortalityhazard_logistic_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = exp(theta[1]);
  double gamma = exp(theta[2]);
  double delta = exp(theta[3]);

  double k = 0.0;
  double num = 0.0;
  double denom = 0.0;

  for(int i=0; i < len; i++) {
    k = exp(beta * z[i]);
    num = alpha * k;
    denom = 1 + (delta*k);

    /*Rprintf("denom is %f\n", denom);*/

    res[i] = (num/denom) + gamma;
  }

  bool anyNA = is_true( any( is_na(res) ) );

  if (anyNA) {
    for(int i=0; i<len; i++) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

}

// [[Rcpp::export]]
NumericVector mortalityhazard_to_prob_logistic_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = exp(theta[1]);
  double gamma = exp(theta[2]);
  double delta = exp(theta[3]);

  double temp = 0.0;


  for(int i=0; i < len; i++) {
    temp = (delta*exp(beta*z[i]) + 1.0) / (delta*exp(beta*(z[i]+1)) + 1.0);
    res[i] = 1.0 - exp(-1.0*gamma) * pow(temp, alpha/(beta*delta));
  }

  //for(int i=0; i < len; i++) {
  //  k0 = exp(beta * z[i]);
  //  k1 = exp(beta * (z[i]+1));
  //  temp = gamma + (alpha/(beta*delta))*(log(beta*delta*k1 + beta) -
  //                                       log(beta*delta*k0 + beta));
  //  res[i] = 1 - exp(-1.0*temp);
  //}

  return(res);

}


