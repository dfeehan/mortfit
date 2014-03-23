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

RcppExport SEXP mortalityhazard_logistic_cpp(SEXP theta, SEXP z)
{

  BEGIN_RCPP

  using namespace Rcpp;

  NumericVector xtheta(theta);
  NumericVector xz(z);

  int len = xz.size();
  NumericVector res(len);

  double alpha = exp(xtheta[0]);
  double beta = exp(xtheta[1]);
  double gamma = exp(xtheta[2]);
  double delta = exp(xtheta[3]);

  double k = 0.0;
  double num = 0.0;
  double denom = 0.0;

  for(int i=0; i < len; i++) {
    k = exp(beta * xz[i]);
    num = alpha * k;
    denom = 1 + (delta*k);
    res[i] = (num/denom) + gamma;
  }

  bool anyNA = is_true( any( is_na(res) ) );

  if (anyNA) {
    for(int i=0; i<len; i++) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

  END_RCPP

}

RcppExport SEXP mortalityhazard_to_prob_logistic_cpp(SEXP theta, SEXP z)
{

  BEGIN_RCPP

  using namespace Rcpp;

  NumericVector xtheta(theta);
  NumericVector xz(z);

  int len = xz.size();
  NumericVector res(len);

  double alpha = exp(xtheta[0]);
  double beta = exp(xtheta[1]);
  double gamma = exp(xtheta[2]);
  double delta = exp(xtheta[3]);

  double k0 = 0.0;
  double k1 = 0.0;
  double temp = 0.0;

  for(int i=0; i < len; i++) {
    k0 = exp(beta * xz[i]);
    k1 = exp(beta * (xz[i]+1));
    temp = gamma + (alpha/(beta*delta))*(log1p(delta*k1) - log1p(delta*k0));
    res[i] = 1 - exp(-1*temp);
  }

  return(res);

  END_RCPP

}


