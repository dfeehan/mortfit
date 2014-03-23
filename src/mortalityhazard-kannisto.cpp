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

RcppExport SEXP mortalityhazard_kannisto_cpp(SEXP theta, SEXP z)
{

  BEGIN_RCPP

  using namespace Rcpp;

  NumericVector xtheta(theta);
  NumericVector xz(z);

  int len = xz.size();
  NumericVector res(len);

  double alpha = exp(xtheta[0]);
  double beta = xtheta[1];

  double num = 0.0;
  double denom = 0.0;

  for(int i=0; i < len; i++) {
    num = alpha*exp(beta*xz[i]);
    denom = 1 + alpha*exp(beta*xz[i]);
    res[i] = num/denom;
    if (res[i] < 0) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

  END_RCPP
}

RcppExport SEXP mortalityhazard_to_prob_kannisto_cpp(SEXP theta, SEXP z)
{

  BEGIN_RCPP

  using namespace Rcpp;

  NumericVector xtheta(theta);
  NumericVector xz(z);

  int len = xz.size();
  NumericVector res(len);

  double alpha = exp(xtheta[0]);
  double beta = xtheta[1];
  
  double temp = 0.0;

  for(int i=0; i < len; i++) {
    temp = (1/beta)*(log1p(alpha*exp(beta*(xz[i]+1)))-log1p(alpha*exp(beta*xz[i])));
    res[i] = 1 - exp(-1*temp);
  }

  return(res);

  END_RCPP

}


