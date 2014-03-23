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

RcppExport SEXP mortalityhazard_gompertz_cpp(SEXP theta, SEXP z)
{

  BEGIN_RCPP

  using namespace Rcpp;

  NumericVector xtheta(theta);
  NumericVector xz(z);

  int len = xz.size();
  NumericVector res(len);

  double alpha = exp(xtheta[0]);
  double beta = xtheta[1];

  for(int i=0; i < len; i++) {
    res[i] = alpha*exp(beta*xz[i]);
  }

  return(res);

  END_RCPP
}

RcppExport SEXP mortalityhazard_to_prob_gompertz_cpp(SEXP theta, SEXP z)
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
    temp = (alpha/beta)*(exp(beta*(xz[i]+1))-exp(beta*xz[i]));
    res[i] = 1 - exp(-1*temp);
  }

  return(res);

  END_RCPP

}


