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

RcppExport SEXP mortalityhazard_makeham_cpp(SEXP theta, SEXP z)
{

  BEGIN_RCPP

  using namespace Rcpp;

  NumericVector xtheta(theta);
  NumericVector xz(z);

  int len = xz.size();
  int parlen = xtheta.size();

  if (parlen != 3) {
    
  }

  NumericVector res(len);

  double alpha = exp(xtheta[0]);
  double beta = xtheta[1];
  double gamma = exp(xtheta[2]);

  for(int i=0; i < len; i++) {

    res[i] = gamma + alpha*exp(beta*xz[i]);

    if (res[i] < 0) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

  END_RCPP
}

RcppExport SEXP mortalityhazard_to_prob_makeham_cpp(SEXP theta, SEXP z)
{

  BEGIN_RCPP

  using namespace Rcpp;

  NumericVector xtheta(theta);
  NumericVector xz(z);

  int len = xz.size();
  NumericVector res(len);

  double alpha = exp(xtheta[0]);
  double beta = xtheta[1];
  double gamma = exp(xtheta[2]);

  double temp = 0.0;

  for(int i=0; i < len; i++) {
    temp = (alpha/beta)*(exp(beta*(xz[i]+1))-exp(beta*xz[i])) + gamma;
    res[i] = 1 - exp(-1*temp);
  }

  return(res);

  END_RCPP
}


