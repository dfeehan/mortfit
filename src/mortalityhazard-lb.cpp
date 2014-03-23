/************************************************
 * mortalityhazard-lb.cpp
 * 
 * 
 * 
 * see "Writing a package that uses Rcpp"
 * by Edelbuettel and Francois, Sept 29 2011
 * 
 * dennis, dec 2011
 ************************************************/

#include <Rcpp.h>

RcppExport SEXP mortalityhazard_lb_cpp(SEXP theta, SEXP z)
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

  for(int i=0; i < len; i++) {

    res[i] = alpha + beta*atan(gamma*(xz[i]-delta));

    if (res[i] < 0) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

  END_RCPP
}

RcppExport SEXP mortalityhazard_to_prob_lb_cpp(SEXP theta, SEXP z)
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

  double res1 = alpha;
  double res2 = (1/(2*gamma))*beta;
  double res3 = 0.0;
  double res4 = 0.0;
  double res5 = 0.0;
  double res6 = 0.0;
  double temp = 0.0;

  for(int i=0; i < len; i++) {

    res3 = 2*gamma*(xz[i]-delta)*atan(gamma*(delta-xz[i]));
    res4 = 2*gamma*(delta-xz[i]-1)*atan(gamma*(delta-xz[i]-1));
    res5 = log1p(pow(gamma,2)*(pow(delta-xz[i],2)));
    res6 = log1p(pow(gamma,2)*(pow(delta-xz[i]-1,2)));

    temp = res1 + res2*(res3 + res4 + res5 - res6);

    res[i] = 1 - exp(-1*temp);
  }

  return(res);

  END_RCPP
}


