/************************************************
 * mortalityhazard-quadratic.cpp
 *
 *
 *
 * see "Writing a package that uses Rcpp"
 * by Edelbuettel and Francois, Sept 29 2011
 *
 * dennis, dec 2011
 ************************************************/

// [[Rcpp::depends(RcppFaddeeva)]]

#include <complex>
#include <Rcpp.h>
#include <RcppFaddeeva.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mortalityhazard_quadratic_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();

  double alpha = theta[0];
  double beta = theta[1];
  double gamma = theta[2];

  NumericVector res(len);

  for(int i=0; i < len; i++) {
    res[i] = exp(alpha + beta*z[i] + gamma*z[i]*z[i]);
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
NumericVector mortalityhazard_to_prob_quadratic_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  const double pi = 3.1415926535897;

  double alpha = theta[0];
  double beta = theta[1];

  std::complex<double> gamma = theta[2];
  std::complex<double> sqrt_minus_gamma = sqrt(-1.0*gamma);

  std::complex<double> k0;
  k0 = exp(alpha - (pow(beta,2.0)/(4.0*gamma)))*(sqrt(pi) / (2.0 * sqrt_minus_gamma));

  // if we're getting Inf or -Inf from this exponential,
  // return NA...
  if (real(k0) == R_PosInf || real(k0) == R_NegInf) {
    for(int i=0; i<len; i++) {
      res[i] = NA_INTEGER;
    }
    return(res);
  }

  std::complex<double> temp;

  std::vector<std::complex<double> > arg1(len), arg2(len);

  for(int i=0; i < len; i++) {

    arg1[i] = (beta + (2.0*gamma*(z[i]+1)))/(2.0*sqrt_minus_gamma);
    arg2[i] = (beta + (2.0*gamma*z[i]))/(2.0*sqrt_minus_gamma);

  }

  arg1 = RcppFaddeeva::erf(arg1);
  arg2 = RcppFaddeeva::erf(arg2);

  //Rprintf("log-quadratic pis:\n");
  for(int i=0; i < len; i++) {

    temp = k0 * (arg2[i] - arg1[i]);

    //Rprintf("%f + i%f / ", real(arg2[i]), imag(arg2[i]));

    temp = 1.0 - exp(-1.0 * temp);

    //Rprintf("%f + i%f / ", real(temp), imag(temp));

    res[i] = real(temp);

    //Rprintf("%f / ", res[i]);
  }

  //Rprintf("\n");

  return(res);

}

