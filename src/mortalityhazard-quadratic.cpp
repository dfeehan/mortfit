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

#include <complex>
#include <Rcpp.h>

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

std::complex<double> erfi_approx_nr1(std::complex<double> x)
{

  std::complex<double> erfiof1 (1.65043, 0.0);
  std::complex<double> expof1 (exp(1.0), 0.0);
  std::complex<double> sqrtofpi (sqrt(3.1415926535897), 0.0);

  std::complex<double> res (0.0, 0.0);

  res = (173465337818777746628.0*expof1 - 161189343277473823463.0*sqrtofpi*erfiof1 + x*(-181236135883415078228.0*expof1 + 58207050515281020997.0*sqrtofpi*erfiof1 +
x*(-32990023275939643620.0*expof1 + 70448356034529178041.0*sqrtofpi*erfiof1 +
x*(53037366273114757940.0*expof1 - 37842691311251513855.0*sqrtofpi*erfiof1 +
x*(-7666360342517050580.0*expof1 - 5516314131408582425.0*sqrtofpi*erfiof1 +
x*(-6665334075150129708.0*expof1 + 7528886304748135803.0*sqrtofpi*erfiof1 +
x*(2371517339152075268.0*expof1 - 1823214457143623537.0*sqrtofpi*erfiof1 +
x*(-316367854022677700.0*expof1 + 146437170211591039.0*sqrtofpi*erfiof1))))))))/
   (-161189343277473823463.0*sqrtofpi +
     x*(58207050515281020997.0*sqrtofpi +
        x*(70448356034529178041.0*sqrtofpi +
           x*(-37842691311251513855.0*sqrtofpi +
              x*(-5516314131408582425.0*sqrtofpi +
                 x*(7528886304748135803.0*sqrtofpi +
                    x*(-1823214457143623537.0*sqrtofpi +
                       146437170211591039.0*sqrtofpi*x)))))));

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
  std::complex<double> sqrt_gamma = sqrt(gamma);

  std::complex<double> k0;
  k0 = exp(alpha - (pow(beta,2.0)/(4.0*gamma)))*sqrt(pi);

  // if we're getting Inf or -Inf from this exponential,
  // return NA...
  if (real(k0) == R_PosInf || real(k0) == R_NegInf) {
    for(int i=0; i<len; i++) {
      res[i] = NA_INTEGER;
    }
    return(res);
  }

  std::complex<double> temp;

  for(int i=0; i < len; i++) {
    temp = (k0/(2.0*sqrt_gamma)) *
      (erfi_approx_nr1((beta + (2.0*gamma*(z[i]+1)))/(2.0*sqrt_gamma)) -
       erfi_approx_nr1((beta + (2.0*gamma*z[i]))/(2.0*sqrt_gamma)));

    // maddeningly, this:
    // 1.0 - exp(-1.0*real(temp))
    // did *not* work!
    res[i] = real(1.0 - exp(-1.0* temp));

  }

  return(res);

}


/*
// use euler's eqn b/c can't find complex-valued exponential
// in Rcpp...
Rcomplex complex_exp(Rcomplex arg)
{

  Rcomplex ival;
  ival.r = 0;
  ival.i = 1;

  //return(exp(arg.r)*(cos(arg.i) + ival*sin(arg.i));
  return ival;

}

// these are
// based on Rcomplex operators from Rcpp package;
Rcomplex operator*( const Rcomplex& lhs, double rhs){
  Rcomplex y ;
  y.r = lhs.r * rhs ;
  return y ;
}

Rcomplex operator*( double lhs, const Rcomplex& rhs){
  return(rhs*lhs);
}

Rcomplex operator+( const Rcomplex& lhs, double rhs){
  Rcomplex y ;
  y.r = lhs.r + rhs ;
  return y ;
}

//Rcomplex operator+( double lhs, const Rcomplex& rhs){
//  return(rhs-lhs);
//}

Rcomplex operator/( double lhs, const Rcomplex& rhs){
  Rcomplex lhs_complex;
  lhs_complex.i = 0;
  lhs_complex.r = lhs;
  return(lhs/rhs);
}


Rcomplex OLDerfi_approx_nr1(Rcomplex x)
{

  Rcomplex erfiof1, expof1, sqrtofpi;

  erfiof1.i = 0.0;
  expof1.i = 0.0;
  sqrtofpi.i = 0.0;

  erfiof1.r = 1.65043;
  expof1.r = exp(1.0);
  sqrtofpi.r = sqrt(3.1415926535897);

  Rcomplex res;
  res.r = res.i = 0.0;

  res = (173465337818777746628.0*expof1 - 161189343277473823463.0*sqrtofpi*erfiof1 + x*(-181236135883415078228.0*expof1 + 58207050515281020997.0*sqrtofpi*erfiof1 +
x*(-32990023275939643620.0*expof1 + 70448356034529178041.0*sqrtofpi*erfiof1 +
x*(53037366273114757940.0*expof1 - 37842691311251513855.0*sqrtofpi*erfiof1 +
x*(-7666360342517050580.0*expof1 - 5516314131408582425.0*sqrtofpi*erfiof1 +
x*(-6665334075150129708.0*expof1 + 7528886304748135803.0*sqrtofpi*erfiof1 +
x*(2371517339152075268.0*expof1 - 1823214457143623537.0*sqrtofpi*erfiof1 +
x*(-316367854022677700.0*expof1 + 146437170211591039.0*sqrtofpi*erfiof1))))))))/
   (-161189343277473823463.0*sqrtofpi +
     x*(58207050515281020997.0*sqrtofpi +
        x*(70448356034529178041.0*sqrtofpi +
           x*(-37842691311251513855.0*sqrtofpi +
              x*(-5516314131408582425.0*sqrtofpi +
                 x*(7528886304748135803.0*sqrtofpi +
                    x*(-1823214457143623537.0*sqrtofpi +
                       146437170211591039.0*sqrtofpi*x)))))));

  return(res);

}

RcppExport SEXP OLDmortalityhazard_to_prob_quadratic_cpp(SEXP theta, SEXP z)
{
  BEGIN_RCPP

  using namespace Rcpp;

  NumericVector xtheta(theta);
  NumericVector xz(z);

  int len = xz.size();
  NumericVector res(len);

  const double pi = 3.1415926535897;

  double alpha = xtheta[0];
  double beta = xtheta[1];
  double gamma = xtheta[2];

  Rcomplex sqrt_gamma;

  if (gamma < 0) {
    sqrt_gamma.r = 0;
    sqrt_gamma.i = sqrt(gamma);
  } else {
    sqrt_gamma.r = sqrt(gamma);
    sqrt_gamma.i = 0;
  }

  Rcomplex k0;
  k0 = alpha - (pow(beta,2)/(4*gamma));

  // TODO LEFT OFF:
  // EXP FUNCTION ISN'T WORKING HERE...

  k0 = complex_exp(k0)*sqrt(pi);

  Rcomplex temp;

  for(int i=0; i < len; i++) {
    temp = (k0/(2*sqrt_gamma)) *
      (erfi_approx_nr1((beta + (2*gamma*(xz[i]+1)))/(2*sqrt_gamma)) -
       erfi_approx_nr1((beta + (2*gamma*xz[i]))/(2*sqrt_gamma)));

    res[i] = 1 - exp(-1* temp.r);

  }

  return(res);

  END_RCPP
}
*/
