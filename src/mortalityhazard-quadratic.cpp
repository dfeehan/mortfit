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

#define ALPHA_IDX 0
#define BETA_IDX 1
#define GAMMA_IDX 2
#define DELTA_IDX 3

using namespace Rcpp;

/*
 * helper fn: calculate erf on a single complex number
 *   (RcppFaddeeva::erf wants a vector of complex numbers)
 */
std::complex<double> vec_erf(std::complex<double> x) {

    std::vector< std::complex<double> > arglist, reslist;

    arglist.push_back(x);

    reslist = RcppFaddeeva::erf(arglist);

    return(reslist[0]);

}

/*
 * CODE BELOW PASTED IN FROM SAGE
 * (but then edited quite a bit so that it could call RcppFaddeeva::erf
 *  and more generally handle complex-valued intermediate quantities)
 */

double log_quadratic_partial_alpha_partI(double a, double b, double g, double x) {

   std::complex<double> alpha(a, 0.0), beta(b, 0.0), gamma(g, 0.0), res;

   res =
      (1.0/2.0)*(-sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma)) /
      sqrt(-gamma) + sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma)) /
      sqrt(-gamma))*exp(-1.0/2.0*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/
      sqrt(-gamma) + (1.0/2.0)*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      sqrt(-gamma))/(exp(-1.0/2.0*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/
      sqrt(-gamma) + (1.0/2.0)*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      sqrt(-gamma)) - 1.0);

     //Rcout << "partial wrt alpha, part I: " << res << std::endl;

    return(res.real());
}

double log_quadratic_partial_alpha_partII(double a, double b, double g, double x) {

   std::complex<double> alpha(a, 0.0), beta(b, 0.0), gamma(g, 0.0), res;

   res =
      -1.0/2.0*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma)) /sqrt(-gamma) +
      (1.0/2.0)*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      sqrt(-gamma);


   // Rcout << "partial wrt alpha, part I: " << res << std::endl;

   return(res.real());

}

double log_quadratic_partial_beta_partI(double a, double b, double g, double x) {

   std::complex<double> alpha(a, 0.0), beta(b, 0.0), gamma(g, 0.0), res;

   res =
      -1.0/4.0*(-sqrt(M_PI)*beta*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/
      (gamma*sqrt(-gamma)) + sqrt(M_PI)*beta*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      (gamma*sqrt(-gamma)) - 2.0*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma +
      (1.0/4.0)*pow(beta + 2.0*gamma*x, 2.0)/gamma)/gamma +
      2.0*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma +
      (1.0/4.0)*pow(beta + 2.0*gamma*x + 2.0*gamma, 2.0)/gamma)/gamma)*
      exp(-1.0/2.0*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/
      sqrt(-gamma) + (1.0/2.0)*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      sqrt(-gamma))/(exp(-1.0/2.0*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/
      sqrt(-gamma) + (1.0/2.0)*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      sqrt(-gamma)) - 1.0);


  // Rcout << "partial wrt beta, part I: " << res << std::endl;

  return(res.real());
}

double log_quadratic_partial_beta_partII(double a, double b, double g, double x) {

   std::complex<double> alpha(a, 0.0), beta(b, 0.0), gamma(g, 0.0), res;

   res =
      (1.0/4.0)*sqrt(M_PI)*beta*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/
      (gamma*sqrt(-gamma)) - 1.0/4.0*sqrt(M_PI)*beta*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      (gamma*sqrt(-gamma)) + (1.0/2.0)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma +
      (1.0/4.0)*pow(beta + 2.0*gamma*x, 2.0)/gamma)/gamma -
      1.0/2.0*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma +
      (1.0/4.0)*pow(beta + 2.0*gamma*x + 2.0*gamma, 2.0)/gamma)/gamma;

   // Rcout << "partial wrt beta, part II: " << res << std::endl;

   return res.real();
}

double log_quadratic_partial_gamma_partI(double a, double b, double g, double x) {

   std::complex<double> alpha(a, 0.0), beta(b, 0.0), gamma(g, 0.0), res;

   res =
      (1.0/8.0)*(-sqrt(M_PI)*pow(beta, 2.0)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf( (1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/
      (pow(gamma, 2.0)*sqrt(-gamma)) +
      sqrt(M_PI)*pow(beta, 2.0)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      (pow(gamma, 2.0)*sqrt(-gamma)) - 2.0*(4.0*x/sqrt(-gamma)
      + (beta + 2.0*gamma*x)/pow(-gamma, 3.0/2.0))*exp(alpha -
      1.0/4.0*pow(beta, 2.0)/gamma + (1.0/4.0)*pow(beta + 2.0*gamma*x,
      2.0)/gamma)/sqrt(-gamma) + 2.0*(4.0*(x + 1.0)/sqrt(-gamma) + (beta + 2.0*gamma*x +
      2.0*gamma)/pow(-gamma, 3.0/2.0))*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma
      + (1.0/4.0)*pow(beta + 2.0*gamma*x + 2.0*gamma, 2.0)/gamma)/sqrt(-gamma) -
      2.0*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/
      pow(-gamma, 3.0/2.0) + 2.0*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      pow(-gamma, 3.0/2.0))*exp(-1.0/2.0*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/
      sqrt(-gamma) + (1.0/2.0)*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      sqrt(-gamma))/(exp(-1.0/2.0*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/
      sqrt(-gamma) + (1.0/2.0)*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      sqrt(-gamma)) - 1.0);

  // Rcout << "partial wrt gamma, part I: " << res << std::endl;

      return res.real();

}

double log_quadratic_partial_gamma_partII(double a, double b, double g, double x) {

   std::complex<double> alpha(a, 0.0), beta(b, 0.0), gamma(g, 0.0), res;

   res =
      -1.0/8.0*sqrt(M_PI)*pow(beta, 2.0)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/ (pow(gamma, 2.0)*sqrt(-gamma)) +
      (1.0/8.0)*sqrt(M_PI)*pow(beta, 2.0)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      (pow(gamma, 2.0)*sqrt(-gamma)) - 1.0/4.0*(4.0*x/sqrt(-gamma) +
      (beta + 2.0*gamma*x)/pow(-gamma, 3.0/2.0))*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma +
      (1.0/4.0)*pow(beta + 2.0*gamma*x, 2.0)/gamma)/sqrt(-gamma) +
      (1.0/4.0)*(4.0*(x + 1.0)/sqrt(-gamma) + (beta + 2.0*gamma*x + 2.0*gamma)/
      pow(-gamma, 3.0/2.0))*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma
      + (1.0/4.0)*pow(beta + 2.0*gamma*x + 2.0*gamma, 2.0)/gamma)/sqrt(-gamma) -
      1.0/4.0*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x)/sqrt(-gamma))/
      pow(-gamma, 3.0/2.0) + (1.0/4.0)*sqrt(M_PI)*exp(alpha - 1.0/4.0*pow(beta, 2.0)/gamma)*
      vec_erf((1.0/2.0)*(beta + 2.0*gamma*x + 2.0*gamma)/sqrt(-gamma))/
      pow(-gamma, 3.0/2.0);

  // Rcout << "partial wrt gamma, part II: " << res << std::endl;

  return res.real();
}

/*
 * END CODE PASTED IN FROM SAGE
 */

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


// [[Rcpp::export]]
NumericVector mortalityhazard_quadratic_binomial_grad_cpp(NumericVector theta,
                                                          NumericVector ages,
                                                          NumericVector Dx,
                                                          NumericVector Nx)
{

  int theta_len = theta.size();
  int ages_len = ages.size();

  NumericVector res(theta_len);

  double alpha = theta[0];
  double beta = theta[1];
  double gamma = theta[2];

  double partI = 0.0, partII =0.0;

  double cur_age = 0.0;

  // initialize result vector to 0
  for(int i=0; i < theta_len; i++) {
    res[i] = 0.0;
  }

  /*
   * alpha
   */

  partI = 0.0;
  partII = 0.0;

  for(int ageidx=0; ageidx < ages_len; ageidx++) {

      cur_age = ages[ageidx];

      // for each age, we want
      // D_z * partI + (N_z - D_z) * partII
      partI += (double)Dx[ageidx] * log_quadratic_partial_alpha_partI(alpha, beta, gamma, cur_age);
      partII += (double)(Nx[ageidx] - Dx[ageidx]) * log_quadratic_partial_alpha_partII(alpha, beta, gamma,cur_age);

  }

  //Rcout << "partial wrt alpha: " << partI + partII << std::endl;

  res[ALPHA_IDX] = partI + partII;

  /*
   * beta
   */

  partI = 0.0;
  partII = 0.0;

  for(int ageidx=0; ageidx < ages_len; ageidx++) {

      cur_age = ages[ageidx];

      // for each age, we want
      // D_z * partI + (N_z - D_z) * partII
      partI += (double)Dx[ageidx] * log_quadratic_partial_beta_partI(alpha, beta, gamma, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * log_quadratic_partial_beta_partII(alpha,
                                                                                 beta,
                                                                                 gamma,
                                                                                 cur_age);

  }

  //Rcout << "partial wrt beta: " << partI + partII << std::endl;

  res[BETA_IDX] = partI + partII;

  /*
   * gamma
   */

  partI = 0.0;
  partII = 0.0;

  for(int ageidx=0; ageidx < ages_len; ageidx++) {

      cur_age = ages[ageidx];

      // for each age, we want
      // D_z * partI + (N_z - D_z) * partII
      partI += (double)Dx[ageidx] * log_quadratic_partial_gamma_partI(alpha, beta, gamma, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * log_quadratic_partial_gamma_partII(alpha,
                                                                                      beta,
                                                                                      gamma,
                                                                                      cur_age);

  }

  //Rcout << "partial wrt gamma: " << partI + partII << std::endl;

  res[GAMMA_IDX] = partI + partII;

  return(res);

}

