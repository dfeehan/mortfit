/************************************************
 * mortalityhazard-weibull.cpp
 *
 *
 *
 * see "Writing a package that uses Rcpp"
 * by Edelbuettel and Francois, Sept 29 2011
 *
 * dennis, dec 2011
 ************************************************/

#include <Rcpp.h>

#define ALPHA_IDX 0
#define BETA_IDX 1
#define GAMMA_IDX 2
#define DELTA_IDX 3

using namespace Rcpp;

/*
 * CODE BELOW PASTED IN FROM SAGE
 */

double weibull_partial_alpha_partI(double alpha, double beta, double x) {

   double weibull_partial_alpha_partI_result;
      weibull_partial_alpha_partI_result = -(-pow(x, beta + 1)/(beta + 1) + pow(x + 1, beta + 1)/(beta + 1))*exp(-alpha*(-pow(x, beta + 1)/(beta + 1) + pow(x + 1, beta + 1)/(beta + 1)))/(-1 + exp(-alpha*(-pow(x, beta + 1)/(beta + 1) + pow(x + 1, beta + 1)/(beta + 1))));
         return weibull_partial_alpha_partI_result;

}

double weibull_partial_alpha_partII(double beta, double x) {

   double weibull_partial_alpha_partII_result;
      weibull_partial_alpha_partII_result = pow(x, beta + 1)/(beta + 1) - pow(x + 1, beta + 1)/(beta + 1);
         return weibull_partial_alpha_partII_result;

}

double weibull_partial_beta_partI(double alpha, double beta, double x) {

   double weibull_partial_beta_partI_result;
      weibull_partial_beta_partI_result = -alpha*(-pow(x, beta + 1)*log(x)/(beta + 1) + pow(x, beta + 1)/pow(beta + 1, 2) + pow(x + 1, beta + 1)*log(x + 1)/(beta + 1) - pow(x + 1, beta + 1)/pow(beta + 1, 2))*exp(-alpha*(-pow(x, beta + 1)/(beta + 1) + pow(x + 1, beta + 1)/(beta + 1)))/(-1 + exp(-alpha*(-pow(x, beta + 1)/(beta + 1) + pow(x + 1, beta + 1)/(beta + 1))));
         return weibull_partial_beta_partI_result;

}

double weibull_partial_beta_partII(double alpha, double beta, double x) {

   double weibull_partial_beta_partII_result;
      weibull_partial_beta_partII_result = -alpha*(-pow(x, beta + 1)*log(x)/(beta + 1) + pow(x, beta + 1)/pow(beta + 1, 2) + pow(x + 1, beta + 1)*log(x + 1)/(beta + 1) - pow(x + 1, beta + 1)/pow(beta + 1, 2));
         return weibull_partial_beta_partII_result;

}


/*
 * END CODE PASTED IN FROM SAGE
 */

// [[Rcpp::export]]
NumericVector mortalityhazard_weibull_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = theta[1];

  for(int i=0; i < len; i++) {

    res[i] = alpha*pow(z[i], beta);

    if (res[i] < 0) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

}

// [[Rcpp::export]]
NumericVector mortalityhazard_to_prob_weibull_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = theta[1];

  double temp1 = 0.0;
  double temp2 = 0.0;

  for(int i=0; i < len; i++) {
    temp1 = (alpha / (1.0 + beta)) * pow(z[i],(1.0+beta));
    temp2 = -1.0 * (alpha / (1.0 + beta)) * pow((z[i]+1),(1.0+beta));
    res[i] = 1 - (exp(temp1)*exp(temp2));
  }

  return(res);

}


// [[Rcpp::export]]
NumericVector mortalityhazard_weibull_binomial_grad_cpp(NumericVector theta,
                                                        NumericVector ages,
                                                        NumericVector Dx,
                                                        NumericVector Nx)
{

  int theta_len = theta.size();
  int ages_len = ages.size();

  NumericVector res(theta_len);

  double alpha = exp(theta[0]);
  double beta = theta[1];

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
      partI += (double)Dx[ageidx] * weibull_partial_alpha_partI(alpha, beta, cur_age);
      partII += (double)(Nx[ageidx] - Dx[ageidx]) * weibull_partial_alpha_partII(beta, cur_age);

  }

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
      partI += (double)Dx[ageidx] * weibull_partial_beta_partI(alpha, beta, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * weibull_partial_beta_partII(alpha, beta, cur_age);

  }

  res[BETA_IDX] = partI + partII;

  return(res);

}

