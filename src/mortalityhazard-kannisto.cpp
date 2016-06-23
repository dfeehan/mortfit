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

#define ALPHA_IDX 0
#define BETA_IDX 1
#define GAMMA_IDX 2
#define DELTA_IDX 3

using namespace Rcpp;

/*
 * CODE BELOW PASTED IN FROM SAGE
 */

double kannisto_partial_alpha_partI(double alpha, double beta, double x) {

   double kannisto_partial_alpha_partI_result;
      kannisto_partial_alpha_partI_result = -(alpha*(exp(beta*x + beta)/(alpha*(alpha*beta*exp(beta*x + beta) + beta)) - exp(beta*x)/(alpha*(alpha*beta*exp(beta*x) + beta)) + log(alpha*beta*exp(beta*x) + beta)/(pow(alpha, 2)*beta) - log(alpha*beta*exp(beta*x + beta) + beta)/(pow(alpha, 2)*beta)) - log(alpha*beta*exp(beta*x) + beta)/(alpha*beta) + log(alpha*beta*exp(beta*x + beta) + beta)/(alpha*beta))*exp(-alpha*(-log(alpha*beta*exp(beta*x) + beta)/(alpha*beta) + log(alpha*beta*exp(beta*x + beta) + beta)/(alpha*beta)))/(-1 + exp(-alpha*(-log(alpha*beta*exp(beta*x) + beta)/(alpha*beta) + log(alpha*beta*exp(beta*x + beta) + beta)/(alpha*beta))));
         return kannisto_partial_alpha_partI_result;

}

double kannisto_partial_alpha_partII(double alpha, double beta, double x) {

   double kannisto_partial_alpha_partII_result;
      kannisto_partial_alpha_partII_result = -alpha*(exp(beta*x + beta)/(alpha*(alpha*beta*exp(beta*x + beta) + beta)) - exp(beta*x)/(alpha*(alpha*beta*exp(beta*x) + beta)) + log(alpha*beta*exp(beta*x) + beta)/(pow(alpha, 2)*beta) - log(alpha*beta*exp(beta*x + beta) + beta)/(pow(alpha, 2)*beta)) + log(alpha*beta*exp(beta*x) + beta)/(alpha*beta) - log(alpha*beta*exp(beta*x + beta) + beta)/(alpha*beta);
         return kannisto_partial_alpha_partII_result;

}

double kannisto_partial_beta_partI(double alpha, double beta, double x) {

   double kannisto_partial_beta_partI_result;
      kannisto_partial_beta_partI_result = -alpha*((alpha*beta*(x + 1)*exp(beta*x + beta) + alpha*exp(beta*x + beta) + 1)/(alpha*beta*(alpha*beta*exp(beta*x + beta) + beta)) - (alpha*beta*x*exp(beta*x) + alpha*exp(beta*x) + 1)/(alpha*beta*(alpha*beta*exp(beta*x) + beta)) + log(alpha*beta*exp(beta*x) + beta)/(alpha*pow(beta, 2)) - log(alpha*beta*exp(beta*x + beta) + beta)/(alpha*pow(beta, 2)))*exp(-alpha*(-log(alpha*beta*exp(beta*x) + beta)/(alpha*beta) + log(alpha*beta*exp(beta*x + beta) + beta)/(alpha*beta)))/(-1 + exp(-alpha*(-log(alpha*beta*exp(beta*x) + beta)/(alpha*beta) + log(alpha*beta*exp(beta*x + beta) + beta)/(alpha*beta))));
         return kannisto_partial_beta_partI_result;

}

double kannisto_partial_beta_partII(double alpha, double beta, double x) {

   double kannisto_partial_beta_partII_result;
      kannisto_partial_beta_partII_result = -alpha*((alpha*beta*(x + 1)*exp(beta*x + beta) + alpha*exp(beta*x + beta) + 1)/(alpha*beta*(alpha*beta*exp(beta*x + beta) + beta)) - (alpha*beta*x*exp(beta*x) + alpha*exp(beta*x) + 1)/(alpha*beta*(alpha*beta*exp(beta*x) + beta)) + log(alpha*beta*exp(beta*x) + beta)/(alpha*pow(beta, 2)) - log(alpha*beta*exp(beta*x + beta) + beta)/(alpha*pow(beta, 2)));
         return kannisto_partial_beta_partII_result;

}


/*
 * END CODE PASTED IN FROM SAGE
 */


// [[Rcpp::export]]
NumericVector mortalityhazard_kannisto_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = exp(theta[1]);
  //double beta = theta[1];

  double num = 0.0;
  double denom = 0.0;

  for(int i=0; i < len; i++) {
    num = alpha*exp(beta*z[i]);
    denom = 1 + alpha*exp(beta*z[i]);
    res[i] = num/denom;
    if (res[i] < 0) {
      res[i] = NA_INTEGER;
    }
  }

  return(res);

}

// [[Rcpp::export]]
NumericVector mortalityhazard_to_prob_kannisto_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = exp(theta[0]);
  double beta = exp(theta[1]);
  //double beta = theta[1];

  double temp = 0.0;

  for(int i=0; i < len; i++) {
    temp = (alpha * exp(beta * z[i]) + 1.0) / (alpha*exp(beta*(z[i]+1.0)) + 1.0);
    res[i] = 1.0 - pow(temp, 1.0/beta);
  }

  return(res);

}


// [[Rcpp::export]]
NumericVector mortalityhazard_kannisto_binomial_grad_cpp(NumericVector theta,
                                                         NumericVector ages,
                                                         NumericVector Dx,
                                                         NumericVector Nx)
{

  int theta_len = theta.size();
  int ages_len = ages.size();

  NumericVector res(theta_len);

  double alpha = exp(theta[0]);
  double beta = exp(theta[1]);

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
      partI += (double)Dx[ageidx] * kannisto_partial_alpha_partI(alpha, beta, cur_age);
      partII += (double)(Nx[ageidx] - Dx[ageidx]) * kannisto_partial_alpha_partII(alpha, beta, cur_age);

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
      partI += (double)Dx[ageidx] * kannisto_partial_beta_partI(alpha, beta, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * kannisto_partial_beta_partII(alpha, beta, cur_age);

  }

  res[BETA_IDX] = partI + partII;

  return(res);

}

