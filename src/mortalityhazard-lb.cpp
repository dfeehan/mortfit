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

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mortalityhazard_lb_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = theta[0];
  double beta = exp(theta[1]);
  double gamma = exp(theta[2]);
  double delta = theta[3];

  //Rprintf("alpha = %f; ", alpha);
  //Rprintf("beta = %f; ", beta);
  //Rprintf("gamma = %f; ", gamma);
  //Rprintf("delta = %f\n", delta);

  for(int i=0; i < len; i++) {

    res[i] = alpha + beta*atan(gamma*(z[i]-delta));

    if (res[i] < 0) {
      //res[i] = NA_INTEGER;
      res[i] = R_NegInf;
    }
    //Rprintf("res[%d] = %f\n", i, res[i]);
  }

  return(res);

}

// [[Rcpp::export]]
NumericVector mortalityhazard_to_prob_lb_cpp(NumericVector theta, NumericVector z)
{

  int len = z.size();
  NumericVector res(len);

  double alpha = theta[0];
  double beta = exp(theta[1]);
  double gamma = exp(theta[2]);
  double delta = theta[3];

  double k0 = 0;
  double k1 = 0;

  double res1 = -alpha;
  double res2 = (-beta/(2*gamma));
  double res3 = 0.0;
  double res4 = 0.0;
  double res5 = 0.0;
  double res6 = 0.0;
  double temp = 0.0;

  //double res1 = alpha;
  //double res2 = (1/(2*gamma))*beta;
  //double res3 = 0.0;
  //double res4 = 0.0;
  //double res5 = 0.0;
  //double res6 = 0.0;
  //double temp = 0.0;

  for(int i=0; i < len; i++) {

    k0 = gamma * (z[i] - delta);
    k1 = gamma * (z[i] - delta + 1);

    res3 = 2.0 * k1 * atan(k1);
    res4 = 2.0 * k0 * atan(k0);
    res5 = log1p(pow(k0,2));
    res6 = log1p(pow(k1,2));

    temp = res1 + res2*(res3 - res4 + res5 - res6);

    //res3 = 2*gamma*(z[i]-delta)*atan(gamma*(delta-z[i]));
    //res4 = 2*gamma*(delta-z[i]-1)*atan(gamma*(delta-z[i]-1));
    //res5 = log1p(pow(gamma,2)*(pow(delta-z[i],2)));
    //res6 = log1p(pow(gamma,2)*(pow(delta-z[i]-1,2)));
    //temp = res1 + res2*(res3 + res4 + res5 - res6);

    res[i] = 1 - exp(temp);

  }

  return(res);

}


