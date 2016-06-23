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

#define ALPHA_IDX 0
#define BETA_IDX 1
#define GAMMA_IDX 2
#define DELTA_IDX 3

using namespace Rcpp;

/*
 * CODE BELOW PASTED IN FROM SAGE
 */

double lb_partial_alpha_partI(double alpha, double beta, double delta, double gamma, double x) {

   double lb_partial_alpha_partI_result;

   lb_partial_alpha_partI_result = -exp(beta*delta*atan(gamma*x - gamma*(delta - 1)) - beta*x*atan(delta*gamma - gamma*x) + (1.0L/2.0L)*beta*log(fabs(pow(gamma, 2)*pow(x, 2) - 2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta + 1) + 1))/gamma)/(exp(beta*delta*atan(gamma*x - gamma*(delta - 1)) - beta*x*atan(delta*gamma - gamma*x) + (1.0L/2.0L)*beta*log(fabs(pow(gamma, 2)*pow(x, 2) - 2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta + 1) + 1))/gamma) - exp(alpha - beta*delta*atan(delta*gamma - gamma*x) + beta*x*atan(gamma*x - gamma*(delta - 1)) + beta*atan(gamma*x - gamma*(delta - 1)) + (1.0L/2.0L)*beta*log(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x + pow(gamma, 2)*pow(x, 2) + 1))/gamma));

    return lb_partial_alpha_partI_result;

}

double lb_partial_alpha_partII() {

    double lb_partial_alpha_partII_result;

    lb_partial_alpha_partII_result = -1;

    return lb_partial_alpha_partII_result;

}

double lb_partial_beta_partI(double alpha, double beta, double delta, double gamma, double x) {

    double lb_partial_beta_partI_result;

    lb_partial_beta_partI_result = (1.0L/2.0L)*(2*((delta*gamma -
                    gamma*x)*exp(beta*delta*atan(gamma*x - gamma*(delta - 1)) -
                        beta*x*atan(delta*gamma - gamma*x))*atan(delta*gamma -
                        gamma*x) + (-gamma*x + gamma*(delta -
                                1))*exp(beta*delta*atan(gamma*x - gamma*(delta
                                        - 1)) - beta*x*atan(delta*gamma -
                                        gamma*x))*atan(gamma*x - gamma*(delta -
                                            1)))*exp((1.0L/2.0L)*beta*log(fabs(pow(gamma,
                                                    2)*pow(x, 2) - 2*pow(gamma,
                                                        2)*x*(delta - 1) +
                                                pow(gamma, 2)*(pow(delta, 2) -
                                                    2*delta + 1) + 1))/gamma) -
            exp(beta*delta*atan(gamma*x - gamma*(delta - 1)) -
                beta*x*atan(delta*gamma - gamma*x) +
                (1.0L/2.0L)*beta*log(fabs(pow(gamma, 2)*pow(x, 2) -
                        2*pow(gamma, 2)*x*(delta - 1) + pow(gamma,
                            2)*(pow(delta, 2) - 2*delta + 1) +
                        1))/gamma)*log(fabs(pow(delta, 2)*pow(gamma, 2) -
                        2*delta*pow(gamma, 2)*x + pow(gamma, 2)*pow(x, 2) + 1))
            + exp(beta*delta*atan(gamma*x - gamma*(delta - 1)) -
                beta*x*atan(delta*gamma - gamma*x) +
                (1.0L/2.0L)*beta*log(fabs(pow(gamma, 2)*pow(x, 2) -
                        2*pow(gamma, 2)*x*(delta - 1) + pow(gamma,
                            2)*(pow(delta, 2) - 2*delta + 1) +
                        1))/gamma)*log(fabs(pow(gamma, 2)*pow(x, 2) -
                        2*pow(gamma, 2)*x*(delta - 1) + pow(gamma,
                            2)*(pow(delta, 2) - 2*delta + 1) +
                        1)))/(gamma*exp(beta*delta*atan(gamma*x - gamma*(delta
                                - 1)) - beta*x*atan(delta*gamma - gamma*x) +
                        (1.0L/2.0L)*beta*log(fabs(pow(gamma, 2)*pow(x, 2) -
                                2*pow(gamma, 2)*x*(delta - 1) + pow(gamma,
                                    2)*(pow(delta, 2) - 2*delta + 1) +
                                1))/gamma) - gamma*exp(alpha -
                            beta*delta*atan(delta*gamma - gamma*x) +
                            beta*x*atan(gamma*x - gamma*(delta - 1)) +
                            beta*atan(gamma*x - gamma*(delta - 1)) +
                            (1.0L/2.0L)*beta*log(fabs(pow(delta, 2)*pow(gamma,
                                        2) - 2*delta*pow(gamma, 2)*x +
                                    pow(gamma, 2)*pow(x, 2) + 1))/gamma));

    return lb_partial_beta_partI_result;

}

double lb_partial_beta_partII(double delta, double gamma, double x) {

    double lb_partial_beta_partII_result;

    lb_partial_beta_partII_result = (1.0L/2.0L)*(2*(delta*gamma - gamma*x)*atan(delta*gamma - gamma*x) + 2*(-gamma*x + gamma*(delta - 1))*atan(gamma*x - gamma*(delta - 1)) - log(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x + pow(gamma, 2)*pow(x, 2) + 1)) + log(fabs(pow(gamma, 2)*pow(x, 2) - 2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta + 1) + 1)))/gamma;

    return lb_partial_beta_partII_result;

}

double lb_partial_gamma_partI(double alpha, double beta, double delta, double gamma, double x) {

   double lb_partial_gamma_partI_result;
      lb_partial_gamma_partI_result = -1.0L/2.0L*(-2*((2*beta*pow(gamma, 2)*x -
      pow(gamma, 2)*(2*beta*delta - beta))*exp(beta*delta*atan(gamma*x -
      gamma*(delta - 1)) - beta*x*atan(delta*gamma -
      gamma*x))*pow(fabs(pow(gamma, 2)*pow(x, 2) - 2*pow(gamma, 2)*x*(delta -
      1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta + 1) + 1), 2) -
      (beta*pow(gamma, 8)*pow(x, 8) - 2*pow(gamma, 8)*pow(x, 7)*(4*beta*delta -
      3*beta) + pow(gamma, 8)*(beta*pow(delta, 8) - 6*beta*pow(delta, 7) +
      15*beta*pow(delta, 6) - 20*beta*pow(delta, 5) + 15*beta*pow(delta, 4) -
      6*beta*pow(delta, 3) + beta*pow(delta, 2)) + pow(gamma,
      6)*(3*beta*pow(delta, 6) - 14*beta*pow(delta, 5) + 27*beta*pow(delta, 4)
      - 28*beta*pow(delta, 3) + 17*beta*pow(delta, 2) - 6*beta*delta + beta) +
      pow(gamma, 4)*(3*beta*pow(delta, 4) - 10*beta*pow(delta, 3) +
      13*beta*pow(delta, 2) - 8*beta*delta + 2*beta) + pow(gamma,
      2)*(beta*pow(delta, 2) - 2*beta*delta + beta) + pow(x,
      6)*(3*beta*pow(gamma, 6) + pow(gamma, 8)*(28*beta*pow(delta, 2) -
      42*beta*delta + 15*beta)) - 2*pow(x, 5)*(pow(gamma,
      8)*(28*beta*pow(delta, 3) - 63*beta*pow(delta, 2) + 45*beta*delta -
      10*beta) + pow(gamma, 6)*(9*beta*delta - 7*beta)) + pow(x,
      4)*(3*beta*pow(gamma, 4) + 5*pow(gamma, 8)*(14*beta*pow(delta, 4) -
      42*beta*pow(delta, 3) + 45*beta*pow(delta, 2) - 20*beta*delta + 3*beta) +
      pow(gamma, 6)*(45*beta*pow(delta, 2) - 70*beta*delta + 27*beta)) -
      2*pow(x, 3)*(pow(gamma, 8)*(28*beta*pow(delta, 5) - 105*beta*pow(delta,
      4) + 150*beta*pow(delta, 3) - 100*beta*pow(delta, 2) + 30*beta*delta -
      3*beta) + 2*pow(gamma, 6)*(15*beta*pow(delta, 3) - 35*beta*pow(delta, 2)
      + 27*beta*delta - 7*beta) + pow(gamma, 4)*(6*beta*delta - 5*beta)) +
      pow(x, 2)*(beta*pow(gamma, 2) + pow(gamma, 8)*(28*beta*pow(delta, 6) -
      126*beta*pow(delta, 5) + 225*beta*pow(delta, 4) - 200*beta*pow(delta, 3)
      + 90*beta*pow(delta, 2) - 18*beta*delta + beta) + pow(gamma,
      6)*(45*beta*pow(delta, 4) - 140*beta*pow(delta, 3) + 162*beta*pow(delta,
      2) - 84*beta*delta + 17*beta) + pow(gamma, 4)*(18*beta*pow(delta, 2) -
      30*beta*delta + 13*beta)) - 2*x*(pow(gamma, 8)*(4*beta*pow(delta, 7) -
      21*beta*pow(delta, 6) + 45*beta*pow(delta, 5) - 50*beta*pow(delta, 4) +
      30*beta*pow(delta, 3) - 9*beta*pow(delta, 2) + beta*delta) + pow(gamma,
      6)*(9*beta*pow(delta, 5) - 35*beta*pow(delta, 4) + 54*beta*pow(delta, 3)
      - 42*beta*pow(delta, 2) + 17*beta*delta - 3*beta) + pow(gamma,
      4)*(6*beta*pow(delta, 3) - 15*beta*pow(delta, 2) + 13*beta*delta -
      4*beta) + pow(gamma, 2)*(beta*delta - beta)))*exp(beta*delta*atan(gamma*x
      - gamma*(delta - 1)) - beta*x*atan(delta*gamma -
      gamma*x)))*exp((1.0L/2.0L)*beta*log(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1))/gamma)*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma,
      2)*x + pow(gamma, 2)*pow(x, 2) + 1), 2) + (beta*pow(gamma, 4)*pow(x, 4) +
      beta - 2*pow(gamma, 4)*pow(x, 3)*(2*beta*delta - beta) + pow(gamma,
      4)*(beta*pow(delta, 4) - 2*beta*pow(delta, 3) + beta*pow(delta, 2)) +
      pow(gamma, 2)*(2*beta*pow(delta, 2) - 2*beta*delta + beta) + pow(x,
      2)*(2*beta*pow(gamma, 2) + pow(gamma, 4)*(6*beta*pow(delta, 2) -
      6*beta*delta + beta)) - 2*x*(pow(gamma, 4)*(2*beta*pow(delta, 3) -
      3*beta*pow(delta, 2) + beta*delta) + pow(gamma, 2)*(2*beta*delta -
      beta)))*exp(beta*delta*atan(gamma*x - gamma*(delta - 1)) -
      beta*x*atan(delta*gamma - gamma*x) + (1.0L/2.0L)*beta*log(fabs(pow(gamma,
      2)*pow(x, 2) - 2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta,
      2) - 2*delta + 1) + 1))/gamma)*log(fabs(pow(delta, 2)*pow(gamma, 2) -
      2*delta*pow(gamma, 2)*x + pow(gamma, 2)*pow(x, 2) +
      1))*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x +
      pow(gamma, 2)*pow(x, 2) + 1), 2)*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2) - (beta*pow(gamma, 4)*pow(x, 4) + beta - 2*pow(gamma,
      4)*pow(x, 3)*(2*beta*delta - beta) + pow(gamma, 4)*(beta*pow(delta, 4) -
      2*beta*pow(delta, 3) + beta*pow(delta, 2)) + pow(gamma,
      2)*(2*beta*pow(delta, 2) - 2*beta*delta + beta) + pow(x,
      2)*(2*beta*pow(gamma, 2) + pow(gamma, 4)*(6*beta*pow(delta, 2) -
      6*beta*delta + beta)) - 2*x*(pow(gamma, 4)*(2*beta*pow(delta, 3) -
      3*beta*pow(delta, 2) + beta*delta) + pow(gamma, 2)*(2*beta*delta -
      beta)))*exp(beta*delta*atan(gamma*x - gamma*(delta - 1)) -
      beta*x*atan(delta*gamma - gamma*x) + (1.0L/2.0L)*beta*log(fabs(pow(gamma,
      2)*pow(x, 2) - 2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta,
      2) - 2*delta + 1) + 1))/gamma)*log(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1))*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x +
      pow(gamma, 2)*pow(x, 2) + 1), 2)*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2) - 2*(beta*pow(delta, 2)*pow(gamma, 2) + beta*pow(gamma,
      8)*pow(x, 8) - 2*pow(gamma, 8)*pow(x, 7)*(4*beta*delta - beta) +
      pow(gamma, 8)*(beta*pow(delta, 8) - 2*beta*pow(delta, 7) +
      beta*pow(delta, 6)) + pow(gamma, 6)*(3*beta*pow(delta, 6) -
      4*beta*pow(delta, 5) + 2*beta*pow(delta, 4)) + pow(gamma,
      4)*(3*beta*pow(delta, 4) - 2*beta*pow(delta, 3) + beta*pow(delta, 2)) +
      pow(x, 6)*(3*beta*pow(gamma, 6) + pow(gamma, 8)*(28*beta*pow(delta, 2) -
      14*beta*delta + beta)) - 2*pow(x, 5)*(pow(gamma, 8)*(28*beta*pow(delta,
      3) - 21*beta*pow(delta, 2) + 3*beta*delta) + pow(gamma, 6)*(9*beta*delta
      - 2*beta)) + pow(x, 4)*(3*beta*pow(gamma, 4) + 5*pow(gamma,
      8)*(14*beta*pow(delta, 4) - 14*beta*pow(delta, 3) + 3*beta*pow(delta, 2))
      + pow(gamma, 6)*(45*beta*pow(delta, 2) - 20*beta*delta + 2*beta)) -
      2*pow(x, 3)*(pow(gamma, 8)*(28*beta*pow(delta, 5) - 35*beta*pow(delta, 4)
      + 10*beta*pow(delta, 3)) + 2*pow(gamma, 6)*(15*beta*pow(delta, 3) -
      10*beta*pow(delta, 2) + 2*beta*delta) + pow(gamma, 4)*(6*beta*delta -
      beta)) + pow(x, 2)*(beta*pow(gamma, 2) + pow(gamma,
      8)*(28*beta*pow(delta, 6) - 42*beta*pow(delta, 5) + 15*beta*pow(delta,
      4)) + pow(gamma, 6)*(45*beta*pow(delta, 4) - 40*beta*pow(delta, 3) +
      12*beta*pow(delta, 2)) + pow(gamma, 4)*(18*beta*pow(delta, 2) -
      6*beta*delta + beta)) - 2*x*(beta*delta*pow(gamma, 2) + pow(gamma,
      8)*(4*beta*pow(delta, 7) - 7*beta*pow(delta, 6) + 3*beta*pow(delta, 5)) +
      pow(gamma, 6)*(9*beta*pow(delta, 5) - 10*beta*pow(delta, 4) +
      4*beta*pow(delta, 3)) + pow(gamma, 4)*(6*beta*pow(delta, 3) -
      3*beta*pow(delta, 2) + beta*delta)))*exp(beta*delta*atan(gamma*x -
      gamma*(delta - 1)) - beta*x*atan(delta*gamma - gamma*x) +
      (1.0L/2.0L)*beta*log(fabs(pow(gamma, 2)*pow(x, 2) - 2*pow(gamma,
      2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta + 1) +
      1))/gamma)*pow(fabs(pow(gamma, 2)*pow(x, 2) - 2*pow(gamma, 2)*x*(delta -
      1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta + 1) + 1), 2))/((-pow(gamma,
      6)*pow(x, 4) + 2*pow(gamma, 6)*pow(x, 3)*(2*delta - 1) - pow(gamma,
      6)*(pow(delta, 4) - 2*pow(delta, 3) + pow(delta, 2)) - pow(gamma,
      4)*(2*pow(delta, 2) - 2*delta + 1) - pow(gamma, 2) - pow(x,
      2)*(pow(gamma, 6)*(6*pow(delta, 2) - 6*delta + 1) + 2*pow(gamma, 4)) +
      2*x*(pow(gamma, 6)*(2*pow(delta, 3) - 3*pow(delta, 2) + delta) +
      pow(gamma, 4)*(2*delta - 1)))*exp(beta*delta*atan(gamma*x - gamma*(delta
      - 1)) - beta*x*atan(delta*gamma - gamma*x) +
      (1.0L/2.0L)*beta*log(fabs(pow(gamma, 2)*pow(x, 2) - 2*pow(gamma,
      2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta + 1) +
      1))/gamma)*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x
      + pow(gamma, 2)*pow(x, 2) + 1), 2)*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2) + (pow(gamma, 6)*pow(x, 4)*exp(alpha) - 2*pow(gamma,
      6)*pow(x, 3)*(2*delta*exp(alpha) - exp(alpha)) + pow(gamma,
      6)*(pow(delta, 4)*exp(alpha) - 2*pow(delta, 3)*exp(alpha) + pow(delta,
      2)*exp(alpha)) + pow(gamma, 4)*(2*pow(delta, 2)*exp(alpha) -
      2*delta*exp(alpha) + exp(alpha)) + pow(gamma, 2)*exp(alpha) + pow(x,
      2)*(pow(gamma, 6)*(6*pow(delta, 2)*exp(alpha) - 6*delta*exp(alpha) +
      exp(alpha)) + 2*pow(gamma, 4)*exp(alpha)) - 2*x*(pow(gamma,
      6)*(2*pow(delta, 3)*exp(alpha) - 3*pow(delta, 2)*exp(alpha) +
      delta*exp(alpha)) + pow(gamma, 4)*(2*delta*exp(alpha) -
      exp(alpha))))*exp(-beta*delta*atan(delta*gamma - gamma*x) +
      beta*x*atan(gamma*x - gamma*(delta - 1)) + beta*atan(gamma*x -
      gamma*(delta - 1)) + (1.0L/2.0L)*beta*log(fabs(pow(delta, 2)*pow(gamma,
      2) - 2*delta*pow(gamma, 2)*x + pow(gamma, 2)*pow(x, 2) +
      1))/gamma)*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x
      + pow(gamma, 2)*pow(x, 2) + 1), 2)*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2));

      return lb_partial_gamma_partI_result;

}

double lb_partial_gamma_partII(double beta, double delta, double gamma, double x) {

   double lb_partial_gamma_partII_result;
      lb_partial_gamma_partII_result = -1.0L/2.0L*((beta*pow(gamma, 4)*pow(x,
      4) + beta - 2*pow(gamma, 4)*pow(x, 3)*(2*beta*delta - beta) + pow(gamma,
      4)*(beta*pow(delta, 4) - 2*beta*pow(delta, 3) + beta*pow(delta, 2)) +
      pow(gamma, 2)*(2*beta*pow(delta, 2) - 2*beta*delta + beta) + pow(x,
      2)*(2*beta*pow(gamma, 2) + pow(gamma, 4)*(6*beta*pow(delta, 2) -
      6*beta*delta + beta)) - 2*x*(pow(gamma, 4)*(2*beta*pow(delta, 3) -
      3*beta*pow(delta, 2) + beta*delta) + pow(gamma, 2)*(2*beta*delta -
      beta)))*log(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x +
      pow(gamma, 2)*pow(x, 2) + 1))*pow(fabs(pow(delta, 2)*pow(gamma, 2) -
      2*delta*pow(gamma, 2)*x + pow(gamma, 2)*pow(x, 2) + 1),
      2)*pow(fabs(pow(gamma, 2)*pow(x, 2) - 2*pow(gamma, 2)*x*(delta - 1) +
      pow(gamma, 2)*(pow(delta, 2) - 2*delta + 1) + 1), 2) - (beta*pow(gamma,
      4)*pow(x, 4) + beta - 2*pow(gamma, 4)*pow(x, 3)*(2*beta*delta - beta) +
      pow(gamma, 4)*(beta*pow(delta, 4) - 2*beta*pow(delta, 3) +
      beta*pow(delta, 2)) + pow(gamma, 2)*(2*beta*pow(delta, 2) - 2*beta*delta
      + beta) + pow(x, 2)*(2*beta*pow(gamma, 2) + pow(gamma,
      4)*(6*beta*pow(delta, 2) - 6*beta*delta + beta)) - 2*x*(pow(gamma,
      4)*(2*beta*pow(delta, 3) - 3*beta*pow(delta, 2) + beta*delta) +
      pow(gamma, 2)*(2*beta*delta - beta)))*log(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1))*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x +
      pow(gamma, 2)*pow(x, 2) + 1), 2)*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2) - 2*(beta*pow(delta, 2)*pow(gamma, 2) + beta*pow(gamma,
      8)*pow(x, 8) - 2*pow(gamma, 8)*pow(x, 7)*(4*beta*delta - beta) +
      pow(gamma, 8)*(beta*pow(delta, 8) - 2*beta*pow(delta, 7) +
      beta*pow(delta, 6)) + pow(gamma, 6)*(3*beta*pow(delta, 6) -
      4*beta*pow(delta, 5) + 2*beta*pow(delta, 4)) + pow(gamma,
      4)*(3*beta*pow(delta, 4) - 2*beta*pow(delta, 3) + beta*pow(delta, 2)) +
      pow(x, 6)*(3*beta*pow(gamma, 6) + pow(gamma, 8)*(28*beta*pow(delta, 2) -
      14*beta*delta + beta)) - 2*pow(x, 5)*(pow(gamma, 8)*(28*beta*pow(delta,
      3) - 21*beta*pow(delta, 2) + 3*beta*delta) + pow(gamma, 6)*(9*beta*delta
      - 2*beta)) + pow(x, 4)*(3*beta*pow(gamma, 4) + 5*pow(gamma,
      8)*(14*beta*pow(delta, 4) - 14*beta*pow(delta, 3) + 3*beta*pow(delta, 2))
      + pow(gamma, 6)*(45*beta*pow(delta, 2) - 20*beta*delta + 2*beta)) -
      2*pow(x, 3)*(pow(gamma, 8)*(28*beta*pow(delta, 5) - 35*beta*pow(delta, 4)
      + 10*beta*pow(delta, 3)) + 2*pow(gamma, 6)*(15*beta*pow(delta, 3) -
      10*beta*pow(delta, 2) + 2*beta*delta) + pow(gamma, 4)*(6*beta*delta -
      beta)) + pow(x, 2)*(beta*pow(gamma, 2) + pow(gamma,
      8)*(28*beta*pow(delta, 6) - 42*beta*pow(delta, 5) + 15*beta*pow(delta,
      4)) + pow(gamma, 6)*(45*beta*pow(delta, 4) - 40*beta*pow(delta, 3) +
      12*beta*pow(delta, 2)) + pow(gamma, 4)*(18*beta*pow(delta, 2) -
      6*beta*delta + beta)) - 2*x*(beta*delta*pow(gamma, 2) + pow(gamma,
      8)*(4*beta*pow(delta, 7) - 7*beta*pow(delta, 6) + 3*beta*pow(delta, 5)) +
      pow(gamma, 6)*(9*beta*pow(delta, 5) - 10*beta*pow(delta, 4) +
      4*beta*pow(delta, 3)) + pow(gamma, 4)*(6*beta*pow(delta, 3) -
      3*beta*pow(delta, 2) + beta*delta)))*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2) + 2*(beta*pow(gamma, 8)*pow(x, 8) - 2*pow(gamma, 8)*pow(x,
      7)*(4*beta*delta - 3*beta) + pow(gamma, 8)*(beta*pow(delta, 8) -
      6*beta*pow(delta, 7) + 15*beta*pow(delta, 6) - 20*beta*pow(delta, 5) +
      15*beta*pow(delta, 4) - 6*beta*pow(delta, 3) + beta*pow(delta, 2)) +
      pow(gamma, 6)*(3*beta*pow(delta, 6) - 14*beta*pow(delta, 5) +
      27*beta*pow(delta, 4) - 28*beta*pow(delta, 3) + 17*beta*pow(delta, 2) -
      6*beta*delta + beta) + pow(gamma, 4)*(3*beta*pow(delta, 4) -
      10*beta*pow(delta, 3) + 13*beta*pow(delta, 2) - 8*beta*delta + 2*beta) +
      pow(gamma, 2)*(beta*pow(delta, 2) - 2*beta*delta + beta) + pow(x,
      6)*(3*beta*pow(gamma, 6) + pow(gamma, 8)*(28*beta*pow(delta, 2) -
      42*beta*delta + 15*beta)) - 2*pow(x, 5)*(pow(gamma,
      8)*(28*beta*pow(delta, 3) - 63*beta*pow(delta, 2) + 45*beta*delta -
      10*beta) + pow(gamma, 6)*(9*beta*delta - 7*beta)) + pow(x,
      4)*(3*beta*pow(gamma, 4) + 5*pow(gamma, 8)*(14*beta*pow(delta, 4) -
      42*beta*pow(delta, 3) + 45*beta*pow(delta, 2) - 20*beta*delta + 3*beta) +
      pow(gamma, 6)*(45*beta*pow(delta, 2) - 70*beta*delta + 27*beta)) -
      2*pow(x, 3)*(pow(gamma, 8)*(28*beta*pow(delta, 5) - 105*beta*pow(delta,
      4) + 150*beta*pow(delta, 3) - 100*beta*pow(delta, 2) + 30*beta*delta -
      3*beta) + 2*pow(gamma, 6)*(15*beta*pow(delta, 3) - 35*beta*pow(delta, 2)
      + 27*beta*delta - 7*beta) + pow(gamma, 4)*(6*beta*delta - 5*beta)) +
      pow(x, 2)*(beta*pow(gamma, 2) + pow(gamma, 8)*(28*beta*pow(delta, 6) -
      126*beta*pow(delta, 5) + 225*beta*pow(delta, 4) - 200*beta*pow(delta, 3)
      + 90*beta*pow(delta, 2) - 18*beta*delta + beta) + pow(gamma,
      6)*(45*beta*pow(delta, 4) - 140*beta*pow(delta, 3) + 162*beta*pow(delta,
      2) - 84*beta*delta + 17*beta) + pow(gamma, 4)*(18*beta*pow(delta, 2) -
      30*beta*delta + 13*beta)) - 2*x*(pow(gamma, 8)*(4*beta*pow(delta, 7) -
      21*beta*pow(delta, 6) + 45*beta*pow(delta, 5) - 50*beta*pow(delta, 4) +
      30*beta*pow(delta, 3) - 9*beta*pow(delta, 2) + beta*delta) + pow(gamma,
      6)*(9*beta*pow(delta, 5) - 35*beta*pow(delta, 4) + 54*beta*pow(delta, 3)
      - 42*beta*pow(delta, 2) + 17*beta*delta - 3*beta) + pow(gamma,
      4)*(6*beta*pow(delta, 3) - 15*beta*pow(delta, 2) + 13*beta*delta -
      4*beta) + pow(gamma, 2)*(beta*delta - beta)) - (2*beta*pow(gamma, 2)*x -
      pow(gamma, 2)*(2*beta*delta - beta))*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2))*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma,
      2)*x + pow(gamma, 2)*pow(x, 2) + 1), 2))/((-pow(gamma, 6)*pow(x, 4) +
      2*pow(gamma, 6)*pow(x, 3)*(2*delta - 1) - pow(gamma, 6)*(pow(delta, 4) -
      2*pow(delta, 3) + pow(delta, 2)) - pow(gamma, 4)*(2*pow(delta, 2) -
      2*delta + 1) - pow(gamma, 2) - pow(x, 2)*(pow(gamma, 6)*(6*pow(delta, 2)
      - 6*delta + 1) + 2*pow(gamma, 4)) + 2*x*(pow(gamma, 6)*(2*pow(delta, 3) -
      3*pow(delta, 2) + delta) + pow(gamma, 4)*(2*delta -
      1)))*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x +
      pow(gamma, 2)*pow(x, 2) + 1), 2)*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2));

       return lb_partial_gamma_partII_result;

}

double lb_partial_delta_partI(double alpha, double beta, double delta, double gamma, double x) {

   double lb_partial_delta_partI_result;
      lb_partial_delta_partI_result = -(((-(beta*pow(gamma, 3)*pow(x, 2) -
      beta*gamma - pow(gamma, 3)*x*(2*beta*delta - beta) + pow(gamma,
      3)*(beta*pow(delta, 2) - beta*delta))*exp(beta*delta*atan(gamma*x -
      gamma*(delta - 1)) - beta*x*atan(delta*gamma - gamma*x)) +
      (beta*pow(gamma, 4)*pow(x, 4) + beta - 2*pow(gamma, 4)*pow(x,
      3)*(2*beta*delta - beta) + pow(gamma, 4)*(beta*pow(delta, 4) -
      2*beta*pow(delta, 3) + beta*pow(delta, 2)) + pow(gamma,
      2)*(2*beta*pow(delta, 2) - 2*beta*delta + beta) + pow(x,
      2)*(2*beta*pow(gamma, 2) + pow(gamma, 4)*(6*beta*pow(delta, 2) -
      6*beta*delta + beta)) - 2*x*(pow(gamma, 4)*(2*beta*pow(delta, 3) -
      3*beta*pow(delta, 2) + beta*delta) + pow(gamma, 2)*(2*beta*delta -
      beta)))*exp(beta*delta*atan(gamma*x - gamma*(delta - 1)) -
      beta*x*atan(delta*gamma - gamma*x))*atan(delta*gamma - gamma*x) +
      (beta*pow(gamma, 4)*pow(x, 4) + beta - 2*pow(gamma, 4)*pow(x,
      3)*(2*beta*delta - beta) + pow(gamma, 4)*(beta*pow(delta, 4) -
      2*beta*pow(delta, 3) + beta*pow(delta, 2)) + pow(gamma,
      2)*(2*beta*pow(delta, 2) - 2*beta*delta + beta) + pow(x,
      2)*(2*beta*pow(gamma, 2) + pow(gamma, 4)*(6*beta*pow(delta, 2) -
      6*beta*delta + beta)) - 2*x*(pow(gamma, 4)*(2*beta*pow(delta, 3) -
      3*beta*pow(delta, 2) + beta*delta) + pow(gamma, 2)*(2*beta*delta -
      beta)))*exp(beta*delta*atan(gamma*x - gamma*(delta - 1)) -
      beta*x*atan(delta*gamma - gamma*x))*atan(gamma*x - gamma*(delta -
      1)))*pow(fabs(pow(gamma, 2)*pow(x, 2) - 2*pow(gamma, 2)*x*(delta - 1) +
      pow(gamma, 2)*(pow(delta, 2) - 2*delta + 1) + 1), 2) - (beta*pow(gamma,
      7)*pow(x, 7) - pow(gamma, 7)*pow(x, 6)*(7*beta*delta - 5*beta) -
      pow(gamma, 7)*(beta*pow(delta, 7) - 5*beta*pow(delta, 6) +
      10*beta*pow(delta, 5) - 10*beta*pow(delta, 4) + 5*beta*pow(delta, 3) -
      beta*pow(delta, 2)) - pow(gamma, 5)*(3*beta*pow(delta, 5) -
      11*beta*pow(delta, 4) + 16*beta*pow(delta, 3) - 12*beta*pow(delta, 2) +
      5*beta*delta - beta) - pow(gamma, 3)*(3*beta*pow(delta, 3) -
      7*beta*pow(delta, 2) + 6*beta*delta - 2*beta) - gamma*(beta*delta - beta)
      + pow(x, 5)*(3*beta*pow(gamma, 5) + pow(gamma, 7)*(21*beta*pow(delta, 2)
      - 30*beta*delta + 10*beta)) - pow(x, 4)*(5*pow(gamma,
      7)*(7*beta*pow(delta, 3) - 15*beta*pow(delta, 2) + 10*beta*delta -
      2*beta) + pow(gamma, 5)*(15*beta*delta - 11*beta)) + pow(x,
      3)*(3*beta*pow(gamma, 3) + 5*pow(gamma, 7)*(7*beta*pow(delta, 4) -
      20*beta*pow(delta, 3) + 20*beta*pow(delta, 2) - 8*beta*delta + beta) +
      2*pow(gamma, 5)*(15*beta*pow(delta, 2) - 22*beta*delta + 8*beta)) -
      pow(x, 2)*(pow(gamma, 7)*(21*beta*pow(delta, 5) - 75*beta*pow(delta, 4) +
      100*beta*pow(delta, 3) - 60*beta*pow(delta, 2) + 15*beta*delta - beta) +
      6*pow(gamma, 5)*(5*beta*pow(delta, 3) - 11*beta*pow(delta, 2) +
      8*beta*delta - 2*beta) + pow(gamma, 3)*(9*beta*delta - 7*beta)) +
      x*(beta*gamma + pow(gamma, 7)*(7*beta*pow(delta, 6) - 30*beta*pow(delta,
      5) + 50*beta*pow(delta, 4) - 40*beta*pow(delta, 3) + 15*beta*pow(delta,
      2) - 2*beta*delta) + pow(gamma, 5)*(15*beta*pow(delta, 4) -
      44*beta*pow(delta, 3) + 48*beta*pow(delta, 2) - 24*beta*delta + 5*beta) +
      pow(gamma, 3)*(9*beta*pow(delta, 2) - 14*beta*delta +
      6*beta)))*exp(beta*delta*atan(gamma*x - gamma*(delta - 1)) -
      beta*x*atan(delta*gamma -
      gamma*x)))*exp((1.0L/2.0L)*beta*log(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1))/gamma)*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma,
      2)*x + pow(gamma, 2)*pow(x, 2) + 1), 2) + (-beta*delta*gamma +
      beta*pow(gamma, 7)*pow(x, 7) - pow(gamma, 7)*pow(x, 6)*(7*beta*delta -
      2*beta) - pow(gamma, 7)*(beta*pow(delta, 7) - 2*beta*pow(delta, 6) +
      beta*pow(delta, 5)) - pow(gamma, 5)*(3*beta*pow(delta, 5) -
      4*beta*pow(delta, 4) + 2*beta*pow(delta, 3)) - pow(gamma,
      3)*(3*beta*pow(delta, 3) - 2*beta*pow(delta, 2) + beta*delta) + pow(x,
      5)*(3*beta*pow(gamma, 5) + pow(gamma, 7)*(21*beta*pow(delta, 2) -
      12*beta*delta + beta)) - pow(x, 4)*(5*pow(gamma, 7)*(7*beta*pow(delta, 3)
      - 6*beta*pow(delta, 2) + beta*delta) + pow(gamma, 5)*(15*beta*delta -
      4*beta)) + pow(x, 3)*(3*beta*pow(gamma, 3) + 5*pow(gamma,
      7)*(7*beta*pow(delta, 4) - 8*beta*pow(delta, 3) + 2*beta*pow(delta, 2)) +
      2*pow(gamma, 5)*(15*beta*pow(delta, 2) - 8*beta*delta + beta)) - pow(x,
      2)*(pow(gamma, 7)*(21*beta*pow(delta, 5) - 30*beta*pow(delta, 4) +
      10*beta*pow(delta, 3)) + 6*pow(gamma, 5)*(5*beta*pow(delta, 3) -
      4*beta*pow(delta, 2) + beta*delta) + pow(gamma, 3)*(9*beta*delta -
      2*beta)) + x*(beta*gamma + pow(gamma, 7)*(7*beta*pow(delta, 6) -
      12*beta*pow(delta, 5) + 5*beta*pow(delta, 4)) + pow(gamma,
      5)*(15*beta*pow(delta, 4) - 16*beta*pow(delta, 3) + 6*beta*pow(delta, 2))
      + pow(gamma, 3)*(9*beta*pow(delta, 2) - 4*beta*delta +
      beta)))*exp(beta*delta*atan(gamma*x - gamma*(delta - 1)) -
      beta*x*atan(delta*gamma - gamma*x) + (1.0L/2.0L)*beta*log(fabs(pow(gamma,
      2)*pow(x, 2) - 2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta,
      2) - 2*delta + 1) + 1))/gamma)*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2))/((-pow(gamma, 4)*pow(x, 4) + 2*pow(gamma, 4)*pow(x,
      3)*(2*delta - 1) - pow(gamma, 4)*(pow(delta, 4) - 2*pow(delta, 3) +
      pow(delta, 2)) - pow(gamma, 2)*(2*pow(delta, 2) - 2*delta + 1) - pow(x,
      2)*(pow(gamma, 4)*(6*pow(delta, 2) - 6*delta + 1) + 2*pow(gamma, 2)) +
      2*x*(pow(gamma, 4)*(2*pow(delta, 3) - 3*pow(delta, 2) + delta) +
      pow(gamma, 2)*(2*delta - 1)) - 1)*exp(beta*delta*atan(gamma*x -
      gamma*(delta - 1)) - beta*x*atan(delta*gamma - gamma*x) +
      (1.0L/2.0L)*beta*log(fabs(pow(gamma, 2)*pow(x, 2) - 2*pow(gamma,
      2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta + 1) +
      1))/gamma)*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x
      + pow(gamma, 2)*pow(x, 2) + 1), 2)*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2) + (pow(gamma, 4)*pow(x, 4)*exp(alpha) - 2*pow(gamma,
      4)*pow(x, 3)*(2*delta*exp(alpha) - exp(alpha)) + pow(gamma,
      4)*(pow(delta, 4)*exp(alpha) - 2*pow(delta, 3)*exp(alpha) + pow(delta,
      2)*exp(alpha)) + pow(gamma, 2)*(2*pow(delta, 2)*exp(alpha) -
      2*delta*exp(alpha) + exp(alpha)) + pow(x, 2)*(pow(gamma, 4)*(6*pow(delta,
      2)*exp(alpha) - 6*delta*exp(alpha) + exp(alpha)) + 2*pow(gamma,
      2)*exp(alpha)) - 2*x*(pow(gamma, 4)*(2*pow(delta, 3)*exp(alpha) -
      3*pow(delta, 2)*exp(alpha) + delta*exp(alpha)) + pow(gamma,
      2)*(2*delta*exp(alpha) - exp(alpha))) +
      exp(alpha))*exp(-beta*delta*atan(delta*gamma - gamma*x) +
      beta*x*atan(gamma*x - gamma*(delta - 1)) + beta*atan(gamma*x -
      gamma*(delta - 1)) + (1.0L/2.0L)*beta*log(fabs(pow(delta, 2)*pow(gamma,
      2) - 2*delta*pow(gamma, 2)*x + pow(gamma, 2)*pow(x, 2) +
      1))/gamma)*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x
      + pow(gamma, 2)*pow(x, 2) + 1), 2)*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2));

      return lb_partial_delta_partI_result;

}

double lb_partial_delta_partII(double beta, double delta, double gamma, double x) {

   double lb_partial_delta_partII_result;
      lb_partial_delta_partII_result = (-(-beta*delta*gamma + beta*pow(gamma,
      7)*pow(x, 7) - pow(gamma, 7)*pow(x, 6)*(7*beta*delta - 2*beta) -
      pow(gamma, 7)*(beta*pow(delta, 7) - 2*beta*pow(delta, 6) +
      beta*pow(delta, 5)) - pow(gamma, 5)*(3*beta*pow(delta, 5) -
      4*beta*pow(delta, 4) + 2*beta*pow(delta, 3)) - pow(gamma,
      3)*(3*beta*pow(delta, 3) - 2*beta*pow(delta, 2) + beta*delta) + pow(x,
      5)*(3*beta*pow(gamma, 5) + pow(gamma, 7)*(21*beta*pow(delta, 2) -
      12*beta*delta + beta)) - pow(x, 4)*(5*pow(gamma, 7)*(7*beta*pow(delta, 3)
      - 6*beta*pow(delta, 2) + beta*delta) + pow(gamma, 5)*(15*beta*delta -
      4*beta)) + pow(x, 3)*(3*beta*pow(gamma, 3) + 5*pow(gamma,
      7)*(7*beta*pow(delta, 4) - 8*beta*pow(delta, 3) + 2*beta*pow(delta, 2)) +
      2*pow(gamma, 5)*(15*beta*pow(delta, 2) - 8*beta*delta + beta)) - pow(x,
      2)*(pow(gamma, 7)*(21*beta*pow(delta, 5) - 30*beta*pow(delta, 4) +
      10*beta*pow(delta, 3)) + 6*pow(gamma, 5)*(5*beta*pow(delta, 3) -
      4*beta*pow(delta, 2) + beta*delta) + pow(gamma, 3)*(9*beta*delta -
      2*beta)) + x*(beta*gamma + pow(gamma, 7)*(7*beta*pow(delta, 6) -
      12*beta*pow(delta, 5) + 5*beta*pow(delta, 4)) + pow(gamma,
      5)*(15*beta*pow(delta, 4) - 16*beta*pow(delta, 3) + 6*beta*pow(delta, 2))
      + pow(gamma, 3)*(9*beta*pow(delta, 2) - 4*beta*delta +
      beta)))*pow(fabs(pow(gamma, 2)*pow(x, 2) - 2*pow(gamma, 2)*x*(delta - 1)
      + pow(gamma, 2)*(pow(delta, 2) - 2*delta + 1) + 1), 2) + (beta*pow(gamma,
      7)*pow(x, 7) - pow(gamma, 7)*pow(x, 6)*(7*beta*delta - 5*beta) -
      pow(gamma, 7)*(beta*pow(delta, 7) - 5*beta*pow(delta, 6) +
      10*beta*pow(delta, 5) - 10*beta*pow(delta, 4) + 5*beta*pow(delta, 3) -
      beta*pow(delta, 2)) - pow(gamma, 5)*(3*beta*pow(delta, 5) -
      11*beta*pow(delta, 4) + 16*beta*pow(delta, 3) - 12*beta*pow(delta, 2) +
      5*beta*delta - beta) - pow(gamma, 3)*(3*beta*pow(delta, 3) -
      7*beta*pow(delta, 2) + 6*beta*delta - 2*beta) - gamma*(beta*delta - beta)
      + pow(x, 5)*(3*beta*pow(gamma, 5) + pow(gamma, 7)*(21*beta*pow(delta, 2)
      - 30*beta*delta + 10*beta)) - pow(x, 4)*(5*pow(gamma,
      7)*(7*beta*pow(delta, 3) - 15*beta*pow(delta, 2) + 10*beta*delta -
      2*beta) + pow(gamma, 5)*(15*beta*delta - 11*beta)) + pow(x,
      3)*(3*beta*pow(gamma, 3) + 5*pow(gamma, 7)*(7*beta*pow(delta, 4) -
      20*beta*pow(delta, 3) + 20*beta*pow(delta, 2) - 8*beta*delta + beta) +
      2*pow(gamma, 5)*(15*beta*pow(delta, 2) - 22*beta*delta + 8*beta)) -
      pow(x, 2)*(pow(gamma, 7)*(21*beta*pow(delta, 5) - 75*beta*pow(delta, 4) +
      100*beta*pow(delta, 3) - 60*beta*pow(delta, 2) + 15*beta*delta - beta) +
      6*pow(gamma, 5)*(5*beta*pow(delta, 3) - 11*beta*pow(delta, 2) +
      8*beta*delta - 2*beta) + pow(gamma, 3)*(9*beta*delta - 7*beta)) +
      x*(beta*gamma + pow(gamma, 7)*(7*beta*pow(delta, 6) - 30*beta*pow(delta,
      5) + 50*beta*pow(delta, 4) - 40*beta*pow(delta, 3) + 15*beta*pow(delta,
      2) - 2*beta*delta) + pow(gamma, 5)*(15*beta*pow(delta, 4) -
      44*beta*pow(delta, 3) + 48*beta*pow(delta, 2) - 24*beta*delta + 5*beta) +
      pow(gamma, 3)*(9*beta*pow(delta, 2) - 14*beta*delta + 6*beta)) +
      (beta*pow(gamma, 3)*pow(x, 2) - beta*gamma - pow(gamma,
      3)*x*(2*beta*delta - beta) + pow(gamma, 3)*(beta*pow(delta, 2) -
      beta*delta) - (beta*pow(gamma, 4)*pow(x, 4) + beta - 2*pow(gamma,
      4)*pow(x, 3)*(2*beta*delta - beta) + pow(gamma, 4)*(beta*pow(delta, 4) -
      2*beta*pow(delta, 3) + beta*pow(delta, 2)) + pow(gamma,
      2)*(2*beta*pow(delta, 2) - 2*beta*delta + beta) + pow(x,
      2)*(2*beta*pow(gamma, 2) + pow(gamma, 4)*(6*beta*pow(delta, 2) -
      6*beta*delta + beta)) - 2*x*(pow(gamma, 4)*(2*beta*pow(delta, 3) -
      3*beta*pow(delta, 2) + beta*delta) + pow(gamma, 2)*(2*beta*delta -
      beta)))*atan(delta*gamma - gamma*x) - (beta*pow(gamma, 4)*pow(x, 4) +
      beta - 2*pow(gamma, 4)*pow(x, 3)*(2*beta*delta - beta) + pow(gamma,
      4)*(beta*pow(delta, 4) - 2*beta*pow(delta, 3) + beta*pow(delta, 2)) +
      pow(gamma, 2)*(2*beta*pow(delta, 2) - 2*beta*delta + beta) + pow(x,
      2)*(2*beta*pow(gamma, 2) + pow(gamma, 4)*(6*beta*pow(delta, 2) -
      6*beta*delta + beta)) - 2*x*(pow(gamma, 4)*(2*beta*pow(delta, 3) -
      3*beta*pow(delta, 2) + beta*delta) + pow(gamma, 2)*(2*beta*delta -
      beta)))*atan(gamma*x - gamma*(delta - 1)))*pow(fabs(pow(gamma, 2)*pow(x,
      2) - 2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) -
      2*delta + 1) + 1), 2))*pow(fabs(pow(delta, 2)*pow(gamma, 2) -
      2*delta*pow(gamma, 2)*x + pow(gamma, 2)*pow(x, 2) + 1), 2))/((-pow(gamma,
      4)*pow(x, 4) + 2*pow(gamma, 4)*pow(x, 3)*(2*delta - 1) - pow(gamma,
      4)*(pow(delta, 4) - 2*pow(delta, 3) + pow(delta, 2)) - pow(gamma,
      2)*(2*pow(delta, 2) - 2*delta + 1) - pow(x, 2)*(pow(gamma,
      4)*(6*pow(delta, 2) - 6*delta + 1) + 2*pow(gamma, 2)) + 2*x*(pow(gamma,
      4)*(2*pow(delta, 3) - 3*pow(delta, 2) + delta) + pow(gamma, 2)*(2*delta -
      1)) - 1)*pow(fabs(pow(delta, 2)*pow(gamma, 2) - 2*delta*pow(gamma, 2)*x +
      pow(gamma, 2)*pow(x, 2) + 1), 2)*pow(fabs(pow(gamma, 2)*pow(x, 2) -
      2*pow(gamma, 2)*x*(delta - 1) + pow(gamma, 2)*(pow(delta, 2) - 2*delta +
      1) + 1), 2));

      return lb_partial_delta_partII_result;

}

/*
 * END OF CODE PASTED IN FROM SAGE
 */

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



// [[Rcpp::export]]
NumericVector mortalityhazard_lb_binomial_grad_cpp(NumericVector theta,
                                                         NumericVector ages,
                                                         NumericVector Dx,
                                                         NumericVector Nx)
{

  int theta_len = theta.size();
  int ages_len = ages.size();

  NumericVector res(theta_len);

  double alpha = theta[0];
  double beta = exp(theta[1]);
  double gamma = exp(theta[2]);
  double delta = theta[3];

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
      partI += (double)Dx[ageidx] * lb_partial_alpha_partI(alpha, beta, delta, gamma, cur_age);
      // this is just -1, but we'll keep the fn call to make maintaining this easier
      partII += (double)(Nx[ageidx] - Dx[ageidx]) * lb_partial_alpha_partII();

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
      partI += (double)Dx[ageidx] * lb_partial_beta_partI(alpha, beta, delta, gamma, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * lb_partial_beta_partII(delta,
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
      partI += (double)Dx[ageidx] * lb_partial_gamma_partI(alpha, beta, delta, gamma, cur_age);
      // partII is always -1, but we keep the function call here for clarity and maintainability
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * lb_partial_gamma_partII(beta, delta, gamma, cur_age);

  }

  //Rcout << "partial wrt gamma: " << partI + partII << std::endl;

  res[GAMMA_IDX] = partI + partII;

  /*
   * delta
   */

  partI = 0.0;
  partII = 0.0;

  for(int ageidx=0; ageidx < ages_len; ageidx++) {

      cur_age = ages[ageidx];

      // for each age, we want
      // D_z * partI + (N_z - D_z) * partII
      partI += (double)Dx[ageidx] * lb_partial_delta_partI(alpha, beta, delta, gamma, cur_age);
      partII+= (double)(Nx[ageidx] - Dx[ageidx]) * lb_partial_delta_partII(beta, delta, gamma, cur_age);

  }

  //Rcout << "partial wrt delta: " << partI + partII << std::endl;

  res[DELTA_IDX] = partI + partII;

  return(res);

}

