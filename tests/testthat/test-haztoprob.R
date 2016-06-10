#########################################################
## test the conversion of hazards to probabilities
##
## TODO -- comment the SAGE part of this analysis and
##         include sourcefiles in package...
##

#############################
## logistic hazard

context("testing haz.to.prob - logistic hazard")


## i obtained these from the sage package
## these are for the default parameter values for the
## logistic function
## (which should be one of the trickiest ones...)
sage.val <- c(0.090244779227302008, 0.099046566325919461, 0.10853442803743796,
0.11870174755925278, 0.12952850516478176, 0.14098000032205571,
0.15300609788046815, 0.1655411485189987, 0.17850469844914119,
0.1918030462157555, 0.20533162954612094, 0.21897814117835945,
0.23262619146477859, 0.24615927051417363, 0.25946472524199482,
0.27243746396819768, 0.28498313382985607, 0.29702057838092255,
0.30848346336233612, 0.319321044366077)

test.logistic <- logistic.haz@haz.to.prob.fn(logistic.haz@theta.default,
                                             0:19)

expect_that(test.logistic,equals(sage.val))

#############################
## log-quadratic hazard

context("testing haz.to.prob - log-quadratic hazard")

#######
## i obtained these numerical estimates for probs of death
## from the sage package

## these are for the parameter values for a fit to denma-m-cohort-1841
## theta = {'alpha':-2.8905, 'beta':0.1262, 'gamma':-0.0023}
lq.sage.val <- c(0.0574430728656160,0.0646226401604976,0.0723434382034114,
              0.0805901232194092,0.0893378715517511,0.0985520398947749,
              0.108188073700139,0.118191683497049,0.128499297052722,
              0.139038782324000,0.149730422948784,0.160488115613563,
              0.171220748002960,0.181833708048967,0.192230470480138,
              0.202314205550720,0.211989357294424,0.221163144352396,
              0.229746944738201,0.237657535982441)

test.theta <- c(-2.8905, 0.1262, -0.0023)

test.log_quadratic <- quad.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.log_quadratic,equals(lq.sage.val))

## these are for the parameter values for a fit to denma-m-cohort-1841
## theta={alpha:-3.4649, beta:0.1039, gamma:-5e-4} 
lq.sage.val <- c(0.0324160203916671,0.0358652529680389,0.0396350632396515,
  0.0437492469642166,0.0482325824635806,0.0531107588024771,
  0.0584102830182595,0.0641583642832294,0.0703827729262864,
  0.0771116723419858,0.0843734219838408,0.0921963498928657,
  0.100608493551380,0.109637308298288,0.119309343093989,
  0.129649884094233,0.140682567283930,0.152428962336345,
  0.164908130889175,0.178136163571156)

test.theta <- c(-3.4649, 0.1039, -5e-4)

test.log_quadratic <- quad.haz@haz.to.prob.fn(test.theta, 0:19)

expect_that(test.log_quadratic,equals(lq.sage.val))


## NB: using the approximation exp(-Mx) does not
##     approximate numerical results especially well,
##     according to this test

## TODO -- run for a few scenarios and be sure that
##         it generally converges, and also that the values
##         are probabilities (ie, 0 <= phat <= 1)

## TODO -- take some simple examples with easy closed-form answers
##         and check them (Gompertz, Weibull, etc)

## TODO -- eventually, use Mathematica to get some external
##         computations for some examples; these can be good
##         validators
