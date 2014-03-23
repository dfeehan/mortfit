#########################################################
## test the conversion of hazards to probabilities
##
## TODO -- comment the SAGE part of this analysis and
##         include sourcefiles in package...
##


##### first test the approximation(s)



##### then test the numerical integration version...

##this.theta <- c(-2.5, .08)

##test1 <- haz.to.prob( gomp.haz@haz.fn,
##                      this.theta,
##                      0:19 )

## now compute analytic result and compare...

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

test.logistic <- haz.to.prob(logistic.haz@haz.fn,
                             logistic.haz@theta.default,
                             0:19)

expect_that(test.logistic,equals(sage.val))

##expect_that(1,equals(0))

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
