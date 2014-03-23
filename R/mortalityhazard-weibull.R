#########################################################
## weibull hazard object

## weibull hazard fn
weibull.haz.fn <- function(theta, z) {
    ##alpha <- exp(theta[1])
    alpha <- theta[1]
    beta <-  theta[2]

    res <- alpha*(z^beta)

    ## ensure only non-negative results are returned
    ## (hazards can't be negative)
    res[res<0] <- NA

    return(res)
}

weibull.haz.to.prob <- function(haz.fn, theta, z) {

    alpha <- theta[1]
    beta <-  theta[2]
  
    res <- (alpha*((z+1)^(1+beta) - z^(1+beta)))/(1+beta)
    res <- 1 - exp(-res)

    return(res)
    
}

## weibull hazard fn implemented in c++
weibull.haz.fn.cpp <- function(theta, z) {
  res <- .Call("mortalityhazard_weibull_cpp", theta, z)
  return(res)
}

## weibull hazard fn implemented in c++
weibull.haz.to.prob.cpp <- function(haz.fn, theta, z) {
  ## note that we don't need the first argument, haz.fn
  res <- .Call("mortalityhazard_to_prob_weibull_cpp", theta, z)
  return(res)
}

## (these starting values have been updated based on preliminary analysis)
weibull.haz <- new("mortalityHazard",
                   name="Weibull",
                   num.param=2L,
                   theta.default=c(0.038, 0.487),
                   theta.range=list(c(0.0212, 0.0548),
                                    c(0.386, 0.623)),
                   theta.start.fn=function(data.obj) {
                     ## choose starting values by getting b from
                     ## a logistic regression...
                     dat <- data.obj@data
                     crude <- dat$Dx/(dat$Nx-0.5*dat$Dx)
                     crude[crude==0] <- 1e-5
                     prelim.t <- coef(lm(log(crude) ~ log(age),data=dat))
                     return(c(exp(prelim.t[1]), prelim.t[2]))
                     ##return(c(prelim.t[1], prelim.t[2]))                     
                   },                   
                   ##NB: was Nelder-Mead
                   optim.default=list(method="BFGS",
                                      control=list(reltol=1e-8,
                                                   ##parscale=c(0.03,
                                                   parscale=c(0.038,
                                                              0.487))),
                   haz.fn=weibull.haz.fn.cpp,
                   haz.to.prob.fn=weibull.haz.to.prob.cpp)
                   ##haz.fn=weibull.haz.fn,
                   ##haz.to.prob.fn=weibull.haz.to.prob)
