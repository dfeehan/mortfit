########################
# makeham hazard 

## makeham hazard fn
mak.haz.fn <- function(theta, z) {

  alpha <- exp(theta[1])
  beta <- theta[2]
  gamma <- exp(theta[3])

  res <- gamma + alpha*exp(beta*z)

  ## ensure only non-negative results are returned
  ## (hazards can't be negative)
  res[res<0] <- NA
  
  return(res)

}

mak.haz.to.prob <- function(haz.fn, theta, z) {

  alpha <- exp(theta[1])
  beta <- theta[2]
  gamma <- exp(theta[3])

  res <- (alpha/beta)*(exp(beta*(z+1))-exp(beta*z)) + gamma
  res <- 1 - exp(-res)

  return(res)

}

## mak hazard fn implemented in c++
mak.haz.fn.cpp <- function(theta, z) {
  res <- .Call("mortalityhazard_makeham_cpp", PACKAGE='mortfit', theta, z)
  return(res)
}

## mak hazard fn implemented in c++
mak.haz.to.prob.cpp <- function(haz.fn, theta, z) {
  ## note that we don't need the first argument, haz.fn
  res <- .Call("mortalityhazard_to_prob_makeham_cpp", PACKAGE='mortfit', theta, z)
  return(res)
}

mak.haz   <- new("mortalityHazard",
                 name="Makeham",
                 num.param=3L,
                 theta.default=c(0.05, 0.1, 1e-3),
                 theta.range=list(c(.01, .1),
                                  c(0.04, .15),
                                  c(1e-6, 0.1)),
                 theta.start.fn=function(data.obj) {

                      ## choose starting values by getting b from
                      ## a regression on the log approximate rates
                      dat <- data.obj@data
                      crude <- dat$Dx/(dat$Nx-0.5*dat$Dx)
                      crude[crude<1e-10] <- 1e-10

                      #offset <- min(crude) 
                      offset <- 1e-10
                      logcrude.shifted <- log(crude) - offset

                      prelim.b <- coef(lm(logcrude.shifted ~ age,data=dat))
                      return(c(exp(prelim.b[1]), prelim.b[2], offset))

                 },                 
                 optim.default=list(method="BFGS",
                                    control=list(parscale=c(0.01, 0.01, 0.0001),
                                                 reltol=1e-12,
                                                 maxit=10000)),
                 haz.fn=mortalityhazard_makeham_cpp,
                 #binomial.grad.fn=NULL,
                 binomial.grad.fn=mortalityhazard_makeham_binomial_grad_cpp,
                 haz.to.prob.fn=mortalityhazard_to_prob_makeham_cpp)
