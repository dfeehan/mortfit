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
                 theta.default=c(-2.93, 7.85e-02, 7.04e-06),
                 theta.range=list(c(-3.35, -2.64),
                                  c(0.065, 9.15e-02),
                                  c(-.022, 1.00e-05)),
                 theta.start.fn=function(data.obj) {
                   ## choose starting values by getting b from
                   ## a logistic regression...
                   dat <- data.obj@data
                   crude <- dat$Dx/(dat$Nx-0.5*dat$Dx)
                   crude[crude==0] <- 1e-5
                   prelim.b <- coef(lm(log(crude) ~ age,data=dat))
                   return(c(prelim.b[1], prelim.b[2], log(1e-5)))
                 },                 
                 optim.default=list(method="BFGS",
                                    control=list(parscale=c(-2.93, 7.85e-02, log(1e-5)),##7.04e-06),
                                                 reltol=1e-10,
                                                 maxit=10000)),
                 haz.fn=mortalityhazard_makeham_cpp,
                 haz.to.prob.fn=mortalityhazard_to_prob_makeham_cpp)
                 #haz.fn=mak.haz.fn.cpp,
                 #haz.to.prob.fn=mak.haz.to.prob.cpp)
                 ##haz.fn=mak.haz.fn,
                 ##haz.to.prob.fn=mak.haz.to.prob)
