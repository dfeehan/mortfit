###############################
# gompertz hazard function

### hazard
gomp.haz.fn <- function(theta, z) {

    alpha <- exp(theta[1])
    beta <- theta[2]

    return(alpha*exp(beta*z))
}

gomp.haz.to.prob <- function(haz.fn, theta, z) {

    alpha <- exp(theta[1])
    beta <- theta[2]
  
    res <- (alpha/beta)*(exp(beta*(z+1))-exp(beta*z))
    res <- 1 - exp(-res)

    return(res)
}

## gomp hazard fn implemented in c++
gomp.haz.fn.cpp <- function(theta, z) {
  res <- .Call("mortalityhazard_gompertz_cpp", PACKAGE='mortfit', theta, z)
  return(res)
}

## gomp hazard fn implemented in c++
gomp.haz.to.prob.cpp <- function(haz.fn, theta, z) {
  ## note that we don't need the first argument, haz.fn
  res <- .Call("mortalityhazard_to_prob_gompertz_cpp", PACKAGE='mortfit', theta, z)
  return(res)
}

## these starting values have been updated based on preliminary analysis
gomp.haz <- new("mortalityHazard",
                name="Gompertz",
                num.param=2L,
                theta.default=c(-3, 0.08),
                theta.range=list(c(-3.42, -2.71),
                                 c(0.07, 0.09)),
                theta.start.fn=function(data.obj) {
                  ## choose starting values by getting b from
                  ## a regression on the log approximate rates
                  dat <- data.obj@data
                  crude <- dat$Dx/(dat$Nx-0.5*dat$Dx)
                  crude[crude<1e-10] <- 1e-10
                  regout <- coef(lm(log(crude)~age,data=dat))
                  return(c(regout[1], regout[2]))
                },
                optim.default=list(method="BFGS",
                                   control=list(parscale=c(0.1, 0.1),
                                                reltol=1e-12)),
                #binomial.grad.fn=NULL,
                binomial.grad.fn=mortalityhazard_gompertz_binomial_grad_cpp,
                haz.fn=mortalityhazard_gompertz_cpp,
                haz.to.prob.fn=mortalityhazard_to_prob_gompertz_cpp)
                #haz.fn=gomp.haz.fn.cpp,
                #haz.to.prob.fn=gomp.haz.to.prob.cpp)
                ##haz.fn=gomp.haz.fn,
                ##haz.to.prob.fn=gomp.haz.to.prob)
