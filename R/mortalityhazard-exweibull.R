########################
# extended weibull hazard object

## extended weibull hazard fn
exweibull.haz.fn <- function(theta, z) {
    alpha <- theta[1]
    beta <-  theta[2]
    gamma <- theta[3]

    res <- alpha*(z-gamma)^beta

    ## ensure only non-negative results are returned
    ## (hazards can't be negative)
    res[res<0] <- NA

    ## also ensure that non numeric results aren't returned
    res[is.nan(res)] <- NA

    return(res)
}

exweibull.haz <- new("mortalityHazard",
                     name="Extended Weibull",
                     num.param=2L,
                     theta.default=c(2.78e-04, 1.879, -14.82),
                     theta.range=list(c(9.22e-05, 0.0023),
                                      c(1.36, 2.13),
                                      c(-1.799e+01, -9.834)),
                     theta.start.fn=function(data.obj) {
                       ## choose starting values by getting b from
                       ## a regression...
                       dat <- data.obj@data
                       crude <- dat$Dx/(dat$Nx-0.5*dat$Dx)
                       crude[crude==0] <- 1e-5
                       crude[crude>1] <- .999999

                       prelim.b <- coef(lm(log(crude) ~ log(age),data=dat))
                       return(c(exp(prelim.b[1]),prelim.b[2],1e-5))
                     },                     
                     optim.default=list(method="BFGS",
                                        control=list(reltol=1e-10,
                                                     parscale=c(2.78e-04, 1.879, .1))),
                     haz.fn=exweibull.haz.fn,
                     haz.to.prob.fn=haz.to.prob)
