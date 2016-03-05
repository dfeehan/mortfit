########################
# constant hazard object

consthaz.haz.fn <- function(theta, z) {
    alpha <- exp(theta[1])
    return(rep(alpha,length(z)))
}

## these starting values have been updated based on preliminary analysis
consthaz.haz <- new("mortalityHazard",
                    name="Constant Hazard",
                    num.param=1L,
                    theta.default=c(-2.48),
                    theta.range=list(c(-2.73, -2.29)),
                    optim.default=list(method="BFGS",
                                       control=list(reltol=1e-10)),
                    haz.fn=consthaz.haz.fn,
                    haz.to.prob.fn=Curry(haz.to.prob,
                                         haz.fn=consthaz.haz.fn))
