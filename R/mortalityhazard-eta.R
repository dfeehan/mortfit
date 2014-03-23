########################
## ETA hazard object
##
## TODO -- double check this fn...
##         should it be -alpha in the numerator?
##         and need alpha be strictly positive?

## hazard function
eta.haz.fn <- function(theta, z) {
    alpha <- exp(theta[1])
    T <- theta[2]

    res <- (alpha / (1-exp(alpha*(z - T))))

    ## ensure only non-negative results are returned
    ## (hazards can't be negative)
    res[res<0] <- NA

    return(res)
}

eta.haz   <- new("mortalityHazard",
                 name="ETA",
                 num.param=2L,
                 ## TODO -- defaults here are tricky; second param is
                 ##  a maximum age...
                 theta.default=c(log(.1), 30),
                 optim.default=list(method="BFGS",
                                    control=list(parscale=c(log(.1), 30),
                                                 reltol=1e-8)),
                 haz.fn=eta.haz.fn,
                 haz.to.prob.fn=haz.to.prob)
