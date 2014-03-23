#######################################################
##
## classdefns-showmethods.R
##
## show methods, which are used for printing out contents
## of objects when 
##
## @include classdefns.R

setMethod("show",
          "mortalityHazard",
          function(object) {
            cat("------------------------------------\n")
            cat("mortalityHazard object: ", object@name, "\n")
            cat("          num.param : ", object@num.param, "\n")
            cat("      theta.default : ", object@theta.default, "\n")            
            cat("------------------------------------\n")            
          })

setMethod("show",
          "mortalityData",
          function(object) {
            cat("------------------------------------\n")
            cat("mortalityData object:", object@name, "\n")
            cat("   data has columns ", colnames(object@data), "\n")            
            cat("   data has dim ", paste(dim(object@data)), "\n")
            cat("   age offset is ", paste(object@age.offset), "\n")
            if (! is.null(object@tags)) {
              cat("   tags are: ",
                  paste(names(object@tags), object@tags,
                        sep="=", collapse="; "), "\n")
            }
            cat("------------------------------------\n")            
          })

setMethod("show",
          "mortalityDataFolded",
          function(object) {
            cat("------------------------------------\n")
            cat("mortalityDataFolded object: ", object@name, "\n")
            cat("   num.folds : ", object@num.folds, "\n")
            cat("   data has columns ", colnames(object@data), "\n")
            cat("   data has dim ", paste(dim(object@data)), "\n")
            cat("------------------------------------\n")            
          })

setMethod("show",
          "fitMethod",
          function(object) {
            cat("------------------------------------\n")
            cat("fitMethod object: ", object@name, "\n")
            cat("------------------------------------\n")            
          })

setMethod("show",
          "fitSummary",
          function(object) {
            cat("------------------------------------\n\n")
            cat("fitSummary object\n")
            cat("                AIC: ", object@AIC, "\n")
            cat("                BIC: ", object@BIC, "\n")
            cat("             SSE.Dx: ", object@SSE.Dx, "\n\n")
            cat("        chi-squared: ", object@chisq, "\n")
            cat("          (p value): ", object@chisq.p, "\n\n")
            cat("           deviance: ", object@deviance, "\n\n")

            chisq.byagemat <- data.frame(object@chisq.byage,
                                         object@chisq.byage.p)
            colnames(chisq.byagemat) <- c("chi-squared by age", "(p values)")
            print(chisq.byagemat, row.names=FALSE)

            cat("\n------------------------------------\n")            
          })

setMethod("show",
          "fitSummaryCV",
          function(object) {
            cat("------------------------------------\n\n")
            cat("fitSummary object\n")
            cat("                AIC: ", object@AIC, "\n")
            cat("                BIC: ", object@BIC, "\n")
            cat("             SSE.Dx: ", object@SSE.Dx, "\n\n")
            cat("        chi-squared: ", object@chisq, "\n")
            cat("          (p value): ", object@chisq.p, "\n\n")
            cat(" chi-squared by age: ", object@chisq.byage, "\n")
            cat("         (p values): ", object@chisq.byage.p, "\n\n")
            cat("           deviance: ", object@deviance, "\n\n")
            cat("         CV rmse.Dx: ", object@cv.rmse.Dx, "\n")
            cat("------------------------------------\n")            
          })

