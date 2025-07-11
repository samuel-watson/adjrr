#' Estimates Relative Risks, Risk Differences, and Marginal Effects from Mixed Models
#'
#' @description
#' Calculates the marginal effect of variable x. There are several options for
#' marginal effect and several types of conditioning or averaging. The type of marginal
#' effect can be the derivative of the mean with respect to x (`dydx`), the expected
#' difference E(y|x=a)-E(y|x=b) (`diff`), or the expected log ratio log(E(y|x=a)/E(y|x=b)) (`ratio`).
#' Other fixed effect variables can be set at specific values (`at`), set at their mean values
#' (`atmeans`), or averaged over (`average`). Averaging over a fixed effects variable here means
#' using all observed values of the variable in the relevant calculation.
#' The random effects can similarly be set at their
#' estimated value (`re="estimated"`), set to zero (`re="zero"`), set to a specific value
#' (`re="at"`), or averaged over (`re="average"`). The standard errors are calculated using the delta method with one
#' of several options for the variance matrix of the fixed effect parameters.
#' Several of the arguments require the names
#' of the variables as given to the model object. Most variables are as specified in the formula,
#' factor variables are specified as the name of the `variable_value`, e.g. `t_1`. To see the names
#' of the stored parameters and data variables see the member function `names()`.
#' @param fit. Either a lme4, glmmTMB, or glmmrBase model fit.
#' @param x String. Name of the variable to calculate the marginal effect for.
#' @param type String. Either `dydx` for derivative, `diff` for difference, or `ratio` for log ratio. See description.
#' @param re String. Either `estimated` to condition on estimated values, `zero` to set to zero, `at` to
#' provide specific values, or `average` to average over the random effects.
#' @param se String. Type of standard error to use, either `GLS` for the GLS standard errors, `KR` for
#' Kenward-Roger estimated standard errors, `KR2` for the improved Kenward-Roger correction (see `small_sample_correction()`),
#'  or `robust` to use a robust sandwich estimator.
#' @param at Optional. A vector of strings naming the fixed effects for which a specified value is given.
#' @param atmeans Optional. A vector of strings naming the fixed effects that will be set at their mean value.
#' @param average Optional. A vector of strings naming the fixed effects which will be averaged over.
#' @param xvals Optional. A vector specifying the values of `a` and `b` for `diff` and `ratio`. The default is (1,0).
#' @param atvals Optional. A vector specifying the values of fixed effects specified in `at` (in the same order).
#' @param revals Optional. If `re="at"` then this argument provides a vector of values for the random effects.
#' @param oim Logical. If TRUE use the observed information matrix, otherwise use the expected information matrix for standard error and related calculations.
#' @param sampling Integer. Number of MCMC samples to use.
#' @return A named vector with elements `margin` specifying the point estimate and `se` giving the standard error.
#' @importFrom glmmrBase lme4_to_glmmr
#' @importFrom glmmrBase Model
#' @examples
#' TBC
#' @export
margin <- function(fit, ..., sampling = 250){

  args <- list(...)

  if(is(fit,"glmmTMB")){
    f1 <- glmmrBase::lme4_to_glmmr(fit$call$formula,names(df))
    if(fit$modelInfo$REML){
      mean_pars <- fit$fit$parfull[grepl("beta",names(fit$fit$parfull))]
    } else {
      mean_pars <- fit$fit$par[which(names(fit$fit$par)=="beta")]
    }
    model <- glmmrBase::Model$new(
      f1,
      data = fit$frame,
      mean = mean_pars,
      covariance = exp(fit$fit$par[which(names(fit$fit$par)=="theta")])^2,
      family = binomial()
    )
  } else if(is(fit,"glmerMod")){
    f1 <- glmmrBase::lme4_to_glmmr(fit@call$formula,names(df))
    model <- glmmrBase::Model$new(
      f1,
      data = fit@frame,#df[!is.na(df[,as.character(f1[[2]])]),],
      mean = fit@beta,
      covariance = fit@theta^2,
      family = binomial()
    )
  } else if(is(fit,"Model")) {
    model <- Model$new(fit)
  } else {
    stop("Fit should be glmerMod, glmmTMB, or Model class")
  }

  model$mcmc_options$samps <- sampling
  suppressMessages(suppressWarnings(model$mcmc_sample()))
  result <- model$marginal(...)
  out <- list(
    result = result,
    formula = f1,
    x = args$x,
    type = args$type,
    se = args$se,
    re = args$re,
    sampling = sampling
  )
  class(out) <- "margin"

  return(out)
}

#' Prints the marginal output
#'
#' @export
print.margin <- function(x, digits = 4){
  cat("Marginal Effects from Mixed Model")
  f1 <- as.character(x$formula)
  cat("\nFormula: ",f1[[2]]," ~ ",f1[[3]])
  res <- round(x$result$margin,digits)
  if(x$type == "ratio"){
    cat("\nLog RR: ",res)
  } else if(x$type == "diff"){
    cat("\nRisk diff.: ",res)
  } else {
    cat("\nMarginal effect (dydx): ",res)
  }
  cat(" SE: ",round(x$result$SE,digits))
}


#' Prints the marginal output
#'
#' @export
summary.margin <- function(x, digits = 3){
  cat("Marginal Effects from Mixed Model")
  f1 <- as.character(x$formula)
  cat("\nFormula: ",f1[[2]]," ~ ",f1[[3]])
  if(x$type == "ratio"){
    name <- "Log RR "
  } else if(x$type == "diff"){
    name <- "Risk diff. "
  } else {
    name <- "Marginal effect (dydx) "
  }
  name <- paste0(name, "(",x$x,")")
  cat("\n\n")
  out <- data.frame(est = round(x$result$margin, digits),
                    se = round(x$result$SE, digits),
                    z = round(x$result$margin/x$result$SE,digits))
  #print(out)
  colnames(out) <- c("Estimate","Std. Err.","z value")
  rownames(out) <- name
  print(out)
  cat("\nRE type: ",x$re,", SE: ",x$se,", MCMC samples: ",x$sampling)
}


