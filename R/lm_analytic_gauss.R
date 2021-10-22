#' Estimate a linear regression model with the analytic gaussian mechanism
#'
#' @param formula A formula for the regression model
#' @param data A data frame with the outcome variable and predictor variables
#' @param epsilon A numeric value for epsilon
#' @param delta A numeric value for delta
#' @param var_names A character vector with the names of the variables. The 
#' dependent variable should be last. 
#' @param var_types A named character vector with the types of variables. The
#' dependent variable should be last. 
#' @param bounds A named list with the variable bounds. 
#' @param reference_classes A named list with the reference classes for 
#' categorical predictors.
#' @param alpha A numeric for the alpha-level of the confidence intervals. 
#' Defaults to 0.05 for a 95% confidence interval.
#'
#' @return A list with two (asymptotic and bootstrap) tidy tibbles in the form 
#' of tidy.lm
#' 
lm_analytic_gauss <- function(formula, data, epsilon, delta, var_names, var_types, bounds, reference_classes, alpha = 0.05) {

  #### Getting DP estimates
  DPfit <- DpBoots_AGM.lm(
    Dataset = data, 
    Names = var_names, 
    Type = var_types, 
    Bounds = bounds, 
    Ref = reference_classes, 
    epsilon = epsilon,
    delta = delta,
    alpha = alpha
  )
  
  asymptotic_coef_table <- DPfit$asymptotic %>%
    rownames_to_column(var = "term") %>%
    as_tibble() %>%
    mutate(
      statistic = as.double(NA),
      p.value = as.double(NA)
    ) %>%
    select(
      term,
      estimate = Beta,
      std.error = se,
      statistic,
      p.value,
      conf.low = CI.l,
      conf.high = CI.u
    )

  bootstrap_coef_table <- DPfit$bootstrap %>%
    rownames_to_column(var = "term") %>%
    as_tibble() %>%
    mutate(
      statistic = as.double(NA),
      p.value = as.double(NA)
    ) %>%
    select(
      term,
      estimate = Beta,
      std.error = se,
      statistic,
      p.value,
      conf.low = CI.l,
      conf.high = CI.u
    )
  
  return(list(asymptotic = asymptotic_coef_table,
              bootstrap = bootstrap_coef_table))

}
