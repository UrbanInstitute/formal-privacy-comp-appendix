#' Calculate confidence interval overlap
#'
#' @param private_lower A numeric lower bound for a private estimate
#' @param private_upper A numeric upper bound for a private estimate
#' @param nonprivate_lower A numeric lower bound for a noisy estimate
#' @param nonprivate_upper A numeric upper bound for a noisy estimate
#'
#' @return A numeric confidence interval overlap
#'
compute_confidence_interval_overlap <- function(private_lower, private_upper, nonprivate_lower, nonprivate_upper) {
  0.5 * 
    (((pmin(private_upper, nonprivate_upper) - pmax(private_lower, nonprivate_lower)) / (nonprivate_upper - nonprivate_lower)) +
       ((pmin(private_upper, nonprivate_upper) - pmax(private_lower, nonprivate_lower)) / (private_upper - private_lower)))
}