#' Compute attribution metrics from bootstrapped analogs
#'
#' Calculates standard extreme event attribution metrics — Probability Ratio (PR),
#' Fraction Attributable Risk (FAR), and mean temperature shift (delta_T) — from
#' the bootstrap simulation output of \code{\link{bs_analogs}}.
#' Confidence intervals are estimated by resampling the bootstrap distributions.
#'
#' @param x Output list from \code{\link{bs_analogs}}, containing
#'   \code{bootstrap_simulation} and \code{observed}.
#' @param threshold Numeric. The threshold value defining the extreme event
#'   (e.g. the observed event magnitude or a percentile value).
#' @param tail Character. Direction of the exceedance: \code{"upper"} for
#'   hot/wet events (\eqn{P(X \ge threshold)}), \code{"lower"} for cold/dry
#'   events (\eqn{P(X \le threshold)}). Default is \code{"upper"}.
#' @param cf_period Character. Label of the counterfactual period as it appears
#'   in \code{analogs$period}. Defaults to the first period in the data.
#' @param f_period Character. Label of the factual period. Defaults to the last
#'   period in the data.
#' @param conf_level Numeric. Confidence level for the intervals. Default \code{0.95}.
#' @param n_boot Integer. Number of resamples for CI estimation. Default \code{1000}.
#'
#' @return A \code{tibble} with one row and columns:
#'   \describe{
#'     \item{counterfactual_period}{Label of the counterfactual period}
#'     \item{factual_period}{Label of the factual period}
#'     \item{threshold}{The threshold used}
#'     \item{tail}{Direction of exceedance}
#'     \item{P_counterfactual}{Exceedance probability in the counterfactual world}
#'     \item{P_factual}{Exceedance probability in the factual world}
#'     \item{PR}{Probability Ratio (\eqn{P_f / P_{cf}})}
#'     \item{PR_low, PR_high}{Lower and upper CI bounds for PR}
#'     \item{FAR}{Fraction Attributable Risk (\eqn{1 - 1/PR})}
#'     \item{FAR_low, FAR_high}{CI bounds for FAR}
#'     \item{delta_T}{Mean shift between factual and counterfactual (\eqn{\mu_f - \mu_{cf}})}
#'   }
#'
#' @importFrom dplyr filter group_by summarise mutate
#' @importFrom tibble tibble
#' @importFrom stats quantile
#'
#' @export
#'
#' @examples
#' \dontrun{
#' metrics <- attribution_metrics(bs_tab, threshold = 35, tail = "upper")
#' print(metrics)
#' }
attribution_metrics <- function(x,
                                threshold,
                                tail       = c("upper", "lower"),
                                cf_period  = NULL,
                                f_period   = NULL,
                                conf_level = 0.95,
                                n_boot     = 1000) {

  tail <- match.arg(tail)

  if (!inherits(x, "climattr_bs")) {
    stop("x must be the output of bs_analogs() (a 'climattr_bs' object)")
  }

  # Aggregate bootstraps to event level (mean over all event days per sim/period)
  bs_event <- x$bootstrap_simulation |>
    dplyr::group_by(sim, period) |>
    dplyr::summarise(var_event = mean(var, na.rm = TRUE), .groups = "drop")

  yr_split <- unique(bs_event$period)
  if (length(yr_split) < 2) {
    stop("attribution_metrics() requires at least 2 periods in the bootstrap output.")
  }

  cf_period <- if (is.null(cf_period)) yr_split[1]                  else cf_period
  f_period  <- if (is.null(f_period))  yr_split[length(yr_split)]   else f_period

  cf_vals <- dplyr::filter(bs_event, period == cf_period)$var_event
  f_vals  <- dplyr::filter(bs_event, period == f_period)$var_event

  exceed_fun <- if (tail == "upper") {
    function(v, thr) mean(v >= thr, na.rm = TRUE)
  } else {
    function(v, thr) mean(v <= thr, na.rm = TRUE)
  }

  p_cf <- exceed_fun(cf_vals, threshold)
  p_f  <- exceed_fun(f_vals,  threshold)

  # Continuity correction: if P_cf = 0, use 1/(2*n) to obtain a finite PR
  # This is standard in attribution studies when the event is not observed
  # in the counterfactual sample. PR is then a lower bound.
  n_sims  <- length(cf_vals)
  p_cf_cc <- if (p_cf == 0) 1 / (2 * n_sims) else p_cf
  p_f_cc  <- if (p_f  == 0) 1 / (2 * n_sims) else p_f

  PR  <- p_f_cc / p_cf_cc
  FAR <- 1 - 1 / PR

  # CI via resampling the bootstrap distributions
  alpha   <- 1 - conf_level

  pr_boot <- replicate(n_boot, {
    cf_s <- sample(cf_vals, n_sims, replace = TRUE)
    f_s  <- sample(f_vals,  n_sims, replace = TRUE)
    p_c  <- exceed_fun(cf_s, threshold)
    p_fs <- exceed_fun(f_s,  threshold)
    p_c_cc  <- if (p_c  == 0) 1 / (2 * n_sims) else p_c
    p_fs_cc <- if (p_fs == 0) 1 / (2 * n_sims) else p_fs
    p_fs_cc / p_c_cc
  })

  pr_ci  <- quantile(pr_boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  far_ci <- 1 - 1 / pr_ci  # FAR is monotone in PR, so CI inverts

  tibble::tibble(
    counterfactual_period = cf_period,
    factual_period        = f_period,
    threshold             = threshold,
    tail                  = tail,
    P_counterfactual      = round(p_cf, 4),
    P_factual             = round(p_f,  4),
    PR                    = round(PR, 3),
    PR_low                = round(pr_ci[1], 3),
    PR_high               = round(pr_ci[2], 3),
    FAR                   = round(FAR, 3),
    FAR_low               = round(min(far_ci), 3),
    FAR_high              = round(max(far_ci), 3),
    delta_T               = round(mean(f_vals, na.rm = TRUE) - mean(cf_vals, na.rm = TRUE), 3)
  )
}
