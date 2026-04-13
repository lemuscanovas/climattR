#' Extreme Value Theory attribution metrics
#'
#' Fits a parametric extreme value distribution (GEV or GPD) to the bootstrap
#' simulation output from \code{\link{bs_analogs}} and computes attribution
#' metrics (PR, FAR, return periods) with confidence intervals. Parametric
#' fitting yields smoother probability estimates than counting exceedances,
#' especially in the distribution tails.
#'
#' @details
#' **GEV mode** (\code{type = "GEV"}): the Generalized Extreme Value
#' distribution is fitted via MLE to all bootstrap samples for each period.
#' Exceedance probabilities are then read from the fitted survival function.
#'
#' **GPD mode** (\code{type = "GP"}): a Generalized Pareto Distribution is
#' fitted to values above the \code{gpd_quantile} threshold of each period's
#' distribution. Exceedance probability is \eqn{P(X > u) \times P(Y > thr - u)}
#' where \eqn{u} is the data threshold and \eqn{Y} follows the fitted GPD.
#'
#' Confidence intervals are computed by resampling the bootstrap vectors and
#' refitting the distribution \code{n_boot} times.
#'
#' @param x A \code{climattr_bs} object from \code{\link{bs_analogs}}.
#' @param threshold Numeric. The threshold defining the extreme event.
#' @param tail Character. \code{"upper"} for hot/wet events
#'   (\eqn{P(X \ge threshold)}), \code{"lower"} for cold/dry events. Default
#'   \code{"upper"}.
#' @param type Character. Distribution type: \code{"GEV"} (default) or
#'   \code{"GP"} (Generalized Pareto).
#' @param gpd_quantile Numeric in (0, 1). Quantile used as the data threshold
#'   for GPD fitting. Only used when \code{type = "GP"}. Default \code{0.75}.
#' @param cf_period Character. Label of the counterfactual period. Defaults to
#'   the first period in the data.
#' @param f_period Character. Label of the factual period. Defaults to the last
#'   period in the data.
#' @param conf_level Numeric. Confidence level for the intervals. Default
#'   \code{0.95}.
#' @param n_boot Integer. Number of resamples for CI estimation. Default
#'   \code{500}. Each resample refits the distribution, so keep this moderate
#'   for speed.
#'
#' @return A \code{tibble} with columns:
#'   \describe{
#'     \item{counterfactual_period, factual_period}{Period labels.}
#'     \item{threshold, tail, dist_type}{Input parameters.}
#'     \item{P_counterfactual, P_factual}{Parametric exceedance probabilities.}
#'     \item{return_period_cf, return_period_f}{Return periods (years) in each
#'       world.}
#'     \item{PR, PR_low, PR_high}{Probability Ratio and 95% CI.}
#'     \item{FAR, FAR_low, FAR_high}{Fraction Attributable Risk and 95% CI.}
#'   }
#'
#' @seealso \code{\link{attribution_metrics}} for the empirical (counting)
#'   equivalent.
#'
#' @references
#' Coles, S. (2001). \emph{An Introduction to Statistical Modeling of Extreme
#' Values}. Springer.
#'
#' @importFrom dplyr group_by summarise filter
#' @importFrom tibble tibble
#' @importFrom stats quantile
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Z500 <- terra::rast(system.file("extdata", "z500_0509_1950_2023_eu.nc", package = "climattR"))
#' TX   <- terra::rast(system.file("extdata", "tx_0608_1950_2023.nc",      package = "climattR"))
#' prep    <- prepare_data(Z500, event_dates = as.Date("2023-07-12"), time_window = 31)
#' analogs <- analogs_searcher(list(prep$ts_wo_event), list(prep$event),
#'                             n = 20, periods = c(1951, 1980, 1993, 2022))
#' ts_data <- as_ts(TX)
#' bs_tab  <- bs_analogs(ts_data, analogs$analogs_subperiods, n = 500,
#'                       anom = TRUE, ref_period = c(1950, 2022))
#'
#' obs_mean <- mean(bs_tab$observed$var, na.rm = TRUE)
#' evt_metrics(bs_tab, threshold = obs_mean, type = "GEV")
#' evt_metrics(bs_tab, threshold = obs_mean, type = "GP")
#' }
evt_metrics <- function(x,
                        threshold,
                        tail         = c("upper", "lower"),
                        type         = c("GEV", "GP"),
                        gpd_quantile = 0.75,
                        cf_period    = NULL,
                        f_period     = NULL,
                        conf_level   = 0.95,
                        n_boot       = 500) {

  tail <- match.arg(tail)
  type <- match.arg(type)

  if (!requireNamespace("extRemes", quietly = TRUE)) {
    stop("Package 'extRemes' is required. Install with install.packages('extRemes').")
  }

  if (!inherits(x, "climattr_bs")) {
    stop("x must be the output of bs_analogs() (a 'climattr_bs' object).")
  }

  # Aggregate bootstraps to event level (mean over days per sim/period)
  bs_event <- x$bootstrap_simulation |>
    dplyr::group_by(sim, period) |>
    dplyr::summarise(var_event = mean(var, na.rm = TRUE), .groups = "drop")

  periods   <- unique(bs_event$period)
  if (length(periods) < 2) {
    stop("evt_metrics() requires at least 2 periods in the bootstrap output.")
  }

  cf_period <- if (is.null(cf_period)) periods[1]              else cf_period
  f_period  <- if (is.null(f_period))  periods[length(periods)] else f_period

  cf_vals <- dplyr::filter(bs_event, period == cf_period)$var_event
  f_vals  <- dplyr::filter(bs_event, period == f_period)$var_event

  # ------------------------------------------------------------------
  # Internal: fit distribution and return exceedance probability
  # ------------------------------------------------------------------
  fit_exceedance <- function(vals, thr, tail, type, gpd_quantile) {
    tryCatch({
      if (type == "GEV") {
        fit <- extRemes::fevd(vals, type = "GEV", method = "MLE", verbose = FALSE)
        p_exceed <- extRemes::pextRemes(fit, thr, lower.tail = FALSE)
        if (tail == "lower") p_exceed <- 1 - p_exceed

      } else {  # GPD
        u      <- quantile(vals, gpd_quantile, na.rm = TRUE)
        p_over <- mean(vals > u, na.rm = TRUE)   # P(X > u) empirically
        exc    <- vals[vals > u] - u

        if (length(exc) < 10) {
          # Not enough exceedances — fall through to empirical
          stop("Too few exceedances for GPD fit.")
        }

        fit <- extRemes::fevd(exc, type = "GP", method = "MLE", verbose = FALSE)
        thr_shifted <- thr - u

        if (thr_shifted <= 0) {
          # threshold below GPD base — use full empirical
          p_exceed <- mean(vals >= thr, na.rm = TRUE)
        } else {
          p_gpd    <- extRemes::pextRemes(fit, thr_shifted, lower.tail = FALSE)
          p_exceed <- p_over * p_gpd
        }

        if (tail == "lower") p_exceed <- 1 - p_exceed
      }

      p_exceed

    }, error = function(e) {
      # Fallback: empirical probability
      if (tail == "upper") mean(vals >= thr, na.rm = TRUE)
      else                 mean(vals <= thr, na.rm = TRUE)
    })
  }

  p_cf <- fit_exceedance(cf_vals, threshold, tail, type, gpd_quantile)
  p_f  <- fit_exceedance(f_vals,  threshold, tail, type, gpd_quantile)

  # Continuity correction (same as attribution_metrics)
  n_sims  <- min(length(cf_vals), length(f_vals))
  cc      <- 1 / (2 * n_sims)
  p_cf_cc <- if (p_cf == 0) cc else p_cf
  p_f_cc  <- if (p_f  == 0) cc else p_f

  PR  <- p_f_cc / p_cf_cc
  FAR <- 1 - 1 / PR

  # Return periods
  rp_cf <- 1 / p_cf_cc
  rp_f  <- 1 / p_f_cc

  # ------------------------------------------------------------------
  # CI: resample bootstrap draws, refit distribution, recompute PR
  # ------------------------------------------------------------------
  alpha   <- 1 - conf_level

  pr_boot <- replicate(n_boot, {
    cf_s <- sample(cf_vals, n_sims, replace = TRUE)
    f_s  <- sample(f_vals,  n_sims, replace = TRUE)
    pc   <- fit_exceedance(cf_s, threshold, tail, type, gpd_quantile)
    pf   <- fit_exceedance(f_s,  threshold, tail, type, gpd_quantile)
    pc_c <- if (pc == 0) cc else pc
    pf_c <- if (pf == 0) cc else pf
    pf_c / pc_c
  })

  pr_ci  <- quantile(pr_boot, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  far_ci <- 1 - 1 / pr_ci   # FAR monotone in PR

  tibble::tibble(
    counterfactual_period = cf_period,
    factual_period        = f_period,
    threshold             = threshold,
    tail                  = tail,
    dist_type             = type,
    P_counterfactual      = round(p_cf, 6),
    P_factual             = round(p_f,  6),
    return_period_cf      = round(rp_cf, 1),
    return_period_f       = round(rp_f,  1),
    PR                    = round(PR, 3),
    PR_low                = round(pr_ci[1], 3),
    PR_high               = round(pr_ci[2], 3),
    FAR                   = round(FAR, 3),
    FAR_low               = round(min(far_ci), 3),
    FAR_high              = round(max(far_ci), 3)
  )
}
