#' Return level curves per climate period
#'
#' Fits a Generalized Extreme Value (GEV) distribution to each sub-period
#' from the output of \code{\link{bs_analogs}} and returns return levels for a
#' grid of return periods with bootstrap confidence intervals. Plotting the two
#' period curves side-by-side shows how the tail of the distribution has shifted
#' due to climate change.
#'
#' @details
#' For each period:
#' \enumerate{
#'   \item Aggregate bootstrap simulations to event level (mean over event days, one value per bootstrap replicate).
#'   \item Fit a GEV distribution via MLE to the \eqn{n} event-level values.
#'   \item Evaluate the fitted CDF at each requested return period \eqn{T}:
#'     \deqn{z_T = \text{GEV quantile at } p = 1 - 1/T}
#'   \item Confidence intervals: resample the bootstrap draws \code{n_boot}
#'     times, refit the GEV, recompute \eqn{z_T} → empirical CI quantiles.
#' }
#' The vertical distance between the two period curves at any given return
#' period is the **attribution signal** at that extremity.
#'
#' @param x A \code{climattr_bs} object from \code{\link{bs_analogs}}.
#' @param return_periods Numeric vector of return periods (in years / number of
#'   bootstrap samples) at which to evaluate return levels. Defaults to a
#'   log-spaced grid from 2 to 200.
#' @param conf_level Numeric. Confidence level for the bands. Default \code{0.95}.
#' @param n_boot Integer. Number of resamples for CI estimation. Default \code{500}.
#' @param tail Character. \code{"upper"} (default) for maxima-type events
#'   (hot/wet), \code{"lower"} for minima-type events (cold/dry). For lower
#'   events the GEV is fitted to \code{-x}.
#'
#' @return A \code{climattr_rl} tibble with columns:
#'   \describe{
#'     \item{period}{Period label from the bootstrap object.}
#'     \item{return_period}{Return period (same units as bootstrap sample size).}
#'     \item{return_level}{Point estimate of the return level.}
#'     \item{rl_low, rl_high}{Lower and upper CI bounds.}
#'   }
#'
#' @seealso \code{\link{plot_return_levels}}, \code{\link{evt_metrics}}
#'
#' @references
#' Coles, S. (2001). \emph{An Introduction to Statistical Modeling of Extreme
#' Values}. Springer.
#'
#' @importFrom dplyr group_by summarise filter bind_rows
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
#' rl <- return_level_curves(bs_tab)
#' plot_return_levels(rl, obs_val = mean(bs_tab$observed$var))
#' }
return_level_curves <- function(x,
                                return_periods = round(exp(seq(log(2), log(200), length.out = 40))),
                                conf_level     = 0.95,
                                n_boot         = 500,
                                tail           = c("upper", "lower")) {

  tail <- match.arg(tail)

  if (!requireNamespace("extRemes", quietly = TRUE)) {
    stop("Package 'extRemes' is required. Install with install.packages('extRemes').")
  }

  if (!inherits(x, "climattr_bs")) {
    stop("x must be the output of bs_analogs() (a 'climattr_bs' object).")
  }

  # Aggregate bootstrap simulations to one value per sim/period
  bs_event <- x$bootstrap_simulation |>
    dplyr::group_by(sim, period) |>
    dplyr::summarise(val = mean(var, na.rm = TRUE), .groups = "drop")

  periods <- unique(bs_event$period)
  alpha   <- 1 - conf_level

  # Probabilities corresponding to return periods
  p_exc <- 1 - 1 / return_periods   # non-exceedance prob for upper tail

  # Internal: fit GEV and extract return levels at p_exc
  gev_return_levels <- function(vals, p_exc, tail) {
    if (tail == "lower") vals <- -vals  # flip for minima

    fit <- tryCatch(
      extRemes::fevd(vals, type = "GEV", method = "MLE", verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(rep(NA_real_, length(p_exc)))

    par   <- fit$results$par
    loc   <- par["location"]
    scale <- par["scale"]
    shape <- par["shape"]

    # GEV quantile: loc + scale/shape * ((-log(p))^(-shape) - 1)  [shape != 0]
    rl <- vapply(p_exc, function(p) {
      y <- -log(p)   # reduced variate
      tryCatch({
        if (abs(shape) < 1e-6) {
          loc - scale * log(y)      # Gumbel limit
        } else {
          loc + scale / shape * (y^(-shape) - 1)
        }
      }, error = function(e) NA_real_)
    }, numeric(1))

    if (tail == "lower") rl <- -rl   # flip back
    rl
  }

  result_list <- lapply(periods, function(per) {
    vals  <- dplyr::filter(bs_event, period == per)$val
    n     <- length(vals)

    # Point estimate
    rl_est <- gev_return_levels(vals, p_exc, tail)

    # CI via resampling
    rl_boot <- replicate(n_boot, {
      s <- sample(vals, n, replace = TRUE)
      gev_return_levels(s, p_exc, tail)
    })
    # rl_boot is a matrix [n_rp x n_boot]
    rl_ci <- apply(rl_boot, 1, quantile,
                   probs = c(alpha / 2, 1 - alpha / 2),
                   na.rm = TRUE)

    tibble::tibble(
      period        = per,
      return_period = return_periods,
      return_level  = round(rl_est, 4),
      rl_low        = round(rl_ci[1, ], 4),
      rl_high       = round(rl_ci[2, ], 4)
    )
  })

  result <- dplyr::bind_rows(result_list)
  structure(result, class = c("climattr_rl", class(result)))
}


#' Plot return level curves for multiple climate periods
#'
#' Produces a return level plot (return period vs. return level) with one curve
#' per period and shaded confidence bands. A horizontal line marks the observed
#' event value, instantly showing its return period in each climate world.
#'
#' @param x A \code{climattr_rl} tibble from \code{\link{return_level_curves}}.
#' @param obs_val Numeric. Observed event value (e.g. mean anomaly). If
#'   supplied, a dashed horizontal line is drawn and the intersection with each
#'   curve is annotated. Default \code{NULL} (no line).
#' @param log_x Logical. Use log10 scale on the x-axis (return period)? Strongly
#'   recommended; default \code{TRUE}.
#' @param ylab Character. Y-axis label. Default \code{"Return level (°C)"}.
#' @param title Character. Plot title.
#' @param palette Named or unnamed character vector of colours for the periods.
#'   Defaults to blue (counterfactual) / red (factual).
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_hline annotate
#'   scale_colour_manual scale_fill_manual scale_x_log10 labs theme_bw theme
#'   element_text guides guide_legend
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rl <- return_level_curves(bs_tab, n_boot = 200)
#' plot_return_levels(rl, obs_val = 4.86)
#' }
plot_return_levels <- function(x,
                               obs_val = NULL,
                               log_x   = TRUE,
                               ylab    = "Return level (\u00b0C anomaly)",
                               title   = "Return level curves by climate period",
                               palette = NULL) {

  if (!inherits(x, "climattr_rl")) {
    stop("x must be the output of return_level_curves().")
  }

  periods <- unique(x$period)

  if (is.null(palette)) {
    # default: generate a blue-to-red ramp for 2 or more periods
    palette <- if (length(periods) == 2) {
      c("#4575b4", "#d73027")
    } else {
      grDevices::colorRampPalette(c("#4575b4", "#d73027"))(length(periods))
    }
    names(palette) <- periods
  }

  # Transparent fills (ribbon)
  fill_pal <- palette
  fill_alpha <- 0.18

  p <- ggplot2::ggplot(
    x,
    ggplot2::aes(x = return_period, colour = period, fill = period)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = rl_low, ymax = rl_high),
      alpha = fill_alpha, colour = NA
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = return_level),
      linewidth = 1.1
    ) +
    ggplot2::scale_colour_manual(values = palette, name = "Period") +
    ggplot2::scale_fill_manual(values = palette, name = "Period") +
    ggplot2::labs(
      x     = "Return period",
      y     = ylab,
      title = title
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )

  if (isTRUE(log_x)) {
    p <- p + ggplot2::scale_x_log10(
      breaks = c(2, 5, 10, 20, 50, 100, 200),
      labels = c("2", "5", "10", "20", "50", "100", "200")
    )
  }

  # Observed event line + annotation
  if (!is.null(obs_val)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = obs_val,
        linetype   = "dashed",
        colour     = "grey30",
        linewidth  = 0.7
      ) +
      ggplot2::annotate(
        "label",
        x     = min(x$return_period, na.rm = TRUE) * 1.1,
        y     = obs_val,
        label = paste0("Observed: ", round(obs_val, 2)),
        size  = 3.2,
        fill  = "white",
        colour = "grey30"
      )
  }

  p
}
