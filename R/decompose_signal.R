#' Dynamic/thermodynamic decomposition of the attribution signal
#'
#' Decomposes the total climate change signal (factual minus counterfactual mean)
#' into a **dynamic** component (driven by changes in atmospheric circulation)
#' and a **thermodynamic** component (driven by background warming at fixed
#' circulation). The decomposition requires two runs of \code{\link{bs_analogs}}
#' on the same input: one with \code{detrend = FALSE} (raw) and one with
#' \code{detrend = TRUE} (detrended, circulation-only).
#'
#' @details
#' The decomposition follows :
#' \describe{
#'   \item{Total}{
#'     \eqn{\Delta X_\text{total} = \bar{X}^\text{raw}_f - \bar{X}^\text{raw}_{cf}}
#'   }
#'   \item{Dynamic}{
#'     \eqn{\Delta X_\text{dyn} = \bar{X}^\text{det}_f - \bar{X}^\text{det}_{cf}}
#'     — same analog dates, temperatures detrended to remove the changing background.
#'   }
#'   \item{Thermodynamic}{
#'     \eqn{\Delta X_\text{thermo} = \Delta X_\text{total} - \Delta X_\text{dyn}}
#'     — residual after removing the circulation contribution.
#'   }
#' }
#' Confidence intervals are estimated by block-bootstrap resampling of the
#' per-simulation means from both objects.
#'
#' @param x_raw A \code{climattr_bs} object from \code{\link{bs_analogs}} with
#'   \code{detrend = FALSE}.
#' @param x_detrended A \code{climattr_bs} object from \code{\link{bs_analogs}}
#'   with \code{detrend = TRUE}. Must have the same periods as \code{x_raw}.
#' @param cf_period Character. Label of the counterfactual period. Defaults to
#'   the first period found in the data.
#' @param f_period Character. Label of the factual period. Defaults to the last
#'   period found in the data.
#' @param conf_level Numeric. Confidence level for the intervals. Default \code{0.95}.
#' @param n_boot Integer. Number of bootstrap resamples for CI. Default \code{1000}.
#'
#' @return A \code{climattr_decomp} tibble with columns:
#'   \describe{
#'     \item{counterfactual_period, factual_period}{Period labels.}
#'     \item{component}{One of \code{"Total"}, \code{"Dynamic"}, \code{"Thermodynamic"}.}
#'     \item{signal}{Point estimate of the signal component (same units as the
#'       variable, e.g. °C).}
#'     \item{ci_low, ci_high}{Lower and upper confidence interval bounds.}
#'   }
#'
#' @seealso \code{\link{bs_analogs}}, \code{\link{plot_decomposition}}
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
#'
#' bs_raw <- bs_analogs(ts_data, analogs$analogs_subperiods, n = 500,
#'                      anom = TRUE, ref_period = c(1950, 2022), detrend = FALSE)
#' bs_det <- bs_analogs(ts_data, analogs$analogs_subperiods, n = 500,
#'                      anom = TRUE, ref_period = c(1950, 2022), detrend = TRUE)
#'
#' dcmp <- decompose_signal(bs_raw, bs_det)
#' print(dcmp)
#' plot_decomposition(dcmp)
#' }
decompose_signal <- function(x_raw,
                             x_detrended,
                             cf_period  = NULL,
                             f_period   = NULL,
                             conf_level = 0.95,
                             n_boot     = 1000) {

  if (!inherits(x_raw, "climattr_bs") || !inherits(x_detrended, "climattr_bs")) {
    stop("x_raw and x_detrended must be 'climattr_bs' objects (output of bs_analogs()).")
  }

  # Aggregate bootstrap to event level (mean over all event days per sim)
  agg_event <- function(obj) {
    obj$bootstrap_simulation |>
      dplyr::group_by(sim, period) |>
      dplyr::summarise(val = mean(var, na.rm = TRUE), .groups = "drop")
  }

  raw <- agg_event(x_raw)
  det <- agg_event(x_detrended)

  periods   <- unique(raw$period)
  cf_period <- if (is.null(cf_period)) periods[1]              else cf_period
  f_period  <- if (is.null(f_period))  periods[length(periods)] else f_period

  get_vals <- function(df, per) dplyr::filter(df, period == per)$val

  cf_raw <- get_vals(raw, cf_period)
  f_raw  <- get_vals(raw, f_period)
  cf_det <- get_vals(det, cf_period)
  f_det  <- get_vals(det, f_period)

  total_signal <- mean(f_raw, na.rm = TRUE) - mean(cf_raw, na.rm = TRUE)
  dyn_signal   <- mean(f_det, na.rm = TRUE) - mean(cf_det, na.rm = TRUE)
  thermo_signal <- total_signal - dyn_signal

  # Bootstrap CI: resample draws and recompute all three components
  alpha <- 1 - conf_level
  n     <- min(length(cf_raw), length(f_raw), length(cf_det), length(f_det))

  boot_mat <- replicate(n_boot, {
    cf_r <- sample(cf_raw, n, replace = TRUE)
    f_r  <- sample(f_raw,  n, replace = TRUE)
    cf_d <- sample(cf_det, n, replace = TRUE)
    f_d  <- sample(f_det,  n, replace = TRUE)
    tot   <- mean(f_r) - mean(cf_r)
    dyn   <- mean(f_d) - mean(cf_d)
    c(Total = tot, Dynamic = dyn, Thermodynamic = tot - dyn)
  })

  ci <- apply(boot_mat, 1, quantile,
              probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

  result <- tibble::tibble(
    counterfactual_period = cf_period,
    factual_period        = f_period,
    component             = c("Total", "Dynamic", "Thermodynamic"),
    signal                = round(c(total_signal, dyn_signal, thermo_signal), 3),
    ci_low                = round(ci[1, ], 3),
    ci_high               = round(ci[2, ], 3)
  )

  structure(result, class = c("climattr_decomp", class(result)))
}


#' Plot a dynamic/thermodynamic decomposition
#'
#' Produces a bar chart of the three signal components (Total, Dynamic,
#' Thermodynamic) from \code{\link{decompose_signal}}, with 95% CI error bars.
#'
#' @param x A \code{climattr_decomp} tibble from \code{\link{decompose_signal}}.
#' @param ylab Character. Y-axis label. Default \code{"Signal (°C)"}.
#' @param title Character. Plot title.
#' @param palette Named character vector of colours for the three components.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_errorbar geom_hline
#'   scale_fill_manual labs theme_bw theme element_text
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dcmp <- decompose_signal(bs_raw, bs_det)
#' plot_decomposition(dcmp)
#' }
plot_decomposition <- function(x,
                               ylab    = "Signal (\u00b0C)",
                               title   = "Dynamic / thermodynamic decomposition",
                               palette = c(Total         = "#2c3e50",
                                           Dynamic       = "#2980b9",
                                           Thermodynamic = "#e74c3c")) {

  if (!inherits(x, "climattr_decomp")) {
    stop("x must be the output of decompose_signal().")
  }

  # factor order: Total first
  x$component <- factor(x$component, levels = c("Total", "Dynamic", "Thermodynamic"))

  ggplot2::ggplot(x, ggplot2::aes(x = component, y = signal, fill = component)) +
    ggplot2::geom_col(width = 0.55, alpha = 0.85, colour = "white") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      width = 0.18, linewidth = 0.85, colour = "grey20"
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_fill_manual(values = palette, guide = "none") +
    ggplot2::labs(
      x     = NULL,
      y     = ylab,
      title = title,
      subtitle = paste0(
        x$counterfactual_period[1], " (CF)  \u2192  ",
        x$factual_period[1], " (F)"
      )
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey40"),
      axis.text.x   = ggplot2::element_text(size = 11)
    )
}
