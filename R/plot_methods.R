#' Print method for climattr_bs objects
#'
#' @param x A \code{climattr_bs} object (output of \code{\link{bs_analogs}}).
#' @param ... Further arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.climattr_bs <- function(x, ...) {
  cat("climattR bootstrap attribution result\n")
  cat("--------------------------------------\n")
  cat("Periods:", paste(unique(x$bootstrap_simulation$period), collapse = " | "), "\n")
  cat("Event days:", nrow(x$observed), "\n")
  n_sims <- length(unique(x$bootstrap_simulation$sim))
  cat("Bootstrap replicates:", n_sims, "\n\n")
  cat("Summary by event day and period:\n")
  print(x$summary_bs, n = Inf)
  invisible(x)
}

#' Plot factual vs. counterfactual bootstrap distributions
#'
#' Generates a density plot comparing the bootstrapped event distributions
#' across periods. The dashed line marks the observed event value.
#' If \code{threshold} is provided, a dotted red line is added.
#'
#' @param object A \code{climattr_bs} object (output of \code{\link{bs_analogs}}).
#' @param threshold Optional numeric. Draws a dotted red line at this value.
#' @param xlab Character. Label for the x axis. Default \code{"Anomaly"}.
#' @param title Character. Plot title. Default \code{"Attribution: factual vs. counterfactual"}.
#' @param ... Further arguments (ignored).
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom dplyr group_by summarise filter
#' @importFrom ggplot2 ggplot aes geom_density geom_vline annotate labs
#'   theme_bw scale_fill_manual scale_colour_manual theme element_text
#'   autoplot
#'
#' @export
autoplot.climattr_bs <- function(object,
                                 threshold = NULL,
                                 xlab  = "Anomaly",
                                 title = "Attribution: factual vs. counterfactual",
                                 ...) {

  # Event-level aggregation: mean over all event days per sim/period
  bs_event <- object$bootstrap_simulation |>
    dplyr::group_by(sim, period) |>
    dplyr::summarise(var_event = mean(var, na.rm = TRUE), .groups = "drop")

  obs_val <- mean(object$observed$var, na.rm = TRUE)

  periods   <- unique(bs_event$period)
  pal_fills <- c("#4575b4", "#d73027")          # blue = cf, red = factual
  pal_cols  <- c("#2c5f8a", "#b02014")

  p <- ggplot2::ggplot(bs_event,
                       ggplot2::aes(x = var_event,
                                    fill = period,
                                    colour = period)) +
    ggplot2::geom_density(alpha = 0.35, linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = obs_val,
                        linetype = "dashed", linewidth = 0.9, colour = "grey20") +
    ggplot2::annotate("text",
                      x = obs_val, y = Inf,
                      label = paste0("observed\n(", round(obs_val, 1), ")"),
                      vjust = 1.4, hjust = -0.05,
                      size = 3.2, colour = "grey20") +
    ggplot2::scale_fill_manual(values   = stats::setNames(pal_fills, periods)) +
    ggplot2::scale_colour_manual(values = stats::setNames(pal_cols,  periods)) +
    ggplot2::labs(x      = xlab,
                  y      = "Density",
                  title  = title,
                  fill   = "Period",
                  colour = "Period") +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"))

  if (!is.null(threshold)) {
    p <- p + ggplot2::geom_vline(xintercept = threshold,
                                 linetype = "dotted",
                                 linewidth = 0.8,
                                 colour = "#d73027") +
      ggplot2::annotate("text",
                        x = threshold, y = Inf,
                        label = paste0("threshold\n(", round(threshold, 1), ")"),
                        vjust = 1.4, hjust = 1.05,
                        size = 3.2, colour = "#d73027")
  }

  p
}
