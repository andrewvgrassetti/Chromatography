# Load required packages
library(R6)
library(ggplot2)
library(signal)
library(minpack.lm)
library(pracma)

# R6 class definition
Chromatogram <- R6Class("Chromatogram",
  public = list(
    time = NULL,
    intensity = NULL,
    smoothed = NULL,
    fit_params = NULL,
    baseline_subtracted = NULL,

    initialize = function(time, intensity) {
      stopifnot(length(time) == length(intensity))
      self$time <- as.numeric(time)
      self$intensity <- as.numeric(intensity)
      message("Chromatogram object created with ", length(time), " points.")
    },

    smooth = function(window = 11, poly = 3) {
      self$smoothed <- signal::sgolayfilt(self$intensity, p = poly, n = window)
      message("Signal smoothed using Savitzky–Golay filter.")
      invisible(self)
    },

    fit_gaussian = function() {
      if (is.null(self$smoothed))
        stop("Smooth the signal first with $smooth().")

      gaussian_model <- function(x, a, b, c)
        a * exp(-((x - b)^2) / (2 * c^2))

      start <- list(a = max(self$smoothed),
                    b = self$time[which.max(self$smoothed)],
                    c = 1)

      fit <- minpack.lm::nlsLM(self$smoothed ~ gaussian_model(self$time, a, b, c),
                               start = start)
      self$fit_params <- coef(fit)
      message("Gaussian peak fitted.")
      invisible(self)
    },

    baseline = function(method = c("min", "median")) {
      method <- match.arg(method)
      base <- if (method == "min") min(self$intensity) else median(self$intensity)
      self$baseline_subtracted <- self$intensity - base
      message(sprintf("Baseline (%s) subtracted.", method))
      invisible(self)
    },

    auc = function(use = c("smoothed", "raw", "baseline")) {
      use <- match.arg(use)
      y <- switch(use,
        smoothed = self$smoothed,
        raw = self$intensity,
        baseline = self$baseline_subtracted
      )
      if (is.null(y)) stop("No data available for ", use)
      area <- pracma::trapz(self$time, y)
      message(sprintf("AUC (%s): %.3f", use, area))
      return(area)
    },

    plot = function(save_path = NULL) {
      df <- data.frame(time = self$time,
                       raw = self$intensity,
                       smoothed = self$smoothed)
      p <- ggplot2::ggplot(df, ggplot2::aes(time)) +
        ggplot2::geom_line(ggplot2::aes(y = raw), alpha = 0.5, color = "gray40") +
        ggplot2::geom_line(ggplot2::aes(y = smoothed), color = "blue") +
        ggplot2::labs(title = "Chromatogram", x = "Time", y = "Intensity")

      if (!is.null(self$fit_params)) {
        fit_curve <- self$fit_params["a"] *
          exp(-((self$time - self$fit_params["b"])^2) /
              (2 * self$fit_params["c"]^2))
        df$fit <- fit_curve
        p <- p + ggplot2::geom_line(ggplot2::aes(y = fit),
                                    color = "green", linetype = "dashed")
      }

      if (is.null(save_path)) {
        print(p)
      } else {
        dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave(save_path, p, width = 7, height = 4, dpi = 150)
        message("Plot saved to: ", save_path)
      }
      invisible(p)
    }
  )
)
