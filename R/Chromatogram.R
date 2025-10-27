# Load required packages
library(R6)
library(ggplot2)
library(signal)
library(minpack.lm)
library(pracma)

`%||%` <- function(a, b) if (!is.null(a)) a else b

# R6 class definition
Chromatogram <- R6Class("Chromatogram",
  public = list(
    time = NULL,
    intensity = NULL,
    smoothed = NULL,
    fit_params = NULL,
    baseline_subtracted = NULL,
    peaks = NULL,

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
        
        # create a data frame for nlsLM to use
        df <- data.frame(
            time = self$time,
            smoothed = self$smoothed
            )
            
        # reasonable starting guesses
        start <- list(
                a = max(df$smoothed, na.rm = TRUE),
                b = df$time[which.max(df$smoothed)],
                c = (max(df$time) - min(df$time)) / 10
                )
                
        fit <- minpack.lm::nlsLM(smoothed ~ gaussian_model(time, a, b, c),
                data = df,
                start = start)
        self$fit_params <- coef(fit)
        message("Gaussian peak fitted successfully.")
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
    
    find_peaks = function(min_height = NULL, min_distance = 5) {
      if (is.null(self$smoothed)) {
        stop("Smooth the signal before finding peaks.")
      }

      peaks <- pracma::findpeaks(
        self$smoothed,
        minpeakheight = min_height %||% (0.1 * max(self$smoothed, na.rm = TRUE)),
        minpeakdistance = min_distance
      )

      if (is.null(peaks)) {
        message("⚠️ No peaks detected.")
        self$peaks <- NULL
      } else {
        colnames(peaks) <- c("height", "position", "start", "end")
        peaks_df <- as.data.frame(peaks)
        peaks_df$time <- self$time[peaks_df$position]
        self$peaks <- peaks_df
        message(nrow(peaks_df), " peaks detected.")
      }
      invisible(self)
    },

    integrate_peaks = function() {
      if (is.null(self$peaks)) stop("No peaks detected. Run $find_peaks() first.")
      results <- self$peaks
      results$area <- mapply(function(start, end) {
        start_idx <- max(1, start)
        end_idx <- min(length(self$time), end)
        pracma::trapz(self$time[start_idx:end_idx], self$smoothed[start_idx:end_idx])
      }, self$peaks$start, self$peaks$end)
      self$peaks <- results
      message("Peak integration complete.")
      invisible(self)
    },

    summarize_peaks = function() {
      if (is.null(self$peaks)) stop("No peak data available.")
      df <- self$peaks[, c("time", "height", "area")]
      print(df)
      return(df)
    },
    
    plot = function(save_path = NULL) {
        df <- data.frame(
            time = self$time,
            raw = self$intensity,
            smoothed = self$smoothed
        )

        p <- ggplot2::ggplot(df, ggplot2::aes(time)) +
            ggplot2::geom_line(ggplot2::aes(y = raw), alpha = 0.5, color = "gray40") +
            ggplot2::geom_line(ggplot2::aes(y = smoothed), color = "blue") +
            ggplot2::labs(title = "Chromatogram", x = "Time", y = "Intensity")

        # Add fitted curve only if fit_params exist
        if (!is.null(self$fit_params)) {
            gaussian_model <- function(x, a, b, c) {
                a * exp(-((x - b)^2) / (2 * c^2))
            }

            df$fit <- gaussian_model(
                df$time,
                self$fit_params["a"],
                self$fit_params["b"],
                self$fit_params["c"]
            )

        # Add fit layer only after 'fit' column exists
        p <- p + ggplot2::geom_line(
            data = df,
            ggplot2::aes(y = fit),
            color = "green",
            linetype = "dashed"
        )
    }

    # Save or print plot
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
