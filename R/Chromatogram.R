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

    auc = function(use = "raw") {
      area <- pracma::trapz(self$time, self$intensity)
      message(sprintf("AUC: %f", area))
      return(area)
      },
    
    find_peaks = function(min_height = NULL, min_distance = 5) {
      # Use raw intensity instead of smoothed data
      if (is.null(self$intensity)) {
        stop("No intensity data available for peak finding.")
      }

      # Use pracma::findpeaks on the raw signal
      peaks <- pracma::findpeaks(
        self$intensity,
        minpeakheight = min_height %||% (0.1 * max(self$intensity, na.rm = TRUE)),
        minpeakdistance = min_distance
      )

      if (is.null(peaks)) {
        message("⚠️ No peaks detected.")
        self$peaks <- NULL
      } else {
        colnames(peaks) <- c("height", "position", "start", "end")
        peaks_df <- as.data.frame(peaks)
        peaks_df$time <- self$time[peaks_df$position]
        # Filter out peaks narrower than 3 points
        peaks_df <- peaks_df[(peaks_df$end - peaks_df$start) >= 3, , drop = FALSE]
        self$peaks <- peaks_df  
        message(nrow(peaks_df), " peaks detected after filtering.")
        }
      
      invisible(self)
      },

    integrate_peaks = function(min_width = 3) {
      if (is.null(self$peaks)) stop("No peaks detected. Run $find_peaks() first.")
      if (is.null(self$intensity)) stop("No intensity data found.")

      results <- self$peaks

      results$area <- mapply(function(start, end) {
        # Ensure numeric indices within bounds
        start_idx <- as.integer(max(1, min(start, length(self$time))))
        end_idx   <- as.integer(max(1, min(end, length(self$time))))

        # Skip peaks that are too narrow
        if ((end_idx - start_idx) < min_width) return(NA_real_)

        # Subset safely
        x <- as.numeric(self$time[start_idx:end_idx])
        y <- as.numeric(self$intensity[start_idx:end_idx])

        # Validate real numeric data
        if (length(x) < 2 || anyNA(x) || anyNA(y)) return(NA_real_)

        # Integrate
        area <- tryCatch(pracma::trapz(x, y), error = function(e) NA_real_)
        area
      }, results$start, results$end)

      # Drop invalid (tiny) peaks
      valid <- which(!is.na(results$area) & results$area > 0)
      results <- results[valid, , drop = FALSE]

      self$peaks <- results
      message("Peak integration complete: ", nrow(results), " valid peaks.")
      invisible(self)
    },

    summarize_peaks = function() {
      if (is.null(self$peaks)) stop("No peak data available.")
      df <- self$peaks[, c("time", "height", "area")]
      print(df)
      return(df)
    },
    
    plot = function(save_path = NULL, show_peaks = TRUE, label_peaks = TRUE) {
  df <- data.frame(
    time = self$time,
    raw = self$intensity,
    smoothed = self$smoothed
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(time)) +
    ggplot2::geom_line(ggplot2::aes(y = raw), alpha = 0.5, color = "gray40") +
    ggplot2::geom_line(ggplot2::aes(y = smoothed), color = "blue") +
    ggplot2::labs(title = "Chromatogram", x = "Time", y = "Intensity")

  # Add fitted Gaussian curve if available
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

    p <- p + ggplot2::geom_line(
      data = df,
      ggplot2::aes(y = fit),
      color = "green",
      linetype = "dashed"
    )
  }

  # Add vertical lines for peaks (if found)
  if (show_peaks && !is.null(self$peaks) && nrow(self$peaks) > 0) {
    peaks_df$peak_id <- factor(seq_len(nrow(peaks_df)))

      p <- p + ggplot2::geom_vline(
        data = peaks_df,
        ggplot2::aes(xintercept = time, color = peak_id),
        linetype = "dotted",
        linewidth = 0.7,
        alpha = 0.9
      ) +
      ggplot2::scale_color_manual(
        values = scales::hue_pal()(nrow(peaks_df)),
        guide = "none"
      )

    if (label_peaks) {
      p <- p + ggplot2::geom_text(
        data = self$peaks,
        ggplot2::aes(
          x = time,
          y = height + 0.05 * max(df$smoothed, na.rm = TRUE),
          label = sprintf("t=%.2f", time)
        ),
        color = "red",
        angle = 90,
        vjust = -0.5,
        size = 3
      )
    }
  }

  # Save or display the plot
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
