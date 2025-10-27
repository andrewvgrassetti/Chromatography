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
    
    find_peaks = function(min_height = NULL,
                          min_distance_min = 0.2,   # minimum spacing between peaks (minutes)
                          min_width_min = 0.1,      # minimum peak width (minutes)
                          min_prom_frac = 0.05) {   # minimum prominence, as fraction of max intensity

      if (is.null(self$intensity) || is.null(self$time)) {
        stop("No data available for peak finding.")
      }

      # sampling interval (minutes per point)
      dt <- suppressWarnings(median(diff(self$time), na.rm = TRUE))
      if (!is.finite(dt) || dt <= 0) dt <- 1

      # convert time-based thresholds to points
      dist_pts  <- max(1L, as.integer(ceiling(min_distance_min / dt)))
      width_pts <- max(1L, as.integer(ceiling(min_width_min   / dt)))

      y <- as.numeric(self$intensity)
      ymax <- max(y, na.rm = TRUE)
      # absolute minima for height and prominence
      min_height_abs <- min_height %||% (0.10 * ymax)
      min_prom_abs   <- min_prom_frac * ymax

      # raw peak candidates
      pk <- pracma::findpeaks(y,
                              minpeakheight  = min_height_abs,
                              minpeakdistance = dist_pts,
                              threshold = min_prom_abs)

      if (is.null(pk)) {
        message("⚠️ No peaks detected with current thresholds.")
        self$peaks <- NULL
        return(invisible(self))
      }

      # pk columns: height, position, start, end
      colnames(pk) <- c("height","position","start","end")
      df <- as.data.frame(pk)

      # clamp indices to valid range & compute time/width/prominence
      n <- length(y)
      df$start    <- pmax(1L, pmin(n, as.integer(df$start)))
      df$end      <- pmax(1L, pmin(n, as.integer(df$end)))
      df$position <- pmax(1L, pmin(n, as.integer(df$position)))

      df$time       <- self$time[df$position]
      df$width_pts  <- pmax(0L, df$end - df$start + 1L)
      df$width_min  <- df$width_pts * dt
      # simple prominence estimate vs bounds
      bound_base    <- pmax(y[df$start], y[df$end])
      df$prominence <- df$height - bound_base

      # strict filtering
      keep <- (df$width_min >= min_width_min) &
              (df$prominence >= min_prom_abs) &
              is.finite(df$height) & is.finite(df$time)

      df <- df[keep, , drop = FALSE]
      df <- df[order(df$time), , drop = FALSE]

      # finalize
      if (!nrow(df)) {
        message("⚠️ All candidate peaks filtered out; relax thresholds.")
        self$peaks <- NULL
      } else {
        self$peaks <- df
        message(nrow(df), " peaks retained after strict filtering.")
      }

      invisible(self)
    },

    integrate_peaks = function(min_width_pts = 3L) {
      if (is.null(self$peaks)) stop("No peaks detected. Run $find_peaks() first.")
      if (is.null(self$intensity)) stop("No intensity data found.")

      res <- self$peaks
      res$area <- mapply(function(s, e) {
        s <- as.integer(max(1, min(s, length(self$time))))
        e <- as.integer(max(1, min(e, length(self$time))))
        if ((e - s) < min_width_pts) return(NA_real_)
        x <- as.numeric(self$time[s:e]); y <- as.numeric(self$intensity[s:e])
        if (length(x) < 2 || anyNA(x) || anyNA(y)) return(NA_real_)
        tryCatch(pracma::trapz(x, y), error = function(.) NA_real_)
      }, res$start, res$end)

      res <- res[is.finite(res$area) & res$area > 0, , drop = FALSE]
      self$peaks <- res
      message("Peak integration complete: ", nrow(res), " valid peaks.")
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
