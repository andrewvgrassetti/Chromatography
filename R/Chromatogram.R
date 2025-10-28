# Load required packages
library(R6)
library(ggplot2)
library(signal)
library(minpack.lm)
library(pracma)

`%||%` <- function(a, b) if (!is.null(a)) a else b

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
      invisible(self)
    },

    smooth = function(window = 11, poly = 3) {
      if (length(self$intensity) >= window && window %% 2 == 1) {
        self$smoothed <- signal::sgolayfilt(self$intensity, p = poly, n = window)
      } else {
        # fallback: no smoothing if params are invalid for short series
        self$smoothed <- self$intensity
      }
      invisible(self)
    },

    fit_gaussian = function() {
      if (is.null(self$smoothed)) stop("Smooth the signal first with $smooth().")

      gaussian_model <- function(x, a, b, c)
        a * exp(-((x - b)^2) / (2 * c^2))

      df <- data.frame(time = self$time, smoothed = self$smoothed)

      start <- list(
        a = max(df$smoothed, na.rm = TRUE),
        b = df$time[which.max(df$smoothed)],
        c = (max(df$time) - min(df$time)) / 10
      )

      fit <- minpack.lm::nlsLM(smoothed ~ gaussian_model(time, a, b, c),
                               data = df, start = start)
      self$fit_params <- coef(fit)
      invisible(self)
    },

    baseline = function(method = c("min", "median")) {
      method <- match.arg(method)
      base <- if (method == "min") min(self$intensity) else median(self$intensity)
      self$baseline_subtracted <- self$intensity - base
      invisible(self)
    },

    auc = function() {
      pracma::trapz(self$time, self$intensity)
    },

    # NOTE: thresholds are in *minutes* (distance/width) and absolute intensity for height
    find_peaks = function(min_height = NULL,
                          min_distance_min = 0.2,
                          min_width_min = 0.1,
                          min_prom_frac = 0.05) {

      if (is.null(self$intensity) || is.null(self$time)) {
        stop("No data available for peak finding.")
      }

      dt <- suppressWarnings(median(diff(self$time), na.rm = TRUE))
      if (!is.finite(dt) || dt <= 0) dt <- 1

      dist_pts  <- max(1L, as.integer(ceiling(min_distance_min / dt)))
      width_pts <- max(1L, as.integer(ceiling(min_width_min   / dt)))

      y <- as.numeric(self$intensity)
      ymax <- max(y, na.rm = TRUE)
      min_height_abs <- min_height %||% (0.10 * ymax)
      min_prom_abs   <- min_prom_frac * ymax

      pk <- pracma::findpeaks(y,
                              minpeakheight   = min_height_abs,
                              minpeakdistance = dist_pts,
                              threshold       = min_prom_abs)

      if (is.null(pk)) {
        self$peaks <- NULL
        return(invisible(self))
      }

      colnames(pk) <- c("height","position","start","end")
      df <- as.data.frame(pk)

      n <- length(y)
      df$start    <- pmax(1L, pmin(n, as.integer(df$start)))
      df$end      <- pmax(1L, pmin(n, as.integer(df$end)))
      df$position <- pmax(1L, pmin(n, as.integer(df$position)))

      df$time       <- self$time[df$position]
      df$width_pts  <- pmax(0L, df$end - df$start + 1L)
      df$width_min  <- df$width_pts * dt
      bound_base    <- pmax(y[df$start], y[df$end])
      df$prominence <- df$height - bound_base

      keep <- (df$width_min >= min_width_min) &
              (df$prominence >= min_prom_abs) &
              is.finite(df$height) & is.finite(df$time)

      df <- df[keep, , drop = FALSE]
      df <- df[order(df$time), , drop = FALSE]

      self$peaks <- if (nrow(df)) df else NULL
      invisible(self)
    },

    integrate_peaks = function(min_width_pts = 3L) {
      if (is.null(self$peaks)) stop("No peaks detected. Run $find_peaks() first.")

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

      # add relative area (% of total)
      if (nrow(res)) {
        total <- sum(res$area, na.rm = TRUE)
        res$rel_area_pct <- if (total > 0) (res$area / total) * 100 else NA_real_
      }

      self$peaks <- res
      invisible(self)
    }
  )
)
