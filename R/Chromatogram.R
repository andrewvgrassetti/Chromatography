# Load required packages
library(R6)
library(ggplot2)
library(minpack.lm)
library(pracma)

`%||%` <- function(a, b) if (!is.null(a)) a else b

Chromatogram <- R6Class("Chromatogram",
  public = list(
    time = NULL,
    intensity = NULL,
    fit_params = NULL,
    peaks = NULL,

    initialize = function(time, intensity) {
      stopifnot(length(time) == length(intensity))
      self$time <- as.numeric(time)
      self$intensity <- as.numeric(intensity)
      invisible(self)
    },

    # Optional: fit single Gaussian to RAW (no smoothing required)
    fit_gaussian = function() {
      y <- self$intensity
      x <- self$time
      if (length(x) < 3 || length(y) < 3) return(invisible(self))

      gaussian_model <- function(x, a, b, c) a * exp(-((x - b)^2) / (2 * c^2))
      start <- list(
        a = max(y, na.rm = TRUE),
        b = x[which.max(y)],
        c = (max(x) - min(x)) / 10
      )

      df <- data.frame(x = x, y = y)
      fit <- try(
        minpack.lm::nlsLM(y ~ gaussian_model(x, a, b, c), data = df, start = start),
        silent = TRUE
      )
      if (!inherits(fit, "try-error")) {
        self$fit_params <- coef(fit)
      }
      invisible(self)
    },

    # Area under the raw curve
    auc = function() {
      pracma::trapz(self$time, self$intensity)
    },

    # -------- Height-only peak picking ----------
    # Keeps every local maximum with intensity >= (min_height_abs),
    # and enforces only a minimum distance in minutes.
    find_peaks = function(min_height = NULL, min_distance_min = 0.2) {
      if (is.null(self$intensity) || is.null(self$time)) {
        stop("No data available for peak finding.")
      }

      y <- as.numeric(self$intensity)
      x <- as.numeric(self$time)

      ymax <- max(y, na.rm = TRUE)
      min_height_abs <- min_height %||% (0.10 * ymax)

      # sampling interval (minutes per point)
      dt <- suppressWarnings(median(diff(x), na.rm = TRUE))
      if (!is.finite(dt) || dt <= 0) dt <- 1

      dist_pts <- max(1L, as.integer(ceiling(min_distance_min / dt)))

      # Use pracma::findpeaks with ONLY height and distance; set threshold = 0
      pk <- pracma::findpeaks(
        y,
        minpeakheight   = min_height_abs,
        minpeakdistance = dist_pts,
        threshold       = 0
      )

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

      df$time      <- x[df$position]
      df$width_pts <- pmax(0L, df$end - df$start + 1L)
      df$width_min <- df$width_pts * dt

      df <- df[order(df$time), , drop = FALSE]
      self$peaks <- df
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
      self$peaks <- res
      invisible(self)
    },

    plot = function(save_path = NULL, show_peaks = TRUE, label_peaks = TRUE) {
      df <- data.frame(time = self$time, raw = self$intensity)
      p <- ggplot(df, aes(time, raw)) +
        geom_line(linewidth = 0.8, color = "#2c3e50") +
        labs(title = "Chromatogram", x = "Time (min)", y = "Intensity") +
        theme_classic(base_size = 13)

      # Fitted Gaussian
      if (!is.null(self$fit_params)) {
        a <- unlist(self$fit_params)[c("a","b","c")]
        if (all(is.finite(a))) {
          fit <- a["a"] * exp(-((df$time - a["b"])^2) / (2 * a["c"]^2))
          p <- p + geom_line(aes(y = fit), linetype = "dashed", linewidth = 0.8, color = "darkgreen")
        }
      }

      # Peaks
      if (show_peaks && !is.null(self$peaks) && nrow(self$peaks) > 0) {
        peaks_df <- self$peaks
        peaks_df$peak_id <- factor(seq_len(nrow(peaks_df)))
        p <- p +
          geom_vline(data = peaks_df, aes(xintercept = time, color = peak_id),
                     linetype = "dotted", linewidth = 0.8, alpha = 0.9) +
          scale_color_manual(values = scales::hue_pal()(nrow(peaks_df)), guide = "none")

        if (label_peaks) {
          p <- p + geom_text(
            data = peaks_df,
            aes(x = time, y = max(df$raw, na.rm = TRUE) * 1.05, label = sprintf("t=%.2f", time), color = peak_id),
            angle = 90, vjust = -0.4, size = 3, show.legend = FALSE
          )
        }
      }

      if (is.null(save_path)) {
        print(p)
      } else {
        dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave(save_path, p, width = 7, height = 4, dpi = 150)
      }
      invisible(p)
    }
  )
)

