# app.R — Multi-file chromatogram explorer with overlay/separate plots and summary CSV

# ---- Packages ----
library(shiny)
library(bslib)
library(readr)
library(ggplot2)
library(colourpicker)
library(scales)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)

suppressWarnings({
  if (!requireNamespace("ragg", quietly = TRUE)) {
    message("Tip: install.packages('ragg') for high-quality TIFF/PNG exports.")
  }
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Load R6 backend ----
source("R/Chromatogram.R")

# ---- Helpers ----

read_two_col_csv <- function(path) {
  suppressWarnings({
    df <- tryCatch(
      readr::read_csv(path, col_names = FALSE, show_col_types = FALSE),
      error = function(e) NULL
    )
  })
  if (is.null(df) || ncol(df) < 2) return(NULL)
  df <- df[, 1:2]
  names(df) <- c("time", "intensity")
  df <- df[complete.cases(df$time, df$intensity), ]
  df$time <- suppressWarnings(as.numeric(df$time))
  df$intensity <- suppressWarnings(as.numeric(df$intensity))
  df <- df[!is.na(df$time) & !is.na(df$intensity), ]
  df
}

build_single_plot <- function(chrom, label, color_line, show_fit, show_peaks, label_peaks, show_area_pct) {
  df <- data.frame(time = chrom$time, intensity = chrom$intensity, label = label)

  p <- ggplot(df, aes(x = time, y = intensity)) +
    geom_line(linewidth = 0.9, color = color_line) +
    labs(x = "Time (min)", y = "Intensity")

  # (Optional) overlay fit
  if (show_fit && !is.null(chrom$fit_params)) {
    gaussian <- function(x, a, b, c) a * exp(-((x - b)^2) / (2 * c^2))
    df$fit <- gaussian(df$time, chrom$fit_params["a"], chrom$fit_params["b"], chrom$fit_params["c"])
    if (any(is.finite(df$fit))) {
      p <- p + geom_line(data = df, aes(y = fit), linewidth = 0.8, linetype = "dashed")
    }
  }

  # peaks & labels
  if (show_peaks && !is.null(chrom$peaks) && nrow(chrom$peaks) > 0) {
    pk <- chrom$peaks
    p <- p + geom_vline(data = pk, aes(xintercept = time), linetype = "dotted", linewidth = 0.6)

    if (label_peaks) {
      lbl <- pk
      ymax <- max(df$intensity, na.rm = TRUE)
      lbl$lab <- if (show_area_pct && "rel_area_pct" %in% names(lbl)) {
        sprintf("t=%.2f, %.1f%%", lbl$time, lbl$rel_area_pct)
      } else {
        sprintf("t=%.2f", lbl$time)
      }
      p <- p + geom_text(
        data = lbl, aes(x = time, y = ymax * 1.04, label = lab),
        angle = 0, vjust = 0, size = 3
      )
    }
  }

  p +
    theme_classic(base_size = 13) +
    theme(plot.title.position = "plot")
}

# ---- UI ----
ui <- page_fillable(
  theme = bs_theme(bootswatch = "flatly"),
  layout_sidebar(
    sidebar = sidebar(
      h4("Data"),
      fileInput("csvs", "Upload CSV files (2 columns: time, intensity)", multiple = TRUE, accept = ".csv"),
      checkboxInput("use_demo", "Use demo data if no files uploaded", TRUE),
      hr(),

      h4("Display"),
      radioButtons("display_mode", "Mode", choices = c("Overlay all" = "overlay", "Separate panels" = "separate"), inline = TRUE),
      checkboxInput("show_fit",   "Show fitted Gaussian", TRUE),
      checkboxInput("show_peaks", "Show peaks", TRUE),
      checkboxInput("label_peaks","Label peaks (horizontal)", TRUE),
      checkboxInput("show_area_pct","Include % area in labels", TRUE),
      hr(),

      h4("Peak Picking"),
      sliderInput("min_height_pct", "Min height (% of max)", min = 0, max = 100, value = 10, step = 1),
      numericInput("min_dist_min", "Min distance (minutes)", value = 0.2, min = 0, step = 0.05),
      numericInput("min_width_min","Min width (minutes)", value = 0.10, min = 0, step = 0.05),
      hr(),

      h4("Appearance"),
      colourInput("base_col", "Base color (overlay lines auto-colored)", value = "#2c3e50"),
      hr(),

      h4("Nicknames"),
      uiOutput("nickname_ui"),
      hr(),

      h4("Export"),
      downloadButton("download_summary", "Download summary CSV")
    ),

    card(
      card_header("Chromatograms"),
      plotOutput("main_plot", height = 500)
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  # Reactive: list of datasets (either uploaded or demo)
  files_df <- reactive({
    if (!is.null(input$csvs) && nrow(input$csvs) > 0) {
      tibble::tibble(
        key = seq_len(nrow(input$csvs)),
        path = input$csvs$datapath,
        filename = input$csvs$name
      )
    } else if (isTRUE(input$use_demo)) {
      # create two demo sets
      demo1 <- tibble::tibble(
        key = 1L,
        path = NA_character_,
        filename = "demo_A.csv"
      )
      demo2 <- tibble::tibble(
        key = 2L,
        path = NA_character_,
        filename = "demo_B.csv"
      )
      dplyr::bind_rows(demo1, demo2)
    } else {
      tibble::tibble(key = integer(), path = character(), filename = character())
    }
  })

  # Nickname inputs (dynamic)
  output$nickname_ui <- renderUI({
    fdf <- files_df()
    if (nrow(fdf) == 0) return(helpText("No files loaded yet."))
    tagList(lapply(seq_len(nrow(fdf)), function(i) {
      textInput(paste0("nick_", fdf$key[i]),
                label = sprintf("Nickname for %s", fdf$filename[i]),
                value = tools::file_path_sans_ext(fdf$filename[i]))
    }))
  })

  # Build Chromatogram objects for each dataset with current controls
  chrom_list <- reactive({
    fdf <- files_df()
    if (nrow(fdf) == 0) return(list())

    res <- vector("list", nrow(fdf))

    for (i in seq_len(nrow(fdf))) {
      # load data
      if (is.na(fdf$path[i])) {
        # demo data
        set.seed(100 + i)
        n <- 800
        time <- seq(0, 20, length.out = n)
        peak1 <- (80 + 10*i) * exp(-((time - (9.5 + 0.2*i))^2) / (2 * 0.25^2))
        peak2 <- (60 + 10*i) * exp(-((time - (10.6 + 0.2*i))^2) / (2 * 0.22^2))
        intensity <- peak1 + peak2 + rnorm(n, sd = 4)
        df <- data.frame(time = time, intensity = intensity)
      } else {
        df <- read_two_col_csv(fdf$path[i])
        if (is.null(df)) next
      }

      # build Chromatogram
      ch <- Chromatogram$new(df$time, df$intensity)
      ch$smooth(window = 11, poly = 3)

      # min_height as % of max of *smoothed* (robuster than raw)
      h_pct <- (input$min_height_pct %||% 10) / 100
      max_ref <- max(ch$smoothed %||% ch$intensity, na.rm = TRUE)
      min_h_abs <- h_pct * max_ref

      # find/integrate peaks
      ch$find_peaks(
        min_height      = min_h_abs,
        min_distance_min = input$min_dist_min %||% 0.2,
        min_width_min    = input$min_width_min %||% 0.10,
        min_prom_frac    = 0.05
      )
      if (!is.null(ch$peaks)) ch$integrate_peaks()

      # store label (nickname or filename)
      nick_val <- input[[paste0("nick_", fdf$key[i])]]
      label <- nick_val %||% fdf$filename[i]

      res[[i]] <- list(
        filename = fdf$filename[i],
        label = label,
        chrom = ch
      )
    }
    Filter(Negate(is.null), res)
  })

  # Assign colors per series (legend)
  series_colors <- reactive({
    lst <- chrom_list()
    n <- length(lst)
    setNames(hue_pal()(n), vapply(lst, function(z) z$label, ""))
  })

  # Main plot: overlay or separate
  output$main_plot <- renderPlot({
    lst <- chrom_list()
    validate(need(length(lst) > 0, "Load at least one CSV (or enable demo)."))

    cols <- series_colors()

    if (identical(input$display_mode, "overlay") || length(lst) == 1) {
      # overlay: bind rows with label and color by label
      df_all <- dplyr::bind_rows(lapply(lst, function(z) {
        data.frame(time = z$chrom$time, intensity = z$chrom$intensity, label = z$label, stringsAsFactors = FALSE)
      }))

      p <- ggplot(df_all, aes(x = time, y = intensity, color = label)) +
        geom_line(linewidth = 0.9) +
        scale_color_manual(values = cols) +
        labs(x = "Time (min)", y = "Intensity", color = "Sample") +
        theme_classic(base_size = 13)

      # peaks + horizontal labels
      if (isTRUE(input$show_peaks)) {
        pk_all <- dplyr::bind_rows(lapply(lst, function(z) {
          if (is.null(z$chrom$peaks) || !nrow(z$chrom$peaks)) return(NULL)
          cbind(z$chrom$peaks, label = z$label, stringsAsFactors = FALSE)
        }))
        if (!is.null(pk_all) && nrow(pk_all) > 0) {
          ymax <- df_all %>% dplyr::group_by(label) %>% dplyr::summarise(ymax = max(intensity, na.rm = TRUE))
          pk_all <- dplyr::left_join(pk_all, ymax, by = "label")
          p <- p +
            geom_vline(data = pk_all, aes(xintercept = time, color = label), linetype = "dotted", linewidth = 0.6)

          if (isTRUE(input$label_peaks)) {
            pk_all$lab <- if (isTRUE(input$show_area_pct) && "rel_area_pct" %in% names(pk_all)) {
              sprintf("t=%.2f, %.1f%%", pk_all$time, pk_all$rel_area_pct)
            } else {
              sprintf("t=%.2f", pk_all$time)
            }
            p <- p + geom_text(
              data = pk_all,
              aes(x = time, y = ymax * 1.04, label = lab, color = label),
              angle = 0, vjust = 0, size = 3, show.legend = FALSE
            )
          }
        }
      }

      p

    } else {
      # separate panels: facet by label
      df_all <- dplyr::bind_rows(lapply(lst, function(z) {
        data.frame(time = z$chrom$time, intensity = z$chrom$intensity, label = z$label, stringsAsFactors = FALSE)
      }))

      p <- ggplot(df_all, aes(x = time, y = intensity)) +
        geom_line(aes(color = label), linewidth = 0.9, show.legend = FALSE) +
        scale_color_manual(values = cols) +
        labs(x = "Time (min)", y = "Intensity") +
        facet_wrap(~ label, scales = "free_y", ncol = min(2, max(1, ceiling(length(lst)/2)))) +
        theme_classic(base_size = 13)

      if (isTRUE(input$show_peaks)) {
        pk_list <- lapply(lst, function(z) {
          if (is.null(z$chrom$peaks) || !nrow(z$chrom$peaks)) return(NULL)
          cbind(z$chrom$peaks, label = z$label, stringsAsFactors = FALSE)
        })
        pk_all <- dplyr::bind_rows(pk_list)
        if (!is.null(pk_all) && nrow(pk_all) > 0) {
          # compute per-panel ymax for text placement
          ymax <- df_all %>% dplyr::group_by(label) %>% dplyr::summarise(ymax = max(intensity, na.rm = TRUE))
          pk_all <- dplyr::left_join(pk_all, ymax, by = "label")

          p <- p +
            geom_vline(data = pk_all, aes(xintercept = time), linetype = "dotted", linewidth = 0.6, color = "gray40")

          if (isTRUE(input$label_peaks)) {
            pk_all$lab <- if (isTRUE(input$show_area_pct) && "rel_area_pct" %in% names(pk_all)) {
              sprintf("t=%.2f, %.1f%%", pk_all$time, pk_all$rel_area_pct)
            } else {
              sprintf("t=%.2f", pk_all$time)
            }
            p <- p + geom_text(
              data = pk_all,
              aes(x = time, y = ymax * 1.04, label = lab),
              angle = 0, vjust = 0, size = 3
            )
          }
        }
      }

      p
    }
  })

  # ---- Summary CSV (one row per file) ----
  output$download_summary <- downloadHandler(
    filename = function() {
      paste0("chromatogram_summary_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv")
    },
    content = function(file) {
      lst <- chrom_list()
      validate(need(length(lst) > 0, "No data to summarize."))

      # Build per-file summaries
      summaries <- lapply(lst, function(z) {
        fn <- z$filename
        nick <- z$label
        pk <- z$chrom$peaks
        if (is.null(pk) || !nrow(pk)) {
          tibble::tibble(filename = fn, nickname = nick, total_peaks = 0)
        } else {
          tibble::tibble(
            filename = fn,
            nickname = nick,
            total_peaks = nrow(pk),
            peak_time = pk$time,
            peak_rel_area = pk$rel_area_pct
          )
        }
      })

      # Long to wide: one row per file, add peak_1_time, peak_1_rel_area, ...
      long <- dplyr::bind_rows(summaries, .id = "row_id") %>%
        dplyr::group_by(filename, nickname, total_peaks) %>%
        dplyr::mutate(k = dplyr::row_number()) %>%
        dplyr::ungroup()

      wide_time <- long %>%
        dplyr::select(filename, nickname, total_peaks, k, peak_time) %>%
        tidyr::pivot_wider(names_from = k, values_from = peak_time, names_prefix = "peak_", names_sep = "") %>%
        dplyr::rename_with(~str_replace(.x, "peak_(\\d+)$", "peak_\\1_time"))

      wide_area <- long %>%
        dplyr::select(filename, nickname, total_peaks, k, peak_rel_area) %>%
        tidyr::pivot_wider(names_from = k, values_from = peak_rel_area, names_prefix = "peak_", names_sep = "") %>%
        dplyr::rename_with(~str_replace(.x, "peak_(\\d+)$", "peak_\\1_rel_area"))

      out <- dplyr::full_join(wide_time, wide_area, by = c("filename","nickname","total_peaks"))

      # order columns: filename, nickname, total_peaks, then paired time/area columns by index
      time_cols <- grep("^peak_\\d+_time$", names(out), value = TRUE)
      area_cols <- sub("_time$", "_rel_area", time_cols)
      keep_cols <- c("filename","nickname","total_peaks", as.vector(rbind(time_cols, area_cols)))

      # ensure missing columns are created as NA if some files have fewer peaks
      missing_cols <- setdiff(keep_cols, names(out))
      if (length(missing_cols)) out[missing_cols] <- NA_real_

      out <- out[, keep_cols, drop = FALSE]
      readr::write_csv(out, file, na = "NA")
    }
  )
}

shinyApp(ui, server)
