# app.R — Interactive chromatogram explorer (raw data only, no smoothing)

# ---- Packages ----
library(shiny)
library(bslib)
library(readr)
library(ggplot2)
library(colourpicker)
library(scales)

suppressWarnings({
  if (!requireNamespace("ragg", quietly = TRUE)) {
    message("Tip: install.packages('ragg') for high-quality TIFF/PNG exports.")
  }
})

# ---- Utilities ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- Plot Builder ----
build_plot <- function(chrom, show_fit = TRUE, show_peaks = TRUE, label_peaks = TRUE,
                       palette = NULL, show_area_pct = TRUE) {

  df <- data.frame(time = chrom$time, intensity = chrom$intensity)

  p <- ggplot(df, aes(x = time, y = intensity)) +
    geom_line(linewidth = 0.8, color = "#2c3e50") +
    labs(title = "Chromatogram", x = "Time (min)", y = "Absorbance") +
    theme_classic(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white"),
      axis.line        = element_line(color = "black")
    ) +
    expand_limits(y = max(df$intensity, na.rm = TRUE) * 1.10)

  # ---- Gaussian fit (optional) ----
  if (isTRUE(show_fit) && !is.null(chrom$fit_params)) {
    pars <- tryCatch(unlist(chrom$fit_params)[c("a", "b", "c")], error = function(e) NULL)
    if (!is.null(pars) && all(is.finite(pars))) {
      df$fit <- pars["a"] * exp(-((df$time - pars["b"])^2) / (2 * pars["c"]^2))
      p <- p + geom_line(
        data = df,
        aes(y = fit),
        linetype = "dashed",
        linewidth = 0.8,
        color = "darkgreen"
      )
    }
  }

  # ---- Peaks ----
  if (isTRUE(show_peaks) && !is.null(chrom$peaks) && nrow(chrom$peaks) > 0) {
    peaks_df <- chrom$peaks[!is.na(chrom$peaks$time), , drop = FALSE]
    peaks_df$peak_id <- factor(seq_len(nrow(peaks_df)))
    pal <- palette %||% hue_pal()(nrow(peaks_df))

    p <- p +
      geom_vline(
        data = peaks_df,
        aes(xintercept = time, color = peak_id),
        linetype = "dotted", linewidth = 0.8, alpha = 0.9
      ) +
      scale_color_manual(values = pal, guide = "none")

    if (isTRUE(show_area_pct) && "area" %in% names(peaks_df)) {
      total_area <- sum(peaks_df$area, na.rm = TRUE)
      peaks_df$area_pct <- 100 * peaks_df$area / total_area
    }

    if (isTRUE(label_peaks)) {
      label_df <- peaks_df
      label_df$label <- if (isTRUE(show_area_pct) && "area_pct" %in% names(label_df)) {
        sprintf("t=%.2f\n%.1f%%", label_df$time, label_df$area_pct)
      } else {
        sprintf("t=%.2f", label_df$time)
      }

      p <- p + geom_text(
        data = label_df,
        aes(
          x = time,
          y = max(df$intensity, na.rm = TRUE) * 1.05,
          label = label,
          color = peak_id
        ),
        angle = 90, vjust = -0.4, size = 3, show.legend = FALSE
      )
    }
  }

  p
}

# ---- Export Helper ----
save_plot <- function(p, file, fmt = c("png","tiff","eps"), width = 7, height = 4, dpi = 150) {
  fmt <- match.arg(fmt)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (fmt == "png") {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(file, width = width, height = height, units = "in", res = dpi); print(p); dev.off()
    } else {
      ggsave(file, p, width = width, height = height, dpi = dpi)
    }
  } else if (fmt == "tiff") {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_tiff(file, width = width, height = height, units = "in", res = dpi); print(p); dev.off()
    } else {
      tiff(file, width = width, height = height, units = "in", res = dpi, compression = "lzw"); print(p); dev.off()
    }
  } else { # eps
    grDevices::cairo_ps(file, width = width, height = height, fallback_resolution = dpi); print(p); dev.off()
  }
}

# ---- Load Chromatogram R6 ----
source("R/Chromatogram.R")

# ---- UI ----
ui <- page_fillable(
  theme = bs_theme(bootswatch = "flatly"),
  layout_sidebar(
    sidebar = sidebar(
      h4("Data"),
      fileInput("csv", "Upload CSV (time, intensity — no header needed)", accept = c(".csv")),
      checkboxInput("use_demo", "Use demo data if no file", TRUE),
      hr(),
      h4("Peak Picking"),
      sliderInput("min_height_frac", "Min height (% of max)", min = 0.01, max = 1, value = 0.15, step = 0.01),
      sliderInput("min_prom_frac",   "Min prominence (% of max)", min = 0, max = 1, value = 0.10, step = 0.01),
      numericInput("min_dist_min",   "Min distance (minutes)", value = 0.25, min = 0, step = 0.05),
      numericInput("min_width_min",  "Min width (minutes)", value = 0.10, min = 0, step = 0.05),
      checkboxInput("show_peaks", "Show peaks", TRUE),
      checkboxInput("label_peaks", "Label peaks (t, % area)", TRUE),
      checkboxInput("show_area_pct", "Show area percentages", TRUE),
      checkboxInput("show_fit", "Show fitted Gaussian curve", TRUE),
      hr(),
      h4("Colors"),
      colourInput("raw_col", "Raw", value = "#2c3e50"),
      selectInput("palette", "Peak palette",
                  choices = c("Auto (hue)" = "auto",
                              "Viridis" = "viridis",
                              "Warm" = "warm",
                              "Cool" = "cool"),
                  selected = "auto"),
      hr(),
      h4("Add to Gallery / Export"),
      actionButton("add_plot", "Add current plot to gallery", class = "btn-primary"),
      actionButton("clear_plots", "Clear gallery", class = "btn-warning"),
      br(), br(),
      radioButtons("fmt", "Export format", c("PNG"="png","TIFF"="tiff","EPS"="eps"), inline = TRUE),
      numericInput("dpi", "DPI", value = 150, min = 72, step = 1),
      numericInput("w", "Width (in)", value = 7, min = 3, step = 0.5),
      numericInput("h", "Height (in)", value = 4, min = 2, step = 0.5),
      downloadButton("download_all", "Export gallery as files")
    ),
    card(
      card_header("Current View"),
      plotOutput("main_plot", height = 420),
      div(style="margin-top:6px;color:#777", "Tip: click 'Add to gallery' to snapshot the current configuration.")
    ),
    card(
      card_header("Plot Gallery"),
      uiOutput("gallery_ui")
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {

  # ---- Reactive data (CSV or demo) ----
  data_reactive <- reactive({
    req(input$use_demo || !is.null(input$csv))
    if (!is.null(input$csv)) {
      suppressWarnings({
        df <- tryCatch(readr::read_csv(input$csv$datapath, col_names = FALSE, show_col_types = FALSE),
                       error = function(e) NULL)
      })
      if (is.null(df) || ncol(df) < 2) {
        showNotification("CSV must have at least two numeric columns (time, intensity).", type = "error")
        return(NULL)
      }
      df <- df[, 1:2]
      names(df) <- c("time", "intensity")
      df <- df[complete.cases(df), ]
      df$time <- as.numeric(df$time)
      df$intensity <- as.numeric(df$intensity)
      df <- df[!is.na(df$time) & !is.na(df$intensity), ]
      return(df)
    }

    # Demo data
    set.seed(123)
    n <- 800
    time <- seq(0, 20, length.out = n)
    peak1 <- 110 * exp(-((time - 9.8)^2) / (2 * 0.25^2))
    peak2 <- 90  * exp(-((time - 10.6)^2) / (2 * 0.22^2))
    intensity <- peak1 + peak2 + rnorm(n, sd = 4)
    data.frame(time, intensity)
  })

  # ---- Chromatogram Object ----
  chrom_reactive <- reactive({
    df <- data_reactive()
    obj <- Chromatogram$new(df$time, df$intensity)
    suppressWarnings(try(obj$fit_gaussian(), silent = TRUE))
    min_h <- input$min_height_frac * max(obj$intensity, na.rm = TRUE)
    obj$find_peaks(min_height = min_h, min_distance = input$min_dist_min)
    obj$integrate_peaks()
    obj
  })

  # ---- Palette helper ----
  compute_palette <- function(n) {
    if (n <= 0) return(NULL)
    switch(input$palette,
           "viridis" = viridis_pal()(n),
           "warm"    = hue_pal(h = c(0,90))(n),
           "cool"    = hue_pal(h = c(180,270))(n),
           hue_pal()(n))
  }

  # ---- Main plot ----
  output$main_plot <- renderPlot({
    chrom <- chrom_reactive()
    pal <- compute_palette(nrow(chrom$peaks %||% data.frame()))
    p <- build_plot(
      chrom,
      show_fit = isTRUE(input$show_fit),
      show_peaks = isTRUE(input$show_peaks),
      label_peaks = isTRUE(input$label_peaks),
      palette = pal,
      show_area_pct = isTRUE(input$show_area_pct)
    ) + theme_minimal(base_size = 12)
    p
  })

  # ---- Gallery management ----
  rv <- reactiveValues(gallery = list())

  observeEvent(input$add_plot, {
    chrom <- chrom_reactive()
    cfg <- list(
      chrom = chrom,
      show_fit = isTRUE(input$show_fit),
      show_peaks = isTRUE(input$show_peaks),
      label_peaks = isTRUE(input$label_peaks),
      show_area_pct = isTRUE(input$show_area_pct),
      palette = compute_palette(nrow(chrom$peaks %||% data.frame()))
    )
    rv$gallery <- append(rv$gallery, list(cfg))
    showNotification("Added plot to gallery.", type = "message")
  })

  observeEvent(input$clear_plots, {
    rv$gallery <- list()
    showNotification("Gallery cleared.", type = "warning")
  })

  output$gallery_ui <- renderUI({
    if (length(rv$gallery) == 0) return(div("No plots yet."))
    tagList(
      lapply(seq_along(rv$gallery), function(i) {
        card(
          card_header(paste("Plot", i)),
          plotOutput(paste0("plot_", i), height = 300),
          fluidRow(
            column(6, downloadButton(paste0("dl_", i, "_png"),  "PNG")),
            column(6, downloadButton(paste0("dl_", i, "_tiff"), "TIFF"))
          ),
          fluidRow(
            column(6, downloadButton(paste0("dl_", i, "_eps"),  "EPS"))
          )
        )
      })
    )
  })

  # ---- Render each gallery plot + downloads ----
  observe({
    lapply(seq_along(rv$gallery), function(i) {
      local({
        idx <- i
        cfg <- rv$gallery[[idx]]
        output[[paste0("plot_", idx)]] <- renderPlot({
          p <- build_plot(cfg$chrom, cfg$show_fit, cfg$show_peaks, cfg$label_peaks,
                          cfg$palette, cfg$show_area_pct) +
            theme_minimal(base_size = 12)
          p
        })
        for (fmt in c("png", "tiff", "eps")) {
          output[[paste0("dl_", idx, "_", fmt)]] <- downloadHandler(
            filename = function() sprintf("plot_%02d.%s", idx, fmt),
            content = function(file) {
              p <- build_plot(cfg$chrom, cfg$show_fit, cfg$show_peaks, cfg$label_peaks,
                              cfg$palette, cfg$show_area_pct) +
                theme_minimal(base_size = 12)
              save_plot(p, file, fmt, width = input$w, height = input$h, dpi = input$dpi)
            }
          )
        }
      })
    })
  })
}

shinyApp(ui, server)
