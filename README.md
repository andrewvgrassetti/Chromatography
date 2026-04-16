# Chromatography Peak Analysis App

An R Shiny app for chromatography peak detection, integration, and visualization from `(time, intensity)` CSV files.

## What it does

- Detects peaks from smoothed signals
- Integrates peak areas (trapezoidal rule)
- Calculates relative area percentages
- Plots multiple chromatograms (overlay or separate panels)
- Exports:
  - Peak summary CSV
  - Plot as PNG, TIFF, or EPS

## Repository structure

```text
Chromatography/
├── app.R                    # Shiny app entrypoint
├── R/
│   └── Chromatogram.R       # Core R6 analysis class
├── scripts/
│   ├── setup_environment.R  # renv-based dependency restore
│   └── test_chrom.R         # Basic analysis test script
├── renv.lock                # Locked dependency snapshot
└── README.md
```

## Prerequisites

- R installed locally
- Internet access for first-time package installation

## Setup

From the repository root:

```bash
cd /path/to/Chromatography
```

### Recommended (reproducible) setup with `renv`

```bash
Rscript scripts/setup_environment.R
```

### Alternative: install packages manually in R

```r
install.packages(c(
  "shiny", "bslib", "readr", "ggplot2", "colourpicker", "scales",
  "dplyr", "purrr", "tidyr", "stringr",
  "R6", "signal", "minpack.lm", "pracma", "viridisLite"
))

# Optional for higher-quality exports:
install.packages(c("ragg", "Cairo"))
```

## Run the app

```bash
Rscript -e "shiny::runApp('.')"
```

Or in an R session:

```r
library(shiny)
shiny::runApp(".")
```

If a browser does not open automatically, use the local URL printed in the console.

## Input file format

Upload one or more CSV files with exactly two numeric columns:

1. `time`
2. `intensity`

Example:

```text
0.00,12.3
0.01,12.9
0.02,13.1
```

Notes:
- No header is required
- Time units can be minutes (or any consistent unit)
- Non-numeric/incomplete rows are ignored

## Testing

Run the included test script from the repository root:

```bash
Rscript scripts/test_chrom.R
```

The script exercises smoothing, fitting, baseline/area calculations, and peak analysis with synthetic data.

## Troubleshooting

- **`there is no package called ...`**  
  Run `Rscript scripts/setup_environment.R` (or install missing packages manually).
- **Plot export quality is poor**  
  Install optional packages `ragg` and `Cairo`.
- **No peaks detected**  
  Lower **Min height (% of max)** in the app sidebar.
