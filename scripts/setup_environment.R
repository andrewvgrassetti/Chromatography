# setup_environment.R
# ------------------------------------------
# Purpose: Ensure the R environment is reproducible across systems.
# It installs renv if needed and restores all packages defined in renv.lock.

if (!requireNamespace("renv", quietly = TRUE)) {
  message("Installing 'renv'...")
  install.packages("renv", repos = "https://cloud.r-project.org")
}

message("Restoring the R environment with renv...")
renv::restore()
message("✅ Environment restored. All required packages are installed and ready.")
