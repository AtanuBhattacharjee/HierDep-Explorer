## Run once from the repository folder, before the first push and after any package update.
## It installs the packages the app uses and writes manifest.json, which Posit Connect Cloud
## reads to rebuild the environment from GitHub.

pkgs <- c("shiny", "bslib", "ggplot2", "DT", "lme4", "lmerTest", "nlme", "pbkrtest", "Matrix", "rsconnect")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)

## check that the app starts
shiny::runApp(".", launch.browser = TRUE)

## write the manifest (run after closing the app)
rsconnect::writeManifest(appDir = ".", appPrimaryDoc = "app.R")
