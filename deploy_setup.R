## Rewrites manifest.json after a package update. Run from the repository folder,
## then commit and push manifest.json; Posit Connect Cloud republishes from GitHub.

pkgs <- c("shiny", "bslib", "ggplot2", "DT", "lme4", "lmerTest", "nlme", "pbkrtest", "Matrix", "rsconnect")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)

rsconnect::writeManifest(appDir = ".", appFiles = "app.R")
