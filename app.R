## =============================================================================
##  HierDep Explorer
##  Companion Shiny application to:
##  "Longitudinal Mental Health and Neurodevelopmental Studies: Guide to Dealing
##   with Non-independent Data Using Hierarchical Dependence"
##  Atanu Bhattacharjee, IoPPN, King's College London
##
##  Layers of dependence: site -> measure -> time.
##  Data sources: simulated ABIDE-format demo, the public ABIDE I phenotypic file,
##  or a user CSV (one row per participant, or per participant-visit).
## =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
  library(DT)
  library(lme4)
  library(lmerTest)
  library(nlme)
  library(pbkrtest)
  library(Matrix)
})

options(shiny.maxRequestSize = 50 * 1024^2)

## plots rendered at screen resolution so text stays legible
renderPlot <- function(expr, ..., res = 110) {
  e <- substitute(expr)
  shiny::renderPlot(e, ..., res = res, env = parent.frame(), quoted = TRUE)
}

## -----------------------------------------------------------------------------
##  Constants
## -----------------------------------------------------------------------------

ABIDE_URL <- paste0("https://s3.amazonaws.com/fcp-indi/data/Projects/",
                    "ABIDE_Initiative/Phenotypic_V1_0b_preprocessed1.csv")

ABIDE_MEASURES <- data.frame(
  var = c("FIQ", "VIQ", "PIQ", "ADOS_total", "ADOS_comm", "ADOS_social",
          "ADOS_stereo", "ADIR_social", "ADIR_verbal", "ADIR_rrb",
          "log_fd", "log_pfd", "dvars"),
  label = c("Full-scale IQ", "Verbal IQ", "Performance IQ", "ADOS total",
            "ADOS communication", "ADOS social", "ADOS stereotyped behaviour",
            "ADI-R social", "ADI-R verbal", "ADI-R RRB", "Mean FD (log)",
            "% volumes FD > 0.2 mm (log1p)", "Functional DVARS"),
  block = c(rep("Cognitive", 3), rep("Clinician-observed (ADOS)", 4),
            rep("Parent interview (ADI-R)", 3), rep("MRI acquisition", 3)),
  both = c(rep(TRUE, 3), rep(FALSE, 7), rep(TRUE, 3)),
  stringsAsFactors = FALSE)

BLOCK_COL <- c("Cognitive" = "#1B9E77",
               "Clinician-observed (ADOS)" = "#D95F02",
               "Parent interview (ADI-R)" = "#7570B3",
               "MRI acquisition" = "#E7298A",
               "Uploaded" = "#1F78B4",
               "LEAP (Table 2)" = "#E6AB02")

MODEL_COL <- c("Naive (no site)" = "#E41A1C",
               "Site as fixed effect" = "#FF7F00",
               "Random site intercept" = "#377EB8",
               "Random intercept + slope" = "#4DAF4A")

## Table 2 reference values (ABIDE I and LEAP) used as planning presets
TABLE2 <- data.frame(
  measure = c("Full-scale IQ", "Verbal IQ", "Performance IQ", "ADOS total",
              "ADOS communication", "ADOS social", "ADOS stereotyped behaviour",
              "ADI-R social", "ADI-R verbal", "ADI-R RRB", "Mean FD (log)",
              "% volumes FD > 0.2 mm (log1p)", "Functional DVARS",
              "False Belief score", "RMET reaction time",
              "SWM between-search errors", "SWM total errors",
              "RMET % correct", "SWM strategy"),
  dataset = c(rep("ABIDE I", 13), rep("LEAP", 6)),
  N = c(813, 723, 737, 279, 259, 260, 224, 278, 279, 279, 873, 873, 873,
        648, 636, 714, 714, 636, 714),
  sites = c(19, 16, 17, 16, 14, 15, 11, 16, 16, 16, 20, 20, 20, 6, 6, 6, 6, 6, 6),
  icc = c(0.050, 0.048, 0.070, 0.100, 0.115, 0.076, 0.101, 0.000, 0.056, 0.096,
          0.190, 0.205, 0.569, 0.035, 0.031, 0.028, 0.019, 0.012, 0.001),
  stringsAsFactors = FALSE)
TABLE2$nbar <- TABLE2$N / TABLE2$sites

## LEAP: reported between-centre eta^2 recovered from the chance-corrected ICC
eta2_from_icc <- function(icc, N, k) {
  nbar <- N / k
  F <- 1 + icc * nbar / (1 - icc)
  F * (k - 1) / (F * (k - 1) + (N - k))
}
LEAP_DEFAULT <- with(subset(TABLE2, dataset == "LEAP"),
                     paste(measure, sites, N,
                           sprintf("%.5f", eta2_from_icc(icc, N, sites)),
                           sep = ", ", collapse = "\n"))

## -----------------------------------------------------------------------------
##  Helpers
## -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a
num <- function(x) {
  x <- suppressWarnings(as.numeric(as.character(x)))
  x[!is.na(x) & x == -9999] <- NA
  x
}
zs <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
f3 <- function(x) { x <- suppressWarnings(as.numeric(x)); ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 3)) }
f2 <- function(x) { x <- suppressWarnings(as.numeric(x)); ifelse(is.na(x), "NA", formatC(x, format = "f", digits = 2)) }
fpt <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format = "f", digits = 3)))
fp <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))

icc_of <- function(m) {
  v <- as.data.frame(VarCorr(m))
  s <- v$vcov[v$grp == "site"]
  s / sum(v$vcov)
}

## ANOVA (F-based) ICC and 95% interval with the effective cluster size n0
icc_F <- function(y, g, nuse = NULL) {
  g <- droplevels(factor(g)); k <- nlevels(g); N <- length(y)
  ns <- table(g); n0 <- (N - sum(ns^2) / N) / (k - 1)
  if (is.null(nuse)) nuse <- n0
  a <- anova(lm(y ~ g)); F <- a[1, 3] / a[2, 3]
  FL <- F / qf(0.975, k - 1, N - k); FU <- F / qf(0.025, k - 1, N - k)
  icc <- function(f) max(0, (f - 1) / (f - 1 + nuse))
  c(icc = icc(F), lo = icc(FL), hi = icc(FU), n0 = n0, F = F)
}

theme_hd <- function() {
  theme_minimal(base_size = 14) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#E3E8F0"),
          axis.text = element_text(colour = "#1E293B", size = 12.5),
          axis.title = element_text(colour = "#0F172A", size = 13.5),
          legend.position = "bottom",
          legend.text = element_text(colour = "#1E293B", size = 12.5),
          legend.title = element_text(colour = "#0F172A", size = 12.5),
          plot.subtitle = element_text(colour = "#334155", size = 12.5),
          plot.title = element_text(face = "bold", size = 14, colour = "#1E3A5F"),
          strip.text = element_text(face = "bold", colour = "#1E3A5F"),
          strip.background = element_rect(fill = "#F1F5FB", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA))
}

## -----------------------------------------------------------------------------
##  Simulated ABIDE-format phenotypic file (demo)
##  Site sizes, groups, sex, ages and eye status follow Table 1; site ICCs follow
##  Table 2; between- and within-site correlations follow Section 4.
## -----------------------------------------------------------------------------

ABIDE_SITES <- data.frame(
  site = c("NYU", "UM_1", "USM", "UCLA_1", "PITT", "MAX_MUN", "TRINITY", "YALE",
           "UM_2", "KKI", "LEUVEN_1", "LEUVEN_2", "OLIN", "SDSU", "SBL", "OHSU",
           "STANFORD", "UCLA_2", "CALTECH", "CMU"),
  n = c(172, 88, 67, 64, 50, 46, 44, 41, 34, 33, 28, 28, 28, 27, 26, 25, 25, 21, 15, 11),
  aut = c(74, 36, 43, 37, 24, 19, 19, 22, 13, 12, 14, 12, 14, 8, 12, 12, 12, 11, 5, 6),
  fem = c(36, 25, 0, 9, 7, 4, 0, 16, 2, 9, 0, 7, 5, 6, 0, 0, 7, 2, 5, 4),
  amed = c(13.7, 13.7, 19.8, 13.5, 16.6, 26.5, 16.4, 13.8, 15.3, 10.2, 22, 14.2, 17,
           14.1, 33.5, 10.5, 9.3, 12.2, 21.2, 27),
  alo = c(6.5, 8.2, 8.8, 8.4, 9.3, 7, 12, 7, 12.8, 8.2, 18, 12.1, 10, 8.7, 20, 8,
          7.5, 9.8, 17, 19),
  ahi = c(39.1, 19.2, 50.2, 17.9, 35.2, 58, 25.7, 17.8, 28.8, 12.8, 32, 16.9, 24,
          17.1, 49, 15.2, 12.9, 16.5, 56.2, 33),
  closed = c(25, 0, 0, 0, 50, 38, 44, 0, 0, 0, 0, 28, 0, 0, 26, 0, 25, 0, 15, 11),
  stringsAsFactors = FALSE)

sim_abide_raw <- function(seed = 2026) {
  set.seed(seed)
  vars <- c("FIQ", "VIQ", "PIQ", "ADOS_TOTAL", "ADOS_COMM", "ADOS_SOCIAL",
            "ADOS_STEREO_BEHAV", "ADI_R_SOCIAL_TOTAL_A", "ADI_R_VERBAL_TOTAL_BV",
            "ADI_RRB_TOTAL_C", "fd", "pfd", "dvars")
  M <- length(vars)
  icc <- c(0.01, 0.01, 0.02, 0.10, 0.115, 0.076, 0.101, 0.001, 0.056, 0.096,
           0.19, 0.205, 0.57)
  pair <- function(R, i, j, r) { R[i, j] <- R[j, i] <- r; R }
  Re <- diag(M)
  for (p in list(c(1, 2, .80), c(1, 3, .80), c(2, 3, .55), c(4, 5, .80),
                 c(4, 6, .90), c(4, 7, .50), c(5, 6, .60), c(5, 7, .30),
                 c(6, 7, .30), c(8, 9, .50), c(8, 10, .40), c(9, 10, .35),
                 c(4, 8, .11), c(6, 8, .10), c(11, 12, .85), c(11, 13, -.09),
                 c(12, 13, -.07), c(1, 11, -.10)))
    Re <- pair(Re, p[1], p[2], p[3])
  Ru <- diag(M)
  for (p in list(c(1, 2, .6), c(1, 3, .6), c(2, 3, .4), c(4, 5, .7), c(4, 6, .8),
                 c(4, 7, .5), c(5, 6, .6), c(8, 9, .5), c(8, 10, .4), c(4, 8, .57),
                 c(11, 12, .9), c(11, 13, .21), c(12, 13, .2)))
    Ru <- pair(Ru, p[1], p[2], p[3])
  Re <- as.matrix(nearPD(Re, corr = TRUE)$mat)
  Ru <- as.matrix(nearPD(Ru, corr = TRUE)$mat)
  Le <- chol(diag(sqrt(1 - icc)) %*% Re %*% diag(sqrt(1 - icc)))
  Lu <- chol(diag(sqrt(icc)) %*% Ru %*% diag(sqrt(icc)))

  instr <- c("WASI", "WISC_IV", "WAIS_III", "DAS_II", "HAWIK_IV", "WISC",
             "STANFORD_BINET", "RAVENS")
  w_instr <- setNames(rnorm(length(instr), 0, 0.2), instr)
  prim <- setNames(sample(instr, nrow(ABIDE_SITES), replace = TRUE,
                          prob = c(.35, .2, .15, .1, .05, .05, .05, .05)),
                   ABIDE_SITES$site)
  second <- sample(ABIDE_SITES$site, 6)
  beta <- c(-0.50, -0.55, -0.35, rep(0, 7), 0.20, 0.20, 0)
  tau <- c(0.41, 0.47, 0.10, rep(0, 7), 0.33, 0.05, 0.05)

  rows <- list()
  for (s in seq_len(nrow(ABIDE_SITES))) {
    S <- ABIDE_SITES[s, ]; n <- S$n
    dx <- sample(c(rep(1, S$aut), rep(2, n - S$aut)))
    sex <- sample(c(rep(2, S$fem), rep(1, n - S$fem)))
    lo <- runif(ceiling(n / 2), S$alo, S$amed); hi <- runif(n - ceiling(n / 2), S$amed, S$ahi)
    age <- sample(c(lo, hi))
    eye <- sample(c(rep(2, S$closed), rep(1, n - S$closed)))
    test <- rep(prim[S$site], n)
    if (S$site %in% second) {
      k <- sample(n, round(0.3 * n)); test[k] <- sample(setdiff(instr, prim[S$site]), 1)
    }
    u <- as.vector(rnorm(M) %*% Lu)
    v <- rnorm(M, 0, tau)
    E <- matrix(rnorm(n * M), n) %*% Le
    aut <- as.numeric(dx == 1)
    Z <- sweep(E, 2, u, "+") + outer(aut, beta + v)
    Z[, 1:3] <- Z[, 1:3] + w_instr[test]
    Z[, 11] <- Z[, 11] - 0.015 * (age - 17)
    raw <- data.frame(
      SITE_ID = S$site, SUB_ID = 50000 + 1000 * s + seq_len(n),
      FILE_ID = paste0(S$site, "_", seq_len(n)), DX_GROUP = dx,
      AGE_AT_SCAN = round(age, 2), SEX = sex,
      FIQ = round(108 + 15 * Z[, 1]), VIQ = round(107 + 15 * Z[, 2]),
      PIQ = round(106 + 15 * Z[, 3]), FIQ_TEST_TYPE = test,
      ADOS_TOTAL = pmax(0, round(11 + 3.8 * Z[, 4])),
      ADOS_COMM = pmax(0, round(3.8 + 1.6 * Z[, 5])),
      ADOS_SOCIAL = pmax(0, round(7.5 + 2.6 * Z[, 6])),
      ADOS_STEREO_BEHAV = pmax(0, round(2.3 + 1.5 * Z[, 7])),
      ADOS_MODULE = ifelse(age < 16, 3, 4),
      ADI_R_SOCIAL_TOTAL_A = pmax(0, round(19 + 5.5 * Z[, 8])),
      ADI_R_VERBAL_TOTAL_BV = pmax(0, round(15.5 + 4.4 * Z[, 9])),
      ADI_RRB_TOTAL_C = pmax(0, round(5.5 + 2.5 * Z[, 10])),
      EYE_STATUS_AT_SCAN = eye,
      func_mean_fd = round(exp(-2.2 + 0.6 * Z[, 11]), 4),
      func_perc_fd = round(pmax(0, expm1(2.3 + 1.1 * Z[, 12])), 3),
      func_dvars = round(1.4 + 0.3 * Z[, 13], 4),
      qc_rater_1 = "OK", qc_anat_rater_2 = "OK", qc_func_rater_2 = "OK",
      qc_anat_rater_3 = "OK", qc_func_rater_3 = "OK",
      stringsAsFactors = FALSE)
    clin <- c("ADOS_TOTAL", "ADOS_COMM", "ADOS_SOCIAL", "ADOS_STEREO_BEHAV",
              "ADOS_MODULE", "ADI_R_SOCIAL_TOTAL_A", "ADI_R_VERBAL_TOTAL_BV",
              "ADI_RRB_TOTAL_C")
    raw[raw$DX_GROUP == 2, clin] <- -9999
    rows[[s]] <- raw
  }
  r <- do.call(rbind, rows)
  miss <- function(cols, sites) r[r$SITE_ID %in% sites, cols] <<- -9999
  ados <- c("ADOS_TOTAL", "ADOS_COMM", "ADOS_SOCIAL", "ADOS_STEREO_BEHAV", "ADOS_MODULE")
  miss(ados, c("UM_1", "UM_2", "LEUVEN_1", "LEUVEN_2"))
  miss("ADOS_COMM", c("CMU", "CALTECH"))
  miss("ADOS_SOCIAL", "CMU")
  miss("ADOS_STEREO_BEHAV", c("CMU", "CALTECH", "SBL", "KKI", "OHSU"))
  miss(c("VIQ", "PIQ"), c("KKI", "OLIN", "OHSU"))
  miss("VIQ", "CMU")
  miss("FIQ", "SBL")
  miss(c("ADI_R_SOCIAL_TOTAL_A", "ADI_R_VERBAL_TOTAL_BV", "ADI_RRB_TOTAL_C"),
       c("CMU", "CALTECH", "SBL", "MAX_MUN"))
  for (cc in c("FIQ", "VIQ", "PIQ", "ADI_R_SOCIAL_TOTAL_A", "ADOS_TOTAL")) {
    k <- sample(nrow(r), round(0.03 * nrow(r))); r[k, cc] <- -9999
  }
  r$FIQ_TEST_TYPE[r$FIQ == -9999] <- ""
  ## participants excluded by quality control or without a functional file
  ex <- r[sample(nrow(r), 60), ]
  ex$SUB_ID <- 90000 + seq_len(60)
  ex$qc_func_rater_2[1:30] <- "fail"
  ex$qc_rater_1[31:45] <- "fail"
  ex$FILE_ID[46:60] <- "no_filename"
  rbind(r, ex)
}

## -----------------------------------------------------------------------------
##  Data preparation
## -----------------------------------------------------------------------------

is_abide <- function(df) all(c("SITE_ID", "SUB_ID", "DX_GROUP", "func_dvars") %in% names(df))

prep_abide <- function(ph, source = "abide") {
  n_in <- nrow(ph)
  if ("FILE_ID" %in% names(ph)) ph <- ph[ph$FILE_ID != "no_filename", ]
  for (q in intersect(c("qc_rater_1", "qc_anat_rater_2", "qc_func_rater_2",
                        "qc_anat_rater_3", "qc_func_rater_3"), names(ph)))
    ph <- ph[is.na(ph[[q]]) | ph[[q]] != "fail", ]
  d <- data.frame(
    id = as.character(ph$SUB_ID), site = factor(ph$SITE_ID),
    dx = ifelse(ph$DX_GROUP == 1, 1, 0),
    age = num(ph$AGE_AT_SCAN), female = ifelse(ph$SEX == 2, 1, 0),
    FIQ = num(ph$FIQ), VIQ = num(ph$VIQ), PIQ = num(ph$PIQ),
    iqtest = ifelse(ph$FIQ_TEST_TYPE %in% c("", NA, "-9999"), NA, ph$FIQ_TEST_TYPE),
    ADOS_total = num(ph$ADOS_TOTAL), ADOS_comm = num(ph$ADOS_COMM),
    ADOS_social = num(ph$ADOS_SOCIAL), ADOS_stereo = num(ph$ADOS_STEREO_BEHAV),
    ADOS_module = num(ph$ADOS_MODULE),
    ADIR_social = num(ph$ADI_R_SOCIAL_TOTAL_A), ADIR_verbal = num(ph$ADI_R_VERBAL_TOTAL_BV),
    ADIR_rrb = num(ph$ADI_RRB_TOTAL_C),
    eyes_closed = ifelse(num(ph$EYE_STATUS_AT_SCAN) == 2, 1, 0),
    fd = num(ph$func_mean_fd), pfd = num(ph$func_perc_fd), dvars = num(ph$func_dvars),
    stringsAsFactors = FALSE)
  d$iqtest[d$iqtest %in% c("WISC_IV_FULL", "WISC_IV_4_SUBTESTS")] <- "WISC_IV"
  d$iqtest[d$iqtest %in% c("WAIS")] <- "WAIS_III"
  d$iqtest[d$iqtest %in% c("GIT")] <- NA
  clin <- grep("^ADOS_|^ADIR_", names(d), value = TRUE)
  d[d$dx == 0, clin] <- NA
  d$log_fd <- log(d$fd)
  d$log_pfd <- log1p(d$pfd)
  d$institution <- factor(sub("_[0-9]+$", "", as.character(d$site)))
  list(d = d, d_base = d, meas = ABIDE_MEASURES, source = source,
       groups = c("Autistic", "Non-autistic"), expo = "eyes_closed",
       expo_labels = c(eyes_closed = "Eyes closed at scan"),
       covs = c("dx", "age", "female", "eyes_closed"),
       instrument = "iqtest", longitudinal = FALSE,
       n_in = n_in, n_out = nrow(d))
}

prep_generic <- function(raw, m) {
  nn <- function(x) make.names(x)
  d <- data.frame(id = as.character(raw[[m$id]]), site = factor(raw[[m$site]]),
                  stringsAsFactors = FALSE)
  d$institution <- d$site
  has_dx <- !is.null(m$group) && m$group != "(none)"
  if (has_dx) d$dx <- as.numeric(as.character(raw[[m$group]]) == m$group_level)
  if (!is.null(m$age) && m$age != "(none)") d$age <- num(raw[[m$age]])
  if (!is.null(m$sex) && m$sex != "(none)")
    d$female <- as.numeric(as.character(raw[[m$sex]]) == m$sex_female)
  has_time <- !is.null(m$time) && m$time != "(none)"
  if (has_time) d$time <- num(raw[[m$time]])
  instrument <- NULL
  if (!is.null(m$instr) && m$instr != "(none)") {
    d$iqtest <- as.character(raw[[m$instr]]); instrument <- "iqtest"
  }
  expo <- character(0)
  for (x in m$sitecov) { d[[nn(x)]] <- num(raw[[x]]); expo <- c(expo, nn(x)) }
  for (x in m$outcomes) d[[nn(x)]] <- num(raw[[x]])
  meas <- data.frame(var = nn(m$outcomes), label = m$outcomes, block = "Uploaded",
                     both = has_dx, stringsAsFactors = FALSE)
  d_base <- d
  if (has_time) {
    d <- d[!is.na(d$time), ]
    o <- order(d$id, d$time)
    d_base <- d[o, ][!duplicated(d$id[o]), ]
  }
  covs <- c(if (has_dx) "dx", intersect(c("age", "female"), names(d)), expo)
  list(d = d, d_base = d_base, meas = meas, source = "upload",
       groups = if (has_dx) c(as.character(m$group_level), "Other") else NULL,
       expo = expo, expo_labels = setNames(expo, expo), covs = covs,
       instrument = instrument, longitudinal = has_time,
       n_in = nrow(raw), n_out = nrow(d_base))
}

## -----------------------------------------------------------------------------
##  Analyses
## -----------------------------------------------------------------------------

has_var <- function(dd, v) v %in% names(dd) && sum(!is.na(dd[[v]])) > 0 && sd(dd[[v]], na.rm = TRUE) > 0

## Site layer: empty model (1), adjusted model (3), design effect (4)
site_layer_one <- function(D, k, ci = "F", nsim = 200) {
  d <- D$d_base; meas <- D$meas
  v <- meas$var[k]
  dd <- d[!is.na(d[[v]]), ]
  if (grepl("^ADOS", v) && "ADOS_module" %in% names(dd)) dd <- dd[!is.na(dd$ADOS_module), ]
  dd$site <- droplevels(dd$site)
  nsite <- nlevels(dd$site); N <- nrow(dd); nbar <- N / nsite
  if (nsite < 3) stop("fewer than three sites")
  ns <- table(dd$site); n0 <- (N - sum(ns^2) / N) / (nsite - 1)
  m0 <- lme4::lmer(as.formula(paste(v, "~ 1 + (1|site)")), data = dd, REML = TRUE)
  icc0 <- icc_of(m0)
  if (ci == "boot") {
    bt <- bootMer(m0, icc_of, nsim = nsim, type = "parametric", use.u = FALSE)
    cil <- quantile(bt$t, c(.025, .975), na.rm = TRUE)
  } else {
    cil <- icc_F(dd[[v]], dd$site)[c("lo", "hi")]
  }
  rhs <- character(0)
  if (meas$both[k] && has_var(dd, "dx")) rhs <- c(rhs, "dx")
  if (has_var(dd, "age")) {
    dd$age_b <- ave(dd$age, dd$site, FUN = function(z) mean(z, na.rm = TRUE))
    dd$age_w <- dd$age - dd$age_b
    rhs <- c(rhs, "age_w", "age_b")
  }
  if (has_var(dd, "female")) rhs <- c(rhs, "female")
  if (grepl("^ADOS", v) && "ADOS_module" %in% names(dd) &&
      length(unique(dd$ADOS_module)) > 1) rhs <- c(rhs, "factor(ADOS_module)")
  icc1 <- NA
  if (length(rhs)) {
    m1 <- lme4::lmer(as.formula(paste(v, "~", paste(rhs, collapse = " + "), "+ (1|site)")),
                     data = dd, REML = TRUE)
    icc1 <- icc_of(m1)
  }
  m0ml <- lme4::lmer(as.formula(paste(v, "~ 1 + (1|site)")), data = dd, REML = FALSE)
  lr <- max(0, as.numeric(2 * (logLik(m0ml) - logLik(lm(as.formula(paste(v, "~ 1")), data = dd)))))
  re <- ranef(m0, condVar = TRUE)$site
  pv <- attr(re, "postVar")[1, 1, ]
  sdY <- sd(dd[[v]])
  list(
    res = data.frame(var = v, label = meas$label[k], block = meas$block[k],
                     N = N, sites = nsite, nbar = nbar, n0 = n0, icc = icc0,
                     lo = unname(cil[1]), hi = unname(cil[2]), icc_adj = icc1,
                     deff = 1 + (nbar - 1) * icc0, neff = N / (1 + (nbar - 1) * icc0),
                     p_site = 0.5 * pchisq(lr, 1, lower.tail = FALSE),
                     singular = isSingular(m0), stringsAsFactors = FALSE),
    blup = data.frame(var = v, label = meas$label[k], site = rownames(re),
                      n = as.vector(ns[rownames(re)]),
                      est = re[, 1] / sdY, se = sqrt(pv) / sdY, stringsAsFactors = FALSE))
}

## Group difference under four analyses, random slope (6), prediction interval (16)
dx_four <- function(D, v) {
  d <- D$d_base
  dd <- d[!is.na(d[[v]]) & !is.na(d$dx), ]
  if (has_var(dd, "age")) dd <- dd[!is.na(dd$age), ]
  dd$site <- droplevels(dd$site); S <- nlevels(dd$site); N <- nrow(dd)
  dd$z <- zs(dd[[v]])
  cn <- "dx"; cm <- "dx"
  if (has_var(dd, "age")) {
    dd$age_b <- ave(dd$age, dd$site); dd$age_w <- dd$age - dd$age_b
    cn <- c(cn, "age"); cm <- c(cm, "age_w", "age_b")
  }
  if (has_var(dd, "female")) { cn <- c(cn, "female"); cm <- c(cm, "female") }
  fn <- as.formula(paste("z ~", paste(cn, collapse = " + ")))
  a <- summary(lm(fn, data = dd))$coef["dx", 1:2]
  b <- summary(lm(update(fn, . ~ . + site), data = dd))$coef["dx", 1:2]
  fri <- as.formula(paste("z ~", paste(cm, collapse = " + "), "+ (1|site)"))
  frs <- as.formula(paste("z ~", paste(cm, collapse = " + "), "+ (1 + dx|site)"))
  ctl <- lmerControl(optimizer = "bobyqa")
  mri <- lmerTest::lmer(fri, data = dd)
  c1 <- summary(mri, ddf = "Kenward-Roger")$coef["dx", 1:3]
  mrs <- lmerTest::lmer(frs, data = dd, control = ctl)
  c2 <- summary(mrs, ddf = "Kenward-Roger")$coef["dx", 1:3]
  vc <- as.data.frame(VarCorr(mrs))
  tau <- sqrt(vc$vcov[vc$grp == "site" & vc$var1 == "dx" & is.na(vc$var2)])
  l0 <- lme4::lmer(fri, data = dd, REML = FALSE)
  l1 <- lme4::lmer(frs, data = dd, REML = FALSE, control = ctl)
  chi <- max(0, as.numeric(2 * (logLik(l1) - logLik(l0))))
  p_slope <- 0.5 * pchisq(chi, 1, lower.tail = FALSE) + 0.5 * pchisq(chi, 2, lower.tail = FALSE)
  tab <- data.frame(model = names(MODEL_COL),
                    est = c(a[1], b[1], c1[1], c2[1]), se = c(a[2], b[2], c1[2], c2[2]),
                    df = c(N - length(cn) - 1, N - length(cn) - S, c1[3], c2[3]))
  tab$lo <- tab$est - qt(.975, tab$df) * tab$se
  tab$hi <- tab$est + qt(.975, tab$df) * tab$se
  pi_hw <- qt(.975, S - 2) * sqrt(c2[2]^2 + tau^2)
  list(tab = tab, tau = tau, p_slope = p_slope, N = N, S = S,
       pi = c(c2[1] - pi_hw, c2[1] + pi_hw), singular = isSingular(mrs))
}

## Site-specific differences with DerSimonian-Laird pooling
persite_pool <- function(D, v, minn = 3) {
  d <- D$d_base
  dd <- d[!is.na(d[[v]]) & !is.na(d$dx), ]; dd$z <- zs(dd[[v]])
  ps <- do.call(rbind, lapply(levels(droplevels(dd$site)), function(s) {
    x <- dd[dd$site == s, ]
    if (sum(x$dx == 1) >= minn && sum(x$dx == 0) >= minn) {
      f <- summary(lm(z ~ dx, data = x))$coef
      data.frame(site = s, n = nrow(x), est = f["dx", 1], se = f["dx", 2])
    }
  }))
  if (is.null(ps) || nrow(ps) < 3) return(NULL)
  w <- 1 / ps$se^2; fe <- sum(w * ps$est) / sum(w)
  Q <- sum(w * (ps$est - fe)^2); k <- nrow(ps)
  t2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (ps$se^2 + t2); re <- sum(ws * ps$est) / sum(ws); se <- sqrt(1 / sum(ws))
  list(ps = ps, pool = data.frame(k = k, re = re, re_se = se, tau = sqrt(t2),
                                  I2 = max(0, (Q - (k - 1)) / Q), Q = Q,
                                  pQ = pchisq(Q, k - 1, lower.tail = FALSE),
                                  pi_lo = re - qt(.975, k - 2) * sqrt(t2 + se^2),
                                  pi_hi = re + qt(.975, k - 2) * sqrt(t2 + se^2)))
}

## Intraclass correlation of a covariate (rho_x)
rho_x <- function(d, x) {
  dd <- d[!is.na(d[[x]]), ]
  if (!has_var(dd, x)) return(NA)
  dd$site <- droplevels(dd$site)
  if (all(tapply(dd[[x]], dd$site, function(z) length(unique(z)) == 1))) return(1)
  tryCatch({
    m <- suppressMessages(lme4::lmer(as.formula(paste(x, "~ 1 + (1|site)")), data = dd))
    icc_of(m)
  }, error = function(e) unname(icc_F(dd[[x]], dd$site)["icc"]))
}

## Site-level exposure: single-level against random-intercept (Kenward-Roger)
expo_test <- function(D, v, x) {
  d <- D$d_base
  dd <- d[!is.na(d[[v]]) & !is.na(d[[x]]), ]
  extra_n <- c(if (has_var(dd, "dx")) "dx", if (has_var(dd, "age")) "age",
               if (has_var(dd, "female")) "female")
  dd <- dd[complete.cases(dd[, c(v, x, extra_n), drop = FALSE]), ]
  dd$site <- droplevels(dd$site); dd$z <- zs(dd[[v]])
  extra_m <- extra_n
  if ("age" %in% extra_n) {
    dd$age_b <- ave(dd$age, dd$site); dd$age_w <- dd$age - dd$age_b
    extra_m <- c(setdiff(extra_n, "age"), "age_w", "age_b")
  }
  fa <- as.formula(paste("z ~", paste(c(x, extra_n), collapse = " + ")))
  a <- summary(lm(fa, data = dd))$coef[x, ]
  m <- lmerTest::lmer(as.formula(paste("z ~", paste(c(x, extra_m), collapse = " + "),
                                       "+ (1|site)")), data = dd)
  c1 <- summary(m, ddf = "Kenward-Roger")$coef[x, ]
  data.frame(analysis = c("Single-level (site ignored)", "Random site intercept (Kenward-Roger)"),
             est = c(a[1], c1[1]), se = c(a[2], c1[2]),
             df = c(nrow(dd) - length(extra_n) - 2, c1[3]),
             p = c(a[4], c1[5]), N = nrow(dd), sites = nlevels(dd$site),
             rho_x = rho_x(dd, x))
}

## Real-data null experiment (Listing 5)
null_exp <- function(D, v, B, progress = NULL) {
  d <- D$d_base
  dd <- d[!is.na(d[[v]]), ]; dd$site <- droplevels(dd$site); S <- levels(dd$site)
  rej <- matrix(NA, B, 4)
  f1 <- as.formula(paste(v, "~ x")); f2 <- as.formula(paste(v, "~ x + (1|site)"))
  g1 <- as.formula(paste(v, "~ w")); g2 <- as.formula(paste(v, "~ w + (1|site)"))
  for (b in seq_len(B)) {
    dd$x <- as.numeric(dd$site %in% sample(S, floor(length(S) / 2)))
    dd$w <- ave(rbinom(nrow(dd), 1, .5), dd$site, FUN = function(z) sample(z))
    rej[b, 1] <- summary(lm(f1, dd))$coef["x", 4] < .05
    m <- suppressMessages(lmerTest::lmer(f2, dd))
    rej[b, 2] <- summary(m, ddf = "Satterthwaite")$coef["x", 5] < .05
    rej[b, 3] <- summary(lm(g1, dd))$coef["w", 4] < .05
    m2 <- suppressMessages(lmerTest::lmer(g2, dd))
    rej[b, 4] <- summary(m2, ddf = "Satterthwaite")$coef["w", 5] < .05
    if (!is.null(progress)) progress()
  }
  r <- colMeans(rej, na.rm = TRUE)
  data.frame(var = v,
             assignment = rep(c("Assigned to whole sites", "Assigned within sites"), each = 2),
             analysis = rep(c("Site ignored", "Random site intercept"), 2),
             rate = r, B = B, sites = length(S))
}

## Crossed instrument (7) and institution/cohort decomposition
crossed_instr <- function(D, v) {
  d <- D$d_base
  dq <- d[!is.na(d[[v]]) & !is.na(d$iqtest), ]
  rhs <- c(if (has_var(dq, "dx")) "dx", if (has_var(dq, "female")) "female")
  if (has_var(dq, "age")) {
    dq <- dq[!is.na(dq$age), ]
    dq$age_b <- ave(dq$age, dq$site); dq$age_w <- dq$age - dq$age_b
    rhs <- c(rhs, "age_w", "age_b")
  }
  rhs <- if (length(rhs)) paste(rhs, collapse = " + ") else "1"
  mx <- lme4::lmer(as.formula(paste(v, "~", rhs, "+ (1|site) + (1|iqtest)")), data = dq)
  mn <- lme4::lmer(as.formula(paste(v, "~", rhs, "+ (1|site)")), data = dq)
  vx <- as.data.frame(VarCorr(mx)); vn <- as.data.frame(VarCorr(mn))
  tx <- sum(vx$vcov); tn <- sum(vn$vcov)
  nest <- table(dq$site, dq$iqtest)
  data.frame(model = c("Site only", "Site + instrument (crossed)"),
             site_icc = c(vn$vcov[vn$grp == "site"] / tn, vx$vcov[vx$grp == "site"] / tx),
             instrument_icc = c(NA, vx$vcov[vx$grp == "iqtest"] / tx),
             p_instrument = c(NA, anova(mn, mx, refit = FALSE)$`Pr(>Chisq)`[2] / 2),
             N = nrow(dq), instruments = length(unique(dq$iqtest)),
             sites = length(unique(dq$site)),
             sites_multi_instr = sum(rowSums(nest > 0) > 1),
             instr_multi_site = sum(colSums(nest > 0) > 1))
}

inst_cohort <- function(D, vars) {
  d <- D$d_base
  do.call(rbind, lapply(vars, function(v) {
    dd <- d[!is.na(d[[v]]), ]
    m <- lme4::lmer(as.formula(paste(v, "~ 1 + (1|institution) + (1|institution:site)")), data = dd)
    vc <- as.data.frame(VarCorr(m)); tot <- sum(vc$vcov)
    data.frame(var = v, institution = vc$vcov[vc$grp == "institution"] / tot,
               cohort_within = vc$vcov[vc$grp == "institution:site"] / tot)
  }))
}

## Measure layer: within, between and pooled correlations (8)-(9)
mv_fit <- function(D, vars, sub_dx = FALSE, model = TRUE) {
  d <- D$d_base
  if (sub_dx) d <- d[!is.na(d$dx) & d$dx == 1, ]
  dd <- d[complete.cases(d[, vars]), ]; dd$site <- droplevels(dd$site)
  for (v in vars) dd[[v]] <- zs(dd[[v]])
  naive <- cor(dd[, vars])
  cen <- dd[, vars] - apply(dd[, vars], 2, function(x) ave(x, dd$site))
  within <- cor(cen)
  smeans <- aggregate(dd[, vars], list(site = dd$site), mean)
  betw_desc <- cor(smeans[, vars])
  pairs <- t(combn(vars, 2))
  idx <- cbind(match(pairs[, 1], vars), match(pairs[, 2], vars))
  iccs <- sapply(vars, function(v) icc_of(lme4::lmer(as.formula(paste(v, "~ 1 + (1|site)")), data = dd)))
  out <- data.frame(m1 = pairs[, 1], m2 = pairs[, 2], N = nrow(dd), sites = nlevels(dd$site),
                    pooled = naive[idx], within_centred = within[idx],
                    between_sitemeans = betw_desc[idx],
                    between_model = NA, within_model = NA,
                    icc1 = iccs[pairs[, 1]], icc2 = iccs[pairs[, 2]],
                    stringsAsFactors = FALSE)
  msg <- NULL
  if (model) {
    fit <- tryCatch({
      long <- do.call(rbind, lapply(seq_along(vars), function(j)
        data.frame(id = dd$id, site = dd$site, measure = vars[j], mi = j, y = dd[[vars[j]]])))
      long$measure <- factor(long$measure, levels = vars)
      long <- long[order(long$site, long$id, long$mi), ]
      nlme::lme(y ~ 0 + measure, random = list(site = pdSymm(~ 0 + measure)),
                correlation = corSymm(form = ~ mi | site / id),
                weights = varIdent(form = ~ 1 | measure), data = long, method = "REML",
                control = lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200))
    }, error = function(e) e)
    if (inherits(fit, "error")) {
      msg <- paste0("The multivariate model was not fitted (", conditionMessage(fit), ").",
                   " Within-site correlations are shown from site-mean centring.")
    } else {
      Ru <- cov2cor(as.matrix(getVarCov(fit)))
      Rp <- corMatrix(fit$modelStruct$corStruct)[[1]]
      out$between_model <- Ru[idx]; out$within_model <- Rp[idx]
    }
  } else {
    msg <- "The multivariate model was not requested; within-site correlations come from site-mean centring and between-site correlations from site means."
  }
  bw <- ifelse(is.na(out$between_model), out$between_sitemeans, out$between_model)
  ww <- ifelse(is.na(out$within_model), out$within_centred, out$within_model)
  out$pooled_eq9 <- sqrt(out$icc1 * out$icc2) * bw + sqrt((1 - out$icc1) * (1 - out$icc2)) * ww
  list(tab = out, msg = msg)
}

## Time layer: correlation between visits (11)
visit_corr <- function(times, b0, b1, b01, se2, phi, su2) {
  k <- length(times)
  C <- matrix(NA, k, k)
  for (i in 1:k) for (j in 1:k) {
    C[i, j] <- su2 + b0 + times[i] * times[j] * b1 + (times[i] + times[j]) * b01 +
      se2 * phi^abs(i - j)
  }
  V <- diag(C)
  list(cov = C, cor = C / sqrt(outer(V, V)))
}

sim_long <- function(S, n, times, b0, b1, b01, se2, phi, su2, bt) {
  D <- matrix(c(b0, b01, b01, b1), 2)
  if (min(eigen(D, only.values = TRUE)$values) <= 0) stop("The participant covariance matrix is not positive definite: reduce |sigma_b01|.")
  L <- chol(D)
  out <- list()
  for (s in 1:S) {
    u <- rnorm(1, 0, sqrt(su2))
    for (i in 1:n) {
      b <- as.vector(rnorm(2) %*% L)
      e <- numeric(length(times)); e[1] <- rnorm(1, 0, sqrt(se2))
      if (length(times) > 1) for (t in 2:length(times))
        e[t] <- phi * e[t - 1] + rnorm(1, 0, sqrt(se2 * (1 - phi^2)))
      out[[length(out) + 1]] <- data.frame(site = paste0("S", s), id = paste0("S", s, "_", i),
                                           time = times, visit = seq_along(times),
                                           y = 50 + bt * times + u + b[1] + b[2] * times + e)
    }
  }
  x <- do.call(rbind, out); x$site <- factor(x$site); x
}

fit_long <- function(dd, y = "y", tm = "time", irregular = FALSE) {
  dd <- dd[!is.na(dd[[y]]) & !is.na(dd[[tm]]), ]
  dd$yy <- dd[[y]]; dd$tt <- dd[[tm]]
  dd <- dd[order(dd$site, dd$id, dd$tt), ]
  dd$site <- droplevels(factor(dd$site)); dd$id <- factor(dd$id)
  ctl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200, returnObject = TRUE)
  m_ind <- nlme::lme(yy ~ tt, random = list(site = ~ 1, id = ~ 1 + tt), data = dd,
                     method = "REML", control = ctl)
  cs <- if (irregular) corCAR1(form = ~ tt | site / id) else corAR1(form = ~ visit | site / id)
  if (!irregular) dd$visit <- ave(dd$tt, dd$id, FUN = seq_along)
  m_ar <- tryCatch(nlme::lme(yy ~ tt, random = list(site = ~ 1, id = ~ 1 + tt), data = dd,
                             correlation = cs, method = "REML", control = ctl),
                   error = function(e) NULL)
  comp <- function(m) {
    pm <- pdMatrix(m$modelStruct$reStruct)
    s2 <- m$sigma^2
    Did <- pm$id * s2
    phi <- if (!is.null(m$modelStruct$corStruct)) coef(m$modelStruct$corStruct, unconstrained = FALSE) else 0
    c(su2 = pm$site[1, 1] * s2, b0 = Did[1, 1], b1 = Did[2, 2], b01 = Did[1, 2],
      se2 = s2, phi = unname(phi), slope = unname(fixef(m)[2]), AIC = AIC(m), BIC = BIC(m))
  }
  list(ind = comp(m_ind), ar = if (!is.null(m_ar)) comp(m_ar) else NULL,
       N = length(unique(dd$id)), S = nlevels(dd$site), obs = nrow(dd))
}

## -----------------------------------------------------------------------------
##  UI
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
##  Decision guide (Figure 9 and the analysis cycle of Section 6.3)
## -----------------------------------------------------------------------------

LAYER_COL <- c("Structure and inference" = "#F28E2B", "Site layer" = "#4E79A7",
               "Measure layer" = "#B07AA1", "Time layer" = "#59A14F")
LAYER_TINT <- c("Structure and inference" = "#FDEBD9", "Site layer" = "#DDE7F2",
                "Measure layer" = "#EFE3EC", "Time layer" = "#E1EFDC")

## Exposure effect and its between-site variation (random slope test)
slope_test <- function(dd, y, x) {
  dd <- dd[!is.na(dd[[y]]) & !is.na(dd[[x]]), ]
  dd$site <- droplevels(dd$site); dd$z <- zs(dd[[y]])
  cov <- character(0)
  if (x != "age" && has_var(dd, "age")) {
    dd <- dd[!is.na(dd$age), ]
    dd$age_b <- ave(dd$age, dd$site); dd$age_w <- dd$age - dd$age_b
    cov <- c(cov, "age_w", "age_b")
  }
  if (x != "female" && has_var(dd, "female")) cov <- c(cov, "female")
  if (x != "dx" && has_var(dd, "dx")) cov <- c(cov, "dx")
  rhs <- paste(c(x, cov), collapse = " + ")
  fri <- as.formula(paste("z ~", rhs, "+ (1|site)"))
  frs <- as.formula(paste0("z ~ ", rhs, " + (1 + ", x, "|site)"))
  ctl <- lmerControl(optimizer = "bobyqa")
  mri <- lmerTest::lmer(fri, data = dd)
  cri <- summary(mri, ddf = "Kenward-Roger")$coef[x, ]
  out <- list(est_ri = cri[1], se_ri = cri[2], df_ri = cri[3], p_ri = cri[5],
              tau = NA, p_tau = NA, est = cri[1], se = cri[2], df = cri[3], singular = FALSE,
              N = nrow(dd), S = nlevels(dd$site))
  mrs <- tryCatch(lmerTest::lmer(frs, data = dd, control = ctl), error = function(e) NULL)
  if (!is.null(mrs)) {
    crs <- summary(mrs, ddf = "Kenward-Roger")$coef[x, ]
    vc <- as.data.frame(VarCorr(mrs))
    tau <- sqrt(vc$vcov[vc$grp == "site" & vc$var1 == x & is.na(vc$var2)])
    l0 <- lme4::lmer(fri, data = dd, REML = FALSE)
    l1 <- lme4::lmer(frs, data = dd, REML = FALSE, control = ctl)
    chi <- max(0, as.numeric(2 * (logLik(l1) - logLik(l0))))
    out$tau <- tau
    out$p_tau <- 0.5 * pchisq(chi, 1, lower.tail = FALSE) + 0.5 * pchisq(chi, 2, lower.tail = FALSE)
    out$est <- crs[1]; out$se <- crs[2]; out$df <- crs[3]; out$singular <- isSingular(mrs)
  }
  out
}

## Facts for the decision guide, taken from the loaded data
dg_facts <- function(D, y, x, joint) {
  d <- D$d_base
  k <- match(y, D$meas$var)
  elig <- d
  if (!is.na(k) && !D$meas$both[k] && "dx" %in% names(d)) elig <- d[!is.na(d$dx) & d$dx == 1, ]
  dd <- d[!is.na(d[[y]]), ]
  if (grepl("^ADOS", y) && "ADOS_module" %in% names(dd)) dd <- dd[!is.na(dd$ADOS_module), ]
  dd$site <- droplevels(dd$site)
  S <- nlevels(dd$site); N <- nrow(dd); nbar <- N / S
  m0 <- lme4::lmer(as.formula(paste(y, "~ 1 + (1|site)")), data = dd)
  icc <- icc_of(m0)
  ci <- icc_F(dd[[y]], dd$site)
  if (is.null(x) || !x %in% names(dd) || !has_var(dd, x)) {
    rhox <- 0; st <- NULL
  } else {
    rhox <- rho_x(dd, x)
    st <- tryCatch(slope_test(dd, y, x), error = function(e) NULL)
  }
  ## instrument crossed with site
  instr <- FALSE
  if (!is.null(D$instrument) && D$instrument %in% names(dd) &&
      !(D$source %in% c("demo", "abide") && !y %in% c("FIQ", "VIQ", "PIQ"))) {
    q <- dd[!is.na(dd[[D$instrument]]), ]
    if (nrow(q)) {
      nest <- table(droplevels(q$site), q[[D$instrument]])
      instr <- ncol(nest) > 1 && (sum(rowSums(nest > 0) > 1) > 0 || sum(colSums(nest > 0) > 1) > 0)
    }
  }
  ## institution contributing several cohorts
  inst_multi <- any(table(unique(dd[, c("institution", "site")])$institution) > 1)
  inst_level <- "unchecked"
  if (inst_multi) {
    ic <- tryCatch(inst_cohort(list(d_base = dd), y), error = function(e) NULL)
    if (!is.null(ic)) {
      inst_level <- if (ic$cohort_within < 0.25 * (ic$institution + ic$cohort_within)) "institution"
                    else if (ic$institution < 0.25 * (ic$institution + ic$cohort_within)) "cohort" else "nested"
    }
  }
  sites_all <- nlevels(droplevels(elig$site))
  in_sites <- elig[elig$site %in% levels(dd$site), ]
  long <- isTRUE(D$longitudinal)
  fu <- NA; irregular <- FALSE; int_differ <- FALSE
  if (long) {
    dl <- D$d[!is.na(D$d[[y]]) & !is.na(D$d$time), ]
    nv <- tapply(dl$time, dl$id, length)
    fu_ids <- names(nv)[nv > 1]
    fu <- length(unique(dl$site[dl$id %in% fu_ids]))
    irregular <- length(unique(round(dl$time, 2))) > max(nv)
    last <- tapply(dl$time, dl$id, max); ls <- dl$site[match(names(last), dl$id)]
    int_differ <- length(unique(ls)) > 1 &&
      tryCatch(kruskal.test(last ~ factor(ls))$p.value < 0.05, error = function(e) FALSE)
  }
  list(y = y, x = x, S = S, N = N, nbar = nbar, icc = icc, icc_lo = ci[["lo"]], icc_hi = ci[["hi"]],
       rhox = rhox, tau = if (!is.null(st)) st$tau else NA, p_tau = if (!is.null(st)) st$p_tau else NA,
       est = if (!is.null(st)) st$est else NA, se = if (!is.null(st)) st$se else NA,
       df_ri = if (!is.null(st)) st$df_ri else NA, p_ri = if (!is.null(st)) st$p_ri else NA,
       instr = instr, inst_multi = inst_multi, inst_level = inst_level,
       n_meas = max(1, length(unique(c(y, joint)))),
       miss_sites = sites_all - S,
       miss_pct = round(100 * mean(is.na(in_sites[[y]])), 1),
       long = long, fu_sites = fu, irregular = irregular, int_differ = int_differ,
       singular = isSingular(m0))
}

## Decision engine: returns the decisions taken, the model and the code
dg_decide <- function(f) {
  nz <- function(v, alt = NA) if (is.null(v) || length(v) == 0 || is.na(v)) alt else v
  S <- nz(f$S, 1); nbar <- nz(f$nbar, 1); icc <- nz(f$icc, 0); rhox <- nz(f$rhox, 0)
  deff_y <- 1 + (nbar - 1) * icc
  deff_x <- 1 + (nbar - 1) * rhox * icc
  N <- S * nbar
  xl <- if (rhox >= 0.8) "site" else if (rhox <= 0.1) "participant" else "both"
  rows <- list()
  add <- function(layer, step, q, ans, act, active, why, warn = "") {
    rows[[length(rows) + 1]] <<- data.frame(layer = layer, step = step, question = q, answer = ans,
                                            action = act, active = active, why = why, warn = warn,
                                            stringsAsFactors = FALSE)
  }
  ## 1 multilevel needed
  ml <- !(S < 2 || (deff_x < 1.1 && icc < 0.01))
  add("Structure and inference", "1-2", "Is a multilevel analysis needed?",
      if (ml) "Yes" else "No",
      if (ml) "Random site intercept (1 | site); report ICC with interval, mean n per site, DEFF and N_eff" else "Single-level analysis; report the site ICC and DEFF that justify it",
      TRUE,
      if (S < 2) "Data come from a single site." else
        sprintf("Site ICC %s%s; %s sites with a mean of %s participants; DEFF for the outcome %s and for this exposure %s; N_eff %s of %s.",
                f3(icc), if (is.na(nz(f$icc_lo)) || is.na(nz(f$icc_hi))) "" else sprintf(" (95%% interval %s to %s)", f3(f$icc_lo), f3(f$icc_hi)),
                S, f2(nbar), f2(deff_y), f2(deff_x), round(N / deff_x), round(N)),
      if (ml) "Do not drop site because its variance test is not significant: an ICC is negligible only when the design effect is close to 1." else "")
  ## 2 exposure level
  act2 <- switch(xl,
    site = sprintf("Enter %s as a site-level term z_s and test it against between-site variation with Kenward-Roger df (at most %d)", f$x, max(1, S - 2)),
    participant = sprintf("Enter %s as a participant-level term; the multilevel model uses its within-site information", f$x),
    both = sprintf("Split %s into a within-site deviation and the site mean (%s_w, %s_b)", f$x, f$x, f$x))
  add("Site layer", "3-4", sprintf("Where does the exposure (%s) sit?", f$x),
      sprintf("%s (rho_x = %s)", c(site = "Site property", participant = "Participant property", both = "Within and between sites")[[xl]], f3(rhox)),
      act2, ml,
      sprintf("The exposure's own ICC is %s. The variance ratio of (5), correct to naive, is about %s.", f3(rhox), f2(deff_x)),
      if (xl == "site") sprintf("Never test a site-level exposure with a single-level model: its %s observations behave like about %s.", round(N), S) else "")
  ## 3 random slope
  slope <- FALSE
  if (xl != "site" && ml) {
    tau <- nz(f$tau); pt <- nz(f$p_tau)
    plan <- identical(f$stage, "plan")
    slope <- if (plan) !is.na(tau) && tau > 0 else !is.na(pt) && pt < 0.05
    pi_txt <- ""
    if (slope && !is.na(nz(f$est)) && !is.na(nz(f$se))) {
      hw <- qt(0.975, max(1, S - 2)) * sqrt(f$se^2 + tau^2)
      pi_txt <- sprintf(" Prediction interval for a new site (16): %s to %s SD.", f2(f$est - hw), f2(f$est + hw))
    }
    if (plan) add("Site layer", "5", sprintf("Does the effect of %s vary between sites?", f$x),
        if (slope) sprintf("Expected: tau = %s SD", f2(tau)) else "Not expected",
        if (slope) sprintf("Pre-specify (1 + %s | site); retain it when the likelihood-ratio test of tau^2 = 0 supports it, and report the prediction interval for a new site", f$x)
        else sprintf("Pre-specify the likelihood-ratio test of tau^2 = 0 for (1 + %s | site); random intercept when it is not supported", f$x),
        TRUE,
        if (slope) sprintf("With tau = %s the standard error falls with the number of sites, not with participants per site (17).", f2(tau)) else "No between-site variation in the effect is assumed in the power calculation.",
        "Fixed site terms and ComBat-type harmonisation adjust site means only; they do not cover between-site variation in the effect.")
    else add("Site layer", "5", sprintf("Does the effect of %s vary between sites?", f$x),
        if (is.na(pt)) "Not tested" else if (slope) "Yes" else "No",
        if (slope) sprintf("Random slope (1 + %s | site); report the average effect with its random-slope SE and the prediction interval", f$x)
        else if (is.na(pt)) sprintf("Fit (1 + %s | site) and test tau^2 = 0 by likelihood ratio", f$x)
        else "Random intercept only; report tau and its test",
        TRUE,
        if (is.na(pt)) "No estimate of tau is available." else
          sprintf("tau = %s SD (likelihood-ratio %s).%s", f2(tau), fpt(pt), pi_txt),
        "Fixed site terms and ComBat-type harmonisation adjust site means only; they do not cover between-site variation in the effect.")
  } else {
    add("Site layer", "5", sprintf("Does the effect of %s vary between sites?", f$x), "Not applicable",
        "No random slope", FALSE,
        if (!ml) "Single-level analysis." else "A site-level exposure cannot vary within a site, so there is no slope to vary.")
  }
  ## 4 site decomposition
  unit <- "site"; crossed <- isTRUE(f$instr)
  lev <- nz(f$inst_level, "unchecked")
  if (isTRUE(f$inst_multi)) unit <- switch(lev, institution = "institution", nested = "institution/site", "site")
  parts <- c(if (crossed) "instrument crossed with site", if (isTRUE(f$inst_multi)) "several cohorts from one institution")
  add("Site layer", "6", "Is site made of identifiable parts?",
      if (length(parts)) { a <- paste(parts, collapse = "; "); paste0(toupper(substr(a, 1, 1)), substring(a, 2)) } else "No",
      paste(c(if (crossed) "Add (1 | instrument) crossed with site",
              if (isTRUE(f$inst_multi)) switch(lev, institution = "Use institution as the site unit",
                                               cohort = "Keep cohort as the site unit",
                                               nested = "Use (1 | institution/site)",
                                               "Compare institution and cohort as the site unit")),
            collapse = "; "),
      ml && length(parts) > 0,
      if (length(parts)) paste0("Site unit: ", unit, ".") else "One instrument per measure and one cohort per institution.")
  ## 5 several measures
  nm <- nz(f$n_meas, 1)
  mvfull <- nm > 1 && S >= 10 && nm <= 4
  add("Measure layer", "7", "Are several measures analysed together?",
      if (nm > 1) sprintf("Yes (%d measures)", nm) else "No",
      if (nm <= 1) "Univariate model" else paste0(
        if (mvfull) "Multivariate model (8) with unstructured Sigma_u and Sigma_e" else "Centre each measure on its site mean (an unstructured Sigma_u would be singular with these sites and measures)",
        "; report within-site correlations",
        if (isTRUE(f$predict)) "; validate the prediction model by leaving whole sites out" else ""),
      nm > 1,
      if (nm > 1) "Pooled correlations mix the relationship in participants with agreement between site practices (9)." else "Only one measure enters the analysis.",
      if (nm > 1) "Do not correlate or combine raw measures across sites." else "")
  ## 6 time
  long <- isTRUE(f$long)
  cs <- if (isTRUE(f$irregular)) "continuous-time AR(1), corCAR1" else "AR(1), corAR1"
  add("Time layer", "8", "Are there repeated visits?",
      if (long) "Yes" else "No",
      if (long) paste0("Participant random intercept and slope; residual ", cs, " against independent residuals by AIC/BIC",
                       if (isTRUE(f$int_differ)) "; site random slope for time" else "",
                       if (!is.na(nz(f$fu_sites)) && f$fu_sites < S) sprintf("; restrict to the %d sites with follow-up", f$fu_sites) else "")
      else "Cross-sectional model",
      long,
      if (long) paste0("Visit intervals ", if (isTRUE(f$int_differ)) "differ" else "do not differ", " between sites; visit times are ",
                       if (isTRUE(f$irregular)) "irregular." else "regular.") else "One assessment per participant.",
      if (long) "Do not analyse change scores or averaged visits without site." else "")
  ## 7 missing
  ms <- nz(f$miss_sites, 0); mp <- nz(f$miss_pct, 0)
  add("Structure and inference", "9", "Are values missing, and at which level?",
      if (ms > 0) sprintf("%d site(s) without the outcome; %s%% missing within sites", ms, mp) else if (mp > 0) sprintf("%s%% missing within sites", mp) else "None",
      if (ms > 0 || mp > 0) paste0(if (ms > 0) "Treat sites without the measure as structural missingness; " else "",
                                   if (mp > 0) "multilevel imputation with site as a random effect (mice 2l.pmm) or full-information likelihood" else "state which sites contribute")
      else "Complete data",
      ms > 0 || mp > 0,
      if (ms > 0) "Complete-case comparisons across measures compare different sets of sites." else "No site-level missingness.",
      if (mp > 0) "An imputation model without site transfers between-site covariance into imputed values." else "")
  ## 8 number of sites
  few <- if (S < 10) "fewer than 10" else if (S < 20) "10 to 19" else "20 or more"
  add("Structure and inference", "10", "How many sites inform the site variance?", sprintf("%d (%s)", S, few),
      if (S < 10) "Report the ICC as a bound with its interval; Kenward-Roger inference; Bayesian fit with a half-normal prior on sigma_u (and tau)"
      else if (S < 20) "Report ICC intervals; Kenward-Roger inference; report singular fits; check within-site comparisons against a fixed-site analysis"
      else "REML with bootstrap ICC intervals; Kenward-Roger inference for site-level or random-slope comparisons",
      ml,
      if (isTRUE(f$singular)) "The empty model is singular: the site variance is estimated at zero and should be reported as such." else "20 to 30 clusters are recommended for a stable between-cluster variance.",
      if (S < 10) "Do not rank measures by site dependence unless their intervals separate." else "")
  tab <- do.call(rbind, rows)
  tab$action <- paste0(toupper(substr(tab$action, 1, 1)), substring(tab$action, 2))
  list(tab = tab, ml = ml, slope = slope, xl = xl, unit = unit, crossed = crossed, long = long,
       mvfull = mvfull, deff_x = deff_x, deff_y = deff_y, neff = N / deff_x, N = N, S = S)
}

## Analysis script for the decisions taken
dg_code <- function(f, dc, covs = c("age", "female"), instr = "iqtest") {
  y <- f$y; x <- f$x; ht <- if (is.null(f$htype)) "group" else f$htype
  L <- c("## Analysis script generated by HierDep Explorer (decision guide)",
         if (!is.null(f$h0)) paste0("## H0: ", f$h0, "    H1: ", f$h1) else NULL,
         "## d: one row per participant (per visit for longitudinal data) with columns site, id and the variables below",
         "library(lme4); library(lmerTest); library(nlme)", "")
  cv <- setdiff(covs, x)
  fx <- character(0)
  if (ht == "change") {
  } else if (dc$xl == "both") {
    L <- c(L, sprintf("d$%s_b <- ave(d$%s, d$site, FUN = function(v) mean(v, na.rm = TRUE)); d$%s_w <- d$%s - d$%s_b", x, x, x, x, x))
    fx <- c(fx, paste0(x, "_w"), paste0(x, "_b"))
  } else fx <- c(fx, x)
  if ("age" %in% cv) {
    L <- c(L, "d$age_b <- ave(d$age, d$site, FUN = function(v) mean(v, na.rm = TRUE)); d$age_w <- d$age - d$age_b")
    fx <- c(fx, "age_w", "age_b")
  }
  fx <- c(fx, setdiff(cv, "age"))
  if (dc$unit != "site" && !grepl("/", dc$unit)) L <- c(L, "## site unit: institution")
  L <- c(L, "", "## 1. site ICC, design effect and effective sample size",
         sprintf("m0 <- lmer(%s ~ 1 + (1 | %s), data = d)", y, if (dc$unit == "institution") "institution" else "site"),
         "vc <- as.data.frame(VarCorr(m0)); icc <- vc$vcov[1] / sum(vc$vcov)",
         "nbar <- nobs(m0) / ngrps(m0)[[1]]",
         sprintf("deff <- 1 + (nbar - 1) * %s * icc; neff <- nobs(m0) / deff   # rho_x = %s", f3(f$rhox), f3(f$rhox)),
         "icc_ci <- quantile(bootMer(m0, function(m) { v <- as.data.frame(VarCorr(m)); v$vcov[1] / sum(v$vcov) }, nsim = 500)$t, c(.025, .975))",
         "")
  rhs <- if (length(fx)) paste(fx, collapse = " + ") else "1"
  rhs_nox <- setdiff(fx, c(x, paste0(x, "_w"), paste0(x, "_b")))
  rhs_nox <- if (length(rhs_nox)) paste(rhs_nox, collapse = " + ") else "1"
  if (ht == "corr") {
    L <- c(L, "## 2. hypothesis test: within-site correlation between the two measures",
           sprintf("cc <- d[complete.cases(d[, c(\"%s\", \"%s\")]), ]", y, x),
           sprintf("w1 <- cc$%s - ave(cc$%s, cc$site); w2 <- cc$%s - ave(cc$%s, cc$site)", y, y, x, x),
           "r_w <- cor(w1, w2); df_w <- nrow(cc) - length(unique(cc$site)) - 3   # site means removed",
           "z <- atanh(r_w) * sqrt(df_w); p <- 2 * pnorm(-abs(z)); ci <- tanh(atanh(r_w) + c(-1, 1) * qnorm(.975) / sqrt(df_w))",
           "")
  }
  if (!dc$ml) {
    L <- c(L, "## 2. single-level analysis (site ICC near zero, DEFF < 1.1 for this exposure)",
           sprintf("fit <- lm(%s ~ %s, data = d)", y, rhs), "summary(fit)")
  } else if (!dc$long) {
    sx <- if (dc$xl == "both") paste0(x, "_w") else x
    re <- if (dc$slope) sprintf("(1 + %s | %s)", sx, dc$unit) else sprintf("(1 | %s)", dc$unit)
    if (dc$crossed) re <- paste(re, sprintf("+ (1 | %s)", instr))
    L <- c(L, "## 2. multilevel model", sprintf("fit <- lmer(%s ~ %s + %s, data = d, control = lmerControl(optimizer = \"bobyqa\"))", y, rhs, re),
           "isSingular(fit)", "summary(fit, ddf = \"Kenward-Roger\")")
    if (dc$xl != "site") L <- c(L, "",
           "## between-site variation of the effect: likelihood-ratio test of tau^2 = 0",
           sprintf("f0 <- lmer(%s ~ %s + (1 | %s), data = d, REML = FALSE)", y, rhs, dc$unit),
           sprintf("f1 <- lmer(%s ~ %s + (1 + %s | %s), data = d, REML = FALSE, control = lmerControl(optimizer = \"bobyqa\"))", y, rhs, sx, dc$unit),
           "chi <- max(0, 2 * (logLik(f1) - logLik(f0))); p_tau <- 0.5 * pchisq(chi, 1, lower.tail = FALSE) + 0.5 * pchisq(chi, 2, lower.tail = FALSE)")
    if (dc$slope) L <- c(L, "",
           "## prediction interval for a new site (16)",
           sprintf("b <- summary(fit, ddf = \"Kenward-Roger\")$coef[\"%s\", ]; vc <- as.data.frame(VarCorr(fit))", sx),
           sprintf("tau <- sqrt(vc$vcov[vc$var1 == \"%s\" & is.na(vc$var2)][1]); S <- ngrps(fit)[[\"%s\"]]", sx, sub("/.*", "", dc$unit)),
           "pi <- b[1] + c(-1, 1) * qt(.975, S - 2) * sqrt(b[2]^2 + tau^2)")
    if (dc$S < 20) L <- c(L, "", "## sensitivity: site as a fixed effect for within-site comparisons",
                          sprintf("summary(lm(%s ~ %s + factor(site), data = d))", y, rhs))
  } else {
    rs <- if (isTRUE(f$int_differ)) "~ 1 + time" else "~ 1"
    cs <- if (isTRUE(f$irregular)) "corCAR1(form = ~ time | site / id)" else "corAR1(form = ~ visit | site / id)"
    L <- c(L, "## 2. longitudinal multilevel model",
           if (!isTRUE(f$irregular)) "d$visit <- ave(d$time, d$id, FUN = seq_along)" else NULL,
           if (!is.na(f$fu_sites) && f$fu_sites < f$S) "d <- d[d$site %in% unique(d$site[duplicated(d$id)]), ]   # sites with follow-up" else NULL,
           sprintf("m_ind <- lme(%s ~ %s, random = list(site = %s, id = ~ 1 + time), data = d, na.action = na.omit)", y,
                   if (ht == "groupchange") paste0("time * ", x, " + ", rhs_nox) else if (ht == "change") paste("time +", rhs_nox) else paste("time +", rhs), rs),
           sprintf("m_ar  <- update(m_ind, correlation = %s)", cs),
           "AIC(m_ind, m_ar); BIC(m_ind, m_ar)",
           "summary(m_ar)$tTable   # hypothesis: the time term (change) or the time-by-group term (group difference in change)")
  }
  if (dc$S < 10 && dc$ml) L <- c(L, "", "## few sites: Bayesian fit with a weakly informative prior on the site SD",
                                 sprintf("# brms::brm(%s ~ %s + (1 | site), data = d, prior = brms::prior(normal(0, 0.5), class = sd))", y, rhs))
  if (nz_num(f$n_meas) > 1) {
    L <- c(L, "", "## 3. several measures: within- and between-site correlations",
           "# vars <- c(...)   # the measures analysed together, standardised",
           "# cen <- d[, vars] - apply(d[, vars], 2, function(v) ave(v, d$site)); cor(cen, use = \"pairwise\")   # within sites")
    if (dc$mvfull) L <- c(L,
           "# long <- reshape to one row per participant and measure (columns id, site, measure, mi, y)",
           "# lme(y ~ 0 + measure, random = list(site = pdSymm(~ 0 + measure)), correlation = corSymm(form = ~ mi | site / id),",
           "#     weights = varIdent(form = ~ 1 | measure), data = long, control = lmeControl(opt = \"optim\"))")
    if (isTRUE(f$predict)) L <- c(L, "# leave-one-site-out validation: for (s in unique(d$site)) { train on d[d$site != s, ]; predict d[d$site == s, ] }")
  }
  if (nz_num(f$miss_pct) > 0) L <- c(L, "", "## 4. missing values: multilevel imputation with site as the cluster",
           "# library(mice); pm <- make.predictorMatrix(d); pm[, \"site\"] <- -2",
           sprintf("# meth <- make.method(d); meth[\"%s\"] <- \"2l.pmm\"   # needs package miceadds", y),
           "# imp <- mice(d, method = meth, predictorMatrix = pm, m = 20)")
  paste(L, collapse = "\n")
}
nz_num <- function(v) if (is.null(v) || length(v) == 0 || is.na(v)) 0 else v

## -----------------------------------------------------------------------------
##  Objective, hypothesis, planning-stage power and analysis-stage test
## -----------------------------------------------------------------------------

HYP_TYPES <- c("Difference between groups (participant-level exposure)" = "group",
               "Effect of a site-level condition" = "sitelevel",
               "Within-site relationship between two measures" = "corr",
               "Change over time" = "change",
               "Group difference in change over time" = "groupchange")

hyp_info <- function(ht, ylab, xlab, side) {
  ht <- ht %||% "group"; side <- side %||% "two"; ylab <- ylab %||% "the outcome"; xlab <- xlab %||% "the exposure"
  par <- switch(ht,
    group = list(sym = "\u03b2\u2081", tex = "\\beta_1",
                 desc = sprintf("the average difference in %s for %s across sites (SD units)", ylab, xlab),
                 eq = "$$Y_{is}=\\beta_0+(\\beta_1+v_s)x_{is}+X_{is}^\\top\\beta+u_s+\\varepsilon_{is},\\quad u_s\\sim N(0,\\sigma^2_u),\\ v_s\\sim N(0,\\tau^2)$$"),
    sitelevel = list(sym = "\u03b3", tex = "\\gamma",
                     desc = sprintf("the effect of %s, a site-level condition, on %s (SD units), estimated from between-site information", xlab, ylab),
                     eq = "$$Y_{is}=\\beta_0+\\gamma z_s+X_{is}^\\top\\beta+u_s+\\varepsilon_{is},\\quad u_s\\sim N(0,\\sigma^2_u)$$"),
    corr = list(sym = "\u03c1w", tex = "\\rho_{\\text{within}}",
                desc = sprintf("the within-site correlation between %s and %s", ylab, xlab),
                eq = "$$Y_{mis}=\\beta_{0m}+u_{ms}+e_{mis},\\quad \\mathbf u_s\\sim N(0,\\Sigma_u),\\ \\mathbf e_{is}\\sim N(0,\\Sigma_e),\\quad \\rho_{\\text{within}}=\\mathrm{corr}_{\\Sigma_e}(Y_1,Y_2)$$"),
    change = list(sym = "\u03b2t", tex = "\\beta_t",
                  desc = sprintf("the average change in %s per unit time (SD units)", ylab),
                  eq = "$$Y_{ist}=\\beta_0+\\beta_t t^*+X_{ist}^\\top\\beta+u_s+b_{0i}+b_{1i}t^*+\\varepsilon_{it},\\quad \\varepsilon_{it}=\\phi\\varepsilon_{i,t-1}+\\omega_{it}$$"),
    groupchange = list(sym = "\u03b2xt", tex = "\\beta_{x\\times t}",
                       desc = sprintf("the difference in change per unit time in %s for %s (SD units)", ylab, xlab),
                       eq = "$$Y_{ist}=\\beta_0+\\beta_t t^*+\\beta_x x_i+\\beta_{x\\times t}x_it^*+u_s+b_{0i}+b_{1i}t^*+\\varepsilon_{it},\\quad \\varepsilon_{it}=\\phi\\varepsilon_{i,t-1}+\\omega_{it}$$"))
  op <- switch(side, two = "\u2260", greater = ">", less = "<")
  c(par, h0 = paste(par$sym, "= 0"), h1 = paste(par$sym, op, "0"))
}

p_side <- function(t, df, side) {
  df <- ifelse(is.na(df) | !is.finite(df), 1e6, df)
  switch(side, two = 2 * pt(-abs(t), df), greater = pt(t, df, lower.tail = FALSE), less = pt(t, df))
}
crit <- function(alpha, side, df) {
  df <- ifelse(is.na(df) | !is.finite(df), 1e6, df)
  qt(1 - if (side == "two") alpha / 2 else alpha, df)
}

## Standard error and degrees of freedom of the tested parameter for a planned design
plan_se <- function(ht, S, n, icc, tau, miss, visits, gap, rel, sb1, taut) {
  n <- n * (1 - miss / 100)
  if (ht %in% c("change", "groupchange")) {
    if (visits < 2) return(list(se = NA, df = NA, msg = "Change over time needs at least two visits."))
    tt <- (seq_len(visits) - 1) * gap; sst <- sum((tt - mean(tt))^2)
    occ <- max(1 - rel, 1e-6)
    if (ht == "change") return(list(se = sqrt((sb1 + occ / sst) / (S * n) + taut^2 / S),
                                    df = if (taut > 0) S - 1 else S * n - 1, msg = NULL))
    return(list(se = sqrt((taut^2 + 4 * (sb1 + occ / sst) / n) / S),
                df = if (taut > 0) S - 2 else S * n - S - 1, msg = NULL))
  }
  switch(ht,
    group = list(se = sqrt((tau^2 + 4 * (1 - icc) / n) / S), df = if (tau > 0) S - 2 else S * n - S - 1, msg = NULL),
    sitelevel = list(se = sqrt(4 * (icc + (1 - icc) / n) / S), df = S - 2, msg = NULL),
    corr = list(se = 1 / sqrt(max(S * n - S - 3, 1)), df = Inf, msg = NULL))
}
plan_power <- function(ht, delta, alpha, side, ...) {
  s <- plan_se(ht, ...)
  if (is.na(s$se)) return(c(s, power = NA))
  d <- if (ht == "corr") atanh(min(abs(delta), 0.99)) else abs(delta)
  c(s, power = pt(d / s$se - crit(alpha, side, s$df), ifelse(is.finite(s$df), max(s$df, 1), 1e6)))
}
plan_required <- function(ht, delta, alpha, side, target, S, n, icc, tau, miss, visits, gap, rel, sb1, taut) {
  pw <- function(S, n) plan_power(ht, delta, alpha, side, S, n, icc, tau, miss, visits, gap, rel, sb1, taut)$power
  rs <- NA; for (s in 3:2000) { p <- pw(s, n); if (!is.na(p) && p >= target) { rs <- s; break } }
  rn <- NA; for (m in 2:5000) { p <- pw(S, m); if (!is.na(p) && p >= target) { rn <- m; break } }
  list(S = rs, n = rn)
}
mde <- function(se, df, alpha, side, power, corr = FALSE) {
  dfx <- ifelse(is.na(df) | !is.finite(df), 1e6, max(df, 1))
  d <- se * (crit(alpha, side, df) + qt(power, dfx))
  if (corr) tanh(d) else d
}

## Test of the stated hypothesis on the collected data, with the model chosen by the decisions
fit_hypothesis <- function(D, ht, y, x, dc, f, alpha, side) {
  out <- list(ok = FALSE)
  d <- D$d_base
  if (ht %in% c("group", "sitelevel")) {
    if (is.null(x) || !x %in% names(d)) { out$why <- "The exposure is not a column in the loaded data, so the hypothesis was not tested."; return(out) }
    dd <- d[!is.na(d[[y]]), ]; dd$site <- droplevels(dd$site)
    st <- slope_test(dd, y, x)
    use_slope <- ht == "group" && dc$slope && !is.na(st$tau)
    est <- if (use_slope) st$est else st$est_ri; se <- if (use_slope) st$se else st$se_ri
    df <- if (use_slope) st$df else st$df_ri
    out <- list(ok = TRUE, est = est, se = se, df = df, p = p_side(est / se, df, side),
                lo = est - qt(1 - alpha / 2, df) * se, hi = est + qt(1 - alpha / 2, df) * se,
                model = if (use_slope) sprintf("random intercept and random slope for %s by site", x) else "random site intercept",
                tau = st$tau, N = st$N, S = st$S, units = "SD units", pi = NULL)
    if (use_slope) {
      hw <- qt(0.975, max(1, st$S - 2)) * sqrt(se^2 + st$tau^2); out$pi <- c(est - hw, est + hw)
    }
    return(out)
  }
  if (ht == "corr") {
    if (is.null(x) || !x %in% names(d)) { out$why <- "The second measure is not a column in the loaded data, so the hypothesis was not tested."; return(out) }
    cc <- d[complete.cases(d[, c(y, x)]), ]; cc$site <- droplevels(cc$site)
    S <- nlevels(cc$site); N <- nrow(cc)
    w1 <- cc[[y]] - ave(cc[[y]], cc$site); w2 <- cc[[x]] - ave(cc[[x]], cc$site)
    r <- cor(w1, w2); dfw <- N - S - 3
    z <- atanh(r) * sqrt(dfw)
    out <- list(ok = TRUE, est = r, se = 1 / sqrt(dfw), df = Inf,
                p = switch(side, two = 2 * pnorm(-abs(z)), greater = pnorm(z, lower.tail = FALSE), less = pnorm(z)),
                lo = tanh(atanh(r) - qnorm(1 - alpha / 2) / sqrt(dfw)), hi = tanh(atanh(r) + qnorm(1 - alpha / 2) / sqrt(dfw)),
                model = "within-site correlation after removing site means", N = N, S = S, units = "correlation",
                pooled = cor(cc[[y]], cc[[x]]),
                between = if (S > 3) cor(tapply(cc[[y]], cc$site, mean), tapply(cc[[x]], cc$site, mean)) else NA)
    return(out)
  }
  ## change over time
  if (!isTRUE(D$longitudinal)) {
    out$why <- "The loaded data have one assessment per participant, so change over time cannot be estimated and the hypothesis was not tested. Map 'Time since baseline' on the Data tab for longitudinal data."
    return(out)
  }
  dl <- D$d[!is.na(D$d[[y]]) & !is.na(D$d$time), ]
  if (ht == "groupchange" && (is.null(x) || !x %in% names(dl))) { out$why <- "No group variable is available, so the group difference in change was not tested."; return(out) }
  base <- dl[order(dl$id, dl$time), ]; base <- base[!duplicated(base$id), ]
  dl$yy <- (dl[[y]] - mean(base[[y]])) / sd(base[[y]])
  dl <- dl[order(dl$site, dl$id, dl$time), ]; dl$site <- droplevels(factor(dl$site)); dl$id <- factor(dl$id)
  fx <- if (ht == "groupchange") as.formula(paste("yy ~ time *", x)) else yy ~ time
  rs <- if (isTRUE(f$int_differ)) ~ 1 + time else ~ 1
  ctl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200, returnObject = TRUE)
  m1 <- nlme::lme(fx, random = list(site = rs, id = ~ 1 + time), data = dl, method = "REML", control = ctl)
  m2 <- tryCatch(update(m1, correlation = corCAR1(form = ~ time | site / id)), error = function(e) NULL)
  m <- if (!is.null(m2) && AIC(m2) < AIC(m1)) m2 else m1
  tt <- summary(m)$tTable
  term <- if (ht == "groupchange") rownames(tt)[grepl(":", rownames(tt))][1] else "time"
  est <- tt[term, 1]; se <- tt[term, 2]; df <- tt[term, 3]
  list(ok = TRUE, est = est, se = se, df = df, p = p_side(est / se, df, side),
       lo = est - qt(1 - alpha / 2, df) * se, hi = est + qt(1 - alpha / 2, df) * se,
       model = paste0("participant random intercept and slope, site random ", if (isTRUE(f$int_differ)) "intercept and slope for time" else "intercept",
                      if (identical(m, m2)) ", continuous-time AR(1) residuals (lower AIC)" else ", independent residuals (lower AIC)"),
       N = length(unique(dl$id)), S = nlevels(dl$site), units = "SD units per unit time")
}

## Decision path as a responsive HTML flowchart (text wraps to the available width)
dg_flow <- function(tab) {
  rows <- lapply(seq_len(nrow(tab)), function(i) {
    r <- tab[i, ]
    col <- if (r$active) LAYER_COL[[r$layer]] else "#9CA3AF"
    tint <- if (r$active) LAYER_TINT[[r$layer]] else "#F1F3F6"
    tagList(
      tags$div(class = "flow-row",
        tags$div(class = "flow-step", paste("Step", r$step)),
        tags$div(class = "flow-q", style = sprintf("border-color:%s", col),
                 tags$span(class = "flow-dia", style = sprintf("background:%s", col)), r$question),
        tags$div(class = "flow-arrow", style = sprintf("color:%s", col), HTML("&#10140;")),
        tags$div(class = paste("flow-a", if (!r$active) "inactive"), style = sprintf("border-color:%s;background:%s", col, tint),
                 tags$div(class = "flow-ans", style = sprintf("color:%s", if (r$active) col else "#4B5563"), r$answer),
                 tags$div(class = "flow-act", if (r$active) r$action else paste(r$action, "\u2014 not needed")),
                 if (nzchar(r$warn) && r$active) tags$div(class = "flow-warn", HTML("&#9888; "), r$warn))),
      if (i < nrow(tab)) tags$div(class = "flow-conn", HTML("&#8595;")))
  })
  tags$div(class = "flow", rows,
           tags$div(class = "flow-legend",
                    lapply(names(LAYER_COL), function(n) tags$span(class = "flow-key",
                      tags$span(class = "flow-dot", style = sprintf("background:%s", LAYER_COL[[n]])), n)),
                    tags$span(class = "flow-key", tags$span(class = "flow-dot", style = "background:#9CA3AF"), "Not needed")))
}

## Decision path figure
dg_plot <- function(tab) {
  n <- nrow(tab); y <- -(seq_len(n) - 1) * 1.6
  wrap <- function(s, w) vapply(s, function(z) paste(strwrap(z, w), collapse = "\n"), "")
  tab$y <- y
  tab$fill <- ifelse(tab$active, LAYER_TINT[tab$layer], "#F1F3F6")
  tab$col <- ifelse(tab$active, LAYER_COL[tab$layer], "#9CA3AF")
  dia <- do.call(rbind, lapply(seq_len(n), function(i)
    data.frame(id = i, x = c(-3.2, 0, 3.2, 0), yy = y[i] + c(0, 0.62, 0, -0.62), col = tab$col[i])))
  ggplot() +
    geom_segment(data = data.frame(y = y[-n], yend = y[-1] + 0.62), aes(x = 0, xend = 0, y = y - 0.62, yend = yend),
                 colour = "grey55", arrow = arrow(length = unit(0.18, "cm"))) +
    geom_segment(data = tab, aes(x = 3.2, xend = 4.6, y = y, yend = y, colour = I(col)), linewidth = 0.9,
                 arrow = arrow(length = unit(0.2, "cm"))) +
    geom_polygon(data = dia, aes(x, yy, group = id, colour = I(col)), fill = "white", linewidth = 0.9) +
    geom_text(data = tab, aes(0, y, label = wrap(question, 22)), size = 3.7, lineheight = 0.92, colour = "#0F172A", fontface = "bold") +
    geom_rect(data = tab, aes(xmin = 4.6, xmax = 12.4, ymin = y - 0.74, ymax = y + 0.74, fill = I(fill), colour = I(col)),
              linewidth = 0.8) +
    geom_text(data = tab, aes(8.5, y + 0.47, label = wrap(answer, 38), colour = I(ifelse(active, col, "#4B5563"))),
              size = 3.4, fontface = "bold", lineheight = 0.9) +
    geom_text(data = tab, aes(8.5, y - 0.05, label = wrap(action, 42)),
              size = 3.15, lineheight = 0.92, colour = ifelse(tab$active, "#0F172A", "#4B5563")) +
    geom_text(data = tab[nzchar(tab$warn) & tab$active, ], aes(8.5, y - 0.6, label = "\u26a0 warning in the decision list"),
              size = 3.0, colour = "#C2414F", fontface = "bold") +
    geom_text(data = tab, aes(-3.5, y, label = paste("Step", step)), hjust = 1, size = 3.6, colour = "#334155", fontface = "bold") +
    scale_x_continuous(limits = c(-5.6, 12.6)) +
    theme_void()
}

eq <- function(...) withMathJax(helpText(paste0(...)))
note <- function(...) tags$p(class = "text-muted small", ...)

## -----------------------------------------------------------------------------
##  Decision panel helpers (shown beside each inspection)
## -----------------------------------------------------------------------------

dec_item <- function(status, title, decision, evidence = NULL, rule = NULL) {
  lab <- c(act = "Action", ok = "Fine", warn = "Check", stop = "Avoid")[[status]]
  tags$div(class = paste("dcard", status),
           tags$div(class = "dhead", tags$span(class = paste("pill", status), lab), tags$span(class = "dtitle", title)),
           tags$div(class = "ddec", decision),
           if (!is.null(evidence)) tags$div(class = "devi", evidence),
           if (!is.null(rule)) tags$div(class = "drule", icon("scale-balanced"), " ", rule))
}
dec_empty <- function(msg) tags$div(class = "dempty", icon("hourglass-half"), " ", msg)
decision_card <- function(id, title = "Decision") {
  card(class = "decision-card", fill = FALSE,
       card_header(tags$span(icon("compass"), " ", title)),
       card_body(fillable = FALSE, uiOutput(id)))
}
split_view <- function(..., dec, title = "Decision") {
  layout_columns(col_widths = c(8, 4), tags$div(...), decision_card(dec, title))
}
vbt <- function(bg) value_box_theme(bg = bg, fg = "#1E293B")

APP_CSS <- "
body{background:#F5F7FB;}
.navbar{background:linear-gradient(90deg,#E9F1FB 0%,#F2ECF8 55%,#E8F5EE 100%)!important;border-bottom:1px solid #DCE3EE;}
.navbar-brand{font-weight:700;color:#1E3A5F!important;letter-spacing:.2px}
.navbar .nav-link{font-weight:500;color:#1E293B!important;padding:.85rem .45rem!important;font-size:.88rem}
.navbar .nav-link.active{color:#2F5FA7!important;box-shadow:inset 0 -3px 0 #2F5FA7}
.card{border:1px solid #E3E9F2;border-radius:14px;box-shadow:0 1px 2px rgba(30,41,59,.05),0 8px 20px rgba(30,41,59,.04);}
.card-header{background:#FFFFFF;border-bottom:1px solid #EEF2F7;font-weight:600;color:#1E3A5F;border-top-left-radius:14px!important;border-top-right-radius:14px!important}
.bslib-sidebar-layout>.sidebar{background:#FBFCFE;border-right:1px solid #E6ECF4}
.bslib-sidebar-layout{border-radius:14px;border:1px solid #E3E9F2;background:#F5F7FB}
.btn-primary{border-radius:10px;font-weight:600;box-shadow:0 2px 6px rgba(47,95,167,.25)}
.btn-default,.btn-outline-primary,.shiny-download-link{border-radius:10px}
.bslib-value-box{border-radius:14px!important;border:1px solid #E3E9F2!important}
.bslib-value-box .value-box-title{font-size:.95rem;color:#1E293B;font-weight:500}
.bslib-value-box .value-box-value{font-weight:700}
.hero{background:linear-gradient(120deg,#E9F1FB 0%,#F3EDF9 55%,#E8F5EE 100%);border:1px solid #DCE6F3;border-radius:20px;padding:28px 34px;margin:6px 0 18px 0}
.hero h1{font-family:'Source Serif 4',Georgia,serif;font-weight:700;color:#1E3A5F;font-size:2rem;margin-bottom:.4rem}
.hero .lead{color:#334155;max-width:70ch}
.case{background:#FFFFFF;border-radius:14px;border:1px solid #E3E9F2;padding:14px 18px;margin-top:14px;max-width:80ch}
.case .num{font-weight:700;color:#C2410C}
.tile{border-radius:16px;padding:16px 18px;height:100%;border:1px solid}
.tile h5{font-weight:700;margin-bottom:.3rem}
.tile .big{font-size:1.7rem;font-weight:700;line-height:1.1;margin:.3rem 0}
.tile p{margin-bottom:0;color:#1E293B;font-size:.98rem}
.t-site{background:#E8F0FA;border-color:#C9DAF0}.t-site h5,.t-site .big{color:#2F5FA7}
.t-measure{background:#F4EAF3;border-color:#E3CFE1}.t-measure h5,.t-measure .big{color:#8E4F87}
.t-time{background:#E8F5E6;border-color:#C9E5C4}.t-time h5,.t-time .big{color:#3F7F36}
.t-infer{background:#FDF0E2;border-color:#F3D9B7}.t-infer h5,.t-infer .big{color:#B5631B}
.stage{border-radius:16px;padding:16px 18px;background:#FFFFFF;border:1px solid #E3E9F2;height:100%}
.stage .n{display:inline-block;width:28px;height:28px;border-radius:50%;background:#E8F0FA;color:#2F5FA7;text-align:center;font-weight:700;line-height:28px;margin-right:8px}
.decision-card{position:sticky;top:72px;align-self:start;height:auto!important;border-top:4px solid #2F5FA7!important}
.decision-card .card-body{max-height:calc(100vh - 140px);overflow-y:auto}
.dcard{border-radius:12px;padding:10px 12px;margin-bottom:10px;border:1px solid}
.dcard.act{background:#EAF2FC;border-color:#C6D9F2}
.dcard.ok{background:#EAF7EF;border-color:#C3E4CF}
.dcard.warn{background:#FDF4E4;border-color:#F2D9AA}
.dcard.stop{background:#FCECEE;border-color:#F1C4CA}
.dhead{display:flex;align-items:center;gap:8px;margin-bottom:4px}
.dtitle{font-weight:600;color:#0F172A;font-size:1rem}
.pill{display:inline-block;padding:2px 10px;border-radius:999px;font-size:.8rem;font-weight:700;white-space:nowrap}
.pill.act{background:#2F5FA7;color:#FFF}.pill.ok{background:#2E8B57;color:#FFF}
.pill.warn{background:#D98E1F;color:#FFF}.pill.stop{background:#C2414F;color:#FFF}
.ddec{color:#0F172A;font-size:.98rem;line-height:1.45}
.devi{color:#334155;font-size:.92rem;margin-top:5px;line-height:1.4}
.drule{color:#334155;font-size:.88rem;margin-top:5px;font-style:italic}
.dempty{color:#334155;font-size:.98rem;padding:10px;background:#F8FAFC;border-radius:10px;border:1px dashed #CBD5E1}
.shiny-notification{border-left:4px solid #2F5FA7;border-radius:10px}
pre{border-radius:10px}
body{font-size:1rem;color:#0F172A}
.text-muted{color:#475569!important}
p.small,.small{font-size:.92rem!important}
.help-block,.shiny-text-output,.form-text{color:#1E293B}
.help-block{font-size:1.02rem}
.MathJax_Display,.MathJax{font-size:112%!important;color:#0F172A}
pre.shiny-text-output,pre{font-size:.9rem;line-height:1.5;color:#0F172A;background:#F8FAFC;border:1px solid #E2E8F0;padding:12px;white-space:pre-wrap;word-break:break-word}
table.table,.table{font-size:.95rem;color:#0F172A}
.table th{color:#1E3A5F;font-weight:600}
table.dataTable{font-size:.93rem!important;color:#0F172A}
.control-label,label{color:#1E293B;font-weight:500}
.card-header{font-size:1rem}
.flow{display:flex;flex-direction:column}
.flow-row{display:grid;grid-template-columns:62px minmax(150px,1fr) 26px minmax(200px,1.55fr);align-items:center;gap:8px}
.flow-step{font-weight:700;color:#334155;font-size:.9rem;text-align:right}
.flow-q{background:#FFFFFF;border:2px solid;border-radius:14px;padding:10px 12px;font-weight:600;color:#0F172A;font-size:.98rem;line-height:1.35;display:flex;gap:8px;align-items:flex-start}
.flow-dia{flex:0 0 auto;width:11px;height:11px;transform:rotate(45deg);margin-top:6px;border-radius:2px}
.flow-arrow{font-size:1.3rem;text-align:center;font-weight:700}
.flow-a{border:2px solid;border-radius:14px;padding:9px 12px}
.flow-a.inactive{border-style:dashed}
.flow-ans{font-weight:700;font-size:.98rem;margin-bottom:3px}
.flow-act{color:#0F172A;font-size:.95rem;line-height:1.4}
.flow-a.inactive .flow-act{color:#4B5563}
.flow-warn{color:#B4232F;font-size:.88rem;margin-top:5px;font-weight:500}
.flow-conn{margin-left:calc(62px + 8px);width:calc((100% - 62px - 26px - 24px) * 0.392);text-align:center;color:#94A3B8;font-size:1.2rem;line-height:1.1;margin-top:2px;margin-bottom:2px}
.flow-legend{display:flex;flex-wrap:wrap;gap:14px;margin-top:14px;padding-top:10px;border-top:1px solid #E2E8F0;font-size:.9rem;color:#1E293B}
.flow-key{display:inline-flex;align-items:center;gap:6px}
.flow-dot{width:12px;height:12px;border-radius:3px;display:inline-block}
"

thm <- bs_theme(version = 5, bg = "#FFFFFF", fg = "#1E293B", primary = "#2F5FA7",
                secondary = "#64748B", success = "#2E8B57", warning = "#D98E1F", danger = "#C2414F",
                base_font = font_link("Inter", "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"),
                heading_font = font_link("Source Serif 4", "https://fonts.googleapis.com/css2?family=Source+Serif+4:wght@600;700&display=swap"),
                code_font = font_collection("JetBrains Mono", "Menlo", "Consolas", "monospace"),
                "border-radius" = "0.6rem", "font-size-base" = "0.95rem")

ui <- page_navbar(
  title = tags$span(icon("diagram-project"), " HierDep"),
  theme = thm,
  fillable = FALSE,
  header = tagList(withMathJax(), tags$style(HTML(APP_CSS))),

  ## ---------------------------------------------------------------- Overview
  nav_panel(
    "Overview", icon = icon("house"),
    tags$div(class = "hero",
      tags$h1("Pooled multi-site data are not independent"),
      tags$p(class = "lead",
             "Participants at the same site share its scanner, examiners, instruments and referral routes; measures on the same participant share that participant; repeated visits share both. ",
             "HierDep Explorer estimates each layer of dependence, shows what it does to inference, and turns the result into a decision at every step, from planning a study to testing its hypothesis."),
      tags$div(class = "case",
        tags$b("A case from ABIDE I. "),
        "Eye status at scan was set by the site in 18 of 20 cohorts. A single-level regression gave eyes closed an effect of ",
        tags$span(class = "num", "0.56 SD on DVARS (p < 0.001)"), ". With a random site intercept the effect was ",
        tags$span(class = "num", "\u22120.01 SD (p = 0.92)"),
        ". The single-level analysis had compared scanners, not eye states: its 873 observations behaved like about 20.")),
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      tags$div(class = "tile t-site", tags$h5(icon("building"), " Site layer"),
               tags$div(class = "big", "ICC 0.00 \u2013 0.57"),
               tags$p("From parent-reported ADI-R social to DVARS. Design effects reach 25.3: 873 DVARS scans carry the information of 35 participants.")),
      tags$div(class = "tile t-measure", tags$h5(icon("layer-group"), " Measure layer"),
               tags$div(class = "big", "\u22120.09 vs +0.21"),
               tags$p("Head motion and DVARS correlate negatively within sites and positively between them; the pooled value reflects neither.")),
      tags$div(class = "tile t-time", tags$h5(icon("clock"), " Time layer"),
               tags$div(class = "big", "0.80 \u00b7 0.75 \u00b7 0.64"),
               tags$p("Correlations between visits built from stable person differences, differences in change and carry-over between occasions.")),
      tags$div(class = "tile t-infer", tags$h5(icon("triangle-exclamation"), " Inference"),
               tags$div(class = "big", "30% \u2013 86%"),
               tags$p("False-positive rates for site-level exposures when site is ignored, against 4% to 6% with a random site intercept."))
    ),
    tags$br(),
    layout_columns(
      col_widths = c(4, 4, 4),
      tags$div(class = "stage", tags$h5(tags$span(class = "n", "1"), "Plan the study"),
               tags$p("State the objective and hypothesis, set the assumed site dependence and the design, and obtain power, the sites needed, the decision path and the analysis plan in the ", tags$b("Decision guide"), ".")),
      tags$div(class = "stage", tags$h5(tags$span(class = "n", "2"), "Inspect the collected data"),
               tags$p("Load the data and work through the layers. Every tab shows the inspection on the left and the decision it leads to on the right.")),
      tags$div(class = "stage", tags$h5(tags$span(class = "n", "3"), "Test and report"),
               tags$p("Fit the model the decisions select, test the hypothesis, compare the plan with what was observed, and download the script, the decision record and the reporting text."))
    ),
    tags$br(),
    layout_columns(
      col_widths = c(7, 5),
      card(card_header("The model behind the decisions"),
           eq("$$Y_{is} = \\mu + u_s + \\varepsilon_{is},\\quad u_s\\sim N(0,\\sigma^2_u),\\ \\varepsilon_{is}\\sim N(0,\\sigma^2_\\varepsilon)\\qquad (1)$$",
              "$$\\mathrm{ICC} = \\frac{\\sigma^2_u}{\\sigma^2_u+\\sigma^2_\\varepsilon}\\qquad (2)\\qquad\\qquad \\mathrm{DEFF} = 1+(\\bar n-1)\\,\\mathrm{ICC},\\quad N_{\\mathrm{eff}} = N/\\mathrm{DEFF}\\qquad (4)$$",
              "$$Y_{mist}=\\beta_{0m}+\\beta_{1m}D_i+\\beta_{2m}(x_{ist}-\\bar x_s)+\\beta_{3m}\\bar x_s+\\gamma_m^\\top z_s+\\beta_{tm}t^*+u_{ms}+v_{ms}D_i+w_{m,k(i)}+b_{0mi}+b_{1mi}t^*+\\varepsilon_{mit}\\qquad(12)$$"),
           note("Each term enters only when the data support it; the tabs and the decision guide show which.")),
      card(card_header("Structure of pooled data"), plotOutput("hier_plot", height = "400px"))
    )
  ),

  ## ---------------------------------------------------------------- Study decision guide
  nav_panel(
    "Decision guide", icon = icon("compass"),
    layout_sidebar(
      sidebar = sidebar(
        width = 390,
        radioButtons("stage", "Stage",
                     c("1. Planning the study (before data are collected)" = "plan",
                       "2. Analysing the collected data" = "analyse")),
        accordion(
          id = "acc_h", open = "Objective and hypothesis",
          accordion_panel("Objective and hypothesis",
            textAreaInput("h_obj", "Primary objective",
                          "To estimate the difference in full-scale IQ between autistic and non-autistic participants in pooled multi-site data.",
                          rows = 3),
            selectInput("h_type", "Hypothesis", HYP_TYPES),
            textInput("h_ylab", "Primary outcome", "Full-scale IQ"),
            textInput("h_xlab", "Exposure, comparison or second measure", "autistic versus non-autistic participants"),
            radioButtons("h_side", "Alternative hypothesis",
                         c("Two-sided" = "two", "Greater than 0" = "greater", "Less than 0" = "less"), inline = TRUE),
            numericInput("h_alpha", "Significance level \u03b1", 0.05, min = 0.001, max = 0.2, step = 0.005),
            numericInput("h_delta", "Effect of interest (SD units; a correlation for the relationship hypothesis; change per unit time for change hypotheses)",
                         0.4, step = 0.05))
        ),
        conditionalPanel("input.stage == 'plan'",
          accordion(
            open = "Planned design",
            accordion_panel("Planned design",
              numericInput("p_power", "Target power", 0.8, min = 0.5, max = 0.99, step = 0.05),
              numericInput("p_S", "Number of sites", 12, min = 2),
              numericInput("p_n", "Participants per site", 40, min = 2),
              numericInput("p_miss", "Expected missing values or dropout (%)", 10, min = 0, max = 90)),
            accordion_panel("Assumed site dependence",
              selectInput("p_preset", "Site ICC from Table 2",
                          c("(enter a value)" = "", setNames(seq_len(nrow(TABLE2)), paste0(TABLE2$dataset, ": ", TABLE2$measure, " (", f3(TABLE2$icc), ")")))),
              numericInput("p_icc", "Site ICC of the outcome", 0.05, min = 0, max = 0.99, step = 0.005),
              numericInput("p_tau", "Between-site SD of the effect \u03c4 (SD units)", 0.2, min = 0, step = 0.05)),
            accordion_panel("Measures and visits",
              numericInput("p_nmeas", "Measures analysed together", 1, min = 1),
              checkboxInput("p_predict", "The measures are combined in a prediction model", FALSE),
              numericInput("p_visits", "Visits per participant", 1, min = 1),
              numericInput("p_gap", "Time between visits", 1, min = 0.1, step = 0.5),
              numericInput("p_rel", "Test-retest reliability of the outcome", 0.75, min = 0.1, max = 0.99, step = 0.05),
              numericInput("p_sb1", "Between-person variance of change (SD units\u00b2)", 0.02, min = 0, step = 0.01),
              numericInput("p_taut", "Between-site SD of change (SD units)", 0, min = 0, step = 0.01),
              numericInput("p_fu", "Sites with follow-up by design", NA, min = 0),
              checkboxInput("p_intdiff", "Visit intervals, protocols or scanners differ between sites", FALSE)),
            accordion_panel("Site composition",
              checkboxInput("p_instr", "Different instruments across sites or participants", FALSE),
              checkboxInput("p_inst", "Several cohorts from the same institution", FALSE),
              numericInput("p_misssites", "Sites not recording the outcome", 0, min = 0)))),
        conditionalPanel("input.stage == 'analyse'",
          tags$p(class = "small", "Load data on the Data tab, choose the variables, fill the answers from the data and test the hypothesis. Every answer can be edited and the decisions update at once."),
          selectizeInput("dg_y", "Outcome", choices = NULL, options = list(create = TRUE)),
          selectizeInput("dg_x", "Exposure, comparison or second measure", choices = NULL, options = list(create = TRUE)),
          selectizeInput("dg_joint", "Other measures analysed together with the outcome", choices = NULL, multiple = TRUE,
                         options = list(create = TRUE)),
          checkboxInput("dg_predict", "The measures are combined in a prediction model", FALSE),
          actionButton("dg_fill", "Fill answers from the loaded data", class = "btn-primary w-100"),
          tags$br(), tags$br(),
          actionButton("dg_test", "Test the hypothesis with the chosen model", class = "btn-primary w-100"),
          tags$br(), tags$br(),
          accordion(
            open = FALSE,
            accordion_panel("Site structure",
              numericInput("dg_S", "Number of sites", 20, min = 1),
              numericInput("dg_nbar", "Mean participants per site", 44, min = 1),
              numericInput("dg_icc", "Site ICC of the outcome", 0.05, min = 0, max = 1, step = 0.005),
              numericInput("dg_icc_lo", "ICC lower 95% limit", NA, min = 0, max = 1, step = 0.005),
              numericInput("dg_icc_hi", "ICC upper 95% limit", NA, min = 0, max = 1, step = 0.005),
              checkboxInput("dg_singular", "Empty model is singular (site variance estimated at zero)", FALSE)),
            accordion_panel("Exposure",
              numericInput("dg_rhox", "\u03c1x: share of the exposure's variance between sites", 0, min = 0, max = 1, step = 0.01),
              note("0.8 or more: site property; 0.1 or less: participant property; between: varies within and between sites."),
              numericInput("dg_tau", "\u03c4: between-site SD of the effect (SD units)", NA, step = 0.01),
              numericInput("dg_ptau", "p for \u03c4\u00b2 = 0 (likelihood ratio)", NA, min = 0, max = 1, step = 0.001),
              numericInput("dg_est", "Average effect (SD units)", NA, step = 0.01),
              numericInput("dg_se", "Standard error of the average effect", NA, step = 0.001)),
            accordion_panel("What site is made of",
              checkboxInput("dg_instr", "The instrument differs between participants (crossed with site)", FALSE),
              checkboxInput("dg_instmulti", "Some institutions contribute more than one cohort", FALSE),
              selectInput("dg_instlevel", "Where the site variance sits",
                          c("Not checked" = "unchecked", "Institution" = "institution",
                            "Cohort" = "cohort", "Both levels" = "nested"))),
            accordion_panel("Measures, time and missing data",
              numericInput("dg_nmeas", "Number of measures analysed together", 1, min = 1),
              checkboxInput("dg_long", "Repeated visits", FALSE),
              numericInput("dg_fu", "Sites with follow-up by design", NA, min = 0),
              checkboxInput("dg_irr", "Visit times irregular", FALSE),
              checkboxInput("dg_intdiff", "Visit intervals, protocols or scanners differ between sites", FALSE),
              numericInput("dg_misssites", "Sites where the outcome was not recorded", 0, min = 0),
              numericInput("dg_misspct", "% missing within contributing sites", 0, min = 0, max = 100))))
      ),
      card(card_header("Objective and hypothesis"), uiOutput("h_card")),

      ## ------------------------------------------------ planning stage
      conditionalPanel("input.stage == 'plan'",
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          value_box("Power of the planned design", textOutput("p_vb_power"), theme = vbt("#E8F0FA")),
          value_box("Sites needed for the target power", textOutput("p_vb_S"), theme = vbt("#F4EAF3")),
          value_box("Participants per site needed", textOutput("p_vb_n"), theme = vbt("#E8F5E6")),
          value_box("Total participants planned", textOutput("p_vb_N"), theme = vbt("#FDF0E2"))),
        layout_columns(
          col_widths = c(7, 5),
          card(card_header("Power against the number of sites"), plotOutput("p_pow_plot", height = "460px")),
          card(card_header("Basis of the calculation"), uiOutput("p_pow_text"))),
        layout_columns(
          col_widths = c(7, 5),
          card(card_header("Decision path for the planned analysis"), full_screen = TRUE, uiOutput("dg_plot_plan")),
          card(class = "decision-card", fill = FALSE, card_header(tags$span(icon("compass"), " Decisions for the planned analysis")), card_body(fillable = FALSE, uiOutput("dg_list_plan")))),
        layout_columns(
          col_widths = c(6, 6),
          card(card_header("Statistical analysis plan"), uiOutput("p_sap"),
               downloadButton("p_dl_sap", "Download analysis plan text")),
          card(card_header("Planned analysis script"), full_screen = TRUE,
               verbatimTextOutput("p_code"),
               downloadButton("p_dl_code", "Download planned analysis script (.R)")))
      ),

      ## ------------------------------------------------ analysis stage
      conditionalPanel("input.stage == 'analyse'",
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          value_box("Analysis", textOutput("dg_vb_ml"), theme = vbt("#E8F0FA")),
          value_box("Design effect for the exposure", textOutput("dg_vb_deff"), theme = vbt("#F4EAF3")),
          value_box("Effective sample size", textOutput("dg_vb_neff"), theme = vbt("#E8F5E6")),
          value_box("Exposure sits at", textOutput("dg_vb_xl"), theme = vbt("#FDF0E2"))),
        layout_columns(
          col_widths = c(7, 5),
          card(card_header("Test of the hypothesis"), uiOutput("h_result")),
          card(card_header("Planned design against the collected data"), tableOutput("h_plan_obs"),
               note("Planned values come from stage 1. The minimum detectable effect uses the target power and \u03b1 of the plan."))),
        layout_columns(
          col_widths = c(7, 5),
          card(card_header("Decision path"), full_screen = TRUE, uiOutput("dg_plot")),
          card(class = "decision-card", fill = FALSE, card_header(tags$span(icon("compass"), " Decisions, evidence and warnings")),
               card_body(fillable = FALSE, uiOutput("dg_list"), downloadButton("dg_dl_rec", "Download decision record")))),
        layout_columns(
          col_widths = c(7, 5),
          card(card_header("Analysis script for these decisions"), full_screen = TRUE,
               verbatimTextOutput("dg_code"),
               downloadButton("dg_dl_code", "Download analysis script (.R)")),
          card(card_header("What to report"), uiOutput("dg_report")))
      )
    )
  ),

  ## ---------------------------------------------------------------- Data
  nav_panel(
    "Data", icon = icon("database"),
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        radioButtons("src", "Data source",
                     c("Simulated ABIDE-format demo" = "demo",
                       "ABIDE I phenotypic file (public download)" = "abide_url",
                       "Upload a CSV" = "upload")),
        conditionalPanel("input.src == 'demo'",
                         numericInput("seed", "Simulation seed", 2026, min = 1)),
        conditionalPanel("input.src == 'abide_url'",
                         note("Reads Phenotypic_V1_0b_preprocessed1.csv from the ABIDE Preprocessed release. If the server has no internet access, download the file and use Upload.")),
        conditionalPanel("input.src == 'upload'",
                         fileInput("file", "CSV file", accept = c(".csv", "text/csv")),
                         uiOutput("map_ui")),
        actionButton("load", "Load data", class = "btn-primary w-100"),
        hr(),
        note("ABIDE files are recognised by their column names and processed as in Listing 1: participants with a functional file and no 'fail' rating from any rater; ADOS and ADI-R in autistic participants; mean FD log-transformed and % high-motion volumes log(1 + x)-transformed.")
      ),
      split_view(
      uiOutput("data_status"),
      layout_columns(
        col_widths = c(12, 12),
        card(card_header("Participants by site"), card_body(fillable = FALSE, DTOutput("site_tab"))),
        card(card_header("Availability of each measure by site"),
             plotOutput("avail_plot", height = "680px"),
             note("Filled circle: recorded for 90% or more of eligible participants; open circle: recorded for some; cross: not recorded at the site. Sites ordered by size."))
      ),
        dec = "dec_data")
    )
  ),

  ## ---------------------------------------------------------------- Site layer
  nav_panel(
    "Site layer", icon = icon("building"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        uiOutput("site_meas_ui"),
        radioButtons("ci", "Interval for the ICC",
                     c("F-based (ANOVA, effective cluster size n0)" = "F",
                       "Parametric bootstrap" = "boot")),
        conditionalPanel("input.ci == 'boot'",
                         sliderInput("nsim", "Bootstrap replicates", 100, 1000, 200, step = 100)),
        checkboxInput("show_leap", "Show LEAP values from Table 2", TRUE),
        actionButton("run_site", "Run site layer", class = "btn-primary w-100")
      ),
      split_view(
      card(card_header("Equations"),
           eq("$$Y_{is}=\\beta_0+\\beta_1\\mathrm{Dx}_i+\\beta_2(\\mathrm{Age}_{is}-\\overline{\\mathrm{Age}}_s)+\\beta_3\\overline{\\mathrm{Age}}_s+\\beta_4\\mathrm{Female}_i+u_s+\\varepsilon_{is}\\qquad(3)$$"),
           note("ICC from the empty model (1); adjusted ICC from (3), with ADOS module added for ADOS measures; p from the likelihood-ratio test of \\(\\sigma^2_u = 0\\) (mixture \\(\\chi^2\\)); DEFF and N_eff from (4).")),
      card(card_header("Site intraclass correlation, design effect and effective sample size"),
           card_body(fillable = FALSE, DTOutput("icc_tab"))),
      layout_columns(
        col_widths = c(12, 12),
        card(card_header("Site ICCs with 95% intervals"), plotOutput("icc_plot", height = "640px")),
        card(card_header("Design effect against participants per site"), plotOutput("deff_plot", height = "640px"))
      ),
      card(card_header("Estimated site effects (BLUPs, SD units, 95% interval)"),
           uiOutput("blup_ui"), plotOutput("blup_plot", height = "680px")),
        dec = "dec_site")
    )
  ),

  ## ---------------------------------------------------------------- Exposure level
  nav_panel(
    "Exposure", icon = icon("crosshairs"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        tags$h6("Site-level exposure"),
        uiOutput("expo_ui"),
        actionButton("run_expo", "Test exposure", class = "btn-primary w-100"),
        hr(),
        tags$h6("Real-data null experiment"),
        uiOutput("null_meas_ui"),
        sliderInput("B", "Random assignments per outcome", 20, 400, 100, step = 20),
        numericInput("null_seed", "Seed", 77),
        actionButton("run_null", "Run null experiment", class = "btn-primary w-100"),
        note("Satterthwaite degrees of freedom are used in the null experiment for speed, as in Listing 5.")
      ),
      split_view(
      card(card_header("Where each variable sits: its own intraclass correlation \\(\\rho_x\\)"),
           plotOutput("rho_plot", height = "360px"),
           note("\\(\\rho_x\\) near 1: site property; near 0: participant property; between: varies within and between sites.")),
      layout_columns(
        col_widths = c(12, 12),
        card(card_header("Exposure effect with and without site"),
             tableOutput("expo_tab"), uiOutput("expo_text")),
        card(card_header("Variance inflation by level of the exposure"),
             eq("$$\\frac{\\mathrm{Var}_{\\text{correct}}(\\hat\\beta)}{\\mathrm{Var}_{\\text{naive}}(\\hat\\beta)}\\approx1+(\\bar n-1)\\,\\rho_x\\,\\mathrm{ICC}_Y\\qquad(5)$$"),
             plotOutput("infl_plot", height = "380px"))
      ),
      card(card_header("False-positive rates under a known null"),
           plotOutput("null_plot", height = "480px"), card_body(fillable = FALSE, DTOutput("null_tab"))),
        dec = "dec_expo")
    )
  ),

  ## ---------------------------------------------------------------- Group difference
  nav_panel(
    "Groups", icon = icon("people-arrows"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        uiOutput("dx_meas_ui"),
        actionButton("run_dx", "Fit the four analyses", class = "btn-primary w-100"),
        hr(),
        note("Single-level, fixed site, random site intercept and random intercept plus random slope for the group, with Kenward-Roger degrees of freedom for the mixed models.")
      ),
      split_view(
      card(card_header("Equations"),
           eq("$$Y_{is}=\\beta_0+(\\beta_1+v_s)\\mathrm{Dx}_i+X_{is}^\\top\\beta+u_s+\\varepsilon_{is},\\quad (u_s,v_s)\\sim N_2(0,\\Sigma),\\ \\mathrm{Var}(v_s)=\\tau^2_{\\mathrm{Dx}}\\qquad(6)$$",
              "$$\\hat\\beta_1\\pm t_{0.975,S-2}\\sqrt{\\mathrm{SE}(\\hat\\beta_1)^2+\\hat\\tau^2}\\qquad(16)$$")),
      layout_columns(
        col_widths = c(12, 12),
        card(card_header("Group difference under four analyses (SD units, 95% CI)"),
             plotOutput("dx_plot", height = "440px"), tableOutput("dx_tab"), uiOutput("dx_text")),
        card(card_header("Site-specific differences, pooled difference and prediction interval"),
             plotOutput("ps_plot", height = "640px"),
             note("Blue line and band: random-effects pooled difference and 95% CI; pink band: 95% prediction interval for a new site; point size proportional to site n."),
             uiOutput("ps_text"))
      ),
        dec = "dec_dx")
    )
  ),

  ## ---------------------------------------------------------------- Site composition
  nav_panel(
    "Site make-up", icon = icon("puzzle-piece"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        uiOutput("cx_meas_ui"),
        actionButton("run_cx", "Decompose site", class = "btn-primary w-100")
      ),
      split_view(
      card(card_header("Instrument crossed with site"),
           eq("$$Y_{is}=X_{is}^\\top\\beta+u_s+w_{\\mathrm{test}(i)}+\\varepsilon_{is},\\quad w_{\\mathrm{test}}\\sim N(0,\\sigma^2_w)\\qquad(7)$$"),
           tableOutput("cx_tab"), uiOutput("cx_text")),
      card(card_header("Institution and cohort within institution"),
           plotOutput("inst_plot", height = "360px"), uiOutput("inst_text")),
        dec = "dec_cx")
    )
  ),

  ## ---------------------------------------------------------------- Measure layer
  nav_panel(
    "Measures", icon = icon("layer-group"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        uiOutput("mv_meas_ui"),
        uiOutput("mv_sub_ui"),
        checkboxInput("mv_model", "Fit the multivariate model (8) with nlme", TRUE),
        actionButton("run_mv", "Estimate correlations", class = "btn-primary w-100")
      ),
      split_view(
      card(card_header("Equations"),
           eq("$$Y_{mis}=\\beta_{0m}+X_{is}^\\top\\beta_m+u_{ms}+e_{mis},\\quad \\mathbf u_s\\sim N(0,\\Sigma_u),\\ \\mathbf e_{is}\\sim N(0,\\Sigma_e)\\qquad(8)$$",
              "$$r_{\\text{pooled}}\\approx\\sqrt{\\mathrm{ICC}_1\\mathrm{ICC}_2}\\,r_{\\text{between}}+\\sqrt{(1-\\mathrm{ICC}_1)(1-\\mathrm{ICC}_2)}\\,r_{\\text{within}}\\qquad(9)$$")),
      card(card_header("Within-site, between-site and pooled correlations"),
           plotOutput("mv_plot", height = "500px"), uiOutput("mv_msg"), card_body(fillable = FALSE, DTOutput("mv_tab"))),
        dec = "dec_mv")
    )
  ),

  ## ---------------------------------------------------------------- Time layer
  nav_panel(
    "Time", icon = icon("clock"),
    layout_sidebar(
      sidebar = sidebar(
        width = 330,
        tags$h6("Variance components"),
        numericInput("sb0", "\u03c3\u00b2 b0 (stable between-person)", 12.5, step = 0.5),
        numericInput("sb1", "\u03c3\u00b2 b1 (between-person change)", 0.8, step = 0.1),
        numericInput("sb01", "\u03c3 b01 (intercept-slope covariance)", -2.1, step = 0.1),
        numericInput("se2", "\u03c3\u00b2 \u03b5 (occasion)", 4.2, step = 0.1),
        sliderInput("phi", "\u03c6 (carry-over between visits)", 0, 0.95, 0.35, step = 0.05),
        numericInput("su2", "\u03c3\u00b2 u (site)", 0, step = 0.5),
        textInput("times", "Time since baseline at each visit", "0, 1, 2"),
        hr(),
        tags$h6("Simulate and fit"),
        numericInput("ls_S", "Sites", 5, min = 3, max = 30),
        numericInput("ls_n", "Participants per site", 60, min = 10, max = 300),
        numericInput("ls_bt", "Average change per unit time", 0.5, step = 0.1),
        actionButton("run_ls", "Simulate and fit", class = "btn-primary w-100"),
        uiOutput("long_up_ui")
      ),
      split_view(
      card(card_header("Correlation between visits"),
           eq("$$Y_{ist}=\\beta_0+\\beta_t t^*+X_{ist}^\\top\\beta+u_s+b_{0i}+b_{1i}t^*+\\varepsilon_{it},\\quad (b_{0i},b_{1i})\\sim N_2(0,D),\\ \\varepsilon_{it}=\\phi\\varepsilon_{i,t-1}+\\omega_{it}\\qquad(10)$$",
              "$$\\mathrm{Corr}(Y_{it},Y_{it'})=\\frac{\\sigma^2_{b0}+t^*t'^*\\sigma^2_{b1}+(t^*+t'^*)\\sigma_{b01}+\\sigma^2_\\varepsilon\\phi^{|t-t'|}+\\sigma^2_u}{\\sqrt{\\mathrm{Var}(Y_{it})\\mathrm{Var}(Y_{it'})}}\\qquad(11)$$"),
           layout_columns(col_widths = c(12, 12),
                          tableOutput("vc_tab"),
                          plotOutput("vc_plot", height = "400px"))),
      card(card_header("Simulated longitudinal data: estimated against true components"),
           tableOutput("ls_tab"), uiOutput("ls_text")),
      card(card_header("Longitudinal model for uploaded data"),
           uiOutput("lu_text"), tableOutput("lu_tab")),
        dec = "dec_time")
    )
  ),

  ## ---------------------------------------------------------------- Planning
  nav_panel(
    "Planning", icon = icon("ruler-combined"),
    layout_sidebar(
      sidebar = sidebar(
        width = 330,
        selectInput("preset", "ICC preset (Table 2)",
                    c("(enter values)" = "", setNames(seq_len(nrow(TABLE2)),
                                                      paste0(TABLE2$dataset, ": ", TABLE2$measure, " (ICC ", f3(TABLE2$icc), ")")))),
        numericInput("pl_icc", "Site ICC of the outcome", 0.05, min = 0, max = 0.99, step = 0.005),
        numericInput("pl_n", "Participants per site", 40, min = 2),
        numericInput("pl_S", "Number of sites", 15, min = 2),
        sliderInput("pl_rho", "\u03c1x of the exposure", 0, 1, 1, step = 0.05),
        hr(),
        tags$h6("Group comparison across sites (17)"),
        numericInput("pl_tau", "\u03c4 (between-site SD of the group difference)", 0.41, step = 0.05),
        numericInput("pl_se2", "\u03c3\u00b2 \u03b5 (SD units)", 0.93, step = 0.01)
      ),
      split_view(
      layout_columns(
        col_widths = c(4, 4, 4),
        value_box("Design effect for the exposure", textOutput("vb_deff"), theme = vbt("#E8F0FA")),
        value_box("Total participants", textOutput("vb_N"), theme = vbt("#F4EAF3")),
        value_box("Effective sample size", textOutput("vb_neff"), theme = vbt("#E8F5E6"))
      ),
      card(card_header("Equations"),
           eq("$$\\mathrm{DEFF}_x=1+(\\bar n-1)\\,\\rho_x\\,\\mathrm{ICC}_{\\text{site}},\\quad N_{\\text{eff}}=N/\\mathrm{DEFF}_x\\qquad(15)\\qquad\\qquad \\mathrm{Var}(\\hat\\beta_1)\\approx\\frac{\\tau^2+4\\sigma^2_\\varepsilon/\\bar n}{S}\\qquad(17)$$")),
      layout_columns(
        col_widths = c(12, 12),
        card(card_header("Effective sample size against participants per site"),
             plotOutput("pl_neff", height = "380px")),
        card(card_header("Standard error of the group difference: adding sites or adding participants"),
             plotOutput("pl_se", height = "380px"), uiOutput("pl_text"))
      ),
        dec = "dec_plan")
    )
  ),

  ## ---------------------------------------------------------------- Few sites
  nav_panel(
    "Few sites", icon = icon("circle-half-stroke"),
    layout_sidebar(
      sidebar = sidebar(
        width = 360,
        textAreaInput("eta_in", "Measure, sites (k), participants (N), reported \u03b7\u00b2 — one per line",
                      LEAP_DEFAULT, rows = 9),
        actionButton("run_eta", "Convert", class = "btn-primary w-100"),
        note("Default rows are the six LEAP baseline social-cognition outcomes.")
      ),
      split_view(
      card(card_header("From between-centre \\(\\eta^2\\) to the intraclass correlation"),
           eq("$$F=\\frac{\\eta^2/(k-1)}{(1-\\eta^2)/(N-k)},\\qquad \\mathrm{ICC}=\\frac{F-1}{F-1+\\bar n},\\qquad E(\\eta^2\\mid\\text{no site differences})=\\frac{k-1}{N-1}$$"),
           note("The interval comes from the F distribution with k − 1 and N − k degrees of freedom; negative limits are set to 0.")),
      layout_columns(
        col_widths = c(12, 12),
        card(card_header("Converted values"), tableOutput("eta_tab")),
        card(card_header("ICC with 95% interval"), plotOutput("eta_plot", height = "380px"))
      ),
        dec = "dec_eta")
    )
  ),

  ## ---------------------------------------------------------------- Report
  nav_panel(
    "Report", icon = icon("file-lines"),
    layout_columns(
      col_widths = c(7, 5),
      card(card_header("Methods and results text"),
           uiOutput("report_text"),
           downloadButton("dl_report", "Download text")),
      card(card_header("Result tables"),
           tags$p("Each table is written as CSV from the analyses run in this session."),
           downloadButton("dl_icc", "Site layer (ICC, DEFF, N_eff)"), tags$br(), tags$br(),
           downloadButton("dl_blup", "Site effects (BLUPs)"), tags$br(), tags$br(),
           downloadButton("dl_dx", "Group difference, four analyses"), tags$br(), tags$br(),
           downloadButton("dl_null", "Null experiment"), tags$br(), tags$br(),
           downloadButton("dl_mv", "Within- and between-site correlations"))
    )
  ),
)

## -----------------------------------------------------------------------------
##  Server
## -----------------------------------------------------------------------------

server <- function(input, output, session) {

  D <- reactiveVal(NULL)
  R <- reactiveValues(site = NULL, blup = NULL, dx = NULL, ps = NULL, null = NULL,
                      expo = NULL, mv = NULL, cx = NULL, inst = NULL, ls = NULL, lu = NULL)

  reset_results <- function() for (n in names(R)) R[[n]] <- NULL

  ## ---------------------------------------------------------- overview figure
  output$hier_plot <- renderPlot({
    nodes <- data.frame(
      x = c(2, 6, 1, 3, 5, 7, 0.4, 1, 1.6, 2.4, 3, 3.6, 4.4, 5, 5.6, 6.4, 7, 7.6),
      y = c(3, 3, 2, 2, 2, 2, rep(1, 12)),
      lev = c("Site", "Site", rep("Participant", 4), rep("Measure", 12)),
      lab = c("Site A", "Site B", "P1", "P2", "P3", "P4",
              rep(c("IQ", "AD", "MRI"), 4)))
    edges <- data.frame(x = c(2, 2, 6, 6, rep(c(1, 3, 5, 7), each = 3)),
                        y = c(3, 3, 3, 3, rep(2, 12)),
                        xend = c(1, 3, 5, 7, nodes$x[7:18]),
                        yend = c(2, 2, 2, 2, rep(1, 12)))
    ggplot() +
      geom_segment(data = edges, aes(x, y, xend = xend, yend = yend), colour = "grey65") +
      annotate("rect", xmin = -0.2, xmax = 8.2, ymin = 2.75, ymax = 3.55, fill = "#E3F0FA", alpha = .6) +
      annotate("text", x = 4, y = 3.42, label = "Site-level features: scanner, eye status, examiners, instrument",
               size = 4.2, colour = "#2B6CB0") +
      geom_point(data = nodes, aes(x, y, colour = lev), size = c(19, 19, 13, 13, 13, 13, rep(11, 12))) +
      geom_text(data = nodes, aes(x, y, label = lab), size = c(3.6, 3.6, rep(3.6, 4), rep(2.7, 12)),
                colour = "white", fontface = "bold") +
      annotate("text", x = 4, y = 0.45, label = "Time: the same measures repeated across visits (longitudinal data)",
               size = 4.2, colour = "#2F7A28") +
      scale_colour_manual(values = c(Site = "#2B6CB0", Participant = "#D95F02", Measure = "#7570B3"),
                          name = NULL) +
      scale_y_continuous(limits = c(0.3, 3.6), breaks = 1:3,
                         labels = c("Level 1\nmeasure", "Level 2\nparticipant", "Level 3\nsite")) +
      theme_void(base_size = 13) +
      theme(axis.text.y = element_text(colour = "#1E293B", size = 13), legend.position = "none")
  })

  ## ---------------------------------------------------------- data loading
  raw_up <- reactive({
    req(input$file)
    read.csv(input$file$datapath, stringsAsFactors = FALSE, check.names = FALSE)
  })

  output$map_ui <- renderUI({
    req(input$file)
    df <- raw_up()
    if (is_abide(df)) return(note("ABIDE phenotypic file recognised; columns are mapped automatically."))
    cn <- names(df); cno <- c("(none)", cn)
    tagList(
      selectInput("m_site", "Site", cn),
      selectInput("m_id", "Participant identifier", cn),
      selectInput("m_group", "Group (optional)", cno),
      uiOutput("m_glevel_ui"),
      selectInput("m_age", "Age (optional)", cno),
      selectInput("m_sex", "Sex (optional)", cno),
      uiOutput("m_sexlevel_ui"),
      selectizeInput("m_out", "Outcomes", cn, multiple = TRUE),
      selectInput("m_time", "Time since baseline (optional)", cno),
      selectInput("m_instr", "Instrument (optional)", cno),
      selectizeInput("m_sitecov", "Site-level exposures (optional)", cn, multiple = TRUE)
    )
  })
  output$m_glevel_ui <- renderUI({
    req(input$m_group, input$m_group != "(none)")
    selectInput("m_glevel", "Value coded as group 1", sort(unique(as.character(raw_up()[[input$m_group]]))))
  })
  output$m_sexlevel_ui <- renderUI({
    req(input$m_sex, input$m_sex != "(none)")
    selectInput("m_sexf", "Value coded as female", sort(unique(as.character(raw_up()[[input$m_sex]]))))
  })

  observeEvent(input$load, {
    res <- tryCatch({
      if (input$src == "demo") {
        prep_abide(sim_abide_raw(input$seed), source = "demo")
      } else if (input$src == "abide_url") {
        withProgress(message = "Downloading ABIDE I phenotypic file", {
          ph <- read.csv(url(ABIDE_URL), stringsAsFactors = FALSE)
        })
        prep_abide(ph, source = "abide")
      } else {
        df <- raw_up()
        if (is_abide(df)) prep_abide(df, source = "abide") else {
          validate(need(length(input$m_out) > 0, "Choose at least one outcome."))
          prep_generic(df, list(site = input$m_site, id = input$m_id, group = input$m_group,
                                group_level = input$m_glevel, age = input$m_age, sex = input$m_sex,
                                sex_female = input$m_sexf, outcomes = input$m_out,
                                time = input$m_time, instr = input$m_instr,
                                sitecov = input$m_sitecov))
        }
      }
    }, error = function(e) e)
    if (inherits(res, "error")) {
      showNotification(paste("Data not loaded:", conditionMessage(res)), type = "error", duration = 10)
    } else {
      reset_results(); D(res)
      showNotification(sprintf("Loaded %d participants from %d sites.", res$n_out,
                               nlevels(droplevels(res$d_base$site))), type = "message")
    }
  })

  output$data_status <- renderUI({
    x <- D()
    if (is.null(x)) return(card(card_body(tags$p("Choose a data source and press ", tags$b("Load data"), "."))))
    d <- x$d_base
    src <- switch(x$source, demo = "Simulated ABIDE-format demo (values are simulated, not ABIDE data)",
                  abide = "ABIDE I phenotypic file", upload = "Uploaded CSV")
    grp <- if ("dx" %in% names(d) && !is.null(x$groups))
      sprintf("; %s %d, %s %d", x$groups[1], sum(d$dx == 1, na.rm = TRUE), x$groups[2], sum(d$dx == 0, na.rm = TRUE)) else ""
    card(card_body(
      tags$p(tags$b(src)),
      tags$p(sprintf("%d rows read, %d participants analysed in %d sites%s.", x$n_in, x$n_out,
                     nlevels(droplevels(d$site)), grp),
             if (x$longitudinal) sprintf(" Longitudinal: %d participant-visits; the site and measure layers use each participant's first visit.", nrow(x$d)) else "")
    ))
  })

  output$site_tab <- renderDT({
    x <- req(D()); d <- x$d_base
    ss <- split(d, droplevels(d$site))
    tb <- do.call(rbind, lapply(names(ss), function(s) {
      z <- ss[[s]]
      data.frame(Site = s, N = nrow(z),
                 G1 = if ("dx" %in% names(z)) sum(z$dx == 1, na.rm = TRUE) else NA,
                 G0 = if ("dx" %in% names(z)) sum(z$dx == 0, na.rm = TRUE) else NA,
                 Female = if ("female" %in% names(z)) sum(z$female == 1, na.rm = TRUE) else NA,
                 `Median age` = if ("age" %in% names(z)) round(median(z$age, na.rm = TRUE), 1) else NA,
                 `Age range` = if ("age" %in% names(z)) paste0(f2(min(z$age, na.rm = TRUE)), "–", f2(max(z$age, na.rm = TRUE))) else NA,
                 Eyes = if ("eyes_closed" %in% names(z)) {
                   o <- sum(z$eyes_closed == 0, na.rm = TRUE); cl <- sum(z$eyes_closed == 1, na.rm = TRUE)
                   if (o > 0 && cl > 0) sprintf("Open (%d) / closed (%d)", o, cl) else if (cl > 0) "Closed" else "Open"
                 } else NA, check.names = FALSE, stringsAsFactors = FALSE)
    }))
    tb <- tb[order(-tb$N), ]
    if (!is.null(x$groups)) names(tb)[3:4] <- x$groups else tb <- tb[, -(3:4)]
    tb <- tb[, colSums(!is.na(tb)) > 0, drop = FALSE]
    datatable(tb, rownames = FALSE, options = list(pageLength = 25, dom = "tp"))
  })

  output$avail_plot <- renderPlot({
    x <- req(D()); d <- x$d_base; meas <- x$meas
    sizes <- sort(table(droplevels(d$site)), decreasing = TRUE)
    av <- do.call(rbind, lapply(seq_len(nrow(meas)), function(k) {
      v <- meas$var[k]
      do.call(rbind, lapply(names(sizes), function(s) {
        z <- d[d$site == s, ]
        if (!meas$both[k] && "dx" %in% names(z)) z <- z[!is.na(z$dx) & z$dx == 1, ]
        p <- if (nrow(z)) mean(!is.na(z[[v]])) else 0
        data.frame(measure = meas$label[k], block = meas$block[k],
                   site = sprintf("%s (%d)", s, sizes[[s]]), p = p)
      }))
    }))
    av$status <- cut(av$p, c(-1, 0, 0.8999, 1.1), labels = c("Not recorded", "Recorded for some", "Recorded for 90% or more"))
    av$site <- factor(av$site, levels = rev(unique(av$site)))
    av$measure <- factor(av$measure, levels = meas$label)
    ggplot(av, aes(measure, site, shape = status, colour = block)) +
      geom_point(size = 3.4, stroke = 1.1) +
      scale_shape_manual(values = c("Not recorded" = 4, "Recorded for some" = 1, "Recorded for 90% or more" = 16), name = NULL, drop = FALSE) +
      scale_colour_manual(values = BLOCK_COL, name = NULL) +
      labs(x = NULL, y = NULL) + theme_hd() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.box = "vertical") +
      guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 3))
  })

  ## ---------------------------------------------------------- site layer
  meas_choices <- reactive({ x <- req(D()); setNames(x$meas$var, x$meas$label) })

  output$site_meas_ui <- renderUI(
    checkboxGroupInput("site_meas", "Measures", meas_choices(), selected = meas_choices()))

  observeEvent(input$run_site, {
    x <- req(D()); ks <- which(x$meas$var %in% input$site_meas)
    validate(need(length(ks) > 0, "Choose measures."))
    out <- list(); fails <- character(0)
    withProgress(message = "Fitting site-layer models", value = 0, {
      for (k in ks) {
        incProgress(1 / length(ks), detail = x$meas$label[k])
        r <- tryCatch(site_layer_one(x, k, input$ci, input$nsim), error = function(e) e)
        if (inherits(r, "error")) fails <- c(fails, paste0(x$meas$label[k], " (", conditionMessage(r), ")")) else out[[length(out) + 1]] <- r
      }
    })
    if (length(fails)) showNotification(paste("Not estimated:", paste(fails, collapse = "; ")), type = "warning", duration = 12)
    if (length(out)) {
      R$site <- do.call(rbind, lapply(out, `[[`, "res"))
      R$blup <- do.call(rbind, lapply(out, `[[`, "blup"))
    }
  })

  output$icc_tab <- renderDT({
    r <- R$site
    validate(need(!is.null(r), "Press 'Run site layer' to estimate the site ICCs."))
    tb <- data.frame(Measure = r$label, Block = r$block, N = r$N, Sites = r$sites,
                     `Mean n per site` = f2(r$nbar), ICC = f3(r$icc), `95% interval` = paste0(f3(r$lo), " to ", f3(r$hi)),
                     `Adjusted ICC` = f3(r$icc_adj), DEFF = f2(r$deff), N_eff = round(r$neff),
                     `p (site variance)` = fp(r$p_site),
                     Singular = ifelse(r$singular, "yes", ""), check.names = FALSE)
    datatable(tb, rownames = FALSE, options = list(dom = "t", pageLength = 50))
  })

  output$icc_plot <- renderPlot({
    r <- R$site; validate(need(!is.null(r), "Run the site layer first."))
    p <- data.frame(label = r$label, block = r$block, icc = r$icc, lo = r$lo, hi = r$hi, set = "This analysis")
    if (isTRUE(input$show_leap)) {
      L <- subset(TABLE2, dataset == "LEAP")
      ci <- t(mapply(function(icc, N, k) {
        nbar <- N / k; F <- 1 + icc * nbar / (1 - icc)
        FL <- F / qf(.975, k - 1, N - k); FU <- F / qf(.025, k - 1, N - k)
        c(max(0, (FL - 1) / (FL - 1 + nbar)), max(0, (FU - 1) / (FU - 1 + nbar)))
      }, L$icc, L$N, L$sites))
      p <- rbind(p, data.frame(label = L$measure, block = "LEAP (Table 2)", icc = L$icc,
                               lo = ci[, 1], hi = ci[, 2], set = "LEAP"))
    }
    p$label <- factor(p$label, levels = unique(p$label[order(p$set != "This analysis", p$icc)]))
    p$set <- factor(p$set, levels = c("This analysis", "LEAP"))
    ggplot(p, aes(icc, label, colour = block)) +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.8) +
      geom_point(size = 3) +
      facet_grid(set ~ ., scales = "free_y", space = "free_y") +
      scale_colour_manual(values = BLOCK_COL, name = NULL) +
      labs(x = "Site intraclass correlation", y = NULL) + theme_hd() +
      guides(colour = guide_legend(ncol = 1))
  })

  output$deff_plot <- renderPlot({
    nb <- seq(1, 200, length.out = 300)
    cur <- expand.grid(nbar = nb, ICC = c(0.01, 0.05, 0.10, 0.20, 0.50))
    cur$deff <- 1 + (cur$nbar - 1) * cur$ICC
    g <- ggplot(cur, aes(nbar, deff)) +
      geom_line(aes(group = ICC, linetype = factor(ICC)), colour = "grey45") +
      scale_linetype_manual(values = c("dotted", "dotdash", "dashed", "longdash", "solid"), name = "ICC") +
      scale_y_log10() + labs(x = "Participants per site", y = "Design effect (log scale)") + theme_hd()
    if (!is.null(R$site)) {
      r <- R$site
      g <- g + geom_point(data = data.frame(nbar = r$nbar, deff = r$deff, block = r$block),
                          aes(colour = block), size = 3) +
        scale_colour_manual(values = BLOCK_COL, name = NULL)
    }
    L <- subset(TABLE2, dataset == "LEAP")
    g + geom_point(data = data.frame(nbar = L$nbar, deff = 1 + (L$nbar - 1) * L$icc),
                   shape = 17, size = 3, colour = BLOCK_COL[["LEAP (Table 2)"]]) +
      guides(colour = guide_legend(ncol = 2), linetype = guide_legend(nrow = 1)) +
      theme(legend.box = "vertical")
  })

  output$blup_ui <- renderUI({
    req(R$blup)
    ch <- unique(R$blup$label)
    pre <- intersect(c("Full-scale IQ", "ADOS total", "Functional DVARS"), ch)
    if (!length(pre)) pre <- head(ch, 3)
    selectizeInput("blup_sel", "Measures", ch, selected = pre, multiple = TRUE)
  })
  output$blup_plot <- renderPlot({
    b <- R$blup; validate(need(!is.null(b), "Run the site layer first."))
    b <- b[b$label %in% input$blup_sel, ]; req(nrow(b) > 0)
    ord_var <- if ("Functional DVARS" %in% b$label) "Functional DVARS" else b$label[1]
    o <- b[b$label == ord_var, ]; lev <- o$site[order(o$est)]
    b$site <- factor(b$site, levels = unique(c(lev, setdiff(b$site, lev))))
    ggplot(b, aes(est, site, colour = label)) +
      geom_vline(xintercept = 0, colour = "grey60") +
      geom_errorbarh(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), height = 0,
                     position = position_dodge(width = 0.6), linewidth = 0.7) +
      geom_point(position = position_dodge(width = 0.6), size = 2.4) +
      scale_colour_brewer(palette = "Dark2", name = NULL) +
      labs(x = "Site effect (SD units)", y = NULL,
           subtitle = paste("Sites ordered by the", ord_var, "effect")) + theme_hd()
  })

  ## ---------------------------------------------------------- exposure level
  output$expo_ui <- renderUI({
    x <- req(D())
    if (!length(x$expo)) return(note("No site-level exposure is available. For uploaded data, map one under 'Site-level exposures'."))
    tagList(selectInput("ex_y", "Outcome", meas_choices(),
                        selected = if ("dvars" %in% x$meas$var) "dvars" else x$meas$var[1]),
            selectInput("ex_x", "Exposure", setNames(x$expo, x$expo_labels[x$expo])))
  })
  observeEvent(input$run_expo, {
    x <- req(D()); req(input$ex_x)
    R$expo <- tryCatch(expo_test(x, input$ex_y, input$ex_x),
                       error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
  })
  output$expo_tab <- renderTable({
    e <- R$expo; validate(need(!is.null(e), "Choose an outcome and exposure and press 'Test exposure'."))
    data.frame(Analysis = e$analysis, Estimate = f3(e$est), SE = f3(e$se), df = f2(e$df), p = fp(e$p))
  })
  output$expo_text <- renderUI({
    e <- R$expo; req(e)
    tags$p(sprintf("N = %d in %d sites; \u03c1x of the exposure = %s. Ignoring site gave %s SD (%s); with a random site intercept the estimate was %s SD (%s, %s degrees of freedom).",
                   e$N[1], e$sites[1], f3(e$rho_x[1]), f2(e$est[1]), fpt(e$p[1]), f2(e$est[2]), fpt(e$p[2]), f2(e$df[2])))
  })

  rho_df <- reactive({
    x <- req(D())
    covs <- intersect(x$covs, names(x$d_base))
    lab <- c(dx = "Group", age = "Age", female = "Female", eyes_closed = "Eyes closed at scan")
    data.frame(var = covs, label = ifelse(covs %in% names(lab), lab[covs], covs),
               rho = sapply(covs, function(v) tryCatch(rho_x(x$d_base, v), error = function(e) NA)))
  })
  output$rho_plot <- renderPlot({
    r <- rho_df(); r <- r[!is.na(r$rho), ]; req(nrow(r) > 0)
    r$cls <- cut(r$rho, c(-1, 0.1, 0.9, 2), labels = c("Participant property", "Within and between sites", "Site property"))
    r$label <- factor(r$label, levels = r$label[order(r$rho)])
    ggplot(r, aes(rho, label, colour = cls)) +
      annotate("rect", xmin = 0.9, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#FDE0DD", alpha = 0.6) +
      annotate("rect", xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = "#E5F5E0", alpha = 0.6) +
      geom_segment(aes(x = 0, xend = rho, yend = label), linewidth = 1) +
      geom_point(size = 4) +
      scale_colour_manual(values = c("Participant property" = "#1B9E77", "Within and between sites" = "#7570B3", "Site property" = "#D95F02"), name = NULL, drop = FALSE) +
      scale_x_continuous(limits = c(0, 1)) + labs(x = expression(rho[x]), y = NULL) + theme_hd() +
      guides(colour = guide_legend(ncol = 2))
  })

  output$infl_plot <- renderPlot({
    r <- R$site
    nb <- if (!is.null(r)) round(mean(r$nbar)) else 44
    ic <- if (!is.null(r)) sort(unique(round(c(min(r$icc), median(r$icc), max(r$icc)), 3))) else c(0.05, 0.19, 0.57)
    g <- expand.grid(rho = seq(0, 1, 0.01), icc = ic)
    g$ratio <- 1 + (nb - 1) * g$rho * g$icc
    ggplot(g, aes(rho, ratio, colour = factor(icc))) + geom_line(linewidth = 1.1) +
      scale_colour_brewer(palette = "Set1", name = expression(ICC[Y])) +
      labs(x = expression(rho[x] ~ "of the exposure"), y = "Variance ratio (correct / naive)",
           subtitle = sprintf("Mean of %d participants per site", nb)) + theme_hd()
  })

  output$null_meas_ui <- renderUI({
    ch <- meas_choices()
    pre <- intersect(c("FIQ", "PIQ", "ADOS_total", "ADIR_social", "log_fd", "dvars"), ch)
    if (!length(pre)) pre <- head(ch, 3)
    selectizeInput("null_meas", "Outcomes", ch, selected = pre, multiple = TRUE)
  })
  observeEvent(input$run_null, {
    x <- req(D()); vs <- req(input$null_meas)
    set.seed(input$null_seed)
    tot <- length(vs) * input$B
    withProgress(message = "Null experiment", value = 0, {
      out <- lapply(vs, function(v) {
        setProgress(detail = x$meas$label[x$meas$var == v])
        tryCatch(null_exp(x, v, input$B, function() incProgress(1 / tot)), error = function(e) NULL)
      })
    })
    o <- do.call(rbind, out)
    if (!is.null(o)) o$label <- x$meas$label[match(o$var, x$meas$var)]
    R$null <- o
  })
  output$null_plot <- renderPlot({
    o <- R$null; validate(need(!is.null(o), "Choose outcomes and press 'Run null experiment'."))
    o$label <- factor(o$label, levels = unique(o$label[order(o$rate * (o$analysis == "Site ignored") * (o$assignment == "Assigned to whole sites"))]))
    w <- reshape(o[, c("label", "assignment", "analysis", "rate")], idvar = c("label", "assignment"),
                 timevar = "analysis", direction = "wide")
    ggplot(o, aes(rate, label)) +
      geom_vline(xintercept = 0.05, linetype = "dashed", colour = "grey40") +
      geom_segment(data = w, aes(x = `rate.Random site intercept`, xend = `rate.Site ignored`, y = label, yend = label),
                   colour = "grey70", linewidth = 1.2) +
      geom_point(aes(colour = analysis), size = 4) +
      facet_wrap(~ assignment) +
      scale_colour_manual(values = c("Site ignored" = "#E41A1C", "Random site intercept" = "#377EB8"), name = NULL) +
      scale_x_continuous(labels = function(z) paste0(round(100 * z), "%")) +
      labs(x = "False-positive rate (nominal 5%)", y = NULL) + theme_hd()
  })
  output$null_tab <- renderDT({
    o <- req(R$null)
    datatable(data.frame(Outcome = o$label, Assignment = o$assignment, Analysis = o$analysis,
                         `False-positive rate` = paste0(round(100 * o$rate, 1), "%"), Replicates = o$B,
                         Sites = o$sites, check.names = FALSE),
              rownames = FALSE, options = list(dom = "t", pageLength = 50))
  })

  ## ---------------------------------------------------------- group difference
  dx_choices <- reactive({
    x <- req(D())
    if (!"dx" %in% names(x$d_base)) return(NULL)
    m <- x$meas[x$meas$both, ]; setNames(m$var, m$label)
  })
  output$dx_meas_ui <- renderUI({
    ch <- dx_choices()
    if (is.null(ch) || !length(ch)) return(note("A group variable recorded in both groups is needed. For uploaded data, map 'Group'."))
    selectInput("dx_v", "Measure", ch, selected = if ("FIQ" %in% ch) "FIQ" else ch[1])
  })
  observeEvent(input$run_dx, {
    x <- req(D()); v <- req(input$dx_v)
    withProgress(message = "Fitting the four analyses", {
      R$dx <- tryCatch(c(dx_four(x, v), list(var = v, label = x$meas$label[x$meas$var == v])),
                       error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
      R$ps <- tryCatch(persite_pool(x, v), error = function(e) NULL)
    })
  })
  output$dx_plot <- renderPlot({
    r <- R$dx; validate(need(!is.null(r), "Choose a measure and press 'Fit the four analyses'."))
    t <- r$tab; t$model <- factor(t$model, levels = rev(names(MODEL_COL)))
    ggplot(t, aes(est, model, colour = model)) +
      geom_vline(xintercept = 0, colour = "grey60") +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 1) +
      geom_point(size = 4) +
      scale_colour_manual(values = MODEL_COL, guide = "none") +
      labs(x = paste0(r$label, ": ", D()$groups[1], " minus ", D()$groups[2], " (SD units)"), y = NULL) + theme_hd()
  })
  output$dx_tab <- renderTable({
    r <- req(R$dx); t <- r$tab
    data.frame(Analysis = t$model, Estimate = f3(t$est), SE = f3(t$se), df = f2(t$df),
               `95% CI` = paste0(f2(t$lo), " to ", f2(t$hi)), check.names = FALSE)
  })
  output$dx_text <- renderUI({
    r <- req(R$dx); t <- r$tab
    tags$p(sprintf("N = %d in %d sites. Between-site SD of the group difference \u03c4 = %s (likelihood-ratio %s). Standard error %s with a random intercept and %s with the random slope. Prediction interval for a new site (16): %s to %s SD.%s",
                   r$N, r$S, f2(r$tau), fpt(r$p_slope), f3(t$se[3]), f3(t$se[4]), f2(r$pi[1]), f2(r$pi[2]),
                   if (r$singular) " The random-slope fit is singular; the slope variance is at the boundary." else ""))
  })
  output$ps_plot <- renderPlot({
    p <- R$ps; validate(need(!is.null(p), "Site-specific differences need at least three sites with three or more participants in each group."))
    ps <- p$ps; pl <- p$pool
    ps$site <- factor(ps$site, levels = ps$site[order(ps$est)])
    ggplot(ps, aes(est, site)) +
      annotate("rect", xmin = pl$pi_lo, xmax = pl$pi_hi, ymin = -Inf, ymax = Inf, fill = "#FBB4C9", alpha = 0.35) +
      annotate("rect", xmin = pl$re - 1.96 * pl$re_se, xmax = pl$re + 1.96 * pl$re_se, ymin = -Inf, ymax = Inf, fill = "#9ECAE1", alpha = 0.6) +
      geom_vline(xintercept = pl$re, colour = "#2171B5", linewidth = 1) +
      geom_vline(xintercept = 0, colour = "grey50", linetype = "dotted") +
      geom_errorbarh(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), height = 0, colour = "grey30") +
      geom_point(aes(size = n), colour = "#08519C") +
      scale_size_area(max_size = 7, name = "Site n") +
      labs(x = "Site-specific difference (SD units)", y = NULL) + theme_hd()
  })
  output$ps_text <- renderUI({
    p <- req(R$ps); pl <- p$pool
    tags$p(sprintf("%d sites; pooled difference %s SD (SE %s); \u03c4 = %s; I\u00b2 = %d%%; Q %s; 95%% prediction interval %s to %s SD. Sites with fewer than three participants in either group are omitted.",
                   pl$k, f2(pl$re), f3(pl$re_se), f2(pl$tau), round(100 * pl$I2), fpt(pl$pQ), f2(pl$pi_lo), f2(pl$pi_hi)))
  })

  ## ---------------------------------------------------------- site composition
  output$cx_meas_ui <- renderUI({
    x <- req(D())
    tagList(
      if (!is.null(x$instrument)) selectInput("cx_v", "Measure for the instrument model", meas_choices(),
                                              selected = if ("FIQ" %in% x$meas$var) "FIQ" else x$meas$var[1])
      else note("No instrument variable is available, so the crossed instrument model is not fitted."),
      selectizeInput("inst_v", "Measures for the institution/cohort model", meas_choices(),
                     selected = intersect(c("FIQ", "log_fd", "dvars"), x$meas$var), multiple = TRUE)
    )
  })
  observeEvent(input$run_cx, {
    x <- req(D())
    R$cx <- if (!is.null(x$instrument) && !is.null(input$cx_v))
      tryCatch(crossed_instr(x, input$cx_v), error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    multi <- any(table(unique(x$d_base[, c("institution", "site")])$institution) > 1)
    R$inst <- if (multi && length(input$inst_v))
      tryCatch(inst_cohort(x, input$inst_v), error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    else "none"
  })
  output$cx_tab <- renderTable({
    c <- R$cx; validate(need(!is.null(c), "Press 'Decompose site' (needs an instrument variable)."))
    data.frame(Model = c$model, `Site ICC` = f3(c$site_icc), `Instrument ICC` = ifelse(is.na(c$instrument_icc), "", f3(c$instrument_icc)),
               `p (instrument)` = fp(c$p_instrument), check.names = FALSE)
  })
  output$cx_text <- renderUI({
    c <- req(R$cx)
    tags$p(sprintf("%d participants with a recorded instrument, %d instruments across %d sites; %d sites used more than one instrument and %d instruments were used at more than one site, so instrument is crossed with site.",
                   c$N[1], c$instruments[1], c$sites[1], c$sites_multi_instr[1], c$instr_multi_site[1]))
  })
  output$inst_plot <- renderPlot({
    i <- R$inst
    validate(need(!is.null(i), "Press 'Decompose site'."),
             need(!identical(i, "none"), "No institution contributed more than one cohort, so the institution and cohort levels coincide and the decomposition was not fitted."))
    x <- D(); i$label <- x$meas$label[match(i$var, x$meas$var)]
    l <- rbind(data.frame(label = i$label, level = "Institution", share = i$institution),
               data.frame(label = i$label, level = "Cohort within institution", share = i$cohort_within))
    w <- data.frame(label = i$label, a = i$institution, b = i$cohort_within)
    ggplot(l, aes(share, label)) +
      geom_segment(data = w, aes(x = b, xend = a, y = label, yend = label), colour = "grey70", linewidth = 1.2) +
      geom_point(aes(colour = level), size = 4) +
      scale_colour_manual(values = c("Institution" = "#2B6CB0", "Cohort within institution" = "#E7298A"), name = NULL) +
      labs(x = "Share of total variance", y = NULL) + theme_hd()
  })
  output$inst_text <- renderUI({
    i <- R$inst; req(!is.null(i), !identical(i, "none"))
    tags$p("Where most of the site variance sits at the institution level, cohorts from the same institution share their scanner or examiner effects, and the institution is the better site unit.")
  })

  ## ---------------------------------------------------------- measure layer
  output$mv_meas_ui <- renderUI({
    ch <- meas_choices()
    pre <- intersect(c("FIQ", "log_fd", "dvars"), ch); if (length(pre) < 2) pre <- head(ch, 3)
    selectizeInput("mv_v", "Measures (2 to 5)", ch, selected = pre, multiple = TRUE, options = list(maxItems = 5))
  })
  output$mv_sub_ui <- renderUI({
    x <- req(D())
    if (!"dx" %in% names(x$d_base) || is.null(x$groups)) return(NULL)
    radioButtons("mv_sub", "Participants", setNames(c("all", "g1"), c("All participants", paste(x$groups[1], "participants only"))))
  })
  observeEvent(input$run_mv, {
    x <- req(D()); v <- req(input$mv_v)
    validate(need(length(v) >= 2, "Choose at least two measures."))
    withProgress(message = "Estimating correlations", {
      R$mv <- tryCatch(mv_fit(x, v, identical(input$mv_sub, "g1"), input$mv_model),
                       error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    })
    if (!is.null(R$mv)) {
      R$mv$tab$l1 <- x$meas$label[match(R$mv$tab$m1, x$meas$var)]
      R$mv$tab$l2 <- x$meas$label[match(R$mv$tab$m2, x$meas$var)]
    }
  })
  output$mv_plot <- renderPlot({
    m <- R$mv; validate(need(!is.null(m), "Choose measures and press 'Estimate correlations'."))
    t <- m$tab; t$pair <- paste(t$l1, "\u2013", t$l2)
    l <- rbind(
      data.frame(pair = t$pair, type = "Pooled (site ignored)", r = t$pooled),
      data.frame(pair = t$pair, type = "Within sites", r = ifelse(is.na(t$within_model), t$within_centred, t$within_model)),
      data.frame(pair = t$pair, type = "Between sites", r = ifelse(is.na(t$between_model), t$between_sitemeans, t$between_model)))
    rng <- aggregate(r ~ pair, l, range)
    rng <- data.frame(pair = rng$pair, lo = rng$r[, 1], hi = rng$r[, 2])
    ggplot(l, aes(r, pair)) +
      geom_vline(xintercept = 0, colour = "grey60") +
      geom_segment(data = rng, aes(x = lo, xend = hi, y = pair, yend = pair), colour = "grey75", linewidth = 1.2) +
      geom_point(aes(colour = type), size = 4.2) +
      scale_colour_manual(values = c("Pooled (site ignored)" = "#999999", "Within sites" = "#1B9E77", "Between sites" = "#D95F02"), name = NULL) +
      scale_x_continuous(limits = c(-1, 1)) +
      labs(x = "Correlation", y = NULL, subtitle = sprintf("N = %d participants in %d sites", t$N[1], t$sites[1])) + theme_hd() +
      guides(colour = guide_legend(ncol = 2))
  })
  output$mv_msg <- renderUI({ m <- req(R$mv); if (!is.null(m$msg)) note(m$msg) })
  output$mv_tab <- renderDT({
    m <- req(R$mv); t <- m$tab
    datatable(data.frame(Pair = paste(t$l1, "\u2013", t$l2), Pooled = f3(t$pooled),
                         `Within (centred)` = f3(t$within_centred), `Within (model Σe)` = f3(t$within_model),
                         `Between (site means)` = f3(t$between_sitemeans), `Between (model Σu)` = f3(t$between_model),
                         `Pooled from (9)` = f3(t$pooled_eq9), check.names = FALSE),
              rownames = FALSE, options = list(dom = "t"))
  })

  ## ---------------------------------------------------------- time layer
  times <- reactive({
    t <- suppressWarnings(as.numeric(strsplit(input$times, "[,; ]+")[[1]]))
    t <- t[!is.na(t)]; validate(need(length(t) >= 2, "Give at least two visit times.")); t
  })
  vc <- reactive(visit_corr(times(), input$sb0, input$sb1, input$sb01, input$se2, input$phi, input$su2))
  output$vc_tab <- renderTable({
    v <- vc(); tt <- times(); k <- length(tt)
    pr <- t(combn(k, 2))
    data.frame(`Visit pair` = paste0(pr[, 1], "–", pr[, 2]), `Gap (t)` = f2(tt[pr[, 2]] - tt[pr[, 1]]),
               Covariance = f2(v$cov[pr]), Correlation = f3(v$cor[pr]), check.names = FALSE)
  })
  output$vc_plot <- renderPlot({
    tt <- times(); k <- length(tt); pr <- t(combn(k, 2))
    comp <- function(i, j) c(`Stable between-person difference` = input$sb0,
                             `Between-person difference in change` = tt[i] * tt[j] * input$sb1 + (tt[i] + tt[j]) * input$sb01,
                             `Occasion carry-over` = input$se2 * input$phi^abs(i - j),
                             `Site` = input$su2)
    l <- do.call(rbind, lapply(seq_len(nrow(pr)), function(r) {
      cc <- comp(pr[r, 1], pr[r, 2])
      data.frame(pair = paste0("Visits ", pr[r, 1], "–", pr[r, 2]), part = names(cc), value = cc)
    }))
    l$part <- factor(l$part, levels = names(comp(1, 2)))
    ggplot(l, aes(value, pair, colour = part)) +
      geom_vline(xintercept = 0, colour = "grey60") +
      geom_point(size = 4, position = position_dodge(width = 0.55)) +
      scale_colour_brewer(palette = "Set2", name = NULL) +
      labs(x = "Contribution to the covariance between visits", y = NULL) + theme_hd() +
      guides(colour = guide_legend(nrow = 2))
  })

  observeEvent(input$run_ls, {
    tt <- times()
    withProgress(message = "Simulating and fitting", {
      R$ls <- tryCatch({
        dd <- sim_long(input$ls_S, input$ls_n, tt, input$sb0, input$sb1, input$sb01, input$se2,
                       input$phi, max(input$su2, 1e-6), input$ls_bt)
        fit_long(dd, "y", "time", irregular = FALSE)
      }, error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    })
  })
  output$ls_tab <- renderTable({
    f <- R$ls; validate(need(!is.null(f), "Press 'Simulate and fit' to generate data from the components on the left and re-estimate them."))
    nm <- c("σ²u (site)", "σ²b0", "σ²b1", "σb01", "σ²ε", "φ", "Average change", "AIC", "BIC")
    truth <- c(input$su2, input$sb0, input$sb1, input$sb01, input$se2, input$phi, input$ls_bt, NA, NA)
    data.frame(Component = nm, True = ifelse(is.na(truth), "", f2(truth)),
               `Independent residuals` = f2(f$ind), `AR(1) residuals` = if (!is.null(f$ar)) f2(f$ar) else "not fitted",
               check.names = FALSE)
  })
  output$ls_text <- renderUI({
    f <- req(R$ls)
    best <- if (!is.null(f$ar) && f$ar["AIC"] < f$ind["AIC"]) "AR(1)" else "independent"
    tags$p(sprintf("%d participants in %d sites, %d observations. AIC prefers the %s residual structure.", f$N, f$S, f$obs, best))
  })

  output$long_up_ui <- renderUI({
    x <- D()
    if (is.null(x) || !x$longitudinal) return(NULL)
    tagList(hr(), tags$h6("Uploaded longitudinal data"),
            selectInput("lu_v", "Outcome", meas_choices()),
            checkboxInput("lu_irr", "Irregular visit times (continuous-time AR(1))", TRUE),
            actionButton("run_lu", "Fit longitudinal model", class = "btn-primary w-100"))
  })
  observeEvent(input$run_lu, {
    x <- req(D())
    withProgress(message = "Fitting longitudinal model", {
      R$lu <- tryCatch(fit_long(x$d, input$lu_v, "time", input$lu_irr),
                       error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL })
    })
  })
  output$lu_text <- renderUI({
    x <- D()
    if (is.null(x)) return(note("Load data first."))
    if (!x$longitudinal)
      return(tags$p(if (x$source %in% c("demo", "abide"))
        "Not fitted: ABIDE I is a single-visit dataset, so participant random slopes and the correlation between visits cannot be estimated from it. The calculator and the simulation above show the time layer; upload longitudinal data with a time variable to fit it."
        else "Not fitted: no time variable was mapped. Map 'Time since baseline' on the Data tab to fit the longitudinal model."))
    f <- R$lu
    if (is.null(f)) return(note("Choose an outcome and press 'Fit longitudinal model'."))
    tags$p(sprintf("%d participants in %d sites, %d observations. Site random intercept, participant random intercept and slope, %s residuals.",
                   f$N, f$S, f$obs, if (isTRUE(input$lu_irr)) "continuous-time AR(1)" else "AR(1)"))
  })
  output$lu_tab <- renderTable({
    f <- req(R$lu)
    nm <- c("σ²u (site)", "σ²b0", "σ²b1", "σb01", "σ²ε", "φ", "Average change", "AIC", "BIC")
    data.frame(Component = nm, `Independent residuals` = f2(f$ind),
               `Autoregressive residuals` = if (!is.null(f$ar)) f2(f$ar) else "not fitted", check.names = FALSE)
  })

  ## ---------------------------------------------------------- planning
  observeEvent(input$preset, {
    req(input$preset != "")
    r <- TABLE2[as.integer(input$preset), ]
    updateNumericInput(session, "pl_icc", value = r$icc)
    updateNumericInput(session, "pl_n", value = round(r$nbar))
    updateNumericInput(session, "pl_S", value = r$sites)
  })
  pl <- reactive({
    deff <- 1 + (input$pl_n - 1) * input$pl_rho * input$pl_icc
    N <- input$pl_n * input$pl_S
    list(deff = deff, N = N, neff = N / deff)
  })
  output$vb_deff <- renderText(f2(pl()$deff))
  output$vb_N <- renderText(format(pl()$N, big.mark = ","))
  output$vb_neff <- renderText(format(round(pl()$neff), big.mark = ","))
  output$pl_neff <- renderPlot({
    nb <- 2:250
    g <- expand.grid(n = nb, rho = sort(unique(c(1, 0.5, input$pl_rho))))
    g$neff <- g$n * input$pl_S / (1 + (g$n - 1) * g$rho * input$pl_icc)
    ggplot(g, aes(n, neff, colour = factor(rho))) + geom_line(linewidth = 1.1) +
      geom_point(data = data.frame(n = input$pl_n, neff = pl()$neff), aes(n, neff), inherit.aes = FALSE,
                 size = 4, colour = "black") +
      geom_hline(yintercept = input$pl_S / input$pl_icc * (input$pl_icc > 0), linetype = "dotted", colour = "grey50") +
      scale_colour_brewer(palette = "Dark2", name = expression(rho[x])) +
      labs(x = "Participants per site", y = "Effective sample size",
           subtitle = sprintf("%d sites, ICC %s; dotted line: limit S/ICC for a site-level exposure", input$pl_S, f3(input$pl_icc))) + theme_hd()
  })
  output$pl_se <- renderPlot({
    S0 <- input$pl_S; n0 <- input$pl_n; tau <- input$pl_tau; s2 <- input$pl_se2
    se <- function(S, n) sqrt((tau^2 + 4 * s2 / n) / S)
    mult <- seq(1, 4, by = 0.05)
    l <- rbind(data.frame(N = S0 * mult * n0, se = se(S0 * mult, n0), way = "Add sites"),
               data.frame(N = S0 * n0 * mult, se = se(S0, n0 * mult), way = "Add participants per site"))
    ggplot(l, aes(N, se, colour = way)) + geom_line(linewidth = 1.2) +
      geom_point(data = data.frame(N = S0 * n0, se = se(S0, n0)), aes(N, se), inherit.aes = FALSE, size = 4) +
      scale_colour_manual(values = c("Add sites" = "#2B6CB0", "Add participants per site" = "#E7298A"), name = NULL) +
      labs(x = "Total participants", y = "Standard error of the group difference (SD units)") + theme_hd()
  })
  output$pl_text <- renderUI({
    se0 <- sqrt((input$pl_tau^2 + 4 * input$pl_se2 / input$pl_n) / input$pl_S)
    tags$p(sprintf("With %d sites of %d participants, \u03c4 = %s and \u03c3\u00b2\u03b5 = %s, (17) gives a standard error of %s. The term \u03c4\u00b2/S does not fall as participants are added within sites.",
                   input$pl_S, input$pl_n, f2(input$pl_tau), f2(input$pl_se2), f3(se0)))
  })

  ## ---------------------------------------------------------- few sites
  eta <- eventReactive(input$run_eta, ignoreNULL = FALSE, {
    lines <- strsplit(input$eta_in, "\n")[[1]]; lines <- lines[nzchar(trimws(lines))]
    do.call(rbind, lapply(lines, function(l) {
      p <- trimws(strsplit(l, ",")[[1]]); n <- length(p)
      if (n < 4) return(NULL)
      k <- as.numeric(p[n - 2]); N <- as.numeric(p[n - 1]); e2 <- as.numeric(p[n])
      nm <- paste(p[1:(n - 3)], collapse = ", ")
      nbar <- N / k
      F <- (e2 / (k - 1)) / ((1 - e2) / (N - k))
      FL <- F / qf(.975, k - 1, N - k); FU <- F / qf(.025, k - 1, N - k)
      icc <- function(f) max(0, (f - 1) / (f - 1 + nbar))
      data.frame(measure = nm, k = k, N = N, nbar = nbar, eta2 = e2, chance = (k - 1) / (N - 1),
                 F = F, icc = icc(F), lo = icc(FL), hi = icc(FU), deff = 1 + (nbar - 1) * icc(F),
                 p = pf(F, k - 1, N - k, lower.tail = FALSE))
    }))
  })
  output$eta_tab <- renderTable({
    e <- eta(); validate(need(!is.null(e) && nrow(e) > 0, "Enter rows as: measure, k, N, eta-squared."))
    data.frame(Measure = e$measure, k = e$k, N = e$N, `η²` = formatC(e$eta2, format = "f", digits = 4),
               `Chance η²` = formatC(e$chance, format = "f", digits = 4), ICC = f3(e$icc),
               `95% interval` = paste0(f3(e$lo), " to ", f3(e$hi)), DEFF = f2(e$deff), p = fp(e$p), check.names = FALSE)
  })
  output$eta_plot <- renderPlot({
    e <- eta(); req(e)
    e$measure <- factor(e$measure, levels = e$measure[order(e$icc)])
    ggplot(e, aes(icc, measure)) +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, colour = BLOCK_COL[["LEAP (Table 2)"]], linewidth = 0.9) +
      geom_point(size = 3.5, colour = "#B8860B") +
      geom_point(aes(x = eta2), shape = 4, size = 3, colour = "#E41A1C") +
      labs(x = "ICC (point and interval); reported η² (cross)", y = NULL) + theme_hd()
  })

  ## ---------------------------------------------------------- report
  report <- reactive({
    x <- D(); out <- character(0)
    if (is.null(x)) return("Load data and run the analyses to build the text.")
    d <- x$d_base
    src <- switch(x$source, demo = "a simulated dataset in the format of the ABIDE I phenotypic file",
                  abide = "the ABIDE I phenotypic file (ABIDE Preprocessed release)", upload = "the uploaded dataset")
    out <- c(out, "METHODS", "",
             sprintf("Data were taken from %s: %d participants in %d sites. Site was treated as a random effect throughout.", src, x$n_out, nlevels(droplevels(d$site))),
             "For each measure the empty random-intercept model Y_is = mu + u_s + e_is, u_s ~ N(0, sigma2_u), e_is ~ N(0, sigma2_e) (1) was fitted by restricted maximum likelihood, and the site intraclass correlation computed as ICC = sigma2_u / (sigma2_u + sigma2_e) (2). The adjusted ICC came from model (3), with participant covariates split into within-site deviations and site means. The site variance was tested with the likelihood-ratio test of sigma2_u = 0 against a 50:50 mixture of chi-square distributions with 0 and 1 degrees of freedom. The design effect was DEFF = 1 + (n_bar - 1) ICC and the effective sample size N_eff = N / DEFF (4).",
             "Group differences were estimated under four analyses (single-level, site as a fixed effect, random site intercept, random site intercept and random slope for group) in SD units, with Kenward-Roger degrees of freedom for the mixed models. The random-slope model Y_is = b0 + (b1 + v_s) Dx_i + X_is' b + u_s + e_is, Var(v_s) = tau2 (6), gave the between-site SD of the group difference, and the prediction interval for a new site was b1_hat +/- t(0.975, S - 2) sqrt(SE(b1_hat)^2 + tau2_hat) (16).",
             "Site-level exposures were tested against between-site variation with a random site intercept and Kenward-Roger degrees of freedom. A real-data null experiment assigned a pseudo-exposure to a random half of the sites, and separately within each site, and recorded the rejection rate at the 5% level with and without the random site intercept.",
             "Correlations between measures were separated into within-site and between-site components with the multivariate multilevel model Y_mis = b0m + X_is' bm + u_ms + e_mis, u_s ~ N(0, Sigma_u), e_is ~ N(0, Sigma_e) (8), and by site-mean centring; the pooled correlation was related to them by r_pooled ~ sqrt(ICC1 ICC2) r_between + sqrt((1 - ICC1)(1 - ICC2)) r_within (9).",
             "", "RESULTS", "")
    r <- R$site
    if (!is.null(r)) {
      i1 <- which.min(r$icc); i2 <- which.max(r$icc)
      out <- c(out, sprintf("Site ICCs ranged from %s (%s) to %s (%s). At %s to %s participants per site, the design effects ran from %s to %s; for %s, %d observations carried the information of %d independent participants.",
                            f3(r$icc[i1]), r$label[i1], f3(r$icc[i2]), r$label[i2], round(min(r$nbar)), round(max(r$nbar)),
                            f2(min(r$deff)), f2(max(r$deff)), r$label[i2], r$N[i2], round(r$neff[i2])))
    } else out <- c(out, "Site layer: not run in this session, so no ICCs are reported.")
    dx <- R$dx
    if (!is.null(dx)) {
      t <- dx$tab
      out <- c(out, sprintf("For %s, the group difference was %s SD (SE %s) ignoring site, %s SD (SE %s) with a random site intercept and %s SD (SE %s) with a random slope; the between-site SD of the difference was %s (%s) and the 95%% prediction interval for a new site ran from %s to %s SD.",
                            dx$label, f2(t$est[1]), f3(t$se[1]), f2(t$est[3]), f3(t$se[3]), f2(t$est[4]), f3(t$se[4]),
                            f2(dx$tau), fpt(dx$p_slope), f2(dx$pi[1]), f2(dx$pi[2])))
    } else out <- c(out, "Group difference: not run in this session.")
    e <- R$expo
    if (!is.null(e)) out <- c(out, sprintf("The site-level exposure (rho_x = %s) gave %s SD (%s) in a single-level analysis and %s SD (%s) with a random site intercept.",
                                           f3(e$rho_x[1]), f2(e$est[1]), fpt(e$p[1]), f2(e$est[2]), fpt(e$p[2])))
    o <- R$null
    if (!is.null(o)) {
      s <- o[o$assignment == "Assigned to whole sites", ]
      out <- c(out, sprintf("In the null experiment (%d assignments per outcome), false-positive rates for a site-level pseudo-exposure were %d%% to %d%% when site was ignored and %d%% to %d%% with a random site intercept.",
                            o$B[1], round(100 * min(s$rate[s$analysis == "Site ignored"])), round(100 * max(s$rate[s$analysis == "Site ignored"])),
                            round(100 * min(s$rate[s$analysis == "Random site intercept"])), round(100 * max(s$rate[s$analysis == "Random site intercept"]))))
    }
    m <- R$mv
    if (!is.null(m)) {
      t <- m$tab
      for (j in seq_len(nrow(t))) {
        w <- ifelse(is.na(t$within_model[j]), t$within_centred[j], t$within_model[j])
        b <- ifelse(is.na(t$between_model[j]), t$between_sitemeans[j], t$between_model[j])
        out <- c(out, sprintf("%s and %s were correlated %s within sites and %s between sites (pooled %s).",
                              t$l1[j], t$l2[j], f2(w), f2(b), f2(t$pooled[j])))
      }
    }
    paste(out, collapse = "\n")
  })
  output$report_text <- renderUI(tags$pre(style = "white-space:pre-wrap;font-family:inherit;font-size:0.95rem;background:#F8FAFC;border:1px solid #E2E8F0;padding:12px", report()))
  output$dl_report <- downloadHandler("hierdep_methods_results.txt", function(f) writeLines(report(), f))
  dlcsv <- function(obj, name) downloadHandler(name, function(f) {
    o <- obj(); if (is.null(o)) o <- data.frame(note = "Not run in this session")
    write.csv(o, f, row.names = FALSE)
  })
  output$dl_icc <- dlcsv(function() R$site, "site_layer_icc.csv")
  output$dl_blup <- dlcsv(function() R$blup, "site_blups.csv")
  output$dl_dx <- dlcsv(function() if (!is.null(R$dx)) cbind(measure = R$dx$label, R$dx$tab, tau = R$dx$tau, p_slope = R$dx$p_slope, pi_lo = R$dx$pi[1], pi_hi = R$dx$pi[2]), "group_difference.csv")
  output$dl_null <- dlcsv(function() R$null, "null_experiment.csv")
  output$dl_mv <- dlcsv(function() if (!is.null(R$mv)) R$mv$tab, "measure_correlations.csv")

  ## ---------------------------------------------------------- study decision guide
  observe({
    x <- D()
    if (is.null(x)) return()
    ch <- meas_choices()
    ex <- c(if ("dx" %in% names(x$d_base) && !is.null(x$groups)) setNames("dx", paste0("Group (", x$groups[1], " vs ", x$groups[2], ")")),
            setNames(x$expo, x$expo_labels[x$expo]),
            if (has_var(x$d_base, "age")) c(Age = "age"), if (has_var(x$d_base, "female")) c(Female = "female"),
            ch)
    updateSelectizeInput(session, "dg_y", choices = ch, selected = if ("FIQ" %in% ch) "FIQ" else ch[1])
    updateSelectizeInput(session, "dg_x", choices = ex, selected = ex[1])
    updateSelectizeInput(session, "dg_joint", choices = ch, selected = character(0))
  })

  ## objective and hypothesis (both stages)
  hyp <- reactive(hyp_info(input$h_type, input$h_ylab, input$h_xlab, input$h_side))
  output$h_card <- renderUI({
    h <- hyp()
    side <- c(two = "two-sided", greater = "one-sided (greater than 0)", less = "one-sided (less than 0)")[[input$h_side]]
    withMathJax(tagList(
      tags$p(tags$b("Objective. "), input$h_obj),
      tags$p(tags$b("Primary outcome: "), input$h_ylab, ". ", tags$b("Exposure or comparison: "), input$h_xlab, "."),
      tags$p(tags$b("Parameter tested: "), h$sym, ", ", h$desc, "."),
      tags$p(tags$b("Null hypothesis H0: "), h$h0, ". ", tags$b("Alternative hypothesis H1: "), h$h1, "."),
      tags$p(sprintf("Tested at \u03b1 = %s, %s; effect of interest %s %s.", input$h_alpha, side, input$h_delta,
                     if (input$h_type == "corr") "(correlation)" else if (input$h_type %in% c("change", "groupchange")) "SD units per unit time" else "SD units")),
      helpText(h$eq)))
  })

  ## ---------------------------------------- stage 1: planning
  observeEvent(input$p_preset, {
    req(input$p_preset != "")
    updateNumericInput(session, "p_icc", value = TABLE2$icc[as.integer(input$p_preset)])
  })
  pargs <- reactive(list(S = input$p_S, n = input$p_n, icc = input$p_icc, tau = input$p_tau, miss = input$p_miss,
                         visits = input$p_visits, gap = input$p_gap, rel = input$p_rel, sb1 = input$p_sb1, taut = input$p_taut))
  ppow <- reactive(do.call(plan_power, c(list(ht = input$h_type, delta = input$h_delta, alpha = input$h_alpha, side = input$h_side), pargs())))
  preq <- reactive(do.call(plan_required, c(list(ht = input$h_type, delta = input$h_delta, alpha = input$h_alpha, side = input$h_side,
                                                 target = input$p_power), pargs())))
  output$p_vb_power <- renderText(if (is.na(ppow()$power)) "Not computed" else paste0(round(100 * ppow()$power), "%"))
  output$p_vb_S <- renderText(if (is.na(preq()$S)) "Not reached by 2,000" else sprintf("%d sites of %d", preq()$S, input$p_n))
  output$p_vb_n <- renderText(if (is.na(preq()$n)) sprintf("Not reachable with %d sites", input$p_S) else sprintf("%d per site with %d sites", preq()$n, input$p_S))
  output$p_vb_N <- renderText(format(input$p_S * input$p_n, big.mark = ","))

  output$p_pow_plot <- renderPlot({
    validate(need(!is.na(ppow()$power), ppow()$msg %||% "Power cannot be computed for these inputs."))
    Smax <- max(40, 2 * input$p_S, if (!is.na(preq()$S)) preq()$S + 5 else 0)
    ns <- unique(round(c(input$p_n / 2, input$p_n, 2 * input$p_n)))
    g <- expand.grid(S = 3:Smax, n = ns)
    g$power <- mapply(function(S, n) do.call(plan_power, c(list(ht = input$h_type, delta = input$h_delta, alpha = input$h_alpha,
                                                                  side = input$h_side), modifyList(pargs(), list(S = S, n = n))))$power,
                      g$S, g$n)
    ggplot(g, aes(S, power, colour = factor(n))) +
      geom_hline(yintercept = input$p_power, linetype = "dashed", colour = "grey40") +
      geom_line(linewidth = 1.1) +
      geom_point(data = data.frame(S = input$p_S, power = ppow()$power), aes(S, power), inherit.aes = FALSE, size = 4) +
      scale_colour_brewer(palette = "Dark2", name = "Participants per site") +
      scale_y_continuous(limits = c(0, 1), labels = function(z) paste0(round(100 * z), "%")) +
      labs(x = "Number of sites", y = "Power", subtitle = "Point: planned design; dashed line: target") + theme_hd()
  })
  output$p_pow_text <- renderUI({
    pp <- ppow(); ht <- input$h_type
    if (is.na(pp$power)) return(tags$p(pp$msg))
    n_e <- input$p_n * (1 - input$p_miss / 100)
    rho <- if (ht == "sitelevel") 1 else 0
    deff <- 1 + (input$p_n - 1) * rho * input$p_icc
    form <- switch(ht,
      group = "\\(\\mathrm{Var}(\\hat\\beta_1)=\\{\\tau^2+4(1-\\mathrm{ICC})/n\\}/S\\), equation (17) in SD units",
      sitelevel = "\\(\\mathrm{Var}(\\hat\\gamma)=4\\{\\mathrm{ICC}+(1-\\mathrm{ICC})/n\\}/S\\), with half of the sites exposed",
      corr = "\\(\\mathrm{SE}\\{\\tanh^{-1}(r_w)\\}=1/\\sqrt{N-S-3}\\), within-site correlation after removing site means",
      change = "\\(\\mathrm{Var}(\\hat\\beta_t)=\\{\\sigma^2_{b1}+(1-R)/SS_t\\}/(Sn)+\\tau_t^2/S\\), with reliability R and \\(SS_t=\\sum_j(t_j-\\bar t)^2\\)",
      groupchange = "\\(\\mathrm{Var}(\\hat\\beta_{x\\times t})=[\\tau_t^2+4\\{\\sigma^2_{b1}+(1-R)/SS_t\\}/n]/S\\)")
    withMathJax(tagList(
      tags$p("Standard error of the tested parameter: ", form, "."),
      tags$p(sprintf("With %d sites of %d participants (%s analysable after %s%% missing), the standard error is %s on %s degrees of freedom, giving power %s%% for an effect of %s at \u03b1 = %s (%s).",
                     input$p_S, input$p_n, f2(n_e), input$p_miss, f3(pp$se), if (is.finite(pp$df)) round(pp$df) else "large",
                     round(100 * pp$power), input$h_delta, input$h_alpha, if (input$h_side == "two") "two-sided" else "one-sided")),
      tags$p(sprintf("Design effect for the exposure %s; effective sample size %s of %s.", f2(deff),
                     round(input$p_S * input$p_n / deff), input$p_S * input$p_n)),
      if (ht %in% c("group", "groupchange") && ((ht == "group" && input$p_tau > 0) || (ht == "groupchange" && input$p_taut > 0)))
        tags$p("Because the effect varies between sites, the term \\(\\tau^2/S\\) sets a floor that adding participants within sites cannot remove; adding sites does."),
      if (!is.na(preq()$S)) tags$p(sprintf("The target power of %s%% is reached with %d sites of %d participants.", round(100 * input$p_power), preq()$S, input$p_n))
    ))
  })

  fplan <- reactive({
    ht <- input$h_type; h <- hyp()
    list(stage = "plan", htype = ht, h0 = h$h0, h1 = h$h1, y = "outcome",
         x = switch(ht, group = "group", sitelevel = "site_condition", corr = "measure2", change = "time", groupchange = "group"),
         S = input$p_S, nbar = input$p_n, icc = input$p_icc, icc_lo = NA, icc_hi = NA,
         rhox = if (ht == "sitelevel") 1 else 0,
         tau = if (ht %in% c("change", "groupchange")) input$p_taut else input$p_tau, p_tau = NA, est = NA, se = NA,
         instr = input$p_instr, inst_multi = input$p_inst, inst_level = "unchecked",
         n_meas = max(input$p_nmeas, if (ht == "corr") 2 else 1), predict = input$p_predict,
         long = input$p_visits > 1 || ht %in% c("change", "groupchange"), fu_sites = input$p_fu,
         irregular = FALSE, int_differ = input$p_intdiff, miss_sites = input$p_misssites,
         miss_pct = input$p_miss, singular = FALSE)
  })
  dg_plan <- reactive(dg_decide(fplan()))
  output$dg_plot_plan <- renderUI(dg_flow(dg_plan()$tab))
  plan_script <- reactive(dg_code(fplan(), dg_plan(), covs = c("age", "female")))
  output$p_code <- renderText(plan_script())
  output$p_dl_code <- downloadHandler("planned_analysis.R", function(f) writeLines(plan_script(), f))

  sap <- reactive({
    h <- hyp(); pp <- ppow(); t <- dg_plan()$tab; ht <- input$h_type
    side <- if (input$h_side == "two") "two-sided" else "one-sided"
    c("Objective",
      input$h_obj, "",
      "Hypotheses",
      sprintf("The parameter tested is %s, %s. H0: %s; H1: %s. The test is %s at alpha = %s.", h$sym, h$desc, h$h0, h$h1, side, input$h_alpha), "",
      "Design",
      sprintf("The study will recruit %d participants at each of %d sites (%d in total). Site is treated as a random effect in every analysis, and the site unit is defined before analysis and kept the same throughout.",
              input$p_n, input$p_S, input$p_S * input$p_n), "",
      "Analysis",
      paste(vapply(which(t$active), function(i) paste0(toupper(substr(t$action[i], 1, 1)), substring(t$action[i], 2), "."), ""), collapse = " "), "",
      "Sample size",
      if (is.na(pp$power)) pp$msg else
        sprintf("Assuming a site ICC of %s%s and %s%% missing values, %d sites of %d participants give a standard error of %s and power of %d%% to detect an effect of %s at %s alpha = %s.",
                f3(input$p_icc),
                if (ht == "group" && input$p_tau > 0) sprintf(", a between-site SD of the effect of %s", f2(input$p_tau))
                else if (ht %in% c("change", "groupchange")) sprintf(", test-retest reliability %s, %d visits %s time units apart", f2(input$p_rel), input$p_visits, input$p_gap) else "",
                input$p_miss, input$p_S, input$p_n, f3(pp$se), round(100 * pp$power), input$h_delta, side, input$h_alpha),
      "",
      "Reporting",
      paste0("For the primary outcome the report gives the number of sites, participants per site, site ICC with interval, design effect and effective sample size; the estimate of ",
             h$sym, " with its confidence interval, degrees of freedom and p value",
             if (ht == "group") ", the between-site SD of the effect and the prediction interval for a new site" else "", "."))
  })
  output$p_sap <- renderUI(tags$pre(style = "white-space:pre-wrap;font-family:inherit;font-size:0.93rem;background:#F8FAFC;border:1px solid #E2E8F0;padding:12px",
                                    paste(sap(), collapse = "\n")))
  output$p_dl_sap <- downloadHandler("statistical_analysis_plan.txt", function(f) writeLines(sap(), f))

  ## ---------------------------------------- stage 2: analysis
  x_for_type <- reactive(switch(input$h_type, change = NULL, input$dg_x))

  observeEvent(input$dg_fill, {
    x <- D()
    if (is.null(x)) { showNotification("Load data on the Data tab first, or enter the answers yourself.", type = "warning"); return() }
    req(input$dg_y)
    xx <- x_for_type()
    if (!is.null(xx) && !xx %in% names(x$d_base)) { showNotification("The exposure is not a column in the loaded data.", type = "error"); return() }
    f <- withProgress(message = "Working through the decisions", tryCatch(dg_facts(x, input$dg_y, xx, input$dg_joint),
                      error = function(e) { showNotification(conditionMessage(e), type = "error"); NULL }))
    req(f)
    un <- function(id, v) updateNumericInput(session, id, value = if (is.null(v) || is.na(v)) NA else round(v, 4))
    un("dg_S", f$S); un("dg_nbar", f$nbar); un("dg_icc", f$icc); un("dg_icc_lo", f$icc_lo); un("dg_icc_hi", f$icc_hi)
    un("dg_rhox", f$rhox); un("dg_tau", f$tau); un("dg_ptau", f$p_tau); un("dg_est", f$est); un("dg_se", f$se)
    un("dg_nmeas", max(f$n_meas, if (input$h_type == "corr") 2 else 1)); un("dg_fu", f$fu_sites)
    un("dg_misssites", f$miss_sites); un("dg_misspct", f$miss_pct)
    updateCheckboxInput(session, "dg_singular", value = f$singular)
    updateCheckboxInput(session, "dg_instr", value = f$instr)
    updateCheckboxInput(session, "dg_instmulti", value = f$inst_multi)
    updateSelectInput(session, "dg_instlevel", selected = f$inst_level)
    updateCheckboxInput(session, "dg_long", value = f$long)
    updateCheckboxInput(session, "dg_irr", value = f$irregular)
    updateCheckboxInput(session, "dg_intdiff", value = f$int_differ)
    showNotification("Answers filled from the data. Press 'Test the hypothesis' to fit the chosen model.", type = "message")
  })

  dg_f <- reactive({
    h <- hyp()
    list(stage = "analyse", htype = input$h_type, h0 = h$h0, h1 = h$h1,
         y = if (nzchar(input$dg_y %||% "")) input$dg_y else "y",
         x = if (input$h_type == "change") "time" else if (nzchar(input$dg_x %||% "")) input$dg_x else "x",
         S = input$dg_S, nbar = input$dg_nbar, icc = input$dg_icc, icc_lo = input$dg_icc_lo, icc_hi = input$dg_icc_hi,
         rhox = input$dg_rhox, tau = input$dg_tau, p_tau = input$dg_ptau, est = input$dg_est, se = input$dg_se,
         instr = input$dg_instr, inst_multi = input$dg_instmulti, inst_level = input$dg_instlevel,
         n_meas = input$dg_nmeas, predict = input$dg_predict,
         long = input$dg_long || input$h_type %in% c("change", "groupchange"), fu_sites = input$dg_fu,
         irregular = input$dg_irr, int_differ = input$dg_intdiff, miss_sites = input$dg_misssites,
         miss_pct = input$dg_misspct, singular = input$dg_singular)
  })
  dg <- reactive(dg_decide(dg_f()))

  observeEvent(input$dg_test, {
    x <- D()
    if (is.null(x)) { showNotification("Load data on the Data tab first.", type = "warning"); return() }
    R$hyp <- withProgress(message = "Testing the hypothesis",
      tryCatch(c(fit_hypothesis(x, input$h_type, input$dg_y, x_for_type(), dg(), dg_f(), input$h_alpha, input$h_side),
                 list(ht = input$h_type, alpha = input$h_alpha, side = input$h_side, delta = input$h_delta)),
               error = function(e) list(ok = FALSE, why = paste("The model was not fitted:", conditionMessage(e)))))
  })

  output$h_result <- renderUI({
    r <- R$hyp
    if (is.null(r)) return(note("Fill the answers from the data, check the decisions, then press 'Test the hypothesis with the chosen model'."))
    if (!isTRUE(r$ok)) return(tags$p(r$why))
    h <- hyp(); rej <- r$p < r$alpha
    dft <- if (is.finite(r$df)) sprintf("%s df", f2(r$df)) else "normal approximation"
    tagList(
      tags$table(class = "table table-sm",
        tags$tr(tags$th("Parameter"), tags$th("Estimate"), tags$th(sprintf("%d%% CI", round(100 * (1 - r$alpha)))), tags$th("Inference"), tags$th("p")),
        tags$tr(tags$td(h$sym), tags$td(f3(r$est)), tags$td(paste(f3(r$lo), "to", f3(r$hi))), tags$td(dft), tags$td(fp(r$p)))),
      tags$p(tags$b(if (rej) "H0 is rejected. " else "H0 is not rejected. "),
             sprintf("At \u03b1 = %s (%s), %s is estimated at %s %s (%d%% CI %s to %s; %s) from %d participants in %d sites, using a model with %s.",
                     r$alpha, if (r$side == "two") "two-sided" else "one-sided", h$sym, f3(r$est), r$units,
                     round(100 * (1 - r$alpha)), f3(r$lo), f3(r$hi), fpt(r$p), r$N, r$S, r$model)),
      if (!is.null(r$pi)) tags$p(sprintf("Between-site SD of the effect \u03c4 = %s; a new site could observe an effect from %s to %s %s (95%% prediction interval).",
                                         f2(r$tau), f2(r$pi[1]), f2(r$pi[2]), r$units)),
      if (r$ht == "corr") tags$p(sprintf("For comparison, the pooled correlation ignoring site is %s and the correlation of site means is %s.", f2(r$pooled), f2(r$between))),
      {
        ds <- sign(r$est) * abs(r$delta)
        tags$p(sprintf("The effect of interest stated at the planning stage (%s %s, in the direction of the estimate) lies %s the confidence interval.",
                       f2(ds), r$units, if (ds >= r$lo && ds <= r$hi) "inside" else "outside"))
      }
    )
  })

  output$h_plan_obs <- renderTable({
    f <- dg_f(); ht <- input$h_type; r <- R$hyp
    pse <- ppow()
    obs_se <- if (!is.null(r) && isTRUE(r$ok)) r$se else NA
    obs_df <- if (!is.null(r) && isTRUE(r$ok)) r$df else NA
    data.frame(
      Quantity = c("Sites", "Participants per site", "Site ICC of the outcome", "Between-site SD of the effect",
                   "Standard error of the tested parameter", "Minimum detectable effect at the target power"),
      Planned = c(input$p_S, input$p_n, f3(input$p_icc), f2(if (ht %in% c("change", "groupchange")) input$p_taut else input$p_tau),
                  f3(pse$se), if (is.na(pse$se)) "NA" else f3(mde(pse$se, pse$df, input$h_alpha, input$h_side, input$p_power, ht == "corr"))),
      Observed = c(f$S, f2(f$nbar), f3(f$icc), f2(f$tau), f3(obs_se),
                   if (is.na(obs_se)) "NA" else f3(mde(obs_se, obs_df, input$h_alpha, input$h_side, input$p_power, ht == "corr"))),
      check.names = FALSE)
  })

  output$dg_vb_ml <- renderText(if (dg()$ml) { if (dg()$long) "Longitudinal multilevel" else if (dg()$slope) "Multilevel, random slope" else "Multilevel, random intercept" } else "Single-level")
  output$dg_vb_deff <- renderText(f2(dg()$deff_x))
  output$dg_vb_neff <- renderText(sprintf("%s of %s", format(round(dg()$neff), big.mark = ","), format(round(dg()$N), big.mark = ",")))
  output$dg_vb_xl <- renderText(c(site = "Site level", participant = "Participant level", both = "Both levels")[[dg()$xl]])
  output$dg_plot <- renderUI(dg_flow(dg()$tab))

  output$dg_list <- renderUI(dec_list_ui(dg()$tab))

  dg_script <- reactive({
    x <- D()
    covs <- if (is.null(x)) c("age", "female") else intersect(c("age", "female", "dx"), names(x$d_base))
    covs <- setdiff(covs, dg_f()$x)
    dg_code(dg_f(), dg(), covs = covs, instr = if (!is.null(x) && !is.null(x$instrument)) x$instrument else "instrument")
  })
  output$dg_code <- renderText(dg_script())
  output$dg_dl_code <- downloadHandler(function() paste0("analysis_", dg_f()$y, ".R"), function(f) writeLines(dg_script(), f))

  dg_items <- reactive({
    f <- dg_f(); d <- dg(); h <- hyp(); it <- character(0)
    it <- c(it, sprintf("Objective, primary outcome and the hypotheses (H0: %s; H1: %s) with \u03b1 = %s, %s", h$h0, h$h1, input$h_alpha,
                        if (input$h_side == "two") "two-sided" else "one-sided"),
            sprintf("Estimate of %s with its confidence interval, degrees of freedom and p value", h$sym),
            "Definition of the site unit, kept the same in every analysis",
            "Site table: participants, groups, age, sex and site-level conditions by site")
    if (d$ml) it <- c(it, sprintf("For %s: %s sites, mean %s per site, site ICC %s%s, DEFF %s, N_eff %s",
                                  f$y, f$S, f2(f$nbar), f3(f$icc),
                                  if (!is.na(f$icc_lo)) sprintf(" (95%% interval %s to %s)", f3(f$icc_lo), f3(f$icc_hi)) else "",
                                  f2(d$deff_y), round(d$N / d$deff_y)))
    it <- c(it, sprintf("Level of the exposure %s (rho_x = %s) and the degrees-of-freedom method used to test it", f$x, f3(f$rhox)))
    if (d$slope) it <- c(it, "Between-site SD of the effect and the prediction interval for a new site")
    else if (d$xl != "site" && d$ml) it <- c(it, "Test of between-site variation in the effect (tau and its p value)")
    if (isTRUE(f$instr) || isTRUE(f$inst_multi)) it <- c(it, "Instrument and institution components of the site variance")
    if (nz_num(f$n_meas) > 1) it <- c(it, "Within-site correlations between measures, with pooled correlations only alongside them")
    if (isTRUE(f$predict)) it <- c(it, "Prediction accuracy from leave-one-site-out validation")
    if (isTRUE(f$long)) it <- c(it, "Residual correlation structure chosen by AIC/BIC, estimated correlation between visits and reliability of individual slopes")
    if (nz_num(f$miss_sites) > 0 || nz_num(f$miss_pct) > 0) it <- c(it, "Measure availability by site and the imputation model, with site as a level")
    if (d$ml) it <- c(it, "Estimation by REML, Kenward-Roger inference, and any singular fits")
    it
  })
  output$dg_report <- renderUI(tags$ul(lapply(dg_items(), tags$li)))

  output$dg_dl_rec <- downloadHandler("decision_record.txt", function(fl) {
    t <- dg()$tab; f <- dg_f(); h <- hyp(); r <- R$hyp
    txt <- c(sprintf("Decision record: outcome %s, exposure %s", f$y, f$x), "",
             paste("Objective:", input$h_obj),
             sprintf("Hypotheses: H0 %s; H1 %s (alpha = %s, %s)", h$h0, h$h1, input$h_alpha, if (input$h_side == "two") "two-sided" else "one-sided"), "",
             unlist(lapply(seq_len(nrow(t)), function(i) c(
               sprintf("Step %s. %s %s", t$step[i], t$question[i], t$answer[i]),
               paste("  Decision:", if (t$active[i]) t$action[i] else paste(t$action[i], "(not needed)")),
               paste("  Evidence:", t$why[i]),
               if (nzchar(t$warn[i]) && t$active[i]) paste("  Warning:", t$warn[i]), ""))),
             if (!is.null(r) && isTRUE(r$ok)) c("Test of the hypothesis:",
               sprintf("  %s = %s (CI %s to %s), p %s; H0 %s", h$sym, f3(r$est), f3(r$lo), f3(r$hi), fp(r$p),
                       if (r$p < r$alpha) "rejected" else "not rejected"), "") else NULL,
             "To report:", paste("-", dg_items()))
    writeLines(txt, fl)
  })

  ## ---------------------------------------------------------- decisions beside each inspection
  dec_list_ui <- function(t) {
    tagList(lapply(seq_len(nrow(t)), function(i) {
      r <- t[i, ]
      st <- if (!r$active) "ok" else if (nzchar(r$warn)) "warn" else "act"
      dec_item(st, paste0("Step ", r$step, ". ", r$question), tags$span(tags$b(r$answer), tags$br(), if (r$active) r$action else "Not needed"),
               evidence = r$why, rule = if (nzchar(r$warn) && r$active) r$warn else NULL)
    }))
  }
  output$dg_list_plan <- renderUI(dec_list_ui(dg_plan()$tab))

  output$dec_data <- renderUI({
    x <- D()
    if (is.null(x)) return(dec_empty("Load data to see the decisions about the data structure."))
    d <- x$d_base; d$site <- droplevels(d$site); S <- nlevels(d$site); ns <- table(d$site)
    it <- list()
    multi <- names(which(table(unique(d[, c("institution", "site")])$institution) > 1))
    it[[length(it) + 1]] <- if (length(multi))
      dec_item("warn", "Site unit", "Keep the acquisition cohort as the site unit and compare it with the institution on 'What site is made of'.",
               paste("Institutions with more than one cohort:", paste(multi, collapse = ", ")), "Cohorts from one institution may share a scanner.")
    else dec_item("ok", "Site unit", "Use site as the site unit in every analysis.", "Each site is a separate institution.")
    it[[length(it) + 1]] <- dec_item(if (S < 10) "warn" else if (S < 20) "act" else "ok",
      sprintf("%d sites", S),
      if (S < 10) "Treat the site variance as a bound: report ICC intervals, Kenward-Roger inference and a Bayesian fit with a weakly informative prior."
      else if (S < 20) "Report ICC intervals and use Kenward-Roger degrees of freedom for site-level and random-slope comparisons."
      else "The site variance can be estimated with useful precision; report ICC intervals and Kenward-Roger inference.",
      sprintf("Participants per site %d to %d (mean %s).", min(ns), max(ns), f2(mean(ns))),
      "20 to 30 clusters are needed for a stable between-cluster variance.")
    if ("dx" %in% names(d) && !is.null(x$groups)) {
      tb <- table(d$site, d$dx); p <- suppressWarnings(chisq.test(tb)$p.value)
      it[[length(it) + 1]] <- dec_item(if (p >= 0.05) "ok" else "warn", "Group balance across sites",
        if (p >= 0.05) "Groups are compared within every site; the group comparison is a within-site comparison."
        else "The group proportion differs between sites: split the group into within-site and site-mean parts.",
        sprintf("\u03c7\u00b2 = %s, df %d, %s.", f2(suppressWarnings(chisq.test(tb)$statistic)), (nrow(tb) - 1) * (ncol(tb) - 1), fpt(p)))
    }
    if ("female" %in% names(d)) {
      one <- sum(tapply(d$female, d$site, function(z) length(unique(na.omit(z))) == 1))
      if (one > 0) it[[length(it) + 1]] <- dec_item("act", "Sex is partly a site property",
        "Enter sex with its site mean as well as the participant value, or adjust within sites.",
        sprintf("%d of %d sites recruited one sex only.", one, S))
    }
    if (length(x$expo)) it[[length(it) + 1]] <- dec_item("act", "Site-level conditions present",
      paste0("Analyse ", paste(x$expo_labels[x$expo], collapse = ", "), " as site-level exposures, tested against between-site variation (Exposure level tab)."))
    miss <- sapply(seq_len(nrow(x$meas)), function(k) {
      v <- x$meas$var[k]; z <- d
      if (!x$meas$both[k] && "dx" %in% names(z)) z <- z[!is.na(z$dx) & z$dx == 1, ]
      sum(tapply(z[[v]], droplevels(z$site), function(q) all(is.na(q))), na.rm = TRUE)
    })
    if (any(miss > 0)) it[[length(it) + 1]] <- dec_item("warn", "Measures absent at whole sites",
      "Treat these as structural missingness: compare measures on the same sites and include site as a level in any imputation model.",
      paste(paste0(x$meas$label[miss > 0], " (", miss[miss > 0], ")"), collapse = "; "))
    else it[[length(it) + 1]] <- dec_item("ok", "Measure availability", "Every measure is recorded at every site.")
    tagList(it)
  })

  output$dec_site <- renderUI({
    r <- R$site
    if (is.null(r)) return(dec_empty("Run the site layer to see which measures need a random site effect."))
    need <- r$deff >= 1.1 | r$icc >= 0.01
    it <- list(dec_item("act", sprintf("Random site effect for %d of %d measures", sum(need), nrow(r)),
      "Fit a random site intercept and report the ICC with its interval, mean participants per site, DEFF and N_eff for these measures.",
      paste(r$label[need], collapse = ", "), "An ICC is negligible only when the design effect is below about 1.1."))
    if (any(!need)) it[[length(it) + 1]] <- dec_item("ok", "Single-level analysis acceptable",
      paste("Site ICC near zero and DEFF below 1.1:", paste(r$label[!need], collapse = ", "), ". Report the ICC that justifies it."))
    k <- which.max(r$deff)
    it[[length(it) + 1]] <- dec_item("stop", paste("Largest dependence:", r$label[k]),
      sprintf("Avoid a single-level analysis: it treats %d observations as independent, but they carry the information of %d.", r$N[k], round(r$neff[k])),
      sprintf("ICC %s, DEFF %s.", f3(r$icc[k]), f2(r$deff[k])))
    meas_v <- !is.na(r$icc_adj) & r$icc > 0.01 & r$icc_adj >= 0.8 * r$icc
    if (any(meas_v)) it[[length(it) + 1]] <- dec_item("act", "Site variance reflects measurement",
      "Adjusting for who was recruited did not remove it: model site even after covariate adjustment.",
      paste(r$label[meas_v], collapse = ", "), "Adjusted ICC at least 80% of the empty-model ICC.")
    comp <- !is.na(r$icc_adj) & r$icc > 0.01 & r$icc_adj < 0.5 * r$icc
    if (any(comp)) it[[length(it) + 1]] <- dec_item("ok", "Site variance reflects recruitment",
      "Covariates explain most of the site variance: keep the within-site and site-mean split of the covariates.",
      paste(r$label[comp], collapse = ", "))
    if (any(r$singular)) it[[length(it) + 1]] <- dec_item("warn", "Singular fits",
      "Report the site variance as estimated at zero; do not drop site silently.", paste(r$label[r$singular], collapse = ", "))
    if (any(r$sites < 20)) it[[length(it) + 1]] <- dec_item("warn", "Fewer than 20 sites behind some measures",
      "Rank measures by site dependence only where their intervals separate.",
      paste(paste0(r$label[r$sites < 20], " (", r$sites[r$sites < 20], ")"), collapse = "; "))
    tagList(it)
  })

  output$dec_expo <- renderUI({
    if (is.null(D())) return(dec_empty("Load data to classify the variables."))
    rr <- rho_df(); rr <- rr[!is.na(rr$rho), ]
    it <- lapply(seq_len(nrow(rr)), function(i) {
      v <- rr[i, ]
      if (v$rho >= 0.8) dec_item("stop", paste(v$label, ": site property"),
                                 "Test it against between-site variation with Kenward-Roger df; never with a single-level model.", sprintf("\u03c1x = %s", f3(v$rho)))
      else if (v$rho > 0.1) dec_item("act", paste(v$label, ": within and between sites"),
                                     "Enter it as a within-site deviation and a site mean.", sprintf("\u03c1x = %s", f3(v$rho)))
      else dec_item("ok", paste(v$label, ": participant property"), "Enter it as a participant-level term.", sprintf("\u03c1x = %s", f3(v$rho)))
    })
    e <- R$expo
    if (!is.null(e)) {
      a <- e$p[1] < 0.05; b <- e$p[2] < 0.05
      it[[length(it) + 1]] <- if (a && !b) dec_item("stop", "Single-level result is a site difference",
        "Report the multilevel estimate; the single-level effect compares sites, not the exposure.",
        sprintf("Single-level %s SD (%s); random intercept %s SD (%s).", f2(e$est[1]), fpt(e$p[1]), f2(e$est[2]), fpt(e$p[2])))
      else if (b) dec_item("act", "Effect holds against between-site variation",
        sprintf("Report %s SD with Kenward-Roger df %s.", f2(e$est[2]), f2(e$df[2])), fpt(e$p[2]))
      else dec_item("ok", "No effect in either analysis", sprintf("Report the multilevel estimate %s SD (%s).", f2(e$est[2]), fpt(e$p[2])))
    }
    o <- R$null
    if (!is.null(o)) {
      s <- o[o$assignment == "Assigned to whole sites", ]
      nv <- s$rate[s$analysis == "Site ignored"]; mv <- s$rate[s$analysis == "Random site intercept"]
      it[[length(it) + 1]] <- dec_item("stop", "Ignoring site under a known null",
        sprintf("False-positive rates of %d%% to %d%% for site-level exposures: a single-level test is not valid.", round(100 * min(nv)), round(100 * max(nv))),
        sprintf("%d assignments per outcome.", o$B[1]))
      it[[length(it) + 1]] <- dec_item(if (max(mv) <= 0.08) "ok" else "warn", "Random site intercept under the same null",
        sprintf("Rates of %d%% to %d%% against the nominal 5%%.", round(100 * min(mv)), round(100 * max(mv))))
    }
    tagList(it)
  })

  output$dec_dx <- renderUI({
    r <- R$dx
    if (is.null(r)) return(dec_empty("Fit the four analyses to decide between a random intercept and a random slope."))
    t <- r$tab; it <- list()
    it[[1]] <- if (r$p_slope < 0.05) dec_item("act", "Retain the random slope",
      sprintf("Report %s SD (SE %s, %s df) with the prediction interval %s to %s SD for a new site.", f2(t$est[4]), f3(t$se[4]), f2(t$df[4]), f2(r$pi[1]), f2(r$pi[2])),
      sprintf("\u03c4 = %s SD, likelihood-ratio %s.", f2(r$tau), fpt(r$p_slope)), "Test of \u03c4\u00b2 = 0 against a 50:50 mixture of \u03c7\u00b2 with 1 and 2 df.")
    else dec_item("ok", "Random intercept is sufficient",
      sprintf("Report %s SD (SE %s) and \u03c4 = %s with its test.", f2(t$est[3]), f3(t$se[3]), f2(r$tau)), fpt(r$p_slope))
    ratio <- t$se[4] / t$se[3]
    if (r$p_slope < 0.05 && ratio > 1.3) it[[length(it) + 1]] <- dec_item("stop", "Site-mean adjustment understates uncertainty",
      sprintf("Fixed site terms, a random intercept or ComBat-type harmonisation give an SE %s times too small here.", f2(ratio)),
      sprintf("SE %s (random intercept) against %s (random slope).", f3(t$se[3]), f3(t$se[4])))
    if (abs(t$se[1] - t$se[3]) / t$se[3] < 0.1) it[[length(it) + 1]] <- dec_item("ok", "Group balanced within sites",
      "The random intercept barely changes the SE of this within-site comparison.", sprintf("SE %s single-level, %s random intercept.", f3(t$se[1]), f3(t$se[3])))
    if (isTRUE(r$singular)) it[[length(it) + 1]] <- dec_item("warn", "Singular random-slope fit", "Report the slope variance at the boundary.")
    p <- R$ps
    if (!is.null(p)) {
      pl <- p$pool
      it[[length(it) + 1]] <- dec_item(if (pl$pi_lo < 0 && pl$pi_hi > 0) "warn" else "ok", "Two-stage check across sites",
        if (pl$pi_lo < 0 && pl$pi_hi > 0) "A new site could find a difference in either direction: report the prediction interval with the average."
        else "Every new site is expected to show a difference in the same direction.",
        sprintf("I\u00b2 = %d%%; prediction interval %s to %s SD.", round(100 * pl$I2), f2(pl$pi_lo), f2(pl$pi_hi)))
    }
    tagList(it)
  })

  output$dec_cx <- renderUI({
    c <- R$cx; i <- R$inst; it <- list()
    if (is.null(c) && is.null(i)) return(dec_empty("Decompose site to choose the site unit and the crossed terms."))
    if (!is.null(c)) it[[length(it) + 1]] <- if (!is.na(c$p_instrument[2]) && c$p_instrument[2] < 0.05)
      dec_item("act", "Add the instrument as a crossed random effect",
               sprintf("Use (1 | site) + (1 | instrument): the site ICC falls from %s to %s and the instrument carries %s.", f3(c$site_icc[1]), f3(c$site_icc[2]), f3(c$instrument_icc[2])),
               fpt(c$p_instrument[2]))
      else dec_item("ok", "Instrument adds no variance", "Keep (1 | site) without the instrument term.", fpt(c$p_instrument[2]))
    if (identical(i, "none")) it[[length(it) + 1]] <- dec_item("ok", "Institution and cohort coincide", "Keep site as the site unit.")
    else if (!is.null(i)) for (k in seq_len(nrow(i))) {
      tot <- i$institution[k] + i$cohort_within[k]
      lab <- D()$meas$label[match(i$var[k], D()$meas$var)]
      it[[length(it) + 1]] <- if (tot > 0 && i$cohort_within[k] < 0.25 * tot)
        dec_item("act", paste(lab, ": institution carries the variance"), "Use institution as the site unit; cohorts from one institution share their effects.",
                 sprintf("Institution %s, cohort within institution %s.", f3(i$institution[k]), f3(i$cohort_within[k])))
      else dec_item("ok", paste(lab, ": cohort carries the variance"), "Keep cohort as the site unit.",
                    sprintf("Institution %s, cohort within institution %s.", f3(i$institution[k]), f3(i$cohort_within[k])))
    }
    tagList(it)
  })

  output$dec_mv <- renderUI({
    m <- R$mv
    if (is.null(m)) return(dec_empty("Estimate the correlations to decide how measures are combined."))
    t <- m$tab
    it <- lapply(seq_len(nrow(t)), function(j) {
      w <- ifelse(is.na(t$within_model[j]), t$within_centred[j], t$within_model[j])
      b <- ifelse(is.na(t$between_model[j]), t$between_sitemeans[j], t$between_model[j])
      ev <- sprintf("Within %s, between %s, pooled %s.", f2(w), f2(b), f2(t$pooled[j]))
      if (sign(w) != sign(b) && abs(b - w) > 0.1) dec_item("stop", paste(t$l1[j], "\u2013", t$l2[j]),
        "The sign reverses between levels: report the within-site correlation; the pooled value must not stand alone.", ev)
      else if (abs(b - w) > 0.2) dec_item("warn", paste(t$l1[j], "\u2013", t$l2[j]),
        "Within and between differ: report the within-site correlation with the pooled one alongside.", ev)
      else dec_item("ok", paste(t$l1[j], "\u2013", t$l2[j]), "The levels agree; report the within-site correlation.", ev)
    })
    it[[length(it) + 1]] <- dec_item("act", "Combining these measures",
      "Centre each measure on its site mean before correlating or predicting, and validate any prediction model by leaving whole sites out.")
    if (!is.null(m$msg)) it[[length(it) + 1]] <- dec_item("warn", "Multivariate model", m$msg)
    tagList(it)
  })

  output$dec_time <- renderUI({
    v <- vc(); tt <- times(); k <- length(tt); it <- list()
    adj <- sapply(seq_len(k - 1), function(i) v$cor[i, i + 1]); far <- v$cor[1, k]
    it[[1]] <- dec_item(if (min(adj) < 0.5) "warn" else "ok", "Correlation between visits",
      if (min(adj) < 0.5) "Much of the observed change is occasion noise: individual slopes will be shrunk strongly; report their reliability."
      else "Adjacent visits are strongly related: participant random intercept and slope are supported.",
      sprintf("Adjacent %s to %s; first to last %s.", f2(min(adj)), f2(max(adj)), f2(far)))
    if (input$su2 > 0) it[[length(it) + 1]] <- dec_item("act", "Keep the site random effect",
      "Site adds \u03c3\u00b2u to every covariance between visits; an analysis of change without site attributes site shifts to participants.")
    if (input$sb01 < 0) it[[length(it) + 1]] <- dec_item("ok", "Negative intercept-slope covariance", "Participants starting higher change more slowly; keep the covariance in D unstructured.")
    for (f in list(list(R$ls, "Simulated data"), list(R$lu, "Uploaded data"))) if (!is.null(f[[1]]) && !is.null(f[[1]]$ar)) {
      a <- f[[1]]
      it[[length(it) + 1]] <- dec_item("act", paste(f[[2]], ": residual structure"),
        if (a$ar["AIC"] < a$ind["AIC"]) sprintf("Use autoregressive residuals (\u03c6 = %s).", f2(a$ar["phi"])) else "Use independent residuals.",
        sprintf("AIC %s independent, %s autoregressive.", f2(a$ind["AIC"]), f2(a$ar["AIC"])))
    }
    x <- D()
    if (!is.null(x) && !x$longitudinal) it[[length(it) + 1]] <- dec_item("warn", "Loaded data are single-visit",
      "The time layer cannot be estimated from these data; plan it with the calculator or upload longitudinal data.")
    tagList(it)
  })

  output$dec_plan <- renderUI({
    p <- pl(); se <- function(S, n) sqrt((input$pl_tau^2 + 4 * input$pl_se2 / n) / S)
    gainS <- se(input$pl_S, input$pl_n) - se(2 * input$pl_S, input$pl_n)
    gainN <- se(input$pl_S, input$pl_n) - se(input$pl_S, 2 * input$pl_n)
    tagList(
      dec_item(if (p$deff >= 1.1) "act" else "ok", "Plan with the effective sample size",
               if (p$deff >= 1.1) sprintf("Base power on N_eff = %s, not N = %s.", round(p$neff), p$N) else "The design effect is close to 1 for this exposure.",
               sprintf("DEFF %s.", f2(p$deff))),
      dec_item(if (gainS > 1.5 * gainN) "act" else "ok", "Where to put extra participants",
               if (gainS > 1.5 * gainN) "Add sites rather than participants per site." else "Adding participants per site is about as efficient as adding sites.",
               sprintf("Doubling sites lowers the SE by %s; doubling participants per site by %s.", f3(gainS), f3(gainN)),
               "Equation (17): \u03c4\u00b2/S does not fall when participants are added within sites."),
      if (input$pl_S < 20) dec_item("warn", sprintf("%d sites planned", input$pl_S), "Plan Kenward-Roger inference and report the ICC as an interval."))
  })

  output$dec_eta <- renderUI({
    e <- eta()
    if (is.null(e) || !nrow(e)) return(dec_empty("Convert reported \u03b7\u00b2 values to see the decisions."))
    ov <- max(e$lo) <= min(e$hi)
    tagList(
      dec_item("act", "Report the ICC, not \u03b7\u00b2",
               sprintf("\u03b7\u00b2 includes a chance component of about %s with these sites.", formatC(mean(e$chance), format = "f", digits = 4))),
      if (any(e$k < 10)) dec_item("warn", sprintf("%d sites", max(e$k)), "Present the ICC as a bound with its interval.",
                                  sprintf("Upper limits reach %s.", f3(max(e$hi)))),
      dec_item(if (ov) "warn" else "ok", "Ranking outcomes by site dependence",
               if (ov) "The intervals overlap: do not rank the outcomes." else "At least two outcomes separate; rank only those."),
      dec_item(if (any(e$deff >= 1.1)) "act" else "ok", "Model site under the harmonised protocol",
               sprintf("Design effects of %s to %s at %s to %s participants per site.", f2(min(e$deff)), f2(max(e$deff)), round(min(e$nbar)), round(max(e$nbar)))))
  })

}

shinyApp(ui, server)
