# ============================================================================
# LOIO (leave-one-item-out) R^2 add-on for the existing simulation pipeline.
# ----------------------------------------------------------------------------
# WHY THIS IS A SEPARATE, CHEAP SCRIPT AND NOT A RE-RUN OF THE FULL PIPELINE
# ----------------------------------------------------------------------------
# `process_one_cell()` in your main script computes the raw item data and each
# method's item partition internally, but never saves them -- only the summary
# row (R^2, ARI, Selected_K, ...) survives to the .rds/.csv. That's confirmed
# by inspecting your saved files directly: 1200 rows x 26 columns, ~210 KB --
# far too small to contain any raw item matrices. So LOIO cannot be computed
# from the files you already have; it needs a re-run.
#
# BUT: LOIO only needs the item partition from the TRUE-k1 fit, which is
# exactly what evaluate_efa() and evaluate_E_ReMI() already produce. The three
# expensive parts of your pipeline exist for OTHER purposes and are skipped
# entirely here:
#   - select_k_efa()          (50 parallel-analysis reference iterations)  -- SKIPPED
#   - select_k_ereMI_permnull() (Q-grid x 50-250 permutations each)        -- SKIPPED
#   - evaluate_efa_boot()      (50-500 full bootstrap refits)              -- SKIPPED
# What's left per replication: one psych::fa() call, and one E-ReMI fit at
# Q = k1 (your existing Nruns random starts, permutations trimmed to the
# minimum the function accepts since we don't need the significance test --
# see the PERMUTATIONS_FOR_FIT note below). That should be a small fraction
# of your original per-cell runtime.
#
# HONESTY NOTE: this generates FRESH data (new random draws), not the exact
# same replications behind your saved summary numbers. Re-creating those
# bit-for-bit would require your exact R/package versions and an unmodified
# future.apply RNG stream, which is too fragile to rely on. Treat this as an
# independent supplementary check run under the same design, not a recomputation
# of the existing rows -- which is a perfectly normal thing to report as a
# robustness check, just describe it that way in the manuscript.
# ============================================================================

# ---------- 0. Setup (same as your main script) -----------------------------
library(future)
library(future.apply)
library(progressr)
library(psych)
library(MASS)
library(mclust)
library(pracma)
library(optimbase)
library(Rfast)
library(reshape)

set.seed(12345)
RNGkind("L'Ecuyer-CMRG")

# Path setup -- uses your R session's current working directory (getwd()),
# which the error you hit confirms is where the source files actually live:
#   "C:/Users/ahmed/Desktop/IfADo Projects/Third Project based on Survey and
#    DVS Data/BOOTEFA-Server-number cruncher"
# The old '~/Documents/BOOTEFA-Server-number cruncher/' broke because on this
# Windows R session '~' ALREADY expands to "Documents" (R's default on
# Windows), so '~/Documents/...' doubled it into ".../Documents/Documents/...".
# getwd() sidesteps that entirely -- just make sure you've either opened this
# script from that folder, or run setwd("C:/Users/ahmed/Desktop/IfADo
# Projects/Third Project based on Survey and DVS Data/BOOTEFA-Server-number
# cruncher") once beforehand.
path <- paste0(getwd(), "/")

required_source_files <- c(
  'Randompartition_function.R', 'Unequal_Randompartition_Function.R',
  'Log_Likelihood_function_REMAXINT.R', 'Log_Likelihood_function_E_ReMI.R',
  'Log_LR_Test_Statistic.R', 'Update_column_clusters.R',
  'Update_row_clusters_REMAXINT.R', 'Update_row_clusters_E_ReMI.R',
  'Update_G_Omega.R', 'REMAXINT.R', 'E_ReMI.R', 'Permutation_Function.R',
  'Data_Gen_Function.R'
)
missing_files <- required_source_files[!file.exists(paste0(path, required_source_files))]
if (length(missing_files) > 0) {
  stop("Could not find these source file(s) in '", path, "':\n  ",
       paste(missing_files, collapse = "\n  "),
       "\n\nRun getwd() and confirm it points at the folder containing these .R files ",
       "(setwd(...) there first if not), then re-source this script.", call. = FALSE)
}

source(paste(path, 'Randompartition_function.R', sep = ''))
source(paste(path, 'Unequal_Randompartition_Function.R', sep = ''))
source(paste(path, 'Log_Likelihood_function_REMAXINT.R', sep = ''))
source(paste(path, 'Log_Likelihood_function_E_ReMI.R', sep = ''))
source(paste(path, 'Log_LR_Test_Statistic.R', sep = ''))
source(paste(path, 'Update_column_clusters.R', sep = ''))
source(paste(path, 'Update_row_clusters_REMAXINT.R', sep = ''))
source(paste(path, 'Update_row_clusters_E_ReMI.R', sep = ''))
source(paste(path, 'Update_G_Omega.R', sep = ''))
source(paste(path, 'REMAXINT.R', sep = ''))
source(paste(path, 'E_ReMI.R', sep = ''))
source(paste(path, 'Permutation_Function.R', sep = ''))
source(paste(path, 'Data_Gen_Function.R', sep = ''))

worker_source_files <- c(
  paste(path, 'Randompartition_function.R', sep = ''),
  paste(path, 'Unequal_Randompartition_Function.R', sep = ''),
  paste(path, 'Log_Likelihood_function_REMAXINT.R', sep = ''),
  paste(path, 'Log_Likelihood_function_E_ReMI.R', sep = ''),
  paste(path, 'Log_LR_Test_Statistic.R', sep = ''),
  paste(path, 'Update_column_clusters.R', sep = ''),
  paste(path, 'Update_row_clusters_REMAXINT.R', sep = ''),
  paste(path, 'Update_row_clusters_E_ReMI.R', sep = ''),
  paste(path, 'Update_G_Omega.R', sep = ''),
  paste(path, 'REMAXINT.R', sep = ''),
  paste(path, 'E_ReMI.R', sep = ''),
  paste(path, 'Permutation_Function.R', sep = ''),
  paste(path, 'Data_Gen_Function.R', sep = '')
)
worker_packages <- c("psych", "MASS", "mclust", "pracma", "optimbase", "Rfast", "reshape")

# Same fix as `path` above -- write outputs to the current working directory
# rather than a hardcoded, environment-dependent '~/Documents/...' path.
output_dir <- getwd()

# ---------- 1. Reused verbatim from your main script -------------------------
generate_true_structure <- function(p, k1, type = "Block") {
  true_Lambda <- matrix(0, nrow = p, ncol = k1)
  if (type == "Block") {
    base_size <- floor(p / k1); remainder <- p %% k1; start <- 1
    true_labels <- numeric(p)
    for (i in 1:k1) {
      block_size <- base_size + ifelse(i <= remainder, 1, 0)
      end <- start + block_size - 1
      true_Lambda[start:end, i] <- runif(block_size, 0.5, 0.7)
      true_labels[start:end] <- i
      start <- end + 1
    }
  } else if (type == "Diagonal") {
    for (i in 1:p) {
      factor_index <- ((i - 1) %% k1) + 1
      true_Lambda[i, factor_index] <- runif(1, 0.5, 0.7)
    }
    true_labels <- ((1:p - 1) %% k1) + 1
  }
  list(true_Lambda = true_Lambda, true_labels = true_labels)
}

simulate_data <- function(n, true_Lambda, non_normal = FALSE,
                          outliers = FALSE, psi_value = 0.3) {
  p <- nrow(true_Lambda)
  psi <- diag(p) * psi_value
  Sigma <- true_Lambda %*% t(true_Lambda) + psi
  latent_data <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  cut_points <- qnorm(seq(0, 1, length.out = 6))
  likert_data <- apply(latent_data, 2, function(x) {
    as.numeric(cut(x, breaks = cut_points, labels = FALSE))
  })
  if (non_normal) {
    n_bias <- ceiling(n * 0.65)
    bias_rows <- sample(1:n, n_bias, replace = FALSE)
    for (i in bias_rows) {
      shift_amounts <- sample(1:3, size = p, replace = TRUE, prob = c(.4, .4, .2))
      likert_data[i, ] <- pmin(likert_data[i, ] + shift_amounts, 5)
    }
  }
  if (outliers) {
    n_out_cells <- ceiling(n * p * 0.10)
    outlier_cells <- sample(1:(n * p), n_out_cells, replace = FALSE)
    likert_data[outlier_cells] <- sample(c(1, 5), n_out_cells, replace = TRUE)
    n_out_rows <- ceiling(n * 0.05)
    outlier_rows <- sample(1:n, n_out_rows, replace = FALSE)
    for (i in outlier_rows) {
      likert_data[i, ] <- sample(c(1, 5), size = 1)
    }
  }
  likert_data <- pmax(pmin(likert_data, 5), 1)
  return(likert_data)
}

# ---------- 2. LOIO core: identical formula applied to BOTH methods ---------
# For item i in group g, correlate item i with the mean of the OTHER items in
# g (excluding i), square it, average across items. No refitting, no model
# selection -- just the item data plus a group-assignment vector.
compute_loio_r2 <- function(item_data, group_assignment) {
  item_data <- as.matrix(item_data)
  groups <- unique(group_assignment)
  loio_r2 <- rep(NA_real_, ncol(item_data))
  for (g in groups) {
    idx <- which(group_assignment == g)
    if (length(idx) < 2) next  # singleton group: no "other items" to average
    for (i in idx) {
      others <- setdiff(idx, i)
      loio_mean <- if (length(others) == 1) item_data[, others] else rowMeans(item_data[, others, drop = FALSE])
      loio_r2[i] <- cor(item_data[, i], loio_mean)^2
    }
  }
  mean(loio_r2, na.rm = TRUE)
}

# ---------- 3. Trimmed per-cell worker ---------------------------------------
# PERMUTATIONS_FOR_FIT: Permutation_Function's `permutations` argument controls
# the size of its null distribution for the significance test -- it does not
# change the OBSERVED-data fit (K_EReMI) we actually need. Set to a small
# value to skip almost all of that overhead. If you'd rather bypass this
# argument entirely, call your core E_ReMI() fitting function directly instead
# of going through Permutation_Function -- share that file's signature and I
# can wire it in for a further speedup.
PERMUTATIONS_FOR_FIT <- 2

process_one_cell_loio <- function(sim, p, k1, cond_idx, conditions,
                                  n, P, Nruns, psi_value) {
  cond <- conditions[cond_idx, ]
  
  obj <- generate_true_structure(p, k1, cond$Structure)
  true_Lambda <- obj$true_Lambda
  data <- simulate_data(n, true_Lambda, cond$NonNormal, cond$Outliers, psi_value)
  
  rows <- list()
  
  # ---- EFA at true k1 ----
  out_efa <- tryCatch({
    poly_cor <- psych::polychoric(data)$rho
    ft <- psych::fa(r = poly_cor, nfactors = k1, n.obs = nrow(data),
                    rotate = "oblimin", fm = "minres", cor = "poly")
    L1 <- apply(as.matrix(ft$loadings), 1, function(x) which.max(abs(x)))
    list(loio_r2 = compute_loio_r2(data, L1))
  }, error = function(e) NULL)
  if (!is.null(out_efa)) {
    rows[[length(rows) + 1]] <- data.frame(
      Sim = sim, SampleSize = n, Items = p, TrueClusters = k1,
      Method = "EFA", LOIO_R_squared = out_efa$loio_r2, stringsAsFactors = FALSE)
  }
  
  # ---- E-ReMI at true k1 (Q = k1), permutation test trimmed to the minimum ----
  out_eremi <- tryCatch({
    pf <- Permutation_Function(data, P, k1, Nruns, PERMUTATIONS_FOR_FIT, 0.05)
    K <- pf$K_EReMI
    if (is.null(K) || !is.matrix(K)) NULL else {
      grp <- apply(K, 1, function(x) which.max(abs(x)))
      list(loio_r2 = compute_loio_r2(data, grp))
    }
  }, error = function(e) NULL)
  if (!is.null(out_eremi)) {
    rows[[length(rows) + 1]] <- data.frame(
      Sim = sim, SampleSize = n, Items = p, TrueClusters = k1,
      Method = "E-ReMI", LOIO_R_squared = out_eremi$loio_r2, stringsAsFactors = FALSE)
  }
  
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

# ---------- 4. Driver (same future.apply pattern as your main script) -------
run_loio_for_single_n <- function(n, session_id, Nruns = 20, n_sim = 100,
                                  ncores = 120, output_dir = output_dir,
                                  worker_source_files = worker_source_files,
                                  worker_packages = worker_packages) {
  item_levels     <- c(30, 60)
  cluster_levels1 <- c(4, 6)
  P <- 6
  psi_value <- 0.3
  conditions <- expand.grid(Structure = "Block", NonNormal = FALSE, Outliers = FALSE,
                            stringsAsFactors = FALSE)
  all_cells <- expand.grid(sim = 1:n_sim, p = item_levels, k1 = cluster_levels1,
                           cond_idx = 1:nrow(conditions))
  message(sprintf("LOIO run: %d cells (n=%d, ncores=%d)", nrow(all_cells), n, ncores))
  
  future::plan(future::multisession, workers = ncores)
  on.exit(future::plan(future::sequential), add = TRUE)
  
  t0 <- Sys.time()
  results_list <- future.apply::future_lapply(
    seq_len(nrow(all_cells)),
    future.seed = TRUE,
    future.packages = worker_packages,
    FUN = function(i) {
      if (!isTRUE(getOption("loio_worker_ready", FALSE))) {
        for (f in worker_source_files) if (file.exists(f)) try(source(f), silent = TRUE)
        options(loio_worker_ready = TRUE)
      }
      row <- all_cells[i, ]
      tryCatch(
        process_one_cell_loio(sim = row$sim, p = row$p, k1 = row$k1,
                              cond_idx = row$cond_idx, conditions = conditions,
                              n = n, P = P, Nruns = Nruns, psi_value = psi_value),
        error = function(e) NULL)
    }
  )
  message("LOIO elapsed: ", format(round(difftime(Sys.time(), t0, units = "mins"), 1)))
  
  results_list <- Filter(Negate(is.null), results_list)
  if (length(results_list) == 0) stop("All LOIO cells failed -- check worker errors.")
  results <- do.call(rbind, results_list)
  
  base_name <- paste0("loio_results_n", n, "_session", session_id)
  write.csv(results, file.path(output_dir, paste0(base_name, ".csv")), row.names = FALSE)
  saveRDS(results, file.path(output_dir, paste0(base_name, ".rds")))
  invisible(results)
}

# ---------- 5. Run it --------------------------------------------------------
# Start small (n_sim=20-30) to confirm it works and to get a runtime estimate
# before committing to matching your pilot's n_sim=100 or production's n_sim=500.
loio_pilot <- run_loio_for_single_n(
  n = 100, session_id = "loio_pilot",
  Nruns = 5, n_sim = 5, ncores = 8,
  output_dir = output_dir,
  worker_source_files = worker_source_files,
  worker_packages = worker_packages
)
print(head(loio_pilot))

# ---------- 6. Merge with your EXISTING summary results for reporting -------
# (Run after you have both the original sim_results_*.rds and a loio_results_*.rds
#  for the sample size(s) you care about. Reads directly from output_dir.)
merge_loio_with_existing <- function(n, session_orig, session_loio, output_dir) {
  orig <- readRDS(file.path(output_dir, sprintf("sim_results_n%d_session%s_reduced.rds", n, session_orig)))
  loio <- readRDS(file.path(output_dir, sprintf("loio_results_n%d_session%s.rds", n, session_loio)))
  
  # Original R^2 summarized per (Items, TrueClusters, Method); LOIO summarized
  # the same way (LOIO uses fresh replications, so join on DESIGN CELL, not Sim)
  orig_summary <- aggregate(R_squared_true ~ Items + TrueClusters + Method, orig, mean)
  loio_summary <- aggregate(LOIO_R_squared ~ Items + TrueClusters + Method, loio, mean)
  merge(orig_summary, loio_summary, by = c("Items", "TrueClusters", "Method"), all = TRUE)
}
# Example:
# comparison <- merge_loio_with_existing(500, "pilot", "loio_pilot", output_dir)
# print(comparison)


# ============================================================================
# LOIO R^2 results: figures 
# ============================================================================
for (pkg in c("ggplot2", "dplyr", "tidyr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required. Install it with: install.packages(\"", pkg, "\")",
         call. = FALSE)
}
suppressPackageStartupMessages(library(ggplot2))
`%>%` <- dplyr::`%>%`

data_dir  <- getwd()
fig_dir   <- file.path(data_dir, "figures_highres")
efa_label <- "PCEFA"   # label shown in figures/legends; set to "EFA" if you prefer
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Internally the methods are always called "EFA" and "E-ReMI" (as in your files);
# only the legend/axis text uses efa_label.
method_cols <- c("EFA" = "#1F78B4", "E-ReMI" = "#D95F02")
method_labs <- c("EFA" = efa_label, "E-ReMI" = "E-ReMI")

# ---------- 1. Helpers to read files safely ------------------------------------
to_num <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))
}

read_any <- function(f) {
  if (grepl("\\.rds$", f, ignore.case = TRUE)) return(as.data.frame(readRDS(f)))
  first <- readLines(f, n = 1, warn = FALSE)
  semi  <- grepl(";", first, fixed = TRUE) && !grepl(",", first, fixed = TRUE)
  rd    <- if (semi) utils::read.csv2 else utils::read.csv
  df    <- rd(f, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", check.names = FALSE)
  names(df) <- trimws(names(df))
  df
}

find_files <- function(strict, relaxed) {
  f <- list.files(data_dir, pattern = strict, full.names = TRUE, ignore.case = TRUE)
  if (length(f) == 0) f <- list.files(data_dir, pattern = relaxed, full.names = TRUE, ignore.case = TRUE)
  f
}

need_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0)
    stop(what, " is missing column(s): ", paste(miss, collapse = ", "),
         ".\nColumns found: ", paste(names(df), collapse = ", "), call. = FALSE)
}

tidy_method <- function(df) {
  df$Method <- trimws(as.character(df$Method))
  df$Method[df$Method == "PCEFA"] <- "EFA"
  df
}

key_cols <- c("SampleSize", "Items", "TrueClusters")

# ---------- 2. Read LOIO results ----------------------------------------------
loio_files <- find_files("loio_results_n(100|300|500)_sessionloio_pilot\\.csv$", "loio_results_n[0-9]+.*\\.csv$")
if (length(loio_files) == 0)
  loio_files <- find_files("loio_results_n(100|300|500)_sessionloio_pilot\\.rds$", "loio_results_n[0-9]+.*\\.rds$")
if (length(loio_files) == 0)
  stop("No loio_results_n*_sessionloio_pilot files found in:\n  ", data_dir,
       "\nFiles in this folder: ", paste(utils::head(list.files(data_dir), 15), collapse = ", "),
       "\nUse setwd() to point R at the folder that contains them.", call. = FALSE)
message("LOIO files used: ", paste(basename(loio_files), collapse = ", "))

read_loio <- function(f) {
  df <- read_any(f)
  need_cols(df, c("Sim", key_cols, "Method", "LOIO_R_squared"), basename(f))
  df <- df[, c("Sim", key_cols, "Method", "LOIO_R_squared")]
  for (v in c("Sim", key_cols, "LOIO_R_squared")) df[[v]] <- to_num(df[[v]])
  tidy_method(df)
}
loio <- dplyr::bind_rows(lapply(loio_files, read_loio))
loio <- loio[!is.na(loio$LOIO_R_squared) & loio$Method %in% c("EFA", "E-ReMI"), ]
loio <- loio[!duplicated(loio[, c("Sim", key_cols, "Method")]), ]
if (nrow(loio) == 0) stop("No usable LOIO rows (need Method = EFA/PCEFA or E-ReMI).", call. = FALSE)

n_levels <- sort(unique(loio$SampleSize))
p_levels <- sort(unique(loio$Items))
k_levels <- sort(unique(loio$TrueClusters))

add_labels <- function(df) {
  df$p_lab <- factor(paste0("p = ", df$Items),        levels = paste0("p = ", p_levels))
  df$k_lab <- factor(paste0("k = ", df$TrueClusters), levels = paste0("k = ", k_levels))
  df$n_lab <- factor(paste0("n = ", df$SampleSize),   levels = paste0("n = ", n_levels))
  df$Method <- factor(df$Method, levels = c("EFA", "E-ReMI"))
  df
}
loio  <- add_labels(loio)
n_rep <- max(table(loio$Method, loio$SampleSize, loio$Items, loio$TrueClusters))

# ---------- 3. Read main-simulation results (optional) -------------------------
sim <- NULL
sim_files <- find_files("sim_results_n(100|300|500)_sessionpilot_reduced\\.csv$", "sim_results_n[0-9]+.*\\.csv$")
if (length(sim_files) == 0)
  sim_files <- find_files("sim_results_n(100|300|500)_sessionpilot_reduced\\.rds$", "sim_results_n[0-9]+.*\\.rds$")
if (length(sim_files) > 0) {
  sim <- tryCatch({
    need <- c(key_cols, "Method", "ARI_true", "R_squared_true", "R_squared_true_global")
    one <- function(f) {
      df <- read_any(f); need_cols(df, need, basename(f))
      keep <- intersect(c(need, "Structure", "NonNormal", "Outliers"), names(df))
      df <- df[, keep, drop = FALSE]
      for (v in c(key_cols, "ARI_true", "R_squared_true", "R_squared_true_global")) df[[v]] <- to_num(df[[v]])
      if ("NonNormal" %in% names(df)) df$NonNormal <- as.character(df$NonNormal)
      if ("Outliers"  %in% names(df)) df$Outliers  <- as.character(df$Outliers)
      tidy_method(df)
    }
    s <- dplyr::bind_rows(lapply(sim_files, one))
    if (all(c("Structure", "NonNormal", "Outliers") %in% names(s)))
      s <- s[s$Structure == "Block" & toupper(s$NonNormal) == "FALSE" & toupper(s$Outliers) == "FALSE", ]
    s <- s[s$Method %in% c("EFA", "E-ReMI"), ]
    if (nrow(s) == 0) stop("no EFA/E-ReMI rows after filtering")
    add_labels(s)
  }, error = function(e) {
    message("Main-simulation files could not be used (", conditionMessage(e),
            "): Figure 3 and table C are skipped."); NULL })
} else {
  message("No sim_results files found: Figure 3 and table C are skipped.")
}

# ---------- 4. Numerical summaries (printed) -----------------------------------
wide <- tidyr::pivot_wider(
  droplevels(loio[, c("Sim", key_cols, "p_lab", "k_lab", "n_lab", "Method", "LOIO_R_squared")]),
  names_from = "Method", values_from = "LOIO_R_squared")
if (!all(c("EFA", "E-ReMI") %in% names(wide)))
  stop("Need LOIO values for BOTH methods (EFA and E-ReMI) to compare them.", call. = FALSE)
wide <- as.data.frame(wide)
wide$delta     <- wide[["E-ReMI"]] - wide[["EFA"]]
wide           <- wide[!is.na(wide$delta), ]
wide$identical <- abs(wide$delta) < 1e-12

cat("\n=== (A) How often are the two methods' LOIO R2 identical? ===\n")
tabA <- wide %>%
  dplyr::group_by(SampleSize, Items, TrueClusters) %>%
  dplyr::summarise(n_rep = dplyr::n(), n_identical = sum(identical),
                   pct_identical = 100 * mean(identical),
                   mean_delta_when_different = if (any(!identical)) mean(delta[!identical]) else NA_real_,
                   .groups = "drop")
if (!is.null(sim)) {
  ari <- sim[sim$Method == "E-ReMI", ] %>%
    dplyr::group_by(SampleSize, Items, TrueClusters) %>%
    dplyr::summarise(pct_EReMI_ARI_below_1_in_main_sim = 100 * mean(ARI_true < 1, na.rm = TRUE),
                     .groups = "drop")
  tabA <- dplyr::left_join(tabA, ari, by = key_cols)
}
print(as.data.frame(tabA), digits = 3)
cat("\nTotal identical pairs:", sum(wide$identical), "of", nrow(wide),
    "| E-ReMI lower in", sum(wide$delta < 0), ", higher in", sum(wide$delta > 0), "replications\n")

cat("\n=== (B) Mean and SD of LOIO R2 per cell and method ===\n")
tabB <- loio %>%
  dplyr::group_by(SampleSize, Items, TrueClusters, Method) %>%
  dplyr::summarise(mean = mean(LOIO_R_squared), sd = stats::sd(LOIO_R_squared), .groups = "drop")
print(as.data.frame(tabB), digits = 3)

if (!is.null(sim)) {
  cat("\n=== (C) Raw R2 columns from the main simulation vs LOIO (cell means) ===\n")
  simsum <- sim %>%
    dplyr::group_by(SampleSize, Items, TrueClusters, Method) %>%
    dplyr::summarise(raw = mean(R_squared_true, na.rm = TRUE),
                     item_avg = mean(R_squared_true_global, na.rm = TRUE), .groups = "drop")
  tabC <- as.data.frame(tidyr::pivot_wider(simsum, names_from = "Method", values_from = c("raw", "item_avg")))
  loio_mean <- loio %>% dplyr::group_by(SampleSize, Items, TrueClusters) %>%
    dplyr::summarise(LOIO_mean = mean(LOIO_R_squared), .groups = "drop")
  tabC <- dplyr::left_join(tabC, loio_mean, by = key_cols)
  print(tabC[, intersect(c(key_cols, "raw_EFA", "raw_E-ReMI", "item_avg_E-ReMI", "LOIO_mean"), names(tabC))],
        digits = 3)
}

# ---------- 5. Figures ----------------------------------------------------------
theme_loio <- ggplot2::theme_light(base_size = 11) +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey85", colour = NA),
    strip.text = ggplot2::element_text(colour = "black",face = "bold.italic"),
    legend.position = "top", panel.grid.minor = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(0.5, "lines"),
    plot.caption = ggplot2::element_text(size = 8, hjust = 0),
    
    # Only modify x- and y-axis tick labels
    axis.text.x = ggplot2::element_text(size = 12, face = "bold"),
    axis.text.y = ggplot2::element_text(size = 12, face = "bold")
  )

save_fig <- function(p, name, w = 11, h = 4.8) {
  png_file <- file.path(fig_dir, paste0(name, ".png"))
  ok <- FALSE
  if (requireNamespace("ragg", quietly = TRUE))
    ok <- tryCatch({ ggplot2::ggsave(png_file, p, width = w, height = h, dpi = 300, device = ragg::agg_png); TRUE },
                   error = function(e) FALSE)
  if (!ok)
    ok <- tryCatch({ ggplot2::ggsave(png_file, p, width = w, height = h, dpi = 300); TRUE },
                   error = function(e) { message("PNG export failed for ", name, ": ", conditionMessage(e)); FALSE })
  pdf_file <- file.path(fig_dir, paste0(name, ".pdf"))
  tryCatch(ggplot2::ggsave(pdf_file, p, width = w, height = h, device = grDevices::cairo_pdf),
           error = function(e) tryCatch(ggplot2::ggsave(pdf_file, p, width = w, height = h, device = "pdf"),
                                        error = function(e2) message("PDF export failed for ", name)))
  if (interactive()) print(p)
  invisible(ok)
}




# Figure 1 -- LOIO R2 by method (style of manuscript Figure 3) -------------------
fig1 <- ggplot2::ggplot(loio, ggplot2::aes(x = n_lab, y = LOIO_R_squared, fill = Method, colour = Method)) +
  ggplot2::geom_boxplot(width = 0.6, alpha = 0.75, outlier.shape = NA, colour = "grey25", linewidth = 0.3,
                        position = ggplot2::position_dodge(0.75)) +
  ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
                      size = 0.5, alpha = 0.35, show.legend = FALSE) +
  ggplot2::facet_grid(. ~ p_lab + k_lab) +
  ggplot2::scale_fill_manual(values = method_cols, labels = method_labs) +
  ggplot2::scale_colour_manual(values = method_cols, labels = method_labs) +
  ggplot2::labs(x = "Sample size", y = expression("LOIO " * R^2)) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = 12, face = "bold.italic"
    )
  ) +
  theme_loio
save_fig(fig1, "fig_loio_r2_by_method")

# Figure 2 -- paired differences (the key diagnostic) ----------------------------
yhi  <- max(0.016, max(wide$delta) * 1.15)
ylo  <- min(-0.02, min(wide$delta) * 1.15)
pct  <- wide %>% dplyr::group_by(n_lab, p_lab, k_lab) %>%
  dplyr::summarise(lab = paste0(round(100 * mean(identical)), "%"), .groups = "drop")
fig2 <- ggplot2::ggplot(wide, ggplot2::aes(x = n_lab, y = delta)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  ggplot2::geom_jitter(ggplot2::aes(colour = !identical), width = 0.15, height = 0, size = 1.1, alpha = 0.7) +
  ggplot2::geom_text(data = pct, ggplot2::aes(x = n_lab, y = 0.72 * yhi, label = lab),
                     size = 2.8, colour = "grey20", inherit.aes = FALSE) +
  ggplot2::facet_grid(. ~ p_lab + k_lab) +
  ggplot2::scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#D95F02"), guide = "none") +
  ggplot2::coord_cartesian(ylim = c(ylo, yhi)) +
  ggplot2::labs(x = "Sample size", y = paste0("\u0394 LOIO R\u00b2  (E-ReMI \u2212 ", efa_label, ")")) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      size = 12, face = "bold.italic"
    )
  )+
  theme_loio
save_fig(fig2, "fig_loio_r2_paired_difference")

# Figure 3 -- raw R2 columns of the main simulation next to LOIO R2 ---------------
if (!is.null(sim)) {
  ser <- c(paste0(efa_label, " raw R\u00b2 (mean communality)"),
           "E-ReMI raw R\u00b2 as computed in the main script (mean over clusters of the best item)",
           "E-ReMI item-averaged R\u00b2 (each item vs. its best-matching cluster mean; column R_squared_true_global)",
           paste0("LOIO R\u00b2 (", efa_label, " and E-ReMI)"))
  cell_mean <- function(df, col, method = NULL, series) {
    if (!is.null(method)) df <- df[df$Method == method, ]
    out <- stats::aggregate(df[[col]], by = list(n_lab = df$n_lab, p_lab = df$p_lab, k_lab = df$k_lab),
                            FUN = mean, na.rm = TRUE)
    names(out)[names(out) == "x"] <- "y"; out$series <- series; out
  }
  d3 <- dplyr::bind_rows(
    cell_mean(sim,  "R_squared_true",        "EFA",    ser[1]),
    cell_mean(sim,  "R_squared_true",        "E-ReMI", ser[2]),
    cell_mean(sim,  "R_squared_true_global", "E-ReMI", ser[3]),
    cell_mean(loio, "LOIO_R_squared",        NULL,     ser[4]))
  cond_levels <- unlist(lapply(p_levels, function(p) paste0("p = ", p, ", k = ", k_levels)))
  d3$series <- factor(d3$series, levels = ser)
  d3$cond   <- factor(paste0(d3$p_lab, ", ", d3$k_lab), levels = cond_levels)
  fig3 <- ggplot2::ggplot(d3, ggplot2::aes(x = cond, y = y, colour = series, shape = series, fill = series)) +
    ggplot2::geom_point(size = 2.8, stroke = 0.9, position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::facet_grid(. ~ n_lab) +
    ggplot2::scale_colour_manual(values = c("#1F78B4", "#D95F02", "#D95F02", "black")) +
    ggplot2::scale_shape_manual(values = c(16, 16, 24, 15)) +
    ggplot2::scale_fill_manual(values = c("#1F78B4", "#D95F02", "white", "black")) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1), shape = ggplot2::guide_legend(ncol = 1),
                    fill = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::labs(x = NULL, y = expression("Mean "*R^2*" per cell"), colour = NULL, shape = NULL, fill = NULL) +
    theme_loio + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
                                legend.text = ggplot2::element_text(size = 8), legend.position = "bottom")
  save_fig(fig3, "fig_raw_vs_loio_decomposition", h = 5.8)
}

message("Done. Figures written to: ", fig_dir)