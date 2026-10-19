# ============================================================================
# New simulation + pilot runner (single file).
# EDIT `path` and `output_dir` below to match your machine, then source/run.
# ============================================================================

# ---------- 0. Setup --------------------------------------------------------
# Load Packages
#install.packages(c("future", "future.apply", "progressr"))
library(future)
library(future.apply)
library(progressr)
library(psych)
library(cluster) 
library(mclust) 
library(MASS)
library(ggplot2) 
library(dplyr) 
library(tictoc)
library(optimbase)          ##### This packages will to create a random data matrix 
library(pracma)
library(pracma)             ##### This packages will to create a random data matrix 
library(reshape)            ##### This package will use to take transpose
library(optimbase)
library(gdata)
library(Rfast)
library(R.utils)            ##### This package is for printf function
library(dplyr)
library(janitor)
library(fossil)
library(lintools)
library(pdfCluster)

set.seed(12345)
RNGkind("L'Ecuyer-CMRG")

# Load source functions
path <- '~/Documents/BOOTEFA-Server-number cruncher/'
source(paste(path,'Randompartition_function.R',sep = ''))
source(paste(path,'Unequal_Randompartition_Function.R',sep = ''))
source(paste(path,'Log_Likelihood_function_REMAXINT.R',sep = ''))
source(paste(path,'Log_Likelihood_function_E_ReMI.R',sep = ''))
source(paste(path,'Log_LR_Test_Statistic.R',sep = ''))
source(paste(path,'Update_column_clusters.R',sep = ''))
source(paste(path,'Update_row_clusters_REMAXINT.R',sep = ''))
source(paste(path,'Update_row_clusters_E_ReMI.R',sep = ''))
source(paste(path,'Update_G_Omega.R',sep = ''))
source(paste(path,'REMAXINT.R',sep = ''))
source(paste(path,'E_ReMI.R',sep = ''))
source(paste(path,'Permutation_Function.R',sep = ''))
source(paste(path,'Data_Gen_Function.R',sep = ''))

# Files each PARALLEL WORKER sources so the E-ReMI stack is guaranteed present
# in every worker process (independent of future's auto-export). The simulation
# functions defined below in THIS file reach workers via future's auto-export,
# so this list must NOT include this combined script.
worker_source_files <- c(
  paste(path,'Randompartition_function.R',sep = ''),
  paste(path,'Unequal_Randompartition_Function.R',sep = ''),
  paste(path,'Log_Likelihood_function_REMAXINT.R',sep = ''),
  paste(path,'Log_Likelihood_function_E_ReMI.R',sep = ''),
  paste(path,'Log_LR_Test_Statistic.R',sep = ''),
  paste(path,'Update_column_clusters.R',sep = ''),
  paste(path,'Update_row_clusters_REMAXINT.R',sep = ''),
  paste(path,'Update_row_clusters_E_ReMI.R',sep = ''),
  paste(path,'Update_G_Omega.R',sep = ''),
  paste(path,'REMAXINT.R',sep = ''),
  paste(path,'E_ReMI.R',sep = ''),
  paste(path,'Permutation_Function.R',sep = ''),
  paste(path,'Data_Gen_Function.R',sep = '')
)
worker_packages <- c("psych","MASS","mclust","pracma","optimbase","Rfast","reshape")

# Output directory
output_dir <- "~/Documents/BOOTEFA-Server-number cruncher"

# ============================================================================
# Parallelized simulation script for the manuscript revision (REDUCED DESIGN).
# ============================================================================
# =============================================================================
# 1. Generate True loading and label matrices (Block structure only is used,
#    but the Diagonal branch is harmless so we keep it.
# =============================================================================
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

# =============================================================================
# 2. Replication-specific R²_max (recovery ceiling). Depends only on Lambda
#    and psi; identical denominator for all methods within a cell.
# =============================================================================
compute_R_squared_max <- function(true_Lambda, psi_value = 0.3) {
  sum_lambda_sq_per_item <- rowSums(true_Lambda^2)
  mean(sum_lambda_sq_per_item / (sum_lambda_sq_per_item + psi_value))
}

# =============================================================================
# 3. Simulate Likert data. Kept general (non_normal/outliers flags retained)
#    but the design always calls it with FALSE/FALSE now. They never impact 
#    R suqared ARI and Structural recovery count too much. 
# =============================================================================
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
      shift_amounts <- sample(1:3, size = p, replace = TRUE,
                              prob = c(0.4, 0.4, 0.2))
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
      extreme_value <- sample(c(1, 5), size = 1)
      likert_data[i, ] <- extreme_value
    }
  }
  
  likert_data <- pmax(pmin(likert_data, 5), 1)
  return(likert_data)
}

# =============================================================================
# 4. Structure-count helpers (unchanged; independent of k2)
# =============================================================================
select_k_efa <- function(likert_data, n_iter = 50) {
  pa <- tryCatch(suppressWarnings(suppressMessages(
    psych::fa.parallel(likert_data, fa = "fa", fm = "minres",
                       cor = "poly", plot = FALSE, n.iter = n_iter))),
    error = function(e) NULL)
  if (is.null(pa) || is.null(pa$nfact)) return(NA_integer_)
  as.integer(pa$nfact)
}

# ============================================================================
# select_k_ereMI_permnull(): choose the number of column clusters Q for
# E-ReMI by an INCREMENTAL permutation-null gain rule.
#
# WHY THIS IS VALID EVEN THOUGH THE NULL HAS NO CLUSTERING
# -------------------------------------------------------------------------
# The permutation null (sample() over all cells, inside Permutation_Function)
# is the no-interaction null: H0 has no clustering. Our data, by construction,
# HAS column clustering (H1). A per-Q significance test (observed LLR vs its
# critical value) is therefore useless for COUNTING clusters: structured data
# beats the no-structure null at almost every Q -- under-, correct- and
# over-extracted alike. It only tells us "is there ANY structure" (used as a
# floor, nothing more).
#
# The QUANTITY that does count clusters is the INCREMENT. The null increment
# = how much LLR grows from Q-1 to Q on no-structure data (pure chance gain
# from one extra column cluster). The observed increment = chance + any REAL
# gain. For Q <= K_true each added cluster captures a real item block, so the
# observed increment sits far above the null increment. For Q > K_true the
# extra cluster only splits a real block into identically-loaded pieces, adds
# no interaction, and the observed increment falls back to the null level.
# Rule: keep adding columns while the observed increment beats the null
# increment; stop when it doesn't. The no-structure null is used ONLY to
# calibrate "what a chance-only increment looks like" -- a question it CAN
# answer.
#
# LIMITATION (report it as such): the null increment is built under TOTAL
# absence of structure, not conditional on "Q-1 real clusters present". It is
# an approximate, slightly conservative (toward fewer clusters) incremental
# test, not an exact conditional permutation test.
#
# USE WITH FIXED P (P_equals_Q = FALSE) so each grid step changes ONLY Q.
# Set the driver's P large enough for rank (e.g. P = 6 >= largest k1).
#
# REQUIRES: Permutation_Function in scope.
# ============================================================================

# ---- Robustly pull the E-ReMI permutation-null vector ----------------------
# Permutation_Function names this element "Perm_Log_LR_EReMI " WITH a trailing
# space, so out$Perm_Log_LR_EReMI is NULL. Match on trimmed names so this
# works whether or not the source is fixed.
get_perm_null_EReMI <- function(out) {
  if (!is.list(out)) return(NULL)
  nm <- names(out)
  if (is.null(nm)) return(NULL)
  hit <- which(trimws(nm) == "Perm_Log_LR_EReMI")
  if (length(hit) == 0) return(NULL)
  v <- suppressWarnings(as.numeric(out[[hit[1]]]))
  v[is.finite(v)]
}

# ---- Upper (1 - alpha) quantile of the null gain (LLR_hi - LLR_lo) ----------
# Built from the two independent per-Q permutation nulls. Exact pairwise
# differences when small (deterministic); sampled pairs otherwise.
null_gain_threshold <- function(perm_lo, perm_hi, alpha, n_pairs = 4000) {
  if (is.null(perm_lo) || is.null(perm_hi)) return(NA_real_)
  perm_lo <- perm_lo[is.finite(perm_lo)]; perm_hi <- perm_hi[is.finite(perm_hi)]
  if (length(perm_lo) < 2L || length(perm_hi) < 2L) return(NA_real_)
  alpha <- max(min(alpha, 0.999), 1e-6)
  diffs <- if (as.numeric(length(perm_lo)) * length(perm_hi) <= 1e6) {
    as.numeric(outer(perm_hi, perm_lo, "-"))
  } else {
    sample(perm_hi, n_pairs, TRUE) - sample(perm_lo, n_pairs, TRUE)
  }
  diffs <- diffs[is.finite(diffs)]
  if (!length(diffs)) return(NA_real_)
  as.numeric(stats::quantile(diffs, 1 - alpha, names = FALSE, na.rm = TRUE))
}

# ---- Main selector ---------------------------------------------------------
select_k_ereMI_permnull <- function(likert_data, P, k1_target, Nruns,
                                    permutations_grid, alpha_level,
                                    Q_grid = NULL,
                                    P_equals_Q = FALSE,   # keep FALSE: fixed P
                                    n_pairs = 4000,
                                    return_details = FALSE) {
  
  fail <- function() {
    if (isTRUE(return_details))
      list(selected = NA_integer_, table = NULL, any_structure = NA)
    else NA_integer_
  }
  
  # --- input guards -------------------------------------------------------
  if (missing(likert_data) || is.null(likert_data)) return(fail())
  D <- tryCatch(as.matrix(likert_data), error = function(e) NULL)
  if (is.null(D) || !is.numeric(D) || nrow(D) < 2L || ncol(D) < 2L) return(fail())
  if (!is.numeric(P) || length(P) != 1L || P < 2L) return(fail())
  if (!is.numeric(k1_target) || length(k1_target) != 1L) return(fail())
  
  if (is.null(Q_grid)) {
    r <- ceiling(k1_target / 2)
    Q_grid <- seq(max(2L, k1_target - r), k1_target + r)
  }
  Q_grid <- sort(unique(as.integer(Q_grid)))
  Q_grid <- Q_grid[is.finite(Q_grid) & Q_grid >= 2L]
  m <- length(Q_grid)
  if (m < 1L) return(fail())
  
  obs  <- rep(NA_real_, m)     # observed E-ReMI LLR per Q
  crit <- rep(NA_real_, m)     # permutation critical value per Q
  perm <- vector("list", m)    # permutation null LLR vector per Q
  
  for (i in seq_len(m)) {
    Q     <- Q_grid[i]
    P_use <- if (isTRUE(P_equals_Q)) Q else P
    if (P_use < 2L || Q < 2L) next
    out <- tryCatch(
      Permutation_Function(D, P_use, Q, Nruns, permutations_grid, alpha_level),
      error = function(e) NULL)
    if (is.null(out) || !is.list(out)) next
    o <- suppressWarnings(as.numeric(out$Obs_Log_LR_EReMI)[1])
    cc <- suppressWarnings(as.numeric(out$Crit_Value_Perm_EReMI)[1])
    if (length(o))  obs[i]  <- o
    if (length(cc)) crit[i] <- cc
    pv <- get_perm_null_EReMI(out)
    perm[[i]] <- if (is.null(pv)) numeric(0) else pv   # NULL would delete the slot
  }
  
  valid <- which(is.finite(obs))
  if (!length(valid)) return(fail())
  
  # --- floor: is there ANY structure at all? -----------------------------
  any_structure <- any(is.finite(obs) & is.finite(crit) & obs > crit)
  
  # --- per-step null-calibrated gain significance ------------------------
  g_obs <- rep(NA_real_, m); g_thr <- rep(NA_real_, m); sig <- rep(NA, m)
  if (m >= 2L) {
    for (k in 2:m) {
      if (!is.finite(obs[k]) || !is.finite(obs[k - 1])) next
      g_obs[k] <- obs[k] - obs[k - 1]
      g_thr[k] <- null_gain_threshold(perm[[k - 1]], perm[[k]],
                                      alpha_level, n_pairs)
      sig[k]   <- is.finite(g_thr[k]) && is.finite(g_obs[k]) &&
        (g_obs[k] > g_thr[k])
    }
  }
  
  # --- walk up from the smallest valid Q while gains beat the null --------
  selected <- Q_grid[valid[1]]
  if (m >= 2L) {
    for (k in 2:m) {
      if (isTRUE(sig[k])) selected <- Q_grid[k] else break
    }
  }
  if (!isTRUE(any_structure)) selected <- NA_integer_
  
  if (isTRUE(return_details)) {
    tab <- data.frame(
      Q = Q_grid,
      P = if (isTRUE(P_equals_Q)) Q_grid else rep(P, m),
      Obs_LLR = obs, Crit = crit,
      gain_obs = g_obs, gain_null_thr = g_thr, gain_sig = sig,
      stringsAsFactors = FALSE)
    return(list(selected = as.integer(selected), table = tab,
                any_structure = any_structure))
  }
  as.integer(selected)
}
categorize_k <- function(selected_k, k1) {
  if (is.na(selected_k)) return(NA_character_)
  if (selected_k == k1) return("Correct")
  if (selected_k <  k1) return("Under")
  if (selected_k >  k1) return("Over")
  NA_character_
}

# =============================================================================
# 5. Method evaluation functions (TRUE / k1 model only)
# =============================================================================
#-----------------------------------------------------
# 5.1 EFA [polychoric, no bootstrap]
#-----------------------------------------------------
evaluate_efa <- function(likert_data, p, k1, true_labels) {
  poly_cor <- psych::polychoric(likert_data)$rho
  ft <- psych::fa(r = poly_cor, nfactors = k1, n.obs = nrow(likert_data),
                  rotate = "oblimin", fm = "minres", cor = "poly")
  R1 <- mean(ft$communality)
  R3 <- cor(as.vector(poly_cor), as.vector(ft$model))^2
  L1 <- apply(as.matrix(ft$loadings), 1, function(x) which.max(abs(x)))
  A1 <- tryCatch(mclust::adjustedRandIndex(true_labels, L1), error = function(e) NA)
  list(R_squared_true = R1,
       ARI_true = A1,
       R_squared_true_global = R3)
}

#--------------------------------------------------------------------
# 5.2 EFA [polychoric, nonparametric bootstrap]
#--------------------------------------------------------------------
evaluate_efa_boot <- function(likert_data, p, k1, true_labels, B = 1000) {
  n <- nrow(likert_data)
  v <- function() rep(NA_real_, B)
  R2t <- v(); R2tg <- v(); At <- v()
  for (b in 1:B) {
    idx <- sample.int(n, n, replace = TRUE)
    bd  <- likert_data[idx, , drop = FALSE]
    o <- tryCatch({
      pc <- suppressWarnings(psych::polychoric(bd)$rho)
      ft <- suppressWarnings(psych::fa(r = pc, nfactors = k1, n.obs = n,
                                       rotate = "oblimin", fm = "minres", cor = "poly"))
      Lt <- apply(as.matrix(ft$loadings), 1, function(x) which.max(abs(x)))
      list(R2t  = mean(ft$communality),
           R2tg = cor(as.vector(pc), as.vector(ft$model))^2,
           At   = mclust::adjustedRandIndex(true_labels, Lt))
    }, error = function(e) list(R2t = NA, R2tg = NA, At = NA))
    R2t[b] <- o$R2t; R2tg[b] <- o$R2tg; At[b] <- o$At
  }
  ag <- function(x) list(mean = mean(x, na.rm = TRUE),
                         lo = unname(quantile(x, 0.025, na.rm = TRUE)),
                         hi = unname(quantile(x, 0.975, na.rm = TRUE)),
                         sd = sd(x, na.rm = TRUE))
  aR2t <- ag(R2t); aR2tg <- ag(R2tg); aAt <- ag(At)
  list(R_squared_true = aR2t$mean,
       R_squared_true_global = aR2tg$mean,
       ARI_true = aAt$mean,
       R_squared_true_lo = aR2t$lo, R_squared_true_hi = aR2t$hi,
       ARI_true_lo = aAt$lo, ARI_true_hi = aAt$hi,
       R_squared_true_sd = aR2t$sd, ARI_true_sd = aAt$sd,
       B_used = B, n_valid = sum(!is.na(R2t)))
}

#---------------------------------
# 5.3 E-ReMI
#---------------------------------
evaluate_E_ReMI <- function(likert_data, P, k1, Nruns, permutations,
                            alpha_level, true_labels) {
  fit_one <- function(Q) {
    out <- Permutation_Function(likert_data, P, Q, Nruns, permutations, alpha_level)
    K <- out$K_EReMI
    if (is.null(K) || !is.matrix(K)) return(NULL)
    Data <- as.data.frame(likert_data); p <- ncol(Data)
    # Force list output and guard empty clusters (bulletproofing).
    cols_with_1s <- lapply(seq_len(ncol(K)), function(j) which(K[, j] == 1))
    valid <- which(vapply(cols_with_1s, length, integer(1)) > 0)
    cols_with_1s <- cols_with_1s[valid]
    if (length(cols_with_1s) == 0) return(NULL)
    for (i in seq_along(cols_with_1s)) {
      sel <- Data[, colnames(Data)[cols_with_1s[[i]]], drop = FALSE]
      Data[[paste0("Cluster", i)]] <- rowMeans(sel, na.rm = TRUE)
    }
    ic <- colnames(Data)[1:p]; cc <- paste0("Cluster", seq_along(cols_with_1s))
    cor_mat <- cor(Data[, c(ic, cc)], method = "pearson")
    r2 <- cor_mat[ic, cc, drop = FALSE]^2
    list(R  = mean(apply(r2, 2, max)),
         Rg = mean(apply(r2, 1, max)),
         A  = tryCatch(mclust::adjustedRandIndex(
           true_labels, apply(K, 1, function(x) which.max(abs(x)))),
           error = function(e) NA))
  }
  t1 <- fit_one(k1); if (is.null(t1)) return(NULL)
  list(R_squared_true = t1$R,
       ARI_true = t1$A,
       R_squared_true_global = t1$Rg)
}

# =============================================================================
# 6. Row-builder. Common schema for all methods; bootstrap-only columns are
#    NA for the non-bootstrap methods. No more "_mis" / MisClusters columns.
# =============================================================================
make_row <- function(method, out, R_squared_max, selected_k, k1,
                     sim, n, p, cond, boot_extras = NULL) {
  if (is.null(boot_extras)) {
    boot_extras <- list(R_squared_true_lo = NA_real_, R_squared_true_hi = NA_real_,
                        ARI_true_lo = NA_real_, ARI_true_hi = NA_real_,
                        R_squared_true_sd = NA_real_, ARI_true_sd = NA_real_,
                        B_used = NA_integer_, n_valid = NA_integer_)
  }
  R1 <- out$R_squared_true
  data.frame(
    Sim = sim, SampleSize = n, Items = p,
    TrueClusters = k1,
    Structure = cond$Structure, NonNormal = cond$NonNormal, Outliers = cond$Outliers,
    Method = method,
    R_squared_true = R1,
    R_squared_true_global = out$R_squared_true_global,
    ARI_true = out$ARI_true,
    R_squared_max = R_squared_max,
    R_squared_true_ratio = R1 / R_squared_max,
    Selected_K = selected_k,
    K_Status = categorize_k(selected_k, k1),
    Under_Factoring   = as.integer(!is.na(selected_k) & selected_k <  k1),
    Correct_Factoring = as.integer(!is.na(selected_k) & selected_k == k1),
    Over_Factoring    = as.integer(!is.na(selected_k) & selected_k >  k1),
    R_squared_true_lo = boot_extras$R_squared_true_lo,
    R_squared_true_hi = boot_extras$R_squared_true_hi,
    ARI_true_lo = boot_extras$ARI_true_lo, ARI_true_hi = boot_extras$ARI_true_hi,
    R_squared_true_sd = boot_extras$R_squared_true_sd,
    ARI_true_sd = boot_extras$ARI_true_sd,
    B_used = boot_extras$B_used, B_valid = boot_extras$n_valid,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 7. Per-cell worker (no k2)
# =============================================================================
process_one_cell <- function(sim, p, k1, cond_idx, conditions,
                             n, P, Nruns, permutations, permutations_grid,
                             alpha_level, B_boot, psi_value, pa_n_iter) {
  cond <- conditions[cond_idx, ]
  
  obj <- generate_true_structure(p, k1, cond$Structure)
  true_Lambda <- obj$true_Lambda; true_labels <- obj$true_labels
  R_squared_max <- compute_R_squared_max(true_Lambda, psi_value)
  data <- simulate_data(n, true_Lambda, cond$NonNormal, cond$Outliers, psi_value)
  
  selected_k_efa   <- select_k_efa(data, n_iter = pa_n_iter)
  
  # FIX: assign the result (was discarded) and pass k1 (k1_target did not exist
  # in this scope). Without the assignment, selected_k_ereMI was undefined below
  # and EVERY cell errored -> all cells dropped -> empty output.
  selected_k_ereMI <- select_k_ereMI_permnull(
    data, P, k1, Nruns,
    permutations_grid, alpha_level,
    Q_grid = NULL,
    P_equals_Q = FALSE,   # keep FALSE: fixed P
    n_pairs = 4000,
    return_details = FALSE)
  
  selected_k_boot  <- selected_k_efa  # bootstrap doesn't change PA selection
  
  rows <- list()
  
  # EFA
  out_efa <- evaluate_efa(data, p, k1, true_labels)
  rows[[1]] <- make_row("EFA", out_efa, R_squared_max, selected_k_efa, k1,
                        sim, n, p, cond)
  
  # E-ReMI
  out_ereMI <- evaluate_E_ReMI(data, P, k1, Nruns, permutations,
                               alpha_level, true_labels)
  if (!is.null(out_ereMI)) {
    rows[[length(rows) + 1]] <- make_row("E-ReMI", out_ereMI, R_squared_max,
                                         selected_k_ereMI, k1,
                                         sim, n, p, cond)
  }
  
  # EFA-boot
  out_boot <- evaluate_efa_boot(data, p, k1, true_labels, B = B_boot)
  boot_extras <- out_boot[c("R_squared_true_lo","R_squared_true_hi",
                            "ARI_true_lo","ARI_true_hi",
                            "R_squared_true_sd","ARI_true_sd",
                            "B_used","n_valid")]
  rows[[length(rows) + 1]] <- make_row("EFA-boot", out_boot, R_squared_max,
                                       selected_k_boot, k1,
                                       sim, n, p, cond,
                                       boot_extras = boot_extras)
  
  do.call(rbind, rows)
}

# =============================================================================
# 8. Main driver — parallel over cells via future.apply
# =============================================================================
run_simulation_for_single_n_parallel <- function(n, session_id, Nruns, n_sim,
                                                 B_boot = 500,
                                                 permutations = 500,
                                                 permutations_grid = 50,
                                                 ncores = 120,
                                                 output_dir = "C:/Users/Zaheer/Desktop/BOOTEFA-Server-number cruncher",
                                                 worker_source_files = NULL,
                                                 worker_packages = NULL) {
  # ---- Design grid (REDUCED) ----
  item_levels     <- c(30, 60)
  cluster_levels1 <- c(4, 6)
  P               <- 6     # FIXED, >= largest k1. The interaction core Gamma is
  # P x Q (rank <= min(P,Q)); P must be >= the column
  # structure for E-ReMI to express it. Do NOT couple
  # P to Q. Keep this constant across the whole grid.
  alpha_level     <- 0.05
  psi_value       <- 0.3
  pa_n_iter       <- 50
  
  # Single condition: Block structure, normal, no outliers.
  conditions <- expand.grid(
    Structure = c("Block"),
    NonNormal = c(FALSE),
    Outliers  = c(FALSE),
    stringsAsFactors = FALSE
  )
  
  all_cells <- expand.grid(
    sim      = 1:n_sim,
    p        = item_levels,
    k1       = cluster_levels1,
    cond_idx = 1:nrow(conditions)
  )
  message(sprintf("Total cells to process: %d  (using %d cores)",
                  nrow(all_cells), ncores))
  
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Install future.apply: install.packages('future.apply')")
  }
  future::plan(future::multisession, workers = ncores)
  on.exit(future::plan(future::sequential), add = TRUE)  # always reset
  
  has_progressr <- requireNamespace("progressr", quietly = TRUE)
  if (has_progressr) progressr::handlers("txtprogressbar")
  
  t0 <- Sys.time()
  message("Started at: ", format(t0))
  
  # Packages every worker must load. Includes E-ReMI source dependencies
  # (pracma supplies ones(), etc.) on top of the EFA/cluster stack.
  pkgs_for_workers <- unique(c("psych", "MASS", "mclust",
                               "pracma", "optimbase", "Rfast", "reshape",
                               worker_packages))
  
  run_block <- function() {
    if (has_progressr) {
      p_handler <- progressr::progressor(steps = nrow(all_cells))
    }
    future.apply::future_lapply(
      seq_len(nrow(all_cells)),
      future.seed = TRUE,
      future.packages = pkgs_for_workers,
      FUN = function(i) {
        # ---- one-time per-worker setup -----------------------------------
        # Source the E-ReMI (+ optionally the simulation) files INSIDE each
        # worker so the functions are guaranteed present, independent of
        # future's global auto-detection. Runs once per worker process.
        if (!isTRUE(getOption("sim_worker_ready", FALSE))) {
          if (length(worker_source_files)) {
            suppressWarnings(suppressMessages(
              for (f in worker_source_files) if (file.exists(f)) try(source(f), silent = TRUE))
            )
          }
          options(sim_worker_ready = TRUE)
        }
        row <- all_cells[i, ]
        # Bulletproof: a single failed cell returns NULL instead of aborting
        # the whole gather. Dropped cells are counted afterwards.
        out <- tryCatch(
          process_one_cell(
            sim = row$sim, p = row$p, k1 = row$k1,
            cond_idx = row$cond_idx, conditions = conditions,
            n = n, P = P, Nruns = Nruns,
            permutations = permutations,
            permutations_grid = permutations_grid,
            alpha_level = alpha_level, B_boot = B_boot,
            psi_value = psi_value, pa_n_iter = pa_n_iter
          ),
          error = function(e) NULL
        )
        if (has_progressr) p_handler()
        out
      }
    )
  }
  
  results_list <- if (has_progressr) {
    progressr::with_progress(run_block())
  } else {
    run_block()
  }
  
  n_failed <- sum(vapply(results_list, is.null, logical(1)))
  results_list <- Filter(Negate(is.null), results_list)
  if (n_failed > 0) {
    message(sprintf("WARNING: %d of %d cells failed and were dropped.",
                    n_failed, nrow(all_cells)))
  }
  
  # Do not silently write a blank file if everything failed. Re-run ONE cell
  # locally (sequentially) to surface the real error message for debugging.
  if (length(results_list) == 0) {
    diag_msg <- tryCatch({
      process_one_cell(
        sim = all_cells$sim[1], p = all_cells$p[1], k1 = all_cells$k1[1],
        cond_idx = all_cells$cond_idx[1], conditions = conditions,
        n = n, P = P, Nruns = Nruns,
        permutations = permutations, permutations_grid = permutations_grid,
        alpha_level = alpha_level, B_boot = B_boot,
        psi_value = psi_value, pa_n_iter = pa_n_iter)
      "cell ran when re-tried sequentially (the failure may be worker/export related)"
    }, error = function(e) conditionMessage(e))
    stop("All cells failed; nothing to write. First-cell error was:\n  ", diag_msg)
  }
  
  results <- do.call(rbind, results_list)
  
  t1 <- Sys.time()
  message("Finished at: ", format(t1))
  message("Elapsed: ", format(round(difftime(t1, t0, units = "hours"), 2)))
  
  base_name <- paste0("sim_results_n", n, "_session", session_id, "_reduced")
  write.csv(results,
            file = file.path(output_dir, paste0(base_name, ".csv")),
            row.names = FALSE)
  saveRDS(results,
          file = file.path(output_dir, paste0(base_name, ".rds")))
  
  invisible(results)
}


# ---------- 1. (Optional but recommended) Pilot run ------------------------
# Confirms the pipeline works end-to-end on a small workload (a few minutes
# to ~1 hour). DO this once before launching the multi-day production run.
RUN_PILOT <- TRUE
if (RUN_PILOT) {
  message("==== Running pilot (n=300, n_sim=100, 120 cores) ====")
  pilot <- run_simulation_for_single_n_parallel(
    n=500, session_id = "pilot",
    Nruns = 20, n_sim = 100,
    B_boot = 50, permutations = 50, permutations_grid = 50,
    ncores = 120, output_dir = output_dir,
    worker_source_files = worker_source_files,
    worker_packages = worker_packages
  )
  # Sanity checks on the pilot
  stopifnot(nrow(pilot) > 0)
  stopifnot(all(c("Sim","Method","R_squared_max","Selected_K","K_Status")
                %in% colnames(pilot)))
  # Confirm R²_max is identical across the 3 methods within each cell
  key_cols <- c("Sim","SampleSize","Items","TrueClusters","Structure")
  chk <- aggregate(R_squared_max ~ ., pilot[, c(key_cols, "R_squared_max")],
                   function(x) length(unique(x)))
  stopifnot(all(chk$R_squared_max == 1))
  message("==== Pilot passed all sanity checks ====")
}

# ---------- 2. Production runs across all sample sizes ---------------------
# Flip to TRUE to launch the multi-day production run AFTER the pilot passes.
# RUN_PRODUCTION <- TRUE
# 
# all_results <- list()   # always defined, so section 3 never errors
# 
# if (RUN_PRODUCTION) {
#   n_values <- c(100, 300, 500)
# 
#   for (n_idx in seq_along(n_values)) {
#     n_val <- n_values[n_idx]
# 
#     message("\n", strrep("=", 76))
#     message(sprintf(" Starting production run %d of %d:  n = %d",
#                     n_idx, length(n_values), n_val))
#     message(sprintf(" Started at: %s", format(Sys.time())))
#     message(strrep("=", 76))
# 
#     # tryCatch keeps the loop going even if one n value errors out
#     res_n <- tryCatch({
#       run_simulation_for_single_n_parallel(
#         n            = n_val,
#         session_id   = paste0("prod_n", n_val),
#         Nruns        = 20,
#         n_sim        = 500,
#         B_boot       = 500,
#         permutations = 100,
#         permutations_grid = 250,
#         ncores       = 120,
#         output_dir   = output_dir,
#         worker_source_files = worker_source_files,
#         worker_packages     = worker_packages
#       )
#     }, error = function(e) {
#       message(sprintf("ERROR for n = %d: %s", n_val, conditionMessage(e)))
#       NULL
#     })
# 
#     if (!is.null(res_n)) {
#       all_results[[as.character(n_val)]] <- res_n
#       combined_so_far <- do.call(rbind, all_results)
#       saveRDS(combined_so_far,
#               file = file.path(output_dir, "sim_results_combined_so_far.rds"))
#     }
#     message(sprintf(" Finished n = %d at: %s", n_val, format(Sys.time())))
#   }
# }
# 
# # ---------- 3. Combine all results into a single long-format file ----------
# if (exists("all_results") && length(all_results) > 0) {
#   combined <- do.call(rbind, all_results)
#   saveRDS(combined,
#           file = file.path(output_dir, "sim_results_all3methods_combined.rds"))
#   write.csv(combined,
#             file = file.path(output_dir, "sim_results_all3methods_combined.csv"),
#             row.names = FALSE)
#   message(sprintf("\nCombined file saved: %d rows across n = {%s}",
#                   nrow(combined), paste(names(all_results), collapse = ", ")))
# }
# 
# message("\nAll runs complete.")