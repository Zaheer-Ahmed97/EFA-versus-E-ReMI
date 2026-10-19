#################################################################################################
#  Case study (Section 6): Trier Inventory for Chronic Stress (TICS) and burnout (MBI)
#  "Jointly Clustering Respondents and Items in Likert Response Data"
#################################################################################################

## ============================================================================================
## 0. SETUP
## ============================================================================================
data_dir    <- "C:/Users/ahmed/Desktop/IfADo Projects/Third Project based on Survey and DVS Data/Folder for submitted paper/Case Studies Section"
data_source <- "raw"   # "raw": import raw files and impute here; "imputed": use TICS_Data.csv
k_impute    <- 10      # neighbours for k-NN imputation
n_reps_cv   <- 100     # fold assignments for the repeated cross-validation (Section 6.4.2)
P_eremi     <- 3       # E-ReMI participant clusters
seed        <- 123

out_dir <- "case_study_output"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
if (!dir.exists(data_dir)) stop("data_dir does not exist: ", data_dir, call. = FALSE)
stopifnot(data_source %in% c("raw", "imputed"), n_reps_cv >= 2, P_eremi >= 2)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

## ---- Packages ------------------------------------------------------------------------------
required_pkgs <- c("psych", "GPArotation", "glmnet", "ggplot2", "dplyr", "tidyr", "mclust", "ragg",
                   if (data_source == "raw") "VIM")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall with: install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))", call. = FALSE)
# GPArotation is essential: without it psych::fa(rotate = "oblimin") silently returns an unrotated solution.
has <- vapply(c("DiagrammeR", "DiagrammeRsvg", "rsvg", "magick"), requireNamespace, logical(1), quietly = TRUE)

suppressPackageStartupMessages({
  library(psych); library(glmnet); library(ggplot2); library(dplyr); library(tidyr); library(mclust)
})
set.seed(seed)

## ---- Helpers --------------------------------------------------------------------------------
log_file <- file.path(out_dir, "case_study_results.txt")
cat("TICS case study results -", format(Sys.time()), "\n\n", file = log_file)
report <- function(..., obj = NULL) {
  txt <- paste0(...)
  cat(txt, "\n"); cat(txt, "\n", file = log_file, append = TRUE)
  if (!is.null(obj)) {
    out <- capture.output(print(obj))
    cat(out, sep = "\n"); cat(out, sep = "\n", file = log_file, append = TRUE)
    cat("\n", file = log_file, append = TRUE)
  }
}
save_table <- function(df, name) write.csv(df, file.path(tab_dir, paste0(name, ".csv")), row.names = FALSE)

# Figure export: vector PDF + 600 dpi TIFF + 600 dpi PNG
save_ggplot_hires <- function(plot, filename, width_mm = 174, height_mm = 120, dpi = 600) {
  base <- file.path(fig_dir, filename)
  ggsave(paste0(base, ".pdf"), plot, device = grDevices::cairo_pdf, width = width_mm, height = height_mm,
         units = "mm", bg = "white")
  ggsave(paste0(base, ".tiff"), plot, device = ragg::agg_tiff, width = width_mm, height = height_mm,
         units = "mm", dpi = dpi, compression = "lzw", bg = "white")
  ggsave(paste0(base, ".png"), plot, device = ragg::agg_png, width = width_mm, height = height_mm,
         units = "mm", dpi = dpi, bg = "white")
  message("Saved: ", base, ".pdf / .tiff / .png")
}
theme_paper <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))

fmt_p <- function(p) ifelse(p < .001, "< .001", sub("^0", "", sprintf("%.3f", p)))
ari   <- mclust::adjustedRandIndex
map_to_vec <- function(map, p = 57) { v <- rep(NA_integer_, p); for (q in seq_along(map)) v[map[[q]]] <- q; v }
all_perms <- function(n) {
  if (n == 1) return(matrix(1L, 1, 1))
  sp <- all_perms(n - 1)
  do.call(rbind, lapply(seq_len(n), function(i) cbind(i, matrix(c(seq_len(n))[-i][sp], nrow(sp)))))
}
# Exact one-to-one matching of groups of partition `a` to groups of partition `b` (max shared items):
# returns m such that group m[q] of `a` is matched to group q of `b`
match_groups <- function(a, b, K) {
  O <- table(factor(b, levels = 1:K), factor(a, levels = 1:K))
  P <- all_perms(K)
  score <- rowSums(matrix(O[cbind(rep(1:K, each = nrow(P)), as.vector(P))], nrow(P)))
  P[which.max(score), ]
}
find_data_file <- function(dir, candidates, pattern) {
  for (f in candidates) if (file.exists(file.path(dir, f))) return(file.path(dir, f))
  hits <- list.files(dir, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
  if (length(hits)) return(hits[1])
  stop("No file matching '", pattern, "' in ", normalizePath(dir, mustWork = FALSE), ". CSV files there: ",
       paste(list.files(dir, pattern = "\\.csv$", ignore.case = TRUE), collapse = ", "), call. = FALSE)
}

## ============================================================================================
## 1. DATA AND MEASURES (Section 6.1)
## ============================================================================================
mbi_file  <- find_data_file(data_dir, c("BDI_converted DKA_Yannick.csv", "BDI_converted_DKA_Yannick.csv"),
                            "^BDI.*Yannick.*\\.csv$")
tics_file <- if (data_source == "raw") find_data_file(data_dir, c("TICS data.csv"), "^TICS[ _]data\\.csv$") else
  find_data_file(data_dir, c("TICS_Data.csv"), "^TICS_Data.*\\.csv$")
eremi_rdata <- file.path(data_dir, "E-ReMI simulation results for TICS data setseed 123.Rdata")

read_raw <- function(f) read.csv(f, sep = ";", header = TRUE, na.strings = c("NULL", "NA", "", " "),
                                 stringsAsFactors = FALSE)
numbered_cols <- function(df, prefix, n_expected) {      # e.g. TICS_1_... to TICS_57_..., in item order
  hit <- grep(paste0("^", prefix, "_([0-9]+)(_|$)"), names(df), value = TRUE)
  num <- as.integer(sub(paste0("^", prefix, "_([0-9]+).*$"), "\\1", hit))
  hit <- hit[order(num)]
  if (length(hit) != n_expected || !identical(sort(num), seq_len(n_expected)))
    stop("Expected ", prefix, "_1 ... ", prefix, "_", n_expected, " columns.", call. = FALSE)
  hit
}
to_num <- function(df) as.data.frame(lapply(df, function(v) suppressWarnings(as.numeric(sub(",", ".", as.character(v), fixed = TRUE)))))

mbi_raw <- read_raw(mbi_file)
mbi_item_cols <- numbered_cols(mbi_raw, "MBI", 14)

if (data_source == "raw") {
  tics_raw <- read_raw(tics_file)
  tics_item_cols <- numbered_cols(tics_raw, "TICS", 57)
  id_t <- names(tics_raw)[1]; id_m <- names(mbi_raw)[1]
  matched <- dplyr::inner_join(tics_raw, mbi_raw, by = setNames(id_m, id_t), suffix = c("", ".mbi"))
  imp_cols <- c(tics_item_cols, mbi_item_cols, intersect("MBI_total", names(matched)))
  imp_in   <- to_num(matched[, imp_cols])
  report("Participants with TICS and MBI data: ", nrow(matched),
         "; imputed TICS responses: ", sum(is.na(imp_in[, tics_item_cols])), " of ", nrow(matched) * 57,
         "; participants with imputed MBI responses: ", sum(!complete.cases(imp_in[, mbi_item_cols])))
  set.seed(seed)
  imp <- VIM::kNN(imp_in, k = k_impute, imp_var = FALSE)   # median of the k nearest neighbours
  X   <- as.matrix(imp[, tics_item_cols])
  EE  <- rowMeans(imp[, mbi_item_cols[1:9]])               # emotional exhaustion: MBI 1-9
  DP  <- rowMeans(imp[, mbi_item_cols[10:14]])             # depersonalization:   MBI 10-14
  age <- suppressWarnings(as.numeric(matched$Alter))
} else {
  tics <- read.csv(tics_file)
  X   <- as.matrix(tics[, 1:57]); EE <- tics$Emotional_Exhaustion; DP <- tics$Depersonalization
  age <- suppressWarnings(as.numeric(mbi_raw$Alter))
}
colnames(X) <- paste0("TICS_", 1:57)
stopifnot(!anyNA(X), all(X >= 0 & X <= 4), !anyNA(EE), !anyNA(DP),
          nrow(X) == length(EE), nrow(X) == length(DP), length(age) == nrow(X))
if (nrow(X) != 563) warning("Sample size is ", nrow(X), "; the manuscript reports 563 participants.")

level <- rowMeans(X)   # overall stress level = participant's mean across all 57 TICS items (0-4)
if (all(is.na(age))) {
  report("Age: column 'Alter' not found or empty; age not reported.")
} else report(sprintf("Age: M = %.1f, SD = %.1f, range %.0f-%.0f (n with age = %d)", mean(age, na.rm = TRUE),
                      sd(age, na.rm = TRUE), min(age, na.rm = TRUE), max(age, na.rm = TRUE), sum(!is.na(age))))
report(sprintf("Overall stress level: M = %.2f, SD = %.2f; emotional exhaustion: M = %.2f, SD = %.2f; depersonalization: M = %.2f, SD = %.2f",
               mean(level), sd(level), mean(EE), sd(EE), mean(DP), sd(DP)))

## ============================================================================================
## 2. ITEM STRUCTURES (Section 6.2)
## ============================================================================================
# Interpretive names of the nine E-ReMI item clusters (as in Figure 8)
domain_labels <- c("Work overload", "Social overload", "Pressure to perform", "Work discontent",
                   "Excessive demands at work", "Lack of social recognition", "Social tensions",
                   "Social isolation", "Chronic worrying")

## ---- 2a. Suitability for factor analysis and parallel analysis (Figure 7) -------------------
r57  <- cor(X)   # Pearson correlations: these reproduce KMO = 0.94 and chi2 = 19,346 reported in the text
kmo  <- psych::KMO(r57)$MSA
bart <- psych::cortest.bartlett(r57, n = nrow(X))
report(sprintf("KMO = %.2f; Bartlett chi2(%d) = %.0f, p %s", kmo, bart$df, bart$chisq, fmt_p(bart$p.value)))

pc57 <- polychoric(X)$rho
make_parallel_plot <- function() fa.parallel(pc57, n.obs = nrow(X), fa = "fa", fm = "minres",
                                             n.iter = 100, main = "")
grDevices::cairo_pdf(file.path(fig_dir, "fig07_parallel_analysis.pdf"), width = 7, height = 5)
set.seed(seed); invisible(capture.output(pa <- make_parallel_plot()))
invisible(dev.off())
ragg::agg_png(file.path(fig_dir, "fig07_parallel_analysis.png"), width = 7, height = 5, units = "in", res = 600)
set.seed(seed); invisible(capture.output(make_parallel_plot()))
invisible(dev.off())
report("Parallel analysis suggests ", pa$nfact, " factors.")

## ---- 2b. PCEFA with iterative refinement ------------------------------------------------------
run_pcefa <- function(Xs, nf = 9, pc = NULL) {
  if (is.null(pc)) pc <- polychoric(Xs)$rho
  fa(r = pc, nfactors = nf, n.obs = nrow(Xs), rotate = "oblimin", fm = "minres")
}
drop_iter1 <- c(47, 13, 37, 14, 40, 43, 30, 12, 57, 29, 18, 9)
drop_iter2 <- c(5, 11, 8)
efa_1 <- run_pcefa(X, pc = pc57)
# Check: every item removed in iteration 1 meets the removal rule (complexity > 2 and/or uniqueness > .50)
flag1 <- which(efa_1$complexity > 2 | efa_1$uniquenesses > .50)
if (!all(drop_iter1 %in% flag1))
  warning("Removed item(s) ", paste(setdiff(drop_iter1, flag1), collapse = ", "),
          " do not meet the complexity/uniqueness rule in iteration 1.")
efa_2 <- run_pcefa(X[, -drop_iter1])
keep  <- setdiff(1:57, c(drop_iter1, drop_iter2))
X42   <- X[, keep]
efa_3 <- run_pcefa(X42)
var_expl <- sapply(list(efa_1, efa_2, efa_3), function(f) 100 * f$Vaccounted["Cumulative Var", 9])
report("PCEFA variance explained (iterations 1-3): ", paste(sprintf("%.2f%%", var_expl), collapse = ", "),
       "; items retained: 57 -> ", 57 - length(drop_iter1), " -> ", length(keep))
pcefa_cl42 <- apply(abs(unclass(efa_3$loadings)), 1, which.max)   # item -> factor (largest |loading|)
pcefa_cl   <- rep(NA_integer_, 57); pcefa_cl[keep] <- pcefa_cl42

## ---- 2c. E-ReMI item partition ---------------------------------------------------------------
eremi_map <- list(c(1,4,17,38,44,50,54), c(7,39,43,49), c(8,14,22,23,30,32), c(10,11,21,41,42),
                  c(2,3,9,13,18,20,24,35,55,57), c(19,27,28), c(5,6,12,15,26,29,31,33,37,40,45,46,52),
                  c(34,48,51,53,56), c(16,25,36,47))
eremi_cl <- map_to_vec(eremi_map)
row_membership_from_file <- NULL
if (file.exists(eremi_rdata)) {
  e_env <- new.env(); load(eremi_rdata, envir = e_env)
  if (exists("K_EReMI_3_9", envir = e_env)) {
    K <- as.matrix(get("K_EReMI_3_9", envir = e_env))
    if (nrow(K) != 57 && ncol(K) == 57) K <- t(K)
    k_cl <- apply(K, 1, which.max)
    eremi_cl <- match(k_cl, match_groups(k_cl, eremi_cl, 9))   # keep the manuscript's cluster numbering
  }
  for (nm in ls(e_env)) {                                       # stored participant memberships, if any
    o <- get(nm, envir = e_env)
    if (is.matrix(o) && nrow(o) == nrow(X) && ncol(o) == P_eremi && all(o %in% 0:1) && all(rowSums(o) == 1)) {
      row_membership_from_file <- apply(o, 1, which.max); break
    }
  }
}
stopifnot(!anyNA(eremi_cl))

## ---- 2d. Reliability and agreement --------------------------------------------------------------
alpha_of <- function(cols) suppressWarnings(psych::alpha(X[, cols, drop = FALSE], check.keys = FALSE)$total$raw_alpha)
alpha_eremi <- sapply(1:9, function(q) alpha_of(which(eremi_cl == q)))
alpha_pcefa <- sapply(1:9, function(f) alpha_of(keep[pcefa_cl42 == f]))
report(sprintf("Cronbach's alpha: E-ReMI clusters %.2f-%.2f; PCEFA factors %.2f-%.2f",
               min(alpha_eremi), max(alpha_eremi), min(alpha_pcefa), max(alpha_pcefa)))

fac_for_cluster <- match_groups(pcefa_cl42, eremi_cl[keep], 9)   # PCEFA factor matched to E-ReMI cluster q
xtab <- table(`PCEFA factor` = factor(ifelse(is.na(pcefa_cl), "Removed", paste0("F", pcefa_cl)),
                                      levels = c(paste0("F", fac_for_cluster), "Removed")),
              `E-ReMI item cluster` = factor(paste0("C", eremi_cl), levels = paste0("C", 1:9)))
report(sprintf("ARI (E-ReMI vs PCEFA, 42 items) = %.3f; items in matched group: %d of %d",
               ari(eremi_cl[keep], pcefa_cl[keep]), sum(diag(xtab[1:9, ])), length(keep)))
removed_row <- xtab["Removed", ]
report(sprintf("Items removed by PCEFA: %d, spread over %d E-ReMI clusters; most in %s (%d items)",
               sum(removed_row), sum(removed_row > 0), names(which.max(removed_row)), max(removed_row)))
report("Cross-tabulation PCEFA x E-ReMI (Figure 9):", obj = xtab)
save_table(cbind(PCEFA = rownames(xtab), as.data.frame.matrix(xtab)), "fig09_crosstab_pcefa_eremi")

## ---- 2e. Figure 8: PCEFA structure (left) vs E-ReMI structure (right) ----------------------------
fig8_labels  <- tools::toTitleCase(domain_labels)
fig8_colours <- c("red", "blue", "forestgreen", "purple", "orange", "brown", "#008080", "darkgoldenrod1", "darkred")
groups_pcefa <- lapply(1:9, function(q) sort(keep[pcefa_cl42 == fac_for_cluster[q]]))
groups_eremi <- lapply(1:9, function(q) sort(which(eremi_cl == q)))
build_structure_dot <- function(groups, centre_label) {
  boxes <- vapply(1:9, function(q) sprintf(
    '  subgraph cluster_%d {\n    label = "Instrument %d: (%s)";\n    color = "%s";\n    %s;\n  }',
    q, q, fig8_labels[q], fig8_colours[q], paste(groups[[q]], collapse = "; ")), character(1))
  edges <- vapply(1:9, function(q) sprintf('  Factors -> {%s} [color = "%s"]',
                                           paste(groups[[q]], collapse = " "), fig8_colours[q]), character(1))
  paste0('digraph structure_mapping {
  graph [layout = dot, rankdir = LR, fontname = "Helvetica-Bold", fontsize = 26, nodesep = 0.3, ranksep = 1.0]
  node  [shape = rectangle, style = filled, fontname = "Helvetica-Bold", fontsize = 22, margin = "0.18,0.10", penwidth = 1.3]
  edge  [penwidth = 1.8, arrowsize = 1.0]
', paste(boxes, collapse = "\n"), '
  Factors [shape = ellipse, label = "', centre_label, '", fillcolor = gray95, fontname = "Helvetica-Bold", fontsize = 30]
', paste(edges, collapse = "\n"), '
}')
}
dot_to_svg <- function(dot_code) {      # via DiagrammeR, or the Graphviz program `dot` if installed
  if (has["DiagrammeR"] && has["DiagrammeRsvg"]) return(DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot_code)))
  if (nzchar(Sys.which("dot"))) {
    f_in <- tempfile(fileext = ".gv"); f_out <- tempfile(fileext = ".svg")
    writeLines(dot_code, f_in); system2(Sys.which("dot"), c("-Tsvg", shQuote(f_in), "-o", shQuote(f_out)))
    return(paste(readLines(f_out, warn = FALSE), collapse = "\n"))
  }
  NULL
}
# Place the two panels side by side in one vector SVG (no panel titles), then export PDF / PNG / TIFF
combine_svgs <- function(svgs, filename, png_width_px = 5000, dpi = 600, gap = 60, margin = 40) {
  base <- file.path(fig_dir, filename)
  dims <- lapply(svgs, function(s) {
    vb <- regmatches(s, regexpr('viewBox="[^"]+"', s, perl = TRUE))
    as.numeric(strsplit(gsub('viewBox="|"', "", vb), "[ ,]+")[[1]])[3:4] })
  W <- sum(sapply(dims, `[`, 1)) + gap * (length(svgs) - 1) + 2 * margin
  H <- max(sapply(dims, `[`, 2)) + 2 * margin
  x <- margin; parts <- character(0)
  for (i in seq_along(svgs)) {
    inner <- sub("^[\\s\\S]*?(<svg\\b)", "\\1", svgs[[i]], perl = TRUE)
    inner <- sub('(<svg[^>]*?\\swidth=")[^"]*"',  sprintf('\\1%.2f"', dims[[i]][1]), inner, perl = TRUE)
    inner <- sub('(<svg[^>]*?\\sheight=")[^"]*"', sprintf('\\1%.2f"', dims[[i]][2]), inner, perl = TRUE)
    parts <- c(parts, sub("<svg ", sprintf('<svg x="%.1f" y="%.1f" ', x, margin), inner))
    x <- x + dims[[i]][1] + gap
  }
  combined <- paste0('<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" ',
                     sprintf('width="%.0fpt" height="%.0fpt" viewBox="0 0 %.0f %.0f">', W, H, W, H),
                     sprintf('<rect width="%.0f" height="%.0f" fill="white"/>', W, H),
                     paste(parts, collapse = "\n"), "</svg>")
  rsvg::rsvg_pdf(charToRaw(combined), paste0(base, ".pdf"))
  rsvg::rsvg_png(charToRaw(combined), paste0(base, ".png"), width = png_width_px)
  if (has["magick"]) tryCatch({
    img <- magick::image_background(magick::image_read(paste0(base, ".png")), "white", flatten = TRUE)
    magick::image_write(img, paste0(base, ".tiff"), format = "tiff", compression = "LZW", density = paste0(dpi, "x", dpi))
  }, error = function(e) message("TIFF export of Figure 8 failed; use the PDF/PNG."))
  message("Saved: ", base, ".pdf / .png / .tiff")
}
svg_pcefa <- dot_to_svg(build_structure_dot(groups_pcefa, "PCEFA Structure\\n(9 Constructs, 42 Items)"))
svg_eremi <- dot_to_svg(build_structure_dot(groups_eremi, "E-ReMI Structure\\n(9 Constructs, 57 Items)"))
if (!is.null(svg_pcefa) && !is.null(svg_eremi) && has["rsvg"]) {
  tryCatch(combine_svgs(list(svg_pcefa, svg_eremi), "fig08_item_structures"),
           error = function(e) message("Figure 8 could not be drawn: ", conditionMessage(e)))
} else message("Figure 8 not drawn: install DiagrammeR, DiagrammeRsvg and rsvg (or Graphviz).")

## ---- 2f. Figure 9: cross-tabulation heatmap -----------------------------------------------------
xt_df <- as.data.frame(xtab, responseName = "n")
names(xt_df)[1:2] <- c("PCEFA", "EReMI")
xt_df$EReMI <- factor(xt_df$EReMI, levels = paste0("C", 1:9), labels = paste0("C", 1:9, ": ", domain_labels))
xt_df$PCEFA <- factor(xt_df$PCEFA, levels = rev(levels(xt_df$PCEFA)))
p_xtab <- ggplot(xt_df, aes(EReMI, PCEFA, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(n > 0, n, "")), size = 3.2) +
  scale_fill_gradient(low = "white", high = "#2166AC", name = "Items") +
  labs(x = "E-ReMI item cluster", y = "PCEFA factor") +
  theme_paper + theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
save_ggplot_hires(p_xtab, "fig09_crosstab_pcefa_eremi", width_mm = 174, height_mm = 110)

## ============================================================================================
## 3. PARTICIPANT CLUSTERS (Section 6.3)
## ============================================================================================
## E-ReMI conditions the likelihood on the row means, so the interaction effects w_pq are means of the
## doubly-centred data within biclusters. Given the item partition, each participant is assigned to the
## cluster that maximises the conditional classification likelihood (unequal cluster sizes allowed).
dc <- function(Y) { Y <- as.matrix(Y)
Y - outer(rowMeans(Y), rep(1, ncol(Y))) - outer(rep(1, nrow(Y)), colMeans(Y)) + mean(Y) }
eremi_rows <- function(Y, ccl, P, nstart = 100, maxit = 100, init = NULL, seed = 1) {
  D <- dc(Y); n <- nrow(D); Q <- max(ccl); nq <- tabulate(ccl, Q)
  Dq   <- sapply(1:Q, function(q) rowMeans(D[, ccl == q, drop = FALSE]))
  base <- sum(D^2) - sum(sweep(Dq^2, 2, nq, "*"))
  fit_from <- function(r) {
    for (it in seq_len(maxit)) {
      pis <- tabulate(r, P) / n
      if (any(pis == 0)) return(list(loss = Inf))
      W   <- t(sapply(1:P, function(p) colMeans(Dq[r == p, , drop = FALSE])))
      sse <- base + sum(sapply(1:P, function(p) sum(sweep(sweep(Dq[r == p, , drop = FALSE], 2, W[p, ])^2, 2, nq, "*"))))
      s2  <- sse / (n * ncol(D))
      cost <- sapply(1:P, function(p) colSums((t(Dq) - W[p, ])^2 * nq)) - 2 * s2 * matrix(log(pis), n, P, byrow = TRUE)
      rn <- max.col(-cost, ties.method = "first")
      if (all(rn == r)) break
      r <- rn
    }
    pis <- tabulate(r, P) / n
    if (any(pis == 0)) return(list(loss = Inf))
    W   <- t(sapply(1:P, function(p) colMeans(Dq[r == p, , drop = FALSE])))
    sse <- base + sum(sapply(1:P, function(p) sum(sweep(sweep(Dq[r == p, , drop = FALSE], 2, W[p, ])^2, 2, nq, "*"))))
    s2  <- sse / (n * ncol(D))
    list(r = r, W = W, loss = n * ncol(D) / 2 * log(s2) - sum(n * pis * log(pis)))
  }
  set.seed(seed)
  starts <- c(if (!is.null(init)) list(init), replicate(nstart, sample(rep(1:P, length.out = n)), simplify = FALSE))
  fits <- lapply(starts, fit_from)
  fits[[which.min(sapply(fits, `[[`, "loss"))]]
}
er_fit <- if (!is.null(row_membership_from_file)) {
  eremi_rows(X, eremi_cl, P_eremi, nstart = 0, init = row_membership_from_file, maxit = 0)
} else eremi_rows(X, eremi_cl, P_eremi, nstart = 100, seed = seed)

ord   <- order(-tabulate(er_fit$r, P_eremi))          # P1 = largest cluster
er_cl <- match(er_fit$r, ord)
W_er  <- er_fit$W[ord, , drop = FALSE]
er_sizes <- tabulate(er_cl, P_eremi)
report("Participant clusters: ", paste(sprintf("n = %d (%.1f%%)", er_sizes, 100 * er_sizes / sum(er_sizes)), collapse = ", "),
       if (is.null(row_membership_from_file)) " [re-estimated: stored E-ReMI solution not found]" else " [stored E-ReMI solution]")
if (!identical(as.integer(er_sizes), c(215L, 195L, 153L)))
  warning("Cluster sizes (", paste(er_sizes, collapse = "/"), ") differ from the manuscript (215/195/153). ",
          "Place the E-ReMI .Rdata file in data_dir to use the reported solution.")
W_tab <- data.frame(Cluster = paste0("P", 1:P_eremi, " (n = ", er_sizes, ")"), round(W_er, 2), check.names = FALSE)
names(W_tab)[-1] <- paste0("C", 1:9, ": ", domain_labels)
report("Interaction effects w_pq (Figure 10):", obj = W_tab)
save_table(W_tab, "fig10_interaction_effects")

## ---- Figure 10: interaction heatmap -------------------------------------------------------------
col_lab <- paste0("C", 1:9, ": ", domain_labels)
hm <- expand.grid(person = factor(W_tab$Cluster, levels = rev(W_tab$Cluster)), item = factor(col_lab, levels = col_lab))
hm$value <- as.vector(W_er[cbind(match(hm$person, W_tab$Cluster), match(hm$item, col_lab))])
lim <- max(abs(hm$value))
p_heat <- ggplot(hm, aes(item, person, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim),
                       name = expression(hat(w)[pq])) +
  labs(x = "E-ReMI item cluster", y = "E-ReMI person cluster") +
  theme_paper + theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank())
save_ggplot_hires(p_heat, "fig10_eremi_person_item_heatmap", width_mm = 174, height_mm = 85)

## ---- One-way ANOVA (eta^2) and Kruskal-Wallis test across the participant clusters --------------
assoc_one <- function(y, name) {
  ok <- !is.na(y); g <- factor(er_cl[ok]); a <- summary(aov(y[ok] ~ g))[[1]]
  data.frame(Variable = name,
             M_SD = paste(sprintf("%.2f (%.2f)", tapply(y[ok], g, mean), tapply(y[ok], g, sd)), collapse = "; "),
             F = round(a[1, "F value"], 2), df = paste0(a[1, "Df"], ", ", a[2, "Df"]), p = fmt_p(a[1, "Pr(>F)"]),
             eta2 = round(a[1, "Sum Sq"] / sum(a[, "Sum Sq"]), 3), KW_p = fmt_p(kruskal.test(y[ok], g)$p.value))
}
## ---- Do the clusters differ in overall stress level? (one-way ANOVA, Tukey HSD) -----------------
## E-ReMI removes each participant's mean before clustering, so cluster membership is invariant to a
## constant shift in a person's responses. The clusters can nevertheless differ in level empirically
## if certain response patterns occur more often at high or low stress levels. This is tested here.
lvl_test <- assoc_one(level, "Overall stress level")
lvl_cl   <- factor(paste0("P", er_cl))
tukey    <- TukeyHSD(aov(level ~ lvl_cl))$lvl_cl
tukey_df <- data.frame(Comparison = rownames(tukey), Diff = round(tukey[, "diff"], 2),
                       CI_low = round(tukey[, "lwr"], 2), CI_high = round(tukey[, "upr"], 2),
                       p_adj = fmt_p(tukey[, "p adj"]))
lvl_means <- tapply(level, lvl_cl, mean)
d_max_pts <- max(abs(outer(lvl_means, lvl_means, "-")))      # in scale points (0-4)
report("Overall stress level by participant cluster (M (SD) for P1; P2; P3):", obj = lvl_test)
report("Tukey HSD pairwise comparisons of overall stress level:", obj = tukey_df)
report(sprintf("Largest difference between cluster means in overall stress level: %.2f scale points; largest |w_pq| = %.2f",
               d_max_pts, max(abs(W_er))))
save_table(rbind(lvl_test), "sec6_3_cluster_stress_level")
save_table(tukey_df, "sec6_3_cluster_stress_level_tukey")

## ---- Incremental validity of cluster membership over overall stress level -----------------------
incr <- do.call(rbind, lapply(list(`Emotional exhaustion` = EE, Depersonalization = DP), function(y) {
  m0 <- lm(y ~ level); m1 <- lm(y ~ level + factor(er_cl)); a <- anova(m0, m1)
  data.frame(dR2 = round(summary(m1)$r.squared - summary(m0)$r.squared, 3), F = round(a$F[2], 2),
             df = paste0(a$Df[2], ", ", a$Res.Df[2]), p = fmt_p(a$`Pr(>F)`[2])) }))
report("Cluster membership added to overall stress level:", obj = incr)
save_table(cbind(Outcome = rownames(incr), incr), "sec6_3_incremental_validity")

## ============================================================================================
## 4. PREDICTING BURNOUT (Section 6.4)
## ============================================================================================
factor_scores <- as.data.frame(factor.scores(X42, efa_3, method = "regression")$scores)
names(factor_scores) <- paste0("Factor", 1:9)
cluster_scores <- as.data.frame(sapply(1:9, function(q) rowMeans(X[, eremi_cl == q, drop = FALSE])))
names(cluster_scores) <- paste0("Cluster", 1:9)
dat <- data.frame(EE = EE, DP = DP, factor_scores, cluster_scores)

## ---- 4a. Elastic net with 10-fold CV (Section 6.4.1) -------------------------------------------
fit_elastic_net_cv <- function(Xp, y, seed = 123, nfolds = 10, alpha_grid = seq(0.1, 0.9, by = 0.1)) {
  Xp <- as.matrix(Xp); set.seed(seed)
  foldid <- sample(rep(1:nfolds, length.out = nrow(Xp)))     # identical folds for every alpha
  r2 <- sapply(alpha_grid, function(a)
    1 - min(cv.glmnet(Xp, y, family = "gaussian", alpha = a, foldid = foldid, type.measure = "mse")$cvm) / var(y))
  max(r2)
}
cv_r2 <- data.frame(Outcome = c("Depersonalization", "Emotional exhaustion"),
                    Factors  = c(fit_elastic_net_cv(factor_scores, DP, seed), fit_elastic_net_cv(factor_scores, EE, seed)),
                    Clusters = c(fit_elastic_net_cv(cluster_scores, DP, seed), fit_elastic_net_cv(cluster_scores, EE, seed)))
report("Cross-validated R^2 (elastic net, single 10-fold split):", obj = mutate(cv_r2, across(-Outcome, ~ round(.x, 3))))
save_table(cv_r2, "sec6_4_1_elastic_net_cv_r2")

## ---- 4b. Repeated cross-validation (Section 6.4.2, Figure 11) -----------------------------------
rep_full <- do.call(rbind, lapply(seq_len(n_reps_cv), function(r) {
  f_dp <- fit_elastic_net_cv(factor_scores, DP, seed = r); c_dp <- fit_elastic_net_cv(cluster_scores, DP, seed = r)
  f_ee <- fit_elastic_net_cv(factor_scores, EE, seed = r); c_ee <- fit_elastic_net_cv(cluster_scores, EE, seed = r)
  data.frame(rep = r, Outcome = c("Depersonalization", "Emotional Exhaustion"),
             R2_factors = c(f_dp, f_ee), R2_clusters = c(c_dp, c_ee))
}))
rep_full$R2_diff <- rep_full$R2_clusters - rep_full$R2_factors
save_table(rep_full, "fig11_repeated_cv_all_repetitions")
rep_long <- rep_full[, c("rep", "Outcome", "R2_diff")]
rep_sum <- rep_long %>% group_by(Outcome) %>%
  summarise(Mean_diff_pp = 100 * mean(R2_diff), SD_pp = 100 * sd(R2_diff),
            CI_lower_pp = 100 * (mean(R2_diff) - qt(.975, n() - 1) * sd(R2_diff) / sqrt(n())),
            CI_upper_pp = 100 * (mean(R2_diff) + qt(.975, n() - 1) * sd(R2_diff) / sqrt(n())),
            Reps_favoring_clusters = sum(R2_diff > 0), .groups = "drop") %>%
  left_join(rep_full %>% group_by(Outcome) %>% summarise(Mean_R2_factors = mean(R2_factors), .groups = "drop"),
            by = "Outcome") %>%
  mutate(Relative_gain_pct = Mean_diff_pp / (100 * Mean_R2_factors) * 100)   # mean gain relative to PCEFA baseline
report(sprintf("Repeated CV over %d fold assignments, clusters minus factors (percentage points; 95%% CI of the mean; relative gain vs. mean PCEFA R^2):", n_reps_cv),
       obj = mutate(as.data.frame(rep_sum), across(-Outcome, ~ round(.x, 2))))
save_table(rep_sum, "fig11_repeated_cv_summary")
p_cv <- ggplot(rep_long, aes(Outcome, R2_diff, fill = Outcome)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.3, size = 1) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), guide = "none") +
  labs(x = NULL, y = expression(Delta * " cross-validated " * R^2 * " (clusters " - " factors)")) + theme_paper
save_ggplot_hires(p_cv, "fig11_repeated_cv_r2_difference", width_mm = 140, height_mm = 110)

## ---- 4c. Stepwise models and Figure 12 (Section 6.4.3) -------------------------------------------
## The predictors below are the ones retained by the stepwise selection reported in Section 6.4.3
## (all p < .05); the models are refitted here to obtain the fit indices and coefficients in Figure 12.
models <- list(
  `Depersonalization (Factors)`     = lm(DP ~ Factor2 + Factor5 + Factor8, data = dat),
  `Depersonalization (Clusters)`    = lm(DP ~ Cluster1 + Cluster2 + Cluster5 + Cluster7 + Cluster8, data = dat),
  `Emotional Exhaustion (Factors)`  = lm(EE ~ Factor1 + Factor3 + Factor4 + Factor6 + Factor8 + Factor9, data = dat),
  `Emotional Exhaustion (Clusters)` = lm(EE ~ Cluster1 + Cluster4 + Cluster5 + Cluster7, data = dat))
fit_indices <- data.frame(Model = names(models), R2 = sapply(models, function(m) summary(m)$r.squared),
                          AdjR2 = sapply(models, function(m) summary(m)$adj.r.squared),
                          AIC = sapply(models, AIC), BIC = sapply(models, BIC), row.names = NULL)
report("Stepwise models, fit indices:", obj = mutate(fit_indices, across(-Model, ~ round(.x, 3))))
pmax_sel <- sapply(models, function(m) max(summary(m)$coefficients[-1, 4]))
if (any(pmax_sel >= .05)) warning("Not all retained predictors have p < .05 in: ",
                                   paste(names(pmax_sel)[pmax_sel >= .05], collapse = "; "))
save_table(fit_indices, "fig12_stepwise_fit")
coef_df <- do.call(rbind, lapply(names(models), function(nm) {
  s <- summary(models[[nm]])$coefficients[-1, , drop = FALSE]
  data.frame(Model = nm, Predictor = rownames(s), Estimate = s[, 1], SE = s[, 2], p = s[, 4]) }))
coef_df$signif <- as.character(cut(coef_df$p, c(-Inf, .001, .01, .05, Inf), labels = c("***", "**", "*", "")))
save_table(coef_df, "fig12_stepwise_coefficients")
coef_df <- coef_df %>% left_join(fit_indices, by = "Model") %>%
  mutate(FacetLabel = sprintf("%s\nR\u00b2 = %.3f | Adj.R\u00b2 = %.3f\nAIC = %.1f | BIC = %.1f", Model, R2, AdjR2, AIC, BIC)) %>%
  arrange(Model, Estimate) %>%
  mutate(PredKey = factor(paste(Model, Predictor, sep = "__"), levels = unique(paste(Model, Predictor, sep = "__"))))
p_coef <- ggplot(coef_df, aes(x = PredKey, y = Estimate)) +
  scale_x_discrete(labels = function(x) sub("^.*__", "", x)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = Estimate - 2 * SE, ymax = Estimate + 2 * SE), width = 0.15, colour = "#1f78b4") +
  geom_point(size = 2.6, colour = "#1f78b4") +
  geom_text(aes(y = Estimate + 2 * SE, label = paste0(sprintf("%.2f", Estimate), signif)), hjust = -0.15, size = 2.6) +
  coord_flip(clip = "off") + facet_wrap(~ FacetLabel, scales = "free_y", ncol = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.24))) +
  labs(x = NULL, y = "Unstandardized coefficient (B \u00b1 2 SE)") +
  theme_paper + theme(strip.text = element_text(size = 7, lineheight = 1.05))
save_ggplot_hires(p_coef, "fig12_coefficient_plot", width_mm = 200, height_mm = 180)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
report("\nDone. Figures in ", fig_dir, ", tables in ", tab_dir, ".")