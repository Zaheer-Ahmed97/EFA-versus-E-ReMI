#########################################################################
## BOOTEFA simulation study: merge results + publication-style graphics
##
## Merges the three per-sample-size result files (n = 100 / 300 / 500)
## and builds four figures, each broken down by:
##    - Method     -> colour (legend)
##    - Sample size  -> x-axis
##    - Items (p)    -> its own grey strip, 30 vs 60 shown side-by-side,
##                      separated by a vertical divider line
##    - K            -> its own grey strip (facet), nested inside each
##                      items block
##
##   (i)   R_squared_true         - recovery of the true factor space
##   (ii)  ARI_true               - recovery of the true item partition
##   (iii) R_squared_true_ratio   - R^2 relative to the max achievable R^2
##   (iv)  Selected_K / K_Status  - under-/correct-/over-factoring rates
##
## Every figure is exported as a true vector PDF (journal-ready, scales to
## any size) plus a 600 dpi PNG (quick previewing).
#########################################################################

## ---- 0. Setup ---------------------------------------------------------

# Folder containing the three sim_results_*.csv files. If this exact folder
# doesn't contain them, the script automatically searches nearby folders
# (see step 1) -- you do not need to get this 100% right by hand.
data_dir <- "C:/Users/ahmed/Desktop/IfADo Projects/Third Project based on Survey and DVS Data/BOOTEFA-Server-number cruncher"

# Where the merged data + figures will be written. Does not need to exist.
out_dir  <- "C:/Users/ahmed/Desktop/IfADo Projects/Third Project based on Survey and DVS Data/BOOTEFA-Server-number cruncher/outputs"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "scales", "forcats", "patchwork")
missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
       ". Install them first, e.g. install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(forcats)
  library(patchwork)
})

cat("R version:", R.version.string, "\n")
cat("ggplot2:", as.character(packageVersion("ggplot2")),
    "| patchwork:", as.character(packageVersion("patchwork")), "\n\n")

## ---- 1. Locate the three input files -------------------------------------
## Extracted zip folders rarely land exactly where you expect (double-nested
## folder names, a different destination, etc.), so instead of failing hard
## on the very first guess, this looks in a few sensible places automatically.

required_files <- c(
  "sim_results_n100_sessionpilot_reduced.csv",
  "sim_results_n300_sessionpilot_reduced.csv",
  "sim_results_n500_sessionpilot_reduced.csv"
)

dir_has_all_files <- function(d) {
  dir.exists(d) && all(file.exists(file.path(d, required_files)))
}

locate_data_dir <- function(primary_dir, search_roots, required_files) {
  if (dir_has_all_files(primary_dir)) return(normalizePath(primary_dir))
  
  cat("`data_dir` did not contain the expected files -- searching nearby folders automatically...\n")
  
  for (root in unique(search_roots)) {
    if (!isTRUE(dir.exists(root))) next
    cat("  Searching under:", root, "(this can take a few seconds)\n")
    
    pattern <- paste0("^(", paste(required_files, collapse = "|"), ")$")
    hits <- tryCatch(
      list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE, ignore.case = TRUE),
      error = function(e) character(0)
    )
    if (length(hits) == 0) next
    
    candidate_dirs <- unique(dirname(hits))
    good_dirs <- candidate_dirs[vapply(candidate_dirs, dir_has_all_files, logical(1))]
    
    if (length(good_dirs) >= 1) {
      if (length(good_dirs) > 1) {
        warning("Multiple folders contain all three files; using the first one found:\n  ",
                paste(good_dirs, collapse = "\n  "),
                "\nSet `data_dir` explicitly at the top of the script if this is the wrong one.")
      }
      return(normalizePath(good_dirs[1]))
    }
  }
  NULL
}

search_roots <- c(
  dirname(data_dir),                                                   # one level up
  dirname(dirname(data_dir)),                                          # two levels up
  "C:/Users/ahmed/Desktop",                                            # broad fallback
  Sys.getenv("USERPROFILE")                                            # whole profile (OneDrive-proof)
)

resolved_dir <- locate_data_dir(data_dir, search_roots, required_files)

if (is.null(resolved_dir)) {
  stop(
    "Could not find sim_results_n100/n300/n500_sessionpilot_reduced.csv anywhere under:\n  ",
    paste(unique(search_roots), collapse = "\n  "),
    "\n\nThis almost always means the zip file has not actually been extracted to disk yet ",
    "(browsing INTO a .zip in Windows Explorer is not the same as extracting it).",
    "\nRight-click BOOTEFA-Server-number cruncher.zip -> Extract All..., then re-run this script.",
    "\nIf you know the exact folder, you can also just set `data_dir` to it directly at the top."
  )
}

data_dir <- resolved_dir
cat("Using data_dir:\n ", data_dir, "\n\n")

## ---- 2. Merge the three files ------------------------------------------

files <- file.path(data_dir, required_files)

combined_data <- tryCatch(
  bind_rows(lapply(files, read.csv, stringsAsFactors = FALSE)),
  error = function(e) stop("Failed to read/merge the input CSV files: ", conditionMessage(e))
)

required_cols <- c("SampleSize", "Items", "TrueClusters", "Method",
                   "R_squared_true", "ARI_true", "R_squared_true_ratio",
                   "Selected_K", "K_Status")
missing_cols <- setdiff(required_cols, names(combined_data))
if (length(missing_cols) > 0) {
  stop("Merged data is missing expected column(s): ", paste(missing_cols, collapse = ", "))
}
if (nrow(combined_data) == 0) stop("Merged data has 0 rows -- check the input files.")

cat("Merged", nrow(combined_data), "rows from", length(files), "files\n")
cat("Rows per sample size:\n")
print(table(combined_data$SampleSize))

write.csv(combined_data, file.path(out_dir, "combined_simulation_results.csv"), row.names = FALSE)
saveRDS(combined_data, file.path(out_dir, "combined_simulation_results.rds"))

## ---- 3. Recode factors for plotting -------------------------------------

method_levels_raw <- c("EFA", "EFA-boot", "E-ReMI")
missing_methods <- setdiff(unique(combined_data$Method), method_levels_raw)
if (length(missing_methods) > 0) {
  stop("Found Method value(s) not in the expected set (EFA, EFA-boot, E-ReMI): ",
       paste(missing_methods, collapse = ", "),
       ". Update `method_levels_raw`/`method_recode`/`method_colors` to include them.")
}

# Raw CSV value -> display label shown in the figures
method_recode <- c("EFA" = "PCEFA", "EFA-boot" = "PCEFA-boot", "E-ReMI" = "E-ReMI")
method_levels <- c("PCEFA", "PCEFA-boot", "E-ReMI")   # display labels, in legend order

method_colors <- c(
  "PCEFA"      = "#0072B2",  # blue
  "PCEFA-boot" = "#009E73",  # green
  "E-ReMI"     = "#D55E00"   # vermillion
)

kstatus_levels <- c("Under", "Correct", "Over")

results <- combined_data %>%
  mutate(
    Method     = factor(method_recode[Method], levels = method_levels),
    # Clean standard text labels matching your target screenshot style: n = 100, n = 300, etc.
    SampleSize = factor(SampleSize, levels = sort(unique(SampleSize))),
    K          = factor(paste0("bolditalic(k) == ", TrueClusters),
                        levels = paste0("bolditalic(k) == ", sort(unique(TrueClusters)))),
    K_Status   = factor(K_Status, levels = kstatus_levels)
  )

items_levels <- sort(unique(results$Items))   # e.g. 30, 60
if (length(items_levels) != 2) {
  stop("Expected exactly 2 distinct values of `Items` (e.g. 30 and 60); found: ",
       paste(items_levels, collapse = ", "),
       ". The side-by-side layout below assumes 2 levels -- adjust the script if not.")
}

## ---- 4. Shared theme & helpers -------------------------------------------

theme_elegant <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = rel(1.12), hjust = 0),
      plot.subtitle    = element_text(color = "grey35", size = rel(0.95), hjust = 0,
                                      margin = margin(b = 8)),
      plot.caption     = element_text(color = "grey50", size = rel(0.75), hjust = 1),
      strip.text       = element_text(face = "bold.italic", size = 12),
      strip.background = element_rect(fill = "grey94", color = NA),
      axis.title       = element_text(face = "bold", size = rel(0.95)),
      axis.text        = element_text(color = "black", size = 12, face = "bold"),
      axis.text.x      = element_text(color = "black", size = 12, face = "bold"),
      axis.text.y      = element_text(color = "black", size = 12, face = "bold"),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(4, "pt"),
      legend.position  = "top",
      legend.title     = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(1, "lines"),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}

items_banner <- function(items_val, base_size = 12) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "grey85") +
    annotate("text", x = 0.5, y = 0.5,
             label = paste0("bolditalic(p) == ", items_val), parse = TRUE,
             size = base_size / ggplot2::.pt, color = "black") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(1, 1, 1, 1))
}

divider_v <- function() {
  ggplot() + geom_vline(xintercept = 0, linewidth = 0.9, color = "grey70") +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}
divider_h <- function() {
  ggplot() + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0),
                          linewidth = 0.9, color = "grey70") +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}

save_figure <- function(fig, basename, width, height, png_dpi = 600) {
  pdf_path <- file.path(out_dir, paste0(basename, ".pdf"))
  png_path <- file.path(out_dir, paste0(basename, ".png"))
  
  pdf_ok <- tryCatch({
    ggsave(pdf_path, fig, width = width, height = height,
           device = grDevices::cairo_pdf, bg = "white")
    TRUE
  }, error = function(e) {
    warning(sprintf("cairo_pdf failed for '%s' (%s); falling back to standard pdf device.",
                    basename, conditionMessage(e)))
    ggsave(pdf_path, fig, width = width, height = height, device = "pdf", bg = "white")
    FALSE
  })
  
  ggsave(png_path, fig, width = width, height = height, dpi = png_dpi, bg = "white")
  
  if (!file.exists(pdf_path) || !file.exists(png_path)) {
    stop("Failed to write one or more output files for '", basename, "'")
  }
  cat(sprintf("  Saved %-32s  PDF (vector%s) + PNG (%d dpi, %.0f x %.0f px)\n",
              basename, if (pdf_ok) "" else ", fallback device",
              png_dpi, width * png_dpi, height * png_dpi))
  invisible(list(pdf = pdf_path, png = png_path))
}

## ---- 5. Generic builder for the three continuous measures ---------------

build_measure_panel <- function(data, items_val, yvar, ylab,
                                y_limits = NULL, y_breaks = NULL,
                                hline = NULL) {
  
  d <- filter(data, Items == items_val)
  
  p <- ggplot(d, aes(x = SampleSize, y = .data[[yvar]], fill = Method)) +
    geom_boxplot(position = position_dodge(0.75), width = 0.65,
                 outlier.shape = NA, alpha = 0.85, linewidth = 0.35) +
    geom_point(aes(color = Method),
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
               alpha = 0.25, size = 0.7, show.legend = FALSE) +
    facet_wrap(vars(K), nrow = 1, labeller = labeller(K = label_parsed)) +
    scale_x_discrete(labels = function(x) paste0("n = ", x)) +
    scale_fill_manual(values = method_colors) +
    scale_color_manual(values = method_colors) +
    labs(x = "Sample size", y = ylab, fill = "Method") +
    theme_elegant()
  
  if (!is.null(y_limits)) p <- p + scale_y_continuous(limits = y_limits, breaks = y_breaks)
  if (!is.null(hline))    p <- p + geom_hline(yintercept = hline, linetype = "dashed",
                                              color = "grey40", linewidth = 0.4)
  p
}

combine_measure_figure <- function(data, yvar, ylab, title = NULL, caption = NULL,
                                   y_limits = NULL, y_breaks = NULL, hline = NULL){
  
  banner1 <- items_banner(items_levels[1])
  banner2 <- items_banner(items_levels[2])
  plot1   <- build_measure_panel(data, items_levels[1], yvar, ylab, y_limits, y_breaks, hline)
  plot2   <- build_measure_panel(data, items_levels[2], yvar, ylab, y_limits, y_breaks, hline)
  
  design <- "
  AEB
  CED
  "
  
  fig <- wrap_plots(A = banner1, B = banner2, C = plot1, D = plot2, E = divider_v(),
                    design = design, heights = c(0.07, 1), widths = c(10, 0.15, 10),
                    guides = "collect") +
    plot_annotation(
      title = title,
      caption = caption,
      theme = theme(
        plot.title   = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.caption = element_text(color = "grey50", size = 9, hjust = 1)
      )
    ) &
    theme(legend.position = "top")
  
  fig
}

## ---- 6. (i) R_squared_true ------------------------------------------------

fig_r2 <- tryCatch(
  combine_measure_figure(
    results,
    yvar   = "R_squared_true",
    ylab   = expression(bold(R^2~"(explained variance)")),
    y_limits = c(0.45, 0.75), y_breaks = seq(0.45, 0.75, 0.05)
  ),
  error = function(e) stop("Failed building figure 01_R_squared_true: ", conditionMessage(e))
)
save_figure(fig_r2, "01_R_squared_true", width = 13, height = 6.7)

## ---- 7. (ii) ARI_true -------------------------------------------------------

fig_ari <- tryCatch(
  combine_measure_figure(
    results,
    yvar   = "ARI_true",
    ylab   = "Adjusted Rand Index (ARI)",
    y_limits = c(0.75, 1.0), y_breaks = seq(0.75, 1.0, 0.05)
  ),
  error = function(e) stop("Failed building figure 02_ARI_true: ", conditionMessage(e))
)
save_figure(fig_ari, "02_ARI_true", width = 13, height = 6.7)

## ---- 8. (iii) R_squared_true_ratio -----------------------------------------

fig_ratio <- tryCatch(
  combine_measure_figure(
    results,
    yvar   = "R_squared_true_ratio",
    ylab   = expression(bold(R^2~"ratio  ("*R[true]^2*" / "*R[max]^2*")")),
    y_limits = c(0.9, 1.35), y_breaks = seq(0.9, 1.35, 0.1),
    hline = 1
  ),
  error = function(e) stop("Failed building figure 03_R_squared_true_ratio: ", conditionMessage(e))
)
save_figure(fig_ratio, "03_R_squared_true_ratio", width = 13, height = 6.7)

## ---- 9. (iv) Selected_K / K_Status -----------------------------------------

k_status_summary <- results %>%
  group_by(Items, K, SampleSize, Method, K_Status, .drop = FALSE) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Items, K, SampleSize, Method) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  complete(Items, K, SampleSize, Method, K_Status, fill = list(n = 0, pct = 0))

build_kstatus_panel <- function(data, items_val) {
  d <- filter(data, Items == items_val)
  
  ggplot(d, aes(x = SampleSize, y = pct, fill = Method)) +
    geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.9) +
    geom_text(aes(label = ifelse(pct > 0, sprintf("%.0f", pct), "")),
              position = position_dodge(0.75), vjust = -0.4, size = 2.6,
              color = "grey25", show.legend = FALSE) +
    facet_grid(rows = vars(K), cols = vars(K_Status),
               labeller = labeller(K = label_parsed)) +
    scale_x_discrete(labels = function(x) paste0("n = ", x)) +
    scale_fill_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Sample size", y = "% of replications", fill = "Method") +
    theme_elegant()
}

banner_k30 <- items_banner(items_levels[1])
banner_k60 <- items_banner(items_levels[2])
p_k30 <- build_kstatus_panel(k_status_summary, items_levels[1])
p_k60 <- build_kstatus_panel(k_status_summary, items_levels[2])

design_k <- "
A
C
H
B
D
"

fig_kstatus <- tryCatch(
  wrap_plots(A = banner_k30, C = p_k30, H = divider_h(), B = banner_k60, D = p_k60,
             design = design_k, heights = c(0.09, 1, 0.06, 0.09, 1),
             guides = "collect") +
    plot_annotation(
      theme = theme(
        plot.title   = element_text(hjust = 0.5, face = "bold", size = 15),
        plot.caption = element_text(color = "grey50", size = 9, hjust = 1)
      )
    ) &
    theme(legend.position = "top"),
  error = function(e) stop("Failed building figure 04_Selected_K_by_K_Status: ", conditionMessage(e))
)

save_figure(fig_kstatus, "04_Selected_K_by_K_Status", width = 13, height = 9.3)

## ---- 10. Final verification -------------------------------------------------

expected_files <- file.path(out_dir, c(
  "01_R_squared_true.pdf", "01_R_squared_true.png",
  "02_ARI_true.pdf", "02_ARI_true.png",
  "03_R_squared_true_ratio.pdf", "03_R_squared_true_ratio.png",
  "04_Selected_K_by_K_Status.pdf", "04_Selected_K_by_K_Status.png",
  "combined_simulation_results.csv", "combined_simulation_results.rds"
))
missing_out <- expected_files[!file.exists(expected_files)]
if (length(missing_out) > 0) {
  stop("Script finished but the following expected output file(s) are missing:\n  ",
       paste(missing_out, collapse = "\n  "))
}

cat("\nAll checks passed. Figures written to:\n ", normalizePath(out_dir), "\n")
