## =====================================================================
## analysis_dash_sparsedossa2.R
##
## Local analysis of the DASH SparseDOSSA2-grid results.
##
## INPUT (in the working directory):
##   combined_raw.csv     <- set `infile` below if your file is named differently
##
## OUTPUT (written to the working directory):
##   fig1_typeI_<PREFIX>.{pdf,png}        per-taxon Type I error under the null
##   fig2_fdr_<PREFIX>_conf<TF>.{pdf,png} realized FDR vs effect (BH, target 0.05)
##   fig3_power_<PREFIX>_conf<TF>.{pdf,png} BH power vs effect
##   fig4_tradeoff_<PREFIX>.{pdf,png}     power vs FDR frontier
##   table_summary_<PREFIX>.csv           summary table (readable)
##   table_summary_<PREFIX>.tex           summary table (LaTeX, booktabs)
##
## REQUIRES (install once):
##   install.packages(c("ggplot2","dplyr","tidyr","scales"))
## =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

## ---------------------------------------------------------------- settings
infile    <- "combined_raw.csv"   # rename here if needed
PREFIX    <- "sparsedossa2"       # output filename tag
CONF_MAIN <- TRUE                 # confounded value used for table + tradeoff
N_TABLE   <- 160                  # sample size used in the summary table
ALPHA     <- 0.05                 # nominal level / BH target

## Focus methods: raw name in the CSV -> pretty label. DASH listed first.
focus <- data.frame(
  raw = c("DASH.BonfminP.adj", "DASH.ACAT.adj", "pkg.ZINQ.depthCov",
          "pkg.ANCOMBC.cov", "pkg.DESeq2.cov", "pkg.edgeR.cov",
          "pkg.metagenomeSeq.cov", "obsTwo.depthCov"),
  lab = c("DASH (Bonf-minP)", "DASH (ACAT)", "ZINQ", "ANCOM-BC2",
          "DESeq2", "edgeR", "metagenomeSeq", "Two-part (obs.)"),
  stringsAsFactors = FALSE
)

## Colour + shape per pretty label (category-10 palette; DASH in warm colours).
pal <- c(
  "DASH (Bonf-minP)" = "#D62728", "DASH (ACAT)" = "#FF7F0E",
  "ZINQ" = "#1F77B4", "ANCOM-BC2" = "#2CA02C", "DESeq2" = "#9467BD",
  "edgeR" = "#8C564B", "metagenomeSeq" = "#E377C2", "Two-part (obs.)" = "#7F7F7F"
)
shp <- c(
  "DASH (Bonf-minP)" = 19, "DASH (ACAT)" = 17, "ZINQ" = 15, "ANCOM-BC2" = 18,
  "DESeq2" = 3, "edgeR" = 4, "metagenomeSeq" = 8, "Two-part (obs.)" = 7
)

scen_levels <- c("null", "abundance", "prevalence", "mixed", "both")
scen_labels <- c("Null", "Abundance", "Prevalence", "Mixed", "Both")
nn_scen     <- c("abundance", "prevalence", "mixed", "both")   # non-null

## ---------------------------------------------------------------- read data
if (!file.exists(infile)) stop("Cannot find '", infile, "' in: ", getwd())
df <- read.csv(infile, stringsAsFactors = FALSE)

need <- c("rep", "scenario", "n", "effect", "confounded", "method",
          "raw_fpr", "bh_power", "bh_fdr")
miss <- setdiff(need, colnames(df))
if (length(miss)) stop("Input is missing columns: ", paste(miss, collapse = ", "))

df$confounded <- as.logical(df$confounded)
df$n          <- as.integer(df$n)
df$effect     <- as.numeric(df$effect)
df$scenario   <- factor(df$scenario, levels = scen_levels)

present <- intersect(focus$raw, unique(df$method))
dropped <- setdiff(focus$raw, present)
if (length(dropped)) {
  message("Note: focus methods absent / all-NA in the data, dropped: ",
          paste(focus$lab[match(dropped, focus$raw)], collapse = ", "))
}
focus <- focus[focus$raw %in% present, , drop = FALSE]
lab_levels <- focus$lab

## ---------------------------------------------------------------- aggregate
se <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x)) }
mn <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else mean(x) }

agg <- df %>%
  filter(method %in% focus$raw) %>%
  group_by(scenario, n, effect, confounded, method) %>%
  summarise(
    nrep     = dplyr::n(),
    power    = mn(bh_power),  power_se = se(bh_power),
    fdr      = mn(bh_fdr),    fdr_se   = se(bh_fdr),
    fpr      = mn(raw_fpr),   fpr_se   = se(raw_fpr),
    .groups  = "drop"
  ) %>%
  left_join(focus, by = c("method" = "raw")) %>%
  mutate(
    lab      = factor(lab, levels = lab_levels),
    scen_lab = factor(scenario, levels = scen_levels, labels = scen_labels),
    n_lab    = factor(paste0("n = ", n), levels = paste0("n = ", sort(unique(n))))
  )

## ---------------------------------------------------------------- theme + IO
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = rel(0.9)),
      axis.title       = element_text(face = "bold"),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.key.width = unit(1.1, "lines"),
      plot.title       = element_text(face = "bold"),
      plot.subtitle    = element_text(colour = "grey30")
    )
}

save_fig <- function(p, name, w, h) {
  ggsave(paste0(name, ".pdf"), p, width = w, height = h)
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = 320, bg = "white")
  message("wrote ", name, ".pdf / .png")
}

## ============================================================ FIGURE 1: Type I
t1e <- agg %>%
  filter(scenario == "null", is.finite(fpr)) %>%
  mutate(
    conf_lab = factor(ifelse(confounded, "Confounded", "No confounder"),
                      levels = c("No confounder", "Confounded")),
    status   = ifelse(fpr <= ALPHA + 0.01, "Controlled", "Inflated")
  )

if (nrow(t1e)) {
  p1 <- ggplot(t1e, aes(fpr, lab)) +
    geom_vline(xintercept = ALPHA, linetype = "dashed", colour = "grey40") +
    geom_errorbarh(aes(xmin = pmax(0, fpr - fpr_se), xmax = fpr + fpr_se),
                   height = 0.25, colour = "grey55") +
    geom_point(aes(colour = status), size = 2.6) +
    scale_colour_manual(values = c(Controlled = "#2CA02C", Inflated = "#D62728")) +
    scale_x_continuous(labels = label_number(accuracy = 0.01)) +
    facet_grid(conf_lab ~ n_lab) +
    labs(title = "Type I error under the global null",
         subtitle = "Per-taxon false positive rate; dashed line = nominal 0.05",
         x = "Per-taxon false positive rate", y = NULL, colour = NULL) +
    theme_pub() +
    scale_y_discrete(limits = rev(lab_levels))
  save_fig(p1, paste0("fig1_typeI_", PREFIX), 9, 6)
}

## ====================================================== FIGURES 2/3: FDR & power
make_curve <- function(d, yvar, yse, ylab, ttl, sub, hline = NA) {
  d <- d %>% rename(yy = !!yvar, yyse = !!yse)
  dD <- d %>% filter(grepl("DASH", lab))
  p <- ggplot(d, aes(effect, yy, colour = lab, shape = lab, group = lab)) +
    { if (is.finite(hline)) geom_hline(yintercept = hline, linetype = "dashed", colour = "grey40") } +
    geom_errorbar(aes(ymin = yy - yyse, ymax = yy + yyse),
                  width = 0, alpha = 0.35, linewidth = 0.4) +
    geom_line(linewidth = 0.65) +
    geom_point(size = 1.8) +
    geom_line(data = dD, linewidth = 1.2) +
    geom_point(data = dD, size = 2.7) +
    facet_grid(scen_lab ~ n_lab) +
    scale_colour_manual(values = pal, breaks = lab_levels) +
    scale_shape_manual(values = shp, breaks = lab_levels) +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    labs(title = ttl, subtitle = sub, x = "Effect size", y = ylab) +
    theme_pub()
  p
}

for (cf in c(FALSE, TRUE)) {
  tag <- if (cf) "confT" else "confF"
  sub <- if (cf) "With a measured confounder (covariate-adjusted methods)"
  else    "No measured confounder (covariate-adjusted methods)"
  
  d_nn <- agg %>% filter(confounded == cf, scenario %in% nn_scen) %>%
    mutate(scen_lab = factor(scenario, levels = nn_scen,
                             labels = scen_labels[match(nn_scen, scen_levels)]))
  
  if (nrow(d_nn)) {
    p3 <- make_curve(d_nn, "power", "power_se", "BH power",
                     "Detection power", sub)
    save_fig(p3, paste0("fig3_power_", PREFIX, "_", tag), 9, 9)
    
    p2 <- make_curve(d_nn, "fdr", "fdr_se", "Realized FDR",
                     "False discovery rate control", sub, hline = ALPHA)
    save_fig(p2, paste0("fig2_fdr_", PREFIX, "_", tag), 9, 9)
  }
}

## ===================================================== FIGURE 4: power vs FDR
trade <- agg %>%
  filter(confounded == CONF_MAIN, scenario %in% nn_scen,
         is.finite(power), is.finite(fdr)) %>%
  mutate(scen_lab = factor(scenario, levels = nn_scen,
                           labels = scen_labels[match(nn_scen, scen_levels)]))

if (nrow(trade)) {
  p4 <- ggplot(trade, aes(fdr, power, colour = lab, shape = lab, group = lab)) +
    geom_vline(xintercept = ALPHA, linetype = "dashed", colour = "grey40") +
    geom_path(linewidth = 0.5, alpha = 0.7) +
    geom_point(size = 2) +
    facet_grid(scen_lab ~ n_lab) +
    scale_colour_manual(values = pal, breaks = lab_levels) +
    scale_shape_manual(values = shp, breaks = lab_levels) +
    scale_x_continuous(labels = label_percent(accuracy = 1)) +
    scale_y_continuous(labels = label_percent(accuracy = 1)) +
    labs(title = "Power versus realized FDR",
         subtitle = paste0("Each path traces the effect grid; vertical line = ",
                           "BH target 0.05 (", ifelse(CONF_MAIN, "confounded", "unconfounded"), ")"),
         x = "Realized FDR", y = "BH power") +
    theme_pub()
  save_fig(p4, paste0("fig4_tradeoff_", PREFIX), 9, 8)
}

## ============================================================ SUMMARY TABLE
f2 <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 2), "--")
f3 <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 3), "--")

tab_src <- agg %>% filter(n == N_TABLE, confounded == CONF_MAIN)

t1col <- tab_src %>% filter(scenario == "null") %>%
  transmute(lab, typeI = fpr)

nn <- tab_src %>% filter(scenario %in% nn_scen) %>%
  group_by(lab, scenario) %>%
  summarise(power = mn(power), fdr = mn(fdr), .groups = "drop") %>%
  mutate(cell = paste0(f2(power), " / ", f2(fdr))) %>%
  select(lab, scenario, cell) %>%
  pivot_wider(names_from = scenario, values_from = cell)

tab <- t1col %>% left_join(nn, by = "lab")
tab$lab <- factor(tab$lab, levels = lab_levels)
tab <- tab[order(tab$lab), , drop = FALSE]

scen_present <- intersect(nn_scen, colnames(tab))
disp <- data.frame(
  Method = as.character(tab$lab),
  `Type I` = f3(tab$typeI),
  check.names = FALSE, stringsAsFactors = FALSE
)
for (s in scen_present) disp[[scen_labels[match(s, scen_levels)]]] <- tab[[s]]

write.csv(disp, paste0("table_summary_", PREFIX, ".csv"), row.names = FALSE)

## LaTeX (booktabs); non-null cells are "Power / FDR".
make_latex <- function(disp, caption, label, bold_row = NA) {
  header <- colnames(disp)
  align  <- paste0("l", paste(rep("c", ncol(disp) - 1), collapse = ""))
  lines <- c(
    "\\begin{table}[t]", "\\centering",
    sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
    sprintf("\\begin{tabular}{%s}", align), "\\toprule",
    paste(paste(header, collapse = " & "), "\\\\"), "\\midrule"
  )
  for (i in seq_len(nrow(disp))) {
    row <- as.character(unlist(disp[i, ]))
    if (!is.na(bold_row) && i == bold_row) row <- paste0("\\textbf{", row, "}")
    lines <- c(lines, paste(paste(row, collapse = " & "), "\\\\"))
  }
  c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

cap <- sprintf(paste0("Type I error (null) and detection power / realized FDR by scenario ",
                      "at $n=%d$%s, averaged over the effect grid and replicates. ",
                      "Non-null cells report Power / FDR; BH target FDR is %.2f."),
               N_TABLE, ifelse(CONF_MAIN, " with a measured confounder", ""), ALPHA)
bold <- which(disp$Method == "DASH (Bonf-minP)")
writeLines(make_latex(disp, cap, paste0("tab:", PREFIX),
                      bold_row = if (length(bold)) bold else NA),
           paste0("table_summary_", PREFIX, ".tex"))

cat("\n==== summary table (n =", N_TABLE,
    ", confounded =", CONF_MAIN, ") ====\n")
print(disp, row.names = FALSE)
cat("\nWrote table_summary_", PREFIX, ".csv and .tex\n", sep = "")
cat("Done.\n")