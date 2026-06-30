## =====================================================================
## analysis_dash_first_grid_reorganized.R
##
## Reorganized DASH first-grid analysis.
##
## Main changes from the previous script:
##   1) Keep 7 methods only: DASH Bonf-minP, DASH Cauchy, ZINQ,
##      ANCOM-BC2, DESeq2, edgeR, metagenomeSeq.
##   2) DASH colors fixed to #E73649 and #ee7D09.
##      Comparator colors fixed to #f8bd00, #96c535, #00a496, #0091d8, #d45f9d.
##   3) Comparator methods are drawn with dashed, thinner lines.
##      DASH methods are solid and slightly emphasized.
##   4) Type-I error is exported as a table, not a figure.
##   5) Power and FDR are each exported as one combined figure:
##      unconfounded panels on the left; confounded panels on the right.
##   6) In the power figure only, n = 320 keeps effect sizes up to and including 2.
##
## INPUT:
##   combined_raw.csv
##
## OUTPUT:
##   table_typeI_firstgrid.csv
##   table_typeI_firstgrid.tex
##   fig_power_combined_firstgrid.pdf / .png
##   fig_fdr_combined_firstgrid.pdf / .png
##
## REQUIRES:
##   install.packages(c("ggplot2","dplyr","tidyr","scales"))
## =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

## ---------------------------------------------------------------- settings
infile    <- "combined_raw.csv"
PREFIX    <- "firstgrid"
ALPHA     <- 0.05

## This matches the old red/green rule in your current Type-I plot:
## controlled if FPR <= 0.05 + 0.01; inflated if FPR > 0.06.
## Change to 0 if you want strict nominal 0.05 flagging.
TYPEI_TOL <- 0.01

## In the power figure only, keep effect sizes up to and including this value for n = 320.
N320_POWER_MAX_EFFECT <- 2

## ---------------------------------------------------------------- methods
## Remove "Two-part (obs.)" from the main display because its label is ambiguous
## and it makes the figure crowded. Keep it only for supplement if needed.
focus <- data.frame(
  raw = c("OMNIBUS(Bonf-minP).adj",
          "OMNIBUS(ACAT).adj",
          "pkg.ZINQ.depthCov",
          "pkg.ANCOMBC.cov",
          "pkg.DESeq2.cov",
          "pkg.edgeR.cov",
          "pkg.metagenomeSeq.cov"),
  lab = c("DASH (Bonf-minP)",
          "DASH (Cauchy)",
          "ZINQ",
          "ANCOM-BC2",
          "DESeq2",
          "edgeR",
          "metagenomeSeq"),
  stringsAsFactors = FALSE
)

lab_levels <- focus$lab

pal <- c(
  "DASH (Bonf-minP)" = "#E73649",
  "DASH (Cauchy)"      = "#ee7D09",
  "ZINQ"             = "#f8bd00",
  "ANCOM-BC2"        = "#96c535",
  "DESeq2"           = "#00a496",
  "edgeR"            = "#0091d8",
  "metagenomeSeq"    = "#d45f9d"
)

## Shapes are retained mainly for black-and-white readability.
shp <- c(
  "DASH (Bonf-minP)" = 19,
  "DASH (Cauchy)"      = 17,
  "ZINQ"             = 15,
  "ANCOM-BC2"        = 18,
  "DESeq2"           = 3,
  "edgeR"            = 4,
  "metagenomeSeq"    = 8
)

ltp <- c(
  "DASH (Bonf-minP)" = "solid",
  "DASH (Cauchy)"      = "solid",
  "ZINQ"             = "dashed",
  "ANCOM-BC2"        = "dashed",
  "DESeq2"           = "dashed",
  "edgeR"            = "dashed",
  "metagenomeSeq"    = "dashed"
)

scen_levels <- c("null", "abundance", "prevalence", "mixed", "both")
scen_labels <- c("Null", "Abundance", "Prevalence", "Mixed", "Both")
nn_scen     <- c("abundance", "prevalence", "mixed", "both")

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

pal <- pal[lab_levels]
shp <- shp[lab_levels]
ltp <- ltp[lab_levels]

## ---------------------------------------------------------------- aggregate
se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) NA_real_ else sd(x) / sqrt(length(x))
}
mn <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

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
    n_lab    = factor(paste0("n = ", n), levels = paste0("n = ", sort(unique(n)))),
    conf_lab = factor(ifelse(confounded, "Confounded", "Unconfounded"),
                      levels = c("Unconfounded", "Confounded")),
    is_dash  = grepl("^DASH", as.character(lab))
  )

## ---------------------------------------------------------------- theme + save
theme_dash <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text       = element_text(face = "bold", size = rel(0.82)),
      axis.title       = element_text(face = "bold"),
      axis.text.x      = element_text(size = rel(0.82)),
      axis.text.y      = element_text(size = rel(0.88)),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.key.width = unit(1.35, "lines"),
      legend.key.height= unit(0.8, "lines"),
      plot.title       = element_text(face = "bold", size = rel(1.2)),
      plot.subtitle    = element_text(colour = "grey30"),
      plot.margin      = margin(6, 8, 6, 8)
    )
}

save_fig <- function(p, name, w, h) {
  ggsave(paste0(name, ".pdf"), p, width = w, height = h, useDingbats = FALSE)
  ggsave(paste0(name, ".png"), p, width = w, height = h, dpi = 320, bg = "white")
  message("wrote ", name, ".pdf / .png")
}

## ---------------------------------------------------------------- Type-I table
f3 <- function(x) ifelse(is.finite(x), formatC(x, format = "f", digits = 3), "--")

typeI <- agg %>%
  filter(scenario == "null") %>%
  mutate(
    inflated = is.finite(fpr) & fpr > ALPHA + TYPEI_TOL,
    col_lab  = paste(conf_lab, n_lab, sep = ": ")
  ) %>%
  select(Method = lab, confounded, conf_lab, n, n_lab, nrep, fpr, fpr_se, inflated, col_lab)

## Wide CSV table: numeric values only.
typeI_csv <- typeI %>%
  mutate(value = round(fpr, 4)) %>%
  select(Method, col_lab, value) %>%
  pivot_wider(names_from = col_lab, values_from = value)

typeI_csv$Method <- factor(typeI_csv$Method, levels = lab_levels)
typeI_csv <- typeI_csv[order(typeI_csv$Method), , drop = FALSE]
typeI_csv$Method <- as.character(typeI_csv$Method)

write.csv(typeI_csv, paste0("table_typeI_", PREFIX, ".csv"), row.names = FALSE)

## LaTeX table with grouped headers. Values above ALPHA + TYPEI_TOL are bolded.
n_values <- sort(unique(typeI$n))

cell_typeI <- function(method_lab, conf_val, n_val, latex = TRUE) {
  z <- typeI %>%
    filter(as.character(Method) == method_lab,
           confounded == conf_val,
           n == n_val)
  if (nrow(z) == 0 || !is.finite(z$fpr[1])) return("--")
  out <- f3(z$fpr[1])
  if (latex && isTRUE(z$inflated[1])) out <- paste0("\\textbf{", out, "}")
  out
}

make_typeI_latex <- function() {
  ncol_block <- length(n_values)
  align <- paste0("l", paste(rep("c", 2 * ncol_block), collapse = ""))
  
  lines <- c(
    "\\begin{table}[t]",
    "\\centering",
    sprintf("\\caption{Per-taxon Type I error under the global null. Values are averaged over simulation replicates. Nominal level is %.2f; bold values exceed %.2f.}", ALPHA, ALPHA + TYPEI_TOL),
    sprintf("\\label{tab:typeI_%s}", PREFIX),
    sprintf("\\begin{tabular}{%s}", align),
    "\\toprule",
    paste0(" & \\multicolumn{", ncol_block, "}{c}{Unconfounded}",
           " & \\multicolumn{", ncol_block, "}{c}{Confounded} \\\\"),
    paste0("\\cmidrule(lr){2-", 1 + ncol_block, "} ",
           "\\cmidrule(lr){", 2 + ncol_block, "-", 1 + 2 * ncol_block, "}"),
    paste(c("Method",
            paste0("$n=", n_values, "$"),
            paste0("$n=", n_values, "$")),
          collapse = " & ") |> paste0(" \\\\"),
    "\\midrule"
  )
  
  for (m in lab_levels) {
    row <- c(m,
             vapply(n_values, function(nn) cell_typeI(m, FALSE, nn), character(1)),
             vapply(n_values, function(nn) cell_typeI(m, TRUE,  nn), character(1)))
    lines <- c(lines, paste(row, collapse = " & ") |> paste0(" \\\\"))
  }
  
  c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

writeLines(make_typeI_latex(), paste0("table_typeI_", PREFIX, ".tex"))

cat("\n==== Type-I table ====\n")
print(typeI_csv, row.names = FALSE)
cat("\nWrote table_typeI_", PREFIX, ".csv and .tex\n", sep = "")

## ---------------------------------------------------------------- combined figures
filter_power_tail <- function(d) {
  ## For n = 320, power saturates quickly, so truncate the x-axis at effect size 2.
  ## For other n values, keep the full effect grid.
  d %>% filter(n != 320 | effect <= N320_POWER_MAX_EFFECT)
}

make_combined_curve <- function(d, yvar, yse, ylab, ttl, sub,
                                hline = NA_real_, free_x = FALSE) {
  d <- d %>%
    mutate(
      yy   = .data[[yvar]],
      yyse = .data[[yse]]
    )
  
  d_comp <- d %>% filter(!is_dash)
  d_dash <- d %>% filter(is_dash)
  
  p <- ggplot(d, aes(effect, yy, colour = lab, shape = lab,
                     linetype = lab, group = lab)) +
    { if (is.finite(hline)) geom_hline(yintercept = hline,
                                       linetype = "dashed",
                                       colour = "grey35",
                                       linewidth = 0.35) } +
    geom_errorbar(aes(ymin = pmax(0, yy - yyse),
                      ymax = pmin(1, yy + yyse)),
                  width = 0,
                  alpha = 0.18,
                  linewidth = 0.25,
                  show.legend = FALSE) +
    geom_line(data = d_comp, linewidth = 0.34, alpha = 0.95) +
    geom_point(data = d_comp, size = 1.25, stroke = 0.55) +
    geom_line(data = d_dash, linewidth = 0.62, alpha = 0.98) +
    geom_point(data = d_dash, size = 1.85, stroke = 0.55) +
    facet_grid(
      rows = vars(scen_lab),
      cols = vars(conf_lab, n_lab),
      scales = if (free_x) "free_x" else "fixed"
    ) +
    scale_colour_manual(values = pal, breaks = lab_levels, drop = FALSE) +
    scale_shape_manual(values = shp, breaks = lab_levels, drop = FALSE) +
    scale_linetype_manual(values = ltp, breaks = lab_levels, drop = FALSE) +
    scale_y_continuous(labels = label_percent(accuracy = 1),
                       limits = c(0, 1),
                       breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0.02, 0.05))) +
    scale_x_continuous(breaks = breaks_pretty(n = 4)) +
    labs(x = "Effect size", y = ylab) +
    theme_dash() +
    guides(
      colour   = guide_legend(nrow = 2, byrow = TRUE),
      shape    = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 2, byrow = TRUE)
    )
  
  p
}

d_nn <- agg %>%
  filter(scenario %in% nn_scen) %>%
  mutate(
    scen_lab = factor(scenario, levels = nn_scen,
                      labels = scen_labels[match(nn_scen, scen_levels)])
  )

## Figure 1: combined power
d_power <- filter_power_tail(d_nn)

## Hard diagnostic: the power figure must not contain effect sizes above 2 when n = 320.
## This prevents accidentally regenerating the old "largest two effect sizes" figure.
if (any(d_power$n == 320 & d_power$effect > N320_POWER_MAX_EFFECT + 1e-10, na.rm = TRUE)) {
  stop("Power filtering failed: n = 320 still contains effect sizes above ",
       N320_POWER_MAX_EFFECT, ".")
}

cat("\n==== effect sizes used in the POWER figure ====\n")
print(
  d_power %>%
    group_by(conf_lab, n_lab) %>%
    summarise(effect_sizes = paste(sort(unique(effect)), collapse = ", "),
              .groups = "drop"),
  n = Inf
)

p_power <- make_combined_curve(
  d_power,
  yvar = "power",
  yse  = "power_se",
  ylab = "BH power",
  ttl  = "Detection power",
  sub  = "Unconfounded panels are on the left and confounded panels are on the right; for n = 320, effect sizes are shown only up to 2.",
  hline = NA_real_,
  free_x = TRUE
)

save_fig(p_power, paste0("fig_power_combined_", PREFIX), 15.5, 8.8)

## Figure 2: combined FDR
p_fdr <- make_combined_curve(
  d_nn,
  yvar = "fdr",
  yse  = "fdr_se",
  ylab = "Realized FDR",
  ttl  = "False discovery rate control",
  sub  = "Unconfounded panels are on the left and confounded panels are on the right; dashed horizontal line = BH target 0.05.",
  hline = ALPHA,
  free_x = FALSE
)

save_fig(p_fdr, paste0("fig_fdr_combined_", PREFIX), 15.5, 8.8)

cat("Done.\n")
