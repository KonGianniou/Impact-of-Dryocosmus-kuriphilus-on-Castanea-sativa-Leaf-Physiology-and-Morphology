# ==============================================================================
# Impact of Dryocosmus kuriphilus on Castanea sativa
# Physiological and morphological leaf trait analysis — 5 provenances, 3 countries
# University Forest of Taxiarchis, Halkidiki, Greece (760 m a.s.l.)
# ==============================================================================
# Author:   Konstantina Gianniou
# Contact:  g.tem2106@gmail.com | github.com/KonGianniou
# DISCLAIMER: This script runs on synthetic data. Real measurements remain the
#             property of the University of the Aegean / ELGO Demeter/ Aristotle
#             University of Thessaloniki.
#             Code rewritten from scratch to avoid copyright issues.
# Reference: DCI formula — Gehring et al. (2018), J. Vis. Exp., 138.
# ==============================================================================

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────

pkgs <- c("ggplot2", "dplyr", "tidyr", "car", "gplots",
          "FactoMineR", "factoextra", "GGally", "lme4", "lmerTest", "patchwork")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(car)
  library(gplots);  library(FactoMineR); library(factoextra)
  library(GGally);  library(lme4); library(lmerTest); library(patchwork)
})

# ── SHARED THEME & PALETTE ────────────────────────────────────────────────────

theme_study <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10, colour = "grey45"),
    axis.title       = element_text(size = 11),
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

prov_pal <- c(
  Coruna    = "#1d6fa4",
  Hortiatis = "#2ca25f",
  Malaga    = "#e6550d",
  Pellice   = "#9467bd",
  Sicily    = "#d62728"
)

country_pal <- c(Spain = "#1d6fa4", Greece = "#2ca25f", Italy = "#e6550d")

# Helper: print a clean section header
hdr <- function(txt) cat(sprintf("\n%s\n%s\n", txt, strrep("─", nchar(txt))))

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

source("generate_synthetic_data.R")
dci        <- dci_synth
physiologia <- phys_synth

# Clean complete-case subset for leaf physiology analyses
phys       <- physiologia[complete.cases(physiologia[, c("PI","Fv.Fm","CCI",
                          "Gs","Trm","A_sat","F_W","T_W","D_W","L_Area","LDMC","LMA")]), ]

hdr("1. Data loaded")
cat(sprintf("Branch data: %d trees across %d provenances\n",
            nrow(dci), length(unique(dci$Variety))))
cat(sprintf("Leaf data:   %d leaf records (%d with complete physiology)\n",
            nrow(physiologia), nrow(phys)))

# ── 2. DCI INFESTATION INDEX ──────────────────────────────────────────────────

hdr("2. DCI Infestation Index")

# Summary by provenance and country
cat("DCI by provenance:\n")
print(dci %>% group_by(Variety) %>%
  summarise(n=n(), Mean=round(mean(DCI),2), SD=round(sd(DCI),2), .groups="drop"))

# ── ANOVA DCI ~ Provenance ──
aov_dci_prov <- aov(DCI ~ Variety, data = dci)
cat(sprintf("\nANOVA DCI ~ Variety: F=%.3f, p=%.3f\n",
  summary(aov_dci_prov)[[1]]$`F value`[1],
  summary(aov_dci_prov)[[1]]$`Pr(>F)`[1]))
cat(sprintf("Levene's test: p=%.3f\n",
  leveneTest(DCI ~ Variety, data = dci)$`Pr(>F)`[1]))
res_dci <- residuals(aov_dci_prov)
cat(sprintf("Shapiro-Wilk residuals: p=%.4f\n", shapiro.test(res_dci)$p.value))
cat(sprintf("Kruskal-Wallis: p=%.3f\n", kruskal.test(DCI ~ Variety, data=dci)$p.value))

# ── ANOVA DCI ~ Country ──
aov_dci_ctry <- aov(DCI ~ Country, data = dci)
cat(sprintf("\nANOVA DCI ~ Country: F=%.3f, p=%.3f\n",
  summary(aov_dci_ctry)[[1]]$`F value`[1],
  summary(aov_dci_ctry)[[1]]$`Pr(>F)`[1]))

# Fig 1: DCI by provenance — violin + boxplot
fig_dci <- ggplot(dci, aes(x = reorder(Variety, DCI, median), y = DCI, fill = Variety)) +
  geom_violin(alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.25, alpha = 0.85, outlier.shape = 21, outlier.size = 1.5) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 4,
               fill = "white", colour = "black") +
  scale_fill_manual(values = prov_pal) +
  labs(title = "Infestation index (DCI) by provenance",
       subtitle = "White dot = mean | ANOVA p > 0.05",
       x = NULL, y = "DCI") +
  theme_study + theme(legend.position = "none")

ggsave("fig1_DCI_provenance.png", fig_dci, width = 7, height = 4.5, dpi = 150)

# ── 3. REACTIVATED DORMANT BUDS (RDB), DEAD SHOOTS, GALLS ────────────────────

hdr("3. RDB, Dead Shoots, Galls")

# Function: run ANOVA + assumption checks + KW for one variable
run_anova_kw <- function(var, group, data, label) {
  f <- as.formula(paste(var, "~", group))
  aov_fit  <- aov(f, data = data)
  lev_p    <- leveneTest(f, data = data)$`Pr(>F)`[1]
  sw_p     <- shapiro.test(residuals(aov_fit))$p.value
  kw_p     <- kruskal.test(f, data = data)$p.value
  aov_p    <- summary(aov_fit)[[1]]$`Pr(>F)`[1]
  cat(sprintf("  %-20s | ANOVA p=%.3f | Levene p=%.3f | SW p=%.4f | KW p=%.3f\n",
              label, aov_p, lev_p, sw_p, kw_p))
}

cat("RDB:\n")
run_anova_kw("RDB", "Variety", dci, "by Provenance")
run_anova_kw("RDB", "Country", dci, "by Country")

cat("Dead Shoots:\n")
run_anova_kw("Dead_Sh", "Variety", dci, "by Provenance")
run_anova_kw("Dead_Sh", "Country", dci, "by Country")

cat("Galls:\n")
run_anova_kw("Galls", "Variety", dci, "by Provenance")
run_anova_kw("Galls", "Country", dci, "by Country")

# Fig 2: RDB by provenance
fig_rdb <- ggplot(dci, aes(x = Variety, y = RDB, fill = Variety)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 1.5) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3.5,
               fill = "#1f77b4", colour = "black") +
  scale_fill_manual(values = prov_pal) +
  labs(title = "Reactivated Dormant Buds (RDB) by provenance",
       subtitle = "Blue dot = mean | KW p = 0.49 (ns)",
       x = NULL, y = "RDB count") +
  theme_study + theme(legend.position = "none")

ggsave("fig2_RDB_provenance.png", fig_rdb, width = 7, height = 4.5, dpi = 150)

# ── 4. CHLOROPHYLL CONTENT INDEX (CCI) ───────────────────────────────────────

hdr("4. Chlorophyll Content Index (CCI)")

cat("CCI by provenance:\n")
print(physiologia %>% group_by(Variety) %>%
  summarise(Mean=round(mean(CCI),2), SD=round(sd(CCI),2), .groups="drop"))

# ANOVA CCI ~ Provenance
aov_cci <- aov(CCI ~ Variety, data = physiologia)
cci_p   <- summary(aov_cci)[[1]]$`Pr(>F)`[1]
cat(sprintf("\nANOVA CCI ~ Variety: p=%.4f %s\n",
            cci_p, ifelse(cci_p < 0.05, "✓ Significant", "(ns)")))
cat(sprintf("Levene's test: p=%.3f\n",
            leveneTest(CCI ~ Variety, data=physiologia)$`Pr(>F)`[1]))
cat(sprintf("Shapiro-Wilk residuals: p=%.3f\n",
            shapiro.test(residuals(aov_cci))$p.value))

if (cci_p < 0.05) {
  cat("\nTukey HSD post-hoc (significant pairs):\n")
  tukey <- TukeyHSD(aov_cci)$Variety
  sig   <- tukey[tukey[, "p adj"] < 0.05, , drop = FALSE]
  print(round(sig, 4))
}

# ANOVA CCI ~ Country
aov_cci_c <- aov(CCI ~ Country, data = physiologia)
cat(sprintf("\nANOVA CCI ~ Country: p=%.4f\n",
            summary(aov_cci_c)[[1]]$`Pr(>F)`[1]))

# Fig 3: CCI violin by provenance + country facet
fig_cci_prov <- ggplot(physiologia, aes(x=Variety, y=CCI, fill=Variety)) +
  geom_violin(alpha = 0.45, colour = NA) +
  geom_boxplot(width = 0.22, alpha = 0.85, outlier.shape = 21, outlier.size = 1.5) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3.5,
               fill = "white", colour = "black") +
  scale_fill_manual(values = prov_pal) +
  labs(title = "Chlorophyll Content Index (CCI) by provenance",
       subtitle = sprintf("ANOVA p = %.4f", cci_p),
       x = NULL, y = "CCI") +
  theme_study + theme(legend.position = "none")

fig_cci_ctry <- ggplot(physiologia, aes(x=Country, y=CCI, fill=Country)) +
  geom_violin(alpha = 0.45, colour = NA) +
  geom_boxplot(width = 0.22, alpha = 0.85) +
  scale_fill_manual(values = country_pal) +
  labs(title = "CCI by country of origin", x = NULL, y = "CCI") +
  theme_study + theme(legend.position = "none")

ggsave("fig3_CCI.png", fig_cci_prov + fig_cci_ctry, width=12, height=5, dpi=150)

# ── 5. PHOTOSYNTHESIS (A_sat) ─────────────────────────────────────────────────

hdr("5. Photosynthesis (A_sat)")

A <- physiologia[!is.na(physiologia$A_sat), ]
cat("A_sat by provenance:\n")
print(A %>% group_by(Variety) %>%
  summarise(Mean=round(mean(A_sat),2), SD=round(sd(A_sat),2), .groups="drop"))

aov_asat  <- aov(A_sat ~ Variety, data = A)
asat_p    <- summary(aov_asat)[[1]]$`Pr(>F)`[1]
cat(sprintf("\nANOVA A_sat ~ Variety: p=%.4f %s\n",
            asat_p, ifelse(asat_p<0.05,"✓ Significant","(ns)")))
cat(sprintf("Welch ANOVA: p=%.4f\n", oneway.test(A_sat~Variety, A)$p.value))
cat(sprintf("Kruskal-Wallis: p=%.4f\n", kruskal.test(A_sat~Variety, A)$p.value))

if (asat_p < 0.05) {
  cat("\nTukey HSD post-hoc (significant pairs):\n")
  tukey_a <- TukeyHSD(aov_asat)$Variety
  print(round(tukey_a[tukey_a[,"p adj"]<0.05, , drop=FALSE], 4))
}

# ── 6. OTHER LEAF TRAITS (Fv/Fm, LMA, LDMC, Weights, Leaf Area) ──────────────

hdr("6. Other leaf traits — ANOVA + KW by provenance")

leaf_vars <- list(
  list(var="Fv.Fm",  label="Fv/Fm fluorescence",   subset_na=TRUE),
  list(var="LMA",    label="LMA",                   subset_na=TRUE),
  list(var="LDMC",   label="LDMC",                  subset_na=TRUE),
  list(var="F_W",    label="Fresh weight (F_W)",    subset_na=TRUE),
  list(var="D_W",    label="Dry weight (D_W)",      subset_na=TRUE),
  list(var="T_W",    label="Saturated weight (T_W)",subset_na=TRUE),
  list(var="L_Area", label="Leaf area",             subset_na=TRUE)
)

for (lv in leaf_vars) {
  df_sub <- physiologia[!is.na(physiologia[[lv$var]]), ]
  run_anova_kw(lv$var, "Variety", df_sub, lv$label)
}

# ── 7. SPEARMAN CORRELATIONS: DCI & RDB with Leaf Traits ─────────────────────

hdr("7. Spearman Correlations — DCI & RDB with Leaf Traits")

cor_vars <- c("Fv.Fm","CCI","A_sat","F_W","T_W","D_W","L_Area","LDMC","LMA")

cat(sprintf("%-12s  %8s  %6s    %8s  %6s\n",
            "Variable", "DCI rho", "p", "RDB rho", "p"))
cat(strrep("-", 55), "\n")

for (v in cor_vars) {
  df_v <- physiologia[!is.na(physiologia[[v]]), ]
  c1 <- cor.test(df_v$DCI, df_v[[v]], method = "spearman", exact = FALSE)
  c2 <- cor.test(df_v$RDB, df_v[[v]], method = "spearman", exact = FALSE)
  cat(sprintf("%-12s  %8.3f  %6.4f    %8.3f  %6.4f  %s\n",
              v, c1$estimate, c1$p.value, c2$estimate, c2$p.value,
              ifelse(c1$p.value < 0.05 | c2$p.value < 0.05, "✓", "")))
}

# Extra: RDB ~ DCI (strong correlation expected, rho~0.75)
c_rdb_dci <- cor.test(dci$DCI, dci$RDB, method = "spearman", exact = FALSE)
cat(sprintf("\nRDB ~ DCI:  rho=%.3f, p=%.4f\n",
            c_rdb_dci$estimate, c_rdb_dci$p.value))

# ── 8. PCA ────────────────────────────────────────────────────────────────────

hdr("8. Principal Component Analysis")

# Aggregate to tree level (mean of 2 leaves per tree)
dat_pca <- physiologia %>%
  group_by(ID, Variety) %>%
  summarise(
    PI    = mean(PI,    na.rm = TRUE),
    Fv.Fm = mean(Fv.Fm, na.rm = TRUE),
    CCI   = mean(CCI,   na.rm = TRUE),
    Gs    = mean(Gs,    na.rm = TRUE),
    Trm   = mean(Trm,   na.rm = TRUE),
    A_sat = mean(A_sat, na.rm = TRUE),
    LMA   = mean(LMA,   na.rm = TRUE) * 1000,  # scale up for readability
    LDMC  = mean(LDMC,  na.rm = TRUE),
    DCI   = mean(DCI,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  na.omit()

cat(sprintf("PCA dataset: %d trees\n", nrow(dat_pca)))

pca1 <- PCA(dat_pca[, 3:11], graph = FALSE)
cat("Variance explained (first 3 PCs):\n")
print(round(pca1$eig[1:3, ], 3))

# Full biplot
fig_pca1 <- fviz_pca_biplot(
  pca1,
  habillage  = dat_pca$Variety,
  col.var    = "black",
  alpha.var  = "cos2",
  label      = "var",
  palette    = prov_pal,
  title      = "PCA biplot — all physiological traits + DCI",
  ggtheme    = theme_study
)
ggsave("fig4_PCA_full.png", fig_pca1, width = 7, height = 5.5, dpi = 150)

# Reduced PCA (PI, CCI, A_sat, LMA, LDMC, DCI)
dat_pca2 <- select(dat_pca, Variety, PI, CCI, A_sat, LMA, LDMC, DCI)
pca2     <- PCA(dat_pca2[, 2:7], graph = FALSE)

fig_pca2 <- fviz_pca_biplot(
  pca2,
  habillage  = dat_pca2$Variety,
  col.var    = "black",
  alpha.var  = "cos2",
  label      = "var",
  palette    = prov_pal,
  title      = "PCA biplot — selected traits (PI, CCI, A_sat, LMA, LDMC, DCI)",
  ggtheme    = theme_study
)
ggsave("fig5_PCA_reduced.png", fig_pca2, width = 7, height = 5.5, dpi = 150)

# ── 9. MIXED-EFFECTS MODELS ───────────────────────────────────────────────────

hdr("9. Linear Mixed-Effects Models (tree as random effect)")

# Each tree contributes 2 leaf measurements → ID as random intercept
me_vars <- c("PI", "CCI", "Fv.Fm", "A_sat")

for (v in me_vars) {
  df_me <- physiologia[!is.na(physiologia[[v]]), ]
  m <- tryCatch(
    lmer(as.formula(paste(v, "~ Variety + DCI + LMA + (1|ID)")),
         data = df_me, REML = TRUE),
    error = function(e) NULL
  )
  if (!is.null(m)) {
    cat(sprintf("\n%s ~ Variety + DCI + LMA + (1|ID):\n", v))
    sm <- summary(m)$coefficients
    print(round(sm[, c("Estimate","Std. Error","Pr(>|t|)")], 4))
  }
}

# ── 10. CORRELATION MATRIX VISUALISATION ──────────────────────────────────────

hdr("10. Correlation matrix plot")

cor_mat <- physiologia[, c("Fv.Fm","CCI","DCI","RDB","F_W","D_W",
                            "T_W","L_Area","LDMC","LMA")]
cor_mat <- cor_mat[complete.cases(cor_mat), ]

fig_corr <- ggcorr(
  cor_mat,
  method      = c("complete.obs","spearman"),
  label       = TRUE,
  label_round = 2,
  hjust       = 0.75,
  size        = 3,
  low         = "#e6f1fb",
  mid         = "white",
  high        = "#1d6fa4",
  layout.exp  = 1
) + labs(title = "Spearman correlation matrix — leaf traits, DCI & RDB") +
  theme(plot.title = element_text(face = "bold", size = 12))

ggsave("fig6_correlation_matrix.png", fig_corr, width = 7.5, height = 6.5, dpi = 150)

# ── 11. DCI VALIDATION (Gehring formula reconstruction) ──────────────────────

hdr("11. DCI formula validation")

dci$DCI_check <- (dci$Sd * 0.479 + dci$Bdor * 0.525 + dci$Gons * 0.120) * 100
dci$DCI_diff  <- abs(dci$DCI - dci$DCI_check)
cat(sprintf("Max deviation between stored and recalculated DCI: %.6f\n",
            max(dci$DCI_diff, na.rm = TRUE)))

DCI_lm <- lm(DCI ~ Sd + Bdor + Gons, data = dci)
cat(sprintf("R² of DCI ~ Sd + Bdor + Gons: %.4f (should be ~1)\n",
            summary(DCI_lm)$r.squared))

cat("\n=== Analysis complete ===\n")
cat("Saved figures: fig1_DCI_provenance.png, fig2_RDB_provenance.png,\n")
cat("               fig3_CCI.png, fig4_PCA_full.png, fig5_PCA_reduced.png,\n")
cat("               fig6_correlation_matrix.png\n")
