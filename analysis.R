# ============================================================
# DISCLAIMER: This script generates synthetic data that preserves the statistical
# structure of the original field study. Real measurements are NOT included — they remain
# the intellectual property of the University of the Aegean / ELGO Demeter/
# Aristotle University of Thessaloniki.
# The code has been rewritten from scratch to avoid copyright issues.
#
# Original study: Gianniou K. (2021), Internship Report, Dept. of Environment,
#                 University of the Aegean.
# ============================================================

### Packages ###

install.packages(c("car", "ggplot2", "gplots", "tidyverse", "FactoMineR","factoextra", "lme4", "lmerTest", "sjPlot", "GGally", "wordcloud2"))

library(car)
library(ggplot2)
library(gplots)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(lme4)
library(sjPlot)
library(lmerTest)
library(GGally)
library(wordcloud2)

### Data ###

physiologia=read_csv("Analisi_Physiologia_synthetic.csv")

phys=na.omit(physiologia)
phys_NA=physiologia[complete.cases(physiologia), ]

dci=read_csv("Analisi_Sfika_synthetic.csv")

### Functions ###

summarize_by_group=function(data, response, group) {
  data %>%
    group_by(.data[[group]]) %>%
    summarise(
      mean   = mean(.data[[response]], na.rm = TRUE),
      median = median(.data[[response]], na.rm = TRUE),
      sd     = sd(.data[[response]], na.rm = TRUE),
      n      = sum(!is.na(.data[[response]])),
      .groups = "drop"
    )
}

plot_distribution=function(data, response, group, title = NULL) {
  ggplot(data, aes(x = .data[[group]], y = .data[[response]], fill = .data[[group]])) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.15, alpha = 0.85, outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +
    scale_fill_brewer(palette = "PuBu") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(
      title = title %||% paste(response, "by", group),
      x = group, y = response
    )
}

analyze_group_differences=function(data, response, group, alpha = 0.05,
                                      posthoc_if_significant = TRUE) {
  fmla  <- as.formula(paste(response, "~", group))
  model <- aov(fmla, data = data)
  
  aov_p     <- summary(model)[[1]][["Pr(>F)"]][1]
  levene_p  <- car::leveneTest(fmla, data = data)[["Pr(>F)"]][1]
  shapiro_p <- tryCatch(shapiro.test(residuals(model))$p.value, error = function(e) NA)
  kw_p      <- kruskal.test(fmla, data = data)$p.value
  
  assumptions_ok <- !is.na(shapiro_p) && shapiro_p > alpha && levene_p > alpha
  recommendation <- if (assumptions_ok) {
    sprintf("ANOVA assumptions met -> trust ANOVA (p = %.3f)", aov_p)
  } else {
    sprintf("Assumptions violated -> trust Kruskal-Wallis (p = %.3f)", kw_p)
  }
  
  tukey <- NULL
  if (posthoc_if_significant && aov_p < alpha) {
    tukey <- TukeyHSD(model)
  }
  
  result <- list(
    response = response, group = group, model = model,
    aov_p = aov_p, levene_p = levene_p, shapiro_p = shapiro_p, kw_p = kw_p,
    recommendation = recommendation, tukey = tukey
  )
  class(result) <- "group_diff_result"
  result
}

print.group_diff_result=function(x, ...) {
  cat(sprintf(
    "\n--- %s ~ %s ---\nANOVA p = %.3f | Levene p = %.3f | Shapiro(resid) p = %s | Kruskal-Wallis p = %.3f\n%s\n",
    x$response, x$group, x$aov_p, x$levene_p,
    ifelse(is.na(x$shapiro_p), "NA", sprintf("%.3f", x$shapiro_p)),
    x$kw_p, x$recommendation
  ))
  if (!is.null(x$tukey)) {
    cat("Significant ANOVA -> Tukey HSD post-hoc:\n")
    print(x$tukey)
  }
  invisible(x)
}


test_correlation=function(data, var1, var2, method = "spearman") {
  d <- data[complete.cases(data[, c(var1, var2)]), ]
  norm1 <- tryCatch(shapiro.test(d[[var1]])$p.value, error = function(e) NA)
  norm2 <- tryCatch(shapiro.test(d[[var2]])$p.value, error = function(e) NA)
  ct <- cor.test(d[[var1]], d[[var2]], method = method)
  
  cat(sprintf(
    "%s vs %s (%s): estimate = %.3f, p = %.3f | Shapiro %s p=%s, %s p=%s\n",
    var1, var2, method, unname(ct$estimate), ct$p.value,
    var1, ifelse(is.na(norm1), "NA", sprintf("%.3f", norm1)),
    var2, ifelse(is.na(norm2), "NA", sprintf("%.3f", norm2))
  ))
  invisible(ct)
}

### DCI ###

dci_response_vars=c("Percent", "RDB", "Dead_Sh", "Galls", "Sd", "Bdor", "Gons")
dci_groups=c("Variety", "Country")

dci_results=list()
for (resp in dci_response_vars) {
  for (grp in dci_groups) {
    key <- paste(resp, grp, sep = "_")
    dci_results[[key]] <- analyze_group_differences(dci, resp, grp)
    print(dci_results[[key]])
    print(plot_distribution(dci, resp, grp,
                            title = paste(resp, "by", grp)))
  }
}


for (resp in dci_response_vars) {
  for (grp in dci_groups) {
    plotmeans(as.formula(paste(resp, "~", grp)), data = dci,
              xlab = grp, ylab = resp,
              main = paste("Means with 95% CI:", resp, "by", grp),
              las = 1, col = "darksalmon")
  }
}


leaf_vars=c("DCI", "Galls", "RDB", "Dead_Sh")
leaf_results=setNames(
  lapply(leaf_vars, function(v) analyze_group_differences(dci, v, "Variety")),
  leaf_vars
)
invisible(lapply(leaf_results, print))
invisible(lapply(leaf_vars, function(v) print(plot_distribution(dci, v, "Variety"))))

### Physiology ###

phys_vars_simple=c("CCI") 
for (grp in c("Variety", "Country")) {
  res <- analyze_group_differences(physiologia, "CCI", grp)
  print(res)
  print(plot_distribution(physiologia, "CCI", grp))
}

Fv_phys=physiologia %>%
  select(Country, Variety, Fv.Fm, PI) %>%
  filter(complete.cases(.)) %>%
  slice(-40)

for (grp in c("Variety", "Country")) {
  print(analyze_group_differences(Fv_phys, "Fv.Fm", grp))
  print(plot_distribution(Fv_phys, "Fv.Fm", grp))
}

Leaf=physiologia %>% select(Variety, Country, L_Area) %>% filter(complete.cases(.))
for (grp in c("Variety", "Country")) {
  res <- analyze_group_differences(Leaf, "L_Area", grp)
  print(res)
  print(plot_distribution(Leaf, "L_Area", grp))
  if (res$levene_p < 0.05) {
    cat("Levene significant -> also reporting Welch's ANOVA:\n")
    print(oneway.test(as.formula(paste("L_Area ~", grp)), data = Leaf))
  }
}

Weight=physiologia %>%
  select(Variety, Country, T_W, D_W, F_W) %>%
  filter(complete.cases(.))

for (resp in c("T_W", "D_W", "F_W")) {
  for (grp in c("Variety", "Country")) {
    print(analyze_group_differences(Weight, resp, grp))
    print(plot_distribution(Weight, resp, grp))
  }
}

L=physiologia %>% select(Variety, Country, LMA, LDMC) %>% filter(complete.cases(.))
for (resp in c("LMA", "LDMC")) {
  for (grp in c("Variety", "Country")) {
    print(analyze_group_differences(L, resp, grp))
    print(plot_distribution(L, resp, grp))
  }
}

A=physiologia %>% select(Variety, Country, A_sat) %>% filter(complete.cases(.))
for (grp in c("Variety", "Country")) {
  res <- analyze_group_differences(A, "A_sat", grp)
  print(res)
  print(plot_distribution(A, "A_sat", grp))
  if (res$levene_p < 0.05) {
    cat("Levene significant -> also reporting Welch's ANOVA:\n")
    print(oneway.test(as.formula(paste("A_sat ~", grp)), data = A))
  }
}

### PCA ###

samples=physiologia %>% distinct(ID, Country, Variety)

dat=physiologia %>%
  group_by(ID, Variety) %>%
  summarise(
    PI    = mean(PI, na.rm = TRUE),
    Fv.Fm = mean(Fv.Fm, na.rm = TRUE),
    CCI   = mean(CCI, na.rm = TRUE),
    Gs    = mean(Gs, na.rm = TRUE),
    Trm   = mean(Trm, na.rm = TRUE),
    A_sat = mean(A_sat, na.rm = TRUE),
    LMA   = 1000 * mean(LMA, na.rm = TRUE),
    LDMC  = mean(LDMC, na.rm = TRUE),
    DCI   = mean(DCI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  na.exclude()

pca1=PCA(dat[, 3:11], graph = FALSE)
dimdesc(pca1)
pca1$eig

fviz_pca_ind(pca1, label = "none", habillage = dat$Variety)
fviz_pca_biplot(pca1, habillage = dat$Variety, col.var = "black",
                alpha.var = "cos2", label = "var")

dat2=dat %>% select(ID, Variety, PI, CCI, A_sat, LMA, LDMC, DCI)
pca2=PCA(dat2[, 3:8], graph = FALSE)
dimdesc(pca2)
pca2$eig
pca2$var$cor

fviz_pca_ind(pca2, label = "none", habillage = dat$Variety)
fviz_pca_biplot(pca2, habillage = dat$Variety, col.var = "black",
                alpha.var = "cos2", label = "var")


ggplot(physiologia, aes(CCI, DCI)) + geom_point() + theme_bw()
ggplot(physiologia, aes(LMA, A_sat)) + geom_point() + theme_bw()
ggplot(physiologia, aes(LDMC, Fv.Fm)) + geom_point() + theme_bw()
ggplot(physiologia, aes(LDMC, CCI)) + geom_point() + theme_bw()


PI_me=lmer(PI    ~ Variety + DCI + LMA + (1 | ID), data = physiologia)
CCI_me=lmer(CCI   ~ Variety + DCI + LMA + (1 | ID), data = physiologia)
FvFm_me=lmer(Fv.Fm ~ Variety + DCI + LMA + (1 | ID), data = physiologia)
Asat_me=lmer(A_sat ~ Variety + DCI + LMA + (1 | ID), data = physiologia)
Gs_me=lmer(Gs    ~ Variety + DCI + LMA + (1 | ID), data = physiologia)

invisible(lapply(list(PI_me, CCI_me, FvFm_me, Asat_me, Gs_me), function(m) print(summary(m))))

summary(lm(A_sat ~ Variety * DCI, dat))
summary(lm(PI    ~ Variety * DCI, dat))
summary(lm(Fv.Fm ~ Variety * DCI, dat))
summary(lm(LMA   ~ Variety * DCI, dat))

correlation_pairs=list(
  list(data = physiologia, v1 = "DCI", v2 = "Fv.Fm", method = "spearman"),
  list(data = physiologia, v1 = "DCI", v2 = "CCI",   method = "spearman"),
  list(data = physiologia, v1 = "DCI", v2 = "A_sat", method = "spearman"),
  list(data = physiologia, v1 = "DCI", v2 = "F_W",   method = "kendall"),
  list(data = physiologia, v1 = "DCI", v2 = "T_W",   method = "kendall"),
  list(data = physiologia, v1 = "DCI", v2 = "D_W",   method = "kendall"),
  list(data = physiologia, v1 = "DCI", v2 = "L_Area", method = "spearman"),
  list(data = physiologia, v1 = "DCI", v2 = "LDMC",  method = "spearman"),
  list(data = physiologia, v1 = "DCI", v2 = "LMA",   method = "kendall"),
  list(data = physiologia, v1 = "RDB", v2 = "Fv.Fm", method = "spearman"),
  list(data = physiologia, v1 = "RDB", v2 = "CCI",   method = "kendall"),
  list(data = physiologia, v1 = "RDB", v2 = "A_sat", method = "kendall"),
  list(data = physiologia, v1 = "RDB", v2 = "L_Area", method = "spearman"),
  list(data = physiologia, v1 = "RDB", v2 = "LDMC",  method = "spearman"),
  list(data = physiologia, v1 = "RDB", v2 = "LMA",   method = "spearman"),
  list(data = dci,         v1 = "RDB", v2 = "Dead_Sh", method = "kendall"),
  list(data = dci,         v1 = "DCI", v2 = "RDB",   method = "kendall"),
  list(data = dci,         v1 = "Galls", v2 = "RDB", method = "spearman"),
  list(data = dci,         v1 = "DCI", v2 = "Sprouts", method = "spearman"),
  list(data = dci,         v1 = "Sprouts", v2 = "Dead_Sh", method = "kendall"),
  list(data = dci,         v1 = "Dead_Sh", v2 = "Galls", method = "kendall")
)

invisible(lapply(correlation_pairs, function(p) {
  test_correlation(p$data, p$v1, p$v2, p$method)
}))


DCI_full=lm(DCI ~ Length + Sprouts + Alive_Sh + Dead_Sh + RDB + Galls +
                 Variety + Country, data = dci)
summary(DCI_full)
DCI_stepwise <- step(DCI_full)
summary(DCI_stepwise)

DCI_reduced=lm(DCI ~ Sd + Bdor + Gons, data = dci)
summary(DCI_reduced)

### Map ###

map=read.table("map.txt", header = TRUE)
ggplot(map, aes(x = x, y = y)) +
  geom_point(shape = 17, color = "darkgreen", fill = "darkgreen") +
  geom_text(label = map$Variety, nudge_x = 0.25, nudge_y = 0.50, check_overlap = TRUE) +
  theme_classic()

### Correlations viz ###

cor_vars=physiologia %>%
  select(Fv.Fm, CCI, DCI, RDB, F_W, D_W, T_W, L_Area, LDMC, LMA) %>%
  filter(complete.cases(.))

ggpairs(cor_vars, title = "Correlogram")
ggcorr(cor_vars, method = c("everything", "pearson"))
cor(cor_vars)

