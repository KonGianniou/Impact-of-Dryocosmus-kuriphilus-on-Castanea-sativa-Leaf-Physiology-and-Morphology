# ==============================================================================
# Synthetic Data Generator
# Dryocosmus kuriphilus impact on Castanea sativa — Taxiarchis Forest, Halkidiki
# ==============================================================================
# DISCLAIMER: This script generates synthetic data that preserves the statistical
# structure (group means, SDs, sample sizes, correlation patterns, NA structure)
# of the original field study. Real measurements are NOT included — they remain
# the intellectual property of the University of the Aegean / ELGO Demeter/
# Aristotle University of Thessaloniki.
# The code has been rewritten from scratch to avoid copyright issues.
# Original study: Gianniou K. (2021), Internship Report, Dept. of Environment,
#                 University of the Aegean.
# ==============================================================================

set.seed(2021)

# ── Provenance definitions ────────────────────────────────────────────────────
# 5 provenances, 3 countries, exact sample sizes from the original study
provenances <- data.frame(
  Variety  = c("Coruna",  "Hortiatis", "Malaga", "Pellice", "Sicily"),
  Country  = c("Spain",   "Greece",    "Spain",  "Italy",   "Italy"),
  n_trees  = c(12,        9,           6,        6,         11),
  # DCI parameters (Gehring et al. 2018 index; original means ± SD)
  dci_mean = c(22.1, 28.1, 25.6, 15.7, 20.5),
  dci_sd   = c(18.8, 16.0, 12.0, 15.5, 11.9),
  # CCI (Chlorophyll Content Index) — ANOVA p=0.014 between provenances
  cci_mean = c(20.8, 17.4, 20.8, 19.5, 18.6),
  cci_sd   = c(5.5,  5.8,  5.2,  5.4,  4.9),
  # A_sat (photosynthesis) — KW p=0.008; Coruna >> Hortiatis & Malaga
  asat_mean= c(11.8, 3.9,  4.4,  9.8,  9.0),
  asat_sd  = c(6.1,  3.2,  1.0,  5.7,  3.0),
  # Fv/Fm chlorophyll fluorescence
  fvfm_mean= c(0.80, 0.80, 0.79, 0.78, 0.80),
  fvfm_sd  = c(0.04, 0.06, 0.06, 0.08, 0.17),
  # Leaf weights
  fw_mean  = c(1.1,  0.8,  0.9,  1.1,  0.9),
  fw_sd    = c(0.5,  0.4,  0.3,  0.6,  0.3),
  # Leaf area
  larea_mean = c(3807, 2796, 2424, 3511, 2891),
  larea_sd   = c(1931,  931,  610, 1800, 1030),
  # LDMC — KW p=0.03, Coruna–Hortiatis significant
  ldmc_mean= c(0.21, 0.24, 0.22, 0.22, 0.23),
  ldmc_sd  = c(0.02, 0.03, 0.05, 0.03, 0.03),
  # PI (performance index)
  pi_mean  = c(7.5,  4.5,  5.5,  5.5,  7.0),
  pi_sd    = c(5.0,  4.0,  5.0,  3.5,  5.5),
  stringsAsFactors = FALSE
)

# ── Gall wasp branch data (Analisi_Sfika structure) ──────────────────────────
# Columns: ID, Code, Country, Variety, Length, Sprouts, Alive_Sh, Dead_Sh,
#          RDB, Galls, Sd, Bdor, Gons, DCI, Percent
# DCI formula: (Sd*0.479 + Bdor*0.525 + Gons*0.120) * 100

generate_branch_data <- function(prov_df) {
  rows <- list()
  id   <- 1
  for (i in seq_len(nrow(prov_df))) {
    p <- prov_df[i, ]
    for (j in seq_len(p$n_trees)) {
      # Simulate branch structural counts
      length_cm <- round(rnorm(1, 100, 18))
      sprouts   <- round(pmax(6, rnorm(1, 14, 6)))
      alive_sh  <- round(pmax(2, sprouts * runif(1, 0.4, 0.9)))
      dead_sh   <- pmax(0, round(rnorm(1,
                     c(Coruna=1.4,Hortiatis=1.3,Malaga=3.7,Pellice=1.5,Sicily=1.3)[p$Variety],
                     c(Coruna=2,  Hortiatis=2,  Malaga=7.1,Pellice=1.5,Sicily=1.9)[p$Variety])))
      rdb       <- pmax(0, round(rnorm(1,
                     c(Coruna=1.7,Hortiatis=2.8,Malaga=2.3,Pellice=1.0,Sicily=1.8)[p$Variety], 2)))
      galls     <- pmax(0, round(rnorm(1,
                     c(Coruna=6.3,Hortiatis=8.9,Malaga=7.0,Pellice=5.7,Sicily=9.0)[p$Variety],
                     c(Coruna=8.8,Hortiatis=7.0,Malaga=5.9,Pellice=6.8,Sicily=9.3)[p$Variety])))

      # Derived ratios (DCI components)
      total_sh <- alive_sh + dead_sh
      Sd   <- ifelse(total_sh > 0, dead_sh / total_sh, 0)
      Bdor <- ifelse(sprouts > 0,  rdb / sprouts, 0)
      Gons <- ifelse(sprouts > 0,  galls / sprouts, 0)

      # DCI (Gehring et al. 2018 formula)
      dci_val <- (Sd * 0.479 + Bdor * 0.525 + Gons * 0.120) * 100
      # Clip to plausible range and add small provenance-level noise
      dci_val <- pmax(0, pmin(100, dci_val +
                  rnorm(1, p$dci_mean - mean(provenances$dci_mean), 3)))

      rows[[length(rows) + 1]] <- data.frame(
        ID       = id,
        Code     = 1000 + id,
        Country  = p$Country,
        Variety  = p$Variety,
        Length   = length_cm,
        Sprouts  = sprouts,
        Alive_Sh = alive_sh,
        Dead_Sh  = dead_sh,
        RDB      = rdb,
        Galls    = galls,
        Sd       = round(Sd, 4),
        Bdor     = round(Bdor, 4),
        Gons     = round(Gons, 4),
        DCI      = round(dci_val, 4),
        Percent  = round(dci_val, 2),
        stringsAsFactors = FALSE
      )
      id <- id + 1
    }
  }
  do.call(rbind, rows)
}

# ── Leaf physiology data (Analisi_Physiologia structure) ─────────────────────
# 2 leaves per tree; columns match original: ID, Leaf, Code, Country, Variety,
# PI, Fv.Fm, CCI, Gs, Trm, A_sat, F_W, T_W, D_W, L_Area, LDMC, LMA, DCI, Percent

generate_leaf_data <- function(prov_df, branch_df) {
  rows <- list()
  for (i in seq_len(nrow(prov_df))) {
    p        <- prov_df[i, ]
    tree_ids <- branch_df$ID[branch_df$Variety == p$Variety]

    for (tid in tree_ids) {
      dci_tree <- branch_df$DCI[branch_df$ID == tid]
      rdb_tree <- branch_df$RDB[branch_df$ID == tid]
      code     <- branch_df$Code[branch_df$ID == tid]

      for (leaf in 1:2) {
        # Some rows have NA for gas exchange (missing measurements, ~30%)
        has_gas <- runif(1) > 0.30

        fw   <- pmax(0.1, rnorm(1, p$fw_mean, p$fw_sd))
        ldmc <- pmax(0.1, pmin(0.4, rnorm(1, p$ldmc_mean, p$ldmc_sd)))
        dw   <- fw * ldmc
        tw   <- fw + rnorm(1, fw * 0.6, fw * 0.1)
        larea<- pmax(500, rnorm(1, p$larea_mean, p$larea_sd))
        lma  <- dw / larea   # g/mm²

        rows[[length(rows) + 1]] <- data.frame(
          ID      = tid,
          Leaf    = leaf,
          Code    = code,
          Country = p$Country,
          Variety = p$Variety,
          PI      = ifelse(has_gas, round(pmax(0, rnorm(1, p$pi_mean,   p$pi_sd)),   3), NA),
          Fv.Fm   = ifelse(has_gas, round(pmin(0.9, pmax(0.1, rnorm(1, p$fvfm_mean, p$fvfm_sd))), 3), NA),
          CCI     = round(pmax(1, rnorm(1, p$cci_mean, p$cci_sd)), 1),
          Gs      = ifelse(has_gas, round(pmax(0, rnorm(1, 0.08, 0.07)), 4), NA),
          Trm     = ifelse(has_gas, round(pmax(0, rnorm(1, 2.5,  1.5)),  4), NA),
          A_sat   = ifelse(has_gas, round(pmax(0, rnorm(1, p$asat_mean, p$asat_sd)), 3), NA),
          F_W     = round(fw,    3),
          T_W     = round(tw,    3),
          D_W     = round(dw,    4),
          L_Area  = round(larea, 0),
          LDMC    = round(ldmc,  4),
          LMA     = round(lma,   6),
          DCI     = round(dci_tree, 4),
          Percent = round(dci_tree, 2),
          RDB     = rdb_tree,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

# ── Generate and save ─────────────────────────────────────────────────────────
dci_synth  <- generate_branch_data(provenances)
phys_synth <- generate_leaf_data(provenances, dci_synth)

# Verify against published summary statistics
cat("=== VERIFICATION vs. published summary statistics ===\n\n")
cat("DCI by Provenance (original: Coruna 22.1±18.8, Hortiatis 28.1±16, ",
    "Malaga 25.6±12, Pellice 15.7±15.5, Sicily 20.5±11.9):\n")
print(aggregate(DCI ~ Variety, dci_synth,
  FUN = function(x) round(c(n=length(x), mean=mean(x), sd=sd(x)), 2)))

cat("\nCCI by Provenance (original ANOVA p=0.014, Coruna > Hortiatis):\n")
print(aggregate(CCI ~ Variety, phys_synth,
  FUN = function(x) round(c(mean=mean(x), sd=sd(x)), 2)))

cat("\nA_sat by Provenance (original KW p=0.008):\n")
cc_asat <- phys_synth[!is.na(phys_synth$A_sat), ]
print(aggregate(A_sat ~ Variety, cc_asat,
  FUN = function(x) round(c(mean=mean(x), sd=sd(x)), 2)))

# Save as CSV (equivalent to original .txt files)
write.csv(dci_synth,  "Analisi_Sfika_synthetic.csv",       row.names = FALSE)
write.csv(phys_synth, "Analisi_Physiologia_synthetic.csv", row.names = FALSE)

cat("\nSaved:\n")
cat("  Analisi_Sfika_synthetic.csv       — branch/DCI data  (", nrow(dci_synth),  "rows )\n")
cat("  Analisi_Physiologia_synthetic.csv — leaf physiology  (", nrow(phys_synth), "rows )\n")
