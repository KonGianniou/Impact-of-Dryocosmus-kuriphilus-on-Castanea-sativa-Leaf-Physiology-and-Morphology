# Impact of Dryocosmus kuriphilus on Castanea sativa Leaf Physiology and Morphology
Author: Konstantina Gianniou  
Institution: Department of Environment, University of the Aegean  
Study site: University Forest of Taxiarchis, Halkidiki, Greece (40°44'35"N, 23°18'12"E, 760 m a.s.l.)  
Year: 2021
---
⚠️ Data & Code Disclaimer
> \\\*\\\*The data in this repository are entirely synthetic.\\\*\\\*
>
> The original field measurements were collected as part of an internship within the University of the Aegean and remain the intellectual property of the University of the Aegean, ELGO Demeter (Forest Research Institute) and Aristotle University of Thessaloniki. They are not included here.
>
> The synthetic dataset (`Analisi\\\_Sfika\\\_synthetic.csv`, `Analisi\\\_Physiologia\\\_synthetic.csv`) was generated using `generate\\\_synthetic\\\_data.R`, which reproduces the \\\*\\\*statistical structure\\\*\\\* of the original data — including group sample sizes, provenance-level means and standard deviations, NA patterns, and inter-variable correlation structure — without exposing any real measurements.
>
> \\\*\\\*The analysis code has also been rewritten from scratch\\\*\\\* to avoid reproducing the original internship report verbatim. The statistical methods, test choices, and interpretation are faithful to the original study; the implementation is new.
---
Overview
The Asian Chestnut Gall Wasp (Dryocosmus kuriphilus Yasumatsu, 1951; Hymenoptera: Cynipidae) is an invasive species endemic to China that reached Europe via Italy in 2002, and was first recorded in Greece in 2014. It induces characteristic galls on chestnut (Castanea sativa Mill.) buds and shoots, causing leaf deformation, reduced photosynthetic capacity, and declines in flowering and fruiting.
This study quantifies the impact of D. kuriphilus on the physiology and morphology of C. sativa leaves across five provenances from three countries (Greece, Italy, Spain) grown together in a common-garden plantation, eliminating confounding environmental variation.
---
Research Questions
Does infestation severity (DCI index) differ significantly between provenances or countries of origin?
Do reactivated dormant buds (RDB), gall counts, and dead shoots differ between provenances?
Are physiological leaf traits (CCI, A_sat, Fv/Fm, PI) affected by provenance or infestation level?
Which morphological traits (LMA, LDMC, leaf area, weight) are associated with infestation severity?
What is the multivariate structure of relationships between leaf traits and DCI?
---
Study Design
Sampling
44 trees randomly selected from a common-garden plantation of 143 individuals
Provenances: Coruna (Spain, n=12), Hortiatis (Greece, n=9), Malaga (Spain, n=6), Pellice (Italy, n=6), Sicily (Italy, n=11)
For each tree: one branch (≥50 cm) + two leaves collected (July–August)
Infestation index
The DCI (Damage by Cynipid Infestation) index was computed following Gehring et al. (2018):
```
DCI = (Sd × 0.479 + Bdor × 0.525 + Gons × 0.120) × 100
```
where:
Sd = dead shoots / total shoots
Bdor = reactivated dormant buds / total green shoots
Gons = galls on sprouts / total green shoots
Leaf measurements
Variable	Description
CCI	Chlorophyll Content Index (field measurement)
A_sat	Saturated photosynthesis rate (μmol m⁻² s⁻¹)
Fv/Fm	Maximum quantum yield of PSII (chlorophyll fluorescence)
PI	Performance Index
Gs	Stomatal conductance
F_W / T_W / D_W	Fresh, saturated, and dry leaf weight (g)
L_Area	Leaf area (mm²) — scanned and digitised
LMA	Leaf Mass per Area = D_W / L_Area (g mm⁻²)
LDMC	Leaf Dry Matter Content = D_W / T_W
---
Statistical Methods
Infestation indices (DCI, RDB, Galls, Dead Shoots)
One-way ANOVA to test for differences between provenances and countries
Levene's test for homogeneity of variances (assumption check)
Shapiro-Wilk test on ANOVA residuals (normality check)
Kruskal-Wallis test as non-parametric backup when residuals are non-normal
Leaf physiological traits (CCI, A_sat, Fv/Fm, etc.)
Same ANOVA + Levene + Shapiro-Wilk + Kruskal-Wallis pipeline
Welch's ANOVA (`oneway.test`) when Levene's test indicates unequal variances
Tukey HSD post-hoc for significant ANOVAs
Correlations
Shapiro-Wilk to select correlation method
Spearman's rho (primary): non-parametric, robust to non-normality and outliers
Kendall's tau: reported alongside Spearman for key pairs with small n
Multivariate analysis
Principal Component Analysis (PCA) via `FactoMineR::PCA()` on tree-level averages
Two PCAs: (1) full trait set including gas exchange; (2) reduced set (PI, CCI, A_sat, LMA, LDMC, DCI)
Biplots produced with `factoextra::fviz\\\_pca\\\_biplot()`, coloured by provenance
Mixed-effects models
`lmer(trait \\\~ Variety + DCI + LMA + (1|ID))` — tree ID as random intercept to account for the two-leaves-per-tree structure (pseudoreplication)
Fit via `lme4`; p-values from `lmerTest` (Satterthwaite approximation)
---
Key Findings (from the original study)
Finding	Result
DCI by provenance	No significant differences (ANOVA p = 0.6; KW p > 0.05)
CCI by provenance	Significant (ANOVA p = 0.014); Hortiatis < Coruna (Tukey HSD)
CCI by country	Significant (ANOVA p = 0.01); Spain > Greece
A_sat by provenance	Significant (KW p = 0.008); Coruna > Hortiatis & Malaga
LDMC by provenance	Marginal (KW p = 0.03); Coruna < Hortiatis
DCI ~ F_W	Significant positive correlation (Spearman rho = 0.23, p = 0.04)
DCI ~ D_W	Significant positive correlation (rho = 0.23, p = 0.03)
DCI ~ LMA	Significant positive correlation (rho = 0.30, p = 0.005)
RDB ~ DCI	Strong positive correlation (rho = 0.75, p < 0.001)
---
Repository Structure
```
dryocosmus-castanea/
├── generate\\\_synthetic\\\_data.R     ← Run this first to create the datasets
├── analysis.R                    ← Full statistical pipeline
├── fig1\\\_DCI\\\_provenance.png       ← DCI violin-boxplot by provenance
├── fig2\\\_RDB\\\_provenance.png       ← RDB boxplot by provenance
├── fig3\\\_CCI.png                  ← CCI by provenance and country
├── fig4\\\_PCA\\\_full.png             ← PCA biplot (all traits)
├── fig5\\\_PCA\\\_reduced.png          ← PCA biplot (reduced trait set)
├── fig6\\\_correlation\\\_matrix.png   ← Spearman correlation matrix
└── README.md
```
---
How to Run
```r
# Step 1: generate the synthetic datasets
source("generate\\\_synthetic\\\_data.R")

# Step 2: run the full analysis (installs missing packages automatically)
source("analysis.R")
```
Requires R ≥ 4.1. All packages install automatically on first run.
---
References
Gehring E., Bellosi B., Quacchia A. & Conedera M. (2018). Evaluating Dryocosmus kuriphilus-induced damage on Castanea sativa. Journal of Pest Science, Springer.
R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.
---
Contact
Konstantina Gianniou  
g.tem2106@gmail.com | LinkedIn |
