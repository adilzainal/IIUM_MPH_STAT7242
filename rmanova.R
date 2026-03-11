library(haven)
data <- read_sav("Desktop/Medical Statistic/R MPH IIUM/stat7242/wideformat.sav")
View(data)

# =============================================================================
#  REPEATED MEASURES ANOVA — COMPLETE ANALYSIS SCRIPT
#  Dataset  : data  (wide format)
#  Within   : time (day1, day2, day3, day4)
#  Between  : diet (0 = Control, 1 = Intervention)
#  ID var   : id
#
#  Sections:
#    0. Install & load packages
#    1. Prepare data
#    2. Descriptive statistics & profile plot
#    3. Assumption checks
#       3a. Normality (Shapiro-Wilk per group × time)
#       3b. Outliers (boxplots, standardised residuals)
#       3c. Homogeneity of variance (Levene's test)
#       3d. Sphericity (Mauchly's — auto in ezANOVA)
#    4. Repeated Measures ANOVA
#       4a. Main effects (time, diet)
#       4b. Interaction (time × diet)
#       4c. Sphericity corrections (G-G / H-F)
#    5. Post-hoc comparisons (emmeans)
#    6. Linear Mixed Model (handles missing data)
#    7. Effect sizes
#    8. Publication-quality figures
# =============================================================================


# ─── 0. INSTALL & LOAD PACKAGES ───────────────────────────────────────────────

# Run this block once to install required packages
required_pkgs <- c("tidyverse", "ez", "emmeans", "rstatix",
                   "ggpubr", "nlme", "performance", "effectsize",
                   "patchwork", "scales")

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(tidyverse)      # data wrangling + ggplot2
library(ez)             # ezANOVA — RM-ANOVA + Mauchly's test
library(emmeans)        # post-hoc pairwise comparisons
library(rstatix)        # Shapiro-Wilk, Levene, outlier helpers
library(ggpubr)         # publication-ready ggplot themes
library(nlme)           # linear mixed model
library(performance)    # ICC
library(effectsize)     # eta_squared
library(patchwork)      # combine ggplot panels
library(scales)         # axis formatting


# ─── 1. PREPARE DATA ──────────────────────────────────────────────────────────

# ── 1a. Factor coding ─────────────────────────────────────────────────────────
data <- data %>%
  mutate(
    id   = factor(id),
    diet = factor(diet, levels = c(0, 1), labels = c("Control", "Intervention"))
  )

# ── 1b. Reshape to LONG format (required for ez, emmeans, rstatix) ────────────
data_long <- data %>%
  pivot_longer(
    cols      = c(day1, day2, day3, day4),
    names_to  = "time",
    values_to = "outcome"
  ) %>%
  mutate(
    time = factor(time,
                  levels = c("day1", "day2", "day3", "day4"),
                  labels = c("Day 1", "Day 2", "Day 3", "Day 4"))
  )

# Quick check
cat("\n── Data structure ──────────────────────────────\n")
glimpse(data_long)
cat("\nSample sizes per group:\n")
data_long %>% count(diet, time) %>% print()
cat("\nMissing values per variable:\n")
colSums(is.na(data_long)) %>% print()


# ─── 2. DESCRIPTIVE STATISTICS & PROFILE PLOT ─────────────────────────────────

desc_stats <- data_long %>%
  group_by(diet, time) %>%
  summarise(
    n    = n(),
    mean = mean(outcome, na.rm = TRUE),
    sd   = sd(outcome,   na.rm = TRUE),
    se   = sd / sqrt(n),
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )

cat("\n── Descriptive statistics (Mean ± SD) ──────────\n")
desc_stats %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  print()

# Profile plot (means ± SE)
p_profile <- ggplot(desc_stats, aes(x = time, y = mean,
                                    colour = diet, group = diet)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.15, linewidth = 0.9) +
  scale_colour_manual(values = c("Control" = "#0E9E8E", "Intervention" = "#F59E0B"),
                      name = "Diet group") +
  labs(title    = "Profile Plot: Outcome over Time by Diet Group",
       subtitle = "Mean ± SE",
       x = "Time point", y = "Outcome") +
  theme_pubr(base_size = 13) +
  theme(legend.position = "top",
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey50"))

print(p_profile)


# ─── 3. ASSUMPTION CHECKS ─────────────────────────────────────────────────────

cat("\n\n═══════════════════════════════════════════════\n")
cat("  SECTION 3: ASSUMPTION CHECKS\n")
cat("═══════════════════════════════════════════════\n")

# ── 3a. Normality — Shapiro-Wilk per group × time ─────────────────────────────
cat("\n── 3a. Shapiro-Wilk Normality Test (per diet × time) ──\n")
cat("  Interpretation: p > 0.05 = normality NOT violated\n\n")

normality_results <- data_long %>%
  group_by(diet, time) %>%
  shapiro_test(outcome) %>%
  mutate(
    interpretation = ifelse(p > 0.05, "Normal ✔", "Non-normal ✗")
  )
print(normality_results)

# Q-Q plots per group × time
p_qq <- ggplot(data_long, aes(sample = outcome, colour = diet)) +
  stat_qq(size = 1.5) +
  stat_qq_line(linewidth = 0.8) +
  facet_grid(diet ~ time) +
  scale_colour_manual(values = c("Control" = "#0E9E8E", "Intervention" = "#F59E0B"),
                      guide = "none") +
  labs(title    = "Q-Q Plots: Normality Check",
       subtitle = "Points should follow the diagonal line",
       x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_pubr(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

print(p_qq)

# Histogram of residuals (from long format means)
data_long <- data_long %>%
  group_by(diet, time) %>%
  mutate(group_mean = mean(outcome, na.rm = TRUE),
         residual   = outcome - group_mean) %>%
  ungroup()

p_hist <- ggplot(data_long, aes(x = residual, fill = diet)) +
  geom_histogram(bins = 20, colour = "white", alpha = 0.8) +
  facet_wrap(~ diet) +
  scale_fill_manual(values = c("Control" = "#0E9E8E", "Intervention" = "#F59E0B"),
                    guide = "none") +
  labs(title    = "Histogram of Residuals",
       subtitle = "Should approximate a bell curve",
       x = "Residual", y = "Count") +
  theme_pubr(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

print(p_hist)


# ── 3b. Outliers ──────────────────────────────────────────────────────────────
cat("\n── 3b. Outlier Detection ──────────────────────────────\n")
cat("  Method: rstatix identify_outliers (|z| > 2 = outlier, |z| > 3.29 = extreme)\n\n")

outlier_check <- data_long %>%
  group_by(diet, time) %>%
  identify_outliers(outcome)

if (nrow(outlier_check) == 0) {
  cat("  No outliers detected across any group × time combination.\n")
} else {
  cat("  Outliers found:\n")
  print(outlier_check)
  extreme <- outlier_check %>% filter(is.extreme == TRUE)
  if (nrow(extreme) > 0) {
    cat("\n  ⚠ EXTREME outliers (investigate these):\n")
    print(extreme)
  }
}

# Boxplots per time × diet
p_box <- ggplot(data_long, aes(x = time, y = outcome, fill = diet)) +
  geom_boxplot(alpha = 0.75, outlier.colour = "red", outlier.size = 2.5,
               position = position_dodge(0.8), width = 0.6) +
  scale_fill_manual(values = c("Control" = "#0E9E8E", "Intervention" = "#F59E0B"),
                    name = "Diet group") +
  labs(title    = "Boxplots: Outcome Distribution by Time and Diet",
       subtitle = "Red points = outliers",
       x = "Time point", y = "Outcome") +
  theme_pubr(base_size = 13) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

print(p_box)


# ── 3c. Homogeneity of variance (Levene's test) ───────────────────────────────
cat("\n── 3c. Levene's Test for Homogeneity of Variance (per time) ──\n")
cat("  Interpretation: p > 0.05 = variances are homogeneous\n\n")

levene_results <- data_long %>%
  group_by(time) %>%
  levene_test(outcome ~ diet) %>%
  mutate(interpretation = ifelse(p > 0.05, "Homogeneous ✔", "Heterogeneous ✗"))

print(levene_results)


# ── 3d. Sphericity note ───────────────────────────────────────────────────────
cat("\n── 3d. Sphericity ─────────────────────────────────────────────\n")
cat("  Mauchly's test will be reported automatically in Section 4 (ezANOVA).\n")
cat("  If violated (p < .05): Greenhouse-Geisser (ε < 0.75) or\n")
cat("                         Huynh-Feldt (ε ≥ 0.75) correction is applied.\n")


# ─── 4. REPEATED MEASURES ANOVA ───────────────────────────────────────────────

cat("\n\n═══════════════════════════════════════════════\n")
cat("  SECTION 4: REPEATED MEASURES ANOVA\n")
cat("═══════════════════════════════════════════════\n")

rm_model <- ezANOVA(
  data       = data_long,
  dv         = outcome,
  wid        = id,
  within     = time,
  between    = diet,
  type       = 3,          # Type III SS (standard for unbalanced designs)
  return_aov = TRUE
)

# ── 4a. ANOVA table ───────────────────────────────────────────────────────────
cat("\n── 4a. ANOVA Table ─────────────────────────────────────────\n")
cat("  F, df, p-value, and generalised eta-squared (ges)\n\n")
print(rm_model$ANOVA)

# ── 4b. Mauchly's Sphericity Test ─────────────────────────────────────────────
cat("\n── 4b. Mauchly's Test for Sphericity ──────────────────────\n")
cat("  Applies to the within-subject factor (time) and time:diet interaction\n")
cat("  p > 0.05 = sphericity NOT violated\n\n")
print(rm_model$`Mauchly's Test for Sphericity`)

# ── 4c. Sphericity corrections ────────────────────────────────────────────────
cat("\n── 4c. Sphericity Corrections ──────────────────────────────\n")
cat("  GGe = Greenhouse-Geisser epsilon\n")
cat("  HFe = Huynh-Feldt epsilon\n")
cat("  p[GG] and p[HF] = corrected p-values\n\n")
print(rm_model$`Sphericity Corrections`)

# Summary guidance
mauchly_p <- rm_model$`Mauchly's Test for Sphericity`$p
cat("\n── Auto-guidance ────────────────────────────────────────────\n")
for (i in seq_along(mauchly_p)) {
  effect <- rm_model$`Mauchly's Test for Sphericity`$Effect[i]
  p_val  <- mauchly_p[i]
  GGe    <- rm_model$`Sphericity Corrections`$GGe[i]
  if (p_val > 0.05) {
    cat(sprintf("  ✔  %s: Sphericity met (p = %.3f) → report uncorrected F\n",
                effect, p_val))
  } else if (GGe < 0.75) {
    cat(sprintf("  ✗  %s: Sphericity VIOLATED (p = %.3f), ε = %.3f < 0.75 → use G-G correction\n",
                effect, p_val, GGe))
  } else {
    cat(sprintf("  ✗  %s: Sphericity VIOLATED (p = %.3f), ε = %.3f ≥ 0.75 → use H-F correction\n",
                effect, p_val, GGe))
  }
}


# ─── 5. POST-HOC COMPARISONS ──────────────────────────────────────────────────

cat("\n\n═══════════════════════════════════════════════\n")
cat("  SECTION 5: POST-HOC COMPARISONS (emmeans)\n")
cat("═══════════════════════════════════════════════\n")

# Fit LMM for emmeans (ezANOVA aov object works too)
lme_fit <- lme(outcome ~ time * diet,
               random  = ~ 1 | id,
               data    = data_long,
               method  = "REML",
               na.action = na.omit)

# ── 5a. Simple effects: pairwise time comparisons within each diet group ───────
cat("\n── 5a. Pairwise Time Comparisons WITHIN Each Diet Group ────\n")
cat("  (Bonferroni adjusted; key for interpreting interaction)\n\n")

emm_time_by_diet <- emmeans(lme_fit, ~ time | diet)
posthoc_time     <- pairs(emm_time_by_diet, adjust = "bonferroni")
cat("  Time contrasts within Control group:\n")
print(summary(posthoc_time)[summary(posthoc_time)$diet == "Control", ])
cat("\n  Time contrasts within Intervention group:\n")
print(summary(posthoc_time)[summary(posthoc_time)$diet == "Intervention", ])

# ── 5b. Between-group (diet) comparisons at each time point ────────────────────
cat("\n── 5b. Diet Group Comparisons AT Each Time Point ──────────\n")
cat("  (Bonferroni adjusted; tests where groups differ at each day)\n\n")

emm_diet_by_time <- emmeans(lme_fit, ~ diet | time)
posthoc_diet     <- pairs(emm_diet_by_time, adjust = "bonferroni")
print(summary(posthoc_diet))

# ── 5c. Overall interaction contrasts ─────────────────────────────────────────
cat("\n── 5c. Interaction Contrasts (time × diet) ─────────────────\n")
emm_interaction <- emmeans(lme_fit, ~ time * diet)
cat("\n  Estimated marginal means for all time × diet cells:\n")
print(emm_interaction)


# ─── 6. LINEAR MIXED MODEL ────────────────────────────────────────────────────

cat("\n\n═══════════════════════════════════════════════\n")
cat("  SECTION 6: LINEAR MIXED MODEL (nlme)\n")
cat("═══════════════════════════════════════════════\n")

cat("\n── Fixed Effects Table ─────────────────────────────────────\n")
print(summary(lme_fit)$tTable)

cat("\n── 95% Confidence Intervals for Fixed Effects ──────────────\n")
print(intervals(lme_fit, which = "fixed"))

cat("\n── Variance Components ──────────────────────────────────────\n")
print(VarCorr(lme_fit))

cat("\n── ICC (Intraclass Correlation) ─────────────────────────────\n")
cat("  ICC = between-subject variance / total variance\n")
icc_val <- performance::icc(lme_fit)
print(icc_val)

cat("\n── Marginal & Conditional R² ────────────────────────────────\n")
cat("  R²m = fixed effects only | R²c = fixed + random\n")
# MuMIn method (install if needed)
if (!requireNamespace("MuMIn", quietly = TRUE)) install.packages("MuMIn")
library(MuMIn)
r2_vals <- r.squaredGLMM(lme_fit)
cat(sprintf("  R²marginal  = %.4f\n", r2_vals[1]))
cat(sprintf("  R²conditional = %.4f\n", r2_vals[2]))

# Covariance structure comparison (model fit)
cat("\n── Comparing Covariance Structures (AIC) ────────────────────\n")
m_cs  <- lme(outcome ~ time * diet, random = ~1|id,
             correlation = corCompSymm(form = ~1|id),
             data = data_long, method = "ML", na.action = na.omit)
m_ar1 <- lme(outcome ~ time * diet, random = ~1|id,
             correlation = corAR1(form = ~1|id),
             data = data_long, method = "ML", na.action = na.omit)
m_base <- lme(outcome ~ time * diet, random = ~1|id,
              data = data_long, method = "ML", na.action = na.omit)

aic_compare <- AIC(m_base, m_cs, m_ar1)
rownames(aic_compare) <- c("Default (no corr)", "Compound Symmetry", "AR(1)")
cat("\n"); print(aic_compare)
cat("  → Use the model with the lowest AIC\n")


# ─── 7. EFFECT SIZES ──────────────────────────────────────────────────────────

cat("\n\n═══════════════════════════════════════════════\n")
cat("  SECTION 7: EFFECT SIZES\n")
cat("═══════════════════════════════════════════════\n")

cat("\n── Generalised Eta-Squared (from ezANOVA) ──────────────────\n")
cat("  Small = 0.01 | Medium = 0.06 | Large = 0.14\n\n")
ges_table <- rm_model$ANOVA %>%
  select(Effect, F, `p<.05`, ges) %>%
  mutate(
    ges_interpretation = case_when(
      ges < 0.01 ~ "Negligible",
      ges < 0.06 ~ "Small",
      ges < 0.14 ~ "Medium",
      TRUE       ~ "Large"
    )
  )
print(ges_table)

cat("\n── Cohen's d for Diet Group Comparisons at Each Time Point ─\n")
cat("  Cohen's d = (mean_intervention - mean_control) / pooled SD\n")
cat("  Small = 0.2 | Medium = 0.5 | Large = 0.8\n\n")

# Compute Cohen's d manually from descriptive stats (robust, no package dependency)
cohens_d_table <- desc_stats %>%
  select(diet, time, mean, sd, n) %>%
  pivot_wider(names_from = diet, values_from = c(mean, sd, n)) %>%
  mutate(
    # Pooled SD (Hedges / Cohen formula)
    pooled_sd = sqrt(
      ((n_Control - 1) * sd_Control^2 + (n_Intervention - 1) * sd_Intervention^2) /
        (n_Control + n_Intervention - 2)
    ),
    cohens_d = (mean_Intervention - mean_Control) / pooled_sd,
    d_interpretation = case_when(
      abs(cohens_d) < 0.2  ~ "Negligible",
      abs(cohens_d) < 0.5  ~ "Small",
      abs(cohens_d) < 0.8  ~ "Medium",
      TRUE                 ~ "Large"
    )
  ) %>%
  select(time, mean_Control, mean_Intervention, pooled_sd, cohens_d, d_interpretation) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

print(cohens_d_table)

cat("\n── Cohen's d for Time Comparisons Within Each Group ────────\n")
# Pairwise Cohen's d within each diet group across time points
cohens_d_time <- data_long %>%
  group_by(diet) %>%
  rstatix::cohens_d(outcome ~ time, paired = TRUE, hedges.correction = FALSE)
print(cohens_d_time)

# ─── 8. PUBLICATION-QUALITY FIGURES ──────────────────────────────────────────

cat("\n\n═══════════════════════════════════════════════\n")
cat("  SECTION 8: FIGURES\n")
cat("═══════════════════════════════════════════════\n")

# Colour palette
pal <- c("Control" = "#0E9E8E", "Intervention" = "#F59E0B")

# ── Fig 1: Profile plot with individual trajectories ─────────────────────────
fig1 <- ggplot() +
  # Individual lines (light, per subject)
  geom_line(data = data_long,
            aes(x = time, y = outcome, group = id, colour = diet),
            alpha = 0.2, linewidth = 0.5) +
  # Group means + SE
  geom_line(data = desc_stats,
            aes(x = time, y = mean, colour = diet, group = diet),
            linewidth = 1.8) +
  geom_point(data = desc_stats,
             aes(x = time, y = mean, colour = diet, shape = diet),
             size = 4) +
  geom_errorbar(data = desc_stats,
                aes(x = time, ymin = mean - se, ymax = mean + se, colour = diet),
                width = 0.12, linewidth = 1.1) +
  scale_colour_manual(values = pal, name = "Diet group") +
  scale_shape_manual(values = c("Control" = 16, "Intervention" = 17), name = "Diet group") +
  labs(
    title    = "Figure 1: Outcome Trajectories Over Time by Diet Group",
    subtitle = "Thin lines = individual subjects  |  Bold line = group mean ± SE",
    x = "Time point",
    y = "Outcome (Mean ± SE)"
  ) +
  theme_pubr(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey50", size = 10)
  )

print(fig1)

# ── Fig 2: Interaction plot with 95% CI ───────────────────────────────────────
fig2 <- ggplot(desc_stats, aes(x = time, y = mean,
                               colour = diet, group = diet,
                               fill   = diet)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1.6) +
  geom_point(aes(shape = diet), size = 4.5) +
  scale_colour_manual(values = pal, name = "Diet group") +
  scale_fill_manual(values = pal, name = "Diet group") +
  scale_shape_manual(values = c("Control" = 16, "Intervention" = 17), name = "Diet group") +
  labs(
    title    = "Figure 2: Time × Diet Interaction Plot",
    subtitle = "Shaded bands = 95% CI  |  Diverging lines = significant interaction",
    x = "Time point",
    y = "Outcome (Mean)"
  ) +
  theme_pubr(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey50", size = 10)
  )

print(fig2)

# ── Fig 3: Bar chart with error bars ──────────────────────────────────────────
fig3 <- ggplot(desc_stats, aes(x = time, y = mean, fill = diet)) +
  geom_col(position = position_dodge(0.75), width = 0.65, colour = "white", alpha = 0.9) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.75), width = 0.22, linewidth = 0.8) +
  scale_fill_manual(values = pal, name = "Diet group") +
  labs(
    title    = "Figure 3: Mean Outcome by Time and Diet Group",
    subtitle = "Bars = Mean  |  Error bars = ± 1 SE",
    x = "Time point",
    y = "Outcome (Mean ± SE)"
  ) +
  theme_pubr(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey50", size = 10)
  )

print(fig3)

# ── Fig 4: Residual diagnostic plots (LMM) ────────────────────────────────────
data_long$fitted_vals <- fitted(lme_fit)
data_long$resid_std   <- resid(lme_fit, type = "normalized")

p_rv_f <- ggplot(data_long, aes(x = fitted_vals, y = resid_std)) +
  geom_point(alpha = 0.5, colour = "#1A2E4A", size = 1.8) +
  geom_hline(yintercept = 0, colour = "#0E9E8E", linewidth = 1, linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, colour = "#DC2626", linewidth = 1) +
  labs(title = "Residuals vs Fitted",
       subtitle = "Should show random scatter around 0",
       x = "Fitted values", y = "Standardised residuals") +
  theme_pubr(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

p_qq2 <- ggplot(data_long, aes(sample = resid_std)) +
  stat_qq(colour = "#1A2E4A", size = 1.8, alpha = 0.6) +
  stat_qq_line(colour = "#0E9E8E", linewidth = 1) +
  labs(title    = "Normal Q-Q (Residuals)",
       subtitle = "Points should follow diagonal",
       x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_pubr(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

p_scale <- ggplot(data_long, aes(x = fitted_vals, y = sqrt(abs(resid_std)))) +
  geom_point(alpha = 0.5, colour = "#1A2E4A", size = 1.8) +
  geom_smooth(method = "loess", se = FALSE, colour = "#DC2626", linewidth = 1) +
  labs(title    = "Scale-Location",
       subtitle = "Flat line = homoscedasticity",
       x = "Fitted values", y = "√|Standardised residuals|") +
  theme_pubr(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

p_hist2 <- ggplot(data_long, aes(x = resid_std)) +
  geom_histogram(bins = 25, fill = "#1A2E4A", colour = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, colour = "#0E9E8E", linewidth = 1, linetype = "dashed") +
  labs(title    = "Histogram of Residuals",
       subtitle = "Should be approximately normal",
       x = "Standardised residuals", y = "Count") +
  theme_pubr(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

fig4 <- (p_rv_f | p_qq2) / (p_scale | p_hist2) +
  plot_annotation(
    title   = "Figure 4: LMM Residual Diagnostics",
    theme   = theme(plot.title = element_text(face = "bold", size = 14))
  )

print(fig4)

# ── Fig 5: Boxplots per time × diet (final clean version) ─────────────────────
fig5 <- ggplot(data_long, aes(x = time, y = outcome, fill = diet)) +
  geom_boxplot(alpha = 0.8, outlier.colour = "#DC2626", outlier.size = 2,
               position = position_dodge(0.78), width = 0.65, colour = "grey30") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5,
               position = position_dodge(0.78),
               aes(group = diet), fill = "white", colour = "grey20") +
  scale_fill_manual(values = pal, name = "Diet group") +
  labs(
    title    = "Figure 5: Distribution of Outcome by Time and Diet Group",
    subtitle = "Boxes = IQR  |  Diamonds = group mean  |  Red = outliers",
    x = "Time point",
    y = "Outcome"
  ) +
  theme_pubr(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(colour = "grey50", size = 10)
  )

print(fig5)

# -- Save all figures ----------------------------------------------------------
# output_dir: folder where the PDF will be saved.
# Default = same folder as this .R script file.
# To save elsewhere, replace with a fixed path, e.g.:
#   output_dir <- "C:/Users/YourName/Documents"   # Windows
#   output_dir <- "/Users/YourName/Documents"      # Mac / Linux

output_dir <- "~/Desktop/Medical Statistic/R MPH IIUM/stat7242"

pdf_path <- file.path(output_dir, "rmanova_figures.pdf")

pdf(pdf_path, width = 10, height = 7, onefile = TRUE)
print(fig1)
print(fig2)
print(fig3)
print(fig4)
print(fig5)
print(p_qq)
print(p_box)
dev.off()

cat(paste0("\n[OK] All figures saved to:\n     ", pdf_path, "\n"))


# --- 9. RESULTS SUMMARY -------------------------------------------------------

cat("\n\n===============================================\n")
cat("  SECTION 9: RESULTS SUMMARY\n")
cat("===============================================\n")

anova_tbl <- rm_model$ANOVA
mauchly   <- rm_model$`Mauchly's Test for Sphericity`
sph_corr  <- rm_model$`Sphericity Corrections`

cat("\n-- Recommended Reporting Sentences -------------\n\n")

for (i in 1:nrow(anova_tbl)) {
  
  eff   <- anova_tbl$Effect[i]
  F_val <- round(anova_tbl$F[i], 2)
  dfn   <- round(anova_tbl$DFn[i], 1)
  dfd   <- round(anova_tbl$DFd[i], 1)
  p_val <- anova_tbl$p[i]
  ges_v <- round(anova_tbl$ges[i], 3)
  sig   <- ifelse(p_val < 0.05,
                  "statistically significant",
                  "not statistically significant")
  
  # Check whether sphericity correction is needed for this effect
  sph_row <- which(mauchly$Effect == eff)
  
  if (length(sph_row) > 0 && mauchly$p[sph_row] < 0.05) {
    GGe   <- round(sph_corr$GGe[which(sph_corr$Effect == eff)], 3)
    p_gg  <- sph_corr$`p[GG]`[which(sph_corr$Effect == eff)]
    note  <- paste0("(Greenhouse-Geisser corrected, epsilon = ", GGe, ")")
    p_rep <- p_gg
  } else {
    note  <- ""
    p_rep <- p_val
  }
  
  p_str <- ifelse(p_rep < 0.001, "p < .001", paste0("p = ", round(p_rep, 3)))
  
  ges_label <- dplyr::case_when(
    ges_v < 0.01 ~ "negligible",
    ges_v < 0.06 ~ "small",
    ges_v < 0.14 ~ "medium",
    TRUE         ~ "large"
  )
  
  cat(paste0(
    "  [", toupper(eff), "]\n",
    "  A ", sig, " effect was found for ", eff, ",\n",
    "  F(", dfn, ", ", dfd, ") = ", F_val, ", ", p_str,
    ", ges = ", ges_v, " (", ges_label, " effect size).",
    if (note != "") paste0("\n  Note: ", note) else "",
    "\n\n"
  ))
}

# Save results summary as a text file in the same output_dir
summary_path <- file.path(output_dir, "rmanova_results_summary.txt")

sink(summary_path)
cat("REPEATED MEASURES ANOVA - RESULTS SUMMARY\n")
cat(paste0("Generated: ", Sys.time(), "\n\n"))

cat("=== ANOVA TABLE ===\n")
print(rm_model$ANOVA)
cat("\n=== MAUCHLY'S TEST FOR SPHERICITY ===\n")
print(rm_model$`Mauchly's Test for Sphericity`)
cat("\n=== SPHERICITY CORRECTIONS ===\n")
print(rm_model$`Sphericity Corrections`)
cat("\n=== DESCRIPTIVE STATISTICS ===\n")
print(as.data.frame(desc_stats))
cat("\n=== POST-HOC: TIME WITHIN EACH DIET GROUP ===\n")
print(summary(posthoc_time))
cat("\n=== POST-HOC: DIET AT EACH TIME POINT ===\n")
print(summary(posthoc_diet))
cat("\n=== EFFECT SIZES (Cohen's d between groups) ===\n")
print(as.data.frame(cohens_d_table))
cat("\n=== LMM FIXED EFFECTS ===\n")
print(summary(lme_fit)$tTable)
cat("\n=== ICC ===\n")
print(performance::icc(lme_fit))
cat("\n=== SESSION INFO ===\n")
print(sessionInfo())
sink()

cat(paste0("[OK] Results summary saved to:\n     ", summary_path, "\n"))
cat("\n[DONE] Analysis complete.\n")
cat("  Check: Console (all results) | Plots pane (figures)\n")
cat(paste0("         ", pdf_path, " (all figures)\n"))
cat(paste0("         ", summary_path, " (results text file)\n"))

