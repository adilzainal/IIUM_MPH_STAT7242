# =================================================================
# Survival Analysis Practical Session
# Topics: Kaplan-Meier, Life Tables, Log-Rank, and Cox Regression
# =================================================================

# 1. Setup and Data Preparation
# -----------------------------------------------------------------
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")
install.packages("survival")
library(survival)
library(survminer)

# Using the NCCTG Lung Cancer Data
df <- survival::lung
lung <- df
data(lung)

# Clean data: recode status (1=censored, 2=dead) to (0=censored, 1=dead)
lung$status <- lung$status - 1

# 2. Life Table (Actuarial Table)
# -----------------------------------------------------------------
# Life tables summarize survival data in fixed time intervals
life_tab <- with(lung, survival::survfit(Surv(time, status) ~ 1))
summary(life_tab, times = seq(0, 1000, by = 100))


# 3. Kaplan-Meier (KM) Estimation
# -----------------------------------------------------------------
# Create survival object (Time, Event)
km_fit <- survfit(Surv(time, status) ~ sex, data = lung)

# Plotting the KM Curves
ggsurvplot(km_fit, 
           data = lung, 
           risk.table = TRUE, 
           pval = TRUE, 
           conf.int = TRUE,
           main = "Survival Curves by Sex",
           legend.labs = c("Male", "Female"))


# 4. Log-Rank Test
# -----------------------------------------------------------------
# Tests the null hypothesis that there is no difference in survival between groups
log_rank <- survdiff(Surv(time, status) ~ sex, data = lung)
print(log_rank)

# Interpretation: If p-value < 0.05, survival distributions are significantly different.


# 5. Cox Proportional Hazards Regression
# -----------------------------------------------------------------

# A. Simple Cox Regression (Univariate)
univ_cox <- coxph(Surv(time, status) ~ sex, data = lung)
summary(univ_cox)

# B. Multiple Cox Regression (Multivariate)
multi_cox <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data = lung)
summary(multi_cox)

# --- INTERPRETATION GUIDE ---
# Hazard Ratio (HR) = exp(coef). 
# HR > 1: Increased risk of death (decreased survival).
# HR < 1: Decreased risk of death (protective effect).
# Example: If HR for 'sex' is 0.6, females (coded 2) have 40% lower risk than males.


# 6. Checking Cox Assumptions (Proportional Hazards)
# -----------------------------------------------------------------
# The main assumption: The hazard ratio is constant over time.

# Schoenfeld Residuals Test,If p value not sig the assumption holds; the hazard ratio is constant over time.
test_ph <- cox.zph(multi_cox)
print(test_ph)

# Visual check: If the line is roughly horizontal, the assumption holds.
ggcoxzph(test_ph)

# 7. Visualizing Hazard Ratios (Forest Plot)
# -----------------------------------------------------------------
ggforest(multi_cox, data = lung)
