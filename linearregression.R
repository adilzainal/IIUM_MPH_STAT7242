#==============================
# Install required packages
#==============================
setwd("~/Desktop/Medical Statistic/R MPH IIUM/stat7242") #change with your cd which contain the dataset
install.packages(c("car","lmtest","ggplot2","performance"))

library(car)
library(lmtest)
library(ggplot2)
library(performance)
library(foreign) 

#==============================
# A.Example Data
#==============================

# Replace with your dataset
data = read.spss("healthstatus.sav")  # most natural way to open data in R
Y <- data$sbp
X1 <- data$hba1c
X2 <- data$hcy
X3 <- data$wt
X4 <- data$smoking

# Example sbp# Example model
model <- lm(Y ~ X1 + X2 + X3, data = data)

summary(model)

confint(model, level = 0.95)
#==============================
# 1. Linearity Assumption
#==============================

# Residuals vs fitted plot
plot(model$fitted.values, model$residuals,
     xlab="Fitted Values",
     ylab="Residuals",
     main="Residuals vs Fitted")
abline(h=0, col="red")

# Component + Residual plot
crPlots(model)

# If pattern appears (non-linear relationship)
# Try transformation

model_log <- lm(log(Y) ~ X1 + X2 + X3, data=data)

summary(model_log)

# use polynomial terms
library(ggplot2)
poly_model <- lm(Y ~ poly(X, 2))
plot(X, fitted(poly_model))

#use splines
library(splines)
spline_model <- lm(Y ~ ns(X, df=5))
plot(X, fitted(spline_model))

#==============================
# 2. Independence of Errors
#==============================

# Durbin-Watson test
durbinWatsonTest(model)
dwtest(model)

# If autocorrelation exists

# Use generalized least squares or Time series model (ARIMA Lecture)
install.packages("nlme")
library(nlme)

gls_model <- gls(Y ~ X1 + X2 + X3, data=data)

summary(gls_model)
#==============================
# 3. Homoscedasticity
#==============================

# Breusch-Pagan Test
bptest(model)

# Scale-location plot
plot(model, which=3)

# If heteroscedasticity detected (p < 0.05)

# Use robust standard errors
library(sandwich)
coeftest(model, vcov = vcovHC(model, type = "HC1"))

#==============================
# 4. Normality of Residuals
#==============================

# Histogram
hist(residuals(model),
     main="Histogram of Residuals",
     xlab="Residuals")

# Q-Q Plot
qqnorm(residuals(model))
qqline(residuals(model), col="red")

# Shapiro-Wilk test
shapiro.test(residuals(model))

#if violated use robust regression or transformation  

library(robustbase)
robust_model <- lmrob(Y ~ X1, data = data)
summary(robust_model)

#==============================
# 5. Multicollinearity
#==============================

vif(model)

# VIF interpretation
# VIF > 5 = moderate multicollinearity
# VIF > 10 = serious multicollinearity

# Remove correlated variables or combine them
model_reduced <- lm(Y ~ X1 + X3, data=data)

summary(model_reduced)

# Use ridge regression
library(glmnet)
ridge_model <- glmnet(as.matrix(cbind(X1, X2, X3)), Y, alpha=0)
summary(ridge_model)
plot(ridge_model)
coef(ridge_model, s = 0.1)


#==============================
# 6. Outliers and Influential Points
#==============================

# Cook's distance
plot(model, which=4)

# Identify influential cases
cooksd <- cooks.distance(model)

influential <- as.numeric(names(cooksd)[(cooksd > 4/length(cooksd))])
influential

# Leverage plot
plot(model, which=5)

# if violated use robust regression or remove influential points and refit the model
library(MASS)
robust_model <- rlm(Y ~ X, data = data)
summary(robust_model)

# Compare models

AIC(model, model_log, robust_model)

#=============================================================
# B. Moderator analysis (categorical and continuous)
# ============================================================
# MODERATION ANALYSIS USING LINEAR REGRESSION (base R only)
# Categorical Predictor (X) × Continuous Moderator (M) → Y
# ============================================================
# Outcome (Y):   sbp      (systolic blood pressure)
# Predictor (X): smoking  (categorical)
# Moderator (M): hcy      (homocysteine, continuous)
# ============================================================

# ---- 1. LOAD DATA ------------------------------------------
df <- data

# Mean-center the continuous moderator
df$hcy_c <- df$hcy - mean(df$hcy, na.rm = TRUE)

# ---- 2. DESCRIPTIVE STATISTICS -----------------------------
cat("===== DESCRIPTIVE STATISTICS =====\n")
cat("\nSBP by smoking group:\n")
desc <- tapply(df$sbp, df$smoking, function(x)
  c(Mean = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE), N = length(x)))
print(round(do.call(rbind, desc), 2))

cat("\nHcy (centered) summary:\n")
print(round(summary(df$hcy_c), 2))

# ---- 3. FIT THREE REGRESSION MODELS ------------------------

# Step 1: Main effect of X only
model1 <- lm(sbp ~ smoking, data = df)

# Step 2: Main effects of X + M
model2 <- lm(sbp ~ smoking + hcy_c, data = df)

# Step 3: Full moderation model (X + M + X*M)
model3 <- lm(sbp ~ smoking * hcy_c, data = df)

# ---- 4. MODEL SUMMARIES ------------------------------------
cat("\n===== STEP 1: X only =====\n")
print(summary(model1))

cat("\n===== STEP 2: X + M (main effects) =====\n")
print(summary(model2))

cat("\n===== STEP 3: X + M + X*M (moderation) =====\n")
print(summary(model3))

# ---- 5. MODEL COMPARISON (R² change) -----------------------
cat("\n===== MODEL COMPARISON =====\n")
print(anova(model2, model3))   # tests whether adding interaction improves fit

r2_1 <- summary(model1)$r.squared
r2_2 <- summary(model2)$r.squared
r2_3 <- summary(model3)$r.squared

cat(sprintf("\nR² Step 1 (X only):       %.4f\n", r2_1))
cat(sprintf("R² Step 2 (X + M):        %.4f\n", r2_2))
cat(sprintf("R² Step 3 (X + M + X*M):  %.4f\n", r2_3))
cat(sprintf("ΔR² (Step 2 → Step 3):    %.4f\n", r2_3 - r2_2))

# ---- 6. INTERPRET INTERACTION ------------------------------
b      <- coef(model3)
cf     <- coef(summary(model3))
# Grab the interaction row (always the last coefficient in X*M model)
int_name <- rownames(cf)[nrow(cf)]
p_int    <- cf[int_name, "Pr(>|t|)"]

cat("\n===== COEFFICIENTS (Step 3) =====\n")
print(round(cf, 4))

cat(sprintf("\nIntercept   = %.3f  (reference smoking group at mean hcy)\n", b[1]))
cat(sprintf("smoking     = %.3f  (smoking group difference at mean hcy)\n", b[2]))
cat(sprintf("hcy_c       = %.3f  (hcy slope for reference smoking group)\n", b[3]))
cat(sprintf("Interaction = %.3f  (difference in hcy slope across smoking groups)\n", b[4]))

cat(sprintf("\nInteraction p = %.4f → ", p_int))
if (p_int < .05) {
  cat("✔ Significant moderation.\n")
  cat("  The effect of training type on performance changes across motivation levels.\n")
} else {
  cat("✘ No significant moderation.\n")
}

# ---- 7. SIMPLE SLOPES (manual, no packages) ----------------
# hcy slope for each smoking group
cat("\n===== SIMPLE SLOPES =====\n")

# b[3] = hcy_c slope (reference group), b[4] = interaction term
slope_ref   <- b[3]   # hcy slope for reference group (e.g., No)
slope_other <- b[3] + b[4]   # hcy slope for other group (e.g., Yes)

smoking_levels <- levels(factor(df$smoking))
cat(sprintf("hcy slope (%s): %.3f\n", smoking_levels[1], unname(slope_ref)))
cat(sprintf("hcy slope (%s): %.3f\n", smoking_levels[2], unname(slope_other)))

# ---- 8. PREDICTED VALUES AT LOW / MEAN / HIGH of M ---------
cat("\n===== PREDICTED SBP AT KEY HCY LEVELS =====\n")
low_m  <- -1 * sd(df$hcy_c, na.rm = TRUE)
high_m <-  1 * sd(df$hcy_c, na.rm = TRUE)

new_data <- data.frame(
  smoking  = factor(rep(smoking_levels, 3), levels = smoking_levels),
  hcy_c    = rep(c(low_m, 0, high_m), each = length(smoking_levels)),
  level    = rep(c("Low (-1SD)", "Mean (0)", "High (+1SD)"), each = length(smoking_levels))
)

new_data$predicted <- round(predict(model3, newdata = new_data), 2)
new_data$hcy_c     <- round(new_data$hcy_c, 2)
print(new_data[, c("smoking", "level", "predicted")])

# ---- 9. INTERACTION PLOT -----------------------------------
plot_data       <- new_data
plot_data$level <- factor(plot_data$level,
                          levels = c("Low (-1SD)", "Mean (0)", "High (+1SD)"))

png("/Users/adilzainal/Desktop/Medical Statistic/R MPH IIUM/stat7242/moderation_lm_plot.png",
    width = 800, height = 550, res = 120)

plot(NULL,
     xlim = c(0.5, 3.5), ylim = range(plot_data$predicted) + c(-5, 5),
     xaxt = "n",
     xlab = "Hcy Level",
     ylab = "Predicted SBP",
     main = "Moderation: Smoking × Hcy → SBP")

axis(1, at = 1:3, labels = c("Low (-1SD)", "Mean", "High (+1SD)"))

colors <- c("#3498DB", "#E74C3C")
for (i in seq_along(smoking_levels)) {
  grp <- smoking_levels[i]
  sub <- plot_data[plot_data$smoking == grp, ]
  sub <- sub[order(sub$level), ]
  lines(1:3, sub$predicted, col = colors[i], lwd = 2)
  points(1:3, sub$predicted, col = colors[i], pch = 19, cex = 1.5)
}

legend("topleft", legend = smoking_levels,
       col = colors, lwd = 2, pch = 19, bty = "n")

dev.off()

# ---- 10. DIAGNOSTIC PLOTS ----------------------------------

png("/Users/adilzainal/Desktop/Medical Statistic/R MPH IIUM/stat7242/moderation_lm_diagnostics.png",
    width = 800, height = 550, res = 120)
par(mfrow = c(2, 2))
plot(model3)
dev.off()
par(mfrow = c(1, 1))

cat("\nPlots saved: moderation_lm_plot.png, moderation_lm_diagnostics.png\n")

#==============================
# C. Mediator (classic Baron & Kenny approach)
#==============================
# ============================================================
# MEDIATION ANALYSIS - BARON & KENNY (1986) APPROACH
# ============================================================
# Outcome (Y):   sbp    (systolic blood pressure)
# Predictor (X): wt (main IV)
# Mediator (M):  hcy    (homocysteine - partial mediator)
#
# Baron & Kenny 4-Step Criteria:
# Step 1: X → Y         (c path)  must be significant
# Step 2: X → M         (a path)  must be significant
# Step 3: M → Y | X     (b path)  must be significant
# Step 4: X → Y | M     (c' path) reduced = partial mediation
#                                  zero    = full mediation
# ============================================================

df <- data

# ============================================================
# STEP 1: X → Y (Total Effect)
# Does wt predict SBP?
# ============================================================
model1 <- lm(sbp ~ wt, data = df)

cat("===== STEP 1: Total Effect (X → Y) =====\n")
cat("Model: sbp ~ wt\n\n")
print(summary(model1))

c_path  <- coef(model1)["wt"]
p_step1 <- coef(summary(model1))["wt", "Pr(>|t|)"]
cat(sprintf("c path (total effect of wt on sbp) = %.4f, p = %.4f\n", c_path, p_step1))
if (p_step1 < .05) {
  cat("✔ Step 1 met: wt significantly predicts sbp.\n")
} else {
  cat("✘ Step 1 NOT met: wt does not significantly predict sbp.\n")
  cat("  Mediation analysis may not be appropriate.\n")
}

# ============================================================
# STEP 2: X → M (a path)
# Does wt predict hcy?
# ============================================================
model2 <- lm(hcy ~ wt, data = df)

cat("\n===== STEP 2: X → M (a path) =====\n")
cat("Model: hcy ~ wt\n\n")
print(summary(model2))

a_path  <- coef(model2)["wt"]
p_step2 <- coef(summary(model2))["wt", "Pr(>|t|)"]
cat(sprintf("a path (effect of wt on hcy) = %.4f, p = %.4f\n", a_path, p_step2))
if (p_step2 < .05) {
  cat("✔ Step 2 met: wt significantly predicts hcy.\n")
} else {
  cat("✘ Step 2 NOT met: wt does not significantly predict hcy.\n")
}

# ============================================================
# STEP 3 & 4: X + M → Y
# Does hcy predict SBP after controlling for wt? (b path)
# Is the effect of wt reduced after adding hcy?  (c' path)
# ============================================================
model3 <- lm(sbp ~ wt + hcy, data = df)

cat("\n===== STEP 3 & 4: X + M → Y =====\n")
cat("Model: sbp ~ wt + hcy\n\n")
print(summary(model3))

b_path    <- coef(model3)["hcy"]
c_prime   <- coef(model3)["wt"]
p_step3   <- coef(summary(model3))["hcy",    "Pr(>|t|)"]
p_step4   <- coef(summary(model3))["wt", "Pr(>|t|)"]

cat(sprintf("b path  (hcy → sbp | wt)  = %.4f, p = %.4f\n", b_path,  p_step3))
cat(sprintf("c' path (wt → sbp | hcy)  = %.4f, p = %.4f\n", c_prime, p_step4))

if (p_step3 < .05) {
  cat("✔ Step 3 met: hcy significantly predicts sbp after controlling for wt.\n")
} else {
  cat("✘ Step 3 NOT met: hcy does not significantly predict sbp.\n")
}

# ============================================================
# INDIRECT EFFECT & MEDIATION TYPE
# ============================================================
indirect <- a_path * b_path           # Sobel product-of-coefficients
direct   <- c_prime
total    <- c_path
prop_med <- indirect / total * 100    # % mediated

cat("\n===== INDIRECT EFFECT =====\n")
cat(sprintf("a × b (indirect effect) = %.4f × %.4f = %.4f\n", a_path, b_path, indirect))
cat(sprintf("Direct effect  (c')     = %.4f\n", direct))
cat(sprintf("Total effect   (c)      = %.4f\n", total))
cat(sprintf("Proportion mediated     = %.1f%%\n", prop_med))

# Mediation conclusion
cat("\n===== MEDIATION CONCLUSION =====\n")
all_steps <- all(c(p_step1, p_step2, p_step3) < .05)
if (all_steps && p_step4 < .05) {
  cat("✔ PARTIAL MEDIATION: hcy partially mediates the wt → sbp relationship.\n")
  cat("  Weight still significantly predicts sbp even after controlling for hcy.\n")
} else if (all_steps && p_step4 >= .05) {
  cat("✔ FULL MEDIATION: hcy fully mediates the wt → sbp relationship.\n")
  cat("  Weight no longer significantly predicts sbp after controlling for hcy.\n")
} else {
  cat("✘ Mediation criteria not fully met. Review steps above.\n")
}

# ============================================================
# SOBEL TEST (manual)
# Tests whether the indirect effect (a*b) is significant
# ============================================================
se_a    <- coef(summary(model2))["wt", "Std. Error"]
se_b    <- coef(summary(model3))["hcy",    "Std. Error"]
se_ab   <- sqrt(b_path^2 * se_a^2 + a_path^2 * se_b^2)   # Sobel SE
z_sobel <- indirect / se_ab
p_sobel <- 2 * (1 - pnorm(abs(z_sobel)))

cat("\n===== SOBEL TEST =====\n")
cat(sprintf("Indirect effect (a*b) = %.4f\n", indirect))
cat(sprintf("Sobel SE              = %.4f\n", se_ab))
cat(sprintf("Sobel Z               = %.4f\n", z_sobel))
cat(sprintf("p-value               = %.4f\n", p_sobel))
if (p_sobel < .05) {
  cat("✔ Sobel test significant: indirect effect is statistically significant.\n")
} else {
  cat("✘ Sobel test not significant: indirect effect is not statistically significant.\n")
}

# ============================================================
# SUMMARY TABLE
# ============================================================
cat("\n===== SUMMARY: PATH COEFFICIENTS =====\n")
summary_table <- data.frame(
  Path        = c("c  (X → Y total)",
                  "a  (X → M)",
                  "b  (M → Y | X)",
                  "c' (X → Y direct)",
                  "ab (indirect)"),
  Coefficient = round(c(c_path, a_path, b_path, c_prime, indirect), 4),
  p_value     = c(round(p_step1, 4),
                  round(p_step2, 4),
                  round(p_step3, 4),
                  round(p_step4, 4),
                  round(p_sobel, 4)),
  Significant = ifelse(c(p_step1, p_step2, p_step3, p_step4, p_sobel) < .05,
                       "Yes", "No")
)
print(summary_table, row.names = FALSE)

# ============================================================
# PATH DIAGRAM (text-based)
# ============================================================
cat("\n===== PATH DIAGRAM =====\n")
cat(sprintf("
         a = %.3f*            b = %.3f%s
wt ─────────────► hcy ─────────────► sbp
  │                                       ▲
  │         c' = %.3f%s                  │
  └───────────────────────────────────────┘
  (total effect c = %.3f%s)
\n",
            a_path, b_path, ifelse(p_step3 < .05, "*", ""),
            c_prime, ifelse(p_step4 < .05, "*", ""),
            c_path,  ifelse(p_step1 < .05, "*", "")
))

# ============================================================
# VISUALIZATION
# ============================================================

# Plot 1: Scatter — wt vs sbp (total effect)
png("/Users/adilzainal/Desktop/Medical Statistic/R MPH IIUM/stat7242/mediation_plot1_total.png",
    width = 700, height = 500, res = 120)
plot(df$wt, df$sbp,
     xlab = "Weight", ylab = "SBP",
     main = "Step 1: Total Effect — Weight → SBP",
     pch = 19, col = "#3498DB88")
abline(model1, col = "#E74C3C", lwd = 2)
legend("topleft", legend = sprintf("b = %.3f, p = %.3f", c_path, p_step1),
       bty = "n", text.col = "#E74C3C")
dev.off()

# Plot 2: Scatter — wt vs hcy (a path)
png("/Users/adilzainal/Desktop/Medical Statistic/R MPH IIUM/stat7242/mediation_plot2_apath.png",
    width = 700, height = 500, res = 120)
plot(df$wt, df$hcy,
     xlab = "Weight", ylab = "Hcy",
     main = "Step 2: a path — Weight → Hcy",
     pch = 19, col = "#2ECC7188")
abline(model2, col = "#E74C3C", lwd = 2)
legend("topleft", legend = sprintf("a = %.3f, p = %.3f", a_path, p_step2),
       bty = "n", text.col = "#E74C3C")
dev.off()

# Plot 3: Diagnostic plots for final model
png("/Users/adilzainal/Desktop/Medical Statistic/R MPH IIUM/stat7242/mediation_diagnostics.png",
    width = 900, height = 700, res = 120)
par(mfrow = c(2, 2))
plot(model3, main = "Model Diagnostics: sbp ~ wt + hcy")
dev.off()
par(mfrow = c(1, 1))

cat("Plots saved to: /Users/adilzainal/Desktop/Medical Statistic/R MPH IIUM/stat7242/\n")
