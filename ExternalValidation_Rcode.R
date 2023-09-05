# Load libraries and data ---

# Libraries ----
library(lubridate)
library(dplyr)
library(ggplot2)
library(qwraps2)
options(qwraps2_markup = "markdown")
library(survival)
library(survminer)
library(riskRegression)
library(grid)
library(gridExtra)
library(geepack)
library(kableExtra)
library(boot)
library(Hmisc)
library(pseudo)
library(pec)
library(margins)
library(rms)
library(lmtest)
library(pmsampsize)

# Load datasets
load("kfre_data.Rda")
k <- kfre_data


k$beta.sum <- (-0.2201 * ((k$age/10)-7.036)) + ((0.2467*k$male) - 0.5642) - (0.5567 * ((k$epi/5)-7.222)) +
  (0.4510 * (log(k$acr_mgmmol/0.113)-5.137))

k$kfre <- 1 - 0.9570^exp(k$beta.sum)

# Create subset of South Asian data
sa <- subset(k, k$white_sa == 1)
# Create subset of White data
w <- subset(k, k$white_sa == 0)

#-------------------------------------------------------------------------------
# Sample size calculation for model development
#-------------------------------------------------------------------------------

# Calculate R squared based on: (https://www.bmj.com/content/bmj/368/bmj.m441.full.pdf)

# Royston's D
D <- 5.5(0.926-0.5)+10.26(0.926-0.5)^3 # using Ruperts C statistic
# = 3.136188
# R^2_D
R_d <- ((pi/8)*D^2) / ((pi^2)/6 + (pi/8)*D)
# = 0.7013226519
# R^2_OQ
R_OQ <- ((-pi^2/6)*R_d) / ((1 - (pi^2/6))*R_d - 1)
# = 0.7943428113
# Likelihood ratio
LR <- -429*(1-R_OQ)
# 678.4826458
# R^2_Cox-Snell
R_CS <- 1-exp(-LR/35539)
# = 0.018910130373

# Run sample size
pmsampsize("s", rsquared=0.0189, parameters=9, shrinkage=0.9, rate=0.009, timepoint=5, meanfup=5.5)
# 4241 min sample size

#-------------------------------------------------------------------------------
# Summary statistics table
#-------------------------------------------------------------------------------

# By ethnicity
summ1 <-
  list("Age" =
         list("min"       = ~ min(age),
              "max"       = ~ max(age),
              "mean (sd)" = ~ qwraps2::mean_sd(age)),
       "eGFR" =
         list("min"       = ~ min(epi),
              "median (IQR)" = ~ qwraps2::median_iqr(epi),
              "max"       = ~ max(epi),
              "mean (sd)" = ~ qwraps2::mean_sd(epi)),
       "ACR" =
         list("min"       = ~ min(acr_mgmmol),
              "median (IQR)" = ~ qwraps2::median_iqr(acr_mgmmol),
              "max"      = ~ max(acr_mgmmol),
              "mean (sd)" = ~ qwraps2::mean_sd(acr_mgmmol)),
       "Sex" =
         list("Male"      = ~ qwraps2::n_perc0(female == 0),
              "Female"    = ~ qwraps2::n_perc0(female == 1)),
       "Diabetes" = 
         list("Yes"       = ~ qwraps2::n_perc0(dm == 1)),
       "Heart failure" = 
         list("Yes"       = ~ qwraps2::n_perc0(hf == 1)),
       "Cardiovascular disease" = 
         list("Yes"       = ~ qwraps2::n_perc0(cvd == 1)),
       "Hypertension" = 
         list("Yes"       = ~ qwraps2::n_perc0(htn == 1)),
       "Follow-up time" = 
         list("mean"      = ~ qwraps2::mean_sd(time),
              "median"    = ~ qwraps2::median_iqr(time)),
       "Time to ESKD" =
         list("mean"      = ~ qwraps2::mean_sd(time[rrt == 1]),
              "median"    = ~ qwraps2::median_iqr(time[rrt == 1]))
  )

by_ethnic <- summary_table(dplyr::group_by(k, white_sa), summ1)
by_ethnic

## Death rates
sum(w$death)/sum(w$time)*1000 # 59.6
sum(sa$death)/sum(sa$time)*1000 # 31.0
# number of deaths
table(w$death) # 7625 deaths
table(sa$death) # 479 deaths
# number of deaths within 5 years
table(w$death[w$time<=5]) # 5421
table(sa$death[sa$time<=5]) # 220
# proportion of deaths within 5 years
sum(w$death[w$time<=5])/n_distinct(w$id) # 20%
sum(sa$death[sa$time<=5])/n_distinct(sa$id) # 8.06%

## KRT rates
sum(w$rrt)/sum(w$time)*1000 # 2.95
sum(sa$rrt)/sum(sa$time)*1000 # 9.20
# number of KRT events
table(w$rrt) # 378 events
table(sa$rrt) # 142 events
# number of KRT events within 5 years
table(w$rrt[w$time<=5]) # 290
table(sa$rrt[sa$time<=5]) # 104
# proportion of KRT events within 5 years
sum(w$rrt[w$time<=5])/n_distinct(w$id) # 1.07%
sum(sa$rrt[sa$time<=5])/n_distinct(sa$id) # 3.81%

## KRT rates for events within 5 years
k_5 <- subset(k, k$time<=5)
table(k_5$rrt)
sum(k_5$rrt)/sum(k_5$time)

## Death rates for events within 5 years
sum(k_5$death)/sum(k_5$time)


# Number censored within 5 years
n_distinct(k$id[k$rrt==0 & k$time<5]) #18554
# how many were deaths
n_distinct(k$id[k$death==1 & k$time<5]) # 6692
n_distinct(k$id[k$death==1 & k$time<5])/
  n_distinct(k$id[k$rrt==0 & k$time<5]) # 36.1%

#number with fu time>5
w_5 <- subset(w, w$time>=5) # 12624
sa_5 <- subset(sa, sa$time>=5) # 1614

rm(by_ethnic, sa_5, w_5, summ1)

#-------------------------------------------------------------------------------
# KM plots by risk group (for each cohort)
#-------------------------------------------------------------------------------

k$riskgrp <- cut(k$kfre, breaks=c(0, 0.03, 0.05, 0.15, 0.25, 0.5, 1), 
                 labels=c("<3%", "3-<5%", "5-<15%", "15-<25%", "25-<50%", ">=50%"),
                 right = FALSE)
sa$riskgrp <- cut(sa$kfre, breaks=c(0, 0.03, 0.05, 0.15, 0.25, 0.5, 1), 
                  labels=c("<3%", "3-<5%", "5-<15%", "15-<25%", "25-<50%", ">=50%"),
                  right = FALSE)
w$riskgrp <- cut(w$kfre, breaks=c(0, 0.03, 0.05, 0.15, 0.25, 0.5, 1), 
                 labels=c("<3%", "3-<5%", "5-<15%", "15-<25%", "25-<50%", ">=50%"),
                 right = FALSE)

surv.sa <- survfit(Surv(time, rrt) ~ riskgrp, data = sa, conf.type="log-log")
g1 <- ggsurvplot(surv.sa, 
                 conf.int=TRUE,
                 censor = TRUE,
                 censor.shape = 124, censor.size = 1,
                 xlab = "Time (years)",
                 xlim = c(0,5),
                 break.x.by=1,
                 ylim = c(0,1),
                 ylab = "ESKD-free probability",          
                 risk.table = TRUE, 
                 risk.table.y.text.col = TRUE, 
                 risk.table.fontsize = 3.4,
                 tables.height = 0.3,
                 legend.labs = c("<3%", "3 to <5%", "5 to <15%", "15 to <25%", "25 to <50%", "\u226550%"),
                 legend.title = "Risk groups",
                 title="ESKD-free probability - South Asian",
                 palette = "Spectral")

surv.w <- survfit(Surv(time, rrt) ~ riskgrp, data = w, conf.type="log-log")
g2 <- ggsurvplot(surv.w, 
                 conf.int=TRUE,
                 censor = TRUE,
                 censor.shape = 124, censor.size = 1,
                 xlab = "Time (years)",
                 xlim = c(0,5),
                 break.x.by=1,
                 ylim = c(0,1),
                 ylab = "ESKD-free probability",
                 risk.table = TRUE, 
                 risk.table.y.text.col = TRUE, 
                 risk.table.fontsize = 3.4,
                 tables.height = 0.3,
                 legend.labs = c("<3%", "3 to <5%", "5 to <15%", "15 to <25%", "25 to <50%", "\u226550%"),
                 legend.title = "Risk groups",
                 title="ESKD-free probability - White",
                 palette = "Spectral")

# Save plots as pdf
g3 <- arrange_ggsurvplots(x=list(g2, g1), ncol=2)
ggsave(file="SFig 2.pdf", g3, width=13, height=7)

## Average risk in each risk group

# KRT risk
obj.k <- summary(survfit(Surv(time, rrt) ~ riskgrp, data = k), times=5, extend=TRUE)
1 - obj.k$surv

# Death risk
obj.k1 <- summary(survfit(Surv(time, death) ~ riskgrp, data = k), times=5, extend=T)
1 - obj.k1$surv

rm(g1, g2, g3, surv.sa, surv.w, obj.k, obj.k1)

#-------------------------------------------------------------------------------
# Validation of NA KFRE:
# Calibration plots
#-------------------------------------------------------------------------------

### Decile plot

# SA cohort
# Split into deciles
sa$risk <- as.factor(ntile(sa$kfre, 10))

# Predicted risks
pred.sa <- sa %>%
  group_by(risk) %>%
  summarise_at(vars(kfre), list(mean = mean))

# KM observed risks and confidence interval
cal.sa <- summary(survfit(Surv(time, rrt) ~ risk, data = sa, conf.type="log-log"),
                  extend = TRUE, times = 5)
risk1 <- 1 - cal.sa$surv
risk.u <- 1 - cal.sa$lower
risk.l <- 1 - cal.sa$upper

# plot predicted versus observed risk 
calib.sa <- data.frame(pred.sa$mean, risk1, risk.l, risk.u)
sa1 <- ggplot(calib.sa, aes(pred.sa.mean, risk1)) + 
  geom_errorbar(aes(ymin = risk.l, ymax = risk.u, width = 0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "South Asian cohort", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.5)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.5)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

# The distribution plot        
sa2 <- ggplot(sa, aes(x = kfre)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.5)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# White cohort

# Split into deciles
w$risk <- as.factor(ntile(w$kfre, 10))

# Predicted risks
pred.w <- w %>%
  group_by(risk) %>%
  summarise_at(vars(kfre), list(mean = mean))

# KM observed risks and conf int
cal.w <- summary(survfit(Surv(time, rrt) ~ risk, data = w, conf.type="log-log"),
                 extend = TRUE, times = 5)
risk.w <- 1 - cal.w$surv
risk.u.w <- 1 - cal.w$lower
risk.l.w <- 1 - cal.w$upper

# plot predicted versus observed risk 
calib.w <- data.frame(pred.w$mean, risk.w, risk.l.w, risk.u.w)
w1 <- ggplot(calib.w, aes(pred.w.mean, risk.w)) + 
  geom_errorbar(aes(ymin = risk.l.w, ymax = risk.u.w, width = 0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "White cohort", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.5)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.5)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

# The distribution plot        
w2 <- ggplot(w, aes(x = kfre)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.5)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# Plots at 50% scale
sa3 <- ggplot(calib.sa, aes(pred.sa.mean, risk1)) + 
  geom_errorbar(aes(ymin = risk.l, ymax = risk.u, width = 0.003), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "South Asian cohort: rescaled axis", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.05)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.05)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

sa4 <- ggplot(sa, aes(x = kfre)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.05)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))


w3 <- ggplot(calib.w, aes(pred.w.mean, risk.w)) + 
  geom_errorbar(aes(ymin = risk.l.w, ymax = risk.u.w, width = 0.003), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "White cohort: rescaled axis", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.05)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.05)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

# The distribution plot        
w4 <- ggplot(w, aes(x = kfre)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.05)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))


# Combine graphs  
cal.plot <- arrangeGrob(w1, sa1, w2, sa2, w3, sa3, w4, sa4, respect=TRUE, heights=c(1,0.3,1,0.3),
                        ncol=2, nrow=4)
grid.newpage()
grid.draw(cal.plot)
ggsave("Figure 1.pdf", cal.plot, width=30, height=30, units="cm")

#-------------------------------------------------------------------------------
# Plot using pseudo observations
#-------------------------------------------------------------------------------

pdf(file="Figure 2.pdf", width=15, height=6, 
    title="Calibration plot using pseudo-observations")
par(mfrow=c(1,2))

# W #
score.w <- Score(list(w$kfre), formula = Surv(time, rrt) ~ 1, data=w, plots="cal", metrics=NULL,
                 cens.model="km", conf.int=TRUE, times=5)
plotCalibration(score.w, cens.method="pseudo", round=FALSE, xlim=c(0,1), 
                ylim=c(0,1), xlab="Predicted risk")
title("White cohort")

# SA #
score.sa <- Score(list(sa$kfre), formula = Surv(time, rrt) ~ 1, data=sa, plots="cal", metrics=NULL,
                  cens.model="km", conf.int=TRUE, times=5)
plotCalibration(score.sa, cens.method="pseudo", round=FALSE, xlim=c(0,1), 
                ylim=c(0,1), xlab="Predicted risk")
title("South Asian cohort")
dev.off()


rm("cal.plot", "cal.sa", "cal.w", "calib.sa", "calib.w", "pred.sa", "pred.w",
   "w1", "sa1", "sa2", "w2", "score.w", "score.sa", "risk.w", "risk.u.w",
   "risk.l.w", "risk1", "risk.u", "risk.l")

#-------------------------------------------------------------------------------
# O/E ratio
#-------------------------------------------------------------------------------

# SA 
# observed risk of rrt at 5 years
obj <- summary(survfit(Surv(time, rrt) ~ 1, data = sa), times = 5, extend = TRUE)
# O/E ratio 
OE <- (1 - obj$surv)/mean(sa$kfre) 

# WHITE 
obj1 <- summary(survfit(Surv(time, rrt) ~ 1, data = w), times = 5, extend = TRUE)
# O/E ratio 
OE1 <- (1 - obj1$surv)/mean(w$kfre) 

# Confidence interval
alpha <- 0.05
OE_summary <- cbind(
  "OE.SA" = OE,
  "Lower .95" = exp(log(OE - qnorm(1 - alpha / 2) * obj$std.err / (1-obj$surv))),
  "Upper .95" = exp(log(OE + qnorm(1 - alpha / 2) * obj$std.err / (1-obj$surv))),
  "OE.White" = OE1,
  "Lower .95" = exp(log(OE1 - qnorm(1 - alpha / 2) * obj1$std.err / (1-obj1$surv))),
  "Upper .95" = exp(log(OE1 + qnorm(1 - alpha / 2) * obj1$std.err / (1-obj1$surv)))
)
OE_summary <- round(OE_summary, 3)

rm("obj", "obj1", "OE_summary", "alpha", "OE", "OE1")

#-------------------------------------------------------------------------------
# Intercept and slope
#-------------------------------------------------------------------------------

# SA 
# log cumulative hazard
sa$cumhaz <- log(-log(1-sa$kfre))
# Pseudo-observations
km.sa <- survfit(Surv(time, rrt) ~ 1, data=sa)
sa$yhat <- 1 - pseudo(km.sa, times=5, type="surv")
summary(sa$yhat)

# model for calibration intercept
sa_cal_int <- geese(yhat ~ offset(cumhaz), data = sa, scale.fix = TRUE, 
                    family = gaussian, mean.link = "cloglog",
                    corstr = "independence", jack = TRUE)
# model for calibration slope
sa_cal_slope <- geese(yhat ~ offset(cumhaz) + cumhaz, data = sa, scale.fix = TRUE, 
                      family = gaussian, mean.link = "cloglog",
                      corstr = "independence", jack = TRUE)

res_cal <- rbind(
  # Value, confidence interval and test for calibration intercept
  "Intercept" = with(
    summary(sa_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  # Value, confidence interval and test for calibration slope
  "Slope" = with(
    summary(sa_cal_slope)$mean["cumhaz", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  )
)
res_cal <- round(res_cal, 3)
kable(res_cal) |>
  kable_styling("striped", position = "center")


# WHITE
# log cumulative hazard
w$cumhaz <- log(-log(1-w$kfre))
# Pseudo-observations
km.w <- survfit(Surv(time, rrt) ~ 1, data=w)
w$yhat <- 1 - pseudo(km.w, times=5, type="surv")

# model for calibration intercept
w_cal_int <- geese(yhat ~ offset(cumhaz), data = w, scale.fix = TRUE, 
                   family = gaussian, mean.link = "cloglog",
                   corstr = "independence", jack = TRUE)
# model for calibration slope
w_cal_slope <- geese(yhat ~ offset(cumhaz) + cumhaz, data = w, scale.fix = TRUE, 
                     family = gaussian, mean.link = "cloglog",
                     corstr = "independence", jack = TRUE)

res_cal.w <- rbind(
  # Value, confidence interval and test for calibration intercept
  "Intercept" = with(
    summary(w_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  # Value, confidence interval and test for calibration slope
  "Slope" = with(
    summary(w_cal_slope)$mean["cumhaz", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  )
)
res_cal.w <- round(res_cal.w, 3)
kable(res_cal.w) |>
  kable_styling("striped", position = "center")

rm("km.sa", "km.w", "res_cal", "res_cal.w", "w_cal_int", "w_cal_slope", "sa_cal_slope", "sa_cal_int")

#-------------------------------------------------------------------------------
# Discrimination - C-index
#-------------------------------------------------------------------------------

# Create truncated dataset (5 years cut off)
sa.trunc <- sa
w.trunc <- w
sa.trunc$rrt <- ifelse(sa.trunc$time>5, 0, sa.trunc$rrt)
w.trunc$rrt <- ifelse(w.trunc$time>5, 0, w.trunc$rrt)
sa.trunc$time <- ifelse(sa.trunc$time>5, 5.001, sa.trunc$time)
w.trunc$time <- ifelse(w.trunc$time>5, 5.001, w.trunc$time)

# c index
rcorr.cens(-1*sa.trunc$beta.sum, Surv(sa.trunc$time, sa.trunc$rrt))[1]
rcorr.cens(-1*w.trunc$beta.sum, Surv(w.trunc$time, w.trunc$rrt))[1]

# Bootstrap c-index to get confidence interval:
set.seed(1234)
C_boot <- function(data, i) {
  val <- data[i,] # select sample with boot
  rcorr.cens(-1*val$beta.sum, Surv(val$time, val$rrt))[1]
}
set.seed(1234)
saboot <- boot(sa.trunc, C_boot, R=1000)
boot.ci(boot.out=saboot, type="norm")
set.seed(1234)
wboot <- boot(w.trunc, C_boot, R=1000)
boot.ci(boot.out=wboot, type="norm")

rm("sa.trunc", "w.trunc", "saboot", "wboot", "modmod.w", "modmod.sa")

#-------------------------------------------------------------------------------
# Overall fit - Brier score
#-------------------------------------------------------------------------------

# Brier Score
brier.sa <- Score(list(sa$kfre), formula = Surv(time, rrt) ~ 1, data=sa, metrics=c("auc", "Brier"),
                  conf.int=TRUE, times=5, summary="ipa")
brier.sa
brier.w <- Score(list(w$kfre), formula = Surv(time, rrt) ~ 1, data=w, metrics=c("auc", "Brier"),
                 conf.int=TRUE, times=5, summary="ipa")
brier.w

# Scaled Brier Score - bootstrap for 95% CI
B_boot <- function(data, i) {
  val <- data[i,]
  brier <- Score(list(val$kfre), formula = Surv(time, rrt) ~ 1, data=val, metrics="Brier",
                 cens.model="km", conf.int=TRUE, times=5, summary="ipa")
  scaled <- brier$Brier$score$IPA
}
bboot.sa <- boot(sa, B_boot, R=1000)
boot.ci(boot.out=bboot.sa, index=2, type="perc")

bboot.w <- boot(w, B_boot, R=1000)
boot.ci(boot.out=bboot.w, index=2, type="perc")

rm("brier.sa", "brier.w", "bboot.sa", "bboot.w")

#-------------------------------------------------------------------------------
# Model updating
#-------------------------------------------------------------------------------

### Model 2 (a and b)
# find baseline risk
coxsa2 <- coxph(Surv(time, rrt) ~ offset(beta.sum), data = sa, x=TRUE, y=TRUE)
coxw2 <- coxph(Surv(time, rrt) ~ offset(beta.sum), data = w, x=TRUE, y=TRUE)

p1 <- predictSurvProb(coxsa2, newdata=data.frame(beta.sum=0), times = 5,
                      type="survival", confint=TRUE, se=TRUE) 
p2 <- predictSurvProb(coxw2, newdata=data.frame(beta.sum=0), times = 5,
                      type="survival", confint=TRUE, se=TRUE)

# Predictions for model 2
sa$kfre2 <- 1 - p1[1,1]^exp(sa$beta.sum)
w$kfre2 <- 1 - p2[1,1]^exp(w$beta.sum)


### Model 3 (a and b)
# fit cox to find slope
coxsa3 <- coxph(Surv(time, rrt) ~ beta.sum, data = sa, x=TRUE, y=TRUE)
coxw3 <- coxph(Surv(time, rrt) ~ beta.sum, data = w, x=TRUE, y=TRUE)
# slope
coxsa3$coefficients
coxw3$coefficients

# set beta.sum3 as new lp using slope
sa$beta.sum3 <- coxsa3$coefficients * sa$beta.sum
w$beta.sum3 <- coxw3$coefficients * w$beta.sum

# find baseline risk
coxsa3.1 <- coxph(Surv(time, rrt) ~ offset(beta.sum3), data=sa, x=T, y=T)
p3 <- predictSurvProb(coxsa3.1, newdata=data.frame(beta.sum3=0), times = 5,
                      type="survival", confint=TRUE, se=TRUE) 

coxw3.1 <- coxph(Surv(time, rrt) ~ offset(beta.sum3), data=w, x=T, y=T)
p4 <- predictSurvProb(coxw3.1, newdata=data.frame(beta.sum3=0), times = 5,
                      type="survival", confint=TRUE, se=TRUE) 

# predicted values for model 3
sa$kfre3 <- 1 - p3[1,1]^exp(sa$beta.sum3)
w$kfre3 <- 1 - p4[1,1]^exp(w$beta.sum3)


### Model 4 

# find ethnicity coefficient
coxk4 <- coxph(Surv(time, rrt) ~ offset(beta.sum) + white_sa, data = k, x=TRUE, y=TRUE)
confint(coxk4)
# save as new beta sum
k$beta.sum_eth <- k$beta.sum + coxk4$coefficients*k$white_sa

# refit cox to find slope
coxk4.1 <- coxph(Surv(time, rrt) ~ beta.sum_eth, data = k, x=TRUE, y=TRUE)
# slope
coxk4.1$coefficients
# set beta.sum4 as new lp using slope
k$beta.sum4 <- coxk4.1$coefficients*k$beta.sum_eth

# find baseline risk
coxk4.2 <- coxph(Surv(time, rrt) ~ offset(beta.sum4), data=k, x=T, y=T)
p5 <- predictSurvProb(coxk4.2, newdata=data.frame(beta.sum4=0), times = 5,
                      type="survival", confint=TRUE, se=TRUE)

# predicted values for model 4
k$kfre4 <- 1 - p5[1,1]^exp(k$beta.sum4)

rm("coxsa2","coxw2","coxk4","coxk4.1","coxk4.2","coxsa3","coxsa3.1","coxw3",
   "coxw3.1", "p1", "p2", "p3", "p4", "p5")

## Add kfre2 and kfre3 (i.e. models 2a, 2b, 3a, 3b to overall dataset k)
k <- select(w, id, kfre2, kfre3) %>% 
  rbind(select(sa, id, kfre2, kfre3)) %>% 
  merge(x=k, y=., by="id")

## Add kfre4 to datasets by ethnicity (w and sa)
w <- k %>% 
  filter(white_sa==0) %>% 
  select(id, kfre4) %>% 
  merge(x=w, y=., by="id")

sa <- k %>% 
  filter(white_sa==1) %>% 
  select(id, kfre4) %>% 
  merge(x=sa, y=., by="id")

#-------------------------------------------------------------------------------
# Developing a new model
#-------------------------------------------------------------------------------

# Note: 'predictor'.c is the centred value of each predictor

# Interaction effects

m1 <- coxph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa, data=k)
m2 <- coxph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa + white_sa:male.c, data=k)
m3 <- coxph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa + white_sa:acr.c, data=k)
m4 <- coxph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa + white_sa:epi.c, data=k)
m5 <- coxph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa + white_sa:age.c, data=k)

# Likelihood ratio test

lrtest(m1, m2)
lrtest(m1, m3) # lowest p value 0.039
lrtest(m1, m4)
lrtest(m1, m5)

m6 <- coxph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa + white_sa:acr.c + white_sa:age.c, data=k)
m7 <- coxph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa + white_sa:acr.c + white_sa:male.c, data=k)
m8 <- coxph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa + white_sa:acr.c + white_sa:epi.c, data=k)
lrtest(m3, m6)
lrtest(m3, m7)
lrtest(m3, m8) # best model
summary(m8)

m9 <- coxph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa + white_sa:acr.c + white_sa:epi.c +
              white_sa:age.c, data=k, x=TRUE)
m10 <- coxph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa + white_sa:acr.c + white_sa:epi.c +
               white_sa:male.c, data=k, x=TRUE)

lrtest(m8, m9)
lrtest(m8, m10) # stop at m8

rm(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)

#-------------------------------------------------------------------------------
# Fit the new model
#-------------------------------------------------------------------------------
m8 <- coxph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa + white_sa:acr.c + white_sa:epi.c, data=k)

fit <- cph(Surv(time, rrt) ~ age.c + epi.c + acr.c + male.c + white_sa + white_sa:acr.c + white_sa:epi.c, 
           data=k, x=TRUE, y=TRUE)

# new prognostic index/linear predictor
k$pi <- fit$linear.predictors
# baseline 5 year survival
fit1 <- coxph(Surv(time, rrt) ~ offset(pi), data=k, x=TRUE, y=TRUE)
p6 <- predictSurvProb(fit1, newdata=data.frame(beta.sum=0, white_sa=0), times = 5,
                      type="survival", confint=TRUE, se=TRUE) # 0.9963
# predicted risks
k$outcome <- 1 - p6[1,1]^exp(k$pi)
# coefficients of apparent model
fit$coefficients
confint(fit)

# Apparent performance measures and optimism
# See 'Model_5_InternalValidation' R script
# Shrinkage factor = 1.015

# Find internally validated model
k$pi5 <- fit$linear.predictors*1.1015

# Baseline 5 year survival re-estimated
mod1 <- coxph(Surv(time, rrt) ~ offset(pi5), data=k, x=TRUE, y=TRUE)
surv5 <- predictSurvProb(mod1, newdata=data.frame(pi5=0), times = 5,
                         type="survival", confint=TRUE, se=TRUE) # 0.9975

# internally validated risks
k$kfre5 <- 1 - surv5[1,1]^exp(k$pi5)

# Remove pre-internal validation variables
k <- subset(k, select=-c(pi, outcome))

# optimism adjusted coefficients
fit$coefficients*1.1019

## Add kfre5 to datasets by ethnicity (w and sa)
w <- k %>% 
  filter(white_sa==0) %>% 
  select(id, kfre5) %>% 
  merge(x=w, y=., by="id")

sa <- k %>% 
  filter(white_sa==1) %>% 
  select(id, kfre5) %>% 
  merge(x=sa, y=., by="id")

rm(fit, fit1, m8, mod1, p6, surv5)

#-------------------------------------------------------------------------------
# Models 2-5 Calibration plots
#-------------------------------------------------------------------------------
### MODEL 2
#----- model 2 WHITE ------

# 10 equal risk groups
w$risk2 <- as.factor(ntile(w$kfre2, 10))

## predicted risk
pred.w <- w %>%
  group_by(risk2) %>%
  summarise_at(vars(kfre2), list(mean = mean))

# observed risk: 1 - survival probability
cal.w <- summary(survfit(Surv(time, rrt) ~ risk2, data = w, conf.type="log-log"),
                 times = 5, extend = TRUE)
obs.risk <- 1 - cal.w$surv
obs.risk.u <- 1 - cal.w$lower
obs.risk.l <- 1 - cal.w$upper

# plot predicted versus observed risk 
calib <- data.frame(pred.w$mean, obs.risk, obs.risk.l, obs.risk.u)
w.m2.1 <- ggplot(calib, aes(pred.w.mean, obs.risk)) + 
  geom_errorbar(aes(ymin = obs.risk.l, ymax = obs.risk.u, width=0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "White cohort", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

# The distribution plot        
w.m2.2 <- ggplot(w, aes(x = kfre2)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.6)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# Combine   
w.m2.plot <- arrangeGrob(w.m2.1, w.m2.2, respect = TRUE, heights = c(1, 0.4), ncol = 1)

#----- SA -----

# Risk groups
sa$risk2 <- as.factor(ntile(sa$kfre2, 10))

# Predicted risk by group
pred.sa <- sa %>%
  group_by(risk2) %>%
  summarise_at(vars(kfre2), list(mean = mean))

# Observed risk by group & CI
cal.sa <- summary(survfit(Surv(time, rrt) ~ risk2, data = sa, conf.type="log-log"),
                  extend = TRUE, times = 5)
risk1 <- 1 - cal.sa$surv
risk.u <- 1 - cal.sa$lower
risk.l <- 1 - cal.sa$upper

# plot predicted versus observed risk 
calib.sa <- data.frame(pred.sa$mean, risk1, risk.l, risk.u)
sa.m2.1 <- ggplot(calib.sa, aes(pred.sa.mean, risk1)) + 
  geom_errorbar(aes(ymin = risk.l, ymax = risk.u, width = 0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "South Asian cohort", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

# The distribution plot        
sa.m2.2 <- ggplot(sa, aes(x = kfre2)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.6)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# Combine  
sa.m2.plot <- arrangeGrob(sa.m2.1, sa.m2.2, respect = TRUE, heights = c(1, 0.4), ncol = 1)


#----- OVERALL -------------

# Risk groups
k$risk2 <- as.factor(ntile(k$kfre2, 10))

## predicted risk
pred.k <- k %>%
  group_by(risk2) %>%
  summarise_at(vars(kfre2), list(mean = mean))

# observed risk
cal.k <- summary(survfit(Surv(time, rrt) ~ risk1, data = k, conf.type="log-log"),
                 times = 5, extend = TRUE)
obs.risk <- 1 - cal.k$surv
obs.risk.u <- 1 - cal.k$lower
obs.risk.l <- 1 - cal.k$upper

# plot predicted versus observed risk 
calib <- data.frame(pred.k$mean, obs.risk, obs.risk.l, obs.risk.u)
k.m2.1 <- ggplot(calib, aes(pred.k.mean, obs.risk)) + 
  geom_errorbar(aes(ymin = obs.risk.l, ymax = obs.risk.u, width=0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "Overall", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::percent, limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::percent, limits=c(0, 0.6))

# The distribution plot        
k.m2.2 <- ggplot(k, aes(x = kfre2)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 0.6)) +
  xlab("Predicted Probability") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank())

# Combine  
k.m2.plot <- arrangeGrob(k.m2.1, k.m2.2, respect = TRUE, heights = c(1, 0.4), ncol = 1)

# Title
tgrob2 <- text_grob("Model 2", size=18)
plot_2 <- as_ggplot(tgrob2)
# Plot
mod2.plot <- ggarrange(plot_2, NULL, NULL, w.m2.plot, sa.m2.plot, k.m2.plot, nrow = 2, ncol = 3, heights=c(1,8))


#----- model 3 WHITE ----
w$risk3 <- as.factor(ntile(w$kfre3, 10))

pred.w <- w %>%
  group_by(risk3) %>%
  summarise_at(vars(kfre3), list(mean = mean))

cal.w <- summary(survfit(Surv(time, rrt) ~ risk3, data = w, conf.type="log-log"),
                 times = 5, extend = TRUE)
obs.risk <- 1 - cal.w$surv
obs.risk.u <- 1 - cal.w$lower
obs.risk.l <- 1 - cal.w$upper

calib <- data.frame(pred.w$mean, obs.risk, obs.risk.l, obs.risk.u)
w.m3.1 <- ggplot(calib, aes(pred.w.mean, obs.risk)) + 
  geom_errorbar(aes(ymin = obs.risk.l, ymax = obs.risk.u, width=0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "White cohort", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

w.m3.2 <- ggplot(w, aes(x = kfre3)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.6)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# Combine   
w.m3.plot <- arrangeGrob(w.m3.1, w.m3.2, respect = TRUE, heights = c(1, 0.3), ncol = 1)

#-----  SA -----
sa$risk3 <- as.factor(ntile(sa$kfre3, 10))

pred.sa <- sa %>%
  group_by(risk3) %>%
  summarise_at(vars(kfre3), list(mean = mean))

cal.sa <- summary(survfit(Surv(time, rrt) ~ risk3, data = sa, conf.type="log-log"),
                  extend = TRUE, times = 5)
risk2 <- 1 - cal.sa$surv
risk.u <- 1 - cal.sa$lower
risk.l <- 1 - cal.sa$upper

calib.sa <- data.frame(pred.sa$mean, risk3, risk.l, risk.u)
sa.m3.1 <- ggplot(calib.sa, aes(pred.sa.mean, risk2)) + 
  geom_errorbar(aes(ymin = risk.l, ymax = risk.u, width = 0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "South Asian cohort", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

sa.m3.2 <- ggplot(sa, aes(x = kfre3)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.6)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# Combine  
sa.m3.plot <- arrangeGrob(sa.m3.1, sa.m3.2, respect = TRUE, heights = c(1, 0.3), ncol = 1)


#----- OVERALL -------------

k$risk3 <- as.factor(ntile(k$kfre3, 10))

pred.k <- k %>%
  group_by(risk3) %>%
  summarise_at(vars(kfre3), list(mean = mean))

cal.k <- summary(survfit(Surv(time, rrt) ~ risk3, data = k, conf.type="log-log"),
                 times = 5, extend = TRUE)
obs.risk <- 1 - cal.k$surv
obs.risk.u <- 1 - cal.k$lower
obs.risk.l <- 1 - cal.k$upper

calib <- data.frame(pred.k$mean, obs.risk, obs.risk.l, obs.risk.u)
k.m3.1 <- ggplot(calib, aes(pred.k.mean, obs.risk)) + 
  geom_errorbar(aes(ymin = obs.risk.l, ymax = obs.risk.u, width=0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "Overall", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::percent, limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::percent, limits=c(0, 0.6))

k.m3.2 <- ggplot(k, aes(x = kfre3)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 0.6)) +
  xlab("Predicted Probability") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank())

# Combine  
k.m3.plot <- arrangeGrob(k.m3.1, k.m3.2, respect = TRUE, heights = c(1, 0.3), ncol = 1)
# Title
tgrob3 <- text_grob("Model 3", size=18)
plot_3 <- as_ggplot(tgrob3)
# Plot
mod3.plot <- ggarrange(plot_3, NULL, NULL, w.m3.plot, sa.m3.plot, k.m3.plot, nrow = 2, ncol = 3, heights=c(1,8))

### MODEL 4

#----- model 4 WHITE ------

w$risk4 <- as.factor(ntile(w$kfre4, 10))

pred.w <- w %>%
  group_by(risk4) %>%
  summarise_at(vars(kfre4), list(mean = mean))

cal.w <- summary(survfit(Surv(time, rrt) ~ risk4, data = w, conf.type="log-log"),
                 times = 5, extend = TRUE)
obs.risk <- 1 - cal.w$surv
obs.risk.u <- 1 - cal.w$lower
obs.risk.l <- 1 - cal.w$upper

calib <- data.frame(pred.w$mean, obs.risk, obs.risk.l, obs.risk.u)
w.m4.1 <- ggplot(calib, aes(pred.w.mean, obs.risk)) + 
  geom_errorbar(aes(ymin = obs.risk.l, ymax = obs.risk.u, width=0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "White cohort", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

w.m4.2 <- ggplot(w, aes(x = kfre4)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.6)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# Combine   
w.m4.plot <- arrangeGrob(w.m4.1, w.m4.2, respect = TRUE, heights = c(1, 0.3), ncol = 1)

#----- SA -----

sa$risk4 <- as.factor(ntile(sa$kfre4, 10))

pred.sa <- sa %>%
  group_by(risk4) %>%
  summarise_at(vars(kfre4), list(mean = mean))

cal.sa <- summary(survfit(Surv(time, rrt) ~ risk4, data = sa, conf.type="log-log"),
                  extend = TRUE, times = 5)
risk3 <- 1 - cal.sa$surv
risk.u <- 1 - cal.sa$lower
risk.l <- 1 - cal.sa$upper

calib.sa <- data.frame(pred.sa$mean, risk3, risk.l, risk.u)
sa.m4.1 <- ggplot(calib.sa, aes(pred.sa.mean, risk3)) + 
  geom_errorbar(aes(ymin = risk.l, ymax = risk.u, width = 0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "South Asian cohort", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

sa.m4.2 <- ggplot(sa, aes(x = kfre4)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.6)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# Combine  
sa.m4.plot <- arrangeGrob(sa.m4.1, sa.m4.2, respect = TRUE, heights = c(1, 0.3), ncol = 1)


#----- OVERALL -------------

k$risk4 <- as.factor(ntile(k$kfre4, 10))

pred.k <- k %>%
  group_by(risk4) %>%
  summarise_at(vars(kfre4), list(mean = mean))

cal.k <- summary(survfit(Surv(time, rrt) ~ risk4, data = k, conf.type="log-log"),
                 times = 5, extend = TRUE)
obs.risk <- 1 - cal.k$surv
obs.risk.u <- 1 - cal.k$lower
obs.risk.l <- 1 - cal.k$upper

calib <- data.frame(pred.k$mean, obs.risk, obs.risk.l, obs.risk.u)
k.m4.1 <- ggplot(calib, aes(pred.k.mean, obs.risk)) + 
  geom_errorbar(aes(ymin = obs.risk.l, ymax = obs.risk.u, width=0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "Overall", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::percent, limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::percent, limits=c(0, 0.6))

k.m4.2 <- ggplot(k, aes(x = kfre4)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 0.6)) +
  xlab("Predicted Probability") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank())

# Combine  
k.m4.plot <- arrangeGrob(k.m4.1, k.m4.2, respect = TRUE, heights = c(1, 0.3), ncol = 1)
# Title
tgrob4 <- text_grob("Model 4", size=18)
plot_4 <- as_ggplot(tgrob4)
# Plot
mod4.plot <- ggarrange(plot_4, NULL, NULL, w.m4.plot, sa.m4.plot, k.m4.plot, nrow = 2, ncol = 3, heights=c(1,8))

### MODEL 5

#----- model 5 WHITE ------

w$risk5 <- as.factor(ntile(w$kfre5, 10))

pred.w <- w %>%
  group_by(risk5) %>%
  summarise_at(vars(kfre5), list(mean = mean))

cal.w <- summary(survfit(Surv(time, rrt) ~ risk5, data = w, conf.type="log-log"),
                 times = 5, extend = TRUE)
obs.risk <- 1 - cal.w$surv
obs.risk.u <- 1 - cal.w$lower
obs.risk.l <- 1 - cal.w$upper

calib <- data.frame(pred.w$mean, obs.risk, obs.risk.l, obs.risk.u)
w.m5.1 <- ggplot(calib, aes(pred.w.mean, obs.risk)) + 
  geom_errorbar(aes(ymin = obs.risk.l, ymax = obs.risk.u, width=0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "White cohort", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

w.m5.2 <- ggplot(w, aes(x = kfre5)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.6)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# Combine   
w.m5.plot <- arrangeGrob(w.m5.1, w.m5.2, respect = TRUE, heights = c(1, 0.3), ncol = 1)

#----- SA -----

sa$risk5 <- as.factor(ntile(sa$kfre5, 10))

pred.sa <- sa %>%
  group_by(risk5) %>%
  summarise_at(vars(kfre5), list(mean = mean))

cal.sa <- summary(survfit(Surv(time, rrt) ~ risk5, data = sa, conf.type="log-log"),
                  extend = TRUE, times = 5)
risk4 <- 1 - cal.sa$surv
risk.u <- 1 - cal.sa$lower
risk.l <- 1 - cal.sa$upper

calib.sa <- data.frame(pred.sa$mean, risk4, risk.l, risk.u)
sa.m5.1 <- ggplot(calib.sa, aes(pred.sa.mean, risk4)) + 
  geom_errorbar(aes(ymin = risk.l, ymax = risk.u, width = 0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "South Asian cohort", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.6)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

sa.m5.2 <- ggplot(sa, aes(x = kfre5)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.6)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# Combine  
sa.m5.plot <- arrangeGrob(sa.m5.1, sa.m5.2, respect = TRUE, heights = c(1, 0.3), ncol = 1)


#----- OVERALL -------------

k$risk5 <- as.factor(ntile(k$kfre5, 10))


pred.k <- k %>%
  group_by(risk5) %>%
  summarise_at(vars(kfre5), list(mean = mean))

cal.k <- summary(survfit(Surv(time, rrt) ~ risk5, data = k, conf.type="log-log"),
                 times = 5, extend = TRUE)
obs.risk <- 1 - cal.k$surv
obs.risk.u <- 1 - cal.k$lower
obs.risk.l <- 1 - cal.k$upper

calib <- data.frame(pred.k$mean, obs.risk, obs.risk.l, obs.risk.u)
k.m5.1 <- ggplot(calib, aes(pred.k.mean, obs.risk)) + 
  geom_errorbar(aes(ymin = obs.risk.l, ymax = obs.risk.u, width=0.02), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "Overall", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::percent, limits=c(0, 0.6)) +
  scale_x_continuous(labels = scales::percent, limits=c(0, 0.6))

k.m5.2 <- ggplot(k, aes(x = kfre5)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 0.6)) +
  xlab("Predicted Probability") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank())

# Combine  
k.m5.plot <- arrangeGrob(k.m5.1, k.m5.2, respect = TRUE, heights = c(1, 0.3), ncol = 1)
# Title
tgrob5 <- text_grob("Model 5", size=18)
plot_5 <- as_ggplot(tgrob5)
# Plot
mod5.plot <- ggarrange(plot_5, NULL, NULL, w.m5.plot, sa.m5.plot, k.m5.plot, nrow = 2, ncol = 3, heights=c(1,8))


# Save plots as pdf
multi.plot <- ggarrange(mod2.plot, mod3.plot, mod4.plot, mod5.plot, nrow=4)
ggsave(file="SFig 3.pdf", multi.plot, width=15, height=24)

rm(obs.risk,obs.risk.l,obs.risk.u,risk.l,risk.u,risk1,risk2,risk3,risk4,k.m2.2,
   cal.k,cal.sa,cal.w,calib,calib.sa,k.m2.1,k.m2.1,k.m3.1,k.m3.2,k.m4.1,k.m4.2,
   k.m5.1,k.m5.2,plot_2,plot_3,plot_4,plot_5,pred.k,pred.sa,pred.w,sa.m2.1,sa.m2.2,
   sa.m3.1,sa.m3.2,sa.m4.1,sa.m4.2,sa.m5.1,sa.m5.2,tgrob2,tgrob3,tgrob4,tgrob5,
   w.m2.1,w.m2.2,w.m3.1,w.m3.2,w.m4.1,w.m4.2,w.m5.1,w.m5.2,w.m2.plot,w.m3.plot,
   w.m4.plot,w.m5.plot,sa.m2.plot,sa.m3.plot,sa.m4.plot,sa.m5.plot,k.m2.plot,
   k.m3.plot,k.m4.plot,k.m5.plot)


#-------------------------------------------------------------------------------
# Model Performance measures for models 2-5

#-------------------------------------------------------------------------------
# O/E ratio
#-------------------------------------------------------------------------------
#--- Model 2 ----
# WHITE 
obj1 <- summary(survfit(Surv(time, rrt) ~ 1, data = w), times = 5, extend = TRUE)
# O/E ratio 
OE1 <- (1 - obj1$surv)/mean(w$kfre2)

# SA 
# observed risk of rrt at 5 years
obj <- summary(survfit(Surv(time, rrt) ~ 1, data = sa), times = 5, extend = TRUE)
# O/E ratio 
OE <- (1 - obj$surv)/mean(sa$kfre2)

# Confidence interval
alpha <- 0.05
OE_summary <- cbind(
  "OE.White" = OE1,
  "Lower .95" = exp(log(OE1 - qnorm(1 - alpha / 2) * obj1$std.err / (1-obj1$surv))),
  "Upper .95" = exp(log(OE1 + qnorm(1 - alpha / 2) * obj1$std.err / (1-obj1$surv))),
  "OE.SA" = OE,
  "Lower .95" = exp(log(OE - qnorm(1 - alpha / 2) * obj$std.err / (1-obj$surv))),
  "Upper .95" = exp(log(OE + qnorm(1 - alpha / 2) * obj$std.err / (1-obj$surv)))
)
OE_summary <- round(OE_summary, 5)

rm("obj", "obj1", "OE_summary", "alpha", "OE", "OE1")

#--- Model 3 ----
# WHITE 
obj1 <- summary(survfit(Surv(time, rrt) ~ 1, data = w), times = 5, extend = TRUE)
# O/E ratio 
OE1 <- (1 - obj1$surv)/mean(w$kfre3)

# SA 
# observed risk of rrt at 5 years
obj <- summary(survfit(Surv(time, rrt) ~ 1, data = sa), times = 5, extend = TRUE)
# O/E ratio 
OE <- (1 - obj$surv)/mean(sa$kfre3)

# Confidence interval
alpha <- 0.05
OE_summary <- cbind(
  "OE.White" = OE1,
  "Lower .95" = exp(log(OE1 - qnorm(1 - alpha / 2) * obj1$std.err / (1-obj1$surv))),
  "Upper .95" = exp(log(OE1 + qnorm(1 - alpha / 2) * obj1$std.err / (1-obj1$surv))),
  "OE.SA" = OE,
  "Lower .95" = exp(log(OE - qnorm(1 - alpha / 2) * obj$std.err / (1-obj$surv))),
  "Upper .95" = exp(log(OE + qnorm(1 - alpha / 2) * obj$std.err / (1-obj$surv)))
)
OE_summary <- round(OE_summary, 5)

rm("obj", "obj1", "OE_summary", "alpha", "OE", "OE1")

#--- Model 4 ----

# WHITE
# observed risk of rrt at 5 years
obj1 <- summary(survfit(Surv(time, rrt) ~ 1, data = w), times = 5, extend = TRUE)
# O/E ratio 
OE1 <- (1 - obj1$surv)/mean(w$kfre4)

# SA
# observed risk of rrt at 5 years
obj2 <- summary(survfit(Surv(time, rrt) ~ 1, data = sa), times = 5, extend = TRUE)
# O/E ratio 
OE2 <- (1 - obj2$surv)/mean(sa$kfre4)

# Confidence interval
alpha <- 0.05
OE_summary1 <- cbind(
  "OE_W" = OE1,
  "Lower .95" = exp(log(OE1 - qnorm(1 - alpha / 2) * obj1$std.err / (1-obj1$surv))),
  "Upper .95" = exp(log(OE1 + qnorm(1 - alpha / 2) * obj1$std.err / (1-obj1$surv))),
  "OE_SA" = OE2,
  "Lower .95" = exp(log(OE2 - qnorm(1 - alpha / 2) * obj2$std.err / (1-obj2$surv))),
  "Upper .95" = exp(log(OE2 + qnorm(1 - alpha / 2) * obj2$std.err / (1-obj2$surv)))
)
OE_summary1 <- round(OE_summary1, 5)

# OVERALL
# observed risk of rrt at 5 years
obj <- summary(survfit(Surv(time, rrt) ~ 1, data = k), times = 5, extend = TRUE)
# O/E ratio 
OE <- (1 - obj$surv)/mean(k$kfre4)

# Confidence interval (Competing risks github code)
alpha <- 0.05
OE_summary <- cbind(
  "OE" = OE,
  "Lower .95" = exp(log(OE - qnorm(1 - alpha / 2) * obj$std.err / (1-obj$surv))),
  "Upper .95" = exp(log(OE + qnorm(1 - alpha / 2) * obj$std.err / (1-obj$surv)))
)
OE_summary <- round(OE_summary, 3)

rm("obj", "OE_summary", "alpha", "OE", "obj1", "OE_summary1", "OE1", "OE2")


#--- Model 5 ----
# WHITE
# observed risk of rrt at 5 years
obj1 <- summary(survfit(Surv(time, rrt) ~ 1, data = w), times = 5, extend = TRUE)
# O/E ratio 
OE1 <- (1 - obj1$surv)/mean(w$kfre5)

# SA
# observed risk of rrt at 5 years
obj2 <- summary(survfit(Surv(time, rrt) ~ 1, data = sa), times = 5, extend = TRUE)
# O/E ratio 
OE2 <- (1 - obj2$surv)/mean(sa$kfre5)

# Confidence interval
alpha <- 0.05
OE_summary1 <- cbind(
  "OE_W" = OE1,
  "Lower .95" = exp(log(OE1 - qnorm(1 - alpha / 2) * obj1$std.err / (1-obj1$surv))),
  "Upper .95" = exp(log(OE1 + qnorm(1 - alpha / 2) * obj1$std.err / (1-obj1$surv))),
  "OE_SA" = OE2,
  "Lower .95" = exp(log(OE2 - qnorm(1 - alpha / 2) * obj2$std.err / (1-obj2$surv))),
  "Upper .95" = exp(log(OE2 + qnorm(1 - alpha / 2) * obj2$std.err / (1-obj2$surv)))
)
OE_summary1 <- round(OE_summary1, 5)

# OVERALL 
# observed risk of rrt at 5 years
obj <- summary(survfit(Surv(time, rrt) ~ 1, data = k), times = 5, extend = TRUE)
# O/E ratio 
OE <- (1 - obj$surv)/mean(k$kfre5)

# Confidence interval (Competing risks github code)
alpha <- 0.05
OE_summary <- cbind(
  "OE" = OE,
  "Lower .95" = exp(log(OE - qnorm(1 - alpha / 2) * obj$std.err / (1-obj$surv))),
  "Upper .95" = exp(log(OE + qnorm(1 - alpha / 2) * obj$std.err / (1-obj$surv)))
)
OE_summary <- round(OE_summary, 3)

rm("obj", "OE_summary", "alpha", "OE", "obj1", "OE_summary1", "OE1", "OE2")


#-------------------------------------------------------------------------------
# Calibration intercept and slope
#-------------------------------------------------------------------------------
#--- Model 2 ----

# SA
# log cumulative hazard
sa$cumhaz2 <- log(-log(1-sa$kfre2))
# Pseudo-observations
km.sa <- survfit(Surv(time, rrt) ~ 1, data=sa)
sa$yhat2 <- 1 - pseudo(km.sa, times=5, type="surv")

# model for calibration intercept
sa_cal_int <- geese(yhat2 ~ offset(cumhaz2), data = sa, scale.fix = TRUE, 
                    family = gaussian, mean.link = "cloglog",
                    corstr = "independence", jack = TRUE)
# model for calibration slope
sa_cal_slope <- geese(yhat2 ~ offset(cumhaz2) + cumhaz2, data = sa, scale.fix = TRUE, 
                      family = gaussian, mean.link = "cloglog",
                      corstr = "independence", jack = TRUE)

# WHITE

# log cumulative hazard
w$cumhaz2 <- log(-log(1-w$kfre2))
# Pseudo-observations
km.w <- survfit(Surv(time, rrt) ~ 1, data=w)
w$yhat2 <- 1 - pseudo(km.w, times=5, type="surv")

# model for calibration intercept
w_cal_int <- geese(yhat2 ~ offset(cumhaz2), data = w, scale.fix = TRUE, 
                   family = gaussian, mean.link = "cloglog",
                   corstr = "independence", jack = TRUE)
# model for calibration slope
w_cal_slope <- geese(yhat2 ~ offset(cumhaz2) + cumhaz2, data = w, scale.fix = TRUE, 
                     family = gaussian, mean.link = "cloglog",
                     corstr = "independence", jack = TRUE)

# Calibration intercept and slope
res_cal <- rbind(  # Value and confidence interval for calibration intercept
  "SA Intercept" = with(
    summary(sa_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  # Value and confidence interval for calibration slope
  "SA Slope" = with(
    summary(sa_cal_slope)$mean["cumhaz2", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  ),
  "White Intercept" = with(
    summary(w_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  "White Slope" = with(
    summary(w_cal_slope)$mean["cumhaz2", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  )
)
res_cal <- round(res_cal, 5)
kable(res_cal) |>
  kable_styling("striped", position = "center")

rm(res_cal, sa_cal_int, sa_cal_slope, w_cal_int, w_cal_slope, km.w, km.sa)


#--- Model 3 ----
# SA
# log cumulative hazard
sa$cumhaz3 <- log(-log(1-sa$kfre3))
# Pseudo-observations
km.sa <- survfit(Surv(time, rrt) ~ 1, data=sa)
sa$yhat3 <- 1 - pseudo(km.sa, times=5, type="surv")

# model for calibration intercept
sa_cal_int <- geese(yhat3 ~ offset(cumhaz3), data = sa, scale.fix = TRUE, 
                    family = gaussian, mean.link = "cloglog",
                    corstr = "independence", jack = TRUE)
# model for calibration slope
sa_cal_slope <- geese(yhat3 ~ offset(cumhaz3) + cumhaz3, data = sa, scale.fix = TRUE, 
                      family = gaussian, mean.link = "cloglog",
                      corstr = "independence", jack = TRUE)

# WHITE

# log cumulative hazard
w$cumhaz3 <- log(-log(1-w$kfre3))
# Pseudo-observations
km.w <- survfit(Surv(time, rrt) ~ 1, data=w)
w$yhat3 <- 1 - pseudo(km.w, times=5, type="surv")

# model for calibration intercept
w_cal_int <- geese(yhat3 ~ offset(cumhaz3), data = w, scale.fix = TRUE, 
                   family = gaussian, mean.link = "cloglog",
                   corstr = "independence", jack = TRUE)
# model for calibration slope
w_cal_slope <- geese(yhat3 ~ offset(cumhaz3) + cumhaz3, data = w, scale.fix = TRUE, 
                     family = gaussian, mean.link = "cloglog",
                     corstr = "independence", jack = TRUE)

# Calibration intercept and slope
res_cal <- rbind(  # Value and confidence interval for calibration intercept
  "SA Intercept" = with(
    summary(sa_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  # Value and confidence interval for calibration slope
  "SA Slope" = with(
    summary(sa_cal_slope)$mean["cumhaz3", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  ),
  "White Intercept" = with(
    summary(w_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  "White Slope" = with(
    summary(w_cal_slope)$mean["cumhaz3", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  )
)
res_cal <- round(res_cal, 5)
kable(res_cal) |>
  kable_styling("striped", position = "center")

rm(res_cal, sa_cal_int, sa_cal_slope, w_cal_int, w_cal_slope, km.w, km.sa)
#--- Model 4 ----
# SA
# log cumulative hazard
sa$cumhaz4 <- log(-log(1-sa$kfre4))
# Pseudo-observations
km.sa <- survfit(Surv(time, rrt) ~ 1, data=sa)
sa$yhat4 <- 1 - pseudo(km.sa, times=5, type="surv")

# model for calibration intercept
sa_cal_int <- geese(yhat4 ~ offset(cumhaz4), data = sa, scale.fix = TRUE, 
                    family = gaussian, mean.link = "cloglog",
                    corstr = "independence", jack = TRUE)
# model for calibration slope
sa_cal_slope <- geese(yhat4 ~ offset(cumhaz4) + cumhaz4, data = sa, scale.fix = TRUE, 
                      family = gaussian, mean.link = "cloglog",
                      corstr = "independence", jack = TRUE)

# WHITE

# log cumulative hazard
w$cumhaz4 <- log(-log(1-w$kfre4))
# Pseudo-observations
km.w <- survfit(Surv(time, rrt) ~ 1, data=w)
w$yhat4 <- 1 - pseudo(km.w, times=5, type="surv")

# model for calibration intercept
w_cal_int <- geese(yhat4 ~ offset(cumhaz4), data = w, scale.fix = TRUE, 
                   family = gaussian, mean.link = "cloglog",
                   corstr = "independence", jack = TRUE)
# model for calibration slope
w_cal_slope <- geese(yhat4 ~ offset(cumhaz4) + cumhaz4, data = w, scale.fix = TRUE, 
                     family = gaussian, mean.link = "cloglog",
                     corstr = "independence", jack = TRUE)

# Calibration intercept and slope
res_cal <- rbind(  # Value and confidence interval for calibration intercept
  "SA Intercept" = with(
    summary(sa_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  # Value and confidence interval for calibration slope
  "SA Slope" = with(
    summary(sa_cal_slope)$mean["cumhaz4", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  ),
  "White Intercept" = with(
    summary(w_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  "White Slope" = with(
    summary(w_cal_slope)$mean["cumhaz4", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  )
)
res_cal <- round(res_cal, 5)
kable(res_cal) |>
  kable_styling("striped", position = "center")

rm(res_cal, sa_cal_int, sa_cal_slope, w_cal_int, w_cal_slope, km.w, km.sa)

# OVERALL
# log cumulative hazard
k$cumhaz4 <- log(-log(1-k$kfre4))
# Pseudo-observations
km.k <- survfit(Surv(time, rrt) ~ 1, data=k)
k$yhat4 <- 1 - pseudo(km.k, times=5, type="surv")

# model for calibration intercept
k_cal_int <- geese(yhat4 ~ offset(cumhaz4), data = k, scale.fix = TRUE, 
                   family = gaussian, mean.link = "cloglog",
                   corstr = "independence", jack = TRUE)
# model for calibration slope
k_cal_slope <- geese(yhat4 ~ offset(cumhaz4) + cumhaz4, data = k, scale.fix = TRUE, 
                     family = gaussian, mean.link = "cloglog",
                     corstr = "independence", jack = TRUE)

# Calibration intercept and slope
res_cal <- rbind(  # Value and confidence interval for calibration intercept
  "Intercept" = with(
    summary(k_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  # Value and confidence interval for calibration slope
  "Slope" = with(
    summary(k_cal_slope)$mean["cumhaz4", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  )
)
res_cal <- round(res_cal, 3)
kable(res_cal) |>
  kable_styling("striped", position = "center")

rm(res_cal, k_cal_int, k_cal_slope, km.k)
#--- Model 5 ----
# SA
# log cumulative hazard
sa$cumhaz5 <- log(-log(1-sa$kfre5))
# Pseudo-observations
km.sa <- survfit(Surv(time, rrt) ~ 1, data=sa)
sa$yhat5 <- 1 - pseudo(km.sa, times=5, type="surv")

# model for calibration intercept
sa_cal_int <- geese(yhat5 ~ offset(cumhaz5), data = sa, scale.fix = TRUE, 
                    family = gaussian, mean.link = "cloglog",
                    corstr = "independence", jack = TRUE)
# model for calibration slope
sa_cal_slope <- geese(yhat5 ~ offset(cumhaz5) + cumhaz5, data = sa, scale.fix = TRUE, 
                      family = gaussian, mean.link = "cloglog",
                      corstr = "independence", jack = TRUE)

# WHITE

# log cumulative hazard
w$cumhaz5 <- log(-log(1-w$kfre5))
# Pseudo-observations
km.w <- survfit(Surv(time, rrt) ~ 1, data=w)
w$yhat5 <- 1 - pseudo(km.w, times=5, type="surv")

# model for calibration intercept
w_cal_int <- geese(yhat5 ~ offset(cumhaz5), data = w, scale.fix = TRUE, 
                   family = gaussian, mean.link = "cloglog",
                   corstr = "independence", jack = TRUE)
# model for calibration slope
w_cal_slope <- geese(yhat5 ~ offset(cumhaz5) + cumhaz5, data = w, scale.fix = TRUE, 
                     family = gaussian, mean.link = "cloglog",
                     corstr = "independence", jack = TRUE)

# Calibration intercept and slope
res_cal <- rbind(  # Value and confidence interval for calibration intercept
  "SA Intercept" = with(
    summary(sa_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  # Value and confidence interval for calibration slope
  "SA Slope" = with(
    summary(sa_cal_slope)$mean["cumhaz5", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  ),
  "White Intercept" = with(
    summary(w_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  "White Slope" = with(
    summary(w_cal_slope)$mean["cumhaz5", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  )
)
res_cal <- round(res_cal, 5)
kable(res_cal) |>
  kable_styling("striped", position = "center")

rm(res_cal, sa_cal_int, sa_cal_slope, w_cal_int, w_cal_slope, km.w, km.sa)

# OVERALL
# log cumulative hazard
k$cumhaz5 <- log(-log(1-k$kfre5))
# Pseudo-observations
km.k <- survfit(Surv(time, rrt) ~ 1, data=k)
k$yhat5 <- 1 - pseudo(km.k, times=5, type="surv")

# model for calibration intercept
k_cal_int <- geese(yhat5 ~ offset(cumhaz5), data = k, scale.fix = TRUE, 
                   family = gaussian, mean.link = "cloglog",
                   corstr = "independence", jack = TRUE)
# model for calibration slope
k_cal_slope <- geese(yhat5 ~ offset(cumhaz5) + cumhaz5, data = k, scale.fix = TRUE, 
                     family = gaussian, mean.link = "cloglog",
                     corstr = "independence", jack = TRUE)

# Calibration intercept and slope
res_cal <- rbind(  # Value and confidence interval for calibration intercept
  "Intercept" = with(
    summary(k_cal_int)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  # Value and confidence interval for calibration slope
  "Slope" = with(
    summary(k_cal_slope)$mean["cumhaz5", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  )
)
res_cal <- round(res_cal, 3)
kable(res_cal) |>
  kable_styling("striped", position = "center")

rm(res_cal, k_cal_int, k_cal_slope, km.k)
#-------------------------------------------------------------------------------
# Discrimination - c-index
#-------------------------------------------------------------------------------

# Truncate at 5 years
w.trunc <- w
# If time>5, change outcome to 0
w.trunc$rrt <- ifelse(w.trunc$time>5, 0, w.trunc$rrt)
# For time>5, change to 5.001
w.trunc$time <- ifelse(w.trunc$time>5, 5.001, w.trunc$time)

# Truncate at 5 years
sa.trunc <- sa
# If time>5, change outcome to 0
sa.trunc$rrt <- ifelse(sa.trunc$time>5, 0, sa.trunc$rrt)
# For time>5, change to 5.001
sa.trunc$time <- ifelse(sa.trunc$time>5, 5.001, sa.trunc$time)

# Truncate at 5 years
k.trunc <- k
# If time>5, change outcome to 0
k.trunc$rrt <- ifelse(k.trunc$time>5, 0, k.trunc$rrt)
# For time>5, change to 5.001
k.trunc$time <- ifelse(k.trunc$time>5, 5.001, k.trunc$time)

#--- Model 2 ----

# Create truncated dataset (5 years cut off)
sa.trunc <- sa
w.trunc <- w
sa.trunc$rrt <- ifelse(sa.trunc$time>5, 0, sa.trunc$rrt)
w.trunc$rrt <- ifelse(w.trunc$time>5, 0, w.trunc$rrt)
sa.trunc$time <- ifelse(sa.trunc$time>5, 5.001, sa.trunc$time)
w.trunc$time <- ifelse(w.trunc$time>5, 5.001, w.trunc$time)

# c index SA
rcorr.cens(-1*sa.trunc$beta.sum, Surv(sa.trunc$time, sa.trunc$rrt))[1] 

# Bootstrap c-index to get confidence interval:
C_boot <- function(data, i) {
  val <- data[i,] # select sample with boot
  rcorr.cens(-1*val$beta.sum, Surv(val$time, val$rrt))[1]
}
set.seed(1234)
cindex.boot <- boot(sa.trunc, C_boot, R=1000)
boot.ci(boot.out=cindex.boot, type="perc")


# c index WHITE
rcorr.cens(-1*w.trunc$beta.sum, Surv(w.trunc$time, w.trunc$rrt))[1] 

# Confidence interval
set.seed(1234)
cindex.boot1 <- boot(w.trunc, C_boot, R=1000)
boot.ci(boot.out=cindex.boot1, type="perc")

rm(C_boot, cindex.boot, cindex.boot1)

#--- Model 3 ----

# c index SA
rcorr.cens(-1*sa.trunc$beta.sum3, Surv(sa.trunc$time, sa.trunc$rrt))[1] 

# Bootstrap c-index to get confidence interval:
C_boot <- function(data, i) {
  val <- data[i,] # select sample with boot
  rcorr.cens(-1*val$beta.sum3, Surv(val$time, val$rrt))[1]
}
set.seed(1234)
cindex.boot <- boot(sa.trunc, C_boot, R=1000)
boot.ci(boot.out=cindex.boot, type="perc")


# c index WHITE
rcorr.cens(-1*w.trunc$beta.sum3, Surv(w.trunc$time, w.trunc$rrt))[1] 

# Confidence interval
set.seed(1234)
cindex.boot1 <- boot(w.trunc, C_boot, R=1000)
boot.ci(boot.out=cindex.boot1, type="perc")

rm(C_boot, cindex.boot, cindex.boot1, sa.trunc, w.trunc)

#--- Model 4 ----

# c index SA
rcorr.cens(-1*sa.trunc$beta.sum4, Surv(sa.trunc$time, sa.trunc$rrt))[1] 

# Bootstrap c-index to get confidence interval:
C_boot <- function(data, i) {
  val <- data[i,] # select sample with boot
  rcorr.cens(-1*val$beta.sum4, Surv(val$time, val$rrt))[1]
}
set.seed(1234)
cindex.boot <- boot(sa.trunc, C_boot, R=1000)
boot.ci(boot.out=cindex.boot, type="perc")


# c index WHITE
rcorr.cens(-1*w.trunc$beta.sum4, Surv(w.trunc$time, w.trunc$rrt))[1] 

# Confidence interval
set.seed(1234)
cindex.boot1 <- boot(w.trunc, C_boot, R=1000)
boot.ci(boot.out=cindex.boot1, type="perc")

rm(C_boot, cindex.boot, cindex.boot1, sa.trunc, w.trunc)

# OVERALL
# Create truncated dataset (5 years cut off)
k.trunc <- k
k.trunc$rrt <- ifelse(k.trunc$time>5, 0, k.trunc$rrt)
k.trunc$time <- ifelse(k.trunc$time>5, 5.001, k.trunc$time)
# c index 
rcorr.cens(-1*k.trunc$beta.sum4, Surv(k.trunc$time, k.trunc$rrt))[1] 

# Bootstrap c-index to get confidence interval:
C_boot <- function(data, i) {
  val <- data[i,] # select sample with boot
  rcorr.cens(-1*val$beta.sum4, Surv(val$time, val$rrt))[1]
}
set.seed(1234)
cindex.boot <- boot(k.trunc, C_boot, R=1000)
boot.ci(boot.out=cindex.boot, type="perc")

rm(C_boot, cindex.boot)

#--- Model 5 ----

# c index SA
rcorr.cens(-1*sa.trunc$pi5, Surv(sa.trunc$time, sa.trunc$rrt))[1] 

# Bootstrap c-index to get confidence interval:
C_boot <- function(data, i) {
  val <- data[i,] # select sample with boot
  rcorr.cens(-1*val$pi5, Surv(val$time, val$rrt))[1]
}
set.seed(1234)
cindex.boot <- boot(sa.trunc, C_boot, R=1000)
boot.ci(boot.out=cindex.boot, type="perc")


# c index WHITE
rcorr.cens(-1*w.trunc$pi5, Surv(w.trunc$time, w.trunc$rrt))[1] 

# Confidence interval
set.seed(1234)
cindex.boot1 <- boot(w.trunc, C_boot, R=1000)
boot.ci(boot.out=cindex.boot1, type="perc")

rm(C_boot, cindex.boot, cindex.boot1, sa.trunc, w.trunc)

# OVERALL
# c index 
rcorr.cens(-1*k.trunc$pi5, Surv(k.trunc$time, k.trunc$rrt))[1] 

# Bootstrap c-index to get confidence interval:
C_boot1 <- function(data, i) {
  val <- data[i,] # select sample with boot
  rcorr.cens(-1*val$pi5, Surv(val$time, val$rrt))[1]
}
set.seed(1234)
cindex.boot <- boot(k.trunc, C_boot1, R=1000)
boot.ci(boot.out=cindex.boot, type="perc")

rm(C_boot, C_boot1, cindex.boot)

#-------------------------------------------------------------------------------
# Overall fit - Brier score
#-------------------------------------------------------------------------------
#--- Model 2 ----

brier.sa <- Score(list(sa$kfre2), formula = Surv(time, rrt) ~ 1, data=sa, metrics="Brier",
                  conf.int=TRUE, times=5, summary="ipa")
brier.w <- Score(list(w$kfre2), formula = Surv(time, rrt) ~ 1, data=w, metrics="Brier",
                 conf.int=TRUE, times=5, summary="ipa")

### Scaled Brier Score - bootstrap for 95% CI
B_boot <- function(data, i) {
  val <- data[i,]
  brier <- Score(list(val$kfre2), formula = Surv(time, rrt) ~ 1, data=val, metrics="Brier",
                 cens.model="km", conf.int=TRUE, times=5, summary="ipa")
  scaled <- brier$Brier$score$IPA
}
bboot.sa <- boot(sa, B_boot, R=100)
boot.ci(boot.out=bboot.sa, index=2, type="perc")

bboot.w <- boot(w, B_boot, R=100)
boot.ci(boot.out=bboot.w, index=2, type="perc")

#--- Model 3 ----

brier.sa <- Score(list(sa$kfre3), formula = Surv(time, rrt) ~ 1, data=sa, metrics="Brier",
                  conf.int=TRUE, times=5, summary="ipa")
brier.w <- Score(list(w$kfre3), formula = Surv(time, rrt) ~ 1, data=w, metrics="Brier",
                 conf.int=TRUE, times=5, summary="ipa")

### Scaled Brier Score - bootstrap for 95% CI
B_boot1 <- function(data, i) {
  val <- data[i,]
  brier <- Score(list(val$kfre3), formula = Surv(time, rrt) ~ 1, data=val, metrics="Brier",
                 cens.model="km", conf.int=TRUE, times=5, summary="ipa")
  scaled <- brier$Brier$score$IPA
}
bboot.sa <- boot(sa, B_boot1, R=100)
boot.ci(boot.out=bboot.sa, index=2, type="perc")

bboot.w <- boot(w, B_boot1, R=100)
boot.ci(boot.out=bboot.w, index=2, type="perc")

#--- Model 4 ----

brier.sa <- Score(list(sa$kfre4), formula = Surv(time, rrt) ~ 1, data=sa, metrics="Brier",
                  conf.int=TRUE, times=5, summary="ipa")
brier.w <- Score(list(w$kfre4), formula = Surv(time, rrt) ~ 1, data=w, metrics="Brier",
                 conf.int=TRUE, times=5, summary="ipa")

### Scaled Brier Score - bootstrap for 95% CI
B_boot1 <- function(data, i) {
  val <- data[i,]
  brier <- Score(list(val$kfre4), formula = Surv(time, rrt) ~ 1, data=val, metrics="Brier",
                 cens.model="km", conf.int=TRUE, times=5, summary="ipa")
  scaled <- brier$Brier$score$IPA
}
bboot.sa <- boot(sa, B_boot1, R=100)
boot.ci(boot.out=bboot.sa, index=2, type="perc")

bboot.w <- boot(w, B_boot1, R=100)
boot.ci(boot.out=bboot.w, index=2, type="perc")

# OVERALL
brier.k <- Score(list(k$kfre4), formula = Surv(time, rrt) ~ 1, data=k, metrics="Brier",
                 conf.int=TRUE, times=5, summary="ipa")
bboot.k <- boot(k, B_boot, R=100)
boot.ci(boot.out=bboot.k, index=2, type="perc")

#--- Model 5 ----

brier.sa <- Score(list(sa$kfre5), formula = Surv(time, rrt) ~ 1, data=sa, metrics="Brier",
                  conf.int=TRUE, times=5, summary="ipa")
brier.w <- Score(list(w$kfre5), formula = Surv(time, rrt) ~ 1, data=w, metrics="Brier",
                 conf.int=TRUE, times=5, summary="ipa")

### Scaled Brier Score - bootstrap for 95% CI
B_boot1 <- function(data, i) {
  val <- data[i,]
  brier <- Score(list(val$kfre5), formula = Surv(time, rrt) ~ 1, data=val, metrics="Brier",
                 cens.model="km", conf.int=TRUE, times=5, summary="ipa")
  scaled <- brier$Brier$score$IPA
}
bboot.sa <- boot(sa, B_boot1, R=100)
boot.ci(boot.out=bboot.sa, index=2, type="perc")

bboot.w <- boot(w, B_boot1, R=100)
boot.ci(boot.out=bboot.w, index=2, type="perc")

# OVERALL
brier.k <- Score(list(k$kfre5), formula = Surv(time, rrt) ~ 1, data=k, metrics="Brier",
                 conf.int=TRUE, times=5, summary="ipa")
bboot.k <- boot(k, B_boot, R=100)
boot.ci(boot.out=bboot.k, index=2, type="perc")

rm(brier.sa, brier.w, brier.k, bboot.sa, bboot.w, bboot.k)

#-------------------------------------------------------------------------------
# Competing Risks
#-------------------------------------------------------------------------------

# packages
library(riskRegression)
library(grid)
library(gridExtra)
library(geepack)
library(kableExtra)
library(boot)
library(Hmisc)
library(pseudo)
library(pec)
library(dplyr)
library(survival)
library(prodlim)

# Create new dataset for competing risks analysis
compk <- k
# Add in interactions as variables
compk$epi.c.white_sa <- compk$epi.c*compk$white_sa
compk$acr.c.white_sa <- compk$acr.c*compk$white_sa

# Change event indicator to three variables - ESKD=1, Death=2
compk$event <- ifelse(compk$death==1, 2, 0)
compk$event <- ifelse(compk$rrt==1, 1, compk$event)

# Risk group of pre-proposed risks
compk$riskgrp <- cut(compk$kfre5, breaks=c(0, 0.03, 0.05, 0.15, 0.25, 0.5, 1), 
                     labels=c("<3%", "3-<5%", "5-<15%", "15-<25%", "25-<50%", ">=50%"),
                     right = FALSE)

compsa <- subset(compk, compk$white_sa==1)
compw <- subset(compk, compk$white_sa==0)

#-------------------------------------------------------------------------------
# Cumulative Incidence Curves
#-------------------------------------------------------------------------------

# AJ for SA & White no risk group
aj1 <- prodlim(Hist(time, event) ~ 1, data = compsa)
aj2 <- prodlim(Hist(time, event) ~ 1, data = compw)

# Plot risk of both causes
pdf(file="cuminc.pdf", width=15, height=6)
par(mfrow=c(1,2))

plot(aj1, cause="stacked", confint=FALSE, percent=FALSE, atrisk=FALSE,
     ylab="Observed risk", xlab="Time (years)", xlim = c(0, 5), 
     legend.title="Competing risk",
     legend.legend=c("ESKD", "Death"), 
     legend.x="topleft",legend.cex=0.9, lwd=2)
#plot(km.3, type="risk", add=TRUE, confint=FALSE, xlim = c(0,5), lty="dashed", lwd=2, col="red")
#plot(km.9, type="risk", add=TRUE, confint=FALSE, xlim = c(0,5), lty="dashed", lwd=2, col="blue")
title("Cumulative Incidence curves:\nSouth Asian cohort")


plot(aj2, cause="stacked", confint=FALSE, percent=FALSE, atrisk=FALSE,
     ylab="Observed risk", xlab="Time (years)", xlim = c(0, 5), 
     legend.title="Competing risk",
     legend.legend=c("ESKD", "Death"), 
     legend.x="topleft",legend.cex=0.9, lwd=2)
#plot(km.4, type="risk", add=TRUE, confint=FALSE, xlim = c(0,5), lty="dashed", lwd=2, col="red")
#plot(km.10, type="risk", add=TRUE, confint=FALSE, xlim = c(0,5), lty="dashed", lwd=2, col="blue")
title("Cumulative Incidence curves:\nWhite cohort")
dev.off()

rm(aj1, aj2)

# AJ for k with risk groups - find cifs for KRT and death for each group
aj3 <- prodlim(Hist(time, event) ~ riskgrp, data = compk)
summary(aj3, times = 5, extend = T)

rm(aj3)


#-------------------------------------------------------------------------------
# Kaplan Meier vs Aalen-Johansen estimates (by risk group)
#-------------------------------------------------------------------------------

## AJ estimator: SA & White
aj1 <- prodlim(Hist(time, event) ~ riskgrp, data = compsa)
aj2 <- prodlim(Hist(time, event) ~ riskgrp, data = compw)
# KM estimator: SA & White
km1 <- prodlim(Hist(time, event==1) ~ riskgrp, data=compsa)
km2 <- prodlim(Hist(time, event==1) ~ riskgrp, data=compw)


pdf(file="SFig 2.pdf", width=15, height=6)
par(mfrow=c(1,2))

# White
plot(aj2, cause=1, confint=FALSE, percent=FALSE, atrisk=FALSE,
     ylab="Observed risk of ESKD", xlab="Time (years)", xlim = c(0,5),
     legend.title="Risk group (%)",
     legend.legend=c("<3", "3-<5", "5-<15", "15-<25", "25-<50", "\u226550"), 
     legend.x="topleft",legend.cex=0.85, lwd=2)
plot(km2, type="risk", add=TRUE, confint=FALSE, xlim = c(0,5), lty="dashed", lwd=2)
title("Kaplan-Meier curve & Cumulative Incidence curve:\nWhite cohort")

# South Asian
plot(aj1, cause=1, confint=FALSE, percent=FALSE, atrisk=FALSE,
     ylab="Observed risk of ESKD", xlab="Time (years)", xlim = c(0,5),
     legend.title="Risk group (%)",
     legend.legend=c("<3", "3-<5", "5-<15", "15-<25", "25-<50", "\u226550"), 
     legend.x="topleft",legend.cex=0.85, lwd=2)
plot(km1, type="risk", add=TRUE, confint=FALSE, xlim = c(0,5), lty="dashed", lwd=2)
title("Kaplan-Meier curve & Cumulative Incidence curve:\nSouth Asian cohort")

dev.off()

rm(aj1, aj2, km1, km2)

#-------------------------------------------------------------------------------
# Risks at 5 years for KM and AJ & relative risk difference
#-------------------------------------------------------------------------------

## AJ estimator: SA, White
# by risk group
aj1 <- prodlim(Hist(time, event) ~ riskgrp, data = compsa)
aj2 <- prodlim(Hist(time, event) ~ riskgrp, data = compw)
# overall
aj3 <- prodlim(Hist(time, event) ~ 1, data=compsa)
aj4 <- prodlim(Hist(time, event) ~ 1, data=compw)

# KM estimator: SA, White
# by risk group
km1 <- prodlim(Hist(time, event==1) ~ riskgrp, data=compsa)
km2 <- prodlim(Hist(time, event==1) ~ riskgrp, data=compw)
# overall
km3 <- prodlim(Hist(time, event==1) ~ 1, data=compsa)
km4 <- prodlim(Hist(time, event==1) ~ 1, data=compw)

### South Asian

# by risk group
ar.km1 <- summary(km1, extend = TRUE, times = 5, surv=FALSE) # KM
ar.aj1 <- summary(aj1, extend = TRUE, times = 5, surv=FALSE, cause=1) # AJ

# overall
ar.km3 <- summary(km3, extend = TRUE, times = 5, surv=FALSE) # KM 4.60
ar.aj3 <- summary(aj3, extend = TRUE, times = 5, surv=FALSE, cause=1) # AJ 4.41

### RRR
rrr.sa <- rep(0, 6)
for(i in 1:6) {
  rrr.sa[i] <- (ar.km1$table[[i]][[1,5]] - ar.aj1$table$`1`[[i]][[1,5]]) / ar.km1$table[[i]][[1,5]]
}
rrr.sa
# overall RRR
(ar.km3$table[[1,5]] - ar.aj3$table$`1`[[1,5]]) / ar.km3$table[[1,5]]


### White

# by risk group
ar.km2 <- summary(km2, extend = TRUE, times = 5, surv=FALSE) # KM
ar.aj2 <- summary(aj2, extend = TRUE, times = 5, surv=FALSE, cause=1) # AJ

# overall
ar.km4 <- summary(km4, extend = TRUE, times = 5, surv=FALSE) # KM
ar.aj4 <- summary(aj4, extend = TRUE, times = 5, surv=FALSE, cause=1) # AJ

### RRR
rrr.w <- rep(0, 6)
for(i in 1:6) {
  rrr.w[i] <- (ar.km2$table[[i]][[1,5]] - ar.aj2$table$`1`[[i]][[1,5]]) / ar.km2$table[[i]][[1,5]]
}
rrr.w
# overall RRR
(ar.km4$table[[1,5]] - ar.aj4$table$`1`[[1,5]]) / ar.km4$table[[1,5]]

prop.table(table(k$rrt))

rm(ar.km1, ar.km2, ar.km3, ar.km4, ar.aj1, ar.aj2, ar.aj3, ar.aj4, aj1, aj2, aj3,
   aj4, km1, km2, km3, km4)

#-------------------------------------------------------------------------------
# Fit Fine and Gray model
#-------------------------------------------------------------------------------

# see file 'Model_6_InternalValidation.R' for model fitting and internal validation

### Fit fgr model
# create matrix with covariate values for each patient
cov <- model.matrix(~ age.c + male.c + epi.c + acr.c + white_sa + white_sa:epi.c + white_sa:acr.c, 
                    data=compk)[,-1]
kint <- cbind(data.frame(cov), event = compk$event, time=compk$time)
# fit model
fgr <- FGR(as.formula(paste("Hist(time,event)~",paste(names(kint[1:7]),collapse = "+"))), 
           data=kint, cause=1, maxiter=100)
# Shrinkage factor = 0.948

compk$pi6 <- 0.948*(fgr$crrFit$coef[1]*cov[,1] + fgr$crrFit$coef[2]*cov[,2] +
                      fgr$crrFit$coef[3]*cov[,3] + fgr$crrFit$coef[4]*cov[,4] +
                      fgr$crrFit$coef[5]*cov[,5] + fgr$crrFit$coef[6]*cov[,6] +
                      fgr$crrFit$coef[7]*cov[,7])

summary(compk$pi6)

## Re-estimate baseline risk

## Baseline cif at 5 years (using STATA) == 0.035009

# predictions of int.val model
compk$pred.new <- 1 - (1-0.035009) ^ exp(compk$pi6)
summary(compk$pred.new)

rm(cov, fgr)

#-------------------------------------------------------------------------------
# Compare models
# Scatter plot of model predictions 
#-------------------------------------------------------------------------------

scatter <- ggplot(compk, aes(kfre5, pred.new)) +
  geom_point(alpha=0.7) +
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "Scatter plot of predicted risk", y = "Model 6", x = "Model 5") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1))
ggsave("SFig 4.pdf", scatter, width=15, height=12, units="cm")

rm(scatter)

#-------------------------------------------------------------------------------
# Compare model performance measures
# Average vs predicted risk
#-------------------------------------------------------------------------------

### Predicted risk

# Cox 
summary(compk$kfre5*100) # 1.61%

# FG
summary(compk$pred.new*100) # 1.55%

### Observed risk

# KM
km1 <- prodlim(Surv(time, event==1) ~ 1, data=compk)
summary(km1, extend = TRUE, times = 5, surv=FALSE) # 1.77%

# AJ
aj1 <- prodlim(Hist(time, event) ~ 1, data = compk)
summary(aj1, extend = TRUE, times = 5, surv=FALSE, cause=1) # 1.56%

rm(km1, aj1)

#-------------------------------------------------------------------------------
# Calibration - plots
#-------------------------------------------------------------------------------
#--- Model 5 -------------

# Split into deciles according to risk group - model 5
k$risk5 <- as.factor(ntile(k$kfre5, 10))
table(k$risk5)

# Predicted risks
pred.cox <- k %>%
  group_by(risk5) %>%
  summarise_at(vars(kfre5), list(mean = mean))

# AJ observed risks and conf int
cal.aj <- summary(prodlim(Hist(time, event) ~ risk5, data = compk),extend=TRUE, 
                  times = 5, cause = 1, surv=TRUE)
ajrisk <- c(cal.aj$cuminc[1], cal.aj$cuminc[2], cal.aj$cuminc[3], 
            cal.aj$cuminc[4], cal.aj$cuminc[5], cal.aj$cuminc[6], 
            cal.aj$cuminc[7], cal.aj$cuminc[8], cal.aj$cuminc[9],
            cal.aj$cuminc[10])
ajrisk.l <- c(cal.aj$lower[1], cal.aj$lower[2], cal.aj$lower[3], 
              cal.aj$lower[4], cal.aj$lower[5], cal.aj$lower[6], 
              cal.aj$lower[7], cal.aj$lower[8], cal.aj$lower[9],
              cal.aj$lower[10])
ajrisk.u <- c(cal.aj$upper[1], cal.aj$upper[2], cal.aj$upper[3], 
              cal.aj$upper[4], cal.aj$upper[5], cal.aj$upper[6], 
              cal.aj$upper[7], cal.aj$upper[8], cal.aj$upper[9],
              cal.aj$upper[10])

# plot predicted versus observed risk 
calib.coxcomp <- data.frame(pred.cox$mean, ajrisk, ajrisk.l, ajrisk.u)
cox1 <- ggplot(calib.coxcomp, aes(pred.cox.mean, ajrisk)) + 
  geom_errorbar(aes(ymin = ajrisk.l, ymax = ajrisk.u, width = 0.01), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "Model 5", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.2)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.2)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

# The distribution plot        
cox2 <- ggplot(k, aes(x = kfre5)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.2)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# Zoomed in plot
cox3 <- ggplot(calib.coxcomp, aes(pred.cox.mean, ajrisk)) + 
  geom_errorbar(aes(ymin = ajrisk.l, ymax = ajrisk.u, width = 0.0005), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "Model 5: rescaled axis", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(breaks = c(0, 0.01, 0.02), labels = scales::label_percent(accuracy=1), 
                     limits=c(0, 0.02)) +
  scale_x_continuous(breaks = c(0, 0.01, 0.02), labels = scales::label_percent(accuracy=1), 
                     limits=c(0, 0.02)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

# The distribution plot        
cox4 <- ggplot(k, aes(x = kfre5)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(breaks = c(0, 0.01, 0.02), labels = scales::label_percent(accuracy=), 
                     limits = c(0, 0.02)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))


#--- Model 6 -------------

# Split into deciles according to risk group - FG model (model 6)
compk$risk6 <- as.factor(ntile(compk$pred.new, 10))
table(compk$risk6)

# Predicted risks
pred.fg <- compk %>%
  group_by(risk6) %>%
  summarise_at(vars(pred.new), list(mean = mean))


# plot predicted versus observed risk 
calib.fg <- data.frame(pred.fg$mean, ajrisk, ajrisk.l, ajrisk.u)
fg1 <- ggplot(calib.fg, aes(pred.fg.mean, ajrisk)) + 
  geom_errorbar(aes(ymin = ajrisk.l, ymax = ajrisk.u, width = 0.01), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "Model 6", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.2)) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits=c(0, 0.2)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

# The distribution plot        
fg2 <- ggplot(compk, aes(x = pred.new)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(labels = scales::label_percent(accuracy=1), limits = c(0, 0.2)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

# Zoomed in plot
fg3 <- ggplot(calib.fg, aes(pred.fg.mean, ajrisk)) + 
  geom_errorbar(aes(ymin = ajrisk.l, ymax = ajrisk.u, width = 0.0005), colour="blue") +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, colour="coral1", linetype="dashed") +
  labs(title = "Model 6: rescaled axis", y = "Observed risk", x = "Predicted risk") +
  scale_y_continuous(breaks = c(0, 0.01, 0.02), labels = scales::label_percent(accuracy=1), 
                     limits=c(0, 0.02)) +
  scale_x_continuous(breaks = c(0, 0.01, 0.02), labels = scales::label_percent(accuracy=1), 
                     limits=c(0, 0.02)) +
  theme(axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"))

# The distribution plot        
fg4 <- ggplot(compk, aes(x = pred.new)) +
  geom_histogram(fill = "coral1", bins = 200) +
  scale_x_continuous(breaks = c(0, 0.01, 0.02), labels = scales::label_percent(accuracy=1), 
                     limits = c(0, 0.02)) +
  xlab("Predicted risk") +
  ylab("") +
  theme_minimal() +
  scale_y_continuous() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"))

sfig3 <- arrangeGrob(cox1, fg1, cox2, fg2, cox3, fg3, cox4, fg4, respect=TRUE, heights=c(1,0.3,1,0.3),
                    ncol=2, nrow=4)
grid.newpage()
grid.draw(sfig3)
ggsave("SFig 3.pdf", sfig3, width=20, height=35, units="cm")

rm(pred.cox, cal.aj, ajrisk,ajrisk.l,ajrisk.u,calib.coxcomp,cox1,cox2,pred.fg,
   calib.fg,fg1,fg2,AppI)

#-------------------------------------------------------------------------------
##### Pseudo-value plots
pdf(file="Appendix I.1.pdf", width=15, height=6)
par(mfrow=c(1,2))

# Model 5
score.m5 <- Score(
  list(compk$kfre5),
  formula = Hist(time, event) ~ 1,
  data=compk, 
  plots="cal",
  metrics=NULL,
  cens.model="km", 
  conf.int=TRUE, 
  times=5, 
  cause = 1)


score.m6 <- Score(
  list(compk$pred.new),
  formula = Hist(time, event) ~ 1,
  plots="cal",
  cens.model = "km",
  data = compk,
  conf.int = TRUE,
  times = 5,
  metrics = NULL,
  cause = 1
)

plotCalibration(score.m5, cens.method="pseudo", round=FALSE, xlim=c(0,1), 
                ylim=c(0,1), xlab="Predicted risk")
title("Model 5")

plotCalibration(score.m6, cens.method="pseudo", round=FALSE, xlim=c(0,1), 
                ylim=c(0,1), xlab="Predicted risk")
title("Model 6")
dev.off()


# Model 5 in a non-competing risks setting
score_nocomp.cox <- Score(
  list(k$kfre5),
  formula = Surv(time, rrt) ~ 1,
  cens.model = "km",
  data = k.miss,
  conf.int = TRUE,
  times = 5,
  metrics = NULL,
  plots="cal"
)

plotCalibration(score_nocomp.cox, cens.method="pseudo", round=FALSE, xlim=c(0,1), 
                ylim=c(0,1), xlab="Predicted risk")
title("Model 5 using KM risks")

#-------------------------------------------------------------------------------
# Calibration - O/E ratio
#-------------------------------------------------------------------------------
#--- Model 5 ------

aj1 <- prodlim(Hist(time, event) ~ 1, data = compk)
obj1 <- summary(aj1, extend = TRUE, times = 5, surv=FALSE, cause=1)
# O/E ratio 
OE1 <- obj1[1,6]/mean(compk$kfre5) 

#--- Model 6 ------

# O/E ratio 
OE2 <- obj1[1,6]/mean(compk$pred.new)

# Confidence interval
alpha <- 0.05
OE1_summary <- cbind(
  "Model 5" = OE1,
  "Lower .95" = exp(log(OE1 - qnorm(1 - alpha / 2) * obj1$table$`1`[1,6] / (obj1$table$`1`[1,5]))),
  "Upper .95" = exp(log(OE1 + qnorm(1 - alpha / 2) * obj1$table$`1`[1,6] / (obj1$table$`1`[1,5]))),
  "Model 6" = OE2,
  "Lower .95" = exp(log(OE2 - qnorm(1 - alpha / 2) * obj1$table$`1`[1,6] / (obj1$table$`1`[1,5]))),
  "Upper .95" = exp(log(OE2 + qnorm(1 - alpha / 2) * obj1$table$`1`[1,6] / (obj1$table$`1`[1,5])))
)
OE1_summary <- round(OE1_summary, 3)
OE1_summary


rm(aj1, obj1, OE1, OE2, alpha, OE1_summary)


#-------------------------------------------------------------------------------
# Calibration intercept and slope
#-------------------------------------------------------------------------------
#--- Model 5 -----
# log cumulative hazard
compk$cumhaz5 <- log(-log(1-compk$kfre5))

# pseudo obs
aj1 <- prodlim(Hist(time, event) ~ 1, data = compk)
compk$yhat5 <- jackknife(aj1, times=5, cause = 1)

# model for calibration intercept
m5_intercept <- geese(yhat5 ~ offset(cumhaz5), data = compk, scale.fix=TRUE,
                      family = gaussian, mean.link="cloglog",
                      corstr = "independence", jack = TRUE)

# model for calibration slope
m5_slope <- geese(yhat5 ~ offset(cumhaz5) + cumhaz5, data = compk, scale.fix = TRUE, 
                  family = gaussian, mean.link = "cloglog",
                  corstr = "independence", jack = TRUE)


#--- Model 6 -----

# log cumulative hazard
compk$cumhaz6 <- log(-log(1-compk$pred.new))
# pseudo obs - as above

# model for calibration intercept
m6_intercept <- geese(yhat5 ~ offset(cumhaz6), data = compk, scale.fix=TRUE,
                      family = gaussian, mean.link="cloglog",
                      corstr = "independence", jack = TRUE)

# model for calibration slope
m6_slope <- geese(yhat5 ~ offset(cumhaz6) + cumhaz6, data = compk, scale.fix = TRUE, 
                  family = gaussian, mean.link = "cloglog",
                  corstr = "independence", jack = TRUE)


res_cal <- rbind(
  "Model 5 intercept" = with(
    summary(m5_intercept)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  # Value, confidence interval and test for calibration slope
  "Model 5 slope" = with(
    summary(m5_slope)$mean["cumhaz5", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  ),
  "Model 6 intercept" = with(
    summary(m6_intercept)$mean,
    c(
      "estimate" = estimate,
      `2.5 %` = estimate - qnorm(0.975) * san.se,
      `97.5 %` = estimate + qnorm(0.975) * san.se
    )
  ),
  "Model 6 slope" = with(
    summary(m6_slope)$mean["cumhaz6", ],
    c(
      "estimate" = 1 + estimate,
      `2.5 %` = 1 + (estimate - qnorm(0.975) * san.se),
      `97.5 %` = 1 + (estimate + qnorm(0.975) * san.se)
    )
  )
)
res_cal <- round(res_cal, 3)
kable(res_cal) |>
  kable_styling("striped", position = "center")

rm(m5_intercept,m5_slope,m6_intercept,m6_slope,res_cal,aj1)

#-------------------------------------------------------------------------------
# Discrimination
#-------------------------------------------------------------------------------
#--- Model 5 ----

# Create truncated dataset (5 years cut off)
compk.trunc <- compk
compk.trunc$event <- ifelse(compk.trunc$time>5, 0, compk.trunc$event)
compk.trunc$time <- ifelse(compk.trunc$time>5, 5.001, compk.trunc$time)

cindex.m5 <- pec::cindex(as.matrix(compk.trunc$kfre5), formula = Hist(time, event) ~ 1, cause = 1,
                         eval.times = 5, data = compk.trunc)
cindex.m5$AppCindex$matrix

C_boot1 <- function(data, i) {
  val <- data[i,]
  pec::cindex(as.matrix(val$pi5), formula = Hist(time, event) ~ 1, cause = 1,
              eval.times = 5, data = val)$AppCindex$matrix
}
set.seed(1234)
cindex.boot.m5 <- boot(compk.trunc, C_boot1, R=1000)
boot.ci(boot.out=cindex.boot.m5, type="perc")


#--- Model 6 ----
cindex.m6 <- pec::cindex(as.matrix(compk.trunc$pi6), formula = Hist(time, event) ~ 1, cause = 1,
                         eval.times = 5, data = compk.trunc)
cindex.m6$AppCindex$matrix 

# Bootstrap c-index to get confidence interval:
C_boot2 <- function(data, i) {
  val <- data[i,]
  cindex <- pec::cindex(as.matrix(val$pi6), formula = Hist(time, event) ~ 1, cause = 1,
                        eval.times = 5, data = val)$AppCindex$matrix
}
set.seed(1234)
cindex.boot.m6 <- boot(compk.trunc, C_boot2, R=1000)
boot.ci(boot.out=cindex.boot.m6, type="perc")

rm(C_boot1, C_boot2, cindex.m5, cindex.m6, compk.trunc)

#-------------------------------------------------------------------------------
# Brier score
#-------------------------------------------------------------------------------
#--- Model 5 -----
score_m5 <- Score(
  list(compk$kfre5),
  formula = Hist(time, event) ~ 1,
  cens.model = "km",
  data = compk,
  conf.int = TRUE,
  times = 5,
  metrics = "brier",
  summary = c("ipa"),
  cause = 1
)

# Scaled Brier Score - bootstrap for 95% CI
B_boot <- function(data, i) {
  val <- data[i,]
  brier <- Score(list(val$kfre5), formula = Surv(time, rrt) ~ 1, data=val, metrics="Brier",
                 cens.model="km", conf.int=TRUE, times=5, summary="ipa")
  scaled <- brier$Brier$score$IPA
}
bboot.m5 <- boot(compk, B_boot, R=1000)
boot.ci(boot.out=bboot.m5, index=2, type="perc")

#--- Model 6 ----
score_m6 <- Score(
  list(compk$pred.new),
  formula = Hist(time, event) ~ 1,
  cens.model = "km",
  data = compk,
  conf.int = TRUE,
  times = 5,
  metrics = "brier",
  summary = c("ipa"),
  cause = 1
)

# Scaled Brier Score - bootstrap for 95% CI
B_boot <- function(data, i) {
  val <- data[i,]
  brier <- Score(list(val$pred.new), formula = Surv(time, rrt) ~ 1, data=val, metrics="Brier",
                 cens.model="km", conf.int=TRUE, times=5, summary="ipa")
  scaled <- brier$Brier$score$IPA
}
bboot.m6 <- boot(compk, B_boot, R=1000)
boot.ci(boot.out=bboot.m6, index=2, type="perc")


rm(score_m5, score_m6, B_boot, bboot.m5, bboot.m6)

#-------------------------------------------------------------------------------
# Clinical Impact Assessment - all models
#-------------------------------------------------------------------------------
#--- SET UP ----
source("stdca.R")

# Assess clinical impact in the eligibility assessment cohort - those who are not known 
# already to secondary care, and are from the Leicester cohort (i.e. psp==0,
# neph_known==0)

# Eligibility assessment cohort
compk1 <- subset(compk, compk$neph_known==0 & compk$psp==0)

# Variable for if rrt developed within 5 years:
compk1$rrt5 <- ifelse(compk1$rrt==1 & compk1$time<=5, 1, 0)


# Referral according to NICE criteria (eGFR<30, ACR>=70)
compk1$nice <- ifelse(compk1$epi<30 | compk1$acr_mgmmol>=70, 1, 0)

# indicator variable for predicted risk >5% and ACR >=70 for models 2-6
compk1$m2 <- ifelse(compk1$kfre2>=0.05 | compk1$acr_mgmmol>=70, 1, 0)
compk1$m3 <- ifelse(compk1$kfre3>=0.05 | compk1$acr_mgmmol>=70, 1, 0)
compk1$m4 <- ifelse(compk1$kfre4>=0.05 | compk1$acr_mgmmol>=70, 1, 0)
compk1$m5 <- ifelse(compk1$kfre5>=0.05 | compk1$acr_mgmmol>=70, 1, 0)
compk1$m6 <- ifelse(compk1$pred.new>=0.05 | compk1$acr_mgmmol>=70, 1, 0)


# eligibility assessment cohort - SOUTH ASIAN
compsa1 <- subset(compk1, compk1$white_sa==1)

# eligibility assessment cohort - WHITE
compw1 <- subset(compk1, compk1$white_sa==0)

#--- Decision Curve Analysis ----

### DCA PLOT - WHITE
netben.white = stdca(data=compw1, outcome="rrt", ttoutcome="time", timepoint=5, xstop=0.12,  
                     predictors=c("nice", "kfre2", "kfre3", "kfre4", "kfre5", "pred.new"), graph=FALSE)
netben.white1 <- data.frame(thresh = netben.white$net.benefit$threshold, 
                            all = netben.white$net.benefit$all, 
                            none = netben.white$net.benefit$none,
                            nice =netben.white$net.benefit$nice, 
                            kfre2 = netben.white$net.benefit$kfre2, 
                            kfre3 = netben.white$net.benefit$kfre3, 
                            kfre4 = netben.white$net.benefit$kfre4,
                            kfre5 = netben.white$net.benefit$kfre5,
                            pred.new = netben.white$net.benefit$pred.new)

colours1 <- c("All"="coral2", "None"="black", "NICE"="pink", "Model 2"="green", 
              "Model 3"="yellow", "Model 4"="orange", "Model 5"="#DC71FA", 
              "Model 6"="skyblue")

p1 <- ggplot(netben.white1, aes(x=thresh)) + 
  geom_line(aes(y=none, col="None"), lwd=1) +
  geom_line(aes(y=all, col="All"), lwd=1) +
  geom_line(aes(y=nice, col="NICE"), lwd=1) +
  geom_line(aes(y=kfre2, col="Model 2"), lwd=1) +
  geom_line(aes(y=kfre3, col="Model 3"), linetype="longdash", lwd=1) +
  geom_line(aes(y=kfre4, col="Model 4"), linetype="dashed", lwd=1) +
  geom_line(aes(y=kfre5, col="Model 5"), lwd=1) +
  geom_line(aes(y=pred.new, col="Model 6"), linetype="dashed", lwd=1) +
  labs(title="White cohort", x="Threshold probability", 
       y="Net benefit", colour="Legend") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), limits=c(-0.02, 0.02)) +
  scale_color_manual(values = colours1)

### DCA PLOT - SOUTH ASIAN
netben.sa = stdca(data=compsa1, outcome="rrt", ttoutcome="time", timepoint=5, xstop=0.12,  
                  predictors=c("nice", "kfre2", "kfre3", "kfre4", "kfre5", "pred.new"), graph=FALSE)
netben.sa1 <- data.frame(thresh = netben.sa$net.benefit$threshold, 
                         all = netben.sa$net.benefit$all, 
                         none = netben.sa$net.benefit$none,
                         nice =netben.sa$net.benefit$nice, 
                         kfre2 = netben.sa$net.benefit$kfre2, 
                         kfre3 = netben.sa$net.benefit$kfre3, 
                         kfre4 = netben.sa$net.benefit$kfre4,
                         kfre5 = netben.sa$net.benefit$kfre5,
                         pred.new = netben.sa$net.benefit$pred.new)

p2 <- ggplot(netben.sa1, aes(x=thresh)) + 
  geom_line(aes(y=none, col="None"), lwd=1) +
  geom_line(aes(y=all, col="All"), lwd=1) +
  geom_line(aes(y=nice, col="NICE"), lwd=1) +
  geom_line(aes(y=kfre2, col="Model 2"), lwd=1) +
  geom_line(aes(y=kfre3, col="Model 3"), linetype="longdash", lwd=1) +
  geom_line(aes(y=kfre4, col="Model 4"), linetype="dashed", lwd=1) +
  geom_line(aes(y=kfre5, col="Model 5"), lwd=1) +
  geom_line(aes(y=pred.new, col="Model 6"), linetype="dashed", lwd=1) +
  labs(title="South Asian cohort", x="Threshold probability", 
       y="Net benefit", colour="Legend") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), limits=c(-0.02, 0.02)) +
  scale_color_manual(values = colours1)

### DCA PLOT - OVERALL
netben.overall = stdca(data=compk1, outcome="rrt", ttoutcome="time", timepoint=5, xstop=0.12,  
                       predictors=c("nice", "kfre2", "kfre3", "kfre4", "kfre5", "pred.new"), graph=FALSE)
netben.overall1 <- data.frame(thresh = netben.overall$net.benefit$threshold, 
                              all = netben.overall$net.benefit$all, 
                              none = netben.overall$net.benefit$none,
                              nice =netben.overall$net.benefit$nice, 
                              kfre2 = netben.overall$net.benefit$kfre2, 
                              kfre3 = netben.overall$net.benefit$kfre3, 
                              kfre4 = netben.overall$net.benefit$kfre4,
                              kfre5 = netben.overall$net.benefit$kfre5,
                              pred.new = netben.overall$net.benefit$pred.new)

p3 <- ggplot(netben.overall1, aes(x=thresh)) + 
  geom_line(aes(y=none, col="None"), lwd=1) +
  geom_line(aes(y=all, col="All"), lwd=1) +
  geom_line(aes(y=nice, col="NICE"), lwd=1) +
  geom_line(aes(y=kfre2, col="Model 2"), lwd=1) +
  geom_line(aes(y=kfre3, col="Model 3"), linetype="longdash", lwd=1) +
  geom_line(aes(y=kfre4, col="Model 4"), linetype="dashed", lwd=1) +
  geom_line(aes(y=kfre5, col="Model 5"), lwd=1) +
  geom_line(aes(y=pred.new, col="Model 6"), linetype="dashed", lwd=1) +
  labs(title="Overall", x="Threshold probability", 
       y="Net benefit", colour="Legend") +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), limits=c(-0.02, 0.02)) +
  scale_color_manual(values = colours1)

# Combine graphs
dca.graph <- ggarrange(p1, p2, p3, common.legend=TRUE, legend = "bottom", ncol = 3)
ggsave("Figure 3.pdf", dca.graph, width=38, height=16, units="cm")

rm(p1,p2,p3, netben.white,netben.white1,netben.sa,netben.sa1,netben.overall,
   netben.overall1,dca.graph,colours1)

#--- Net benefit tables ----
#--- NICE ----
# White
# no. of referrals
sum(compw1$nice)
100*sum(compw1$nice)/length(compw1$id)
# correct no.
compw1$sens2 <- ifelse(compw1$nice==1 & compw1$rrt5==1, 1, 0)
sum(compw1$sens2)
100*sum(compw1$sens2)/sum(compw1$nice)
# incorrect no.
compw1$unnec2 <- ifelse(compw1$nice==1 & compw1$rrt5==0, 1, 0)
sum(compw1$unnec2)
100*sum(compw1$unnec2)/sum(compw1$nice)
# no. of missed cases
compw1$false2 <- ifelse(compw1$nice==0 & compw1$rrt5==1, 1, 0)
sum(compw1$false2)
100*sum(compw1$false2)/(sum(compw1$rrt5==1))

# SA
# no. of referrals
sum(compsa1$nice)
100*sum(compsa1$nice)/length(compsa1$id)
# correct no.
compsa1$sens2 <- ifelse(compsa1$nice==1 & compsa1$rrt5==1, 1, 0)
sum(compsa1$sens2)
100*sum(compsa1$sens2)/sum(compsa1$nice)
# incorrect no.
compsa1$unnec2 <- ifelse(compsa1$nice==1 & compsa1$rrt5==0, 1, 0)
sum(compsa1$unnec2)
100*sum(compsa1$unnec2)/sum(compsa1$nice)
# no. of missed cases
compsa1$false2 <- ifelse(compsa1$nice==0 & compsa1$rrt5==1, 1, 0)
sum(compsa1$false2)
100*sum(compsa1$false2)/(sum(compsa1$rrt5==1))

# Overall
# no. of referrals
sum(compk1$nice)
100*sum(compk1$nice)/length(compk1$id)
# correct no.
compk1$sens2 <- ifelse(compk1$nice==1 & compk1$rrt5==1, 1, 0)
sum(compk1$sens2)
100*sum(compk1$sens2)/sum(compk1$nice)
# incorrect no.
compk1$unnec2 <- ifelse(compk1$nice==1 & compk1$rrt5==0, 1, 0)
sum(compk1$unnec2)
100*sum(compk1$unnec2)/sum(compk1$nice)
# no. of missed cases
compk1$false2 <- ifelse(compk1$nice==0 & compk1$rrt5==1, 1, 0)
sum(compk1$false2)
100*sum(compk1$false2)/(sum(compk1$rrt5==1))

#--- Model 2 ----

# White
# no. of referrals
sum(compw1$m2)
100*sum(compw1$m2)/length(compw1$id)
# correct no.
compw1$sens2 <- ifelse(compw1$m2==1 & compw1$rrt5==1, 1, 0)
sum(compw1$sens2)
100*sum(compw1$sens2)/sum(compw1$m2)
# incorrect no.
compw1$unnec2 <- ifelse(compw1$m2==1 & compw1$rrt5==0, 1, 0)
sum(compw1$unnec2)
100*sum(compw1$unnec2)/sum(compw1$m2)
# no. of missed cases
compw1$false2 <- ifelse(compw1$m2==0 & compw1$rrt5==1, 1, 0)
sum(compw1$false2)
100*sum(compw1$false2)/(sum(compw1$rrt5==1))

# SA
# no. of referrals
sum(compsa1$m2)
100*sum(compsa1$m2)/length(compsa1$id)
# correct no.
compsa1$sens2 <- ifelse(compsa1$m2==1 & compsa1$rrt5==1, 1, 0)
sum(compsa1$sens2)
100*sum(compsa1$sens2)/sum(compsa1$m2)
# incorrect no.
compsa1$unnec2 <- ifelse(compsa1$m2==1 & compsa1$rrt5==0, 1, 0)
sum(compsa1$unnec2)
100*sum(compsa1$unnec2)/sum(compsa1$m2)
# no. of missed cases
compsa1$false2 <- ifelse(compsa1$m2==0 & compsa1$rrt5==1, 1, 0)
sum(compsa1$false2)
100*sum(compsa1$false2)/(sum(compsa1$rrt5==1))

# Overall
# no. of referrals
sum(compk1$m2)
100*sum(compk1$m2)/length(compk1$id)
# correct no.
compk1$sens2 <- ifelse(compk1$m2==1 & compk1$rrt5==1, 1, 0)
sum(compk1$sens2)
100*sum(compk1$sens2)/sum(compk1$m2)
# incorrect no.
compk1$unnec2 <- ifelse(compk1$m2==1 & compk1$rrt5==0, 1, 0)
sum(compk1$unnec2)
100*sum(compk1$unnec2)/sum(compk1$m2)
# no. of missed cases
compk1$false2 <- ifelse(compk1$m2==0 & compk1$rrt5==1, 1, 0)
sum(compk1$false2)
100*sum(compk1$false2)/(sum(compk1$rrt5==1))

#--- Model 3 ----
# White
# no. of referrals
sum(compw1$m3)
100*sum(compw1$m3)/length(compw1$id)
# correct no.
compw1$sens3 <- ifelse(compw1$m3==1 & compw1$rrt5==1, 1, 0)
sum(compw1$sens3)
100*sum(compw1$sens3)/sum(compw1$m3)
# incorrect no.
compw1$unnec3 <- ifelse(compw1$m3==1 & compw1$rrt5==0, 1, 0)
sum(compw1$unnec3)
100*sum(compw1$unnec3)/sum(compw1$m3)
# no. of missed cases
compw1$false3 <- ifelse(compw1$m3==0 & compw1$rrt5==1, 1, 0)
sum(compw1$false3)
100*sum(compw1$false3)/(sum(compw1$rrt5==1))

# SA
# no. of referrals
sum(compsa1$m3)
100*sum(compsa1$m3)/length(compsa1$id)
# correct no.
compsa1$sens3 <- ifelse(compsa1$m3==1 & compsa1$rrt5==1, 1, 0)
sum(compsa1$sens3)
100*sum(compsa1$sens3)/sum(compsa1$m3)
# incorrect no.
compsa1$unnec3 <- ifelse(compsa1$m3==1 & compsa1$rrt5==0, 1, 0)
sum(compsa1$unnec3)
100*sum(compsa1$unnec3)/sum(compsa1$m3)
# no. of missed cases
compsa1$false3 <- ifelse(compsa1$m3==0 & compsa1$rrt5==1, 1, 0)
sum(compsa1$false3)
100*sum(compsa1$false3)/(sum(compsa1$rrt5==1))

# Overall
# no. of referrals
sum(compk1$m3)
100*sum(compk1$m3)/length(compk1$id)
# correct no.
compk1$sens3 <- ifelse(compk1$m3==1 & compk1$rrt5==1, 1, 0)
sum(compk1$sens3)
100*sum(compk1$sens3)/sum(compk1$m3)
# incorrect no.
compk1$unnec3 <- ifelse(compk1$m3==1 & compk1$rrt5==0, 1, 0)
sum(compk1$unnec3)
100*sum(compk1$unnec3)/sum(compk1$m3)
# no. of missed cases
compk1$false3 <- ifelse(compk1$m3==0 & compk1$rrt5==1, 1, 0)
sum(compk1$false3)
100*sum(compk1$false3)/(sum(compk1$rrt5==1))
#--- Model 4 ----
# White
# no. of referrals
sum(compw1$m4)
100*sum(compw1$m4)/length(compw1$id)
# correct no.
compw1$sens4 <- ifelse(compw1$m4==1 & compw1$rrt5==1, 1, 0)
sum(compw1$sens4)
100*sum(compw1$sens4)/sum(compw1$m4)
# incorrect no.
compw1$unnec4 <- ifelse(compw1$m4==1 & compw1$rrt5==0, 1, 0)
sum(compw1$unnec4)
100*sum(compw1$unnec4)/sum(compw1$m4)
# no. of missed cases
compw1$false4 <- ifelse(compw1$m4==0 & compw1$rrt5==1, 1, 0)
sum(compw1$false4)
100*sum(compw1$false4)/(sum(compw1$rrt5==1))

# SA
# no. of referrals
sum(compsa1$m4)
100*sum(compsa1$m4)/length(compsa1$id)
# correct no.
compsa1$sens4 <- ifelse(compsa1$m4==1 & compsa1$rrt5==1, 1, 0)
sum(compsa1$sens4)
100*sum(compsa1$sens4)/sum(compsa1$m4)
# incorrect no.
compsa1$unnec4 <- ifelse(compsa1$m4==1 & compsa1$rrt5==0, 1, 0)
sum(compsa1$unnec4)
100*sum(compsa1$unnec4)/sum(compsa1$m4)
# no. of missed cases
compsa1$false4 <- ifelse(compsa1$m4==0 & compsa1$rrt5==1, 1, 0)
sum(compsa1$false4)
100*sum(compsa1$false4)/(sum(compsa1$rrt5==1))

# Overall
# no. of referrals
sum(compk1$m4)
100*sum(compk1$m4)/length(compk1$id)
# correct no.
compk1$sens4 <- ifelse(compk1$m4==1 & compk1$rrt5==1, 1, 0)
sum(compk1$sens4)
100*sum(compk1$sens4)/sum(compk1$m4)
# incorrect no.
compk1$unnec4 <- ifelse(compk1$m4==1 & compk1$rrt5==0, 1, 0)
sum(compk1$unnec4)
100*sum(compk1$unnec4)/sum(compk1$m4)
# no. of missed cases
compk1$false4 <- ifelse(compk1$m4==0 & compk1$rrt5==1, 1, 0)
sum(compk1$false4)
100*sum(compk1$false4)/(sum(compk1$rrt5==1))

#--- Model 5 ----
# White
# no. of referrals
sum(compw1$m5)
100*sum(compw1$m5)/length(compw1$id)
# correct no.
compw1$sens5 <- ifelse(compw1$m5==1 & compw1$rrt5==1, 1, 0)
sum(compw1$sens5)
100*sum(compw1$sens5)/sum(compw1$m5)
# incorrect no.
compw1$unnec5 <- ifelse(compw1$m5==1 & compw1$rrt5==0, 1, 0)
sum(compw1$unnec5)
100*sum(compw1$unnec5)/sum(compw1$m5)
# no. of missed cases
compw1$false5 <- ifelse(compw1$m5==0 & compw1$rrt5==1, 1, 0)
sum(compw1$false5)
100*sum(compw1$false5)/(sum(compw1$rrt5==1))

# SA
# no. of referrals
sum(compsa1$m5)
100*sum(compsa1$m5)/length(compsa1$id)
# correct no.
compsa1$sens5 <- ifelse(compsa1$m5==1 & compsa1$rrt5==1, 1, 0)
sum(compsa1$sens5)
100*sum(compsa1$sens5)/sum(compsa1$m5)
# incorrect no.
compsa1$unnec5 <- ifelse(compsa1$m5==1 & compsa1$rrt5==0, 1, 0)
sum(compsa1$unnec5)
100*sum(compsa1$unnec5)/sum(compsa1$m5)
# no. of missed cases
compsa1$false5 <- ifelse(compsa1$m5==0 & compsa1$rrt5==1, 1, 0)
sum(compsa1$false5)
100*sum(compsa1$false5)/(sum(compsa1$rrt5==1))

# Overall
# no. of referrals
sum(compk1$m5)
100*sum(compk1$m5)/length(compk1$id)
# correct no.
compk1$sens5 <- ifelse(compk1$m5==1 & compk1$rrt5==1, 1, 0)
sum(compk1$sens5)
100*sum(compk1$sens5)/sum(compk1$m5)
# incorrect no.
compk1$unnec5 <- ifelse(compk1$m5==1 & compk1$rrt5==0, 1, 0)
sum(compk1$unnec5)
100*sum(compk1$unnec5)/sum(compk1$m5)
# no. of missed cases
compk1$false5 <- ifelse(compk1$m5==0 & compk1$rrt5==1, 1, 0)
sum(compk1$false5)
100*sum(compk1$false5)/(sum(compk1$rrt5==1))


#--- Model 6 ----

# White
# no. of referrals
sum(compw1$m6)
100*sum(compw1$m6)/length(compw1$id)
# correct no.
compw1$sens6 <- ifelse(compw1$m6==1 & compw1$rrt5==1, 1, 0)
sum(compw1$sens6)
100*sum(compw1$sens6)/sum(compw1$m6)
# incorrect no.
compw1$unnec6 <- ifelse(compw1$m6==1 & compw1$rrt5==0, 1, 0)
sum(compw1$unnec6)
100*sum(compw1$unnec6)/sum(compw1$m6)
# no. of missed cases
compw1$false6 <- ifelse(compw1$m6==0 & compw1$rrt5==1, 1, 0)
sum(compw1$false6)
100*sum(compw1$false6)/(sum(compw1$rrt5==1))

# SA
# no. of referrals
sum(compsa1$m6)
100*sum(compsa1$m6)/length(compsa1$id)
# correct no.
compsa1$sens6 <- ifelse(compsa1$m6==1 & compsa1$rrt5==1, 1, 0)
sum(compsa1$sens6)
100*sum(compsa1$sens6)/sum(compsa1$m6)
# incorrect no.
compsa1$unnec6 <- ifelse(compsa1$m6==1 & compsa1$rrt5==0, 1, 0)
sum(compsa1$unnec6)
100*sum(compsa1$unnec6)/sum(compsa1$m6)
# no. of missed cases
compsa1$false6 <- ifelse(compsa1$m6==0 & compsa1$rrt5==1, 1, 0)
sum(compsa1$false6)
100*sum(compsa1$false6)/(sum(compsa1$rrt5==1))

# Overall
# no. of referrals
sum(compk1$m6)
100*sum(compk1$m6)/length(compk1$id)
# correct no.
compk1$sens6 <- ifelse(compk1$m6==1 & compk1$rrt5==1, 1, 0)
sum(compk1$sens6)
100*sum(compk1$sens6)/sum(compk1$m6)
# incorrect no.
compk1$unnec6 <- ifelse(compk1$m6==1 & compk1$rrt5==0, 1, 0)
sum(compk1$unnec6)
100*sum(compk1$unnec6)/sum(compk1$m6)
# no. of missed cases
compk1$false6 <- ifelse(compk1$m6==0 & compk1$rrt5==1, 1, 0)
sum(compk1$false6)
100*sum(compk1$false6)/(sum(compk1$rrt5==1))

# Percentage decrease in unnec. referrals from NICE to Model 6
((736 - 509)/736)*100 # 30.84

# Perc dec. in missed ESKD cases
((37 - 29)/37)*100 # 21.6

# Percentage increase in unnecessary referrals from Model 5 to Model 6 (overall)
((509 - 506)/506)*100 # 0.59%

# Perc dec in missed ESKD cases
((34-32)/34)*100 # 5.88

# Perc inc in correct referrals
((53-52)/52)*100 # 1.92



