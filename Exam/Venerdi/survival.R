library("survival")
library("splines")
library("lattice")
library("JM")
################################################################################
pbcseq %<>%
  mutate(id = as.integer(id),
         age = as.double(age+(day/365.25)),
         ascites = as.factor(ascites),
         hepato = as.factor(hepato),
         spiders = as.factor(spiders),
         edema = ifelse(edema==1,'serious',ifelse(edema==0.5,'evidence','no')))
################################################################################
df1 <- mice(pbcseq[,-c(8,9,10)], m=10, maxit = 50, method = 'pmm', seed = 1234)
df2 <- mice(pbcseq[,c(8,9,10)], m=10, maxit = 50, method = 'logreg', seed = 1234)

densityplot(df1)

pbcseq2 <- complete(df1,1)
pbcseq2 <- cbind(pbcseq2,complete(df2,1))
################################################################################
pbcseq.survival.head <- pbcseq2 %>%
  mutate(day = as.integer(day),
         edema = as.factor(edema),
         sex = as.factor(ifelse(sex=='m','male','female')),
         stage = as.factor(stage),
         years = futime/365.25,
         status = ifelse(status=='2','dead',ifelse(status=='1','transplanted','alive')),
         status2 = ifelse(status=='dead',1,0),
         trt = as.factor(ifelse(trt==1,'D-penicil','Placebo'))) %>%
  group_by(id) %>%
  do(head(., n = 1)) %>%
  ungroup() %>%
  as.data.frame(.)

pbcseq.survival.tail <- pbcseq2 %>%
  mutate(day = as.integer(day),
         edema = as.factor(edema),
         sex = as.factor(ifelse(sex=='m','male','female')),
         stage = as.factor(stage),
         years = futime/365.25,
         status = ifelse(status=='2','dead',ifelse(status=='1','transplanted','alive')),
         status2 = ifelse(status=='dead',1,0),
         trt = as.factor(ifelse(trt==1,'D-penicil','Placebo'))) %>%
  group_by(id) %>%
  do(tail(., n = 1)) %>%
  ungroup() %>%
  as.data.frame(.)
################################################################################
KM_fit <- survfit(Surv(years, status2) ~ sex, data = pbcseq.survival.head)
ggsurvplot(KM_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

ggsurvplot(
  KM_fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in years",   # customize X axis label.
  break.time.by = 2,     # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs =
    c("Male", "Female"),    # change legend labels.
  palette =
    c("#E7B800", "#2E9FDF"))


KM_fit <- survfit(Surv(years, status2) ~ trt, data = pbcseq.survival.head)
ggsurvplot(KM_fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

fit_weibull <- survreg(Surv(years, status2) ~ trt + sex + age, data = pbcseq.survival.head)
summary(fit_weibull)

fit_exp <- survreg(Surv(years, status2) ~ trt + sex + age, data = pbcseq.survival.head,
                   dist = "exponential")
summary(fit_exp)

fit_lnorm <- survreg(Surv(years, status2) ~ trt + sex + age, data = pbcseq.survival.head,
                     dist = "lognormal")
summary(fit_lnorm)

fit_llogis <- survreg(Surv(years, status2) ~ trt + sex + age, data = pbcseq.survival.head,
                      dist = "loglogistic")
summary(fit_llogis)

fit.quad <- survreg(Surv(years, status2) ~ (trt + sex) * (age + I(age^2)),
               data = pbcseq.survival.head)
summary(fit.quad)

anova(fit_weibull, fit.quad) # no quad

fit_weib <- survreg(Surv(years, status2) ~ trt * sex, data = pbcseq.survival.head)
fitted_values <- fit_weib$linear.predictors
resids <- (log(fit_weib$y[, 1]) - fitted_values) / fit_weib$scale
resKM <- survfit(Surv(resids, status2) ~ 1, data = pbcseq.survival.head)
par(mfrow=c(1,3))
plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability")
xx <- seq(min(resids), max(resids), length.out = 35)
yy <- exp(- exp(xx)) #good
#yy <- pnorm(xx, lower.tail = FALSE) #bad
#yy <- plogis(xx, lower.tail = FALSE) #bad
lines(xx, yy, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate",
                       "Survival function of Extreme Value distribution"),
       lty = c(1,2,1), col = c(1,1,2), bty = "n")
plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability")
xx <- seq(min(resids), max(resids), length.out = 35)
#yy <- exp(- exp(xx)) #good
yy <- pnorm(xx, lower.tail = FALSE) #bad
#yy <- plogis(xx, lower.tail = FALSE) #bad
lines(xx, yy, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate",
                       "Survival function of Extreme Value distribution"),
       lty = c(1,2,1), col = c(1,1,2), bty = "n")
plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability")
xx <- seq(min(resids), max(resids), length.out = 35)
#yy <- exp(- exp(xx)) #good
#yy <- pnorm(xx, lower.tail = FALSE) #bad
yy <- plogis(xx, lower.tail = FALSE) #bad
lines(xx, yy, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate",
                       "Survival function of Extreme Value distribution"),
       lty = c(1,2,1), col = c(1,1,2), bty = "n")

fit <- coxph(Surv(years, status2) ~ trt + sex + age, data = pbcseq.survival.head)
summary(fit) # no trt

fit_null <- coxph(Surv(years, status2) ~ sex, data = pbcseq.survival.head)
fit_alt <- coxph(Surv(years, status2) ~ trt * sex, data = pbcseq.survival.head)

anova(fit_null, fit_alt) # null

fit <- coxph(Surv(years, status2) ~ sex * age, data = pbcseq.survival.head)
check_PH <- cox.zph(fit, transform = "km")
check_PH
ggcoxzph(check_PH)
ggcoxdiagnostics(fit, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxdiagnostics(fit, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())

pbcseq.survival.head$biliGrp <- cut(pbcseq.survival.head$bili, breaks = c(min(pbcseq.survival.head$bili), .8, 1.3, 3.4,max(pbcseq.survival.head$bili)))
pbcseq.survival.head$ageGrp <- cut(pbcseq.survival.head$age, breaks = c(min(pbcseq.survival.head$age),50,max(pbcseq.survival.head$age)))


fit2 <- survfit(Surv(years, status2) ~ ageGrp + biliGrp,
                 data = pbcseq.survival.head)
ggsurv <- ggsurvplot(fit2, fun = "event", conf.int = TRUE,
                     ggtheme = theme_bw())

ggsurv$plot +theme_bw() +
  theme (legend.position = "right")+
  facet_grid(biliGrp ~ ageGrp)

pbcseq.survival.tail$biliGrp <- cut(pbcseq.survival.tail$bili, breaks = c(min(pbcseq.survival.tail$bili), .8, 1.3, 3.4,max(pbcseq.survival.tail$bili)))
pbcseq.survival.tail$ageGrp <- cut(pbcseq.survival.tail$age, breaks = c(min(pbcseq.survival.tail$age),50,max(pbcseq.survival.tail$age)))


fit2 <- survfit(Surv(years, status2) ~ ageGrp + biliGrp,
                data = pbcseq.survival.tail)
ggsurv <- ggsurvplot(fit2, fun = "event", conf.int = TRUE,
                     ggtheme = theme_bw())

ggsurv$plot +theme_bw() +
  theme (legend.position = "right")+
  facet_grid(biliGrp ~ ageGrp)


par(mfrow=c(1,1))
pbcseq.survival.head$status3 <- as.numeric(pbcseq.survival.head$status != "alive")
fit <- survfit(Surv(years, status3) ~ 1, data = pbcseq.survival.head, etype = status)
plot(fit, fun = "event", col = c("black", "red"), lwd = 2, ylim = c(0, 0.7),
     xlab = "Follow-Up (years)", ylab = "Cumulative Incidence")
legend("topleft", levels(as.factor(pbcseq.survival.head$status))[2:3], lty = 1, lwd = 2, col = c("black", "red"),
       bty = "n")



fit1 <- coxph(Surv(years, status == "transplanted") ~ trt + age,
              data = pbcseq.survival.head)
summary(fit1)
fit2 <- coxph(Surv(years, status == "dead") ~ trt + age,
              data = pbcseq.survival.head)
summary(fit2)
