library(survival)
library(survminer)
library(stringr)
data = read.csv("C:/Users/91735/Documents/BVP material/7_cancer genomics/Breast_Cancer.csv")
data[data$Status == "Alive","Status"]<-1
data[data$Status == "Dead","Status"] <- 0

head(data)

colnames(data)[15] = 'OS'
fit=survfit(Surv(OS,Status)~Grade,data=data)
fit
ggsurvplot(fit,data=data)
ggsurvplot(fit,data=data,surv.median.line='hv')
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T)
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T,risk.table = T)


fit=survfit(Surv(OS,Status)~Race + Grade,data=data)
fit
ggsurvplot(fit,data=data)
ggsurvplot(fit,data=data,surv.median.line='hv')
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T)
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T,risk.table = T)
