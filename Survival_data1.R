library(survival)
library(survminer)
library(stringr)
data = read.csv("C:/Users/91735/Documents/BVP material/7_cancer genomics/pancancerInfo.csv")
colnames(data)[10] = 'OS'

fit=survfit(Surv(OS,Event)~gender,data=data)
fit
ggsurvplot(fit,data=data)
ggsurvplot(fit,data=data,surv.median.line='hv')
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T)
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T,risk.table = T)




fit=survfit(Surv(OS,Event)~race + gender,data=data)
fit
ggsurvplot(fit,data=data)
ggsurvplot(fit,data=data,surv.median.line='hv')
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T)
ggsurvplot(fit,data=data,surv.median.line='hv',pval = T,risk.table = T)