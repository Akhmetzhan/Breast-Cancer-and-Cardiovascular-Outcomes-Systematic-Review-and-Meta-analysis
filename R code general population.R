library(metafor)
data_af_0 <- read.csv2("atrial_fibrillation.csv", header = TRUE, sep = ",")
data_af_0$log_estimate <- log(data_af_0$hr_0)
data_af_0$log_lower <- log(data_af_0$hr_0_cilb)
data_af_0$log_upper <- log(data_af_0$hr_0_ciub)
data_af_0$se <- (data_af_0$log_upper - data_af_0$log_lower)/(2*qnorm(0.975))
data_af_0$var <- data_af_0$se^2
#meta-analysis for the risk of AF as compared to general population in a period 0 to 3 months from the diagnosis of breast cancer
meta_af_0 <- rma(log_estimate, vi = var, data=data_af_0, method="REML", slab = paste(data_af_0$Author, data_af_0$Year, sep = ", "))
predict(meta_af_0, transf = exp)
#sensitivity analysis with KNHA adjustment
meta_af_0_knha <- rma(log_estimate, vi = var, data=data_af_0, method="REML", test = "knha", slab = paste(data_af_0$Author, data_af_0$Year, sep = ", "))
predict(meta_af_0_knha, transf = exp)


#maximum likilehood estimation
meta_af_0_ML <- rma(log_estimate, vi = var, data=data_af_0, method="ML", slab = paste(data_af_0$Author, data_af_0$Year, sep = ", "))
predict(meta_af_0_ML, transf = exp)
#sensitivity analysis with KNHA adjustment
meta_af_0_ML_knha <- rma(log_estimate, vi = var, data=data_af_0, method="ML", test = "knha", slab = paste(data_af_0$Author, data_af_0$Year, sep = ", "))
predict(meta_af_0_ML_knha, transf = exp)



#leave-out-sensitivity analyses
Sens_REML_af_0 <- leave1out(meta_af_0, transf = exp, digits = 2)
Sens_ML_af_0 <- leave1out(meta_af_0_ML, transf = exp, digits = 2)
Sens_REML_af_0
Sens_ML_af_0



data_af_1 <- read.csv2("atrial_fibrillation.csv", header = TRUE, sep = ",")

data_af_1$log_estimate <- log(data_af_1$hr_1)
data_af_1$log_lower <- log(data_af_1$hr_1_cilb)
data_af_1$log_upper <- log(data_af_1$hr_1_ciub)
data_af_1$se <- (data_af_1$log_upper - data_af_1$log_lower)/(2*qnorm(0.975))
data_af_1$var <- data_af_1$se^2
#meta-analysis for the risk of AF as compared to general population in a period 2 to 3 years from the diagnosis of breast cancer

meta_af_1 <- rma(log_estimate, vi = var, data=data_af_1, method="REML", slab = paste(data_af_1$Author, data_af_1$Year, sep = ", "))
predict(meta_af_1, transf = exp)

#sensitivity analysis with KNHA adjustment
meta_af_1_knha <- rma(log_estimate, vi = var, data=data_af_1, method="REML", test = "knha", slab = paste(data_af_1$Author, data_af_1$Year, sep = ", "))
predict(meta_af_1_knha, transf = exp)


#maximum likelihood estimation
meta_af_1_ML <- rma(log_estimate, vi = var, data=data_af_1, method="ML", slab = paste(data_af_1$Author, data_af_1$Year, sep = ", "))
predict(meta_af_1_ML, transf = exp)

#sensitivity analysis with KNHA adjustment
meta_af_1_ML_knha <- rma(log_estimate, vi = var, data=data_af_1, test = "knha", method="ML", slab = paste(data_af_1$Author, data_af_1$Year, sep = ", "))
predict(meta_af_1_ML_knha, transf = exp)



Sens_ML_af_1 <- leave1out(meta_af_1, transf = exp, digits = 2)
Sens_ML_af_1_ML <- leave1out(meta_af_1_ML, transf = exp, digits = 2)
#leave-one-out meta-analysis
Sens_ML_af_1
Sens_ML_af_1_ML


#forest plot for AF, comparison with general population - Figure 4.
par(mfrow = c(2,1), mai = c(0.2, 0.8,0.2, 0.4), omi = c(0.2,0.2,0.2,0.2))

forest(meta_af_0, atransf = exp, xlab= text(""), slab = paste(data_af_0$Author, data_af_0$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_af_0$pmid), ilab.xpos=c(-1.75), cex=0.75, ylim=c(-2, meta_af_0$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_af_0$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_af_0$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_af_0$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_af_0$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.75), meta_af_0$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_af_0_knha, atransf = exp, mlab = paste("RE, REML, KNHA", " (Q = ", formatC(meta_af_0_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_af_0_knha$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_af_0_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_af_0_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -2)
addpoly(meta_af_0_ML, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_af_0_ML$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_af_0_ML$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_af_0_ML$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_af_0_ML$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)
addpoly(meta_af_0_ML_knha, atransf = exp, mlab = paste("RE, ML, KNHA", " (Q = ", formatC(meta_af_0_ML_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_af_0_ML_knha$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_af_0_ML_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_af_0_ML_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -4)


title(outer=TRUE,adj=0.5, main="A. Atrial fibrillation, 0-3 months after breast cancer diagnosis.", cex=0.75, col="black",font=2,line=-1.0)

forest(meta_af_1, atransf = exp, xlab= text(""), slab = paste(data_af_1$Author, data_af_1$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_af_1$pmid), ilab.xpos=c(-1.75), cex=1.0, ylim=c(-4, meta_af_1$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_af_1$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_af_1$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_af_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_af_1$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.75), meta_af_1$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_af_1_knha, atransf = exp, mlab = paste("RE, REML, KNHA", " (Q = ", formatC(meta_af_1_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_af_1_knha$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_af_1_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_af_1_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -2)
addpoly(meta_af_1_ML, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_af_1_ML$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_af_1_ML$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_af_1_ML$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_af_1_ML$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)
addpoly(meta_af_1_ML_knha, atransf = exp, mlab = paste("RE, ML, KNHA", " (Q = ", formatC(meta_af_1_ML_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_af_1_ML_knha$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_af_1_ML_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_af_1_ML_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -4)


title(outer=TRUE,adj=0.5, main="B. Atrial fibrillation, 3 months - 3 years after breast cancer diagnosis.", cex=0.75, col="black",font=2,line=-24.0)













#meta-analyses for coronary artery disease as an outcome
data_cad <- read.csv2("cad.csv", header = TRUE, sep = ",")
data_cad$log_estimate <- log(data_cad$hr_0)
data_cad$log_lower <- log(data_cad$hr_0_cilb)
data_cad$log_upper <- log(data_cad$hr_0_ciub)
data_cad$se <- (data_cad$log_upper - data_cad$log_lower)/(2*qnorm(0.975))
data_cad$var <- data_cad$se^2
#meta-analysis for the risk of CAD as compared to general population in a period 0 to 5 years from the diagnosis of breast cancer

meta_cad <- rma(log_estimate, vi = var, data=data_cad, method="REML", slab = paste(data_cad$Author, data_cad$Year, sep = ", "))
predict(meta_cad, transf = exp)
#maximum likilehood estimation
meta_cad_ML <- rma(log_estimate, vi = var, data=data_cad, method="ML", slab = paste(data_cad$Author, data_cad$Year, sep = ", "))
predict(meta_cad_ML, transf = exp)

data_cad$log_estimate_1 <- log(data_cad$hr_1)
data_cad$log_lower_1 <- log(data_cad$hr_1_cilb)
data_cad$log_upper_1 <- log(data_cad$hr_1_ciub)
data_cad$se_1 <- (data_cad$log_upper_1 - data_cad$log_lower_1)/(2*qnorm(0.975))
data_cad$var_1 <- data_cad$se_1^2
#meta-analysis for the risk of CAD as compared to general population in a period 5 to 8 years from the diagnosis of breast cancer

meta_cad_1 <- rma(log_estimate_1, vi = var_1, data=data_cad, method="REML", slab = paste(data_cad$Author, data_cad$Year, sep = ", "))
predict(meta_cad_1, transf = exp)
#maximum likilehood estimation
meta_cad_ML_1 <- rma(log_estimate_1, vi = var_1, data=data_cad, method="ML", slab = paste(data_cad$Author, data_cad$Year, sep = ", "))
predict(meta_cad_ML_1, transf = exp)


#forest plot for CAD, comparison with general population - Supplementary Figure S1

par(mfrow = c(2,1), mai = c(0.2, 0.8,0.2, 0.4), omi = c(0.2,0.2,0.2,0.2))

forest(meta_cad, atransf = exp, xlab= text(""), slab = paste(data_cad$Author, data_cad$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_cad$pmid), ilab.xpos=c(-1.75), cex=0.75, ylim=c(-3, meta_cad$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_cad$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_cad$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cad$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cad$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.75), meta_cad$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_cad_ML, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_cad_ML$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_cad_ML$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cad_ML$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cad_ML$tau2, digits=2, format="f"), ")"),  cex=1.0, -2)

title(outer=TRUE,adj=0.5, main="A. Coronary artery diseas, 0-5 years after breast cancer diagnosis.", cex=0.75, col="black",font=2,line=-1.0)

forest(meta_cad_1, atransf = exp, xlab= text(""), slab = paste(data_cad$Author, data_cad$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_cad$pmid), ilab.xpos=c(-1.75), cex=1.0, ylim=c(-3, meta_cad_1$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_cad_1$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_cad_1$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cad_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cad_1$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.75), meta_cad_1$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_cad_ML_1, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_cad_ML_1$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_cad_ML_1$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cad_ML_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cad_ML_1$tau2, digits=2, format="f"), ")"),  cex=1.0, -2)


title(outer=TRUE,adj=0.5, main="B. Coronary artery disease, 5-8 years after breast cancer diagnosis.", cex=0.75, col="black",font=2,line=-24.0)












#meta-analyses for stroke as an outcome
data_stroke <- read.csv2("stroke_with_navi.csv", header = TRUE, sep = ",")
data_stroke$log_estimate_5 <- log(data_stroke$hr_5)
data_stroke$log_lower_5 <- log(data_stroke$hr_5_cilb)
data_stroke$log_upper_5 <- log(data_stroke$hr_5_ciub)
data_stroke$se_5 <- (data_stroke$log_upper_5 - data_stroke$log_lower_5)/(2*qnorm(0.975))
data_stroke$var_5 <- data_stroke$se_5^2
#meta-analysis for the risk of stroke as compared to general population in a period 0 to 8 years from the diagnosis of breast cancer

meta_stroke_5 <- rma(log_estimate_5, vi = var_5, data=data_stroke, method="REML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_5, transf = exp)
#maximum likilehood estimation
meta_stroke_ML_5 <- rma(log_estimate_5, vi = var_5, data=data_stroke, method="ML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_ML_5, transf = exp)


data_stroke$log_estimate_4 <- log(data_stroke$hr_4)
data_stroke$log_lower_4 <- log(data_stroke$hr_4_cilb)
data_stroke$log_upper_4 <- log(data_stroke$hr_4_ciub)
data_stroke$se_4 <- (data_stroke$log_upper_4 - data_stroke$log_lower_4)/(2*qnorm(0.975))
data_stroke$var_4 <- data_stroke$se_4^2
#meta-analysis for the risk of stroke as compared to general population in a period 9 to 12 months from the diagnosis of breast cancer

meta_stroke_4 <- rma(log_estimate_4, vi = var_4, data=data_stroke, method="REML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_4, transf = exp)
#maximum likilehood estimation
meta_stroke_ML_4 <- rma(log_estimate_4, vi = var_4, data=data_stroke, method="ML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_ML_4, transf = exp)


data_stroke$log_estimate_3 <- log(data_stroke$hr_3)
data_stroke$log_lower_3 <- log(data_stroke$hr_3_cilb)
data_stroke$log_upper_3 <- log(data_stroke$hr_3_ciub)
data_stroke$se_3 <- (data_stroke$log_upper_3 - data_stroke$log_lower_3)/(2*qnorm(0.975))
data_stroke$var_3 <- data_stroke$se_3^2
#meta-analysis for the risk of stroke as compared to general population in a period 6 to 9 months from the diagnosis of breast cancer

meta_stroke_3 <- rma(log_estimate_3, vi = var_3, data=data_stroke, method="REML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_3, transf = exp)
#maximum likilehood estimation
meta_stroke_ML_3 <- rma(log_estimate_3, vi = var_3, data=data_stroke, method="ML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_ML_3, transf = exp)



data_stroke$log_estimate_2 <- log(data_stroke$hr_2)
data_stroke$log_lower_2 <- log(data_stroke$hr_2_cilb)
data_stroke$log_upper_2 <- log(data_stroke$hr_2_ciub)
data_stroke$se_2 <- (data_stroke$log_upper_2 - data_stroke$log_lower_2)/(2*qnorm(0.975))
data_stroke$var_2 <- data_stroke$se_2^2
#meta-analysis for the risk of stroke as compared to general population in a period 3 to 6 months from the diagnosis of breast cancer

meta_stroke_2 <- rma(log_estimate_2, vi = var_2, data=data_stroke, method="REML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_2, transf = exp)
#maximum likilehood estimation
meta_stroke_ML_2 <- rma(log_estimate_2, vi = var_2, data=data_stroke, method="ML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_ML_2, transf = exp)



data_stroke$log_estimate_1 <- log(data_stroke$hr_1)
data_stroke$log_lower_1 <- log(data_stroke$hr_1_cilb)
data_stroke$log_upper_1 <- log(data_stroke$hr_1_ciub)
data_stroke$se_1 <- (data_stroke$log_upper_1 - data_stroke$log_lower_1)/(2*qnorm(0.975))
data_stroke$var_1 <- data_stroke$se_1^2
#meta-analysis for the risk of stroke as compared to general population in a period 1 to 3 months from the diagnosis of breast cancer

meta_stroke_1 <- rma(log_estimate_1, vi = var_1, data=data_stroke, method="REML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_1, transf = exp)
#maximum likilehood estimation
meta_stroke_ML_1 <- rma(log_estimate_1, vi = var_1, data=data_stroke, method="ML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_ML_1, transf = exp)




data_stroke$log_estimate_0 <- log(data_stroke$hr_0)
data_stroke$log_lower_0 <- log(data_stroke$hr_0_cilb)
data_stroke$log_upper_0 <- log(data_stroke$hr_0_ciub)
data_stroke$se_0 <- (data_stroke$log_upper_0 - data_stroke$log_lower_0)/(2*qnorm(0.975))
data_stroke$var_0 <- data_stroke$se_0^2
#meta-analysis for the risk of stroke as compared to general population in a period 0 to 1 month from the diagnosis of breast cancer

meta_stroke_0 <- rma(log_estimate_0, vi = var_0, data=data_stroke, method="REML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_0, transf = exp)
#maximum likilehood estimation
meta_stroke_ML_0 <- rma(log_estimate_0, vi = var_0, data=data_stroke, method="ML", slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "))
predict(meta_stroke_ML_0, transf = exp)



#forest plot for stroke, comparison with general population - Supplementary Figure S3
par(mfrow = c(3,2), mai = c(0.2, 0.8,0.2, 0.4), omi = c(0.2,0.2,0.2,0.2))

forest(meta_stroke_0, atransf = exp, xlab= text(""), slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_stroke$pmid), ilab.xpos=c(-1.4), cex=1.17, ylim=c(-3.5, meta_stroke_0$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_stroke_0$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_stroke_0$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_0$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_0$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_stroke_0$k + 2, c("PMID"))
par(cex = 0.75, font = 2)
#text for ML RE model
addpoly(meta_stroke_ML_0, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_stroke_ML_0$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_stroke_ML_0$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_ML_0$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_ML_0$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)


title(outer=TRUE,adj=0.25, main="A.", cex=0.75, col="black",font=2,line=-1.0)

forest(meta_stroke_1, atransf = exp, xlab= text(""), slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_stroke$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_stroke_1$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_stroke_1$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_stroke_1$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_1$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_stroke_1$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_stroke_ML_1, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_stroke_ML_1$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_stroke_ML_1$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_ML_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_ML_1$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)


title(outer=TRUE,adj=0.75, main="B.", cex=0.75, col="black",font=2,line=-1.0)

forest(meta_stroke_2, atransf = exp, xlab= text(""), slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_stroke$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_stroke_2$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_stroke_2$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_stroke_2$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_2$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_2$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_stroke_2$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_stroke_ML_2, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_stroke_ML_2$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_stroke_ML_2$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_ML_2$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_ML_2$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)


title(outer=TRUE,adj=0.25, main="C.", cex=0.75, col="black",font=2,line=-16.0)


forest(meta_stroke_3, atransf = exp, xlab= text(""), slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_stroke$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_stroke_3$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_stroke_3$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_stroke_3$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_3$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_3$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_stroke_3$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_stroke_ML_3, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_stroke_ML_3$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_stroke_ML_3$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_ML_3$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_ML_3$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)

title(outer=TRUE,adj=0.75, main="D.", cex=0.75, col="black",font=2,line=-16.0)


forest(meta_stroke_4, atransf = exp, xlab= text(""), slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_stroke$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_stroke_4$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_stroke_4$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_stroke_4$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_4$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_4$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_stroke_4$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_stroke_ML_4, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_stroke_ML_4$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_stroke_ML_4$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_ML_4$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_ML_4$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)

title(outer=TRUE,adj=0.25, main="E.", cex=0.75, col="black",font=2,line=-30.0)



forest(meta_stroke_5, atransf = exp, xlab= text(""), slab = paste(data_stroke$Author, data_stroke$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_stroke$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_stroke_5$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_stroke_5$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_stroke_5$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_5$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_5$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_stroke_5$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_stroke_ML_5, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_stroke_ML_5$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_stroke_ML_5$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_stroke_ML_5$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_stroke_ML_5$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)

title(outer=TRUE,adj=0.75, main="F.", cex=0.75, col="black",font=2,line=-30.0)










#meta-analyses for heart failure as an outcome
data_hf <- read.csv2("hf.csv", header = TRUE, sep = ",")
data_hf$log_estimate <- log(data_hf$hr_0)
data_hf$log_lower <- log(data_hf$hr_0_cilb)
data_hf$log_upper <- log(data_hf$hr_0_ciub)
data_hf$se <- (data_hf$log_upper - data_hf$log_lower)/(2*qnorm(0.975))
data_hf$var <- data_hf$se^2
#meta-analysis for the risk of heart failure as compared to general population in a period 0 to 1 year from the diagnosis of breast cancer

meta_hf <- rma(log_estimate, vi = var, data=data_hf, method="REML", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf, transf = exp)
meta_hf_knha <- rma(log_estimate, vi = var, data=data_hf, method="REML", test="knha", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_knha, transf = exp)
#maximum likilehood estimation
meta_hf_ML <- rma(log_estimate, vi = var, data=data_hf, method="ML", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_ML, transf = exp)
#leave-one-out sensitivity analyses
Sens_ML_hf <- leave1out(meta_hf, transf = exp, digits = 2)
Sens_ML_hf_ML <- leave1out(meta_hf_ML, transf = exp, digits = 2)
#sensitivity analysis with KNHA as a test estimator
meta_hf_ML_knha <- rma(log_estimate, vi = var, data=data_hf, method="ML", test="knha", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_ML_knha, transf = exp)
Sens_ML_hf
Sens_ML_hf_ML



data_hf$log_estimate_1 <- log(data_hf$hr_1)
data_hf$log_lower_1 <- log(data_hf$hr_1_cilb)
data_hf$log_upper_1 <- log(data_hf$hr_1_ciub)
data_hf$se_1 <- (data_hf$log_upper_1 - data_hf$log_lower_1)/(2*qnorm(0.975))
data_hf$var_1 <- data_hf$se_1^2
#meta-analysis for the risk of heart failure as compared to general population in a period 1 to 2 years from the diagnosis of breast cancer

meta_hf_1 <- rma(log_estimate_1, vi = var_1, data=data_hf, method="REML", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_1, transf = exp)
#sensitivity analysis with KNHA as a test estimator
meta_hf_1_knha <- rma(log_estimate_1, vi = var_1, data=data_hf, method="REML", test="knha",  slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_1_knha, transf = exp)
#maximum likilehood estimation
meta_hf_ML_1 <- rma(log_estimate_1, vi = var_1, data=data_hf, method="ML", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_ML_1, transf = exp)
#sensitivity analysis with KNHA as a test estimator
meta_hf_ML_1_knha <- rma(log_estimate_1, vi = var_1, data=data_hf, method="ML", test="knha", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_ML_1_knha, transf = exp)
#leave-one-out sensitivity analyses
Sens_ML_hf_1 <- leave1out(meta_hf_1, transf = exp, digits = 2)
Sens_ML_hf_ML_1 <- leave1out(meta_hf_ML_1, transf = exp, digits = 2)
Sens_ML_hf_1
Sens_ML_hf_ML_1 




data_hf$log_estimate_2 <- log(data_hf$hr_2)
data_hf$log_lower_2 <- log(data_hf$hr_2_cilb)
data_hf$log_upper_2 <- log(data_hf$hr_2_ciub)
data_hf$se_2 <- (data_hf$log_upper_2 - data_hf$log_lower_2)/(2*qnorm(0.975))
data_hf$var_2 <- data_hf$se_2^2
#meta-analysis for the risk of heart failure as compared to general population in a period 2 to 5 years from the diagnosis of breast cancer

meta_hf_2 <- rma(log_estimate_2, vi = var_2, data=data_hf, method="REML", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_2, transf = exp)
#sensitivity analysis with KNHA as a test estimator
meta_hf_2_knha <- rma(log_estimate_2, vi = var_2, data=data_hf, method="REML", test="knha", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_2_knha, transf = exp)

#maximum likilehood estimation
meta_hf_ML_2 <- rma(log_estimate_2, vi = var_2, data=data_hf, method="ML", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_ML_2, transf = exp)
#sensitivity analysis with KNHA as a test estimator
meta_hf_ML_2_knha <- rma(log_estimate_2, vi = var_2, data=data_hf, method="ML", test="knha",  slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_ML_2_knha, transf = exp)
#leave-one-out sensitivity analyses
Sens_ML_hf_2 <- leave1out(meta_hf_2, transf = exp, digits = 2)
Sens_ML_hf_ML_2 <- leave1out(meta_hf_ML_2, transf = exp, digits = 2)
Sens_ML_hf_2
Sens_ML_hf_ML_2


data_hf$log_estimate_3 <- log(data_hf$hr_3)
data_hf$log_lower_3 <- log(data_hf$hr_3_cilb)
data_hf$log_upper_3 <- log(data_hf$hr_3_ciub)
data_hf$se_3 <- (data_hf$log_upper_3 - data_hf$log_lower_3)/(2*qnorm(0.975))
data_hf$var_3 <- data_hf$se_3^2
#meta-analysis for the risk of heart failure as compared to general population in a period 5 to 10 years from the diagnosis of breast cancer

meta_hf_3 <- rma(log_estimate_3, vi = var_3, data=data_hf, method="REML", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_3, transf = exp)
#sensitivity analysis with KNHA as a test estimator
meta_hf_3_knha <- rma(log_estimate_3, vi = var_3, data=data_hf, method="REML", test="knha", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_3_knha, transf = exp)

#maximum likilehood estimation
meta_hf_ML_3 <- rma(log_estimate_3, vi = var_3, data=data_hf, method="ML", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_ML_3, transf = exp)
#sensitivity analysis with KNHA as a test estimator
meta_hf_ML_3_knha <- rma(log_estimate_3, vi = var_3, data=data_hf, method="ML", test="knha", slab = paste(data_hf$Author, data_hf$Year, sep = ", "))
predict(meta_hf_ML_3_knha, transf = exp)
#leave-one-out sensitivity analyses
Sens_ML_hf_3 <- leave1out(meta_hf_3, transf = exp, digits = 2)
Sens_ML_hf_ML_3 <- leave1out(meta_hf_ML_3, transf = exp, digits = 2)
Sens_ML_hf_3
Sens_ML_hf_ML_3


#publication bias assessment - less than 10 studies
rtf_meta_hf_1 <- trimfill(meta_hf_1)
regtest_meta_hf_1 <- regtest(meta_hf_1, model = "rma")
ranktest_meta_hf_1 <- ranktest(meta_hf_1)

rtf_meta_hf_2 <- trimfill(meta_hf_2)
regtest_meta_hf_2 <- regtest(meta_hf_2, model = "rma")
ranktest_meta_hf_2 <- ranktest(meta_hf_2)

rtf_meta_hf_3 <- trimfill(meta_hf_3)
regtest_meta_hf_3 <- regtest(meta_hf_3, model = "rma")
ranktest_meta_hf_3 <- ranktest(meta_hf_3)




#forest plot for HF, comparison with general population - Figure 2
par(mfrow = c(2,2), mai = c(0.2, 0.8,0.2, 0.4), omi = c(0.2,0.2,0.2,0.2))

forest(meta_hf, atransf = exp, xlab= text(""), slab = paste(data_hf$Author, data_hf$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_hf$pmid), ilab.xpos=c(-1.5), cex=1.0, ylim=c(-6, meta_hf$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_hf$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_hf$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_hf$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.5), meta_hf$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_hf_knha, atransf = exp, mlab = paste("RE, REML, KNHA", " (Q = ", formatC(meta_hf_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_knha$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -2.5)
addpoly(meta_hf_ML, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_hf_ML$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_ML$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_ML$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_ML$tau2, digits=2, format="f"), ")"),  cex=1.0, -4)
addpoly(meta_hf_ML_knha, atransf = exp, mlab = paste("RE, ML, KNHA", " (Q = ", formatC(meta_hf_ML_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_ML_knha$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_ML_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_ML_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -5.7)


title(outer=TRUE,adj=0.25, main="A.", cex=0.75, col="black",font=2,line=-1.0)

forest(meta_hf_1, atransf = exp, xlab= text(""), slab = paste(data_hf$Author, data_hf$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_hf$pmid), ilab.xpos=c(-1.5), cex=1.0, ylim=c(-6, meta_hf_1$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_hf_1$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_hf_1$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_hf_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_1$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.5), meta_hf_1$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_hf_1_knha, atransf = exp, mlab = paste("RE, REML, KNHA", " (Q = ", formatC(meta_hf_1_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_1_knha$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_1_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_1_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -2.5)
addpoly(meta_hf_ML_1, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_hf_ML_1$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_ML_1$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_ML_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_ML_1$tau2, digits=2, format="f"), ")"),  cex=1.0, -4)
addpoly(meta_hf_ML_1_knha, atransf = exp, mlab = paste("RE, ML, KNHA", " (Q = ", formatC(meta_hf_ML_1_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_ML_1_knha$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_ML_1_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_ML_1_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -5.7)


title(outer=TRUE,adj=0.75, main="B.", cex=0.75, col="black",font=2,line=-1.0)



forest(meta_hf_2, atransf = exp, xlab= text(""), slab = paste(data_hf$Author, data_hf$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_hf$pmid), ilab.xpos=c(-1.5), cex=1.0, ylim=c(-6, meta_hf_2$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_hf_2$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_hf_2$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_hf_2$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_2$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.5), meta_hf_2$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_hf_2_knha, atransf = exp, mlab = paste("RE, REML, KNHA", " (Q = ", formatC(meta_hf_2_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_2_knha$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_2_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_2_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -2.5)
addpoly(meta_hf_ML_2, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_hf_ML_2$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_ML_2$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_ML_2$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_ML_2$tau2, digits=2, format="f"), ")"),  cex=1.0, -4)
addpoly(meta_hf_ML_2_knha, atransf = exp, mlab = paste("RE, ML, KNHA", " (Q = ", formatC(meta_hf_ML_2_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_ML_2_knha$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_ML_2_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_ML_2_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -5.7)


title(outer=TRUE,adj=0.25, main="C.", cex=0.75, col="black",font=2,line=-24.0)


forest(meta_hf_3, atransf = exp, xlab= text(""), slab = paste(data_hf$Author, data_hf$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_hf$pmid), ilab.xpos=c(-1.5), cex=1.0, ylim=c(-6, meta_hf_3$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_hf_3$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_hf_3$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_hf_3$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_3$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.5), meta_hf_3$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_hf_3_knha, atransf = exp, mlab = paste("RE, REML, KNHA", " (Q = ", formatC(meta_hf_3_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_3_knha$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_3_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_3_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -2.5)
addpoly(meta_hf_ML_3, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_hf_ML_3$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_ML_3$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_ML_3$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_ML_3$tau2, digits=2, format="f"), ")"),  cex=1.0, -4)
addpoly(meta_hf_ML_3_knha, atransf = exp, mlab = paste("RE, ML, KNHA", " (Q = ", formatC(meta_hf_ML_3_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_hf_ML_3_knha$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_hf_ML_3_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_hf_ML_3_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -5.7)


title(outer=TRUE,adj=0.75, main="D.", cex=0.75, col="black",font=2,line=-24.0)














#meta-analyses for cardiovascular death as an outcome
data_cvdeath_hr <- read.csv2("cvdeath.csv", header = TRUE, sep = ",")

data_cvdeath_hr$log_estimate_0 <- log(data_cvdeath_hr$hr_0)
data_cvdeath_hr$log_lower_0 <- log(data_cvdeath_hr$hr_0_cilb)
data_cvdeath_hr$log_upper_0 <- log(data_cvdeath_hr$hr_0_ciub)
data_cvdeath_hr$se_0 <- (data_cvdeath_hr$log_upper_0 - data_cvdeath_hr$log_lower_0)/(2*qnorm(0.975))
data_cvdeath_hr$var_0 <- data_cvdeath_hr$se_0^2
#meta-analysis for the risk of CV death as compared to general population in a period 0 to 5 years from the diagnosis of breast cancer

meta_cvdeath_hr_0 <- rma(log_estimate_0, vi = var_0, data=data_cvdeath_hr, method="REML", slab = paste(data_cvdeath_hr$Author, data_cvdeath_hr$Year, sep = ", "))
predict(meta_cvdeath_hr_0, transf = exp)
#sensitivity analysis with KNHA as a test estimator
meta_cvdeath_hr_0_knha <- rma(log_estimate_0, vi = var_0, data=data_cvdeath_hr, method="REML", test="knha", slab = paste(data_cvdeath_hr$Author, data_cvdeath_hr$Year, 
                                                                                                                         
                                                                                                                         sep = ", "))
predict(meta_cvdeath_hr_0_knha, transf = exp)
#maximum likilehood estimation
meta_cvdeath_hr_0_ML <- rma(log_estimate_0, vi = var_0, data=data_cvdeath_hr, method="ML", slab = paste(data_cvdeath_hr$Author, data_cvdeath_hr$Year, sep = ", "))
predict(meta_cvdeath_hr_0_ML, transf = exp)
#sensitivity analysis with KNHA as a test estimator
meta_cvdeath_hr_0_ML_knha <- rma(log_estimate_0, vi = var_0, data=data_cvdeath_hr, method="ML", test="knha", slab = paste(data_cvdeath_hr$Author, data_cvdeath_hr$Year, 
                                                                                                                          
                                                                                                                          sep = ", "))
predict(meta_cvdeath_hr_0_ML_knha, transf = exp)
#leave-one-out sensitivity analyses
Sens_cvdeath_hr_0 <- leave1out(meta_cvdeath_hr_0, transf = exp, digits = 2)
Sens_cvdeath_hr_0_ML <- leave1out(meta_cvdeath_hr_0_ML, transf = exp, digits = 2)
Sens_cvdeath_hr_0
Sens_cvdeath_hr_0_ML

data_cvdeath_hr$log_estimate_1 <- log(data_cvdeath_hr$hr_1)
data_cvdeath_hr$log_lower_1 <- log(data_cvdeath_hr$hr_1_cilb)
data_cvdeath_hr$log_upper_1 <- log(data_cvdeath_hr$hr_1_ciub)
data_cvdeath_hr$se_1 <- (data_cvdeath_hr$log_upper_1 - data_cvdeath_hr$log_lower_1)/(2*qnorm(0.975))
data_cvdeath_hr$var_1 <- data_cvdeath_hr$se_1^2
#meta-analysis for the risk of CV death as compared to general population in a period 8 to 11 years from the diagnosis of breast cancer

meta_cvdeath_hr_1 <- rma(log_estimate_1, vi = var_1, data=data_cvdeath_hr, method="REML", slab = paste(data_cvdeath_hr$Author, data_cvdeath_hr$Year, sep = ", "))
predict(meta_cvdeath_hr_1, transf = exp)
#sensitivity analysis with KNHA as a test estimator
meta_cvdeath_hr_1_knha <- rma(log_estimate_1, vi = var_1, data=data_cvdeath_hr, method="REML", test="knha", slab = paste(data_cvdeath_hr$Author, data_cvdeath_hr$Year, 
                                                                                                                         
                                                                                                                         sep = ", "))
predict(meta_cvdeath_hr_1_knha, transf = exp)

#maximum likilehood estimation
meta_cvdeath_hr_1_ML <- rma(log_estimate_1, vi = var_1, data=data_cvdeath_hr, method="ML", slab = paste(data_cvdeath_hr$Author, data_cvdeath_hr$Year, sep = ", "))
predict(meta_cvdeath_hr_1_ML, transf = exp)
#sensitivity analysis with KNHA as a test estimator
meta_cvdeath_hr_1_ML_knha <- rma(log_estimate_1, vi = var_1, data=data_cvdeath_hr, method="ML", test="knha", slab = paste(data_cvdeath_hr$Author, data_cvdeath_hr$Year, 
                                                                                                                          
                                                                                                                          sep = ", "))
predict(meta_cvdeath_hr_1_ML_knha, transf = exp)
#leave-one-out sensitivity analyses
Sens_cvdeath_hr_1 <- leave1out(meta_cvdeath_hr_1, transf = exp, digits = 2)
Sens_cvdeath_hr_1_ML <- leave1out(meta_cvdeath_hr_1_ML, transf = exp, digits = 2)
Sens_cvdeath_hr_1
Sens_cvdeath_hr_1_ML


#forest plot for CV death, comparison with general population - Figure 2
par(mfrow = c(2,1), mai = c(0.2, 0.8,0.2, 0.4), omi = c(0.2,0.2,0.2,0.2))

forest(meta_cvdeath_hr_0, atransf = exp, xlab= text(""), slab = paste(data_cvdeath_hr$Author, data_cvdeath_hr$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_cvdeath_hr$pmid), ilab.xpos=c(-1.75), cex=0.75, ylim=c(-4, meta_cvdeath_hr_0$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_cvdeath_hr_0$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_cvdeath_hr_0$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cvdeath_hr_0$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cvdeath_hr_0$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.75), meta_cvdeath_hr_0$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_cvdeath_hr_0_knha, atransf = exp, mlab = paste("RE, REML, KNHA", " (Q = ", formatC(meta_cvdeath_hr_0_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_cvdeath_hr_0_knha$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cvdeath_hr_0_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cvdeath_hr_0_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -2)
addpoly(meta_cvdeath_hr_0_ML, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_cvdeath_hr_0_ML$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_cvdeath_hr_0_ML$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cvdeath_hr_0_ML$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cvdeath_hr_0_ML$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)
addpoly(meta_cvdeath_hr_0_ML_knha, atransf = exp, mlab = paste("RE, ML, KNHA", " (Q = ", formatC(meta_cvdeath_hr_0_ML_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_cvdeath_hr_0_ML_knha$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cvdeath_hr_0_ML_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cvdeath_hr_0_ML_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -4)


title(outer=TRUE,adj=0.5, main="A. Cardiovascular death, 0-5 years after breast cancer diagnosis.", cex=0.75, col="black",font=2,line=-1.0)

forest(meta_cvdeath_hr_1, atransf = exp, xlab= text(""), slab = paste(data_cvdeath_hr$Author, data_cvdeath_hr$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_cvdeath_hr$pmid), ilab.xpos=c(-1.75), cex=1.0, ylim=c(-4, meta_cvdeath_hr_1$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_cvdeath_hr_1$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_cvdeath_hr_1$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cvdeath_hr_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cvdeath_hr_1$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.75), meta_cvdeath_hr_1$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_cvdeath_hr_1_knha, atransf = exp, mlab = paste("RE, REML, KNHA", " (Q = ", formatC(meta_cvdeath_hr_1_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_cvdeath_hr_1_knha$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cvdeath_hr_1_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cvdeath_hr_1_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -2)
addpoly(meta_cvdeath_hr_1_ML, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_cvdeath_hr_1_ML$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_cvdeath_hr_1_ML$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cvdeath_hr_1_ML$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cvdeath_hr_1_ML$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)
addpoly(meta_cvdeath_hr_1_ML_knha, atransf = exp, mlab = paste("RE, ML, KNHA", " (Q = ", formatC(meta_cvdeath_hr_1_ML_knha$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_cvdeath_hr_1_ML_knha$QEp, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_cvdeath_hr_1_ML_knha$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_cvdeath_hr_1_ML_knha$tau2, digits=2, format="f"), ")"),  cex=1.0, -4)


title(outer=TRUE,adj=0.5, main="B. Cardiovascular death, 8-11 years after breast cancer diagnosis.", cex=0.75, col="black",font=2,line=-24.0)













#meta-analyses with ischemic stroke as an outcome
data_ischemic_stroke <- read.csv2("ischemic stroke.csv", header = TRUE, sep = ",")

data_ischemic_stroke$log_estimate_0_1 <- log(data_ischemic_stroke$hr_0_1)
data_ischemic_stroke$log_lower_0_1 <- log(data_ischemic_stroke$hr_0_1_cilb)
data_ischemic_stroke$log_upper_0_1 <- log(data_ischemic_stroke$hr_0_1_ciub)
data_ischemic_stroke$se_0_1 <- (data_ischemic_stroke$log_upper_0_1 - data_ischemic_stroke$log_lower_0_1)/(2*qnorm(0.975))
data_ischemic_stroke$var_0_1 <- data_ischemic_stroke$se_0_1^2
#meta-analysis for the risk of ischemic stroke as compared to general population in a period 0 to 1 month from the diagnosis of breast cancer

meta_ischemic_stroke_0_1 <- rma(log_estimate_0_1, vi = var_0_1, data=data_ischemic_stroke, method="REML", slab = paste(data_ischemic_stroke$Author, 
                                                                                                                       
                                                                                                                       data_ischemic_stroke$Year, sep = ", "))
predict(meta_ischemic_stroke_0_1, transf = exp)
#maximum likelihood
meta_ischemic_stroke_0_1_ML <- rma(log_estimate_0_1, vi = var_0_1, data=data_ischemic_stroke, method="ML", slab = paste(data_ischemic_stroke$Author, 
                                                                                                                        
                                                                                                                        data_ischemic_stroke$Year, sep = ", "))
predict(meta_ischemic_stroke_0_1_ML, transf = exp)
#leave-one-out analyses
Sens_ischemic_stroke_0_1 <- leave1out(meta_ischemic_stroke_0_1, transf = exp, digits = 2)
Sens_ischemic_stroke_0_1_ML <- leave1out(meta_ischemic_stroke_0_1_ML, transf = exp, digits = 2)



data_ischemic_stroke$log_estimate_0 <- log(data_ischemic_stroke$hr_0)
data_ischemic_stroke$log_lower_0 <- log(data_ischemic_stroke$hr_0_cilb)
data_ischemic_stroke$log_upper_0 <- log(data_ischemic_stroke$hr_0_ciub)
data_ischemic_stroke$se_0 <- (data_ischemic_stroke$log_upper_0 - data_ischemic_stroke$log_lower_0)/(2*qnorm(0.975))
data_ischemic_stroke$var_0 <- data_ischemic_stroke$se_0^2
#meta-analysis for the risk of ischemic stroke as compared to general population in a period 1 to 3 months from the diagnosis of breast cancer

meta_ischemic_stroke_0 <- rma(log_estimate_0, vi = var_0, data=data_ischemic_stroke, method="REML", slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke
                                                                                                                 
                                                                                                                 $Year, sep = ", "))
predict(meta_ischemic_stroke_0, transf = exp)
#maximum likelihood
meta_ischemic_stroke_0_ML <- rma(log_estimate_0, vi = var_0, data=data_ischemic_stroke, method="ML", slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke
                                                                                                                  
                                                                                                                  $Year, sep = ", "))
predict(meta_ischemic_stroke_0_ML, transf = exp)


data_ischemic_stroke$log_estimate_1 <- log(data_ischemic_stroke$hr_1)
data_ischemic_stroke$log_lower_1 <- log(data_ischemic_stroke$hr_1_cilb)
data_ischemic_stroke$log_upper_1 <- log(data_ischemic_stroke$hr_1_ciub)
data_ischemic_stroke$se_1 <- (data_ischemic_stroke$log_upper_1 - data_ischemic_stroke$log_lower_1)/(2*qnorm(0.975))
data_ischemic_stroke$var_1 <- data_ischemic_stroke$se_1^2
#meta-analysis for the risk of ischemic stroke as compared to general population in a period 3 to 6 months from the diagnosis of breast cancer

meta_ischemic_stroke_1 <- rma(log_estimate_1, vi = var_1, data=data_ischemic_stroke, method="REML", slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke
                                                                                                                 
                                                                                                                 $Year, sep = ", "))
predict(meta_ischemic_stroke_1, transf = exp)
#maximum likelihood
meta_ischemic_stroke_1_ML <- rma(log_estimate_1, vi = var_1, data=data_ischemic_stroke, method="ML", slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke
                                                                                                                  
                                                                                                                  $Year, sep = ", "))
predict(meta_ischemic_stroke_1_ML, transf = exp)

data_ischemic_stroke$log_estimate_2 <- log(data_ischemic_stroke$hr_2)
data_ischemic_stroke$log_lower_2 <- log(data_ischemic_stroke$hr_2_cilb)
data_ischemic_stroke$log_upper_2 <- log(data_ischemic_stroke$hr_2_ciub)
data_ischemic_stroke$se_2 <- (data_ischemic_stroke$log_upper_2 - data_ischemic_stroke$log_lower_2)/(2*qnorm(0.975))
data_ischemic_stroke$var_2 <- data_ischemic_stroke$se_2^2
#meta-analysis for the risk of ischemic stroke as compared to general population in a period 6 to 12 months from the diagnosis of breast cancer

meta_ischemic_stroke_2 <- rma(log_estimate_2, vi = var_2, data=data_ischemic_stroke, method="REML", slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke
                                                                                                                 
                                                                                                                 $Year, sep = ", "))
predict(meta_ischemic_stroke_2, transf = exp)
#maximum likelihood
meta_ischemic_stroke_2_ML <- rma(log_estimate_2, vi = var_2, data=data_ischemic_stroke, method="ML", slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke
                                                                                                                  
                                                                                                                  $Year, sep = ", "))
predict(meta_ischemic_stroke_2_ML, transf = exp)

data_ischemic_stroke$log_estimate_4 <- log(data_ischemic_stroke$hr_4)
data_ischemic_stroke$log_lower_4 <- log(data_ischemic_stroke$hr_4_cilb)
data_ischemic_stroke$log_upper_4 <- log(data_ischemic_stroke$hr_4_ciub)
data_ischemic_stroke$se_4 <- (data_ischemic_stroke$log_upper_4 - data_ischemic_stroke$log_lower_4)/(2*qnorm(0.975))
data_ischemic_stroke$var_4 <- data_ischemic_stroke$se_4^2
#meta-analysis for the risk of ischemic stroke as compared to general population in a period 2 to 3 years from the diagnosis of breast cancer

meta_ischemic_stroke_4 <- rma(log_estimate_4, vi = var_4, data=data_ischemic_stroke, method="REML", slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke
                                                                                                                 
                                                                                                                 $Year, sep = ", "))
predict(meta_ischemic_stroke_4, transf = exp)
#maximum likelihood
meta_ischemic_stroke_4_ML <- rma(log_estimate_4, vi = var_4, data=data_ischemic_stroke, method="ML", slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke
                                                                                                                  
                                                                                                                  $Year, sep = ", "))
predict(meta_ischemic_stroke_4_ML, transf = exp)



#forest plot for ischemic stroke, comparison with general population - Supplementary Figure S4
par(mfrow = c(3,2), mai = c(0.2, 0.8,0.2, 0.4), omi = c(0.2,0.2,0.2,0.2))

forest(meta_ischemic_stroke_0_1, atransf = exp, xlab= text(""), slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_ischemic_stroke$pmid), ilab.xpos=c(-1.4), cex=1.17, ylim=c(-3.5, meta_ischemic_stroke_0_1$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_ischemic_stroke_0_1$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_ischemic_stroke_0_1$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_ischemic_stroke_0_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_ischemic_stroke_0_1$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_ischemic_stroke_0_1$k + 2, c("PMID"))
par(cex = 0.75, font = 2)
#text for ML RE model
addpoly(meta_ischemic_stroke_0_1_ML, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_ischemic_stroke_0_1_ML$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_ischemic_stroke_0_1_ML$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_ischemic_stroke_0_1_ML$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_ischemic_stroke_0_1_ML$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)


title(outer=TRUE,adj=0.25, main="A.", cex=0.75, col="black",font=2,line=-1.0)

forest(meta_ischemic_stroke_0, atransf = exp, xlab= text(""), slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_ischemic_stroke$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_ischemic_stroke_0$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_ischemic_stroke_0$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_ischemic_stroke_0$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_ischemic_stroke_0$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_ischemic_stroke_0$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_ischemic_stroke_0$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_ischemic_stroke_0_ML, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_ischemic_stroke_0_ML$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_ischemic_stroke_0_ML$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_ischemic_stroke_0_ML$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_ischemic_stroke_0_ML$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)


title(outer=TRUE,adj=0.75, main="B.", cex=0.75, col="black",font=2,line=-1.0)

forest(meta_ischemic_stroke_1, atransf = exp, xlab= text(""), slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_ischemic_stroke$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_ischemic_stroke_1$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_ischemic_stroke_1$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_ischemic_stroke_1$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_ischemic_stroke_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_ischemic_stroke_1$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_ischemic_stroke_1$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_ischemic_stroke_1_ML, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_ischemic_stroke_1_ML$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_ischemic_stroke_1_ML$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_ischemic_stroke_1_ML$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_ischemic_stroke_1_ML$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)


title(outer=TRUE,adj=0.25, main="C.", cex=0.75, col="black",font=2,line=-16.0)


forest(meta_ischemic_stroke_2, atransf = exp, xlab= text(""), slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_ischemic_stroke$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_ischemic_stroke_2$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_ischemic_stroke_2$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_ischemic_stroke_2$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_ischemic_stroke_2$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_ischemic_stroke_2$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_ischemic_stroke_2$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_ischemic_stroke_2_ML, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_ischemic_stroke_2_ML$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_ischemic_stroke_2_ML$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_ischemic_stroke_2_ML$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_ischemic_stroke_2_ML$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)

title(outer=TRUE,adj=0.75, main="D.", cex=0.75, col="black",font=2,line=-16.0)


forest(meta_ischemic_stroke_4, atransf = exp, xlab= text(""), slab = paste(data_ischemic_stroke$Author, data_ischemic_stroke$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_ischemic_stroke$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_ischemic_stroke_4$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_ischemic_stroke_4$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_ischemic_stroke_4$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_ischemic_stroke_4$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_ischemic_stroke_4$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_ischemic_stroke_4$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_ischemic_stroke_4_ML, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_ischemic_stroke_4_ML$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_ischemic_stroke_4_ML$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_ischemic_stroke_4_ML$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_ischemic_stroke_4_ML$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)

title(outer=TRUE,adj=0.25, main="E.", cex=0.75, col="black",font=2,line=-30.0)














#meta-analyses for myocardial infarction as an outcome
data_mi_hr <- read.csv2("mi.csv", header = TRUE, sep = ",")

data_mi_hr$log_estimate_0 <- log(data_mi_hr$hr_0)
data_mi_hr$log_lower_0 <- log(data_mi_hr$hr_0_cilb)
data_mi_hr$log_upper_0 <- log(data_mi_hr$hr_0_ciub)
data_mi_hr$se_0 <- (data_mi_hr$log_upper_0 - data_mi_hr$log_lower_0)/(2*qnorm(0.975))
data_mi_hr$var_0 <- data_mi_hr$se_0^2
#meta-analysis for the risk of myocardial infarction as compared to general population in a period 0 to 3 months from the diagnosis of breast cancer

meta_mi_hr_0 <- rma(log_estimate_0, vi = var_0, data=data_mi_hr, method="REML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_0, transf = exp)
#maximum likelihood
meta_mi_hr_ML_0 <- rma(log_estimate_0, vi = var_0, data=data_mi_hr, method="ML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_ML_0, transf = exp)
Sens_mi_hr <- leave1out(meta_mi_hr_0, transf = exp, digits = 2)
Sens_mi_hr_ML <- leave1out(meta_mi_hr_ML_0, transf = exp, digits = 2)


data_mi_hr$log_estimate_1 <- log(data_mi_hr$hr_1)
data_mi_hr$log_lower_1 <- log(data_mi_hr$hr_1_cilb)
data_mi_hr$log_upper_1 <- log(data_mi_hr$hr_1_ciub)
data_mi_hr$se_1 <- (data_mi_hr$log_upper_1 - data_mi_hr$log_lower_1)/(2*qnorm(0.975))
data_mi_hr$var_1 <- data_mi_hr$se_1^2
#meta-analysis for the risk of myocardial infarction as compared to general population in a period 3 to 6 month from the diagnosis of breast cancer

meta_mi_hr_1 <- rma(log_estimate_1, vi = var_1, data=data_mi_hr, method="REML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_1, transf = exp)
#maximum likelihood
meta_mi_hr_ML_1 <- rma(log_estimate_1, vi = var_1, data=data_mi_hr, method="ML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_ML_1, transf = exp)
Sens_mi_hr_1 <- leave1out(meta_mi_hr_1, transf = exp, digits = 2)
Sens_mi_hr_ML_1 <- leave1out(meta_mi_hr_ML_1, transf = exp, digits = 2)


data_mi_hr$log_estimate_2 <- log(data_mi_hr$hr_2)
data_mi_hr$log_lower_2 <- log(data_mi_hr$hr_2_cilb)
data_mi_hr$log_upper_2 <- log(data_mi_hr$hr_2_ciub)
data_mi_hr$se_2 <- (data_mi_hr$log_upper_2 - data_mi_hr$log_lower_2)/(2*qnorm(0.975))
data_mi_hr$var_2 <- data_mi_hr$se_2^2
#meta-analysis for the risk of myocardial infarction as compared to general population in a period 6 to 12 month from the diagnosis of breast cancer

meta_mi_hr_2 <- rma(log_estimate_2, vi = var_2, data=data_mi_hr, method="REML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_2, transf = exp)
#maximum likelihood
meta_mi_hr_ML_2 <- rma(log_estimate_2, vi = var_2, data=data_mi_hr, method="ML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_ML_2, transf = exp)
Sens_mi_hr_2 <- leave1out(meta_mi_hr_2, transf = exp, digits = 2)
Sens_mi_hr_ML_2 <- leave1out(meta_mi_hr_ML_2, transf = exp, digits = 2)



data_mi_hr$log_estimate_3 <- log(data_mi_hr$hr_3)
data_mi_hr$log_lower_3 <- log(data_mi_hr$hr_3_cilb)
data_mi_hr$log_upper_3 <- log(data_mi_hr$hr_3_ciub)
data_mi_hr$se_3 <- (data_mi_hr$log_upper_3 - data_mi_hr$log_lower_3)/(2*qnorm(0.975))
data_mi_hr$var_3 <- data_mi_hr$se_3^2
#meta-analysis for the risk of myocardial infarction as compared to general population in a period 1 to 2 years from the diagnosis of breast cancer

meta_mi_hr_3 <- rma(log_estimate_3, vi = var_3, data=data_mi_hr, method="REML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_3, transf = exp)
#maximum likelihood
meta_mi_hr_ML_3 <- rma(log_estimate_3, vi = var_3, data=data_mi_hr, method="ML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_ML_3, transf = exp)
Sens_mi_hr_3 <- leave1out(meta_mi_hr_3, transf = exp, digits = 2)
Sens_mi_hr_ML_3 <- leave1out(meta_mi_hr_ML_3, transf = exp, digits = 2)

data_mi_hr$log_estimate_4 <- log(data_mi_hr$hr_4)
data_mi_hr$log_lower_4 <- log(data_mi_hr$hr_4_cilb)
data_mi_hr$log_upper_4 <- log(data_mi_hr$hr_4_ciub)
data_mi_hr$se_4 <- (data_mi_hr$log_upper_4 - data_mi_hr$log_lower_4)/(2*qnorm(0.975))
data_mi_hr$var_4 <- data_mi_hr$se_4^2
#meta-analysis for the risk of myocardial infarction as compared to general population in a period 2 to 3 years from the diagnosis of breast cancer

meta_mi_hr_4 <- rma(log_estimate_4, vi = var_4, data=data_mi_hr, method="REML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_4, transf = exp)
#maximum likelihood
meta_mi_hr_ML_4 <- rma(log_estimate_4, vi = var_4, data=data_mi_hr, method="ML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_ML_4, transf = exp)
Sens_mi_hr_4 <- leave1out(meta_mi_hr_4, transf = exp, digits = 2)
Sens_mi_hr_ML_4 <- leave1out(meta_mi_hr_ML_4, transf = exp, digits = 2)





data_mi_hr$log_estimate_5 <- log(data_mi_hr$hr_5)
data_mi_hr$log_lower_5 <- log(data_mi_hr$hr_5_cilb)
data_mi_hr$log_upper_5 <- log(data_mi_hr$hr_5_ciub)
data_mi_hr$se_5 <- (data_mi_hr$log_upper_5 - data_mi_hr$log_lower_5)/(2*qnorm(0.975))
data_mi_hr$var_5 <- data_mi_hr$se_5^2
#meta-analysis for the risk of myocardial infarction as compared to general population in a period 5 to 7 years from the diagnosis of breast cancer

meta_mi_hr_5 <- rma(log_estimate_5, vi = var_5, data=data_mi_hr, method="REML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_5, transf = exp)
#maximum likelihood
meta_mi_hr_ML_5 <- rma(log_estimate_5, vi = var_5, data=data_mi_hr, method="ML", slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "))
predict(meta_mi_hr_ML_5, transf = exp)
Sens_mi_hr_5 <- leave1out(meta_mi_hr_5, transf = exp, digits = 2)
Sens_mi_hr_ML_5 <- leave1out(meta_mi_hr_ML_5, transf = exp, digits = 2)



#forest plot for MI, comparison with general population - Supplementary Figure S2
par(mfrow = c(3,2), mai = c(0.2, 0.8,0.2, 0.4), omi = c(0.2,0.2,0.2,0.2))

forest(meta_mi_hr_0, atransf = exp, xlab= text(""), slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_mi_hr$pmid), ilab.xpos=c(-1.4), cex=1.17, ylim=c(-3.5, meta_mi_hr_0$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_0$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_mi_hr_0$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_0$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_0$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_mi_hr_0$k + 2, c("PMID"))
par(cex = 0.75, font = 2)
#text for ML RE model
addpoly(meta_mi_hr_ML_0, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_ML_0$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_mi_hr_ML_0$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_ML_0$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_ML_0$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)


title(outer=TRUE,adj=0.25, main="A.", cex=0.75, col="black",font=2,line=-1.0)

forest(meta_mi_hr_1, atransf = exp, xlab= text(""), slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_mi_hr$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_mi_hr_1$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_1$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_mi_hr_1$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_1$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_mi_hr_1$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_mi_hr_ML_1, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_ML_1$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_mi_hr_ML_1$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_ML_1$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_ML_1$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)


title(outer=TRUE,adj=0.75, main="B.", cex=0.75, col="black",font=2,line=-1.0)

forest(meta_mi_hr_2, atransf = exp, xlab= text(""), slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_mi_hr$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_mi_hr_2$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_2$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_mi_hr_2$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_2$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_2$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_mi_hr_2$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_mi_hr_ML_2, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_ML_2$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_mi_hr_ML_2$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_ML_2$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_ML_2$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)


title(outer=TRUE,adj=0.25, main="C.", cex=0.75, col="black",font=2,line=-16.0)


forest(meta_mi_hr_3, atransf = exp, xlab= text(""), slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_mi_hr$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_mi_hr_3$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_3$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_mi_hr_3$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_3$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_3$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_mi_hr_3$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_mi_hr_ML_3, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_ML_3$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_mi_hr_ML_3$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_ML_3$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_ML_3$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)

title(outer=TRUE,adj=0.75, main="D.", cex=0.75, col="black",font=2,line=-16.0)


forest(meta_mi_hr_4, atransf = exp, xlab= text(""), slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_mi_hr$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_mi_hr_4$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_4$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_mi_hr_4$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_4$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_4$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_mi_hr_4$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_mi_hr_ML_4, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_ML_4$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_mi_hr_ML_4$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_ML_4$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_ML_4$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)

title(outer=TRUE,adj=0.25, main="E.", cex=0.75, col="black",font=2,line=-30.0)



forest(meta_mi_hr_5, atransf = exp, xlab= text(""), slab = paste(data_mi_hr$Author, data_mi_hr$Year, sep = ", "), xlim = c(-3,3), ilab = cbind(data_mi_hr$pmid), ilab.xpos=c(-1.4), cex=1.0, ylim=c(-3.5, meta_mi_hr_5$k + 3),  header = c("Author(s) and Year", "Hazard ratios [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_5$QE, digits=2, format="f"), ", p", (metafor:::.pval(meta_mi_hr_5$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_5$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_5$tau2, digits=2, format="f"), ")"), showweights = TRUE, psize = 0.75)
abline(h=0, col = "white")
op <- par(cex=0.75, font=2)
text(c(-1.4), meta_mi_hr_5$k + 2, c("PMID"))
par(cex = 0.75, font = 1)
#text for ML RE model
addpoly(meta_mi_hr_ML_5, atransf = exp, mlab = paste("RE, ML, Wald-type CI", " (Q = ", formatC(meta_mi_hr_ML_5$QE, digits=2, format="f"), ", p ", (metafor:::.pval(meta_mi_hr_ML_5$QEp, digits=2, showeq=TRUE, sep=" ")), ";
I^2 = ", formatC(meta_mi_hr_ML_5$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_mi_hr_ML_5$tau2, digits=2, format="f"), ")"),  cex=1.0, -3)

title(outer=TRUE,adj=0.75, main="F.", cex=0.75, col="black",font=2,line=-30.0)





