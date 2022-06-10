library(metafor)
library(lme4)

#CV death as an outcome. This analysis include only nationwide SEER databases. Regional SEER based cohorts were excluded
data_cvdeath_main_nation_us <- read.csv2("cvdeath_incidence_nationwide_seer.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
#per 1000 person-years of follow-up
data_cvdeath_main_nation_us$PY_FU <- data_cvdeath_main_nation_us$PY_FU/1000
#The median year of enrollment date was calculated
data_cvdeath_main_nation_us$enrollment_median <- round((data_cvdeath_main_nation_us$enrol_end-data_cvdeath_main_nation_us$enrol_begin)/2 + data_cvdeath_main_nation_us$enrol_begin)
#The studies were ordered by their dates of publication
data_cvdeath_main_nation_us <- data_cvdeath_main_nation_us[order(data_cvdeath_main_nation_us$Year),]
#meta-analysis for incidence data. Poisson model for incidence rates
meta_data_cvdeath_main_nation_us <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))
forest(meta_data_cvdeath_main_nation_us, transf = exp)
predict(meta_data_cvdeath_main_nation_us, transf = exp)
#publication bias assessment
ranktest(meta_data_cvdeath_main_nation_us)
funnel_data_cvdeath_main_nation_us <- funnel(meta_data_cvdeath_main_nation_us)


#forest plot - Supplementary Figure S5
forest(meta_data_cvdeath_main_nation_us, transf = exp, xlim = c(-70, 40), slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "), ilab = cbind(data_cvdeath_main_nation_us$country, data_cvdeath_main_nation_us$pmid, data_cvdeath_main_nation_us$n_events, data_cvdeath_main_nation_us$PY_FU*1000), ilab.xpos=c(-50, -35, -20, -10),  header = c("Author(s) and Year", "Incidence rates [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_data_cvdeath_main_nation_us$QE.Wld, digits=2, format="f"), ", p", (metafor:::.pval(meta_data_cvdeath_main_nation_us$QEp.Wld, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_data_cvdeath_main_nation_us$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_data_cvdeath_main_nation_us$tau2, digits=2, format="f"), ")"),  psize = 0.75, cex = 0.67)
abline(h=0, col = "white")
op <- par(cex=0.67, font=2)
text(c(-50, -35, -20, -10), meta_data_cvdeath_main_nation_us$k + 2, c("Country", "PMID", "Events, n", "FU"))



#leave-one-out sensitivity analyses

l1o_data_cvdeath_main_nation_us <- lapply(unique(data_cvdeath_main_nation_us$pmid), function(i)
  rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=subset(data_cvdeath_main_nation_us, pmid != i))
)
sapply(l1o_data_cvdeath_main_nation_us, predict, transf = exp)



#meta-regression analyses. "Mods" argument indicate on an independent variable
metaregr_data_cvdeath_main_nation_us_Year <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~Year, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_enrol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~enrollment_median, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_median_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~median_FU, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_mean_age <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~mean_age, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_males <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~males, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_postmenopause <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~postmenopause, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_race_white <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_white, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_race_hispanic <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_hispanic, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_race_black <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_black, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_race_asian <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_asian, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_history_ihd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_ihd, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_history_cardio <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_cardio, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_history_diabetes <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_diabetes, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_hypertension <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hypertension, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_smoking <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~smoking, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_alcohol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~alcohol, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_dyslipidemia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~dyslipidemia, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_atrial_fibrillation <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~atrial_fibrillation, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_prior_mi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~prior_mi, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_history_hf <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_hf, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_stroke <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_stroke, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_vhd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~vhd, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_ckd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~ckd, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_hormone_replacement <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hormone_replacement, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_bmi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_overweight <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_stage_0 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_0, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_stage_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_1, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_stage_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_2, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_stage_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_3, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_stage_4 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_4, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_grade_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_1, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_grade_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_2, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_grade_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_3, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_estrogen_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~estrogen_positive, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))


metaregr_data_cvdeath_main_nation_us_progesteron_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~progesteron_positive, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))


metaregr_data_cvdeath_main_nation_us_her2_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~her2_positive, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_triple_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))


metaregr_data_cvdeath_main_nation_us_hormone_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_node_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_negative, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_node_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_positive, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_side_left <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_left, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_side_right <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_right, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_size_less_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_less_2, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_size_2_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_2_5, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_size_more_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_more_5, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_surgery_total <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_total, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_surgery_mastectomia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_mastectomia, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_surgery_breast_conserving <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_breast_conserving, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_treatment_adjuvant <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_adjuvant, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_treatment_any_hormone <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_any_hormone, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_treatment_tamoxifen <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_tamoxifen, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_treatment_other_aromatase_inhibitor <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_other_aromatase_inhibitor, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))


metaregr_data_cvdeath_main_nation_us_treatment_chemotherapy <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_chemotherapy, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_treatment_antracycline <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_antracycline, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_treatment_anti_her <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_anti_her, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_treatment_trastuzumab <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_trastuzumab, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

metaregr_data_cvdeath_main_nation_us_radiation_percent <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~radiation_percent, data=data_cvdeath_main_nation_us, slab = paste(data_cvdeath_main_nation_us$Author, data_cvdeath_main_nation_us$Year, sep = ", "))

#As the study by Wildiers et al provided outlied results, we conducted a sensitivity analysis without this study
data_cvdeath_main_nation_us_without_wildiers <- subset(data_cvdeath_main_nation_us, Author != "Wildiers ")

meta_data_cvdeath_main_nation_us_without_wildiers <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_main_nation_us_without_wildiers, slab = paste(data_cvdeath_main_nation_us_without_wildiers$Author, data_cvdeath_main_nation_us_without_wildiers$Year, sep = ", "))
forest(meta_data_cvdeath_main_nation_us_without_wildiers, transf = exp)
predict(meta_data_cvdeath_main_nation_us_without_wildiers, transf = exp)




#CV death as an outcome. This analysis include only regional SEER databases. Nationalwide SEER based cohorts were excluded.
#The comments are similar to the analysis presented previously.
data_cvdeath_main_region_us <- read.csv2("cvdeath_incidence_regional_seer.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cvdeath_main_region_us$PY_FU <- data_cvdeath_main_region_us$PY_FU/1000
data_cvdeath_main_region_us$enrollment_median <- round((data_cvdeath_main_region_us$enrol_end-data_cvdeath_main_region_us$enrol_begin)/2 + data_cvdeath_main_region_us$enrol_begin)
data_cvdeath_main_region_us <- data_cvdeath_main_region_us[order(data_cvdeath_main_region_us$Year),]
meta_data_cvdeath_main_region_us <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
forest(meta_data_cvdeath_main_region_us, transf = exp)
predict(meta_data_cvdeath_main_region_us, transf = exp)



#Forest plot - Figure 5
forest(meta_data_cvdeath_main_region_us, transf = exp, xlim = c(-70, 40), slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "), ilab = cbind(data_cvdeath_main_region_us$country, data_cvdeath_main_region_us$pmid, data_cvdeath_main_region_us$n_events, data_cvdeath_main_region_us$PY_FU*1000), ilab.xpos=c(-50, -35, -20, -10),  header = c("Author(s) and Year", "Incidence rates [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_data_cvdeath_main_region_us$QE.Wld, digits=2, format="f"), ", p", (metafor:::.pval(meta_data_cvdeath_main_region_us$QEp.Wld, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_data_cvdeath_main_region_us$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_data_cvdeath_main_region_us$tau2, digits=2, format="f"), ")"),  psize = 0.75, cex = 0.67)
abline(h=0, col = "white")
op <- par(cex=0.67, font=2)
text(c(-50, -35, -20, -10), meta_data_cvdeath_main_region_us$k + 2, c("Country", "PMID", "Events, n", "FU"))





l1o_data_cvdeath_main_region_us <- lapply(unique(data_cvdeath_main_region_us$pmid), function(i)
  rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=subset(data_cvdeath_main_region_us, pmid != i))
)
sapply(l1o_data_cvdeath_main_region_us, predict, transf = exp)








#subgroup analyses
meta_data_cvdeath_main_region_us_asia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (Asia == "1"), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
meta_data_cvdeath_main_region_us_nonasia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (Asia == "0"), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
metaregr_data_cvdeath_main_region_us_asia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(Asia), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

meta_data_cvdeath_main_region_us_rct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (RCT == "1"), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
meta_data_cvdeath_main_region_us_nonrct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (RCT == "0"), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
metaregr_data_cvdeath_main_region_us_rct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(RCT), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))


#subgroup analyses for DCIS vs invasive BC
data_cvdeath_main_region_us$stage_category <- ifelse(data_cvdeath_main_region_us$stage_0 == "0", "invasive", ifelse(data_cvdeath_main_region_us$stage_0 == "100", "DCIS", ifelse(data_cvdeath_main_region_us$stage_0 < 100 & data_cvdeath_main_region_us$stage_0 > 0, "mixed", NA)))
meta_data_cvdeath_main_region_us_invasive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "invasive"), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
meta_data_cvdeath_main_region_us_DCIS <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "DCIS"), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
meta_data_cvdeath_main_region_us_mixed <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "mixed"), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
#test for subgroup difference
metaregr_data_cvdeath_main_region_us_stage <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(stage_category), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

#subgroup analyses for older vs young BC
meta_data_cvdeath_main_region_us_older <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (older == "1"), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
meta_data_cvdeath_main_region_us_young <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (older == "0"), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
#test for subgroup difference
metaregr_data_cvdeath_main_region_us_older <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(older), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

#subgroup analyses for person-years of follow-up <10000 vs >10000
data_cvdeath_main_region_us$PY_FU_category <- ifelse(data_cvdeath_main_region_us$PY_FU >= "10000", "long_FU", ifelse(data_cvdeath_main_region_us$PY_FU < "10000", "short_FU", NA))
meta_data_cvdeath_main_region_us_short_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (PY_FU_category == "short_FU"), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
meta_data_cvdeath_main_region_us_long_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (PY_FU_category == "long_FU"), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
#test for subgroup difference
metaregr_data_cvdeath_main_region_us_PY_FU_category <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(PY_FU_category), data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))






metaregr_data_cvdeath_main_region_us_Year <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~Year, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_enrol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~enrollment_median, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_median_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~median_FU, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_mean_age <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~mean_age, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_males <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~males, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_postmenopause <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~postmenopause, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_race_white <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_white, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_race_hispanic <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_hispanic, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_race_black <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_black, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_race_asian <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_asian, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_history_ihd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_ihd, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_history_cardio <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_cardio, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_history_diabetes <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_diabetes, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_hypertension <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hypertension, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_smoking <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~smoking, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_alcohol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~alcohol, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_dyslipidemia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~dyslipidemia, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_atrial_fibrillation <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~atrial_fibrillation, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_prior_mi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~prior_mi, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_history_hf <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_hf, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_stroke <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_stroke, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

#metaregr_data_cvdeath_main_region_us_pad <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~pad, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))
#Number of studies that reported PAD is too small

metaregr_data_cvdeath_main_region_us_vhd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~vhd, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_ckd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~ckd, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_hormone_replacement <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hormone_replacement, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_bmi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_overweight <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_stage_0 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_0, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_stage_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_1, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_stage_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_2, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_stage_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_3, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_stage_4 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_4, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_grade_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_1, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_grade_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_2, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_grade_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_3, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_estrogen_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~estrogen_positive, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))


metaregr_data_cvdeath_main_region_us_progesteron_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~progesteron_positive, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))


metaregr_data_cvdeath_main_region_us_her2_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~her2_positive, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_triple_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))


metaregr_data_cvdeath_main_region_us_hormone_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_node_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_negative, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_node_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_positive, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_side_left <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_left, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_side_right <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_right, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_size_less_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_less_2, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_size_2_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_2_5, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_size_more_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_more_5, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_surgery_total <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_total, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_surgery_mastectomia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_mastectomia, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_surgery_breast_conserving <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_breast_conserving, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_treatment_adjuvant <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_adjuvant, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_treatment_any_hormone <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_any_hormone, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_treatment_tamoxifen <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_tamoxifen, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_treatment_other_aromatase_inhibitor <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_other_aromatase_inhibitor, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))


metaregr_data_cvdeath_main_region_us_treatment_chemotherapy <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_chemotherapy, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_treatment_antracycline <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_antracycline, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_treatment_anti_her <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_anti_her, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_treatment_trastuzumab <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_trastuzumab, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))

metaregr_data_cvdeath_main_region_us_radiation_percent <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~radiation_percent, data=data_cvdeath_main_region_us, slab = paste(data_cvdeath_main_region_us$Author, data_cvdeath_main_region_us$Year, sep = ", "))



#The next analyses are sensitivity analyses for CV death. The studies that were conducted on the same cohort of patients were included one after another.
#The name of the first author and PubMed number were incorporated in the names of the R objects 

data_cvdeath_sensitivity_Boekel <- read.csv2("cvdeath_incidence_sensitivity_Boekel_25128694_27026313.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cvdeath_sensitivity_Boekel$PY_FU <- data_cvdeath_sensitivity_Boekel$PY_FU/1000


data_cvdeath_sensitivity_Boekel <- data_cvdeath_sensitivity_Boekel[order(data_cvdeath_sensitivity_Boekel$Year),]


meta_data_cvdeath_sensitivity_Boekel <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_sensitivity_Boekel, slab = paste(data_cvdeath_sensitivity_Boekel$Author, data_cvdeath_sensitivity_Boekel$Year, sep = ", "))
forest(meta_data_cvdeath_sensitivity_Boekel, transf = exp)
predict(meta_data_cvdeath_sensitivity_Boekel, transf = exp)


data_cvdeath_sensitivity_Deen <- read.csv2("cvdeath_incidence_sensitivity_Deen_30121599.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cvdeath_sensitivity_Deen$PY_FU <- data_cvdeath_sensitivity_Deen$PY_FU/1000


data_cvdeath_sensitivity_Deen <- data_cvdeath_sensitivity_Deen[order(data_cvdeath_sensitivity_Deen$Year),]


meta_data_cvdeath_sensitivity_Deen <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_sensitivity_Deen, slab = paste(data_cvdeath_sensitivity_Deen$Author, data_cvdeath_sensitivity_Deen$Year, sep = ", "))
forest(meta_data_cvdeath_sensitivity_Deen, transf = exp)
predict(meta_data_cvdeath_sensitivity_Deen, transf = exp)

data_cvdeath_sensitivity_Elshof <- read.csv2("cvdeath_incidence_sensitivity_Elshof_28375855.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cvdeath_sensitivity_Elshof$PY_FU <- data_cvdeath_sensitivity_Elshof$PY_FU/1000


data_cvdeath_sensitivity_Elshof <- data_cvdeath_sensitivity_Elshof[order(data_cvdeath_sensitivity_Elshof$Year),]


meta_data_cvdeath_sensitivity_Elshof <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_sensitivity_Elshof, slab = paste(data_cvdeath_sensitivity_Elshof$Author, data_cvdeath_sensitivity_Elshof$Year, sep = ", "))
forest(meta_data_cvdeath_sensitivity_Elshof, transf = exp)
predict(meta_data_cvdeath_sensitivity_Elshof, transf = exp)

data_cvdeath_sensitivity_Henson <- read.csv2("cvdeath_incidence_sensitivity_Henson_23257897.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cvdeath_sensitivity_Henson$PY_FU <- data_cvdeath_sensitivity_Henson$PY_FU/1000


data_cvdeath_sensitivity_Henson <- data_cvdeath_sensitivity_Henson[order(data_cvdeath_sensitivity_Henson$Year),]


meta_data_cvdeath_sensitivity_Henson <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_sensitivity_Henson, slab = paste(data_cvdeath_sensitivity_Henson$Author, data_cvdeath_sensitivity_Henson$Year, sep = ", "))
forest(meta_data_cvdeath_sensitivity_Henson, transf = exp)
predict(meta_data_cvdeath_sensitivity_Henson, transf = exp)

data_cvdeath_sensitivity_Khosrow <- read.csv2("cvdeath_incidence_sensitivity_Khosrow_Khavar_32065766.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cvdeath_sensitivity_Khosrow$PY_FU <- data_cvdeath_sensitivity_Khosrow$PY_FU/1000


data_cvdeath_sensitivity_Khosrow <- data_cvdeath_sensitivity_Khosrow[order(data_cvdeath_sensitivity_Khosrow$Year),]


meta_data_cvdeath_sensitivity_Khosrow <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_sensitivity_Khosrow, slab = paste(data_cvdeath_sensitivity_Khosrow$Author, data_cvdeath_sensitivity_Khosrow$Year, sep = ", "))
forest(meta_data_cvdeath_sensitivity_Khosrow, transf = exp)
predict(meta_data_cvdeath_sensitivity_Khosrow, transf = exp)


data_cvdeath_sensitivity_Leung <- read.csv2("cvdeath_incidence_sensitivity_Leung_26359709.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cvdeath_sensitivity_Leung$PY_FU <- data_cvdeath_sensitivity_Leung$PY_FU/1000


data_cvdeath_sensitivity_Leung <- data_cvdeath_sensitivity_Leung[order(data_cvdeath_sensitivity_Leung$Year),]


meta_data_cvdeath_sensitivity_Leung <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_sensitivity_Leung, slab = paste(data_cvdeath_sensitivity_Leung$Author, data_cvdeath_sensitivity_Leung$Year, sep = ", "))
forest(meta_data_cvdeath_sensitivity_Leung, transf = exp)
predict(meta_data_cvdeath_sensitivity_Leung, transf = exp)


data_cvdeath_sensitivity_Lim <- read.csv2("cvdeath_incidence_sensitivity_Lim_33531527.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cvdeath_sensitivity_Lim$PY_FU <- data_cvdeath_sensitivity_Lim$PY_FU/1000


data_cvdeath_sensitivity_Lim <- data_cvdeath_sensitivity_Lim[order(data_cvdeath_sensitivity_Lim$Year),]


meta_data_cvdeath_sensitivity_Lim <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_sensitivity_Lim, slab = paste(data_cvdeath_sensitivity_Leung$Author, data_cvdeath_sensitivity_Lim$Year, sep = ", "))
forest(meta_data_cvdeath_sensitivity_Lim, transf = exp)
predict(meta_data_cvdeath_sensitivity_Lim, transf = exp)


data_cvdeath_sensitivity_Rushton <- read.csv2("cvdeath_incidence_sensitivity_Rushton_32343801.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cvdeath_sensitivity_Rushton$PY_FU <- data_cvdeath_sensitivity_Rushton$PY_FU/1000


data_cvdeath_sensitivity_Rushton <- data_cvdeath_sensitivity_Rushton[order(data_cvdeath_sensitivity_Rushton$Year),]


meta_data_cvdeath_sensitivity_Rushton <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_sensitivity_Rushton, slab = paste(data_cvdeath_sensitivity_Rushton$Author, data_cvdeath_sensitivity_Rushton$Year, sep = ", "))
forest(meta_data_cvdeath_sensitivity_Rushton, transf = exp)
predict(meta_data_cvdeath_sensitivity_Rushton, transf = exp)


data_cvdeath_sensitivity_Ye <- read.csv2("cvdeath_incidence_sensitivity_Ye_25223278.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cvdeath_sensitivity_Ye$PY_FU <- data_cvdeath_sensitivity_Ye$PY_FU/1000


data_cvdeath_sensitivity_Ye <- data_cvdeath_sensitivity_Ye[order(data_cvdeath_sensitivity_Ye$Year),]


meta_data_cvdeath_sensitivity_Ye <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cvdeath_sensitivity_Ye, slab = paste(data_cvdeath_sensitivity_Ye$Author, data_cvdeath_sensitivity_Ye$Year, sep = ", "))
forest(meta_data_cvdeath_sensitivity_Ye, transf = exp)
predict(meta_data_cvdeath_sensitivity_Ye, transf = exp)



#Outcome - heart failure. This analysis include only nationwide SEER databases. Regional SEER based cohorts were excluded
data_hf_main_nation_us <- read.csv2("hf_incidence_main_national_seer.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_main_nation_us$PY_FU <- data_hf_main_nation_us$PY_FU/1000

data_hf_main_nation_us$enrollment_median <- round((data_hf_main_nation_us$enrol_end-data_hf_main_nation_us$enrol_begin)/2 + data_hf_main_nation_us$enrol_begin)
data_hf_main_nation_us <- data_hf_main_nation_us[order(data_hf_main_nation_us$Year),]

meta_data_hf_main_nation_us <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))
forest(meta_data_hf_main_nation_us, transf = exp)
predict(meta_data_hf_main_nation_us, transf = exp)
ranktest(meta_data_hf_main_nation_us)
funnel_data_hf_main_nation_us <- funnel(meta_data_hf_main_nation_us)


#FOREST PLOT - supplementary figure S8
forest(meta_data_hf_main_nation_us, transf = exp, xlim = c(-70, 100), slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "), ilab = cbind(data_hf_main_nation_us$country, data_hf_main_nation_us$pmid, data_hf_main_nation_us$n_events, data_hf_main_nation_us$PY_FU*1000), ilab.xpos=c(-50, -35, -20, -10),  header = c("Author(s) and Year", "Incidence rates [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_data_hf_main_nation_us$QE.Wld, digits=2, format="f"), ", p", (metafor:::.pval(meta_data_hf_main_nation_us$QEp.Wld, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_data_hf_main_nation_us$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_data_hf_main_nation_us$tau2, digits=2, format="f"), ")"),  psize = 0.75, cex = 0.6)
abline(h=0, col = "white")
op <- par(cex=0.6, font=2)
text(c(-50, -35, -20, -10), meta_data_hf_main_nation_us$k + 2, c("Country", "PMID", "Events, n", "FU"))




l1o_data_hf_main_nation_us <- lapply(unique(data_hf_main_nation_us$pmid), function(i)
  rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=subset(data_hf_main_nation_us, pmid != i))
)
sapply(l1o_data_hf_main_nation_us, predict, transf = exp)

















metaregr_data_hf_main_nation_us_Year <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~Year, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_enrol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~enrollment_median, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_median_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~median_FU, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_mean_age <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~mean_age, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_males <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~males, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_postmenopause <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~postmenopause, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_race_white <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_white, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_race_hispanic <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_hispanic, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_race_black <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_black, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_race_asian <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_asian, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_history_ihd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_ihd, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_history_cardio <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_cardio, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_history_diabetes <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_diabetes, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_hypertension <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hypertension, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_smoking <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~smoking, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_alcohol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~alcohol, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_dyslipidemia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~dyslipidemia, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_atrial_fibrillation <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~atrial_fibrillation, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_prior_mi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~prior_mi, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_history_hf <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_hf, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_stroke <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_stroke, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_pad <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~pad, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_vhd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~vhd, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_ckd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~ckd, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

#metaregr_data_hf_main_nation_us_hormone_replacement <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hormone_replacement, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))
#N of studies is too small
metaregr_data_hf_main_nation_us_bmi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_overweight <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_stage_0 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_0, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_stage_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_1, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_stage_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_2, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_stage_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_3, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_stage_4 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_4, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_grade_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_1, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_grade_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_2, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_grade_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_3, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_estrogen_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~estrogen_positive, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))


metaregr_data_hf_main_nation_us_progesteron_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~progesteron_positive, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))


metaregr_data_hf_main_nation_us_her2_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~her2_positive, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_triple_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))


metaregr_data_hf_main_nation_us_hormone_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_node_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_negative, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_node_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_positive, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_side_left <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_left, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_side_right <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_right, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_size_less_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_less_2, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_size_2_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_2_5, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_size_more_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_more_5, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_surgery_total <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_total, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_surgery_mastectomia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_mastectomia, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_surgery_breast_conserving <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_breast_conserving, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_treatment_adjuvant <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_adjuvant, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_treatment_any_hormone <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_any_hormone, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_treatment_tamoxifen <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_tamoxifen, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_treatment_other_aromatase_inhibitor <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_other_aromatase_inhibitor, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))


metaregr_data_hf_main_nation_us_treatment_chemotherapy <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_chemotherapy, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_treatment_antracycline <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_antracycline, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_treatment_anti_her <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_anti_her, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_treatment_trastuzumab <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_trastuzumab, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))

metaregr_data_hf_main_nation_us_radiation_percent <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~radiation_percent, data=data_hf_main_nation_us, slab = paste(data_hf_main_nation_us$Author, data_hf_main_nation_us$Year, sep = ", "))


data_hf_main_nation_us_without_wildiers <- subset(data_hf_main_nation_us, Author != "Wildiers ")

meta_data_hf_main_nation_us_without_wildiers <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_main_nation_us_without_wildiers, slab = paste(data_hf_main_nation_us_without_wildiers$Author, data_hf_main_nation_us_without_wildiers$Year, sep = ", "))
forest(meta_data_hf_main_nation_us_without_wildiers, transf = exp)
predict(meta_data_hf_main_nation_us_without_wildiers, transf = exp)











#Outcome - heart failure. This analysis include only regional SEER databases. Nationwide SEER based cohorts were excluded
data_hf_main_region_us <- read.csv2("hf_incidence_main_regional_seer.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_main_region_us$PY_FU <- data_hf_main_region_us$PY_FU/1000


data_hf_main_region_us$enrollment_median <- round((data_hf_main_region_us$enrol_end-data_hf_main_region_us$enrol_begin)/2 + data_hf_main_region_us$enrol_begin)
data_hf_main_region_us <- data_hf_main_region_us[order(data_hf_main_region_us$Year),]


meta_data_hf_main_region_us <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
forest(meta_data_hf_main_region_us, transf = exp)
predict(meta_data_hf_main_region_us, transf = exp)
ranktest(meta_data_hf_main_region_us)
funnel_data_hf_main_region_us <- funnel(meta_data_hf_main_region_us)


# forest plot - Supplementary figure S7
forest(meta_data_hf_main_region_us, transf = exp, xlim = c(-70, 100), slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "), ilab = cbind(data_hf_main_region_us$country, data_hf_main_region_us$pmid, data_hf_main_region_us$n_events, data_hf_main_region_us$PY_FU*1000), ilab.xpos=c(-50, -35, -20, -10),  header = c("Author(s) and Year", "Incidence rates [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_data_hf_main_region_us$QE.Wld, digits=2, format="f"), ", p", (metafor:::.pval(meta_data_hf_main_region_us$QEp.Wld, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_data_hf_main_region_us$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_data_hf_main_region_us$tau2, digits=2, format="f"), ")"),  psize = 0.75, cex = 0.6)
abline(h=0, col = "white")
op <- par(cex=0.6, font=2)
text(c(-50, -35, -20, -10), meta_data_hf_main_region_us$k + 2, c("Country", "PMID", "Events, n", "FU"))






l1o_data_hf_main_region_us <- lapply(unique(data_hf_main_region_us$pmid), function(i)
  rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=subset(data_hf_main_region_us, pmid != i))
)
sapply(l1o_data_hf_main_region_us, predict, transf = exp)


#subgroup analyses
meta_data_hf_main_region_us_asia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (Asia == "1"), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
meta_data_hf_main_region_us_nonasia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (Asia == "0"), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
metaregr_data_hf_main_region_us_asia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(Asia), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

meta_data_hf_main_region_us_rct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (RCT == "1"), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
meta_data_hf_main_region_us_nonrct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (RCT == "0"), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
metaregr_data_hf_main_region_us_rct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(RCT), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))


#subgroup analyses for DCIS vs invasive BC
#data_hf_main_region_us$stage_category <- ifelse(data_hf_main_region_us$stage_0 == "0", "invasive", ifelse(data_hf_main_region_us$stage_0 == "100", "DCIS", ifelse(data_hf_main_region_us$stage_0 < 100 & data_hf_main_region_us$stage_0 > 0, "mixed", NA)))
#meta_data_hf_main_region_us_invasive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "invasive"), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
#meta_data_hf_main_region_us_DCIS <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "DCIS"), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
#meta_data_hf_main_region_us_mixed <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "mixed"), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
#test for subgroup difference
#metaregr_data_hf_main_region_us_stage <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(stage_category), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

##Error in rma.glmm(measure = "IRLN", xi = n_events, ti = PY_FU, subset = (stage_category ==  : 
## Number of parameters to be estimated is larger than the number of observations.


#subgroup analyses for older vs young BC
meta_data_hf_main_region_us_older <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (older == "1"), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
meta_data_hf_main_region_us_young <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (older == "0"), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
#test for subgroup difference
metaregr_data_hf_main_region_us_older <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(older), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

#subgroup analyses for person-years of follow-up <10000 vs >10000
data_hf_main_region_us$PY_FU_category <- ifelse(data_hf_main_region_us$PY_FU >= "10000", "long_FU", ifelse(data_hf_main_region_us$PY_FU < "10000", "short_FU", NA))
meta_data_hf_main_region_us_short_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (PY_FU_category == "short_FU"), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
meta_data_hf_main_region_us_long_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (PY_FU_category == "long_FU"), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
#test for subgroup difference
metaregr_data_hf_main_region_us_PY_FU_category <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(PY_FU_category), data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

















#meta-regression analyses
metaregr_data_hf_main_region_us_Year <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~Year, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_enrol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~enrollment_median, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_median_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~median_FU, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_mean_age <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~mean_age, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_males <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~males, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_postmenopause <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~postmenopause, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_race_white <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_white, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_race_hispanic <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_hispanic, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_race_black <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_black, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_race_asian <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_asian, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_history_ihd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_ihd, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_history_cardio <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_cardio, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_history_diabetes <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_diabetes, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_hypertension <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hypertension, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_smoking <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~smoking, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_alcohol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~alcohol, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_dyslipidemia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~dyslipidemia, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_atrial_fibrillation <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~atrial_fibrillation, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_prior_mi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~prior_mi, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_history_hf <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_hf, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_stroke <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_stroke, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_pad <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~pad, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_vhd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~vhd, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_ckd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~ckd, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

#metaregr_data_hf_main_region_us_hormone_replacement <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hormone_replacement, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
#N of studies is too small
metaregr_data_hf_main_region_us_bmi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_overweight <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_stage_0 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_0, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_stage_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_1, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_stage_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_2, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_stage_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_3, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_stage_4 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_4, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_grade_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_1, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_grade_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_2, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_grade_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_3, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_estrogen_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~estrogen_positive, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))


metaregr_data_hf_main_region_us_progesteron_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~progesteron_positive, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))


metaregr_data_hf_main_region_us_her2_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~her2_positive, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_triple_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))


#metaregr_data_hf_main_region_us_hormone_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))
#N of studies is too small
metaregr_data_hf_main_region_us_node_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_negative, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_node_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_positive, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_side_left <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_left, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_side_right <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_right, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_size_less_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_less_2, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_size_2_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_2_5, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_size_more_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_more_5, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_surgery_total <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_total, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_surgery_mastectomia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_mastectomia, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_surgery_breast_conserving <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_breast_conserving, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_treatment_adjuvant <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_adjuvant, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_treatment_any_hormone <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_any_hormone, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_treatment_tamoxifen <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_tamoxifen, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_treatment_other_aromatase_inhibitor <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_other_aromatase_inhibitor, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))


metaregr_data_hf_main_region_us_treatment_chemotherapy <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_chemotherapy, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_treatment_antracycline <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_antracycline, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_treatment_anti_her <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_anti_her, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_treatment_trastuzumab <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_trastuzumab, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

metaregr_data_hf_main_region_us_radiation_percent <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~radiation_percent, data=data_hf_main_region_us, slab = paste(data_hf_main_region_us$Author, data_hf_main_region_us$Year, sep = ", "))

data_hf_main_region_us_without_wildiers <- subset(data_hf_main_region_us, Author != "Wildiers ")

meta_data_hf_main_region_us_without_wildiers <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_main_region_us_without_wildiers, slab = paste(data_hf_main_region_us_without_wildiers$Author, data_hf_main_region_us_without_wildiers$Year, sep = ", "))
forest(meta_data_hf_main_region_us_without_wildiers, transf = exp)
predict(meta_data_hf_main_region_us_without_wildiers, transf = exp)

#sensitivity analyses The studies that were conducted on the same cohort of patients were incorporated one after another
data_hf_sensitivity_Wang <- read.csv2("hf_incidence_nation_sensitivity_Wang_24951268.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Wang$PY_FU <- data_hf_sensitivity_Wang$PY_FU/1000


data_hf_sensitivity_Wang <- data_hf_sensitivity_Wang[order(data_hf_sensitivity_Wang$Year),]


meta_data_hf_sensitivity_Wang <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Wang, slab = paste(data_hf_sensitivity_Wang$Author, data_hf_sensitivity_Wang$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Wang, transf = exp)
predict(meta_data_hf_sensitivity_Wang, transf = exp)


data_hf_sensitivity_Cespedes <- read.csv2("hf_incidence_region_sensitivity_Cespedes_28176174.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Cespedes$PY_FU <- data_hf_sensitivity_Cespedes$PY_FU/1000


data_hf_sensitivity_Cespedes <- data_hf_sensitivity_Cespedes[order(data_hf_sensitivity_Cespedes$Year),]


meta_data_hf_sensitivity_Cespedes <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Cespedes, slab = paste(data_hf_sensitivity_Cespedes$Author, data_hf_sensitivity_Cespedes$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Cespedes, transf = exp)
predict(meta_data_hf_sensitivity_Cespedes, transf = exp)





data_hf_sensitivity_Cespedes_2 <- read.csv2("hf_incidence_region_sensitivity_Cespedes_31369302.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Cespedes_2$PY_FU <- data_hf_sensitivity_Cespedes_2$PY_FU/1000


data_hf_sensitivity_Cespedes_2 <- data_hf_sensitivity_Cespedes_2[order(data_hf_sensitivity_Cespedes_2$Year),]


meta_data_hf_sensitivity_Cespedes_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Cespedes_2, slab = paste(data_hf_sensitivity_Cespedes_2$Author, data_hf_sensitivity_Cespedes_2$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Cespedes_2, transf = exp)
predict(meta_data_hf_sensitivity_Cespedes_2, transf = exp)






data_hf_sensitivity_Chou <- read.csv2("hf_incidence_sensitivity_Chou_31996695.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Chou$PY_FU <- data_hf_sensitivity_Chou$PY_FU/1000


data_hf_sensitivity_Chou <- data_hf_sensitivity_Chou[order(data_hf_sensitivity_Chou$Year),]


meta_data_hf_sensitivity_Chou <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Chou, slab = paste(data_hf_sensitivity_Chou$Author, data_hf_sensitivity_Chou$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Chou, transf = exp)
predict(meta_data_hf_sensitivity_Chou, transf = exp)





data_hf_sensitivity_Chung <- read.csv2("hf_incidence_sensitivity_Chung_32771950.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Chung$PY_FU <- data_hf_sensitivity_Chung$PY_FU/1000


data_hf_sensitivity_Chung <- data_hf_sensitivity_Chung[order(data_hf_sensitivity_Chung$Year),]


meta_data_hf_sensitivity_Chung <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Chung, slab = paste(data_hf_sensitivity_Chung$Author, data_hf_sensitivity_Chung$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Chung, transf = exp)
predict(meta_data_hf_sensitivity_Chung, transf = exp)




data_hf_sensitivity_Deen <- read.csv2("hf_incidence_sensitivity_Deen_30121599.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Deen$PY_FU <- data_hf_sensitivity_Deen$PY_FU/1000


data_hf_sensitivity_Deen <- data_hf_sensitivity_Deen[order(data_hf_sensitivity_Deen$Year),]


meta_data_hf_sensitivity_Deen <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Deen, slab = paste(data_hf_sensitivity_Deen$Author, data_hf_sensitivity_Deen$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Deen, transf = exp)
predict(meta_data_hf_sensitivity_Deen, transf = exp)



data_hf_sensitivity_Gong <- read.csv2("hf_incidence_sensitivity_Gong_32236828.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Gong$PY_FU <- data_hf_sensitivity_Gong$PY_FU/1000


data_hf_sensitivity_Gong <- data_hf_sensitivity_Gong[order(data_hf_sensitivity_Gong$Year),]


meta_data_hf_sensitivity_Gong <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Gong, slab = paste(data_hf_sensitivity_Gong$Author, data_hf_sensitivity_Gong$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Gong, transf = exp)
predict(meta_data_hf_sensitivity_Gong, transf = exp)



data_hf_sensitivity_Khan <- read.csv2("hf_incidence_sensitivity_Khan_22048030.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Khan$PY_FU <- data_hf_sensitivity_Khan$PY_FU/1000


data_hf_sensitivity_Khan <- data_hf_sensitivity_Khan[order(data_hf_sensitivity_Khan$Year),]


meta_data_hf_sensitivity_Khan <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Khan, slab = paste(data_hf_sensitivity_Khan$Author, data_hf_sensitivity_Khan$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Khan, transf = exp)
predict(meta_data_hf_sensitivity_Khan, transf = exp)



data_hf_sensitivity_Lee <- read.csv2("hf_incidence_sensitivity_Lee_30690687.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Lee$PY_FU <- data_hf_sensitivity_Lee$PY_FU/1000


data_hf_sensitivity_Lee <- data_hf_sensitivity_Lee[order(data_hf_sensitivity_Lee$Year),]


meta_data_hf_sensitivity_Lee <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Lee, slab = paste(data_hf_sensitivity_Lee$Author, data_hf_sensitivity_Lee$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Lee, transf = exp)
predict(meta_data_hf_sensitivity_Lee, transf = exp)




data_hf_sensitivity_Lee_2 <- read.csv2("hf_incidence_sensitivity_Lee_31454422.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Lee_2$PY_FU <- data_hf_sensitivity_Lee_2$PY_FU/1000


data_hf_sensitivity_Lee_2 <- data_hf_sensitivity_Lee_2[order(data_hf_sensitivity_Lee_2$Year),]


meta_data_hf_sensitivity_Lee_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Lee_2, slab = paste(data_hf_sensitivity_Lee_2$Author, data_hf_sensitivity_Lee_2$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Lee_2, transf = exp)
predict(meta_data_hf_sensitivity_Lee_2, transf = exp)



data_hf_sensitivity_Leung <- read.csv2("hf_incidence_sensitivity_Leung_26359709.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Leung$PY_FU <- data_hf_sensitivity_Leung$PY_FU/1000


data_hf_sensitivity_Leung <- data_hf_sensitivity_Leung[order(data_hf_sensitivity_Leung$Year),]


meta_data_hf_sensitivity_Leung <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Leung, slab = paste(data_hf_sensitivity_Leung$Author, data_hf_sensitivity_Leung$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Leung, transf = exp)
predict(meta_data_hf_sensitivity_Leung, transf = exp)


data_hf_sensitivity_Thavendiranathan <- read.csv2("hf_incidence_sensitivity_Thavendiranathan_27091709.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Thavendiranathan$PY_FU <- data_hf_sensitivity_Thavendiranathan$PY_FU/1000


data_hf_sensitivity_Thavendiranathan <- data_hf_sensitivity_Thavendiranathan[order(data_hf_sensitivity_Thavendiranathan$Year),]


meta_data_hf_sensitivity_Thavendiranathan <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Thavendiranathan, slab = paste(data_hf_sensitivity_Thavendiranathan$Author, data_hf_sensitivity_Thavendiranathan$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Thavendiranathan, transf = exp)
predict(meta_data_hf_sensitivity_Thavendiranathan, transf = exp)



data_hf_sensitivity_Torres <- read.csv2("hf_incidence_sensitivity_Torres_26026467.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Torres$PY_FU <- data_hf_sensitivity_Torres$PY_FU/1000


data_hf_sensitivity_Torres <- data_hf_sensitivity_Torres[order(data_hf_sensitivity_Torres$Year),]


meta_data_hf_sensitivity_Torres <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Torres, slab = paste(data_hf_sensitivity_Torres$Author, data_hf_sensitivity_Torres$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Torres, transf = exp)
predict(meta_data_hf_sensitivity_Torres, transf = exp)




data_hf_sensitivity_Yang <- read.csv2("hf_incidence_sensitivity_Yang_35293856.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_hf_sensitivity_Yang$PY_FU <- data_hf_sensitivity_Yang$PY_FU/1000


data_hf_sensitivity_Yang <- data_hf_sensitivity_Yang[order(data_hf_sensitivity_Yang$Year),]


meta_data_hf_sensitivity_Yang <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_hf_sensitivity_Yang, slab = paste(data_hf_sensitivity_Yang$Author, data_hf_sensitivity_Yang$Year, sep = ", "))
forest(meta_data_hf_sensitivity_Yang, transf = exp)
predict(meta_data_hf_sensitivity_Yang, transf = exp)







#Outcome - coronary artery disease
data_cad_main_nation_us <- read.csv2("cad_incidence_main.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cad_main_nation_us$PY_FU <- data_cad_main_nation_us$PY_FU/1000


data_cad_main_nation_us$enrollment_median <- round((data_cad_main_nation_us$enrol_end-data_cad_main_nation_us$enrol_begin)/2 + data_cad_main_nation_us$enrol_begin)
data_cad_main_nation_us <- data_cad_main_nation_us[order(data_cad_main_nation_us$Year),]


meta_data_cad_main_nation_us <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
forest(meta_data_cad_main_nation_us, transf = exp)
predict(meta_data_cad_main_nation_us, transf = exp)
ranktest(meta_data_cad_main_nation_us)
funnel_data_cad_main_nation_us <- funnel(meta_data_cad_main_nation_us)



#forest plot - Supplementary Figure S9
forest(meta_data_cad_main_nation_us, transf = exp, xlim = c(-70, 40), slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "), ilab = cbind(data_cad_main_nation_us$country, data_cad_main_nation_us$pmid, data_cad_main_nation_us$n_events, data_cad_main_nation_us$PY_FU*1000), ilab.xpos=c(-50, -35, -20, -10),  header = c("Author(s) and Year", "Incidence rates [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_data_cad_main_nation_us$QE.Wld, digits=2, format="f"), ", p", (metafor:::.pval(meta_data_cad_main_nation_us$QEp.Wld, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_data_cad_main_nation_us$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_data_cad_main_nation_us$tau2, digits=2, format="f"), ")"),  psize = 0.75, cex = 0.67)
abline(h=0, col = "white")
op <- par(cex=0.67, font=2)
text(c(-50, -35, -20, -10), meta_data_cad_main_nation_us$k + 2, c("Country", "PMID", "Events, n", "FU"))







l1o_data_cad_main_nation_us <- lapply(unique(data_cad_main_nation_us$pmid), function(i)
  rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=subset(data_cad_main_nation_us, pmid != i))
)
sapply(l1o_data_cad_main_nation_us, predict, transf = exp)



#subgroup analyses
meta_data_cad_main_nation_us_asia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (Asia == "1"), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
meta_data_cad_main_nation_us_nonasia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (Asia == "0"), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
metaregr_data_cad_main_nation_us_asia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(Asia), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

meta_data_cad_main_nation_us_rct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (RCT == "1"), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
meta_data_cad_main_nation_us_nonrct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (RCT == "0"), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
metaregr_data_cad_main_nation_us_rct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(RCT), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))


#subgroup analyses for DCIS vs invasive BC
data_cad_main_nation_us$stage_category <- ifelse(data_cad_main_nation_us$stage_0 == "0", "invasive", ifelse(data_cad_main_nation_us$stage_0 == "100", "DCIS", ifelse(data_cad_main_nation_us$stage_0 < 100 & data_cad_main_nation_us$stage_0 > 0, "mixed", NA)))
meta_data_cad_main_nation_us_invasive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "invasive"), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
meta_data_cad_main_nation_us_DCIS <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "DCIS"), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
meta_data_cad_main_nation_us_mixed <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "mixed"), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
#test for subgroup difference
metaregr_data_cad_main_nation_us_stage <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(stage_category), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))


#subgroup analyses for older vs young BC
meta_data_cad_main_nation_us_older <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (older == "1"), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
meta_data_cad_main_nation_us_young <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (older == "0"), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
#test for subgroup difference
metaregr_data_cad_main_nation_us_older <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(older), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

#subgroup analyses for person-years of follow-up <10000 vs >10000
data_cad_main_nation_us$PY_FU_category <- ifelse(data_cad_main_nation_us$PY_FU >= "10000", "long_FU", ifelse(data_cad_main_nation_us$PY_FU < "10000", "short_FU", NA))
meta_data_cad_main_nation_us_short_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (PY_FU_category == "short_FU"), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
meta_data_cad_main_nation_us_long_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (PY_FU_category == "long_FU"), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
#test for subgroup difference
metaregr_data_cad_main_nation_us_PY_FU_category <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(PY_FU_category), data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))












#meta-regression analyses
metaregr_data_cad_main_nation_us_Year <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~Year, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_enrol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~enrollment_median, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_median_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~median_FU, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_mean_age <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~mean_age, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_males <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~males, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_postmenopause <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~postmenopause, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_race_white <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_white, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_race_hispanic <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_hispanic, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_race_black <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_black, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_race_asian <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_asian, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_history_ihd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_ihd, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_history_cardio <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_cardio, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_history_diabetes <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_diabetes, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_hypertension <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hypertension, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_smoking <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~smoking, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_alcohol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~alcohol, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_dyslipidemia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~dyslipidemia, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_atrial_fibrillation <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~atrial_fibrillation, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_prior_mi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~prior_mi, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_history_hf <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_hf, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_stroke <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_stroke, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

#metaregr_data_cad_main_nation_us_pad <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~pad, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
#N of studies with PAD is too small

metaregr_data_cad_main_nation_us_vhd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~vhd, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

#metaregr_data_cad_main_nation_us_ckd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~ckd, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))
#N of studies is too small
metaregr_data_cad_main_nation_us_hormone_replacement <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hormone_replacement, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_bmi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_overweight <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_stage_0 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_0, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_stage_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_1, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_stage_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_2, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_stage_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_3, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_stage_4 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_4, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_grade_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_1, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_grade_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_2, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_grade_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_3, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_estrogen_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~estrogen_positive, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))


metaregr_data_cad_main_nation_us_progesteron_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~progesteron_positive, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))


metaregr_data_cad_main_nation_us_her2_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~her2_positive, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))



#metaregr_data_cad_main_nation_us_hormone_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_node_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_negative, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_node_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_positive, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_side_left <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_left, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_side_right <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_right, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_size_less_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_less_2, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_size_2_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_2_5, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_size_more_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_more_5, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_surgery_total <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_total, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_surgery_mastectomia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_mastectomia, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_surgery_breast_conserving <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_breast_conserving, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

#metaregr_data_cad_main_nation_us_treatment_adjuvant <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_adjuvant, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_treatment_any_hormone <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_any_hormone, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_treatment_tamoxifen <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_tamoxifen, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_treatment_other_aromatase_inhibitor <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_other_aromatase_inhibitor, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))


metaregr_data_cad_main_nation_us_treatment_chemotherapy <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_chemotherapy, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_treatment_antracycline <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_antracycline, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_treatment_anti_her <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_anti_her, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_treatment_trastuzumab <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_trastuzumab, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))

metaregr_data_cad_main_nation_us_radiation_percent <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~radiation_percent, data=data_cad_main_nation_us, slab = paste(data_cad_main_nation_us$Author, data_cad_main_nation_us$Year, sep = ", "))


data_cad_sensitivity_Bradbury <- read.csv2("cad_incidence_sensitivity_Bradbury_15712362.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cad_sensitivity_Bradbury$PY_FU <- data_cad_sensitivity_Bradbury$PY_FU/1000


data_cad_sensitivity_Bradbury <- data_cad_sensitivity_Bradbury[order(data_cad_sensitivity_Bradbury$Year),]


meta_data_cad_sensitivity_Bradbury <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cad_sensitivity_Bradbury, slab = paste(data_cad_sensitivity_Bradbury$Author, data_cad_sensitivity_Bradbury$Year, sep = ", "))
forest(meta_data_cad_sensitivity_Bradbury, transf = exp)
predict(meta_data_cad_sensitivity_Bradbury, transf = exp)




data_cad_sensitivity_Khan <- read.csv2("cad_incidence_sensitivity_Khan_22048030.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cad_sensitivity_Khan$PY_FU <- data_cad_sensitivity_Khan$PY_FU/1000


data_cad_sensitivity_Khan <- data_cad_sensitivity_Khan[order(data_cad_sensitivity_Khan$Year),]


meta_data_cad_sensitivity_Khan <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cad_sensitivity_Khan, slab = paste(data_cad_sensitivity_Khan$Author, data_cad_sensitivity_Khan$Year, sep = ", "))
forest(meta_data_cad_sensitivity_Khan, transf = exp)
predict(meta_data_cad_sensitivity_Khan, transf = exp)



data_cad_sensitivity_Lee <- read.csv2("cad_incidence_sensitivity_Lee_30690687.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cad_sensitivity_Lee$PY_FU <- data_cad_sensitivity_Lee$PY_FU/1000


data_cad_sensitivity_Lee <- data_cad_sensitivity_Lee[order(data_cad_sensitivity_Lee$Year),]


meta_data_cad_sensitivity_Lee <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cad_sensitivity_Lee, slab = paste(data_cad_sensitivity_Lee$Author, data_cad_sensitivity_Lee$Year, sep = ", "))
forest(meta_data_cad_sensitivity_Lee, transf = exp)
predict(meta_data_cad_sensitivity_Lee, transf = exp)





data_cad_sensitivity_Leung <- read.csv2("cad_incidence_sensitivity_Leung_26359709.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cad_sensitivity_Leung$PY_FU <- data_cad_sensitivity_Leung$PY_FU/1000


data_cad_sensitivity_Leung <- data_cad_sensitivity_Leung[order(data_cad_sensitivity_Leung$Year),]


meta_data_cad_sensitivity_Leung <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cad_sensitivity_Leung, slab = paste(data_cad_sensitivity_Leung$Author, data_cad_sensitivity_Leung$Year, sep = ", "))
forest(meta_data_cad_sensitivity_Leung, transf = exp)
predict(meta_data_cad_sensitivity_Leung, transf = exp)



data_cad_sensitivity_Matthews <- read.csv2("cad_incidence_sensitivity_Matthews_33177117.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_cad_sensitivity_Matthews$PY_FU <- data_cad_sensitivity_Matthews$PY_FU/1000


data_cad_sensitivity_Matthews <- data_cad_sensitivity_Matthews[order(data_cad_sensitivity_Matthews$Year),]


meta_data_cad_sensitivity_Matthews <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_cad_sensitivity_Matthews, slab = paste(data_cad_sensitivity_Matthews$Author, data_cad_sensitivity_Matthews$Year, sep = ", "))
forest(meta_data_cad_sensitivity_Matthews, transf = exp)
predict(meta_data_cad_sensitivity_Matthews, transf = exp)









#Myocardial infarction as an outcome. This analysis include only nationwide SEER databases. Regional SEER based cohorts were excluded
data_mi_main_nation_us <- read.csv2("mi_incidence_national_seer.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_mi_main_nation_us$PY_FU <- data_mi_main_nation_us$PY_FU/1000


data_mi_main_nation_us$enrollment_median <- round((data_mi_main_nation_us$enrol_end-data_mi_main_nation_us$enrol_begin)/2 + data_mi_main_nation_us$enrol_begin)
data_mi_main_nation_us <- data_mi_main_nation_us[order(data_mi_main_nation_us$Year),]


meta_data_mi_main_nation_us <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))
forest(meta_data_mi_main_nation_us, transf = exp)
predict(meta_data_mi_main_nation_us, transf = exp)
ranktest(meta_data_mi_main_nation_us)
funnel_data_mi_main_nation_us <- funnel(meta_data_mi_main_nation_us)



#forest plot - Supplementary Figure S11
forest(meta_data_mi_main_nation_us, transf = exp, xlim = c(-70, 40), slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "), ilab = cbind(data_mi_main_nation_us$country, data_mi_main_nation_us$pmid, data_mi_main_nation_us$n_events, data_mi_main_nation_us$PY_FU*1000), ilab.xpos=c(-50, -35, -20, -10),  header = c("Author(s) and Year", "Incidence rates [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_data_mi_main_nation_us$QE.Wld, digits=2, format="f"), ", p", (metafor:::.pval(meta_data_mi_main_nation_us$QEp.Wld, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_data_mi_main_nation_us$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_data_mi_main_nation_us$tau2, digits=2, format="f"), ")"),  psize = 0.75, cex = 0.9)
abline(h=0, col = "white")
op <- par(cex=0.9, font=2)
text(c(-50, -35, -20, -10), meta_data_mi_main_nation_us$k + 2, c("Country", "PMID", "Events, n", "FU"))






l1o_data_mi_main_nation_us <- lapply(unique(data_mi_main_nation_us$pmid), function(i)
  rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=subset(data_mi_main_nation_us, pmid != i))
)
sapply(l1o_data_mi_main_nation_us, predict, transf = exp)



#meta-regression analyses
metaregr_data_mi_main_nation_us_Year <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~Year, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_enrol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~enrollment_median, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_median_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~median_FU, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_mean_age <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~mean_age, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_males <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~males, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_postmenopause <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~postmenopause, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_race_white <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_white, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_race_hispanic <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_hispanic, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_race_black <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_black, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_race_asian <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_asian, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_history_ihd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_ihd, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_history_cardio <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_cardio, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_history_diabetes <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_diabetes, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_hypertension <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hypertension, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_smoking <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~smoking, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_alcohol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~alcohol, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_dyslipidemia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~dyslipidemia, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_atrial_fibrillation <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~atrial_fibrillation, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_prior_mi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~prior_mi, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_history_hf <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_hf, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_stroke <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_stroke, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_pad <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~pad, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_vhd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~vhd, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_ckd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~ckd, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_hormone_replacement <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hormone_replacement, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_bmi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_overweight <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_stage_0 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_0, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_stage_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_1, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_stage_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_2, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_stage_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_3, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_stage_4 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_4, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_grade_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_1, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_grade_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_2, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_grade_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_3, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_estrogen_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~estrogen_positive, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))


#metaregr_data_mi_main_nation_us_progesteron_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~progesteron_positive, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))


metaregr_data_mi_main_nation_us_her2_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~her2_positive, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_triple_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))


#metaregr_data_mi_main_nation_us_hormone_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_node_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_negative, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_node_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_positive, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_side_left <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_left, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_side_right <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_right, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_size_less_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_less_2, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_size_2_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_2_5, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_size_more_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_more_5, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_surgery_total <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_total, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_surgery_mastectomia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_mastectomia, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_surgery_breast_conserving <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_breast_conserving, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

#metaregr_data_mi_main_nation_us_treatment_adjuvant <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_adjuvant, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_treatment_any_hormone <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_any_hormone, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_treatment_tamoxifen <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_tamoxifen, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_treatment_other_aromatase_inhibitor <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_other_aromatase_inhibitor, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))


metaregr_data_mi_main_nation_us_treatment_chemotherapy <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_chemotherapy, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_treatment_antracycline <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_antracycline, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_treatment_anti_her <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_anti_her, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_treatment_trastuzumab <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_trastuzumab, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))

metaregr_data_mi_main_nation_us_radiation_percent <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~radiation_percent, data=data_mi_main_nation_us, slab = paste(data_mi_main_nation_us$Author, data_mi_main_nation_us$Year, sep = ", "))







#Myocardial infarction as an outcome. This analysis include only nationwide SEER databases. Regional SEER based cohorts were excluded
data_mi_main_region_us <- read.csv2("mi_incidence_regional_seer.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_mi_main_region_us$PY_FU <- data_mi_main_region_us$PY_FU/1000

data_mi_main_region_us$enrollment_median <- round((data_mi_main_region_us$enrol_end-data_mi_main_region_us$enrol_begin)/2 + data_mi_main_region_us$enrol_begin)
data_mi_main_region_us <- data_mi_main_region_us[order(data_mi_main_region_us$Year),]

meta_data_mi_main_region_us <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
forest(meta_data_mi_main_region_us, transf = exp)
predict(meta_data_mi_main_region_us, transf = exp)
ranktest(meta_data_mi_main_region_us)



#forest plot - supplementary figure S9
forest(meta_data_mi_main_region_us, transf = exp, xlim = c(-70, 40), slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "), ilab = cbind(data_mi_main_region_us$country, data_mi_main_region_us$pmid, data_mi_main_region_us$n_events, data_mi_main_region_us$PY_FU*1000), ilab.xpos=c(-50, -35, -20, -10),  header = c("Author(s) and Year", "Incidence rates [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_data_mi_main_region_us$QE.Wld, digits=2, format="f"), ", p", (metafor:::.pval(meta_data_mi_main_region_us$QEp.Wld, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_data_mi_main_region_us$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_data_mi_main_region_us$tau2, digits=2, format="f"), ")"),  psize = 0.75, cex = 0.9)
abline(h=0, col = "white")
op <- par(cex=0.9, font=2)
text(c(-50, -35, -20, -10), meta_data_mi_main_region_us$k + 2, c("Country", "PMID", "Events, n", "FU"))









l1o_data_mi_main_region_us <- lapply(unique(data_mi_main_region_us$pmid), function(i)
  rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=subset(data_mi_main_region_us, pmid != i))
)
sapply(l1o_data_mi_main_region_us, predict, transf = exp)



#subgroup analyses
meta_data_mi_main_region_us_asia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (Asia == "1"), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
meta_data_mi_main_region_us_nonasia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (Asia == "0"), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
metaregr_data_mi_main_region_us_asia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(Asia), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

meta_data_mi_main_region_us_rct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (RCT == "1"), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
meta_data_mi_main_region_us_nonrct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (RCT == "0"), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
metaregr_data_mi_main_region_us_rct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(RCT), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))



#subgroup analyses for DCIS vs invasive BC
#data_mi_main_region_us$stage_category <- ifelse(data_mi_main_region_us$stage_0 == "0", "invasive", ifelse(data_mi_main_region_us$stage_0 == "100", "DCIS", ifelse(data_mi_main_region_us$stage_0 < 100 & data_mi_main_region_us$stage_0 > 0, "mixed", NA)))
#meta_data_mi_main_region_us_invasive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "invasive"), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
#meta_data_mi_main_region_us_DCIS <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "DCIS"), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
#meta_data_mi_main_region_us_mixed <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "mixed"), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
#test for subgroup difference
#metaregr_data_mi_main_region_us_stage <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(stage_category), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
#
##Error: Stopped because k = 0 after subsetting.
#
#subgroup analyses for older vs young BC
meta_data_mi_main_region_us_older <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (older == "1"), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
meta_data_mi_main_region_us_young <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (older == "0"), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
#test for subgroup difference
metaregr_data_mi_main_region_us_older <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(older), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

#subgroup analyses for person-years of follow-up <10000 vs >10000
data_mi_main_region_us$PY_FU_category <- ifelse(data_mi_main_region_us$PY_FU >= "10000", "long_FU", ifelse(data_mi_main_region_us$PY_FU < "10000", "short_FU", NA))
meta_data_mi_main_region_us_short_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (PY_FU_category == "short_FU"), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
meta_data_mi_main_region_us_long_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (PY_FU_category == "long_FU"), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
#test for subgroup difference
metaregr_data_mi_main_region_us_PY_FU_category <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(PY_FU_category), data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))
























#meta-regression analyses
metaregr_data_mi_main_region_us_Year <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~Year, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_enrol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~enrollment_median, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_median_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~median_FU, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_mean_age <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~mean_age, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_males <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~males, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_postmenopause <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~postmenopause, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_race_white <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_white, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_race_hispanic <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_hispanic, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_race_black <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_black, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_race_asian <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_asian, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_history_ihd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_ihd, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_history_cardio <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_cardio, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_history_diabetes <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_diabetes, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_hypertension <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hypertension, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_smoking <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~smoking, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_alcohol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~alcohol, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_dyslipidemia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~dyslipidemia, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

#metaregr_data_mi_main_region_us_atrial_fibrillation <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~atrial_fibrillation, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_prior_mi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~prior_mi, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_history_hf <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_hf, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_stroke <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_stroke, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

#metaregr_data_mi_main_region_us_pad <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~pad, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

#metaregr_data_mi_main_region_us_vhd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~vhd, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_ckd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~ckd, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

#metaregr_data_mi_main_region_us_hormone_replacement <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hormone_replacement, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_bmi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_overweight <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_stage_0 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_0, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_stage_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_1, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_stage_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_2, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_stage_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_3, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_stage_4 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_4, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_grade_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_1, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_grade_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_2, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_grade_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_3, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_estrogen_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~estrogen_positive, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))


metaregr_data_mi_main_region_us_progesteron_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~progesteron_positive, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))


metaregr_data_mi_main_region_us_her2_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~her2_positive, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

#metaregr_data_mi_main_region_us_triple_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))


metaregr_data_mi_main_region_us_hormone_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hormone_positive, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_node_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_negative, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_node_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_positive, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_side_left <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_left, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_side_right <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_right, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

#metaregr_data_mi_main_region_us_size_less_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_less_2, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

#metaregr_data_mi_main_region_us_size_2_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_2_5, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

#metaregr_data_mi_main_region_us_size_more_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_more_5, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_surgery_total <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_total, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_surgery_mastectomia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_mastectomia, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_surgery_breast_conserving <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_breast_conserving, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_treatment_adjuvant <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_adjuvant, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_treatment_any_hormone <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_any_hormone, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_treatment_tamoxifen <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_tamoxifen, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

#metaregr_data_mi_main_region_us_treatment_other_aromatase_inhibitor <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_other_aromatase_inhibitor, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))


metaregr_data_mi_main_region_us_treatment_chemotherapy <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_chemotherapy, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_treatment_antracycline <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_antracycline, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_treatment_anti_her <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_anti_her, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_treatment_trastuzumab <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_trastuzumab, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))

metaregr_data_mi_main_region_us_radiation_percent <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~radiation_percent, data=data_mi_main_region_us, slab = paste(data_mi_main_region_us$Author, data_mi_main_region_us$Year, sep = ", "))


data_mi_sensitivity_Boekel <- read.csv2("mi_incidence_sensitivity_boekel_30065254_25128694.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_mi_sensitivity_Boekel$PY_FU <- data_mi_sensitivity_Boekel$PY_FU/1000


data_mi_sensitivity_Boekel <- data_mi_sensitivity_Boekel[order(data_mi_sensitivity_Boekel$Year),]


meta_data_mi_sensitivity_Boekel <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_mi_sensitivity_Boekel, slab = paste(data_mi_sensitivity_Boekel$Author, data_mi_sensitivity_Boekel$Year, sep = ", "))
forest(meta_data_mi_sensitivity_Boekel, transf = exp)
predict(meta_data_mi_sensitivity_Boekel, transf = exp)

data_mi_sensitivity_Matthews <- read.csv2("mi_incidence_sensitivity_Matthews_33177117.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_mi_sensitivity_Matthews$PY_FU <- data_mi_sensitivity_Matthews$PY_FU/1000


data_mi_sensitivity_Matthews <- data_mi_sensitivity_Matthews[order(data_mi_sensitivity_Matthews$Year),]


meta_data_mi_sensitivity_Matthews <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_mi_sensitivity_Matthews, slab = paste(data_mi_sensitivity_Matthews$Author, data_mi_sensitivity_Matthews$Year, sep = ", "))
forest(meta_data_mi_sensitivity_Matthews, transf = exp)
predict(meta_data_mi_sensitivity_Matthews, transf = exp)










#Stroke as an outcome
data_stroke_main <- read.csv2("stroke_incidence_main.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_stroke_main$PY_FU <- data_stroke_main$PY_FU/1000


data_stroke_main$enrollment_median <- round((data_stroke_main$enrol_end-data_stroke_main$enrol_begin)/2 + data_stroke_main$enrol_begin)
data_stroke_main <- data_stroke_main[order(data_stroke_main$Year),]


meta_data_stroke_main <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
forest(meta_data_stroke_main, transf = exp)
predict(meta_data_stroke_main, transf = exp)
ranktest(meta_data_stroke_main)
funnel_data_stroke_main <- funnel(meta_data_stroke_main)



#forest plot - supplementary Figure S12
forest(meta_data_stroke_main, transf = exp, xlim = c(-70, 40), slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "), ilab = cbind(data_stroke_main$country, data_stroke_main$pmid, data_stroke_main$n_events, data_stroke_main$PY_FU*1000), ilab.xpos=c(-50, -35, -20, -10),  header = c("Author(s) and Year", "Incidence rates [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_data_stroke_main$QE.Wld, digits=2, format="f"), ", p", (metafor:::.pval(meta_data_stroke_main$QEp.Wld, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_data_stroke_main$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_data_stroke_main$tau2, digits=2, format="f"), ")"),  psize = 0.75, cex = 0.9)
abline(h=0, col = "white")
op <- par(cex=0.9, font=2)
text(c(-50, -35, -20, -10), meta_data_stroke_main$k + 2, c("Country", "PMID", "Events, n", "FU"))








l1o_data_stroke_main <- lapply(unique(data_stroke_main$pmid), function(i)
  rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=subset(data_stroke_main, pmid != i))
)
sapply(l1o_data_stroke_main, predict, transf = exp)


#subgroup analyses
meta_data_stroke_main_asia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (Asia == "1"), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
meta_data_stroke_main_nonasia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (Asia == "0"), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
metaregr_data_stroke_main_asia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(Asia), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

meta_data_stroke_main_rct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (RCT == "1"), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
meta_data_stroke_main_nonrct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (RCT == "0"), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
metaregr_data_stroke_main_rct <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(RCT), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

#
#subgroup analyses for DCIS vs invasive BC
#data_stroke_main$stage_category <- ifelse(data_stroke_main$stage_0 == "0", "invasive", ifelse(data_stroke_main$stage_0 == "100", "DCIS", ifelse(data_stroke_main$stage_0 < 100 & data_stroke_main$stage_0 > 0, "mixed", NA)))
#meta_data_stroke_main_invasive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "invasive"), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
#meta_data_stroke_main_DCIS <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "DCIS"), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
#meta_data_stroke_main_mixed <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (stage_category == "mixed"), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
#test for subgroup difference
#metaregr_data_stroke_main_stage <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(stage_category), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
#
####Error: Stopped because k = 0 after subsetting.
#
#subgroup analyses for older vs young BC
meta_data_stroke_main_older <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (older == "1"), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
meta_data_stroke_main_young <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (older == "0"), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
#test for subgroup difference
metaregr_data_stroke_main_older <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(older), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

#subgroup analyses for person-years of follow-up <10000 vs >10000
data_stroke_main$PY_FU_category <- ifelse(data_stroke_main$PY_FU >= "10000", "long_FU", ifelse(data_stroke_main$PY_FU < "10000", "short_FU", NA))
meta_data_stroke_main_short_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (PY_FU_category == "short_FU"), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
meta_data_stroke_main_long_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, subset = (PY_FU_category == "long_FU"), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
#test for subgroup difference
metaregr_data_stroke_main_PY_FU_category <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~factor(PY_FU_category), data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
















#meta-regression analyses
metaregr_data_stroke_main_Year <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~Year, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_enrol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~enrollment_median, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_median_FU <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~median_FU, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_mean_age <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~mean_age, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_males <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~males, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_postmenopause <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~postmenopause, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_race_white <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_white, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_race_hispanic <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_hispanic, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_race_black <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_black, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_race_asian <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~race_asian, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_history_ihd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_ihd, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_history_cardio <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_cardio, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_history_diabetes <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_diabetes, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_hypertension <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hypertension, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_smoking <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~smoking, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_alcohol <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~alcohol, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_dyslipidemia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~dyslipidemia, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

#metaregr_data_stroke_main_atrial_fibrillation <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~atrial_fibrillation, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_prior_mi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~prior_mi, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_history_hf <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_hf, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_stroke <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~history_stroke, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

#metaregr_data_stroke_main_pad <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~pad, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))
#N of studies is too small

#metaregr_data_stroke_main_ckd <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~ckd, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

#metaregr_data_stroke_main_hormone_replacement <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hormone_replacement, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_bmi <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_overweight <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~bmi, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_stage_0 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_0, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_stage_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_1, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_stage_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_2, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_stage_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_3, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_stage_4 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~stage_4, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_grade_1 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_1, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_grade_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_2, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_grade_3 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~grade_3, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_estrogen_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~estrogen_positive, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))


#metaregr_data_stroke_main_progesteron_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~progesteron_positive, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))


metaregr_data_stroke_main_her2_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~her2_positive, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

#metaregr_data_stroke_main_triple_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~triple_negative, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))


metaregr_data_stroke_main_hormone_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~hormone_positive, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_node_negative <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_negative, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_node_positive <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~node_positive, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_side_left <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_left, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_side_right <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~side_right, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

#metaregr_data_stroke_main_size_less_2 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_less_2, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

#metaregr_data_stroke_main_size_more_5 <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~size_more_5, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_surgery_total <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_total, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_surgery_mastectomia <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_mastectomia, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_surgery_breast_conserving <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~surgery_breast_conserving, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

#metaregr_data_stroke_main_treatment_any_hormone <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_any_hormone, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_treatment_tamoxifen <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_tamoxifen, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_treatment_other_aromatase_inhibitor <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_other_aromatase_inhibitor, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))


metaregr_data_stroke_main_treatment_chemotherapy <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_chemotherapy, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_treatment_antracycline <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_antracycline, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_treatment_anti_her <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_anti_her, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

#metaregr_data_stroke_main_treatment_trastuzumab <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~treatment_trastuzumab, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))

metaregr_data_stroke_main_radiation_percent <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, mods =~radiation_percent, data=data_stroke_main, slab = paste(data_stroke_main$Author, data_stroke_main$Year, sep = ", "))


#sensitivity analyses
data_stroke_sensitivity_Cespedes <- read.csv2("stroke_incidence_sensitivity_Cespedes_31369302.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_stroke_sensitivity_Cespedes$PY_FU <- data_stroke_sensitivity_Cespedes$PY_FU/1000


data_stroke_sensitivity_Cespedes <- data_stroke_sensitivity_Cespedes[order(data_stroke_sensitivity_Cespedes$Year),]


meta_data_stroke_sensitivity_Cespedes <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_stroke_sensitivity_Cespedes, slab = paste(data_stroke_sensitivity_Cespedes$Author, data_stroke_sensitivity_Cespedes$Year, sep = ", "))
forest(meta_data_stroke_sensitivity_Cespedes, transf = exp)
predict(meta_data_stroke_sensitivity_Cespedes, transf = exp)

data_stroke_sensitivity_Kim <- read.csv2("stroke_incidence_sensitivity_Kim_34369199.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_stroke_sensitivity_Kim$PY_FU <- data_stroke_sensitivity_Kim$PY_FU/1000


data_stroke_sensitivity_Kim <- data_stroke_sensitivity_Kim[order(data_stroke_sensitivity_Kim$Year),]


meta_data_stroke_sensitivity_Kim <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_stroke_sensitivity_Kim, slab = paste(data_stroke_sensitivity_Kim$Author, data_stroke_sensitivity_Kim$Year, sep = ", "))
forest(meta_data_stroke_sensitivity_Kim, transf = exp)
predict(meta_data_stroke_sensitivity_Kim, transf = exp)








#Ischemic stroke as an outcome
data_ischemic_stroke_main <- read.csv2("ischemic_stroke_incidence_main.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_ischemic_stroke_main$PY_FU <- data_ischemic_stroke_main$PY_FU/1000


data_ischemic_stroke_main$enrollment_median <- round((data_ischemic_stroke_main$enrol_end-data_ischemic_stroke_main$enrol_begin)/2 + data_ischemic_stroke_main$enrol_begin)
data_ischemic_stroke_main <- data_ischemic_stroke_main[order(data_ischemic_stroke_main$Year),]


meta_data_ischemic_stroke_main <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_ischemic_stroke_main, slab = paste(data_ischemic_stroke_main$Author, data_ischemic_stroke_main$Year, sep = ", "))
forest(meta_data_ischemic_stroke_main, transf = exp)
predict(meta_data_ischemic_stroke_main, transf = exp)
ranktest(meta_data_ischemic_stroke_main)
funnel_data_ischemic_stroke_main <- funnel(meta_data_ischemic_stroke_main)

#forest plot - supplementary figure S13
forest(meta_data_ischemic_stroke_main, transf = exp, xlim = c(-70, 40), slab = paste(data_ischemic_stroke_main$Author, data_ischemic_stroke_main$Year, sep = ", "), ilab = cbind(data_ischemic_stroke_main$country, data_ischemic_stroke_main$pmid, data_ischemic_stroke_main$n_events, data_ischemic_stroke_main$PY_FU*1000), ilab.xpos=c(-50, -35, -20, -10),  header = c("Author(s) and Year", "Incidence rates [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_data_ischemic_stroke_main$QE.Wld, digits=2, format="f"), ", p", (metafor:::.pval(meta_data_ischemic_stroke_main$QEp.Wld, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_data_ischemic_stroke_main$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_data_ischemic_stroke_main$tau2, digits=2, format="f"), ")"),  psize = 0.75, cex = 1.0)
abline(h=0, col = "white")
op <- par(cex=1.0, font=2)
text(c(-50, -35, -20, -10), meta_data_ischemic_stroke_main$k + 2, c("Country", "PMID", "Events, n", "FU"))










l1o_data_ischemic_stroke_main <- lapply(unique(data_ischemic_stroke_main$pmid), function(i)
  rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=subset(data_ischemic_stroke_main, pmid != i))
)
sapply(l1o_data_ischemic_stroke_main, predict, transf = exp)





#atrial fibrillation as an outcome
data_af_main <- read.csv2("af_incidence_main.csv", header = TRUE, sep = ",", dec = ".", stringsAsFactors=FALSE)
data_af_main$PY_FU <- data_af_main$PY_FU/1000


data_af_main$enrollment_median <- round((data_af_main$enrol_end-data_af_main$enrol_begin)/2 + data_af_main$enrol_begin)
data_af_main <- data_af_main[order(data_af_main$Year),]


meta_data_af_main <- rma.glmm(measure="IRLN", xi=n_events, ti=PY_FU, data=data_af_main, slab = paste(data_af_main$Author, data_af_main$Year, sep = ", "))
forest(meta_data_af_main, transf = exp)
predict(meta_data_af_main, transf = exp)



#forest plot - supplementary figure S14
forest(meta_data_af_main, transf = exp, xlim = c(-70, 40), slab = paste(data_af_main$Author, data_af_main$Year, sep = ", "), ilab = cbind(data_af_main$country, data_af_main$pmid, data_af_main$n_events, data_af_main$PY_FU*1000), ilab.xpos=c(-50, -35, -20, -10),  header = c("Author(s) and Year", "Incidence rates [95% CI]"), mlab = paste("RE, REML, Wald-type CI", " (Q = ", formatC(meta_data_af_main$QE.Wld, digits=2, format="f"), ", p", (metafor:::.pval(meta_data_af_main$QEp.Wld, digits=2, showeq=TRUE, sep=" ")), "; I^2 = ", formatC(meta_data_af_main$I2, digits=1, format="f"), "%, tau^2 = ", formatC(meta_data_af_main$tau2, digits=2, format="f"), ")"),  psize = 0.75, cex = 1.0)
abline(h=0, col = "white")
op <- par(cex=1.0, font=2)
text(c(-50, -35, -20, -10), meta_data_af_main$k + 2, c("Country", "PMID", "Events, n", "FU"))





#Funnel plots for IR analyses. Supplementary Figure S6

par(mfrow = c(3,2), mai = c(0.2, 0.8,0.2, 0.4), omi = c(0.2,0.2,0.2,0.2))
funnel(meta_data_cvdeath_main_region_us)
test_data_cvdeath_main_region_us <- ranktest(meta_data_cvdeath_main_region_us, digits = 2)
title(outer=TRUE,adj=0.25, main=paste("A. Cardiovascular death. Rank test p value = ", round(test_data_cvdeath_main_region_us$pval, 2)), cex=0.75, col="black",font=2,line=-1.0)

funnel(meta_data_hf_main_region_us)
test_data_hf_main_region_us <- ranktest(meta_data_hf_main_region_us, digits = 2)
title(outer=TRUE,adj=0.8, main=paste("B. Heart failure. Rank test p value = ", round(test_data_hf_main_region_us$pval, 2)), cex=0.75, col="black",font=2,line=-1.0)

funnel(meta_data_cad_main_nation_us)
test_data_cad_main_nation_us <- ranktest(meta_data_cad_main_nation_us, digits = 2)
title(outer=TRUE,adj=0.25, main=paste("C. Coronary artery disease. Rank test p value = ", round(test_data_cad_main_nation_us$pval, 2)), cex=0.75, col="black",font=2,line=-17.5)

funnel(meta_data_mi_main_region_us)
test_data_mi_main_region_us <- ranktest(meta_data_mi_main_region_us, digits = 2)
title(outer=TRUE,adj=0.85, main=paste("C. Myocardial infarction. Rank test p value = ", round(test_data_mi_main_region_us$pval, 2)), cex=0.75, col="black",font=2,line=-17.5)

funnel(meta_data_stroke_main)
test_data_stroke_main <- ranktest(meta_data_stroke_main, digits = 2)
title(outer=TRUE,adj=0.25, main=paste("C. Stroke. Rank test p value = ", round(test_data_stroke_main$pval, 2)), cex=0.75, col="black",font=2,line=-33.5)


#The results of meta-regression analyses in one table
list_obj <- ls()
names <- as.vector(grep("metaregr*", list_obj, value = TRUE))

n_studies <- sapply(lapply(ls(pattern="metaregr*"), get), function(x) x$k)
beta_for_variable <- sapply(lapply(ls(pattern="metaregr*"), get), function(x) x$beta[2])
p_for_variable <- sapply(lapply(ls(pattern="metaregr*"), get), function(x) x$pval[2])
I2_for_variable <- sapply(lapply(ls(pattern="metaregr*"), get), function(x) x$I2)
metaregr_results <- data.frame (n_studies, beta_for_variable, p_for_variable,I2_for_variable, rownames = names)
metaregr_results







