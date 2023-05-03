library(readxl)
library(lme4)
library(emmeans)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(lm.beta)
library(dplyr)
library(arsenal)
library(knitr)
library(kableExtra)
library(epikit)
library(reshape2)
library(ggpubr)
library(flextable)
library(r2glmm)
library(mice)
library(finalfit)



# select volume data
ROIs <- as.data.frame(read_csv("[path]"))

# Calculate mean or copy unilateral volume data if other hemisphere did not pass quality control
mean_lh_rh <- function (lh, rh){
  if (! is.na(lh)){
    if (! is.na(rh)) {
      (lh + rh ) /2
    }
    else if(is.na(rh)){
      (lh + lh ) /2
    }
  }
  else if (is.na(lh)){
    if (! is.na(rh)){
      (rh + rh ) /2
    }
  } else {NA}
}

ROIs <- ROIs %>%
  filter(MRI_Visit == "1")%>%
  transmute(Pseudonym, eTIV, 
            AD_mask = mapply(mean_lh_rh, AD_mask_lh, AD_mask_rh), 
            Precentral = mapply(mean_lh_rh, precentral_lh, precentral_rh),
            #Amygdala
            Whole_amygdala = mapply(mean_lh_rh, Whole_amygdala_lh, Whole_amygdala_rh),
            Lat_nucl = mapply(mean_lh_rh, Lat_nucl_lh, Lat_nucl_rh),
            Bas_nucl = mapply(mean_lh_rh, Bas_nucl_lh, Bas_nucl_rh),
            Acces_Bas_nucl = mapply(mean_lh_rh, Acces_Bas_nucl_lh, Acces_Bas_nucl_rh),
            AAA = mapply(mean_lh_rh, AAA_lh, AAA_rh),
            Centr_nucl = mapply(mean_lh_rh, Centr_nucl_lh, Centr_nucl_rh),
            Med_nucl = mapply(mean_lh_rh, Med_nucl_lh, Med_nucl_rh),
            Cort_nucl = mapply(mean_lh_rh, Cort_nucl_lh, Cort_nucl_rh),
            CATA = mapply(mean_lh_rh, CATA_lh, CATA_rh),
            paral_nucl = mapply(mean_lh_rh, paral_nucl_lh, paral_nucl_rh),
            
            #Hippocampus
            Whole_Hippocampus = mapply(mean_lh_rh, Whole_Hippocampus_lh, Whole_Hippocampus_rh),
            CA1 = mapply(mean_lh_rh, CA1_lh, CA1_rh),
            CA3 = mapply(mean_lh_rh, CA3_lh, CA3_rh),
            CA4 = mapply(mean_lh_rh, CA4_lh, CA4_rh),
            Subiculum = mapply(mean_lh_rh, Subiculum_lh, Subiculum_rh),
            molec_layer = mapply(mean_lh_rh, molec_layer_lh, molec_layer_rh),
            GC_ML_DG = mapply(mean_lh_rh, GC_ML_DG_lh, GC_ML_DG_rh),
            HATA = mapply(mean_lh_rh, HATA_lh, HATA_rh))
ROIs[ROIs == "NULL"] <- NA

# add diagnosis to volume data
setwd("[path]")
Vol <- as.data.frame(read_excel("[path]"))
ROIs <- Vol %>%
  filter(Visit == "1")%>%
  select('Pseudonym', 'Group_0HC_1SCD_2_MCI_3AD')%>%
  full_join(ROIs, by = "Pseudonym")

# adjust volumes to eTIV  
ROI_list <- names(ROIs[4:23])
ROIs_HC <- ROIs %>%
  filter (Group_0HC_1SCD_2_MCI_3AD=="0")
meanTIV_HC = mean(ROIs_HC$eTIV, na.rm = TRUE)
adj_ROIs <- sapply(ROI_list, function (ROI) {
  ROIs_HC$ROI_HC <- unlist(ROIs_HC[, ROI])
  regr <- lm(ROI_HC ~ eTIV, ROIs_HC)
  slope <- summary(regr)$coefficients[2,1]
  adjVol = unlist(ROIs[, ROI]) - slope * (ROIs$eTIV - meanTIV_HC)
})
adj_ROIs <- data.frame(Pseudonym = ROIs$Pseudonym, as.data.frame(adj_ROIs))

# merge volumes to main data set. Baseline Volumes will be inserted also for follow-up visits (necessary for linear mixed models)
Vol <- Vol %>%
  full_join(adj_ROIs, by = "Pseudonym")

#Check ROIs
#Vol %>% select(28:37)%>%  boxplot()  # Amygdala
#Vol %>% select(38:45)%>%  boxplot()  # Hippocampus
#Vol %>% select(26,27)%>%  boxplot()  # AD-mask, precentral gyrus

##### prepare main data set #####
Vol$Diagnose <- factor(Vol$Group_0HC_1SCD_2_MCI_3AD, labels=c('HC', 'SCD', 'MCI', 'AD'))
Vol$ApoE <- factor(Vol$ApoE, labels=c('22','23','24','33','34','44'))
Vol$APOEe4[Vol$ApoE=="22" | Vol$ApoE=="23" | Vol$ApoE=="33"] <- "0"
Vol$APOEe4[Vol$ApoE=="24" | Vol$ApoE=="34" | Vol$ApoE=="44"] <- "1"
Vol$APOEe4 <- factor(Vol$APOEe4, labels=c('non-carrier','carrier'))
Vol$Sex <- factor(Vol$Sex_0man1woman, labels=c('male','female'))
Vol$Age <- Vol$Age_Years
Vol$Education <- Vol$Education_CEARD_Years
Vol$Visit <- factor(Vol$Visit, level = c("1", "2", "3", "4"))
Diagnose_1 <- Vol %>% filter(Visit == '1') %>% mutate(Diagnose_Baseline= factor(Group_0HC_1SCD_2_MCI_3AD, labels=c('HC', 'SCD', 'MCI', 'AD'))) %>% select(Pseudonym, Diagnose_Baseline)
Vol <- Vol %>%
  full_join(Diagnose_1, by = "Pseudonym")%>%
  filter(is.na(`exclusion reason`), no_MRI_at_Baseline=="0", !Visit == '4')%>%
  select(-`exclusion reason`, -Sex_0man1woman, -Education_CEARD_Years, -Age_Years, -Group_0HC_1SCD_2_MCI_3AD, -`...9`, -`...14`, -no_MRI_at_Baseline)

# Calulate years
Vol <- Vol %>%
  group_by(Pseudonym)%>%
  mutate(year = as.numeric(as.Date(Date, format = "%m/%d/%Y") - 
                             lag(as.Date(Date, format = "%m/%d/%Y"), default = as.Date(Date, format = "%m/%d/%Y")[1])
  ))%>%
  mutate(year = cumsum(year))%>%
  mutate(year = round(year/365,2))%>%
  ungroup()%>%
  filter (!is.na(ID))

# Calculate Memory # tested an alternative replacing DS_back by Logical Memory_delayed but crutial missing values -> longitudinal memory overestimated
sum <- Vol %>%
  filter(Diagnose == 'HC'& Visit  == '1')%>%
  dplyr::summarise(mean_AVLT1to5 = mean(`AVLT_trial1:5`, na.rm = TRUE), sd_AVLT1to5 = sd(`AVLT_trial1:5`, na.rm = TRUE),
            mean_AVLT7 = mean(AVLT_trial7, na.rm = TRUE), sd_AVLT7 = sd(AVLT_trial7, na.rm = TRUE),
            mean_AVLTrecogn = mean(AVLT_recogn, na.rm = TRUE), sd_AVLTrecogn = sd(AVLT_recogn, na.rm = TRUE),
            mean_DS_back = mean(DigitSpan_backwards, na.rm = TRUE), sd_DS_back = sd(DigitSpan_backwards, na.rm = TRUE),
  )

z_score_AVLT1to5 <- (Vol$`AVLT_trial1:5` - sum$mean_AVLT1to5) / sum$sd_AVLT1to5
z_score_AVLT7 <- (Vol$AVLT_trial7 - sum$mean_AVLT7) / sum$sd_AVLT7   
z_score_AVLTrecogn <- (Vol$AVLT_recogn - sum$mean_AVLTrecogn) / sum$sd_AVLTrecogn 
z_score_DS_back <- (Vol$DigitSpan_backwards - sum$mean_DS_back) / sum$sd_DS_back

Vol$Memory_vdRest <- (z_score_AVLT1to5 + z_score_AVLT7 + z_score_AVLTrecogn + z_score_DS_back)/4

# Explore missing data and exclude subjects with missing Memory at Baseline
Vol%>%
  missing_plot()

#Identify who doesn't have Memory Score at V1. All 3 visits have to be excluded
excl <- Vol %>% filter(Vol$Visit == '1' & is.na(Vol$Memory_vdRest)) %>% select(Pseudonym)
excl$Pseudonym

Vol <- Vol %>% filter(!(Pseudonym == excl$Pseudonym))

Vol %>%
  select(Memory_vdRest, Whole_amygdala, Lat_nucl, Bas_nucl, Acces_Bas_nucl, AAA, Centr_nucl,
         Med_nucl, Cort_nucl, CATA, paral_nucl, Whole_Hippocampus, 
         CA1, CA3, CA4, molec_layer, Subiculum, GC_ML_DG, HATA, Precentral)%>%
  is.na() %>%
  colSums()

# Define Data at Visits
Vol_1 <- Vol %>% 
  filter(Visit == '1')
Vol_1 %>% missing_plot()

Vol_2 <- Vol %>% 
  filter(Visit == '2')
Vol_2 %>% missing_plot()

Vol_3 <- Vol %>% 
  filter(Visit == '3')
Vol_3 %>% missing_plot()

# Plot Memory at Baseline
Vol_1 %>%
  ggplot(aes(x = factor(Diagnose, labels=c('HC', 'SCD', 'MCI', 'AD')), y = Memory_vdRest))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(position=position_jitter(0.2))

# std. ß function
std.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  std.beta.lmer <- b*sd.x/sd.y
}

##### Figure A2: Boxplots of all ROI volumes ####
#Calculating breaks
means_HC <- Vol_1 %>%
  filter(Diagnose_Baseline =="HC")%>%
  select(Diagnose_Baseline, Whole_amygdala, Lat_nucl, Bas_nucl, Acces_Bas_nucl, AAA, Centr_nucl,
         Med_nucl, Cort_nucl, CATA, paral_nucl, Whole_Hippocampus, 
         CA1, CA3, CA4, molec_layer, Subiculum, GC_ML_DG, HATA, Precentral)%>%
  group_by(Diagnose_Baseline)%>%
  summarise_at (vars(Whole_amygdala:Precentral),mean, na.rm = TRUE)%>%
  mutate()

breaks <- function(mean){
  c(mean, mean-0.1*mean, mean-0.2*mean, mean-0.3*mean)
}# use as breaks(means_HC$Precentral)

breaks_Vol_unscaled <- lapply(as.list(means_HC[2:20]), breaks)
ROI_levels_Figure_A2 <- c("Whole_amygdala", "Lat_nucl", "Bas_nucl", "Acces_Bas_nucl", "AAA", "Centr_nucl",
  "Med_nucl", "Cort_nucl", "CATA", "paral_nucl", "Whole_Hippocampus", 
  "CA1", "CA3", "CA4", "molec_layer", "Subiculum", "GC_ML_DG", "HATA", "Precentral")
ROIs_names_Figure_A2 <- c("Amygdala", "Lateral Nucleus", "Basal Nucleus", "Acc. Basal Nucleus", "AAA", "Central Nucleus",
                "Medial Nucleus", "Cortical Nucleus", "CATA", "Paralaminar Nucleus","Hippocampus", "CA1", "CA3", 
                "CA4", "Molecular Layer", "Subiculum","GC-ML-DG","HATA", "Precentral Gyrus")


FigureA2 <- Vol_1 %>%
  filter(!is.na(Diagnose_Baseline))%>%
  select(Pseudonym, Diagnose_Baseline, 36:54)%>%
  pivot_longer(cols = 3:21, names_to = "ROI", values_drop_na = TRUE)%>%
  mutate(ROI = factor(ROI, levels = ROI_levels_Figure_A2, labels = ROIs_names_Figure_A2)) %>%
  #group_by(Diagnose_Baseline, ROI)%>%
  #ungroup () %>%
  group_by(Diagnose_Baseline)%>%
  ggplot()+ 
  aes(x= Diagnose_Baseline, y= value, fill = Diagnose_Baseline)+
  geom_boxplot(outlier.shape = NA, alpha = 0.95, lwd = 0.3)+
  geom_jitter(position=position_jitter(0.2), size = 0.5, alpha = 0.5)+
  xlab(NULL)+
  ylab("Volume [mm³]") +
  #geom_text_repel(aes(label = outlier), na.rm = TRUE, max.overlaps = 95, size = 1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+ ####remove background grid
  scale_fill_manual(values=c("gray50", "#ABDDA4","#41B6C4" , "#225EA8"))+
  facet_wrap(~ROI, scales = "free", ncol = 4)+
  theme(legend.position = "none")+ 
  geom_hline(data = as.data.frame(breaks_Vol_unscaled) %>% 
                        filter(row_number()==1) %>% 
                        pivot_longer(cols = everything(), names_to = "ROI", values_to = "mean") %>% 
                        mutate(ROI = factor(ROI, levels = ROI_levels_Figure_A2, labels = ROIs_names_Figure_A2)), 
                      aes(yintercept = mean), color = "dark grey", size = 0.25)+
  geom_hline(data = as.data.frame(breaks_Vol_unscaled) %>% 
               filter(row_number()==2) %>% 
               pivot_longer(cols = everything(), names_to = "ROI", values_to = "ten") %>% 
               mutate(ROI = factor(ROI, levels = ROI_levels_Figure_A2, labels = ROIs_names_Figure_A2)), 
             aes(yintercept = ten), linetype = "longdash", color = "dark grey", size = 0.25)+
  geom_hline(data = as.data.frame(breaks_Vol_unscaled) %>% 
               filter(row_number()==3) %>% 
               pivot_longer(cols = everything(), names_to = "ROI", values_to = "twenty") %>% 
               mutate(ROI = factor(ROI, levels = ROI_levels_Figure_A2, labels = ROIs_names_Figure_A2)), 
             aes(yintercept = twenty), linetype = "dashed", color = "dark grey", size = 0.25)
FigureA2
#ggsave("FigureA2.pdf", height = 11, width = 7)



##### Table 1:  Participant's characteristics at baseline.#####
table_1 <- tableby(Diagnose_Baseline ~ 
                     Sex +
                     anova(Age, "meansd", digits=0) + 
                     anova(Education, "meansd", digits=1) +
                     APOEe4,
                   data = Vol_1,
                   cat.simplify = TRUE, numeric.simplify = TRUE)
#write2word(summary(table_1), "[path]", title = "Participant's characteristics at baseline.")

######## Post-Hoc Demographics ########
ANOVA_Age <- lm(Age ~ Diagnose, Vol_1)
# tidy(ANOVA_Age, conf.int = TRUE) #(library(broom))
em_ANOVA_Age <- emmeans(ANOVA_Age, pairwise ~ Diagnose)
em_ANOVA_Age

ANOVA_Education <- lm(Education ~ Diagnose, Vol_1)
em_ANOVA_Education <- emmeans(ANOVA_Education, pairwise ~ Diagnose)
em_ANOVA_Education

##### Table 2: See below(after Model Memory)#####
##### Table A1: Adjusted mean volumes of Amygdala and Hippocampus ROIs across the study groups. Group differences were tested by ANCOVA adjusting for age, sex and education. The resulting adjusted means [95 % CI] and p-values are reported. #####

ROIs_names <- c("Amygdala", "Lateral Nucleus", "Basal Nucleus", "Accessory basal Nucleus", "Anterior Amygdaloid Area", "Central Nucleus",
                "Medial Nucleus", "Cortical Nucleus", "CATA", "Paralaminar Nucleus","Hippocampus", "CA1", "CA3", 
                "CA4", "Molecular Layer", "Subiculum","GC-ML-DG","HATA", "Precentral Gyrus")
ROIs <- Vol_1 %>%
  select(Whole_amygdala, Lat_nucl, Bas_nucl, Acces_Bas_nucl, AAA, Centr_nucl,
         Med_nucl, Cort_nucl, CATA, paral_nucl, Whole_Hippocampus, 
         CA1, CA3, CA4, molec_layer, Subiculum, GC_ML_DG, HATA, Precentral)
ROIs <-  as.list(ROIs)

my_anovas_vol <- lapply(ROIs, function(roi){
  model <- lm(roi ~ Vol_1$Diagnose + Vol_1$Age + Vol_1$Education + Vol_1$Sex, Vol_1)
  attr(model, 'call')$formula <- model
})
my_anovas_emmeans <- lapply(my_anovas_vol, function(model)(emmeans(model, pairwise ~ "Diagnose")))
my_anovas_emmeans_confint <- lapply(my_anovas_emmeans, confint)

coefficients_tableA1 <- function(anovas_vol, anovas_emmeans){
  c(n = nobs(anovas_vol), 
    emmean_HC = round(as.data.frame(anovas_emmeans$emmeans)[1,2],0),
    emmean_HC_confint_low = round(as.data.frame(anovas_emmeans$emmeans)[1,5],0),
    emmean_HC_confint_up = round(as.data.frame(anovas_emmeans$emmeans)[1,6],0),
    emmean_SCD = round(as.data.frame(anovas_emmeans$emmeans)[2,2],0),
    emmean_SCD_confint_low = round(as.data.frame(anovas_emmeans$emmeans)[2,5],0),
    emmean_SCD_confint_up = round(as.data.frame(anovas_emmeans$emmeans)[2,6],0),
    emmean_MCI = round(as.data.frame(anovas_emmeans$emmeans)[3,2],0),
    emmean_MCI_confint_low = round(as.data.frame(anovas_emmeans$emmeans)[3,5],0),
    emmean_MCI_confint_up = round(as.data.frame(anovas_emmeans$emmeans)[3,6],0),
    emmean_AD = round(as.data.frame(anovas_emmeans$emmeans)[4,2],0),
    emmean_AD_confint_low = round(as.data.frame(anovas_emmeans$emmeans)[4,5],0),
    emmean_AD_confint_up = round(as.data.frame(anovas_emmeans$emmeans)[4,6],0),
    p_value = round((anova(anovas_vol)$`Pr(>F)`)[1],3))
}

m <- mapply (coefficients_tableA1, my_anovas_vol, my_anovas_emmeans)
df <- data.frame(ROIs_names, as.data.frame(t(m)))
df$p_value[df$p_value <  0.001] <- "<0.001"

set_flextable_defaults(big.mark = " ", 
                       font.size = 8, digits = 0)
ft <- flextable(
  data = df, 
  col_keys = c("ROIs_names","n", "dummy_HC", "dummy_SCD", "dummy_MCI", "dummy_AD", "p_value")) %>%
  compose(j = "dummy_HC", value = as_paragraph(as_chunk(emmean_HC), as_bracket(emmean_HC_confint_low, emmean_HC_confint_up, sep = "; ", p = " [", s = "]")))%>%
  compose(j = "dummy_SCD", value = as_paragraph(as_chunk(emmean_SCD), as_bracket(emmean_SCD_confint_low, emmean_SCD_confint_up, sep = "; ", p = " [", s = "]")))%>%
  compose(j = "dummy_MCI", value = as_paragraph(as_chunk(emmean_MCI), as_bracket(emmean_MCI_confint_low, emmean_MCI_confint_up, sep = "; ", p = " [", s = "]")))%>%
  compose(j = "dummy_AD", value = as_paragraph(as_chunk(emmean_AD), as_bracket(emmean_AD_confint_low, emmean_AD_confint_up, sep = "; ", p = " [", s = "]")))%>%
  add_footer_lines("Abbreviations: AAA = anterior amygdaloidarea, AD = Alzheimer’s disease, ANCOVA = analysis of covariance, CA = cornu ammonis, CATA = corticoamygdaloid transition area, CI = confidence interval, GC-ML-DG = Granule cell and molecular layer of the dentate gyrus, HATA = hippocampus-amygdala-transition-area, HC = healthy control, MCI = mild cognitive impairment, ROI = region of interest, SCD = subjective cognitive decline, TIV = total intracranial volume")%>%
  set_table_properties(width = 1, layout = "autofit")%>%
  set_header_labels(values = list(ROIs_names = "ROI",
                                  dummy_HC = "HC",
                                  dummy_SCD = "SCD",
                                  dummy_MCI = "MCI",
                                  dummy_AD = "AD",
                                  p_value = "p-value"))%>%
  autofit()
ft

#save_as_docx("TableA1" = ft, path = "[path]")

##### Data for Figure 1(% differences and 95%CIs) ####
coefficients_text <- function (anovas_emmeans, anovas_emmeans_confint){
  c('HC vs. SCD' = round(((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame((anovas_emmeans_confint)$contrasts)[1,2])*-1,2),
    'HC vs. SCD 95%CI lower' = round(((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame((anovas_emmeans_confint)$contrasts)[1,5])*-1,0),
    'HC vs. SCD 95%CI upper' = round(((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame((anovas_emmeans_confint)$contrasts)[1,6])*-1,0),
    'HC vs. SCD p' = round(as.data.frame((anovas_emmeans)$contrasts)[1,6],3),
    'HC vs. MCI' = round(((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame((anovas_emmeans_confint)$contrasts)[2,2])*-1,2),
    'HC vs. MCI 95%CI lower' = round(((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame((anovas_emmeans_confint)$contrasts)[2,5])*-1,0),
    'HC vs. MCI 95%CI upper' = round(((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame((anovas_emmeans_confint)$contrasts)[2,6])*-1,0),
    'HC vs. MCI p' = round(as.data.frame((anovas_emmeans)$contrasts)[2,6],3),
    'HC vs. AD' = round(((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame((anovas_emmeans_confint)$contrasts)[3,2])*-1,2),
    'HC vs. AD 95%CI lower' = round(((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame((anovas_emmeans_confint)$contrasts)[3,5])*-1,0),
    'HC vs. AD 95%CI upper' = round(((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame((anovas_emmeans_confint)$contrasts)[3,6])*-1,0),
    'HC vs. AD p' = round(as.data.frame((anovas_emmeans)$contrasts)[3,6],3)
    )
}
m <- mapply (coefficients_text, my_anovas_emmeans, my_anovas_emmeans_confint)
df <- data.frame(ROIs_names, as.data.frame(t(m)))
df$HC.vs..SCD.p[df$HC.vs..SCD.p <  0.001] <- "<0.001"
df$HC.vs..MCI.p[df$HC.vs..MCI.p <  0.001] <- "<0.001"
df$HC.vs..AD.p[df$HC.vs..AD.p <  0.001] <- "<0.001"

set_flextable_defaults(font.size = 9, digits = 0)
ft <- flextable(
  data = df, 
  col_keys = c("ROIs_names", "dummy_HCvsSCD", "HC.vs..SCD.p", "dummy_HCvsMCI", "HC.vs..MCI.p", "dummy_HCvsAD", "HC.vs..AD.p")) %>%
  compose(j = "dummy_HCvsSCD", value = as_paragraph(as_chunk(HC.vs..SCD), as_bracket(HC.vs..SCD.95.CI.upper, HC.vs..SCD.95.CI.lower, sep = "%; ", p = "% [", s = "%]")))%>%
  compose(j = "dummy_HCvsMCI", value = as_paragraph(as_chunk(HC.vs..MCI), as_bracket(HC.vs..MCI.95.CI.upper, HC.vs..MCI.95.CI.lower, sep = "%; ", p = "% [", s = "%]")))%>%
  compose(j = "dummy_HCvsAD", value = as_paragraph(as_chunk(HC.vs..AD), as_bracket(HC.vs..AD.95.CI.upper, HC.vs..AD.95.CI.lower, sep = "%; ", p = "% [", s = "%]")))%>%
  set_table_properties(width = 1, layout = "autofit")%>%
  set_header_labels(values = list(ROIs_names = "ROI",
                                  dummy_HCvsSCD = "HC vs. SCD",
                                  HC.vs..SCD.p = "p-value",
                                  dummy_HCvsMCI = "HC vs. MCI",
                                  HC.vs..MCI.p = "p-value",
                                  dummy_HCvsAD = "HC vs. AD",
                                  HC.vs..AD.p = "p-value"))%>%
  autofit()
ft

#save_as_docx("Figure1_data" = ft, path = "[path]")

##### Figure 1: Mean group differences of ROI volumes in percent (%). Percentages were calculated from pairwise comparisons of ANCOVAs adjusted to the covariates age, sex and education and to multiple comparison (Tukey). Percentages are highlighted in red for a relative smaller and blue for relative bigger volume. #####

coefficients_figure1 <- function(anovas_emmeans){
  c('HC vs. SCD' = round((100-((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame(anovas_emmeans)[2,3]))*-1,0),
    'HC vs. MCI' = round((100-((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame(anovas_emmeans)[3,3]))*-1,0),
    'HC vs. AD' = round((100-((100/as.data.frame(anovas_emmeans)[1,3])*as.data.frame(anovas_emmeans)[4,3]))*-1,0),
    'SCD vs. MCI' = round((100-((100/as.data.frame(anovas_emmeans)[2,3])*as.data.frame(anovas_emmeans)[3,3]))*-1,0),
    'MCI vs. AD' = round((100-((100/as.data.frame(anovas_emmeans)[3,3])*as.data.frame(anovas_emmeans)[4,3]))*-1,0))
}

m <- mapply (coefficients_figure1, my_anovas_emmeans)
df <- as.data.frame(m)%>%
  setNames(ROIs_names)%>%
  add_rownames()

colourer <- scales::col_numeric(
  palette = c("red", "transparent"),
  domain = c(-34, 0))

set_flextable_defaults(font.size = 8, digits = 0)
ft <- flextable(df, col_keys = c("rowname", "Amygdala", "Lateral Nucleus", "Basal Nucleus", "Accessory basal Nucleus", "Anterior Amygdaloid Area", "Central Nucleus",
                                 "Medial Nucleus", "Cortical Nucleus", "CATA", "Paralaminar Nucleus","Hippocampus", "CA1", "CA2/CA3", 
                                 "CA4GD", "Molecular Layer", "Subikulum", "GC_ML_DG", "HATA"))%>%
  set_table_properties(width = 1, layout = "autofit")%>%
  bg(j = 2:19, bg = colourer)%>%
  add_header_row(colwidths = c(1,10,8), values= c("","Amygdala ROIs",  "Hippocampus ROIs"))%>%
  #rotate(j = 2:17, rotation ="btlr", part="header")%>%
  border_outer(part="all", border = fp_border_default(color = "white") ) %>%
  border_inner(border = fp_border_default(color = "white"), part="all")%>%
  set_header_labels(values = list(rowname = "",
                                  'Amygdala' = "A1",
                                  'Lateral Nucleus' = "A2", 
                                  'Basal Nucleus' = "A3", 
                                  'Accessory basal Nucleus' = "A4", 
                                  'Anterior Amygdaloid Area' = "A5",
                                  'Central Nucleus' = "A6",
                                  'Medial Nucleus' = "A7", 
                                  'Cortical Nucleus' = "A8", 
                                  'CATA' = "A9",
                                  'Paralaminar Nucleus' = "A10",
                                  'Hippocampus' = "H1",
                                  'CA1' = "H2",
                                  'CA2/CA3' = "H3",
                                  'CA4GD' ="H4",
                                  'Molecular Layer' = "H5", 
                                  'Subikulum' = "H6",
                                  'GC_ML_DG' = "H7",
                                  'HATA' = "H8"))%>%
  add_footer_lines("Abbreviations: A1 = whole Amygdala, A2 = Lateral Nucleus, A3 = Basal Nucleus, A4 = Accessory Basal Nucleus, A5 = Anterior Amygdaloid Area (AAA), A6 = Central Nucleus, A7 = Medial Nucleus, A8 = Cortical Nucleus, A9 = Corticoamygdaloid Transition Area (CATA), A10 = Paralaminar Nucleus, H1 = whole Hippocampus, H2 = CA1, H3 = CA2/CA3, H4 = CA4GD, H5 = Molecular Layer, H6 = Subikulum, H7 = GC_ML_DG,
                                  H8 = HATA")%>%
  italic(italic = TRUE, part = "footer" )%>%
  bold(i = 1, bold = TRUE, part = "header") %>%
  vline(j = 11, border = fp_border_default(color = "white", style = "solid", width = 5), part = "all")%>%
  vline(j = 19, border = fp_border_default(color = "white", style = "solid", width = 5), part = "all")
ft <- autofit(ft)
ft

#save_as_docx("Figure1" = ft, path = "[path]")


##### Table A2: Association between baseline memory and baseline volume of ROIs. All volumes were adjusted to individual TIVs and linear regression models included age, sex and education. For comparability of the effects, ? was standardized according to their standard deviations (std. ?).######

ROIs <- Vol_1 %>%
  select(Whole_amygdala, Lat_nucl, Bas_nucl, Acces_Bas_nucl, AAA, Centr_nucl,
         Med_nucl, Cort_nucl, CATA, paral_nucl, Whole_Hippocampus, 
         CA1, CA3, CA4, molec_layer, Subiculum, GC_ML_DG, HATA, Precentral)

mean_ten_percent<- means_HC[2:20]*0.10 #means_HC includes data from T1 only
ROIs_ten_percent <- mapply(function (roi, mean_percent)(roi / mean_percent), ROIs, mean_ten_percent)
ROIs_ten_percent <- as.list(as.data.frame(ROIs_ten_percent))
ROIs_ten_percent <- as.data.frame(ROIs_ten_percent)
colnames(ROIs_ten_percent) <- paste(colnames(ROIs_ten_percent),"perc",sep="_")
ROIs_ten_percent <- data.frame("Pseudonym" = Vol_1$Pseudonym, ROIs_ten_percent)
Vol <- Vol %>% full_join(ROIs_ten_percent, by = "Pseudonym")


my_lms <- lapply(ROIs_ten_percent[2:20], function(roi)lm(Memory_vdRest ~ roi + Age + Education + Sex, Vol_1))

#nobs per group
lapply(Vol_1[36:54], function(roi)Vol_1%>%group_by(Diagnose_Baseline) %>% summarise(n = sum(!is.na(Memory_vdRest) & !is.na(roi))))
Vol_1%>% group_by(Diagnose_Baseline) %>% summarise(n = sum(!is.na(Memory_vdRest) & !is.na(Precentral)))

#control_my_lms <- lapply(my_lms, plot)

coefficients_tableA2 <- function(lm, lm.beta){
  c( n = nobs(lm), 
     beta = round(summary(lm)$coefficients[2],3),   
     confint_low = round(confint(lm)[2,1],2),       
     confint_up = round(confint(lm)[2,2],2),        
     p_value = round(summary(lm)$coefficients[2,4],3),
     #eta_squared = round(etasq(lm)[1,1],2),
     R_squared = round(summary(lm)$r.squared,2),
     adj_R_squared = round(summary(lm)$adj.r.squared,2))
}

m <- sapply (my_lms, coefficients_tableA2)
df <- data.frame(ROIs_names, as.data.frame(t(m)))
df$p_value[df$p_value <  0.001] <- "<0.001"

set_flextable_defaults(big.mark = " ", 
                       font.size = 10, digits = 2)
ft <- flextable(
  data = df, 
  col_keys = c( "ROIs_names","n", "dummy", "p_value", "eta_squared", "R_squared","adj_R_squared")) %>%
  compose(j = "dummy", value = as_paragraph(as_chunk(beta), as_bracket(confint_low, confint_up, sep = "; ", p = " [", s = "]")))%>%
  add_footer_lines("hello note")%>%
  set_table_properties(width = 1, layout = "autofit")%>%
  set_header_labels(values = list(ROIs_names = "ROI",
                                  dummy = "std. ß [95%CI]",
                                  eta_squared = "part. eta?",
                                  p_value = "p-value",
                                  R_squared = "R?",
                                  adj_R_squared = "adj. R?"))%>%
  autofit()
ft

#save_as_docx("TableA2" = ft, path = "[path])

##### Model Memory: lmer (Memory_vdRest ~  Diagnose_Baseline*year...####

Model2a <- lmer (Memory_vdRest ~  Diagnose_Baseline*year + Age + Education + Sex +  (1|ID), Vol, REML = TRUE, na.action = na.omit)
nobs(Model2a)
summary(Model2a)
emmeans_2a1 <- emmeans(Model2a, ~Diagnose_Baseline|year,at=list(year=c(0,1,3)), adjust = "none")

##### Table 2: Memory performance and number of participants at baseline and follow-up visits.#####
n <- Vol %>% filter (!is.na (Diagnose_Baseline)) %>% group_by(Visit, Diagnose_Baseline) %>% summarise(Memory_vdRest_n = sum(!is.na(Memory_vdRest)))

time_diff_V2 <- Vol_1%>% filter(!is.na (Diagnose_Baseline)) %>% inner_join(Vol_2, by = "Pseudonym")%>% transmute(Pseudonym, Diagnose_Baseline.x , V2_A_Time_difference = Date.y - Date.x)%>%group_by(Diagnose_Baseline.x) %>%summarise(V2_A_Time_difference_mean = mean(V2_A_Time_difference,na.rm=TRUE), V2_A_Time_difference_sd = sd(V2_A_Time_difference,na.rm=TRUE))%>% mutate(across(2:3, round, 0))%>%
  mutate("V2_A_Time_difference" = paste0(V2_A_Time_difference_mean, "(",V2_A_Time_difference_sd, ")"))%>%select(Diagnose_Baseline.x, V2_A_Time_difference)
time_diff_V2 <- time_diff_V2 %>% group_by(Diagnose_Baseline.x) %>% gather(column, value, -Diagnose_Baseline.x) %>% spread(Diagnose_Baseline.x, value) %>% mutate_all(as.factor)

time_diff_V3 <- Vol_1 %>% filter (!is.na (Diagnose_Baseline)) %>% inner_join(Vol_3, by = "Pseudonym")%>% transmute(Pseudonym, Diagnose_Baseline.x , V3_A_Time_difference = Date.y - Date.x)%>%group_by(Diagnose_Baseline.x) %>%summarise(V3_A_Time_difference_mean = mean(V3_A_Time_difference,na.rm=TRUE), V3_A_Time_difference_sd = sd(V3_A_Time_difference,na.rm=TRUE))%>% mutate(across(2:3, round, 0))%>%
  mutate("V3_A_Time_difference" = paste0(V3_A_Time_difference_mean, "(",V3_A_Time_difference_sd, ")"))%>%select(Diagnose_Baseline.x, V3_A_Time_difference)
time_diff_V3 <- time_diff_V3 %>% group_by(Diagnose_Baseline.x) %>% gather(column, value, -Diagnose_Baseline.x) %>% spread(Diagnose_Baseline.x, value) %>% mutate_all(as.factor)

em_table2 <- data.frame(emmeans_2a1)%>%
  select (Diagnose_Baseline, year, Memory = emmean, lower.CL, upper.CL)%>%
  cbind(n%>%select(Memory_vdRest_n), data.frame(test(emmeans_2a1))%>%select(p.value))
em_table2<- em_table2%>%
  mutate(Memory = paste0(round(Memory,2), "[", round(lower.CL,2), "; ", round(upper.CL,2), "]"), p= round(p.value, 3))%>%
  select( -lower.CL, -upper.CL, -p.value, -Visit)
options(knitr.kable.NA = "")
table_2 <- kable (em_table2) %>%
  kable_classic() %>%
  save_kable("[path]")
#manually added time_diff_V2 and time_diff_V3 in the table

##### Figure 2: Change of memory performance during three years assessed by a linear mixed model adjusted for age, sex and education. Participant's memory was assessed by a composite score of z-transformed memory subtests at baseline, one and three years after baseline. Abbreviations: AD = Alzheimer Dementia, HC = Healthy Control, MCI = Mild Cognitive Impairment, SCD = Subjective Cognitive Decline. ######

#for figure 2
em2a <- emmeans(Model2a, ~Diagnose_Baseline|year,at=list(year=c(0,1, 3.5)))
emtrends(Model2a, "Diagnose_Baseline", var= "year")                 #yearly memory change
test(emtrends(Model2a, "Diagnose_Baseline", var= "year"))           #p-values
r2beta(Model2a) #Rsq

EM2a <- as.data.frame(em2a)

figure2 <- ggplot (EM2a,aes(x=year, y = emmean, group=Diagnose_Baseline))+
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = Diagnose_Baseline), alpha = .50) +
  geom_line(aes(colour = Diagnose_Baseline), size = 2)+
  labs(x="Time [years]", y="Memory")+
  geom_line(data= Vol %>% filter(!is.na(Diagnose_Baseline) & !is.na(Memory_vdRest)), mapping=aes(x = year, y = Memory_vdRest, group= ID, col=Diagnose_Baseline), alpha= 0.5)+
  geom_point(data=Vol %>% filter(!is.na(Diagnose_Baseline) & !is.na(Memory_vdRest)), mapping=aes(x = year, y = Memory_vdRest, group= ID, col=Diagnose_Baseline), size=1, alpha= 0.5)+
  facet_wrap(~Diagnose_Baseline,nrow = 1)+
  scale_fill_manual(values =c("gray50", "seagreen", "deepskyblue3" ,"#225EA8"))+
  scale_color_manual(values =c("gray50", "seagreen", "deepskyblue3" ,"#225EA8"))+
  theme_bw(base_size=16)+
  theme(legend.position = 'none')
#save as .pdf device size 5 x 11 inches

##### Model 1: lmer(Memory_vdRest ~ roi*year + Age + Education + Sex + (1|ID)....####
library(lmerTest)


ROIs_perc <- as.list(Vol[63:81])  ####use of the scales ROIs because like this the interaction effects are comparable

#calc models
my_lmer <- lapply(ROIs_perc, function(roi)lmer(Memory_vdRest ~ roi*year + Age + Education + Sex + (1|ID), Vol, REML = TRUE, na.action = na.omit))
sum_my_lmer <- lapply(my_lmer, function(model)as.data.frame(summary(model)$coefficients))

##### Table A3: Model 1. Estimated rate of change of memory per year across the entire cohort predicted by baseline ROI volume. Linear mixed models were adjusted to age, sex and education. #### #################

my_Trends <- mapply(function(model)emmeans(model,~roi|year, at=list(roi=c(8,9,10), year=c(0,1,3))), my_lmer) #roi = 80%, 90% and 100%
my_Trends <- lapply(my_Trends, as.data.frame)

my_Trends_m1_yearly <- mapply(function(model)emtrends(model, "roi", var= "year", at = list(roi=c(8,9,10))), my_lmer) #yearly memory change
my_Trends_m1_p <- lapply(my_Trends_m1_yearly, test)#p-values
 

coefficients_tableA3 <- function(model_lmer,sum_model_lmer, Trends, Trends_p){
  std.Coef <- std.beta.lmer(model_lmer)
  c(beta_roiyear = round(sum_model_lmer[7,1],2),
    confint_low = round(confint(model_lmer)[9,1],2),
    confint_up = round(confint(model_lmer)[9,2],2),
    #std_beta = round(std.Coef[6],3),
    #p_value = round(sum_model_lmer[7,5],3),
    "20Vol" = round((as.data.frame(Trends)[1,2]),3),
    "20Vol_lower" = round((as.data.frame(Trends)[1,5]),2),
    "20Vol_upper"  = round((as.data.frame(Trends)[1,6]),2),
    "20Volp" = round((as.data.frame(Trends_p))[1,6],3),
    "10Vol" = round((as.data.frame(Trends)[2,2]),3),
    "10Vol_lower" = round((as.data.frame(Trends)[2,5]),2),
    "10Vol_upper" = round((as.data.frame(Trends)[2,6]),2),
    "10Volp" = round((as.data.frame(Trends_p))[2,6],3),
    meanHCVol = round((as.data.frame(Trends)[3,2]),3),
    meanHC_lower = round((as.data.frame(Trends)[3,5]),2),
    meanHC_upper = round((as.data.frame(Trends)[3,6]),2),
    meanHCp = round((as.data.frame(Trends_p))[3,6],3),
    Rsquared_roi = round(as.data.frame(r2beta(model_lmer))[6,6]+as.data.frame(r2beta(model_lmer))[2,6],2),
    Rsquared_model = round(as.data.frame(r2beta(model_lmer))[1,6],2))
}

m <- mapply(coefficients_tableA3, my_lmer, sum_my_lmer, my_Trends_m1_yearly,my_Trends_m1_p)
df <- data.frame(ROIs_names, as.data.frame(t(m)))
df$p_value[df$p_value <  0.001] <- "<0.001"
df$meanHCp[df$meanHCp <  0.001] <- "<0.001"
df$X10Volp[df$X10Volp <  0.001] <- "<0.001"
df$X20Volp[df$X20Volp <  0.001] <- "<0.001"
df_model1 <- df

set_flextable_defaults(font.size = 10, digits = 2)
ft <- flextable(
  data = df, 
  col_keys = c("ROIs_names", "dummy", "p_value", "dummy_meanHC", "meanHCp", "dummy_10", "X10Volp", "dummy_20", "X20Volp", "Rsquared_roi", "Rsquared_model")) %>%
  compose(j = "dummy", value = as_paragraph(as_chunk(beta_roiyear), as_bracket(confint_low, confint_up, sep = "; ", p = " [", s = "]")))%>%
  compose(j = "dummy_meanHC", value = as_paragraph(as_chunk(meanHCVol), as_bracket(meanHC_lower, meanHC_upper, sep = "; ", p = " [", s = "]")))%>%
  compose(j = "dummy_10", value = as_paragraph(as_chunk(`X10Vol`), as_bracket(`X10Vol_lower`, `X10Vol_upper`, sep = "; ", p = " [", s = "]")))%>%
  compose(j = "dummy_20", value = as_paragraph(as_chunk(`X20Vol`), as_bracket(`X20Vol_lower`, `X20Vol_upper`, sep = "; ", p = " [", s = "]")))%>%
  add_footer_lines("hello note")%>%
  set_table_properties(width = 1, layout = "autofit")%>%
  set_header_labels(values = list(ROIs_names = "ROI",
                                  dummy = "ß volume x time [95%CI]",
                                  #std_beta.roi.year = "std. ?",
                                  p_value = "p-value",
                                  dummy_meanHC = "mean HC",
                                  meanHCp = "p",
                                  "dummy_10" = "-10 %",
                                  "X10Volp" = "p",
                                  "dummy_20" = "-20 %",
                                  "X20Volp" = "p",
                                  Rsquared_roi = "R² volume + volume x time",
                                  Rsquared_model = "R² model"))%>%
  autofit()
ft

#save_as_docx("TableA3" = ft, path = "[path]")

##### Model 2: lmer(Memory_vdRest ~ roi*year*Diagnose + Age + Education + Sex + (1|ID)... ####
library(lmerTest)

my_lmer <- lapply(ROIs_perc, function(roi)lmer(Memory_vdRest ~ roi*year*Diagnose_Baseline + Age + Education + Sex + (1|ID), Vol, REML = TRUE, na.action = na.omit))
sum_my_lmer <- lapply(my_lmer, function(model)as.data.frame(summary(model)$coefficients))
nobs <- as.data.frame(lapply(my_lmer, nobs))


##### Table A4: Model 2. #####
my_Trends_m2_yearly <- mapply(function(model)emtrends(model, ~ Diagnose_Baseline | roi, 
                                                              var= "year", 
                                                              at = list(roi=c(8,9,10))), 
                                                              my_lmer) #yearly memory change
my_Trends_m2_p <- lapply(my_Trends_m2_yearly, test)#p-values

coefficients_table3 <- function(Trends, Trends_p){
  c(HCmeanHCVol = round((as.data.frame(Trends)[9,3]),2),
    HCmeanHC_lower = round((as.data.frame(Trends)[9,6]),2),
    HCmeanHC_upper = round((as.data.frame(Trends)[9,7]),2),
    HCmeanHCp = round((as.data.frame(Trends_p))[9,7],3),
    "HC10Vol" = round((as.data.frame(Trends)[5,3]),2),
    "HC10Vol_lower" = round((as.data.frame(Trends)[5,6]),2),
    "HC10Vol_upper" = round((as.data.frame(Trends)[5,7]),2),
    "HC10Volp" = round((as.data.frame(Trends_p))[5,7],3),
    "HC20Vol" = round((as.data.frame(Trends)[1,3]),2),
    "HC20Vol_lower" = round((as.data.frame(Trends)[1,6]),2),
    "HC20Vol_upper" = round((as.data.frame(Trends)[1,7]),2),
    "HC20Volp" = round((as.data.frame(Trends_p))[1,7],3),
    SCDmeanHCVol = round((as.data.frame(Trends)[10,3]),2),
    SCDmeanHC_lower = round((as.data.frame(Trends)[10,6]),2),
    SCDmeanHC_upper = round((as.data.frame(Trends)[10,7]),2),
    SCDmeanHCp = round((as.data.frame(Trends_p))[10,7],3),
    "SCD10Vol" = round((as.data.frame(Trends)[6,3]),2),
    "SCD10Vol_lower" = round((as.data.frame(Trends)[6,6]),2),
    "SCD10Vol_upper" = round((as.data.frame(Trends)[6,7]),2),
    "SCD10Volp" = round((as.data.frame(Trends_p))[6,7],3),
    "SCD20Vol" = round((as.data.frame(Trends)[2,3]),2),
    "SCD20Vol_lower" = round((as.data.frame(Trends)[2,6]),2),
    "SCD20Vol_upper" = round((as.data.frame(Trends)[2,7]),2),
    "SCD20Volp" = round((as.data.frame(Trends_p))[2,7],3),
    MCImeanHCVol = round((as.data.frame(Trends)[11,3]),2),
    MCImeanHC_lower = round((as.data.frame(Trends)[11,6]),2),
    MCImeanHC_upper = round((as.data.frame(Trends)[11,7]),2),
    MCImeanHCp = round((as.data.frame(Trends_p))[11,7],3),
    "MCI10Vol" = round((as.data.frame(Trends)[7,3]),2),
    "MCI10Vol_lower" = round((as.data.frame(Trends)[7,6]),2),
    "MCI10Vol_upper" = round((as.data.frame(Trends)[7,7]),2),
    "MCI10Volp" = round((as.data.frame(Trends_p))[7,7],3),
    "MCI20Vol" = round((as.data.frame(Trends)[3,3]),2),
    "MCI20Vol_lower" = round((as.data.frame(Trends)[3,6]),2),
    "MCI20Vol_upper" = round((as.data.frame(Trends)[3,7]),2),
    "MCI20Volp" = round((as.data.frame(Trends_p))[3,7],3),
    ADmeanHCVol = round((as.data.frame(Trends)[12,3]),2),
    ADmeanHC_lower = round((as.data.frame(Trends)[12,6]),2),
    ADmeanHC_upper = round((as.data.frame(Trends)[12,7]),2),
    ADmeanHCp = round((as.data.frame(Trends_p))[12,7],3),
    "AD10Vol" = round((as.data.frame(Trends)[8,3]),2),
    "AD10Vol_lower" = round((as.data.frame(Trends)[8,6]),2),
    "AD10Vol_upper" = round((as.data.frame(Trends)[8,7]),2),
    "AD10Volp" = round((as.data.frame(Trends_p))[8,7],3),
    "AD20Vol" = round((as.data.frame(Trends)[4,3]),2),
    "AD20Vol_lower" = round((as.data.frame(Trends)[4,6]),2),
    "AD20Vol_upper" = round((as.data.frame(Trends)[4,7]),2),
    "AD20Volp" = round((as.data.frame(Trends_p))[4,7],3)
    )
}

m <- mapply(coefficients_table3, my_Trends_m2_yearly, my_Trends_m2_p)
df <- data.frame(ROIs_names, as.data.frame(t(m)))%>% mutate_if(is.numeric, round, digits=3)
df$HCmeanHCp[df$HCmeanHCp <  0.001] <- "<0.001"
df$HC10Volp[df$HC10Volp <  0.001] <- "<0.001"
df$HC20Volp[df$HC20Volp <  0.001] <- "<0.001"
df$SCDmeanHCp[df$SCDmeanHCp <  0.001] <- "<0.001"
df$SCD10Volp[df$SCD10Volp <  0.001] <- "<0.001"
df$SCD20Volp[df$SCD20Volp <  0.001] <- "<0.001"
df$MCImeanHCp[df$MCImeanHCp <  0.001] <- "<0.001"
df$MCI10Volp[df$MCI10Volp <  0.001] <- "<0.001"
df$MCI20Volp[df$MCI20Volp <  0.001] <- "<0.001"
df$ADmeanHCp[df$ADmeanHCp <  0.001] <- "<0.001"
df$AD10Volp[df$AD10Volp <  0.001] <- "<0.001"
df$AD20Volp[df$AD20Volp <  0.001] <- "<0.001"

df_model2 <- df

set_flextable_defaults(font.size = 10, digits = 2)
ft <- flextable(
  data = df, 
  col_keys = c("ROIs_names", "HCdummy_meanHC", "HCmeanHCp", "HCdummy_10", "HC10Volp", "HCdummy_20", "HC20Volp",
               "SCDdummy_meanHC", "SCDmeanHCp", "SCDdummy_10", "SCD10Volp", "SCDdummy_20", "SCD20Volp",
               "MCIdummy_meanHC", "MCImeanHCp", "MCIdummy_10", "MCI10Volp", "MCIdummy_20", "MCI20Volp",
               "ADdummy_meanHC", "ADmeanHCp", "ADdummy_10", "AD10Volp", "ADdummy_20", "AD20Volp"
               )) %>%
  add_footer_lines("hello note")%>%
  compose(j = "HCdummy_meanHC", value = as_paragraph(as_chunk(HCmeanHCVol ), as_bracket(HCmeanHC_lower , HCmeanHC_upper , sep = "; ", p = " [", s = "]")))%>%
  compose(j = "HCdummy_10", value = as_paragraph(as_chunk(HC10Vol), as_bracket(HC10Vol_lower  , HC10Vol_upper , sep = "; ", p = " [", s = "]")))%>%
  compose(j = "HCdummy_20", value = as_paragraph(as_chunk(HC20Vol), as_bracket(HC20Vol_lower  , HC20Vol_upper , sep = "; ", p = " [", s = "]")))%>%
  compose(j = "SCDdummy_meanHC", value = as_paragraph(as_chunk(SCDmeanHCVol ), as_bracket(SCDmeanHC_lower , SCDmeanHC_upper , sep = "; ", p = " [", s = "]")))%>%
  compose(j = "SCDdummy_10", value = as_paragraph(as_chunk(SCD10Vol), as_bracket(SCD10Vol_lower  , SCD10Vol_upper , sep = "; ", p = " [", s = "]")))%>%
  compose(j = "SCDdummy_20", value = as_paragraph(as_chunk(SCD20Vol), as_bracket(SCD20Vol_lower  , SCD20Vol_upper , sep = "; ", p = " [", s = "]")))%>%   
  compose(j = "MCIdummy_meanHC", value = as_paragraph(as_chunk(MCImeanHCVol ), as_bracket(MCImeanHC_lower , MCImeanHC_upper , sep = "; ", p = " [", s = "]")))%>%
  compose(j = "MCIdummy_10", value = as_paragraph(as_chunk(MCI10Vol), as_bracket(MCI10Vol_lower  , MCI10Vol_upper , sep = "; ", p = " [", s = "]")))%>%
  compose(j = "MCIdummy_20", value = as_paragraph(as_chunk(MCI20Vol), as_bracket(MCI20Vol_lower  , MCI20Vol_upper , sep = "; ", p = " [", s = "]")))%>%    
  compose(j = "ADdummy_meanHC", value = as_paragraph(as_chunk(ADmeanHCVol ), as_bracket(ADmeanHC_lower , ADmeanHC_upper , sep = "; ", p = " [", s = "]")))%>%
  compose(j = "ADdummy_10", value = as_paragraph(as_chunk(AD10Vol), as_bracket(AD10Vol_lower  , AD10Vol_upper , sep = "; ", p = " [", s = "]")))%>%
  compose(j = "ADdummy_20", value = as_paragraph(as_chunk(AD20Vol), as_bracket(AD20Vol_lower  , AD20Vol_upper , sep = "; ", p = " [", s = "]")))%>%
  set_table_properties(width = 1, layout = "autofit")%>%
  add_header_row(colwidths = c(1,6,6,6,6), values= c("","HC", "SCD", "MCI", "AD"))%>%
  set_header_labels(values = list(ROIs_names = "ROI",
                                  HCdummy_meanHC = "mean",
                                  HCmeanHCp = "p",
                                  HCdummy_10 = "-10%",
                                  HC10Volp = "p",
                                  HCdummy_20 = "-20%",
                                  HC20Volp = "p",
                                  SCDdummy_meanHC = "mean",
                                  SCDmeanHCp = "p",
                                  SCDdummy_10 = "-10%",
                                  SCD10Volp = "p",
                                  SCDdummy_20 = "-20%",
                                  SCD20Volp = "p",
                                  MCIdummy_meanHC = "mean",
                                  MCImeanHCp = "p",
                                  MCIdummy_10 = "-10%",
                                  MCI10Volp = "p",
                                  MCIdummy_20 = "-20%",
                                  MCI20Volp = "p",
                                  ADdummy_meanHC = "mean",
                                  ADmeanHCp = "p",
                                  ADdummy_10 = "-10%",
                                  AD10Volp = "p",
                                  ADdummy_20 = "-20%",
                                  AD20Volp = "p"
                                  ))%>%
  autofit()
ft

#save_as_docx("TableA4" = ft, path = "[path]")

##### Model 2 without scaling ##### 
#scaling of volumes to 10% is not needed, because we don't want to interpret interaction term

ROIs <- Vol %>%
  select(Whole_amygdala, Lat_nucl, Bas_nucl, Acces_Bas_nucl, AAA, Centr_nucl,
         Med_nucl, Cort_nucl, CATA, paral_nucl, Whole_Hippocampus, 
         CA1, CA3, CA4, molec_layer, Subiculum, GC_ML_DG, HATA, Precentral)
ROIs <-  as.list(ROIs[1:19])

means_HC <- Vol %>%
  filter(Visit=="1" & Diagnose_Baseline =="HC")%>%
  select(Diagnose_Baseline, Whole_amygdala, Lat_nucl, Bas_nucl, Acces_Bas_nucl, AAA, Centr_nucl,
         Med_nucl, Cort_nucl, CATA, paral_nucl, Whole_Hippocampus, 
         CA1, CA3, CA4, molec_layer, Subiculum, GC_ML_DG, HATA, Precentral)%>%
  group_by(Diagnose_Baseline)%>%
  summarise_at (vars(Whole_amygdala:Precentral),mean, na.rm = TRUE)%>%
  mutate()

breaks <- function(mean){
  c(mean, mean-0.1*mean, mean-0.2*mean)
}# use as breaks(means_HC$Precentral)
breaks_Vol <- lapply(as.list(means_HC[2:20]), breaks)

my_lmer <- lapply(ROIs, function(roi)lmer(Memory_vdRest ~ roi*year*Diagnose_Baseline + Age + Education + Sex + (1|ID), Vol, REML = TRUE, na.action = na.omit))
sum_my_lmer <- lapply(my_lmer, function(model)as.data.frame(summary(model)$coefficients))

##### Figure 3: Model 2 without scaling #####
my_Trends <- mapply(function(model, breaks)(emtrends(model,specs = pairwise ~roi*Diagnose_Baseline, var = "year",at=list(roi=breaks)))$emtrends, my_lmer, breaks_Vol)
Trends_Figure3 <- as.data.frame(my_Trends)%>%
  select(Volume = Whole_amygdala.roi, Diagnose_Baseline = Whole_amygdala.Diagnose_Baseline, ends_with("trend"))%>%
  pivot_longer(cols = -c(1,2), names_to =  "ROI", values_to = "Memory_change_per_year")%>%
  mutate (ROI = factor (ROI, levels = names(as.data.frame(my_Trends)%>%select(ends_with("trend"))), labels = ROIs_names), 
          Diagnose_Baseline = factor(Diagnose_Baseline, levels = c("HC", "SCD", "MCI", "AD")),
          Volume = factor(Volume, labels = c("-20%", "-10%", "HC mean")))

lCL_Figure3 <- as.data.frame(my_Trends)%>%
  select(Volume = Whole_amygdala.roi, Diagnose_Baseline = Whole_amygdala.Diagnose_Baseline, ends_with("lower.CL"))%>%
  pivot_longer(cols = -c(1,2), names_to = "ROI", values_to = "lower_CL")%>%
  mutate (ROI = factor (ROI, levels = names(as.data.frame(my_Trends)%>%select(ends_with("lower.CL"))), labels = ROIs_names), 
          Diagnose_Baseline = factor(Diagnose_Baseline, levels = c("HC", "SCD", "MCI", "AD")),
          Volume = factor(Volume, labels = c("-20%", "-10%", "HC mean")))

uCL_Figure3 <- as.data.frame(my_Trends)%>%
  select(Volume = Whole_amygdala.roi, Diagnose_Baseline = Whole_amygdala.Diagnose_Baseline, ends_with("upper.CL"))%>%
  pivot_longer(cols = -c(1,2), names_to = "ROI", values_to = "upper_CL")%>%
  mutate (ROI = factor (ROI, levels = names(as.data.frame(my_Trends)%>%select(ends_with("upper.CL"))), labels = ROIs_names), 
          Diagnose_Baseline = factor(Diagnose_Baseline, levels = c("HC", "SCD", "MCI", "AD")),
          Volume = factor(Volume, labels = c("-20%", "-10%", "HC mean")))

df_Figure3 <- Trends_Figure3 %>% 
  inner_join(lCL_Figure3, by = c("Volume", "Diagnose_Baseline", "ROI"))%>% 
  inner_join(uCL_Figure3, by = c("Volume", "Diagnose_Baseline", "ROI"))%>%
  filter(Diagnose_Baseline != "AD")

#Add data from model_1
df_Figure3_model1_HCmean <- df_model1%>%
  select(ROIs_names, "Memory_change_per_year" = meanHCVol, lower_CL = meanHC_lower, upper_CL = meanHC_upper)%>%
  mutate(Volume = "HC mean", Diagnose_Baseline = "Total")%>%
  select(Volume, Diagnose_Baseline, ROI = ROIs_names, "Memory_change_per_year", lower_CL, upper_CL)
df_Figure3_model1_10 <- df_model1%>%
  select(ROIs_names, "Memory_change_per_year" = X10Vol, lower_CL = X10Vol_lower, upper_CL = X10Vol_upper)%>%
  mutate(Volume = "-10%", Diagnose_Baseline = "Total")%>%
  select(Volume, Diagnose_Baseline, ROI = ROIs_names, "Memory_change_per_year", lower_CL, upper_CL)
df_Figure3_model1_20 <- df_model1%>%
  select(ROIs_names, "Memory_change_per_year" = X20Vol, lower_CL = X20Vol_lower, upper_CL = X20Vol_upper)%>%
  mutate(Volume = "-20%", Diagnose_Baseline = "Total")%>%
  select(Volume, Diagnose_Baseline, ROI = ROIs_names, "Memory_change_per_year", lower_CL, upper_CL)

df_Figure3 <- df_Figure3 %>% bind_rows(df_Figure3_model1_HCmean,df_Figure3_model1_10, df_Figure3_model1_20)
df_Figure3$Diagnose_Baseline<- factor(df_Figure3$Diagnose_Baseline, levels = c("Total", "HC", "SCD","MCI"))
ROIs_names[ROIs_names == "Anterior Amygdaloid Area"] <- "AAA"
ROIs_names <- ROIs_names[1:19]
df_Figure3$ROI<- factor(df_Figure3$ROI, levels = ROIs_names)
df_Figure3$ROI[is.na(df_Figure3$ROI)] <- "AAA"
df_Figure3$Volume<- factor(df_Figure3$Volume, levels = c("HC mean", "-10%", "-20%"))

Figure3 <- df_Figure3 %>%
  ggplot(mapping= aes(x = Memory_change_per_year,
                      y = ROI,
                      color = Volume))+
  geom_point(position = position_dodge(width =0.5), size = 0.8)+
  geom_linerange(aes(xmin=lower_CL, xmax=upper_CL), position = position_dodge(width =0.5), size= 0.5)+
  xlab(label = "yearly memory change")+
  ylab(NULL)+
  #guides(color = guide_legend(title = "ROI volume"))+
  facet_grid(~Diagnose_Baseline)+#, scales = "free_x", space = "free_x")+
  ylim(rev(levels(df_Figure3$ROI)))+
  geom_vline(xintercept = 0, color = "dark grey")+
  #scale_color_manual(name="ROI volume")+
  theme_minimal()+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
        legend.position = "top",
        legend.text=element_text(size=9),
        legend.key.size =  unit(0.5, 'cm'),
        legend.title = element_text(size  = 9),
        #axis.text.y = element_text(face = c('bold', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain','bold', 'plain', 'plain',
        #                                       'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'bold'))
        axis.text.y= element_blank(),
        axis.text.x= element_text(size  = 7),
        axis.title.x = element_text(vjust = -2))+
  theme(panel.spacing = unit(1.2, "lines"))
Figure3
#ggsave("Figure3.pdf") # (safe inches 4.7 x 7.2)

#Add nobs in figure 3:
ROI_levels_Figure_A3 <- c("ROI","Amygdala", "Lateral Nucleus", "Basal Nucleus", "Accessory basal Nucleus", "AAA", "Central Nucleus",
                          "Medial Nucleus", "Cortical Nucleus", "CATA", "Paralaminar Nucleus","Hippocampus", "CA1", "CA3", 
                          "CA4", "Molecular Layer", "Subiculum","GC-ML-DG","HATA", "Precentral Gyrus")

df_Figure3_nobsdata <- ROIs_names %>% cbind(nobs%>% select(nparticipants, nobs))
df_Figure3_nobsdata <- rename(df_Figure3_nobsdata, ROIs_names = `.`)
p_names <- df_Figure3_nobsdata %>%
  mutate(ROIs_names = factor(ROIs_names, levels = rev(ROI_levels_Figure_A3)))%>%
  ggplot(aes(y = ROIs_names))+
  geom_text(aes(x = 0, label = ROIs_names), hjust = "right", size = 3, fontface = c('bold', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain','bold',
                                                                    'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'plain', 'bold')
   )+
  theme_void()
p_part<-df_Figure3_nobsdata %>%
  mutate(ROIs_names = factor(ROIs_names, levels = rev(ROI_levels_Figure_A3)))%>%
  ggplot()+
  geom_text(
    aes(x = 0, 
        y = ROIs_names, label = nparticipants),
    hjust = 0, size = 3,
    fontface = "plain")+
  theme_void() +
  ggtitle("n \n part.")+ theme(plot.title = element_text(size = 8, face = "bold", hjust =0.7, vjust = -30))
p_obs<-df_Figure3_nobsdata %>%
  mutate(ROIs_names = factor(ROIs_names, levels = rev(ROI_levels_Figure_A3)))%>%
  ggplot()+
  geom_text(
    aes(x = 0, y = ROIs_names, label = nobs),
    hjust = 0, size = 3,
    fontface = "plain")+
  theme_void() +
  ggtitle("n  \n obs.")+ theme(plot.title = element_text(size = 8, face = "bold", hjust =0.75, vjust = -30))
layout <- c(
  area(t = 0, l = 0, b = 28, r = 8), 
  area(t = 0, l = 5, b = 28, r = 6),
  area(t = 0, l = 6, b = 28, r = 7), 
  area(t = 0, l = 8, b = 28, r = 21)
  )
# final plot arrangement
Figure3b <- p_names + p_part + p_obs + Figure3 + plot_layout(design = layout)
ggsave("Figure3b.png", height = 4.7, width = 7.2,bg="white", dpi = 600) # (safe inches 4.7 x 7.2)

