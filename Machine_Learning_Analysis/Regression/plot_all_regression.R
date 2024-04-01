
setwd('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/')
load("//volume/projects/JF_ECNP/Analysis/HC_SCZ_manuscript/Models_new/project_17/regression_age_unif_dist_age_sel_PID.RData")
data0 <- read.csv('/volume/projects/JF_ECNP/Analysis/HC_SCZ_manuscript/Models_new/project_17/project_17/data_classif_clin.csv')

allatlas_HCnorm=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/allatlas_HCnorm_formanuscript.csv')
allatlas_HCrem=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/allatlas_HCremaining_formanuscript.csv')
allatlas_SCZ=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/allatlas_SCZ_formanuscript.csv')

schaefer_HCnorm=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/schaefer_HCnorm_formanuscript.csv')
schaefer_HCrem=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/schaefer_HCremaining_formanuscript.csv')
schaefer_SCZ=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/schaefer_SCZ_formanuscript.csv')

aal3_HCnorm=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/aal3_HCnorm_formanuscript.csv')
aal3_HCrem=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/aal3_HCremaining_formanuscript.csv')
aal3_SCZ=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/aal3_SCZ_formanuscript.csv')

hammers_HCnorm=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/hammers_HCnorm_formanuscript.csv')
hammers_HCrem=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/hammers_HCremaining_formanuscript.csv')
hammers_SCZ=read.csv('/volume/projects/JF_ECNP/Analysis/Analysis/R_models_150224/Regression/hammers_SCZ_formanuscript.csv')

#histogram
data_unif<-data0[data0$PID %in% sel_PID, ]

png("histograms.png")  # Specify the file name and resolution

par(mfrow = c(1, 2))  # Set up the plotting layout to have 1 row and 2 columns

hist(data_unif$AGE, 
     main = "Uniform-like sample",
     xlab = "Age",
     ylab = "Frequency",
     col = "cornflowerblue",  # Fill bars with blue color
     border = "black",        # Border color of bars
     cex.lab = 1.5,
     cex.main = 1.5,
     font.lab = 2) # Border color of bars

hist(data0$AGE, 
     main = "Complete sample",
     xlab = "Age",
     ylab = "Frequency",
     col = "cornflowerblue",  # Fill bars with blue color
     border = "black",        # Border color of bars
     cex.lab = 1.5,
     cex.main = 1.5,
     font.lab = 2) # Border color of bars

dev.off()  # Close the PNG device
## all atlas
original_ages=allatlas_HCnorm$original_age
average_predictions=allatlas_HCnorm$predicted_age_corr

original_ages_rem=allatlas_HCrem$original_age
average_predictions_rem=allatlas_HCrem$predicted_age_corr

original_age_scz=allatlas_SCZ$original_age
average_predictions_patients_scz=allatlas_SCZ$predicted_age_corr

#plot
#HC vs scz_allsites
Combined_data <- data.frame(
  Age=c(original_ages,original_ages_rem,original_age_scz),
  Age_pred = c(average_predictions,average_predictions_rem,average_predictions_patients_scz),
  Group = factor(
    c(rep("HC normative", length(original_ages)),
      rep("HC remaining", length(original_ages_rem)),
      rep("Schizophrenia", length(original_age_scz))),
    levels = c("HC normative","HC remaining" ,"Schizophrenia")))

ggplot(Combined_data, aes(x = Age, y = Age_pred, color = Group, fill = Group)) +
  geom_point(aes(alpha=0.55,color=Group)) +                            # Scatter plot
  geom_smooth(method = "lm", se = TRUE, level = 0.95, aes(group=Group,fill = Group)) + # Regression line (linear model) with fill aesthetics
  
  labs(x = "Chronological age", y = "Predicted age",title ="All atlas") +  # Add a title
  scale_color_manual(values = c("HC normative" = "royalblue", "HC remaining" = "gray37", "Schizophrenia" = "indianred2")) +
  scale_fill_manual(values = c("HC normative" = "royalblue", "HC remaining" = "gray37", "Schizophrenia" = "indianred2")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'),
        legend.text = element_text(size = 10, face = 'bold'),
        plot.title = element_text(size = 10, face = 'bold'),
        legend.title = element_text(size = 12, face = 'bold'))+
  guides(alpha="none")

ggsave('regression_unif_dist_allatlas_groupoffset_hc_remaining_scz.tiff', width = 7, height = 4, dpi = 500)
ggsave('regression_unif_dist_allatlas_groupoffset_hc_remaining_scz.png', width = 7, height = 4, dpi = 300)

## schaefer
original_ages=schaefer_HCnorm$original_age
average_predictions=schaefer_HCnorm$predicted_age_corr

original_ages_rem=schaefer_HCrem$original_age
average_predictions_rem=schaefer_HCrem$predicted_age_corr

original_age_scz=schaefer_SCZ$original_age
average_predictions_patients_scz=schaefer_SCZ$predicted_age_corr

#plot
#HC vs scz_allsites
Combined_data <- data.frame(
  Age=c(original_ages,original_ages_rem,original_age_scz),
  Age_pred = c(average_predictions,average_predictions_rem,average_predictions_patients_scz),
  Group = factor(
    c(rep("HC normative", length(original_ages)),
      rep("HC remaining", length(original_ages_rem)),
      rep("Schizophrenia", length(original_age_scz))),
    levels = c("HC normative","HC remaining" ,"Schizophrenia")))

ggplot(Combined_data, aes(x = Age, y = Age_pred, color = Group, fill = Group)) +
  geom_point(aes(alpha=0.55,color=Group)) +                            # Scatter plot
  geom_smooth(method = "lm", se = TRUE, level = 0.95, aes(group=Group,fill = Group)) + # Regression line (linear model) with fill aesthetics
  
  labs(x = "Chronological age", y = "Predicted age",title ="Schaefer atlas") +  # Add a title
  scale_color_manual(values = c("HC normative" = "royalblue", "HC remaining" = "gray37", "Schizophrenia" = "indianred2")) +
  scale_fill_manual(values = c("HC normative" = "royalblue", "HC remaining" = "gray37", "Schizophrenia" = "indianred2")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'),
        legend.text = element_text(size = 10, face = 'bold'),
        plot.title = element_text(size = 10, face = 'bold'),
        legend.title = element_text(size = 12, face = 'bold'))+
  guides(alpha="none")

ggsave('regression_unif_dist_schaefer_groupoffset_hc_remaining_scz.png', width = 7, height = 4, dpi = 300)

## aal3
original_ages=aal3_HCnorm$original_age
average_predictions=aal3_HCnorm$predicted_age_corr

original_ages_rem=aal3_HCrem$original_age
average_predictions_rem=aal3_HCrem$predicted_age_corr

original_age_scz=aal3_SCZ$original_age
average_predictions_patients_scz=aal3_SCZ$predicted_age_corr

#plot
#HC vs scz_allsites
Combined_data <- data.frame(
  Age=c(original_ages,original_ages_rem,original_age_scz),
  Age_pred = c(average_predictions,average_predictions_rem,average_predictions_patients_scz),
  Group = factor(
    c(rep("HC normative", length(original_ages)),
      rep("HC remaining", length(original_ages_rem)),
      rep("Schizophrenia", length(original_age_scz))),
    levels = c("HC normative","HC remaining" ,"Schizophrenia")))

ggplot(Combined_data, aes(x = Age, y = Age_pred, color = Group, fill = Group)) +
  geom_point(aes(alpha=0.55,color=Group)) +                            # Scatter plot
  geom_smooth(method = "lm", se = TRUE, level = 0.95, aes(group=Group,fill = Group)) + # Regression line (linear model) with fill aesthetics
  
  labs(x = "Chronological age", y = "Predicted age",title ="AAL3 atlas") +  # Add a title
  scale_color_manual(values = c("HC normative" = "royalblue", "HC remaining" = "gray37", "Schizophrenia" = "indianred2")) +
  scale_fill_manual(values = c("HC normative" = "royalblue", "HC remaining" = "gray37", "Schizophrenia" = "indianred2")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'),
        legend.text = element_text(size = 10, face = 'bold'),
        plot.title = element_text(size = 10, face = 'bold'),
        legend.title = element_text(size = 12, face = 'bold'))+
  guides(alpha="none")

ggsave('regression_unif_dist_aal3_groupoffset_hc_remaining_scz.png', width = 7, height = 4, dpi = 300)


## hammers
original_ages=hammers_HCnorm$original_age
average_predictions=hammers_HCnorm$predicted_age_corr

original_ages_rem=hammers_HCrem$original_age
average_predictions_rem=hammers_HCrem$predicted_age_corr

original_age_scz=hammers_SCZ$original_age
average_predictions_patients_scz=hammers_SCZ$predicted_age_corr

#plot
#HC vs scz_allsites
Combined_data <- data.frame(
  Age=c(original_ages,original_ages_rem,original_age_scz),
  Age_pred = c(average_predictions,average_predictions_rem,average_predictions_patients_scz),
  Group = factor(
    c(rep("HC normative", length(original_ages)),
      rep("HC remaining", length(original_ages_rem)),
      rep("Schizophrenia", length(original_age_scz))),
    levels = c("HC normative","HC remaining" ,"Schizophrenia")))

ggplot(Combined_data, aes(x = Age, y = Age_pred, color = Group, fill = Group)) +
  geom_point(aes(alpha=0.55,color=Group)) +                            # Scatter plot
  geom_smooth(method = "lm", se = TRUE, level = 0.95, aes(group=Group,fill = Group)) + # Regression line (linear model) with fill aesthetics
  
  labs(x = "Chronological age", y = "Predicted age",title ="Hammers atlas") +  # Add a title
  scale_color_manual(values = c("HC normative" = "royalblue", "HC remaining" = "gray37", "Schizophrenia" = "indianred2")) +
  scale_fill_manual(values = c("HC normative" = "royalblue", "HC remaining" = "gray37", "Schizophrenia" = "indianred2")) +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1, face = 'bold'),
        axis.text = element_text(size = 10, face = 'bold'),
        axis.title = element_text(size = 12, face = 'bold'),
        legend.text = element_text(size = 10, face = 'bold'),
        plot.title = element_text(size = 10, face = 'bold'),
        legend.title = element_text(size = 12, face = 'bold'))+
  guides(alpha="none")

ggsave('regression_unif_dist_hammers_groupoffset_hc_remaining_scz.png', width = 7, height = 4, dpi = 300)