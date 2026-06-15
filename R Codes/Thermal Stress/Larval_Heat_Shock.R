library(ggplot2)
#library(ggsignif)
library(readxl)
#library(FSA)
library(stringr)
library(patchwork)
library(tidyverse)
library(car)
library(lattice)
library(rstatix)
library(ggpubr)
library(broom)
library(flextable)
library(emmeans)
library(MASS)
library(dplyr)


Colors = c(
  "STD" ="#00b4d8",
  "3RP"="#f25c54",
  "3RP_HET_1"="#ffbf69",
  "3RP_HET_2"="#ffdab9",
  "2Lt"="#52b788",
  "2Lt_HET_1"="#80ed99",
  "2Lt_HET_2"="#c7f9cc",
  "3RK"="#8367c7",
  "3RK_HET_1"="#e0aaff",
  "3RK_HET_2"="#ebd9fc"
)

plot_theme <- theme_bw()+
  theme(axis.line = element_line(color = "black", linewidth = 0.8),
        axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12)
  )

y_axis <- scale_y_continuous(
  limits = c(0, 0.8),
  breaks = c(0, 0.2, 0.4, 0.6, 0.8)
)


# Import Data

DataSLHS <- read.csv2(file.choose())
SLHS <- DataSLHS


SLHS$Total.survival <- SLHS$Adults.after.experiment + SLHS$Pupae.after.experiment
SLHS$Survival.rate <- SLHS$Total.survival/ SLHS$Total.eggs
SLHS$Adult.Survival.rate <- SLHS$Adults.after.experiment/SLHS$Total.eggs

# Assumptions

M_Survival_Rate <- lm(Survival.rate~Karyotype, SLHS)
shapiro.test(M_Survival_Rate$residuals)

par(mfrow=c(1,2))
hist(aov(M_Survival_Rate)$residuals)
qqPlot(aov(M_Survival_Rate)$residuals, id=FALSE)

leveneTest(Survival.rate~Karyotype, SLHS)

# Analysis -- Non-parametric test

kruskal.test(Survival.rate ~ Karyotype, SLHS)
dunn_test(SLHS, Survival.rate ~ Karyotype)

D_SLHS <- as.data.frame(dunn_test(SLHS, Survival.rate ~ Karyotype))
S_SLHS <- D_SLHS[D_SLHS$p.adj<0.05,]

# Contrast Test

SLHS_3RP <- SLHS[c(which(str_detect(SLHS$Karyotype, "3RP|STD"))),]
SLHS_2Lt <- SLHS[c(which(str_detect(SLHS$Karyotype, "2Lt|STD"))),]
SLHS_3RK <- SLHS[c(which(str_detect(SLHS$Karyotype, "3RK|STD"))),]


a.18 = 1-(1-0.05)^(1/18)

All.comparisons.3RP.SLHS <- dunn_test(SLHS_3RP, Survival.rate ~ Karyotype, p.adjust.method = "none")
All.comparisons.2Lt.SLHS <- dunn_test(SLHS_2Lt, Survival.rate ~ Karyotype, p.adjust.method = "none")
All.comparisons.3RK.SLHS <- dunn_test(SLHS_3RK, Survival.rate ~ Karyotype, p.adjust.method = "none")

Contrast.SLHS <- rbind(All.comparisons.3RP.SLHS, All.comparisons.2Lt.SLHS, All.comparisons.3RK.SLHS)
S_Contrast.SLHS <- Contrast.SLHS[Contrast.SLHS$p< a.18,]
S_Contrast.SLHS <- mutate(S_Contrast.SLHS, p.signif =case_when(p< a.18 ~ "*", TRUE ~ "ns"))

# Add p-value

SLHS_3RP <- SLHS[c(which(str_detect(SLHS$Karyotype, "3RP|STD"))),]
SLHS_2Lt <- SLHS[c(which(str_detect(SLHS$Karyotype, "2Lt|STD"))),]
SLHS_3RK <- SLHS[c(which(str_detect(SLHS$Karyotype, "3RK|STD"))),]
SLHS_Homoz <- SLHS[-c(which(str_detect(SLHS$Karyotype, "HET"))),]


p_SLHS_3RP <- S_Contrast.SLHS[S_Contrast.SLHS$group1 %in% c("STD","3RP","3RP_HET_1", "3RP_HET_2") & S_Contrast.SLHS$group2 %in% c("STD","3RP","3RP_HET_1", "3RP_HET_2"),]
p_SLHS_2Lt <- S_Contrast.SLHS[S_Contrast.SLHS$group1 %in% c("STD","2Lt","2Lt_HET_1", "2Lt_HET_2") & S_Contrast.SLHS$group2 %in% c("STD","2Lt","2Lt_HET_1", "2Lt_HET_2"),]
p_SLHS_3RK <- S_Contrast.SLHS[S_Contrast.SLHS$group1 %in% c("STD","3RK","3RK_HET_1", "3RK_HET_2") & S_Contrast.SLHS$group2 %in% c("STD","3RK","3RK_HET_1", "3RK_HET_2"),]
p_SLHS_Homoz <- S_Contrast.SLHS[S_Contrast.SLHS$group1 %in% c("STD","3RP","2Lt","3RK") & S_Contrast.SLHS$group2 %in% c("STD","3RP","2Lt","3RK"),]


# Box plot 

SLHS_b1 <- ggplot(SLHS_3RP, aes(x = Karyotype, y = Survival.rate, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Survival rate in %" )+
  stat_pvalue_manual(p_SLHS_3RP, label="p.signif", hide.ns = TRUE, y.position = 0.7, step.increase = 0.08)


SLHS_b2 <- ggplot(SLHS_2Lt, aes(x = Karyotype, y = Survival.rate, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  theme_bw()+
  plot_theme+
  ylab("Survival rate in %" )+
  stat_pvalue_manual(p_SLHS_2Lt, label="p.signif", hide.ns = TRUE, y.position=0.80, step.increase = 0.08)


SLHS_b3 <- ggplot(SLHS_3RK, aes(x = Karyotype, y = Survival.rate, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Survival rate in %" )+
  stat_pvalue_manual(p_SLHS_3RK, label="p.signif", hide.ns = TRUE, y.position = 0.8,  step.increase = 0.08)


SLHS_b4 <- ggplot(SLHS_Homoz, aes(x = Karyotype, y = Survival.rate, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Survival rate in %" )+
  stat_pvalue_manual(p_SLHS_Homoz, label="p.signif", hide.ns = TRUE, y.position = 0.65, step.increase = 0.08)


Final <-(SLHS_b1|SLHS_b2)/(SLHS_b3|SLHS_b4) + plot_annotation(tag_levels = "A",
                                    theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),plot.caption = element_text(hjust=0.5, size =12)))
                                    #+ plot_annotation(title="Larvae Survival after a Heat Shock")

# Export plot as png file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/Stress_Resistance_Experiment/Analysis/2026_Stress_Resistance/Graphs _and_Tables/Word/Larval_Heat_Shock.png",
  plot = Final,          
  width = 40,                    
  height = 20,                    
  units = "cm",                   
  dpi = 300                        
)





