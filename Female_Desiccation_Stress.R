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
  limits = c(0, 1.25),
  breaks = c(0.2, 0.4, 0.6, 0.8,1)
)


# Import & Prepare Data

DataSADS.F <- read.csv2(file.choose())
SADS.F <- DataSADS.F
SADS.F <- DataSADS.F[-c(which(is.na(DataSADS.F$Total.Adults))),]


SADS.F$F.Survival <- SADS.F$ALIVE.F.after.exp/SADS.F$TOTAL.F


# Assumptions

M_SADS.F <- lm(F.Survival ~ Karyotype, SADS.F)

shapiro.test(aov(M_SADS.F)$residuals)

par(mfrow=c(1,2))
hist(aov(M_SADS.F)$residuals)
qqPlot(aov(M_SADS.F)$residuals, id=FALSE)

leveneTest(F.Survival ~ Karyotype, SADS.F)

# Statistical Analysis -- Non-parametric test
kruskal.test(F.Survival ~ Karyotype, SADS.F)
dunn_test(SADS.F, F.Survival ~ Karyotype)

D_SADS.F <- as.data.frame(dunn_test(SADS.F, F.Survival ~ Karyotype))
S_SADS.F <- D_SADS.F[D_SADS.F$p.adj<0.05,]


# Contrast Test

SADS.F_3RP <- SADS.F[c(which(str_detect(SADS.F$Karyotype, "3RP|STD"))),]
SADS.F_2Lt <- SADS.F[c(which(str_detect(SADS.F$Karyotype, "2Lt|STD"))),]
SADS.F_3RK <- SADS.F[c(which(str_detect(SADS.F$Karyotype, "3RK|STD"))),]
SADS.F_Homoz <- SADS.F[-c(which(str_detect(SADS.F$Karyotype, "HET"))),]


a.18 = 1-(1-0.05)^(1/18)

All.comparisons.3RP.SADS.F <- dunn_test(SADS.F_3RP, F.Survival ~ Karyotype, p.adjust.method = "none")
All.comparisons.2Lt.SADS.F <- dunn_test(SADS.F_2Lt, F.Survival ~ Karyotype, p.adjust.method = "none")
All.comparisons.3RK.SADS.F <- dunn_test(SADS.F_3RK, F.Survival ~ Karyotype, p.adjust.method = "none")


Contrast.SADS.F <- rbind(All.comparisons.3RP.SADS.F, All.comparisons.2Lt.SADS.F, All.comparisons.3RK.SADS.F)

S_Contrast.SADS.F <- Contrast.SADS.F[Contrast.SADS.F$p< a.18,]
S_Contrast.SADS.F <- mutate(S_Contrast.SADS.F, p.signif =case_when(p < a.18 ~ "*", TRUE ~ "ns"))


# Add p-value

p_SADS.F_3RP <- S_Contrast.SADS.F[S_Contrast.SADS.F$group1 %in% c("STD","3RP","3RP_HET_1", "3RP_HET_2") & S_Contrast.SADS.F$group2 %in% c("STD","3RP","3RP_HET_1", "3RP_HET_2"),]
p_SADS.F_2Lt <- S_Contrast.SADS.F[S_Contrast.SADS.F$group1 %in% c("STD","2Lt","2Lt_HET_1", "2Lt_HET_2") & S_Contrast.SADS.F$group2 %in% c("STD","2Lt","2Lt_HET_1", "2Lt_HET_2"),]
p_SADS.F_3RK <- S_Contrast.SADS.F[S_Contrast.SADS.F$group1 %in% c("STD","3RK","3RK_HET_1", "3RK_HET_2") & S_Contrast.SADS.F$group2 %in% c("STD","3RK","3RK_HET_1", "3RK_HET_2"),]
p_SADS.F_Homoz <- S_Contrast.SADS.F[S_Contrast.SADS.F$group1 %in% c("STD","3RP","2Lt","3RK") & S_Contrast.SADS.F$group2 %in% c("STD","3RP","2Lt","3RK"),]


# Box plot

SADS.F_b1 <- ggplot(SADS.F_3RP, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %")+
  stat_pvalue_manual(p_SADS.F_3RP, label="p.signif", hide.ns = TRUE, y.position = 1.05,  step.increase = 0.08)

SADS.F_b2 <- ggplot(SADS.F_2Lt, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %")+
  stat_pvalue_manual(p_SADS.F_2Lt, label="p.signif", hide.ns = TRUE, y.position= 1.05, step.increase = 0.08)


SADS.F_b3 <- ggplot(SADS.F_3RK, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %")+
  stat_pvalue_manual(p_SADS.F_3RK, label="p.signif", hide.ns = TRUE, y.position = 1.05, step.increase = 0.08)


SADS.F_b4 <- ggplot(SADS.F_Homoz, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %")+
  stat_pvalue_manual(p_SADS.F_Homoz, label="p.signif", hide.ns = TRUE, y.position = 1.05, step.increase = 0.08)


Final <- (SADS.F_b1|SADS.F_b2)/(SADS.F_b3|SADS.F_b4) + plot_annotation(tag_levels='A',
                                                                       theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),
                                                                          plot.caption = element_text(hjust=0.5, size =12)))
                                            #+ plot_annotation(title="Female Survival after a Desiccation Stress", 



# Export plot as png file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/Stress_Resistance_Experiment/Analysis/2026_Stress_Resistance/Graphs _and_Tables/Word/Female_Desiccation_Stress.png",
  plot = Final,          
  width = 40,                    
  height = 20,                    
  units = "cm",                   
  dpi = 300                        
)








