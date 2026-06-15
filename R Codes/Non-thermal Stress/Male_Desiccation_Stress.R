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
  limits = c(0, 1),
  breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
)


# Import Data
DataSADS.M <- read.csv2(file.choose())
SADS.M <- DataSADS.M
SADS.M <- DataSADS.M[-c(which(is.na(DataSADS.M$Total.Adults))),]

SADS.M$M.Survival <- SADS.M$ALIVE.M.after.exp/SADS.M$TOTAL.M


# Assumption

M_SADS.M <- lm(M.Survival ~ Karyotype, SADS.M)

shapiro.test(aov(M_SADS.M)$residuals)

par(mfrow=c(1,2))
hist(aov(M_SADS.M)$residuals)
qqPlot(aov(M_SADS.M)$residuals, id=FALSE)

leveneTest(M.Survival ~ Karyotype, SADS.M)

# Analysis -- Non-parametric test

kruskal.test(M.Survival ~ Karyotype, SADS.M)
dunn_test(SADS.M, M.Survival ~ Karyotype)

D_SADS.M <- as.data.frame(dunn_test(SADS.M, M.Survival ~ Karyotype))
S_SADS.M <- D_SADS.M[D_SADS.M$p.adj<0.05,]


# Contrast Test

SADS.M_3RP <- SADS.M[c(which(str_detect(SADS.M$Karyotype, "3RP|STD"))),]
SADS.M_2Lt <- SADS.M[c(which(str_detect(SADS.M$Karyotype, "2Lt|STD"))),]
SADS.M_3RK <- SADS.M[c(which(str_detect(SADS.M$Karyotype, "3RK|STD"))),]
SADS.M_Homoz <- SADS.M[-c(which(str_detect(SADS.M$Karyotype, "HET"))),]

a.18 = 1-(1-0.05)^(1/18)

All.comparisons.3RP.SADS.M <- dunn_test(SADS.M_3RP, M.Survival ~ Karyotype, p.adjust.method = "none")
All.comparisons.2Lt.SADS.M <- dunn_test(SADS.M_2Lt, M.Survival ~ Karyotype, p.adjust.method = "none")
All.comparisons.3RK.SADS.M <- dunn_test(SADS.M_3RK, M.Survival ~ Karyotype, p.adjust.method = "none")

Contrast.SADS.M <- rbind(All.comparisons.3RP.SADS.M, All.comparisons.2Lt.SADS.M, All.comparisons.3RK.SADS.M)

S_Contrast.SADS.M <- Contrast.SADS.M[Contrast.SADS.M$p< a.18,]
S_Contrast.SADS.M <- mutate(S_Contrast.SADS.M, p.signif =case_when(p < a.18 ~ "*", TRUE ~ "ns"))  


# Add p-value

p_SADS.M_3RP <- S_Contrast.SADS.M[S_Contrast.SADS.M$group1 %in% c("STD","3RP","3RP_HET_1", "3RP_HET_2") & S_Contrast.SADS.M$group2 %in% c("STD","3RP","3RP_HET_1", "3RP_HET_2"),]
p_SADS.M_2Lt <- S_Contrast.SADS.M[S_Contrast.SADS.M$group1 %in% c("STD","2Lt","2Lt_HET_1", "2Lt_HET_2") & S_Contrast.SADS.M$group2 %in% c("STD","2Lt","2Lt_HET_1", "2Lt_HET_2"),]
p_SADS.M_3RK <- S_Contrast.SADS.M[S_Contrast.SADS.M$group1 %in% c("STD","3RK","3RK_HET_1", "3RK_HET_2") & S_Contrast.SADS.M$group2 %in% c("STD","3RK","3RK_HET_1", "3RK_HET_2"),]
p_SADS.M_Homoz <- S_Contrast.SADS.M[S_Contrast.SADS.M$group1 %in% c("STD","3RP","2Lt","3RK") & S_Contrast.SADS.M$group2 %in% c("STD","3RP","2Lt","3RK"),]

# Box plot

SADS.M_b1 <- ggplot(SADS.M_3RP, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")
#stat_pvalue_manual(p_SADS.M_3RP, label="p.signif", hide.ns = TRUE, y.position = 1.05,  step.increase = 0.08)

SADS.M_b2 <- ggplot(SADS.M_2Lt, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")
#stat_pvalue_manual(p_SADS.M_2Lt, label="p.signif", hide.ns = TRUE, y.position= 1.05, step.increase = 0.08)


SADS.M_b3 <- ggplot(SADS.M_3RK, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")+
  stat_pvalue_manual(p_SADS.M_3RK, label="p.signif", hide.ns = TRUE, y.position = 0.9, step.increase = 0.08)


SADS.M_b4 <- ggplot(SADS.M_Homoz, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")+
  stat_pvalue_manual(p_SADS.M_Homoz, label="p.signif", hide.ns = TRUE, y.position = 0.9, step.increase = 0.08)


Final <- (SADS.M_b1|SADS.M_b2)/(SADS.M_b3|SADS.M_b4) + plot_annotation(tag_levels='A',
                                                                       theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),
                                                                          plot.caption = element_text(hjust=0.5, size =12)))
                                                      #plot_annotation(title="Male Survival after a Desiccation Stress"


# Export plot as png file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/Stress_Resistance_Experiment/Analysis/2026_Stress_Resistance/Graphs _and_Tables/Word/Male_Desiccation_Stress.png",
  plot = Final,          
  width = 40,                    
  height = 20,                    
  units = "cm",                   
  dpi = 300                        
)
