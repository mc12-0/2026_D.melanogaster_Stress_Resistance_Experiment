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

sadsfsd
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
  limits = c(0.5, 1.20),
  breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
)


# Import Data

DataSAHS.F <- read.csv2(file.choose())
SAHS.F <- DataSAHS.F

SAHS.F$F.Survival <- SAHS.F$ALIVE.F/SAHS.F$TOTAL.F
SAHS.F <- mutate(SAHS.F, Borders=case_when(str_detect(Remarks, "BORD")~ "YES",
                                           TRUE ~ "NO"))

boxplot(F.Survival~Borders, SAHS.F)
t.test(F.Survival ~ Borders, SAHS.F, var.equal=TRUE)
SAHS.F <- SAHS.F[SAHS.F$Borders == "NO",]


# Assumption

M_SAHS.F <- lm(F.Survival ~ Karyotype, SAHS.F)

shapiro.test(aov(M_SAHS.F)$residuals)

par(mfrow=c(1,2))
hist(aov(M_SAHS.F)$residuals)
qqPlot(aov(M_SAHS.F)$residuals, id=FALSE)

leveneTest(F.Survival ~ Karyotype, SAHS.F)


# Statistical analysis -- Non-parametric test

M_SAHS.F <- lm(F.Survival ~ Karyotype, SAHS.F)

kruskal.test(F.Survival ~ Karyotype, SAHS.F)
dunn_test(SAHS.F, F.Survival ~ Karyotype)

D_SAHS.F <- as.data.frame(dunn_test(SAHS.F, F.Survival ~ Karyotype))
S_SAHS.F <- D_SAHS.F[D_SAHS.F$p.adj<0.05,]

# Contrast Test

SAHS.F_3RP <- SAHS.F[c(which(str_detect(SAHS.F$Karyotype, "3RP|STD"))),]
SAHS.F_2Lt <- SAHS.F[c(which(str_detect(SAHS.F$Karyotype, "2Lt|STD"))),]
SAHS.F_3RK <- SAHS.F[c(which(str_detect(SAHS.F$Karyotype, "3RK|STD"))),]

a.18 = 1-(1-0.05)^(1/18)

All.comparisons.3RP.SAHS.F <- dunn_test(SAHS.F_3RP, F.Survival ~ Karyotype, p.adjust.method = "none")
All.comparisons.2Lt.SAHS.F <- dunn_test(SAHS.F_2Lt, F.Survival ~ Karyotype, p.adjust.method = "none")
All.comparisons.3RK.SAHS.F <- dunn_test(SAHS.F_3RK, F.Survival ~ Karyotype, p.adjust.method = "none")

Contrast.SAHS.F <- rbind(All.comparisons.3RP.SAHS.F, All.comparisons.2Lt.SAHS.F, All.comparisons.3RK.SAHS.F)
S_Contrast.SAHS.F <- Contrast.SAHS.F[Contrast.SAHS.F$p< a.18,]
S_Contrast.SAHS.F <- mutate(S_Contrast.SAHS.F, p.signif =case_when(p< a.18 ~ "*", TRUE ~ "ns"))


# Add p-value

SAHS.F_3RP <- SAHS.F[c(which(str_detect(SAHS.F$Karyotype, "3RP|STD"))),]
SAHS.F_2Lt <- SAHS.F[c(which(str_detect(SAHS.F$Karyotype, "2Lt|STD"))),]
SAHS.F_3RK <- SAHS.F[c(which(str_detect(SAHS.F$Karyotype, "3RK|STD"))),]
SAHS.F_Homoz <- SAHS.F[-c(which(str_detect(SAHS.F$Karyotype, "HET"))),]


p_SAHS.F_3RP <- S_Contrast.SAHS.F[S_Contrast.SAHS.F$group1 %in% c("STD","3RP","3RP_HET_1", "3RP_HET_2") & S_Contrast.SAHS.F$group2 %in% c("STD","3RP","3RP_HET_1","3RP_HET_2"),]
p_SAHS.F_2Lt <- S_Contrast.SAHS.F[S_Contrast.SAHS.F$group1 %in% c("STD","2Lt","2Lt_HET_1", "2Lt_HET_2") & S_Contrast.SAHS.F$group2 %in% c("STD","2Lt","2Lt_HET_1", "2Lt_HET_2"),]
p_SAHS.F_3RK <- S_Contrast.SAHS.F[S_Contrast.SAHS.F$group1 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2") & S_Contrast.SAHS.F$group2 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2"),]
p_SAHS.F_Homoz <- S_Contrast.SAHS.F[S_Contrast.SAHS.F$group1 %in% c("STD","3RP","2Lt","3RK") & S_Contrast.SAHS.F$group2 %in% c("STD","3RP","2Lt","3RK"),]


# Box plot

SAHS.F_b1 <- ggplot(SAHS.F_3RP, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %" )+
  stat_pvalue_manual(p_SAHS.F_3RP, label="p.signif", hide.ns = TRUE, y.position = 1.05, step.increase = 0.08)


SAHS.F_b2 <- ggplot(SAHS.F_2Lt, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+ 
  y_axis+
  plot_theme+
  ylab("Female survival rate in %" )+
  stat_pvalue_manual(p_SAHS.F_2Lt, label="p.signif", hide.ns = TRUE, y.position=1.05, step.increase = 0.095)


SAHS.F_b3 <- ggplot(SAHS.F_3RK, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+  
  y_axis+
  plot_theme+
  ylab("Female survival rate in %" )+
  stat_pvalue_manual(p_SAHS.F_3RK, label="p.signif", hide.ns = TRUE, y.position = 1.05,  step.increase = 0.08)


SAHS.F_b4 <- ggplot(SAHS.F_Homoz, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %" )+
  stat_pvalue_manual(p_SAHS.F_Homoz, label="p.signif", hide.ns = TRUE, y.position = 1.05, step.increase = 0.08)


Final <- (SAHS.F_b1|SAHS.F_b2)/(SAHS.F_b3|SAHS.F_b4) + plot_annotation(tag_levels = "A", 
                                                               theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),plot.caption = element_text(hjust=0.5, size =12)))
                                            #plot_annotation(title="Female Survival after an Adult Heat Shock")
                                                              

# Export plot as png file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/Stress_Resistance_Experiment/Analysis/2026_Stress_Resistance/Graphs _and_Tables/Word/Female_Adult_Heat_Shock.png",
  plot = Final,          
  width = 40,                    
  height = 20,                    
  units = "cm",                   
  dpi = 300                        
)
                                                              
                                                              