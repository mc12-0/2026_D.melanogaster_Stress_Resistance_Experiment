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
  limits = c(0, 0.5),
  breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
)


# Import & Prepare Data

DataSAHS.M <- read.csv2(file.choose())
SAHS.M <- DataSAHS.M


SAHS.M$M.Survival <- SAHS.M$ALIVE.M/SAHS.M$TOTAL.M
SAHS.M <- mutate(SAHS.M, Borders=case_when(str_detect(Remarks, "BORD")~ "YES",
                                           TRUE ~ "NO"))

boxplot(M.Survival~Borders, SAHS.M)
t.test(M.Survival ~ Borders, SAHS.M, var.equal=TRUE)


# Assumptions

M_SAHS.M <- lm(M.Survival ~ Karyotype, SAHS.M)

shapiro.test(aov(M_SAHS.M)$residuals)

par(mfrow=c(1,2))
hist(aov(M_SAHS.M)$residuals)
qqPlot(aov(M_SAHS.M)$residuals, id=FALSE)

leveneTest(M.Survival ~ Karyotype, SAHS.M)

# Data Analysis -- Non-parametric test

kruskal.test(M.Survival ~ Karyotype, SAHS.M)
dunn_test(SAHS.M, M.Survival ~ Karyotype)

D_SAHS.M <- as.data.frame(dunn_test(SAHS.M, M.Survival ~ Karyotype))
S_SAHS.M <- D_SAHS.M[D_SAHS.M$p.adj<0.05,]


# Add p-value

SAHS.M_3RP <- SAHS.M[c(which(str_detect(SAHS.M$Karyotype, "3RP|STD"))),]
SAHS.M_2Lt <- SAHS.M[c(which(str_detect(SAHS.M$Karyotype, "2Lt|STD"))),]
SAHS.M_3RK <- SAHS.M[c(which(str_detect(SAHS.M$Karyotype, "3RK|STD"))),]
SAHS.M_Homoz <- SAHS.M[-c(which(str_detect(SAHS.M$Karyotype, "HET"))),]

# Box plots

SAHS.M_b1 <- ggplot(SAHS.M_3RP, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")


SAHS.M_b2 <- ggplot(SAHS.M_2Lt, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")

SAHS.M_b3 <- ggplot(SAHS.M_3RK, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")


SAHS.M_b4 <- ggplot(SAHS.M_Homoz, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")

Final <-(SAHS.M_b1|SAHS.M_b2)/(SAHS.M_b3|SAHS.M_b4) + plot_annotation(tag_levels='A',
                                                                      theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),
                                                                          plot.caption = element_text(hjust=0.5, size =12)))
                                                    #+ plot_annotation(title="Male Survival after an Adult Heat Shock", 

# Export plot as png file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/Stress_Resistance_Experiment/Analysis/2026_Stress_Resistance/Graphs _and_Tables/Word/Male_Adult_Heat_Shock.png",
  plot = Final,          
  width = 40,                    
  height = 20,                    
  units = "cm",                   
  dpi = 300                        
)





