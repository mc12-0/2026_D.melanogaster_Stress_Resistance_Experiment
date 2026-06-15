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
  limits = c(0, 1.20),
  breaks = c(0.2, 0.4, 0.6, 0.8,1)
)

# Import & Prepare Data
DataFASS <- read.csv2(file.choose())
FASS <- DataFASS[-c(which(is.na(DataFASS$Adults))),] # 379 entry
FASS$Fertility <- FASS$Adults + FASS$Pupae
FASS$Number.of.adults <-FASS$Fertility

# Assumption

M_FASS <- lm(Fertility ~ Karyotype, FASS)

shapiro.test(aov(M_FASS)$residuals)

par(mfrow=c(1,2))
hist(aov(M_FASS)$residuals)
qqPlot(aov(M_FASS)$residuals, id=FALSE)

leveneTest(Fertility ~ Karyotype, FASS)


# Statistical Analysis -- Non-parametric Test

kruskal.test(Fertility ~ Karyotype, FASS)
dunn_test(FASS, Fertility ~ Karyotype)

D_FASS <- as.data.frame(dunn_test(FASS, Fertility ~ Karyotype))
S_FASS <- D_FASS[D_FASS$p.adj<0.05,]


# Box plot

FASS_3RP <- FASS[c(which(str_detect(FASS$Karyotype, "3RP|STD"))),]
FASS_2Lt <- FASS[c(which(str_detect(FASS$Karyotype, "2Lt|STD"))),]
FASS_3RK <- FASS[c(which(str_detect(FASS$Karyotype, "3RK|STD"))),]
FASS_Homoz <- FASS[-c(which(str_detect(FASS$Karyotype, "HET"))),]


FASS_b1 <- ggplot(FASS_3RP, aes(x = Karyotype, y = Number.of.adults, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170) +
  plot_theme+
  ylab("Number of offspring")
  
FASS_b2 <- ggplot(FASS_2Lt, aes(x = Karyotype, y = Number.of.adults, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170) +
  plot_theme+
  ylab("Number of offspring")
  

FASS_b3 <- ggplot(FASS_3RK, aes(x = Karyotype, y = Number.of.adults, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170) +
  plot_theme+
  ylab("Number of offspring")
  
FASS_b4 <- ggplot(FASS_Homoz, aes(x = Karyotype, y = Number.of.adults, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170) +
  plot_theme+
  ylab("Number of offspring")
  
Final <- (FASS_b1|FASS_b2)/(FASS_b3|FASS_b4) + plot_annotation(tag_levels='A', 
                                                      theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),plot.caption = element_text(hjust=0.5, size =12)))
                                    #plot_annotation(title="Fecundity after a starvation shock"


# Export plot as png file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/Stress_Resistance_Experiment/Analysis/2026_Stress_Resistance/Graphs _and_Tables/Word/Adult_Starvation_Stress_Fertility.png",
  plot = Final,          
  width = 40,                    
  height = 20,                    
  units = "cm",                   
  dpi = 300                        
)



