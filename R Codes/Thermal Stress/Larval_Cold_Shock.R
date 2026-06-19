library(ggpubr) 
library(dplyr) 
library(ggplot2)
library(rstatix)
library(car) 
library(stringr)  
library(patchwork)
library(svglite)


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
        axis.title.x = element_text(face = "bold", size = 22),
        axis.title.y = element_text(face = "bold", size = 22),
        axis.text.x = element_text(face = "bold", size = 16),
        axis.text.y = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 20),
        legend.text = element_text(face = "bold", size = 18),
        plot.tag = element_text(size = 20, face="bold")
  )

y_axis <- scale_y_continuous(
  limits = c(0, 0.5),
  breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
)


# Import Data

DataSLCS <- read.csv2(file.choose())
SLCS <- DataSLCS


SLCS$Total.survival <- SLCS$Adults +SLCS$Pupae
SLCS$Survival.rate <- SLCS$Total.survival/SLCS$Total_eggs

# Assumptions

M_Survival_Rate <- lm(Survival.rate~Karyotype, SLCS)
shapiro.test(M_Survival_Rate$residuals)

par(mfrow=c(1,2))
hist(aov(M_Survival_Rate)$residuals,
     main="A. Histogramm of residuals",
     xlab="Residuals value",
     ylab="Frequency")

qqPlot(aov(M_Survival_Rate)$residuals, id=FALSE,
       main="B. Q-Q Plot",
       xlab="Theoretical Quantiles",
       ylab="Residuals value")

leveneTest(Survival.rate~Karyotype, SLCS)


# Analysis -- Non-parametric test

kruskal.test(Survival.rate ~ Karyotype, SLCS)

# Contrast Test

SLCS_3RP <- SLCS[c(which(str_detect(SLCS$Karyotype, "3RP|STD"))),]
SLCS_2Lt <- SLCS[c(which(str_detect(SLCS$Karyotype, "2Lt|STD"))),]
SLCS_3RK <- SLCS[c(which(str_detect(SLCS$Karyotype, "3RK|STD"))),]


a.18 = 1-(1-0.05)^(1/18)

All.comparisons.3RP.SLCS <- dunn_test(SLCS_3RP, Survival.rate ~ Karyotype, p.adjust.method = "none")
All.comparisons.2Lt.SLCS <- dunn_test(SLCS_2Lt, Survival.rate ~ Karyotype, p.adjust.method = "none")
All.comparisons.3RK.SLCS <- dunn_test(SLCS_3RK, Survival.rate ~ Karyotype, p.adjust.method = "none")

Contrast.SLCS <- rbind(All.comparisons.3RP.SLCS, All.comparisons.2Lt.SLCS, All.comparisons.3RK.SLCS)
S_Contrast.SLCS <- Contrast.SLCS[Contrast.SLCS$p< a.18,]
S_Contrast.SLCS <- mutate(S_Contrast.SLCS, p.signif =case_when(p< a.18  ~ "*", TRUE ~ "ns"))

# Add p-value

SLCS_3RP <- SLCS[c(which(str_detect(SLCS$Karyotype, "3RP|STD"))),]
SLCS_2Lt <- SLCS[c(which(str_detect(SLCS$Karyotype, "2Lt|STD"))),]
SLCS_3RK <- SLCS[c(which(str_detect(SLCS$Karyotype, "3RK|STD"))),]
SLCS_Homoz <- SLCS[-c(which(str_detect(SLCS$Karyotype, "HET"))),]


p_SLCS_3RP <- S_Contrast.SLCS[S_Contrast.SLCS$group1 %in% c("STD","3RP","3RP_HET_1", "3RP_HET_2") & S_Contrast.SLCS$group2 %in% c("STD","3RP","3RP_HET_1", "3RP_HET_2"),]
p_SLCS_2Lt <- S_Contrast.SLCS[S_Contrast.SLCS$group1 %in% c("STD","2Lt","2Lt_HET_1", "2Lt_HET_2") & S_Contrast.SLCS$group2 %in% c("STD","2Lt","2Lt_HET_1", "2Lt_HET_2"),]
p_SLCS_3RK <- S_Contrast.SLCS[S_Contrast.SLCS$group1 %in% c("STD","3RK","3RK_HET_1", "3RK_HET_2") & S_Contrast.SLCS$group2 %in% c("STD","3RK","3RK_HET_1", "3RK_HET_2"),]
p_SLCS_Homoz <- S_Contrast.SLCS[S_Contrast.SLCS$group1 %in% c("STD","3RP","2Lt","3RK") & S_Contrast.SLCS$group2 %in% c("STD","3RP","2Lt","3RK"),]


# Box plot 

SLCS_b1 <- ggplot(SLCS_3RP, aes(x = Karyotype, y = Survival.rate, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Survival rate in %" )+
  labs(tag="A")+
  stat_pvalue_manual(p_SLCS_3RP, label="p.signif", hide.ns = TRUE, y.position = 0.40, step.increase = 0.08)


SLCS_b2 <- ggplot(SLCS_2Lt, aes(x = Karyotype, y = Survival.rate, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  theme_bw()+
  plot_theme+
  ylab("Survival rate in %" )+
  labs(tag="B")+
  stat_pvalue_manual(p_SLCS_2Lt, label="p.signif", hide.ns = TRUE, y.position=0.35, step.increase = 0.08)


SLCS_b3 <- ggplot(SLCS_3RK, aes(x = Karyotype, y = Survival.rate, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Survival rate in %" )+
  labs(tag="C")+
  stat_pvalue_manual(p_SLCS_3RK, label="p.signif", hide.ns = TRUE, y.position = 0.35,  step.increase = 0.08)


SLCS_b4 <- ggplot(SLCS_Homoz, aes(x = Karyotype, y = Survival.rate, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Survival rate in %" )+
  labs(tag="D")+
  stat_pvalue_manual(p_SLCS_Homoz, label="p.signif", hide.ns = TRUE, y.position = 0.40, step.increase = 0.08)


Final <- (SLCS_b1|SLCS_b2)/(SLCS_b3|SLCS_b4) 


# Export plot as svg file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/2026_D.melanogaster_Stress_Experiment/Plots/Larval_Cold_Shock.svg",
  plot = Final,          
  width = 45,                    
  height = 25,                    
  units = "cm",                   
  dpi = 300                        
)





