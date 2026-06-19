library(ggpubr)
library(dplyr) 
library(ggplot2)
library(rstatix) 
library(car)
library(stringr)
library(patchwork)
library(emmeans)
library(tidyr)



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
  limits = c(0, 1.1),
  breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
)


# Import & Prepare Data

DataSACS.M <- read.csv2(file.choose())
SACS.M <- DataSACS.M

SACS.M$M.Survival <- SACS.M$ALIVE_M/SACS.M$TOTAL_M


# Assumptions

M_SACS.M <- lm(M.Survival ~ Karyotype, SACS.M)

shapiro.test(aov(M_SACS.M)$residuals)

par(mfrow=c(1,2))
hist(aov(M_SACS.M)$residuals,
     main="A. Histogramm of residuals",
     xlab="Residuals value",
     ylab="Frequency")

qqPlot(aov(M_SACS.M)$residuals, id=FALSE,
       main="B. Q-Q Plot",
       xlab="Theoretical Quantiles",
       ylab="Residuals value")

leveneTest(M.Survival ~ Karyotype, SACS.M)

# Statistical Analysis -- ANOVA

anova(M_SACS.M)

# Contrast Test

SACS.M_3RP <- SACS.M[c(which(str_detect(SACS.M$Karyotype, "3RP|STD"))),]
SACS.M_2Lt <- SACS.M[c(which(str_detect(SACS.M$Karyotype, "2Lt|STD"))),]
SACS.M_3RK <- SACS.M[c(which(str_detect(SACS.M$Karyotype, "3RK|STD"))),]
SACS.M_Homoz <- SACS.M[-c(which(str_detect(SACS.M$Karyotype, "HET"))),]

a.18 = 1-(1-0.05)^(1/18)

M_SACS.M.3RP <- lm(M.Survival ~ Karyotype, SACS.M_3RP)
means.3RP.SACS.M <- emmeans(M_SACS.M.3RP, "Karyotype")
All.comparisons.3RP.SACS.M <- pairs(means.3RP.SACS.M, adjust="none") 
All.comparisons.3RP.SACS.M <- as.data.frame(All.comparisons.3RP.SACS.M)

M_SACS.M.2Lt <- lm(M.Survival ~ Karyotype, SACS.M_2Lt)
means.2Lt.SACS.M <- emmeans(M_SACS.M.2Lt, "Karyotype")
All.comparisons.2Lt.SACS.M <- pairs(means.2Lt.SACS.M, adjust="none")
All.comparisons.2Lt.SACS.M <- as.data.frame(All.comparisons.2Lt.SACS.M)

M_SACS.M.3RK <- lm(M.Survival ~ Karyotype, SACS.M_3RK)
means.3RK.SACS.M <- emmeans(M_SACS.M.3RK, "Karyotype")
All.comparisons.3RK.SACS.M <- pairs(means.3RK.SACS.M, adjust="none")
All.comparisons.3RK.SACS.M <- as.data.frame(All.comparisons.3RK.SACS.M)


Contrast.SACS.M <- rbind(All.comparisons.3RP.SACS.M, All.comparisons.2Lt.SACS.M, All.comparisons.3RK.SACS.M)

S_Contrast.SACS.M <- Contrast.SACS.M[Contrast.SACS.M$p.value < a.18,]
S_Contrast.SACS.M <- mutate(S_Contrast.SACS.M, p.signif =case_when(p.value< a.18 ~ "*", TRUE ~ "ns"))


# Add p-value
S_Contrast.SACS.M <- separate(S_Contrast.SACS.M, contrast, into= c("group1", "group2"), sep = " - ")

p_SACS.M_3RP <- S_Contrast.SACS.M[S_Contrast.SACS.M$group1 %in% c("STD","3RP","3RP_HET_1","3RP_HET_2") & S_Contrast.SACS.M$group2 %in% c("STD","3RP","3RP_HET_1","3RP_HET_2"),]
p_SACS.M_2Lt <- S_Contrast.SACS.M[S_Contrast.SACS.M$group1 %in% c("STD","2Lt","2Lt_HET_1","2Lt_HET_2") & S_Contrast.SACS.M$group2 %in% c("STD","2Lt","2Lt_HET_1","2Lt_HET_2"),]
p_SACS.M_3RK <- S_Contrast.SACS.M[S_Contrast.SACS.M$group1 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2") & S_Contrast.SACS.M$group2 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2"),]
p_SACS.M_Homoz <- S_Contrast.SACS.M[S_Contrast.SACS.M$group1 %in% c("STD","3RP","2Lt","3RK") & S_Contrast.SACS.M$group2 %in% c("STD","3RP","2Lt","3RK"),]


# Box plot

SACS.M_b1 <- ggplot(SACS.M_3RP, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")+
  labs(tag="A")+
  stat_pvalue_manual(p_SACS.M_3RP, label="p.signif", hide.ns = TRUE, y.position = 0.85, step.increase = 0.08)

SACS.M_b2 <- ggplot(SACS.M_2Lt, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")+
  labs(tag="B")+
  stat_pvalue_manual(p_SACS.M_2Lt, label="p.signif", hide.ns = TRUE, y.position=0.97, step.increase = 0.08)


SACS.M_b3 <- ggplot(SACS.M_3RK, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")+
  labs(tag="C")


SACS.M_b4 <- ggplot(SACS.M_Homoz, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")+
  labs(tag="D")+
  stat_pvalue_manual(p_SACS.M_Homoz, label="p.signif", hide.ns = TRUE, y.position = 0.95, step.increase = 0.08)


Final <- (SACS.M_b1|SACS.M_b2)/(SACS.M_b3|SACS.M_b4)

# Export plot as png file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/2026_D.melanogaster_Stress_Experiment/Plots/Male_Adult_Cold_Shock.png",
  plot = Final,          
  width = 45,                    
  height = 25,                    
  units = "cm",                   
  dpi = 300                        
)


