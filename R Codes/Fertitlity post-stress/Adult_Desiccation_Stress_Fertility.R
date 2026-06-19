library(dplyr) 
library(ggplot2)
library(rstatix)
library(ggpubr) 
library(car) 
library(stringr)  
library(patchwork)

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


# Import Data

DataFADS <- read.csv2(file.choose())
FADS <- DataFADS

FADS$Fertility <- FADS$Adults + FADS$Pupae_.dark.dead.


# Assumption

M_FADS <- lm(Fertility ~ Karyotype, FADS)
shapiro.test(aov(M_FADS)$residuals)

par(mfrow=c(1,2))
hist(aov(M_FADS)$residuals,
     main="A. Histogramm of residuals",
     xlab="Residuals value",
     ylab="Frequency")
qqPlot(aov(M_FADS)$residuals, id=FALSE,
       main="B. Q-Q Plot",
       xlab="Theoretical Quantiles",
       ylab="Residuals value")

leveneTest(Fertility ~ Karyotype, FADS)


# Analysis -- Non-parametric Test

kruskal.test(Fertility ~ Karyotype, FADS)

# Contrast Test

FADS_3RP <- FADS[c(which(str_detect(FADS$Karyotype, "3RP|STD"))),]
FADS_2Lt <- FADS[c(which(str_detect(FADS$Karyotype, "2Lt|STD"))),]
FADS_3RK <- FADS[c(which(str_detect(FADS$Karyotype, "3RK|STD"))),]

a.18 = 1-(1-0.05)^(1/18)

All.comparisons.3RP.FADS <- dunn_test(FADS_3RP, Fertility ~ Karyotype, p.adjust.method = "none")
All.comparisons.2Lt.FADS <- dunn_test(FADS_2Lt, Fertility ~ Karyotype, p.adjust.method = "none")
All.comparisons.3RK.FADS <- dunn_test(FADS_3RK, Fertility ~ Karyotype, p.adjust.method = "none")

Contrast.FADS <- rbind(All.comparisons.3RP.FADS, All.comparisons.2Lt.FADS, All.comparisons.3RK.FADS)

S_Contrast.FADS <- Contrast.FADS[Contrast.FADS$p< a.18,]
S_Contrast.FADS <- mutate(S_Contrast.FADS, p.signif =case_when(p < a.18 ~ "*", TRUE ~ "ns"))   



# Add p-value

FADS_3RP <- FADS[c(which(str_detect(FADS$Karyotype, "3RP|STD"))),]
FADS_2Lt <- FADS[c(which(str_detect(FADS$Karyotype, "2Lt|STD"))),]
FADS_3RK <- FADS[c(which(str_detect(FADS$Karyotype, "3RK|STD"))),]
FADS_Homoz <- FADS[-c(which(str_detect(FADS$Karyotype, "HET"))),]

# S_Contrast.FADS <- separate(S_Contrast.FADS, contrast, into= c("group1", "group2"), sep = " - ")

p_FADS_3RP <- S_Contrast.FADS[S_Contrast.FADS$group1 %in% c("STD","3RP","3RP_HET_1","3RP_HET_2") & S_Contrast.FADS$group2 %in% c("STD","3RP","3RP_HET_1","3RP_HET_2"),]
# x %in% y checks if each element in x is present in y
p_FADS_2Lt <- S_Contrast.FADS[S_Contrast.FADS$group1 %in% c("STD","2Lt","2Lt_HET_1","2Lt_HET_2") & S_Contrast.FADS$group2 %in% c("STD","2Lt","2Lt_HET_1","2Lt_HET_2"),]
p_FADS_3RK <- S_Contrast.FADS[S_Contrast.FADS$group1 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2") & S_Contrast.FADS$group2 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2"),]
p_FADS_Homoz <- S_Contrast.FADS[S_Contrast.FADS$group1 %in% c("STD","3RP","2Lt","3RK") & S_Contrast.FADS$group2 %in% c("STD","3RP","2Lt","3RK"),]




# Box plot

FADS_b1 <- ggplot(FADS_3RP, aes(x = Karyotype, y = Fertility, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170) +
  plot_theme+
  ylab("Number of offspring")+
  labs(tag="A")+
  stat_pvalue_manual(p_FADS_3RP, label="p.signif", hide.ns = TRUE, y.position = 130,  step.increase = 0.08)

FADS_b2 <- ggplot(FADS_2Lt, aes(x = Karyotype, y = Fertility, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170) +
  plot_theme+
  ylab("Number of offspring")+
  labs(tag="B")+
  stat_pvalue_manual(p_FADS_2Lt, label="p.signif", hide.ns = TRUE, y.position= 130, step.increase = 0.08)


FADS_b3 <- ggplot(FADS_3RK, aes(x = Karyotype, y = Fertility, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170)+
  plot_theme+
  ylab("Number of offspring")+
  labs(tag="C")+
  stat_pvalue_manual(p_FADS_3RK, label="p.signif", hide.ns = TRUE, y.position = 150, step.increase = 0.08)


FADS_b4 <- ggplot(FADS_Homoz, aes(x = Karyotype, y = Fertility, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170)+
  plot_theme+
  ylab("Number of offspring")+
  labs(tag="D")+
  stat_pvalue_manual(p_FADS_Homoz, label="p.signif", hide.ns = TRUE, y.position = 150, step.increase = 0.08)


Final <- (FADS_b1|FADS_b2)/(FADS_b3|FADS_b4)



# Export plot as png file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/2026_D.melanogaster_Stress_Experiment/Plots/Adult_Desiccation_Stress_Fertility.png",
  plot = Final,          
  width = 45,                    
  height = 25,                    
  units = "cm",                   
  dpi = 300                        
)





