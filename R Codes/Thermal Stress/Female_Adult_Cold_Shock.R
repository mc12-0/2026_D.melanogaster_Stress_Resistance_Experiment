library(dplyr) 
library(ggplot2) 
library(rstatix) 
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

y_axis <- scale_y_continuous(
  limits = c(0.2, 1),
  breaks = c(0.2, 0.4, 0.6, 0.8,1)
)


# Import & Prepare Data

DataSACS.F <- read.csv2(file.choose())
SACS.F <- DataSACS.F

SACS.F$F.Survival <- SACS.F$ALIVE.F/SACS.F$TOTAL.F


# Assumptions

M_SACS.F <- lm(F.Survival ~ Karyotype, SACS.F)
shapiro.test(aov(M_SACS.F)$residuals)

par(mfrow=c(1,2))
hist(aov(M_SACS.F)$residuals,
     main="A. Histogramm of residuals",
     xlab="Residuals value",
     ylab="Frequency")

qqPlot(aov(M_SACS.F)$residuals, id=FALSE,
       main="B. Q-Q Plot",
       xlab="Theoretical Quantiles",
       ylab="Residuals value")

leveneTest(F.Survival ~ Karyotype, SACS.F)

# Statistical Analysis -- Non-parametric test

kruskal.test(F.Survival ~ Karyotype, SACS.F)
dunn_test(SACS.F, F.Survival ~ Karyotype)

D_SACS.F <- as.data.frame(dunn_test(SACS.F, F.Survival ~ Karyotype))
S_SACS.F <- D_SACS.F[D_SACS.F$p.adj<0.05,]

# Contrast Test

SACS.F_3RP <- SACS.F[c(which(str_detect(SACS.F$Karyotype, "3RP|STD"))),]
SACS.F_2Lt <- SACS.F[c(which(str_detect(SACS.F$Karyotype, "2Lt|STD"))),]
SACS.F_3RK <- SACS.F[c(which(str_detect(SACS.F$Karyotype, "3RK|STD"))),]
SACS.F_Homoz <- SACS.F[-c(which(str_detect(SACS.F$Karyotype, "HET"))),]


a.18 = 1-(1-0.05)^(1/18)

All.comparisons.3RP.SACS.F <- dunn_test(SACS.F_3RP, F.Survival ~ Karyotype, p.adjust.method = "none")
All.comparisons.2Lt.SACS.F <- dunn_test(SACS.F_2Lt, F.Survival ~ Karyotype, p.adjust.method = "none")
All.comparisons.3RK.SACS.F <- dunn_test(SACS.F_3RK, F.Survival ~ Karyotype, p.adjust.method = "none")

Contrast.SACS.F <- rbind(All.comparisons.3RP.SACS.F, All.comparisons.2Lt.SACS.F, All.comparisons.3RK.SACS.F)
S_Contrast.SACS.F <- Contrast.SACS.F[Contrast.SACS.F$p< a.18,]
S_Contrast.SACS.F <- mutate(S_Contrast.SACS.F, p.signif =case_when(p< a.18 ~ "*", TRUE ~ "ns"))


# Box plot

SACS.F_b1 <- ggplot(SACS.F_3RP, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %" )+
  labs(tag="A")


SACS.F_b2 <- ggplot(SACS.F_2Lt, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %" )+
  labs(tag="B")


SACS.F_b3 <- ggplot(SACS.F_3RK, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %" )+
  labs(tag="C")



SACS.F_b4 <- ggplot(SACS.F_Homoz, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %" )+
  labs(tag="D")



Final <- (SACS.F_b1|SACS.F_b2)/(SACS.F_b3|SACS.F_b4)



# Export plot as png file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/2026_D.melanogaster_Stress_Experiment/Plots/Female_Adult_Cold_Shock.png",
  plot = Final,          
  width = 45,                    
  height = 25,                    
  units = "cm",                   
  dpi = 300                        
)

