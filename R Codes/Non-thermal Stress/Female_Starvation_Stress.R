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
  limits = c(0, 1),
  breaks = c(0.2, 0.4, 0.6, 0.8,1)
)


# Import & Prepare Data

DataSASS.F <-read.csv2(file.choose())
SASS.F <- DataSASS.F
SASS.F$F.Survival <- SASS.F$ALIVE_F/SASS.F$Total_F

boxplot(F.Survival ~ GAS, SASS.F)
t.test(F.Survival ~ GAS, SASS.F, var.equal=TRUE)

boxplot(F.Survival ~ BORDERS, SASS.F)
t.test(F.Survival ~ BORDERS, SASS.F, var.equal=TRUE)

SASS.F <- SASS.F[SASS.F$BORDERS == "NO",]
SASS.F<- SASS.F[SASS.F$GAS == "YES",]


# Assumptions

M_SASS.F <- lm(F.Survival ~ Karyotype, SASS.F)

shapiro.test(aov(M_SASS.F)$residuals)

par(mfrow=c(1,2))
hist(aov(M_SASS.F)$residuals,
     main="A. Histogramm of residuals",
     xlab="Residuals value",
     ylab="Frequency")

qqPlot(aov(M_SASS.F)$residuals, id=FALSE,
       main="B. Q-Q Plot",
       xlab="Theoretical Quantiles",
       ylab="Residuals value")

leveneTest(F.Survival ~ Karyotype, SASS.F)


# Statistical Analysis -- Non-parametric Test

kruskal.test(F.Survival ~ Karyotype, SASS.F)

# Contrast Test

SASS.F_3RP <- SASS.F[c(which(str_detect(SASS.F$Karyotype, "3RP|STD"))),]
SASS.F_2Lt <- SASS.F[c(which(str_detect(SASS.F$Karyotype, "2Lt|STD"))),]
SASS.F_3RK <- SASS.F[c(which(str_detect(SASS.F$Karyotype, "3RK|STD"))),]
SASS.F_Homoz <- SASS.F[-c(which(str_detect(SASS.F$Karyotype, "HET"))),]

# Contrast with non-parametric test
a.18 = 1-(1-0.05)^(1/18)

All.comparisons.3RP.SASS.F <- dunn_test(SASS.F_3RP, F.Survival ~ Karyotype, p.adjust.method = "none")
All.comparisons.2Lt.SASS.F <- dunn_test(SASS.F_2Lt, F.Survival ~ Karyotype, p.adjust.method = "none")
All.comparisons.3RK.SASS.F <- dunn_test(SASS.F_3RK, F.Survival ~ Karyotype, p.adjust.method = "none")

Contrast.SASS.F <- rbind(All.comparisons.3RP.SASS.F, All.comparisons.2Lt.SASS.F, All.comparisons.3RK.SASS.F)

S_Contrast.SASS.F <- Contrast.SASS.F[Contrast.SASS.F$p< a.18,]
S_Contrast.SASS.F <- mutate(S_Contrast.SASS.F, p.signif =case_when(p < a.18 ~ "*", TRUE ~ "ns")) 


# Box plot
SASS.F_b1 <- ggplot(SASS.F_3RP, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %")+
  labs(tag="A")

SASS.F_b2 <- ggplot(SASS.F_2Lt, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %")+
  labs(tag="B")


SASS.F_b3 <- ggplot(SASS.F_3RK, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %")+
  labs(tag="C")

SASS.F_b4 <- ggplot(SASS.F_Homoz, aes(x = Karyotype, y = F.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Female survival rate in %")+
  labs(tag="D")


Final <- (SASS.F_b1|SASS.F_b2)/(SASS.F_b3|SASS.F_b4)  


# Export plot as svg file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/2026_D.melanogaster_Stress_Experiment/Plots/Female_Starvation_Stress.svg",
  plot = Final,          
  width = 45,                    
  height = 25,                    
  units = "cm",                   
  dpi = 300                        
)
