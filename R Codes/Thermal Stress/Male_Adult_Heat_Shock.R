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
  limits = c(0, 0.5),
  breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
)


# Import & Prepare Data

DataSAHS.M <- read.csv2(file.choose())
SAHS.M <- DataSAHS.M


SAHS.M$M.Survival <- SAHS.M$ALIVE.M/SAHS.M$TOTAL.M


boxplot(M.Survival~BORDERS, SAHS.M)
t.test(M.Survival ~ BORDERS, SAHS.M, var.equal=TRUE)


# Assumptions

M_SAHS.M <- lm(M.Survival ~ Karyotype, SAHS.M)

shapiro.test(aov(M_SAHS.M)$residuals)

par(mfrow=c(1,2))

hist(aov(M_SAHS.M)$residuals,
     main="A. Histogramm of residuals",
     xlab="Residuals value",
     ylab="Frequency")

qqPlot(aov(M_SAHS.M)$residuals, id=FALSE,
       main="B. Q-Q Plot",
       xlab="Theoretical Quantiles",
       ylab="Residuals value")

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
  ylab("Male survival rate in %")+
  labs(tag="A")


SAHS.M_b2 <- ggplot(SAHS.M_2Lt, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")+
  labs(tag="B")

SAHS.M_b3 <- ggplot(SAHS.M_3RK, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")+
  labs(tag="C")


SAHS.M_b4 <- ggplot(SAHS.M_Homoz, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")+
  labs(tag="D")

Final <-(SAHS.M_b1|SAHS.M_b2)/(SAHS.M_b3|SAHS.M_b4) 

# Export plot as png file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/2026_D.melanogaster_Stress_Experiment/Plots/Male_Adult_Heat_Shock.png",
  plot = Final,          
  width = 45,                    
  height = 25,                    
  units = "cm",                   
  dpi = 300                        
)





