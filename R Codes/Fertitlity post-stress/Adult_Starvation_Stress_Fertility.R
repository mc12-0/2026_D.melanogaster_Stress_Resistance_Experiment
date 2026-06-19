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
  limits = c(0, 1.20),
  breaks = c(0.2, 0.4, 0.6, 0.8,1)
)

# Import & Prepare Data
DataFASS <- read.csv2(file.choose())
FASS <- DataFASS

FASS$Fertility <- FASS$Adults + FASS$Pupae

boxplot(Fertility ~ Gas, FASS)
t.test(Fertility ~ Gas, FASS, var.equal=TRUE)


# Assumption

M_FASS <- lm(Fertility ~ Karyotype, FASS)

shapiro.test(aov(M_FASS)$residuals)

par(mfrow=c(1,2))
hist(aov(M_FASS)$residuals,
     main="A. Histogramm of residuals",
     xlab="Residuals value",
     ylab="Frequency")

qqPlot(aov(M_FASS)$residuals, id=FALSE,
       main="B. Q-Q Plot",
       xlab="Theoretical Quantiles",
       ylab="Residuals value")

leveneTest(Fertility ~ Karyotype, FASS)


# Statistical Analysis -- Non-parametric Test

kruskal.test(Fertility ~ Karyotype, FASS)

# Box plot

FASS_3RP <- FASS[c(which(str_detect(FASS$Karyotype, "3RP|STD"))),]
FASS_2Lt <- FASS[c(which(str_detect(FASS$Karyotype, "2Lt|STD"))),]
FASS_3RK <- FASS[c(which(str_detect(FASS$Karyotype, "3RK|STD"))),]
FASS_Homoz <- FASS[-c(which(str_detect(FASS$Karyotype, "HET"))),]


FASS_b1 <- ggplot(FASS_3RP, aes(x = Karyotype, y = Fertility, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170) +
  plot_theme+
  ylab("Number of offspring")+
  labs(tag="A")

FASS_b2 <- ggplot(FASS_2Lt, aes(x = Karyotype, y = Fertility, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170) +
  plot_theme+
  ylab("Number of offspring")+
  labs(tag="B")


FASS_b3 <- ggplot(FASS_3RK, aes(x = Karyotype, y = Fertility, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170) +
  plot_theme+
  ylab("Number of offspring")+
  labs(tag="C")

FASS_b4 <- ggplot(FASS_Homoz, aes(x = Karyotype, y = Fertility, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0,170) +
  plot_theme+
  ylab("Number of offspring")+
  labs(tag="D")

Final <- (FASS_b1|FASS_b2)/(FASS_b3|FASS_b4)

# Export plot as svg file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/2026_D.melanogaster_Stress_Experiment/Plots/Adult_Starvation_Stress_Fertility.svg",
  plot = Final,          
  width = 45,                    
  height = 25,                    
  units = "cm",                   
  dpi = 300                        
)



