library(dplyr) 
library(ggplot2) 
library(rstatix) 
library(car)
library(stringr)  
library(patchwork)
library(emmeans)



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
  limits = c(0, 1.15),
  breaks = c(0.2, 0.4, 0.6, 0.8,1)
)


# Import & Prepare Data

DataSASS.M <-read.csv2(file.choose())
SASS.M <- DataSASS.M

SASS.M$M.Survival <- SASS.M$ALIVE_M/SASS.M$Total_M

boxplot(M.Survival ~ GAS, SASS.M)
t.test(M.Survival ~ GAS, SASS.M, var.equal=TRUE)

 
boxplot(M.Survival ~ BORDERS, SASS.M)
t.test(M.Survival ~ BORDERS, SASS.M, var.equal=TRUE)

SASS.M<- SASS.M[SASS.M$GAS == "YES",]
SASS.M <- SASS.M[SASS.M$BORDERS == "NO",]

# Assumptions

M_SASS.M <- lm(M.Survival ~ Karyotype, SASS.M)
shapiro.test(aov(M_SASS.M)$residuals)

par(mfrow=c(1,2))
hist(aov(M_SASS.M)$residuals,
     main="A. Histogramm of residuals",
     xlab="Residuals value",
     ylab="Frequency")

qqPlot(aov(M_SASS.M)$residuals, id=FALSE,
       main="B. Q-Q Plot",
       xlab="Theoretical Quantiles",
       ylab="Residuals value")

leveneTest(M.Survival ~ Karyotype, SASS.M)

# Statistical Analysis -- ANOVA

anova(M_SASS.M)

# Contrast Test

SASS.M_3RP <- SASS.M[c(which(str_detect(SASS.M$Karyotype, "3RP|STD"))),]
SASS.M_2Lt <- SASS.M[c(which(str_detect(SASS.M$Karyotype, "2Lt|STD"))),]
SASS.M_3RK <- SASS.M[c(which(str_detect(SASS.M$Karyotype, "3RK|STD"))),]
SASS.M_Homoz <- SASS.M[-c(which(str_detect(SASS.M$Karyotype, "HET"))),]


a.18 = 1-(1-0.05)^(1/18)

M_SASS.M.3RP <- lm(M.Survival ~ Karyotype, SASS.M_3RP)
means.SASS.M.3RP <- emmeans(M_SASS.M.3RP, "Karyotype")
All.comparisons.3RP.SASS.M <- pairs(means.SASS.M.3RP, adjust="none") 
All.comparisons.3RP.SASS.M <- as.data.frame(All.comparisons.3RP.SASS.M)

M_SASS.M.2Lt <- lm(M.Survival ~ Karyotype, SASS.M_2Lt)
means.SASS.2Lt <- emmeans(M_SASS.M.2Lt, "Karyotype")
All.comparisons.2Lt.SASS.M <- pairs(means.SASS.2Lt, adjust="none")
All.comparisons.2Lt.SASS.M <- as.data.frame(All.comparisons.2Lt.SASS.M)

M_SASS.M.3RK <- lm(M.Survival ~ Karyotype, SASS.M_3RK)
means.SASS.3RK <- emmeans(M_SASS.M.3RK, "Karyotype")
All.comparisons.3RK.SASS.M <- pairs(means.SASS.3RK, adjust="none")
All.comparisons.3RK.SASS.M <- as.data.frame(All.comparisons.3RK.SASS.M)


Contrast.SASS.M <- rbind(All.comparisons.3RP.SASS.M, All.comparisons.2Lt.SASS.M, All.comparisons.3RK.SASS.M)

S_Contrast.SASS.M <- Contrast.SASS.M[Contrast.SASS.M$p.value < a.18,]
S_Contrast.SASS.M <- mutate(S_Contrast.SASS.M, p.signif =case_when(p.value< a.18 ~ "*", TRUE ~ "ns"))



# Box plot

SASS.M_b1 <- ggplot(SASS.M_3RP, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %")+
  labs(tag="A")

SASS.M_b2 <- ggplot(SASS.M_2Lt, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %") +
  labs(tag="B")


SASS.M_b3 <- ggplot(SASS.M_3RK, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %") +
  labs(tag="C")


SASS.M_b4 <- ggplot(SASS.M_Homoz, aes(x = Karyotype, y = M.Survival, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors,labels = function(x) gsub("_", " ", x)) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  y_axis+
  plot_theme+
  ylab("Male survival rate in %") +
  labs(tag="D")



Final <- (SASS.M_b1|SASS.M_b2)/(SASS.M_b3|SASS.M_b4) 


# Export plot as png file
ggsave(
  filename = "D:/Users/Marine Caussignac/UniFR/Cours/25-26/Travail de Bachelor/2026_D.melanogaster_Stress_Experiment/Plots/Male_Starvation_Stress.png",
  plot = Final,          
  width = 45,                    
  height = 25,                    
  units = "cm",                   
  dpi = 300                        
)

