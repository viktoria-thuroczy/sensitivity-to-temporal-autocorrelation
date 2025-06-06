rm(list = ls())
#load packages
library(dplyr)
library(magrittr)
library(ape)
library(nlme)
library(geiger)
library(ggplot2)


setwd("C:/Users/vikto/OneDrive/Documents/stochastic simulation paper/code/data")

tree <- read.tree("Tree.tre")
sim_data <- read.csv("LHT_SPGR.csv")
sen_data <- read.csv("LHT_sensitivity.csv", sep = ';')

tree$tip.label <- gsub("_", " ", tree$tip.label)


### Effect of autocorrelation on stochastic population growth rate 
spp <- sim_data$SpeciesAccepted
corBM <- corBrownian(phy = tree, form = ~ spp)

pgls1 <- gls(Stoch_GR ~ Autocorrelation, data = sim_data, correlation = corBM)
summary(pgls1)

# effect of species and autocorrelation on stochastic population growth rate 
pgls2 <- gls(Stoch_GR ~ SpeciesAccepted * Autocorrelation, data = sim_data, correlation = corBM)
anova(pgls2)


### Effect of Life history traits on Sensitivity
data_transform <- sen_data %>%
  mutate(Gen_T = log(Gen_T), 
         Net_Rep_Rate = sqrt(Net_Rep_Rate), 
         Age_Maturity = log(Age_Maturity),
         Sensitivity = sqrt(Sensitivity))

# Generation Time
spp <- data_transform$SpeciesAccepted
corBM <- corBrownian(phy = tree, form = ~ spp)

gen <- gls(Gen_T ~ Sensitivity, data = data_transform, correlation = corBM)
summary(gen)


# Net reproductive rate
spp <- data_transform$SpeciesAccepted
corBM <- corBrownian(phy = tree, form = ~ spp)

nrr <- gls(Net_Rep_Rate ~ Sensitivity, data = data_transform, correlation = corBM)
summary(nrr)


# Degree of iteroparity
doi_data <- data_transform %>%
  filter(!is.na(Do_Iteroparity)) %>%
  filter(!is.infinite(Do_Iteroparity))

spp <- doi_data$SpeciesAccepted
corBM <- corBrownian(phy = tree, form = ~ spp)

doi <- gls(Do_Iteroparity ~ Sensitivity, data = doi_data, correlation = corBM)
summary(doi)


# Rate of senescense 
ros_data <- sen_data %>%
  filter(!is.na(Ro_Senescence)) %>%
  filter(!is.infinite((Ro_Senescence)))

spp <- ros_data$SpeciesAccepted
corBM <- corBrownian(phy = tree, form = ~ spp)

ros <- gls(Ro_Senescence ~ sqrt(Sensitivity), data = ros_data, correlation = corBM)
summary(ros)


# Age at sexual maturity
spp <- data_transform$SpeciesAccepted
corBM <- corBrownian(phy = tree, form = ~ spp)

asm <- gls(Age_Maturity ~ Sensitivity, data = data_transform, correlation = corBM)
summary(asm)


## Without outliers 

# Generation time
data_wo <- sen_data %>%
  filter(Gen_T > quantile(Gen_T, 0.25, na.rm = FALSE) - 1.5*IQR(Gen_T) &
           Gen_T < quantile(Gen_T, 0.75, na.rm = FALSE) + 1.5*IQR(Gen_T))

spp <- data_wo$SpeciesAccepted
corBM <- corBrownian(phy = tree, form = ~ spp)

gen_wo <- gls(log(Gen_T) ~ sqrt(Sensitivity), data = data_wo, correlation = corBM)
summary(gen_wo)


# Net reproductive rate
data_wo <- sen_data %>%
  filter(Net_Rep_Rate > quantile(Net_Rep_Rate, 0.25, na.rm = FALSE) - 1.5*IQR(Net_Rep_Rate) &
           Net_Rep_Rate < quantile(Net_Rep_Rate, 0.75, na.rm = FALSE) + 1.5*IQR(Net_Rep_Rate))

spp <- data_wo$SpeciesAccepted
corBM <- corBrownian(phy = tree, form = ~ spp)

nrr_wo <- gls(sqrt(Net_Rep_Rate) ~ sqrt(Sensitivity), data = data_wo, correlation = corBM)
summary(nrr_wo)


# Degree of iteroparity
data_wo <- sen_data %>%
  filter(!is.na(Do_Iteroparity)) %>%
  filter(!is.infinite(Do_Iteroparity)) %>%
  filter(Do_Iteroparity > quantile(Do_Iteroparity, 0.25, na.rm = FALSE) - 1.5*IQR(Do_Iteroparity) &
           Do_Iteroparity < quantile(Do_Iteroparity, 0.75, na.rm = FALSE) + 1.5*IQR(Do_Iteroparity))

spp <- data_wo$SpeciesAccepted
corBM <- corBrownian(phy = tree, form = ~ spp)

doi_wo <- gls(Do_Iteroparity ~ sqrt(Sensitivity), data = data_wo, correlation = corBM)
summary(doi_wo)


# Rate of senescence 
data_wo <- sen_data %>%
  filter(!is.na(Ro_Senescence)) %>%
  filter(!is.infinite(Ro_Senescence)) %>%
  filter(Ro_Senescence > quantile(Ro_Senescence, 0.25, na.rm = FALSE) - 1.5*IQR(Ro_Senescence) &
           Ro_Senescence < quantile(Ro_Senescence, 0.75, na.rm = FALSE) + 1.5*IQR(Ro_Senescence))

spp <- data_wo$SpeciesAccepted
corBM <- corBrownian(phy = tree, form = ~ spp)

ros_wo <- gls(Ro_Senescence ~ sqrt(Sensitivity), data = data_wo, correlation = corBM)
summary(ros_wo)


# Age at maturity
data_wo <- sen_data %>%
  filter(!is.na(Age_Maturity)) %>%
  filter(!is.infinite(Age_Maturity)) %>%
  filter(Age_Maturity > quantile(Age_Maturity, 0.25, na.rm = FALSE) - 1.5*IQR(Age_Maturity) &
           Age_Maturity < quantile(Age_Maturity, 0.75, na.rm = FALSE) + 1.5*IQR(Age_Maturity))

spp <- data_wo$SpeciesAccepted
corBM <- corBrownian(phy = tree, form = ~ spp)

asm_wo <- gls(log(Age_Maturity) ~ sqrt(Sensitivity), data = data_wo, correlation = corBM)
summary(asm_wo)


# make data frame to store p-values
w <- c(summary(gen)$tTable[2,4], 
       summary(nrr)$tTable[2,4], 
       summary(doi)$tTable[2,4], 
       summary(ros)$tTable[2,4],
       summary(asm)$tTable[2,4])

w_o <- c(summary(gen_wo)$tTable[2,4], 
         summary(nrr_wo)$tTable[2,4], 
         summary(doi_wo)$tTable[2,4], 
         summary(ros_wo)$tTable[2,4],
         summary(asm_wo)$tTable[2,4])

p <- data.frame(w, w_o)

rownames(p) <- c("Gen_T", "Net_Rep_Rate", "Do_Iteroparity", "Ro_Senescence", "Age_Maturity")
colnames(p) <- c("With_Outliers", "Without_Outliers")

# Disable scientific notation
options(scipen = 999)

print(p, digits = 5)
