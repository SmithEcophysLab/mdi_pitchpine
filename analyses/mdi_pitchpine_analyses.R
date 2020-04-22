# script to analyze mdi pitch pine data

library(tidyverse)
library(emmeans)

## read in cleaned data
data = read.csv('../data/mdi_all_clean.csv')
data$CN_foliar = data$C_foliar/data$N_foliar
data$CN_soil = data$C_soil/data$N_soil
head(data)

## fit models and explore results

### height
height_lm = lm(log(height) ~ Name, data = data)
#plot(resid(height_lm) ~ fitted(height_lm))
anova(height_lm)
cld(emmeans(height_lm, ~Name))

ggplot(data = data, aes(x = Name, y = log(height))) +
  geom_boxplot()

### canopy
canopy_lm = lm(log(canopy) ~ Name, data = data)
#plot(resid(canopy_lm) ~ fitted(canopy_lm))
anova(canopy_lm)
cld(emmeans(canopy_lm, ~Name))

ggplot(data = data, aes(x = Name, y = log(canopy))) +
  geom_boxplot()

### diam
diam_lm = lm(log(diam) ~ Name, data = data)
#plot(resid(diam_lm) ~ fitted(diam_lm))
anova(diam_lm)
cld(emmeans(diam_lm, ~Name))

ggplot(data = data, aes(x = Name, y = log(diam))) +
  geom_boxplot()

### d13C
d13C_lm = lm((d13C) ~ Name, data = data)
#plot(resid(d13C_lm) ~ fitted(d13C_lm))
anova(d13C_lm)
cld(emmeans(d13C_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (d13C))) +
  geom_boxplot()

### d15N
d15N_lm = lm((d15N) ~ Name, data = data)
#plot(resid(d15N_lm) ~ fitted(d15N_lm))
anova(d15N_lm)
cld(emmeans(d15N_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (d15N))) +
  geom_boxplot()

### C_foliar
C_foliar_lm = lm((C_foliar) ~ Name, data = data)
#plot(resid(C_foliar_lm) ~ fitted(C_foliar_lm))
anova(C_foliar_lm)
cld(emmeans(C_foliar_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (C_foliar))) +
  geom_boxplot()

### N_foliar
N_foliar_lm = lm((N_foliar) ~ Name, data = subset(data, N_foliar < 5))
#plot(resid(N_foliar_lm) ~ fitted(N_foliar_lm))
anova(N_foliar_lm)
cld(emmeans(N_foliar_lm, ~Name))

ggplot(data = subset(data, N_foliar < 5), aes(x = Name, y = (N_foliar))) +
  geom_boxplot()

### CN_foliar
CN_foliar_lm = lm((CN_foliar) ~ Name, data = subset(data, N_foliar < 5))
#plot(resid(CN_foliar_lm) ~ fitted(CN_foliar_lm))
anova(CN_foliar_lm)
cld(emmeans(CN_foliar_lm, ~Name))

ggplot(data = subset(data, N_foliar < 5), aes(x = Name, y = (CN_foliar))) +
  geom_boxplot()

### Ca_foliar
Ca_foliar_lm = lm((Ca_foliar) ~ Name, data = data)
#plot(resid(Ca_foliar_lm) ~ fitted(Ca_foliar_lm))
anova(Ca_foliar_lm)
cld(emmeans(Ca_foliar_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (Ca_foliar))) +
  geom_boxplot()

### P_foliar
P_foliar_lm = lm(log(P_foliar) ~ Name, data = data)
#plot(resid(P_foliar_lm) ~ fitted(P_foliar_lm))
anova(P_foliar_lm)
cld(emmeans(P_foliar_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (P_foliar))) +
  geom_boxplot()

### K_foliar
K_foliar_lm = lm(log(K_foliar) ~ Name, data = data)
#plot(resid(K_foliar_lm) ~ fitted(K_foliar_lm))
anova(K_foliar_lm)
cld(emmeans(K_foliar_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (K_foliar))) +
  geom_boxplot()

### Mg_foliar
Mg_foliar_lm = lm((Mg_foliar) ~ Name, data = data)
#plot(resid(Mg_foliar_lm) ~ fitted(Mg_foliar_lm))
anova(Mg_foliar_lm)
cld(emmeans(Mg_foliar_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (Mg_foliar))) +
  geom_boxplot()

### Al_foliar
Al_foliar_lm = lm((Al_foliar) ~ Name, data = data)
#plot(resid(Al_foliar_lm) ~ fitted(Al_foliar_lm))
anova(Al_foliar_lm)
cld(emmeans(Al_foliar_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (Al_foliar))) +
  geom_boxplot()

### Zn_foliar
Zn_foliar_lm = lm(log(Zn_foliar) ~ Name, data = data)
#plot(resid(Zn_foliar_lm) ~ fitted(Zn_foliar_lm))
anova(Zn_foliar_lm)
cld(emmeans(Zn_foliar_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (Zn_foliar))) +
  geom_boxplot()

### Ca_soil
Ca_soil_lm = lm((Ca_soil) ~ Name, data = data)
#plot(resid(Ca_soil_lm) ~ fitted(Ca_soil_lm))
anova(Ca_soil_lm)
cld(emmeans(Ca_soil_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (Ca_soil))) +
  geom_boxplot()

### P_soil
P_soil_lm = lm(log(P_soil) ~ Name, data = data)
#plot(resid(P_soil_lm) ~ fitted(P_soil_lm))
anova(P_soil_lm)
cld(emmeans(P_soil_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (P_soil))) +
  geom_boxplot()

### K_soil
K_soil_lm = lm((K_soil) ~ Name, data = data)
#plot(resid(K_soil_lm) ~ fitted(K_soil_lm))
anova(K_soil_lm)
cld(emmeans(K_soil_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (K_soil))) +
  geom_boxplot()

### Mg_soil
Mg_soil_lm = lm((Mg_soil) ~ Name, data = data)
#plot(resid(Mg_soil_lm) ~ fitted(Mg_soil_lm))
anova(Mg_soil_lm)
cld(emmeans(Mg_soil_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (Mg_soil))) +
  geom_boxplot()

### Al_soil
Al_soil_lm = lm(log(Al_soil) ~ Name, data = data)
#plot(resid(Al_soil_lm) ~ fitted(Al_soil_lm))
anova(Al_soil_lm)
cld(emmeans(Al_soil_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (Al_soil))) +
  geom_boxplot()

### Zn_soil
Zn_soil_lm = lm(log(Zn_soil) ~ Name, data = data)
#plot(resid(Zn_soil_lm) ~ fitted(Zn_soil_lm))
anova(Zn_soil_lm)
cld(emmeans(Zn_soil_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (Zn_soil))) +
  geom_boxplot()

### pH
pH_lm = lm((pH) ~ Name, data = data)
#plot(resid(pH_lm) ~ fitted(pH_lm))
anova(pH_lm)
cld(emmeans(pH_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (pH))) +
  geom_boxplot()

### CEC
CEC_lm = lm((CEC) ~ Name, data = data)
#plot(resid(CEC_lm) ~ fitted(CEC_lm))
anova(CEC_lm)
cld(emmeans(CEC_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (CEC))) +
  geom_boxplot()

### C_soil
C_soil_lm = lm((C_soil) ~ Name, data = data)
#plot(resid(C_soil_lm) ~ fitted(C_soil_lm))
anova(C_soil_lm)
cld(emmeans(C_soil_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (C_soil))) +
  geom_boxplot()

### N_soil
N_soil_lm = lm((N_soil) ~ Name, data = data)
#plot(resid(N_soil_lm) ~ fitted(N_soil_lm))
anova(N_soil_lm)
cld(emmeans(N_soil_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (N_soil))) +
  geom_boxplot()

### CN_soil
CN_soil_lm = lm(log(CN_soil) ~ Name, data = data)
#plot(resid(CN_soil_lm) ~ fitted(CN_soil_lm))
anova(CN_soil_lm)
cld(emmeans(CN_soil_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (CN_soil))) +
  geom_boxplot()

### retention
retention_lm = lm(asin(sqrt(0.01 * retention)) ~ Name, data = data)
#plot(resid(retention_lm) ~ fitted(retention_lm))
anova(retention_lm)
cld(emmeans(retention_lm, ~Name))

ggplot(data = data, aes(x = Name, y = (retention))) +
  geom_boxplot()


