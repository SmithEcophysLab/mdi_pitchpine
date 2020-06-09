# script to analyze mdi pitch pine data

library(tidyverse)
library(emmeans)
library(lme4)
library(car)

multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}

## read in cleaned data
data = read.csv('../data/mdi_all_clean.csv')
data$CN_foliar = data$C_foliar/data$N_foliar
data$CN_soil = data$C_soil/data$N_soil
data$fire[data$Name == 'CAD'] = 'fire' 
data$fire[data$Name == 'CADCLIFFS'] = 'fire'
data$fire[data$Name == 'STSAUV'] = 'no fire'
data$fire[data$Name == 'WOND'] = 'no fire'
head(data)

## site means
data_group_by_Name = group_by(data, Name)
data_Name_means = summarise(data_group_by_Name,
                            Elevation_mean = mean(Elevation, na.rm = T),
                            Slope_mean = mean(Slope, na.rm = T),
                            Aspect_mean = mean(Aspect, na.rm = T))

## fit models and explore results

### elevation
Elevation_lm = lm(log(Elevation) ~ Name, data = data)
#plot(resid(Elevation_lm) ~ fitted(Elevation_lm))
Anova(Elevation_lm)
cld(emmeans(Elevation_lm, ~Name), alpha = 0.1)

### height
height_lm = lm(log(height) ~ Name, data = data)
#plot(resid(height_lm) ~ fitted(height_lm))
Anova(height_lm)
cld(emmeans(height_lm, ~Name), alpha = 0.1)

height_lmer_cont = lmer(log(height) ~ Elevation * fire + (1|Name), data = data)
Anova(height_lmer_cont)
test(emtrends(height_lmer_cont, ~fire, var = 'Elevation'))

height_plot = ggplot(data = data, aes(x = Name, y = log(height), col = fire)) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_boxplot(outlier.color = NA, fill = 'white') +
  geom_dotplot(binaxis = 'y', binwidth = 0.07, stackdir = 'center', alpha = 0.5) +
  # scale_x_discrete(labels = c('Ambient', 'Added N')) +
  xlab('Site') +
  ylab(expression('ln(Height)'))

height_plot_elevation = ggplot(data = data, aes(x = Elevation, y = log(height), col = fire)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_point(size = 6) +
  ylab(expression('ln(Height)'))

jpeg(filename = "plots/height_plot.jpeg", width = 1000, height = 600, units = 'px')
multiplot(height_plot, height_plot_elevation, cols = 2)
dev.off()

### canopy
canopy_lm = lm(log(canopy) ~ Name, data = data)
#plot(resid(canopy_lm) ~ fitted(canopy_lm))
anova(canopy_lm)
cld(emmeans(canopy_lm, ~Name), alpha = 0.1)

canopy_lmer_cont = lmer(log(canopy) ~ Elevation * fire + (1|Name), data = data)
Anova(canopy_lmer_cont)
test(emtrends(canopy_lmer_cont, ~fire, var = 'Elevation'))

canopy_plot = ggplot(data = data, aes(x = Name, y = log(canopy), col = fire)) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_boxplot(outlier.color = NA, fill = 'white') +
  geom_dotplot(binaxis = 'y', binwidth = 0.07, stackdir = 'center', alpha = 0.5) +
  # scale_x_discrete(labels = c('Ambient', 'Added N')) +
  xlab('Site') +
  ylab(expression('ln(Canopy)'))

canopy_plot_elevation = ggplot(data = data, aes(x = Elevation, y = log(canopy), col = fire)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_point(size = 6) +
  ylab(expression('ln(Canopy)'))

jpeg(filename = "plots/canopy_plot.jpeg", width = 1000, height = 600, units = 'px')
multiplot(canopy_plot, canopy_plot_elevation, cols = 2)
dev.off()

### diam
diam_lm = lm(log(diam) ~ Name, data = data)
#plot(resid(diam_lm) ~ fitted(diam_lm))
anova(diam_lm)
cld(emmeans(diam_lm, ~Name), alpha = 0.1)

diam_lmer_cont = lmer(log(diam) ~ Elevation * fire + (1|Name), data = data)
Anova(diam_lmer_cont)
test(emtrends(diam_lmer_cont, ~fire, var = 'Elevation'))

diam_plot = ggplot(data = data, aes(x = Name, y = log(diam), col = fire)) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_boxplot(outlier.color = NA, fill = 'white') +
  geom_dotplot(binaxis = 'y', binwidth = 0.07, stackdir = 'center', alpha = 0.5) +
  # scale_x_discrete(labels = c('Ambient', 'Added N')) +
  xlab('Site') +
  ylab(expression('ln(Diameter)'))

diam_plot_elevation = ggplot(data = data, aes(x = Elevation, y = log(diam), col = fire)) +
  theme(legend.position = "right", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_point(size = 6) +
  ylab(expression('ln(Diameter)'))

jpeg(filename = "plots/diam_plot.jpeg", width = 1000, height = 600, units = 'px')
multiplot(diam_plot, diam_plot_elevation, cols = 2)
dev.off()

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


