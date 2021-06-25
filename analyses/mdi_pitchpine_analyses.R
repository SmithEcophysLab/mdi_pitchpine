# script to analyze mdi pitch pine data

# load packages
library(tidyverse)
library(emmeans)
library(lme4)
library(car)
library(circular)
library(multcompView)
library(ggthemes)
library(agricolae)
library(patchwork)

#### read in cleaned data ####

data <- read.csv('../data/mdi_all_clean.csv')
data$CN_foliar <- data$C_foliar/data$N_foliar
data$CN_soil <- data$C_soil/data$N_soil

data_density <- read.csv('../data/mdi_stand_density.csv')

# add density data to all data
data <- left_join(data, data_density, by = "ID")

## assign fire history status to each site
data$Fire[data$Name == 'CAD'] <- 'fire' 
data$Fire[data$Name == 'CADCLIFFS'] <- 'fire'
data$Fire[data$Name == 'STSAUV'] <- 'no fire'
data$Fire[data$Name == 'WOND'] <- 'no fire'

## create an elevation factor
data$Elevation_fac[data$Name == 'CAD' | data$Name == 'STSAUV'] <- 'high'
data$Elevation_fac[data$Name == 'CADCLIFFS' | data$Name == 'WOND'] <- 'low'

## reorder levels for elevation factor from "low" to "high"
data$Elevation_fac <- factor(data$Elevation_fac, levels = c("low", "high"))

## rename site labels to match manuscript
data$Site[data$Name == "CAD"] <- "SCT"
data$Site[data$Name == "CADCLIFFS"] <- "GOR"
data$Site[data$Name == "STSAUV"] <- "STS"
data$Site[data$Name == "WOND"] <- "WON"

#### fit models and explore results ####

### aspect
#### turn aspect data into circular data that maps onto a compass
aspect_SCT <- circular(data$Aspect[data$Name == 'CAD'],
                      units = "degrees", template = "geographics")
aspect_GOR <- circular(data$Aspect[data$Name == 'CADCLIFFS'],
                            units = "degrees", template = "geographics")
aspect_STS <- circular(data$Aspect[data$Name == 'STSAUV'],
                         units = "degrees", template = "geographics")
aspect_WON <- circular(data$Aspect[data$Name == 'WOND'],
                       units = "degrees", template = "geographics")

#### Watson's Two Sample Test of Homogeneity (https://bigdata.duke.edu/sites/bigdata.duke.edu/files/site-images/FullLesson.html)
#### tests whether the ‘north’ and ‘south’ groups orient in different directions
watson.two.test(aspect_WON, aspect_GOR) # 0.01 < P < 0.05
watson.two.test(aspect_WON, aspect_STS) # 0.001 < P < 0.01
watson.two.test(aspect_WON, aspect_SCT) # 0.01 < P < 0.05
watson.two.test(aspect_GOR, aspect_STS) # 0.001 < P < 0.01
watson.two.test(aspect_GOR, aspect_SCT) # 0.05 < P < 0.1
watson.two.test(aspect_STS, aspect_SCT) # P < 0.001

#### circular plots for each site
jpeg(filename = "plots/plots_aspect.jpeg", width = 3000, 
     height = 3000, units = 'px')
par(mfrow = c(2, 2), cex.lab = 6, cex.main = 6, mar = c(5.5, 8.5, 5.5, 2.5))
plot_aspect_GOR <- plot.circular(aspect_GOR, main = 'Gorham Cliffs (a)', 
                                      ylab = "Fire", 
                                      cex = 8, col = "red", pch = 16)
plot_aspect_SCT <- plot.circular(aspect_SCT, main = 'South Cadillac (a)',
                                cex = 8, col = "red", pch = 17)
plot_aspect_WON <- plot.circular(aspect_WON, main = 'Wonderland (b)', 
                                 ylab = "No Fire", xlab = "Low Elevation", 
                                 cex = 8, col = "blue", pch = 16)
plot_aspect_STS <- plot.circular(aspect_STS, main = 'St. Sauveur (c)', 
                                   xlab = "High Elevation", 
                                   cex = 8, col = "blue", pch = 17)
dev.off()

## slope
slope_lm <- lmer(Slope ~ Elevation * Fire + (1 | Site), data = data)
# plot(resid(slope_lm) ~ fitted(slope_lm))
Anova(slope_lm)
summary(slope_lm) # N = 60

(plot_slope <- ggplot(data = data, aes(x = Elevation, y = Slope)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1500)) +
    scale_y_continuous(name = "Slope (◦)") +
    guides(color = guide_legend("Fire History")))

## allometry
### height
height_lm <- lmer(log(Height) ~ Elevation * Fire + (1 | Site), data = data)
# plot(resid(height_lm) ~ fitted(height_lm))
Anova(height_lm)
summary(height_lm) # N = 40

emtrends(height_lm, ~ Fire, var = "Elevation") # fire = -0.000905 | no fire = 0.000405
emmeans(height_lm, ~ Fire, at = list(Elevation = 0)) # fire = 1.67 | no fire = 1.25

height_f_slope <- summary(emtrends(height_lm, ~ Fire, var = "Elevation"))[1, 2] 
height_f_intercept <- summary(emmeans(height_lm, ~ Fire, at = list(Elevation = 0)))[1, 2] 
height_f_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
height_f_trend <- height_f_intercept + height_f_seq * height_f_slope
height_f_trend <- as.data.frame(cbind(height_f_seq, height_f_trend))

height_nf_slope <- summary(emtrends(height_lm, ~ Fire, var = "Elevation"))[2, 2] 
height_nf_intercept <- summary(emmeans(height_lm, ~ Fire, at = list(Elevation = 0)))[2, 2] 
height_nf_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
height_nf_trend <- height_nf_intercept + height_nf_seq * height_nf_slope
height_nf_trend <- as.data.frame(cbind(height_nf_seq, height_nf_trend))

(plot_height <- ggplot(data = data, aes(x = Elevation, y = log(Height))) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = height_f_trend, aes(x = height_f_seq, y = height_f_trend), 
              col = 'red', lwd = 2, alpha = 0.8) +
    geom_line(data = height_nf_trend, aes(x = height_nf_seq, y = height_nf_trend), 
              col = 'blue', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    scale_y_continuous(name = "Height (m)") +
    guides(color = guide_legend("Fire History")))

### canopy
canopy_lm <- lmer(log(Canopy) ~ Elevation * Fire + (1 | Site), data = data)
# plot(resid(canopy_lm) ~ fitted(canopy_lm))
Anova(canopy_lm)
summary(canopy_lm) # N = 40

emtrends(canopy_lm, ~ Elevation, var = "Elevation") # slope = -0.000726
emmeans(canopy_lm, ~ Elevation, at = list(Elevation = 0)) # intercept = 1.27

canopy_slope <- summary(emtrends(canopy_lm, ~ Elevation, var = "Elevation"))[1, 2] 
canopy_intercept <- summary(emmeans(canopy_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
canopy_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
canopy_trend <- canopy_intercept + canopy_seq * canopy_slope
canopy_trend <- as.data.frame(cbind(canopy_seq, canopy_trend))

(plot_canopy <- ggplot(data = data, aes(x = Elevation, y = log(Canopy))) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = canopy_trend, aes(x = canopy_seq, y = canopy_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    scale_y_continuous(name = "Canopy Spread (m)") +
    guides(color = guide_legend("Fire History")))

### diam
diam_lm <- lmer(log(Diam) ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(diam_lm) ~ fitted(diam_lm))
Anova(diam_lm)
summary(diam_lm) # N = 40

emtrends(diam_lm, ~ Elevation, var = "Elevation") # slope = -0.000886
emmeans(diam_lm, ~ Elevation, at = list(Elevation = 0)) # intercept = 3.21

diam_slope <- summary(emtrends(diam_lm, ~ Elevation, var = "Elevation"))[1, 2] 
diam_intercept <- summary(emmeans(diam_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
diam_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
diam_trend <- diam_intercept + diam_seq * diam_slope
diam_trend <- as.data.frame(cbind(diam_seq, diam_trend))

(plot_diam <- ggplot(data = data, aes(x = Elevation, y = log(Diam))) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = diam_trend, aes(x = diam_seq, y = diam_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) +
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    scale_y_continuous(name = "DBH (cm)", limits = c(1, 5)) +
    guides(color = guide_legend("Fire History")))

### density
density_lm <- lmer(mean_distance ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(density_lm) ~ fitted(density_lm))
Anova(density_lm)
summary(density_lm) # N = 60

emtrends(density_lm, ~ Fire, var = "Elevation") # fire = 0.000222 | no fire = 0.006379
emmeans(density_lm, ~ Fire, at = list(Elevation = 0)) # fire = 3.36 | no fire = 1.03

density_f_slope <- summary(emtrends(density_lm, ~ Fire, var = "Elevation"))[1, 2] 
density_f_intercept <- summary(emmeans(density_lm, ~ Fire, at = list(Elevation = 0)))[1, 2] 
density_f_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
density_f_trend <- density_f_intercept + density_f_seq * density_f_slope
density_f_trend <- as.data.frame(cbind(density_f_seq, density_f_trend))

density_nf_slope <- summary(emtrends(density_lm, ~ Fire, var = "Elevation"))[2, 2] 
density_nf_intercept <- summary(emmeans(density_lm, ~ Fire, at = list(Elevation = 0)))[2, 2] 
density_nf_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
density_nf_trend <- density_nf_intercept + density_nf_seq * density_nf_slope
density_nf_trend <- as.data.frame(cbind(density_nf_seq, density_nf_trend))

(plot_density <- ggplot(data = data, aes(x = Elevation, y = mean_distance)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = density_f_trend, aes(x = density_f_seq, y = density_f_trend), 
              col = 'red', lwd = 2, alpha = 0.8) +
    geom_line(data = density_nf_trend, aes(x = density_nf_seq, y = density_nf_trend), 
              col = 'blue', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) +
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1500)) +
    scale_y_continuous(name = "Distance Between Neighbors (m)", limits = c(0, 6)) + 
    guides(color = guide_legend("Fire History")))

(plots_allometry <- plot_canopy + plot_diam + plot_density + plot_height + 
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## foliar isotopes
### d13C
d13C_lm <- lmer(d13C ~ Elevation * Fire + (1 | Site), data = data)
plot(resid(d13C_lm) ~ fitted(d13C_lm))
Anova(d13C_lm)
summary(d13C_lm) # N = 55

emtrends(d13C_lm, ~ Elevation, var = "Elevation") # slope = -0.00556
emmeans(d13C_lm, ~ Elevation, at = list(Elevation = 0)) # intercept = -26.5

d13C_slope <- summary(emtrends(d13C_lm, ~ Elevation, var = "Elevation"))[1, 2] 
d13C_intercept <- summary(emmeans(d13C_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
d13C_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
d13C_trend <- d13C_intercept + d13C_seq * d13C_slope
d13C_trend <- as.data.frame(cbind(d13C_seq, d13C_trend))

(plot_d13C <- ggplot(data = data, aes(x = Elevation, y = d13C)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = d13C_trend, aes(x = d13C_seq, y = d13C_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression(delta^{"13"}*"C (‰)")) +
    guides(color = guide_legend("Fire History")))

### d15N
d15N_lm <- lmer(d15N ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(d15N_lm) ~ fitted(d15N_lm))
Anova(d15N_lm)
summary(d15N_lm) # N = 55

## foliar organics
### C_foliar
C_foliar_lm <- lmer(C_foliar ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(C_foliar_lm) ~ fitted(C_foliar_lm))
Anova(C_foliar_lm)
summary(C_foliar_lm) # N = 60

(plot_C_foliar <- ggplot(data = data, aes(x = Fire, y = C_foliar)) +
    geom_jitter(height = 0, aes(color = Fire), size = 2) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    scale_color_manual(values = c('red', 'blue')) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Fire History") +
    ylab(expression("Foliar Carbon (g g"^{-1}*")")) +
    guides(color = "none"))

### N_foliar
N_foliar_lm <- lmer(N_foliar ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(N_foliar_lm) ~ fitted(N_foliar_lm))
Anova(N_foliar_lm)
summary(N_foliar_lm) # N = 56

### CN_foliar
CN_foliar_lm <- lmer(CN_foliar ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(CN_foliar_lm) ~ fitted(CN_foliar_lm))
Anova(CN_foliar_lm)
summary(CN_foliar_lm) # N = 56

## foliar inorganics
### Ca_foliar
Ca_foliar_lm <- lmer(Ca_foliar ~ Elevation * Fire + (1 | Site), data = data)
# plot(resid(Ca_foliar_lm) ~ fitted(Ca_foliar_lm))
Anova(Ca_foliar_lm)
summary(Ca_foliar_lm) # N = 40

emtrends(Ca_foliar_lm, ~ Elevation, var = "Elevation") # slope = -0.00129
emmeans(Ca_foliar_lm, ~ Elevation, at = list(Elevation = 0)) # intercept = 1.82

Ca_foliar_slope <- summary(emtrends(Ca_foliar_lm, ~ Elevation, var = "Elevation"))[1, 2] 
Ca_foliar_intercept <- summary(emmeans(Ca_foliar_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
Ca_foliar_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
Ca_foliar_trend <- Ca_foliar_intercept + Ca_foliar_seq * Ca_foliar_slope
Ca_foliar_trend <- as.data.frame(cbind(Ca_foliar_seq, Ca_foliar_trend))

(plot_Ca_foliar <- ggplot(data = data, aes(x = Elevation, y = Ca_foliar)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = Ca_foliar_trend, aes(x = Ca_foliar_seq, y = Ca_foliar_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Foliar Calcium (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### P_foliar
P_foliar_lm <- lmer(log(P_foliar) ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(P_foliar_lm) ~ fitted(P_foliar_lm))
Anova(P_foliar_lm)
summary(P_foliar_lm) # N = 40

### K_foliar
K_foliar_lm <- lmer(log(K_foliar) ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(K_foliar_lm) ~ fitted(K_foliar_lm))
Anova(K_foliar_lm)
summary(K_foliar_lm) # N = 40

emtrends(K_foliar_lm, ~ Fire, var = "Elevation") # fire = -0.001582 | no fire = 0.000584
emmeans(K_foliar_lm, ~ Fire, at = list(Elevation = 0)) # fire = 0.904 | no fire = 0.731

K_foliar_f_slope <- summary(emtrends(K_foliar_lm, ~ Fire, var = "Elevation"))[1, 2] 
K_foliar_f_intercept <- summary(emmeans(K_foliar_lm, ~ Fire, at = list(Elevation = 0)))[1, 2] 
K_foliar_f_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
K_foliar_f_trend <- K_foliar_f_intercept + K_foliar_f_seq * K_foliar_f_slope
K_foliar_f_trend <- as.data.frame(cbind(K_foliar_f_seq, K_foliar_f_trend))

K_foliar_nf_slope <- summary(emtrends(K_foliar_lm, ~ Fire, var = "Elevation"))[2, 2] 
K_foliar_nf_intercept <- summary(emmeans(K_foliar_lm, ~ Fire, at = list(Elevation = 0)))[2, 2] 
K_foliar_nf_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
K_foliar_nf_trend <- K_foliar_nf_intercept + K_foliar_nf_seq * K_foliar_nf_slope
K_foliar_nf_trend <- as.data.frame(cbind(K_foliar_nf_seq, K_foliar_nf_trend))

(plot_K_foliar <- ggplot(data = data, aes(x = Elevation, y = log(K_foliar))) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = K_foliar_f_trend, aes(x = K_foliar_f_seq, y = K_foliar_f_trend), 
              col = 'red', lwd = 2, alpha = 0.8) +
    geom_line(data = K_foliar_nf_trend, aes(x = K_foliar_nf_seq, y = K_foliar_nf_trend), 
              col = 'blue', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Foliar Potassium (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### Mg_foliar
Mg_foliar_lm <- lmer(Mg_foliar ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(Mg_foliar_lm) ~ fitted(Mg_foliar_lm))
Anova(Mg_foliar_lm)
summary(Mg_foliar_lm) # N = 40

### Al_foliar
Al_foliar_lm <- lmer(Al_foliar ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(Al_foliar_lm) ~ fitted(Al_foliar_lm))
Anova(Al_foliar_lm)
summary(Al_foliar_lm) # N = 40

### Zn_foliar
Zn_foliar_lm <- lmer(log(Zn_foliar) ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(Zn_foliar_lm) ~ fitted(Zn_foliar_lm))
Anova(Zn_foliar_lm)
summary(Zn_foliar_lm) # N = 40

(plots_foliar <- (plot_d13C + plot_C_foliar + plot_Ca_foliar + plot_K_foliar) + 
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## soil organics
### C_soil
C_soil_lm <- lmer(C_soil ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(C_soil_lm) ~ fitted(C_soil_lm))
Anova(C_soil_lm)
summary(C_soil_lm) # N = 31

emtrends(C_soil_lm, ~ Elevation, var = "Elevation") # slope = -0.0141
emmeans(C_soil_lm, ~ Elevation, at = list(Elevation = 0)) # intercept = 24.4

C_soil_slope <- summary(emtrends(C_soil_lm, ~ Elevation, var = "Elevation"))[1, 2] 
C_soil_intercept <- summary(emmeans(C_soil_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
C_soil_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
C_soil_trend <- C_soil_intercept + C_soil_seq * C_soil_slope
C_soil_trend <- as.data.frame(cbind(C_soil_seq, C_soil_trend))

(plot_C_soil <- ggplot(data = data, aes(x = Elevation, y = C_soil)) +
    geom_jitter(height = 0, aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = C_soil_trend, aes(x = C_soil_seq, y = C_soil_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Carbon (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### N_soil
N_soil_lm <- lmer(N_soil ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(N_soil_lm) ~ fitted(N_soil_lm))
Anova(N_soil_lm)
summary(N_soil_lm) # N = 26

### CN_soil
CN_soil_lm <- lmer(log(CN_soil) ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(CN_soil_lm) ~ fitted(CN_soil_lm))
Anova(CN_soil_lm)
summary(CN_soil_lm) # N = 26

emtrends(CN_soil_lm, ~ Elevation, var = "Elevation") # slope = -0.00106
emmeans(CN_soil_lm, ~ Elevation, at = list(Elevation = 0)) # intercept = 3.78

CN_soil_slope <- summary(emtrends(CN_soil_lm, ~ Elevation, var = "Elevation"))[1, 2] 
CN_soil_intercept <- summary(emmeans(CN_soil_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
CN_soil_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
CN_soil_trend <- CN_soil_intercept + CN_soil_seq * CN_soil_slope
CN_soil_trend <- as.data.frame(cbind(CN_soil_seq, CN_soil_trend))

(plot_CN_soil <- ggplot(data = data, aes(x = Elevation, y = log(CN_soil))) +
    geom_jitter(height = 0, aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = CN_soil_trend, aes(x = CN_soil_seq, y = CN_soil_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab("Soil Carbon/Nitrogen") +
    guides(color = guide_legend("Fire History")))

## soil inorganics
### Ca_soil
Ca_soil_lm <- lmer(Ca_soil ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(Ca_soil_lm) ~ fitted(Ca_soil_lm))
Anova(Ca_soil_lm)
summary(Ca_soil_lm) # N = 31

emtrends(Ca_soil_lm, ~ Elevation, var = "Elevation") # slope = -0.000615
emmeans(Ca_soil_lm, ~ Elevation, at = list(Elevation = 0)) # intercept = 0.785

Ca_soil_slope <- summary(emtrends(Ca_soil_lm, ~ Elevation, var = "Elevation"))[1, 2] 
Ca_soil_intercept <- summary(emmeans(Ca_soil_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
Ca_soil_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
Ca_soil_trend <- Ca_soil_intercept + Ca_soil_seq * Ca_soil_slope
Ca_soil_trend <- as.data.frame(cbind(Ca_soil_seq, Ca_soil_trend))

(plot_Ca_soil <- ggplot(data = data, aes(x = Elevation, y = Ca_soil)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = Ca_soil_trend, aes(x = Ca_soil_seq, y = Ca_soil_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Calcium (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### P_soil
P_soil_lm <- lmer(log(P_soil) ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(P_soil_lm) ~ fitted(P_soil_lm))
Anova(P_soil_lm)
summary(P_soil_lm) # N = 31

emtrends(P_soil_lm, ~ Elevation, var = "Elevation") # slope = -0.00107 
emmeans(P_soil_lm, ~ Elevation, at = list(Elevation = 0)) # intercept = -5.4

P_soil_slope <- summary(emtrends(P_soil_lm, ~ Elevation, var = "Elevation"))[1, 2] 
P_soil_intercept <- summary(emmeans(P_soil_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
P_soil_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
P_soil_trend <- P_soil_intercept + P_soil_seq * P_soil_slope
P_soil_trend <- as.data.frame(cbind(P_soil_seq, P_soil_trend))

(plot_P_soil <- ggplot(data = data, aes(x = Elevation, y = log(P_soil))) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = P_soil_trend, aes(x = P_soil_seq, y = P_soil_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Phosphorus (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### K_soil
K_soil_lm <- lmer(K_soil ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(K_soil_lm) ~ fitted(K_soil_lm))
Anova(K_soil_lm)
summary(K_soil_lm) # N = 31

emtrends(K_soil_lm, ~ Elevation, var = "Elevation") # slope = -0.000128
emmeans(K_soil_lm, ~ Elevation, at = list(Elevation = 0)) # intercept = 0.276

K_soil_slope <- summary(emtrends(K_soil_lm, ~ Elevation, var = "Elevation"))[1, 2] 
K_soil_intercept <- summary(emmeans(K_soil_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
K_soil_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
K_soil_trend <- K_soil_intercept + K_soil_seq * K_soil_slope
K_soil_trend <- as.data.frame(cbind(K_soil_seq, K_soil_trend))

(plot_K_soil <- ggplot(data = data, aes(x = Elevation, y = K_soil)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = K_soil_trend, aes(x = K_soil_seq, y = K_soil_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Potassium (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### Mg_soil
Mg_soil_lm <- lmer(Mg_soil ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(Mg_soil_lm) ~ fitted(Mg_soil_lm))
Anova(Mg_soil_lm)
summary(Mg_soil_lm) # N = 31

### Al_soil
Al_soil_lm <- lmer(log(Al_soil) ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(Al_soil_lm) ~ fitted(Al_soil_lm))
Anova(Al_soil_lm)
summary(Al_soil_lm) # N = 31

emtrends(Al_soil_lm, ~ Fire, var = "Elevation") # fire = -0.00101 | no fire = 0.001
emmeans(Al_soil_lm, ~ Fire, at = list(Elevation = 0)) # fire = -1.7 | no fire = -2.42

Al_soil_f_slope <- summary(emtrends(Al_soil_lm, ~ Fire, var = "Elevation"))[1, 2] 
Al_soil_f_intercept <- summary(emmeans(Al_soil_lm, ~ Fire, at = list(Elevation = 0)))[1, 2] 
Al_soil_f_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
Al_soil_f_trend <- Al_soil_f_intercept + Al_soil_f_seq * Al_soil_f_slope
Al_soil_f_trend <- as.data.frame(cbind(Al_soil_f_seq, Al_soil_f_trend))

Al_soil_nf_slope <- summary(emtrends(Al_soil_lm, ~ Fire, var = "Elevation"))[2, 2] 
Al_soil_nf_intercept <- summary(emmeans(Al_soil_lm, ~ Fire, at = list(Elevation = 0)))[2, 2] 
Al_soil_nf_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
Al_soil_nf_trend <- Al_soil_nf_intercept + Al_soil_nf_seq * Al_soil_nf_slope
Al_soil_nf_trend <- as.data.frame(cbind(Al_soil_nf_seq, Al_soil_nf_trend))

(plot_Al_soil <- ggplot(data = data, aes(x = Elevation, y = log(Al_soil))) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = Al_soil_f_trend, aes(x = Al_soil_f_seq, y = Al_soil_f_trend), 
              col = 'red', lwd = 2, alpha = 0.8) +
    geom_line(data = Al_soil_nf_trend, aes(x = Al_soil_nf_seq, y = Al_soil_nf_trend), 
              col = 'blue', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Aluminum (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### Zn_soil
Zn_soil_lm <- lmer(log(Zn_soil) ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(Zn_soil_lm) ~ fitted(Zn_soil_lm))
Anova(Zn_soil_lm)
summary(Zn_soil_lm) # N = 31

emtrends(Zn_soil_lm, ~ Fire, var = "Elevation") # fire = -0.001267 | no fire = 0.000532
emmeans(Zn_soil_lm, ~ Fire, at = list(Elevation = 0)) # fire = -5.34 | no fire = -6.03

Zn_soil_f_slope <- summary(emtrends(Zn_soil_lm, ~ Fire, var = "Elevation"))[1, 2] 
Zn_soil_f_intercept <- summary(emmeans(Zn_soil_lm, ~ Fire, at = list(Elevation = 0)))[1, 2] 
Zn_soil_f_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
Zn_soil_f_trend <- Zn_soil_f_intercept + Zn_soil_f_seq * Zn_soil_f_slope
Zn_soil_f_trend <- as.data.frame(cbind(Zn_soil_f_seq, Zn_soil_f_trend))

Zn_soil_nf_slope <- summary(emtrends(Zn_soil_lm, ~ Fire, var = "Elevation"))[2, 2] 
Zn_soil_nf_intercept <- summary(emmeans(Zn_soil_lm, ~ Fire, at = list(Elevation = 0)))[2, 2] 
Zn_soil_nf_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
Zn_soil_nf_trend <- Zn_soil_nf_intercept + Zn_soil_nf_seq * Zn_soil_nf_slope
Zn_soil_nf_trend <- as.data.frame(cbind(Zn_soil_nf_seq, Zn_soil_nf_trend))

(plot_Zn_soil <- ggplot(data = data, aes(x = Elevation, y = log(Zn_soil))) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = Zn_soil_f_trend, aes(x = Zn_soil_f_seq, y = Zn_soil_f_trend), 
              col = 'red', lwd = 2, alpha = 0.8) +
    geom_line(data = Zn_soil_nf_trend, aes(x = Zn_soil_nf_seq, y = Zn_soil_nf_trend), 
              col = 'blue', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Zinc (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

(plots_soil_CN <- plot_C_soil + plot_CN_soil +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

(plots_soil <- (plot_Ca_soil + plot_K_soil + plot_P_soil) / (plot_Al_soil + plot_Zn_soil) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## soil characteristics
### pH
pH_lm <- lmer(pH ~ Elevation * Fire + (1 | Site), data = data)
# plot(resid(pH_lm) ~ fitted(pH_lm))
Anova(pH_lm)

### CEC
CEC_lm <- lmer(CEC ~ Elevation * Fire + (1 | Site), data = data)
#plot(resid(CEC_lm) ~ fitted(CEC_lm))
Anova(CEC_lm)

## soil characteristics
### retention
retention_lm <- lmer(asin(sqrt(0.01 * Retention)) ~ Elevation * Fire + (1 | Site), data = data)
# plot(resid(retention_lm) ~ fitted(retention_lm))
Anova(retention_lm)
summary(retention_lm) # N = 40

#### tables and posthoc ####

### topography
#### create table with mean latitude, longitude, elevation, slope, and aspect for each site
topography <- data %>% group_by(Site) %>% summarise_at(vars(Latitude, Longitude, Elevation, Slope, Aspect), mean, na.rm = TRUE)
write.csv(topography, "tables/topography.csv")

### allometry
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(Anova(slope_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(height_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(canopy_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(diam_lm)[, c(1, 2, 3)]),
                as.matrix(Anova(density_lm)[, c(1, 2, 3)])),
          'tables/allometry.csv')

### foliar organics
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(Anova(C_foliar_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(N_foliar_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(CN_foliar_lm)[, c(1, 2, 3)])),
          'tables/foliar_cn.csv')

### foliar inorganics
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(Anova(Ca_foliar_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(P_foliar_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(K_foliar_lm)[, c(1, 2, 3)]),
                as.matrix(Anova(Mg_foliar_lm)[, c(1, 2, 3)]),
                as.matrix(Anova(Al_foliar_lm)[, c(1, 2, 3)]),
                as.matrix(Anova(Zn_foliar_lm)[, c(1, 2, 3)])),
          'tables/foliar_inorganics.csv')

### foliar isotopes
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(Anova(d13C_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(d15N_lm)[, c(1, 2, 3)])),
          'tables/foliar_isotope.csv')

### soil organics
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(Anova(C_soil_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(N_soil_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(CN_soil_lm)[, c(1, 2, 3)]),
                as.matrix(Anova(retention_lm)[, c(1, 2, 3)])),
          'tables/soil_organics.csv')

(summary(emmeans(canopy_lm, ~elevation_fac))[1,2] - summary(emmeans(canopy_lm, ~elevation_fac))[2,2])/ summary(emmeans(canopy_lm, ~elevation_fac))[2,2]

(summary(emmeans(P_foliar_lm, ~fire))[1,2] - summary(emmeans(P_foliar_lm, ~fire))[2,2])/ summary(emmeans(P_foliar_lm, ~fire))[2,2]
(summary(emmeans(K_foliar_lm, ~fire))[1,2] - summary(emmeans(K_foliar_lm, ~fire))[2,2])/ summary(emmeans(K_foliar_lm, ~fire))[2,2]
(summary(emmeans(Ca_foliar_lm, ~elevation_fac))[1,2] - summary(emmeans(Ca_foliar_lm, ~elevation_fac))[2,2])/ summary(emmeans(Ca_foliar_lm, ~elevation_fac))[2,2]
(summary(emmeans(Zn_foliar_lm, ~elevation_fac))[1,2] - summary(emmeans(Zn_foliar_lm, ~elevation_fac))[2,2])/ summary(emmeans(Zn_foliar_lm, ~elevation_fac))[2,2]

(summary(emmeans(d13C_lm, ~elevation_fac))[1,2] - summary(emmeans(d13C_lm, ~elevation_fac))[2,2])/ summary(emmeans(d13C_lm, ~elevation_fac))[2,2]

(summary(emmeans(C_soil_lm, ~fire))[1,2] - summary(emmeans(C_soil_lm, ~fire))[2,2])/ summary(emmeans(C_soil_lm, ~fire))[2,2]
(summary(emmeans(C_soil_lm, ~elevation_fac))[1,2] - summary(emmeans(C_soil_lm, ~elevation_fac))[2,2])/ summary(emmeans(C_soil_lm, ~elevation_fac))[2,2]
(summary(emmeans(CN_soil_lm, ~elevation_fac))[1,2] - summary(emmeans(CN_soil_lm, ~elevation_fac))[2,2])/ summary(emmeans(CN_soil_lm, ~elevation_fac))[2,2]

### soil inorganics
write.csv(cbind(as.matrix(Anova(Ca_soil_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(P_soil_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(K_soil_lm)[, c(1, 2, 3)]),
                as.matrix(Anova(Mg_soil_lm)[, c(1, 2, 3)]),
                as.matrix(Anova(Al_soil_lm)[, c(1, 2, 3)]),
                as.matrix(Anova(Zn_soil_lm)[, c(1, 2, 3)])),
          'tables/soil_inorganics.csv')

(summary(emmeans(K_soil_lm, ~fire))[1,2] - summary(emmeans(K_soil_lm, ~fire))[2,2])/ summary(emmeans(K_soil_lm, ~fire))[2,2]
(summary(emmeans(Ca_soil_lm, ~elevation_fac))[1,2] - summary(emmeans(Ca_soil_lm, ~elevation_fac))[2,2])/ summary(emmeans(Ca_soil_lm, ~elevation_fac))[2,2]

## soil characteristics
write.csv(cbind(as.matrix(Anova(retention_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(pH_lm)[, c(1, 2, 3)]), 
                as.matrix(Anova(CEC_lm)[, c(1, 2, 3)])),
          'tables/soil_char.csv')

#### save graphs ####

## slope
ggsave("plots/plot_slope.jpeg", plot = plot_slope,
       width = 28, height = 18, units = "cm", dpi = 600) # 1 panel

## allometry
ggsave("plots/plots_allometry.jpeg", plot = plots_allometry,
       width = 28, height = 18, units = "cm", dpi = 600) # 4 panels

## foliar nutrients
ggsave("plots/plots_foliar.jpeg", plot = plots_foliar,
       width = 45, height = 25, units = "cm", dpi = 600) # 4 panels

## soil CN
ggsave("plots/plots_soil_CN.jpeg", plot = plots_soil_CN,
       width = 42, height = 25, units = "cm", dpi = 600) # 2 panels

## soil nutrients
ggsave("plots/plots_soil.jpeg", plot = plots_soil,
       width = 42, height = 25, units = "cm", dpi = 600) # 5 panels

