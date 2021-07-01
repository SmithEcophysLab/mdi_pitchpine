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
library(car)

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
slope_lm <- lm(Slope ~ Elevation * Fire, data = data)
# plot(resid(slope_lm) ~ fitted(slope_lm))
Anova(slope_lm)

slope_f_slope <- summary(emtrends(slope_lm, ~ Fire, var = "Elevation"))[1, 2] 
slope_f_intercept <- summary(emmeans(slope_lm, ~ Fire, at = list(Elevation = 0)))[1, 2] 
slope_f_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
slope_f_trend <- slope_f_intercept + slope_f_seq * slope_f_slope
slope_f_trend <- as.data.frame(cbind(slope_f_seq, slope_f_trend))

slope_nf_slope <- summary(emtrends(slope_lm, ~ Fire, var = "Elevation"))[2, 2] 
slope_nf_intercept <- summary(emmeans(slope_lm, ~ Fire, at = list(Elevation = 0)))[2, 2] 
slope_nf_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
slope_nf_trend <- slope_nf_intercept + slope_nf_seq * slope_nf_slope
slope_nf_trend <- as.data.frame(cbind(slope_nf_seq, slope_nf_trend))

(plot_slope <- ggplot(data = data, aes(x = Elevation, y = Slope)) +
        geom_jitter(aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        geom_line(data = slope_f_trend, aes(x = slope_f_seq, y = slope_f_trend), 
                  col = 'red', lwd = 2, alpha = 0.8) +
        geom_line(data = slope_nf_trend, aes(x = slope_nf_seq, y = slope_nf_trend), 
                  col = 'blue', lwd = 2, alpha = 0.8) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1500)) +
        scale_y_continuous(name = "Slope (◦)") +
        guides(color = guide_legend("Fire History")))

## allometry
### height
height_lm <- lm(log(Height) ~ Elevation * Fire, data = data)
# plot(resid(height_lm) ~ fitted(height_lm))
Anova(height_lm)
summary(height_lm) # N = 40

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
canopy_lm <- lm(log(Canopy) ~ Elevation * Fire , data = data)
# plot(resid(canopy_lm) ~ fitted(canopy_lm))
Anova(canopy_lm)

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
diam_lm <- lm(log(Diam) ~ Elevation * Fire , data = data)
#plot(resid(diam_lm) ~ fitted(diam_lm))
Anova(diam_lm)

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
    scale_y_continuous(name = "DBH (cm)") +
    guides(color = guide_legend("Fire History")))

### density
density_lm <- lm(mean_distance ~ Elevation * Fire , data = data)
#plot(resid(density_lm) ~ fitted(density_lm))
Anova(density_lm)

density_slope <- summary(emtrends(density_lm, ~ Fire, var = "Elevation"))[1, 2] 
density_intercept <- summary(emmeans(density_lm, ~ Fire, at = list(Elevation = 0)))[1, 2] 
density_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
density_trend <- density_intercept + density_seq * density_slope
density_trend <- as.data.frame(cbind(density_seq, density_trend))

(plot_density <- ggplot(data = data, aes(x = Elevation, y = mean_distance)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = density_trend, aes(x = density_seq, y = density_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
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
d13C_lm <- lm(d13C ~ Elevation * Fire , data = data)
# plot(resid(d13C_lm) ~ fitted(d13C_lm))
Anova(d13C_lm)

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
d15N_lm <- lm(d15N ~ Elevation * Fire , data = data)
#plot(resid(d15N_lm) ~ fitted(d15N_lm))
Anova(d15N_lm)

(plot_d15N <- ggplot(data = data, aes(x = Elevation, y = d15N)) +
        geom_jitter(aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
        ylab(expression(delta^{"15"}*"N (‰)")) +
        guides(color = guide_legend("Fire History")))

## foliar organics
### C_foliar
C_foliar_lm <- lm(C_foliar ~ Elevation * Fire , data = data)
#plot(resid(C_foliar_lm) ~ fitted(C_foliar_lm))
Anova(C_foliar_lm)

(plot_C_foliar <- ggplot(data = data, aes(x = Elevation, y = C_foliar)) +
        geom_jitter(height = 0, aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1500)) +
        ylab(expression("Foliar Carbon (g g"^{-1}*")")) +
        guides(color = guide_legend("Fire History")))

### N_foliar
N_foliar_lm <- lm(N_foliar ~ Elevation * Fire , data = data)
#plot(resid(N_foliar_lm) ~ fitted(N_foliar_lm))
Anova(N_foliar_lm)

(plot_N_foliar <- ggplot(data = data, aes(x = Elevation, y = N_foliar)) +
        geom_jitter(aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1500)) +
        ylab(expression("Foliar Nitrogen (g g"^{-1}*")")) +
        guides(color = guide_legend("Fire History")))

### CN_foliar
CN_foliar_lm <- lm(CN_foliar ~ Elevation * Fire , data = data)
#plot(resid(CN_foliar_lm) ~ fitted(CN_foliar_lm))
Anova(CN_foliar_lm)

(plot_CN_foliar <- ggplot(data = data, aes(x = Elevation, y = N_foliar)) +
        geom_jitter(aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1500)) +
        ylab("Foliar Carbon/Nitrogen") +
        guides(color = guide_legend("Fire History")))

## foliar inorganics
### Ca_foliar
Ca_foliar_lm <- lm(Ca_foliar ~ Elevation * Fire , data = data)
# plot(resid(Ca_foliar_lm) ~ fitted(Ca_foliar_lm))
Anova(Ca_foliar_lm)

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
P_foliar_lm <- lm(log(P_foliar) ~ Elevation * Fire , data = data)
#plot(resid(P_foliar_lm) ~ fitted(P_foliar_lm))
Anova(P_foliar_lm)

(plot_P_foliar <- ggplot(data = data, aes(x = Elevation, y = P_foliar)) +
        geom_jitter(height = 0, aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
        ylab(expression("Foliar Phosphorus (g g"^{-1}*")")) +
        guides(color = guide_legend("Fire History")))

### K_foliar
K_foliar_lm <- lm(log(K_foliar) ~ Elevation * Fire , data = data)
#plot(resid(K_foliar_lm) ~ fitted(K_foliar_lm))
Anova(K_foliar_lm)

# K_foliar_f_slope <- summary(emtrends(K_foliar_lm, ~ Fire, var = "Elevation"))[1, 2] 
# K_foliar_f_intercept <- summary(emmeans(K_foliar_lm, ~ Fire, at = list(Elevation = 0)))[1, 2] 
# K_foliar_f_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
# K_foliar_f_trend <- K_foliar_f_intercept + K_foliar_f_seq * K_foliar_f_slope
# K_foliar_f_trend <- as.data.frame(cbind(K_foliar_f_seq, K_foliar_f_trend))
# 
# K_foliar_nf_slope <- summary(emtrends(K_foliar_lm, ~ Fire, var = "Elevation"))[2, 2] 
# K_foliar_nf_intercept <- summary(emmeans(K_foliar_lm, ~ Fire, at = list(Elevation = 0)))[2, 2] 
# K_foliar_nf_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
# K_foliar_nf_trend <- K_foliar_nf_intercept + K_foliar_nf_seq * K_foliar_nf_slope
# K_foliar_nf_trend <- as.data.frame(cbind(K_foliar_nf_seq, K_foliar_nf_trend))

(plot_K_foliar <- ggplot(data = data, aes(x = Elevation, y = log(K_foliar))) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    # geom_line(data = K_foliar_f_trend, aes(x = K_foliar_f_seq, y = K_foliar_f_trend), 
    #           col = 'red', lwd = 2, alpha = 0.8) +
    # geom_line(data = K_foliar_nf_trend, aes(x = K_foliar_nf_seq, y = K_foliar_nf_trend), 
    #           col = 'blue', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Foliar Potassium (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### Mg_foliar
Mg_foliar_lm <- lm(Mg_foliar ~ Elevation * Fire , data = data)
#plot(resid(Mg_foliar_lm) ~ fitted(Mg_foliar_lm))
Anova(Mg_foliar_lm)

(plot_Mg_foliar <- ggplot(data = data, aes(x = Elevation, y = Mg_foliar)) +
        geom_jitter(aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
        ylab(expression("Foliar Magnesium (g g"^{-1}*")")) +
        guides(color = guide_legend("Fire History")))

### Al_foliar
Al_foliar_lm <- lm(Al_foliar ~ Elevation * Fire , data = data)
#plot(resid(Al_foliar_lm) ~ fitted(Al_foliar_lm))
Anova(Al_foliar_lm)

(plot_Al_foliar <- ggplot(data = data, aes(x = Elevation, y = Al_foliar)) +
        geom_jitter(aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
        ylab(expression("Foliar Aluminum (g g"^{-1}*")")) +
        guides(color = guide_legend("Fire History")))

### Zn_foliar
Zn_foliar_lm <- lm(log(Zn_foliar) ~ Elevation * Fire , data = data)
#plot(resid(Zn_foliar_lm) ~ fitted(Zn_foliar_lm))
Anova(Zn_foliar_lm)

Zn_foliar_slope <- summary(emtrends(Zn_foliar_lm, ~ Elevation, var = "Elevation"))[1, 2] 
Zn_foliar_intercept <- summary(emmeans(Zn_foliar_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
Zn_foliar_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
Zn_foliar_trend <- Zn_foliar_intercept + Zn_foliar_seq * Zn_foliar_slope
Zn_foliar_trend <- as.data.frame(cbind(Zn_foliar_seq, Zn_foliar_trend))

(plot_Zn_foliar <- ggplot(data = data, aes(x = Elevation, y = log(Zn_foliar))) +
        geom_jitter(aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        geom_line(data = Zn_foliar_trend, aes(x = Zn_foliar_seq, y = Zn_foliar_trend), 
                  col = 'black', lwd = 2, alpha = 0.8) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
        ylab(expression("Foliar Zinc (mg g"^{-1}*")")) +
        guides(color = guide_legend("Fire History")))

(plots_foliar_organic <- (plot_d13C + plot_d15N) / (plot_C_foliar + plot_N_foliar + plot_CN_foliar) + 
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

(plots_foliar_inorganic <- plot_Al_foliar + plot_Ca_foliar + plot_K_foliar + plot_Mg_foliar +
        plot_P_foliar + plot_Zn_foliar + 
        plot_layout(guides = 'collect') +
        plot_annotation(tag_levels = 'A') & 
        theme(plot.tag = element_text(size = 16)))

## soil organics
### C_soil
C_soil_lm <- lm(C_soil ~ Elevation * Fire , data = data)
#plot(resid(C_soil_lm) ~ fitted(C_soil_lm))
Anova(C_soil_lm)

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
N_soil_lm <- lm(N_soil ~ Elevation * Fire , data = data)
#plot(resid(N_soil_lm) ~ fitted(N_soil_lm))
Anova(N_soil_lm)

(plot_N_soil <- ggplot(data = data, aes(x = Elevation, y = N_soil)) +
        geom_jitter(height = 0, aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
        ylab(expression("Soil Nitrogen (g g"^{-1}*")")) +
        guides(color = guide_legend("Fire History")))

### CN_soil
CN_soil_lm <- lm(log(CN_soil) ~ Elevation * Fire , data = data)
#plot(resid(CN_soil_lm) ~ fitted(CN_soil_lm))
Anova(CN_soil_lm)

# CN_soil_slope <- summary(emtrends(CN_soil_lm, ~ Elevation, var = "Elevation"))[1, 2] 
# CN_soil_intercept <- summary(emmeans(CN_soil_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
# CN_soil_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
# CN_soil_trend <- CN_soil_intercept + CN_soil_seq * CN_soil_slope
# CN_soil_trend <- as.data.frame(cbind(CN_soil_seq, CN_soil_trend))

(plot_CN_soil <- ggplot(data = data, aes(x = Elevation, y = log(CN_soil))) +
    geom_jitter(height = 0, aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    # geom_line(data = CN_soil_trend, aes(x = CN_soil_seq, y = CN_soil_trend), 
    #           col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab("Soil Carbon/Nitrogen") +
    guides(color = guide_legend("Fire History")))

## soil inorganics
### Ca_soil
Ca_soil_lm <- lm(Ca_soil ~ Elevation * Fire , data = data)
#plot(resid(Ca_soil_lm) ~ fitted(Ca_soil_lm))
Anova(Ca_soil_lm)

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
P_soil_lm <- lm(log(P_soil) ~ Elevation * Fire , data = data)
#plot(resid(P_soil_lm) ~ fitted(P_soil_lm))
Anova(P_soil_lm)

# P_soil_slope <- summary(emtrends(P_soil_lm, ~ Elevation, var = "Elevation"))[1, 2] 
# P_soil_intercept <- summary(emmeans(P_soil_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
# P_soil_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
# P_soil_trend <- P_soil_intercept + P_soil_seq * P_soil_slope
# P_soil_trend <- as.data.frame(cbind(P_soil_seq, P_soil_trend))

(plot_P_soil <- ggplot(data = data, aes(x = Elevation, y = log(P_soil))) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    # geom_line(data = P_soil_trend, aes(x = P_soil_seq, y = P_soil_trend), 
    #           col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Phosphorus (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### K_soil
K_soil_lm <- lm(K_soil ~ Elevation * Fire , data = data)
#plot(resid(K_soil_lm) ~ fitted(K_soil_lm))
Anova(K_soil_lm)

# K_soil_slope <- summary(emtrends(K_soil_lm, ~ Elevation, var = "Elevation"))[1, 2] 
# K_soil_intercept <- summary(emmeans(K_soil_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
# K_soil_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
# K_soil_trend <- K_soil_intercept + K_soil_seq * K_soil_slope
# K_soil_trend <- as.data.frame(cbind(K_soil_seq, K_soil_trend))

(plot_K_soil <- ggplot(data = data, aes(x = Elevation, y = K_soil)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    # geom_line(data = K_soil_trend, aes(x = K_soil_seq, y = K_soil_trend), 
    #           col = 'black', lwd = 2, alpha = 0.8) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Potassium (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### Mg_soil
Mg_soil_lm <- lm(Mg_soil ~ Elevation * Fire , data = data)
#plot(resid(Mg_soil_lm) ~ fitted(Mg_soil_lm))
Anova(Mg_soil_lm)

# Mg_soil_slope <- summary(emtrends(Mg_soil_lm, ~ Elevation, var = "Elevation"))[1, 2] 
# Mg_soil_intercept <- summary(emmeans(Mg_soil_lm, ~ Elevation, at = list(Elevation = 0)))[1, 2] 
# Mg_soil_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
# Mg_soil_trend <- Mg_soil_intercept + Mg_soil_seq * Mg_soil_slope
# Mg_soil_trend <- as.data.frame(cbind(Mg_soil_seq, Mg_soil_trend))

(plot_Mg_soil <- ggplot(data = data, aes(x = Elevation, y = Mg_soil)) +
        geom_jitter(aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        # geom_line(data = Mg_soil_trend, aes(x = Mg_soil_seq, y = Mg_soil_trend), 
        #           col = 'black', lwd = 2, alpha = 0.8) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
        ylab(expression("Soil Magnesium (mg g"^{-1}*")")) +
        guides(color = guide_legend("Fire History")))

### Al_soil
Al_soil_lm <- lm(log(Al_soil) ~ Elevation * Fire , data = data)
#plot(resid(Al_soil_lm) ~ fitted(Al_soil_lm))
Anova(Al_soil_lm)

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
Zn_soil_lm <- lm(log(Zn_soil) ~ Elevation * Fire , data = data)
#plot(resid(Zn_soil_lm) ~ fitted(Zn_soil_lm))
Anova(Zn_soil_lm)

(plot_Zn_soil <- ggplot(data = data, aes(x = Elevation, y = log(Zn_soil))) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Zinc (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### retention
retention_lm <- lm(asin(sqrt(0.01 * Retention)) ~ Elevation * Fire , data = data)
# plot(resid(retention_lm) ~ fitted(retention_lm))
Anova(retention_lm)

retention_f_slope <- summary(emtrends(retention_lm, ~ Fire, var = "Elevation"))[1, 2] 
retention_f_intercept <- summary(emmeans(retention_lm, ~ Fire, at = list(Elevation = 0)))[1, 2] 
retention_f_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
retention_f_trend <- retention_f_intercept + retention_f_seq * retention_f_slope
retention_f_trend <- as.data.frame(cbind(retention_f_seq, retention_f_trend))

retention_nf_slope <- summary(emtrends(retention_lm, ~ Fire, var = "Elevation"))[2, 2] 
retention_nf_intercept <- summary(emmeans(retention_lm, ~ Fire, at = list(Elevation = 0)))[2, 2] 
retention_nf_seq <- seq(min(data$Elevation, na.rm = T), max(data$Elevation, na.rm = T), 0.01)
retention_nf_trend <- retention_nf_intercept + retention_nf_seq * retention_nf_slope
retention_nf_trend <- as.data.frame(cbind(retention_nf_seq, retention_nf_trend))

(plot_retention <- ggplot(data = data, aes(x = Elevation, y = asin(sqrt(0.01 * Retention)))) +
        geom_jitter(aes(color = Fire), size = 2) +
        scale_color_manual(values = c('red', 'blue')) +
        geom_line(data = retention_f_trend, aes(x = retention_f_seq, y = retention_f_trend), 
                  col = 'red', lwd = 2, alpha = 0.8) +
        geom_line(data = retention_nf_trend, aes(x = retention_nf_seq, y = retention_nf_trend), 
                  col = 'blue', lwd = 2, alpha = 0.8) +
        theme_few(base_size = 16) + 
        scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
        ylab("Soil Water Retention (%)") +
        guides(color = guide_legend("Fire History")))

(plots_soil_organics <- plot_C_soil + plot_N_soil + plot_CN_soil + plot_retention +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

(plots_soil_inorganics <- plot_Al_soil + plot_Ca_soil + plot_K_soil + plot_Mg_soil + plot_P_soil + plot_Zn_soil +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))


#### tables and posthoc ####

### topography
#### create table with mean latitude, longitude, elevation, slope, and aspect for each site
topography <- data %>% group_by(Site) %>% summarise_at(vars(Latitude, Longitude, Elevation, Slope, Aspect), mean, na.rm = TRUE)
# write.csv(topography, "tables/topography.csv")

### allometry
#### create table with degrees of f reedom, f-value, p-value results from linear models
# write.csv(cbind(as.matrix(as.data.frame(Anova(canopy_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(diam_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(density_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(slope_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(height_lm))[, c(2:4)])),
#           'tables/allometry.csv')

### foliar organics
#### create table with degrees of f reedom, f-value, p-value results from linear models
# write.csv(cbind(as.matrix(as.data.frame(Anova(d13C_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(d15N_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(C_foliar_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(N_foliar_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(CN_foliar_lm))[, c(2:4)])),
#           'tables/foliar_organics.csv')

### foliar inorganics
#### create table with degrees of f reedom, f-value, p-value results from linear models
# write.csv(cbind(as.matrix(as.data.frame(Anova(Al_foliar_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(Ca_foliar_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(K_foliar_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(Mg_foliar_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(P_foliar_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(Zn_foliar_lm))[, c(2:4)])),
#           'tables/foliar_inorganics.csv')

### soil organics
#### create table with degrees of f reedom, f-value, p-value results from linear models
# write.csv(cbind(as.matrix(as.data.frame(Anova(C_soil_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(N_soil_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(CN_soil_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(retention_lm))[, c(2:4)])),
#           'tables/soil_organics.csv')

# (summary(emmeans(canopy_lm, ~elevation_fac))[1,2] - summary(emmeans(canopy_lm, ~elevation_fac))[2,2])/ summary(emmeans(canopy_lm, ~elevation_fac))[2,2]
# 
# (summary(emmeans(P_foliar_lm, ~fire))[1,2] - summary(emmeans(P_foliar_lm, ~fire))[2,2])/ summary(emmeans(P_foliar_lm, ~fire))[2,2]
# (summary(emmeans(K_foliar_lm, ~fire))[1,2] - summary(emmeans(K_foliar_lm, ~fire))[2,2])/ summary(emmeans(K_foliar_lm, ~fire))[2,2]
# (summary(emmeans(Ca_foliar_lm, ~elevation_fac))[1,2] - summary(emmeans(Ca_foliar_lm, ~elevation_fac))[2,2])/ summary(emmeans(Ca_foliar_lm, ~elevation_fac))[2,2]
# (summary(emmeans(Zn_foliar_lm, ~elevation_fac))[1,2] - summary(emmeans(Zn_foliar_lm, ~elevation_fac))[2,2])/ summary(emmeans(Zn_foliar_lm, ~elevation_fac))[2,2]
# 
# (summary(emmeans(d13C_lm, ~elevation_fac))[1,2] - summary(emmeans(d13C_lm, ~elevation_fac))[2,2])/ summary(emmeans(d13C_lm, ~elevation_fac))[2,2]
# 
# (summary(emmeans(C_soil_lm, ~fire))[1,2] - summary(emmeans(C_soil_lm, ~fire))[2,2])/ summary(emmeans(C_soil_lm, ~fire))[2,2]
# (summary(emmeans(C_soil_lm, ~elevation_fac))[1,2] - summary(emmeans(C_soil_lm, ~elevation_fac))[2,2])/ summary(emmeans(C_soil_lm, ~elevation_fac))[2,2]
# (summary(emmeans(CN_soil_lm, ~elevation_fac))[1,2] - summary(emmeans(CN_soil_lm, ~elevation_fac))[2,2])/ summary(emmeans(CN_soil_lm, ~elevation_fac))[2,2]

### soil inorganics
# write.csv(cbind(as.matrix(as.data.frame(Anova(Al_soil_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(Ca_soil_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(K_soil_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(Mg_soil_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(P_soil_lm))[, c(2:4)]),
#                 as.matrix(as.data.frame(Anova(Zn_soil_lm))[, c(2:4)])),
#           'tables/soil_inorganics.csv')

# (summary(emmeans(K_soil_lm, ~fire))[1,2] - summary(emmeans(K_soil_lm, ~fire))[2,2])/ summary(emmeans(K_soil_lm, ~fire))[2,2]
# (summary(emmeans(Ca_soil_lm, ~elevation_fac))[1,2] - summary(emmeans(Ca_soil_lm, ~elevation_fac))[2,2])/ summary(emmeans(Ca_soil_lm, ~elevation_fac))[2,2]

#### save graphs ####

## slope
ggsave("plots/plot_slope.jpeg", plot = plot_slope,
       width = 28, height = 18, units = "cm", dpi = 600) # 1 panel

## allometry
ggsave("plots/plots_allometry.jpeg", plot = plots_allometry,
       width = 28, height = 18, units = "cm", dpi = 600) # 4 panels

## foliar nutrients
ggsave("plots/plots_foliar_organics.jpeg", plot = plots_foliar_organic,
       width = 45, height = 25, units = "cm", dpi = 600) # 5 panels
ggsave("plots/plots_foliar_inorganics.jpeg", plot = plots_foliar_inorganic,
       width = 45, height = 25, units = "cm", dpi = 600) # 6 panels

## soil nutrients
ggsave("plots/plots_soil_organics.jpeg", plot = plots_soil_organics,
       width = 42, height = 25, units = "cm", dpi = 600) # 4 panels
ggsave("plots/plots_soil_inorganics.jpeg", plot = plots_soil_inorganics,
       width = 42, height = 25, units = "cm", dpi = 600) # 6 panels

