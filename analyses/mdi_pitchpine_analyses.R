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

# function to get pairwise letters from Tukey's HSD for boxplots
letters <- function(df = data, x = "Site", dfy, y){
  .checkvar <- sym(y) #convert string to variable
  .groupvar <- sym(x) #convert string to variable
  
  # find highest value of y data
  abs_max <- max(dfy, na.rm = TRUE) 
  # get the highest point for each species
  max_y <- df %>% dplyr::group_by(!! .groupvar) %>% 
    dplyr::summarise(yaxis = max(!! .checkvar, na.rm = TRUE) + 0.05 * abs_max)
    #"!!" unquotes the variable to allow for evaluation within dplyr
  # get Tukey HSD results
  hsd <- HSD.test(aov(as.formula(paste(y, x, sep = "~")), df), x, group = TRUE) 
  # add Tukey HSD results to dataframe containing graphing positions
  group <- as.data.frame(hsd$groups)
  group <- tibble::rownames_to_column(group, "Site")
  letters <- dplyr::full_join(max_y, group, by = "Site")
  return(letters)
}

# function to get pairwise letters from Tukey's HSD with a transformed y variable for boxplots
letters_adj <- function(df = data, x = "Site", dfy, y, adjy){
  .checkvar <- sym(y)
  .groupvar <- sym(x)
  
  abs_max <- max(dfy, na.rm = TRUE)
  max_y <- df %>% dplyr::group_by(!! .groupvar) %>% 
    dplyr::summarise(yaxis = max(!! .checkvar, na.rm = TRUE) + 0.05 * abs_max) 
  hsd <- HSD.test(aov(as.formula(paste(adjy, x, sep = "~")), df), x, group = TRUE) 
  group <- as.data.frame(hsd$groups)
  group <- tibble::rownames_to_column(group, "Site")
  letters <- dplyr::full_join(max_y, group, by = "Site")
  return(letters)
}

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

## create a generic variable set to pass to formula argument
ind_variables <- c('Elevation', 'Slope', 'Fire')
dep_variables <- c("log(Height)", "log(Canopy)", "log(Diam)", "mean_distance",
                  "d13C", "d15N", "C_foliar", "N_foliar", "CN_foliar", "Ca_foliar", "log(P_foliar)",
                  "log(K_foliar)", "Mg_foliar", "Al_foliar", "log(Zn_foliar)", 
                  "Ca_soil", "log(P_soil)", "K_soil", "Mg_soil", "log(Al_soil)", "log(Zn_soil)", 
                  "pH", "CEC", "C_soil", "N_soil", "log(CN_soil)", "asin(sqrt(0.01 * Retention))")

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

## allometry
### height
height_lm <- lm(as.formula(paste(dep_variables[1],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
#plot(resid(height_lm) ~ fitted(height_lm))
Anova(height_lm)

### canopy
canopy_lm <- lm(as.formula(paste(dep_variables[2],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
# plot(resid(canopy_lm) ~ fitted(canopy_lm))
anova(canopy_lm)

(plot_canopy <- ggplot(data = data, aes(x = Elevation, y = Canopy)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    theme(legend.position = "bottom") +
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    scale_y_continuous(name = "Canopy Spread (m)") +
    guides(color = "none"))

### diam
diam_lm <- lm(as.formula(paste(dep_variables[3],
                              paste(ind_variables, collapse = "*"),
                              sep = "~")), data = data)
#plot(resid(diam_lm) ~ fitted(diam_lm))
anova(diam_lm)

(plot_diam <- ggplot(data = data, aes(x = Elevation, y = Diam)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) +
    theme(legend.position = "bottom") +
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    scale_y_continuous(name = "DBH (cm)", limits = c(0, 60)) +
    guides(color = guide_legend("Fire History")))

### density
density_lm <- lm(as.formula(paste(dep_variables[4],
                                  paste(ind_variables, collapse = "*"),
                                  sep = "~")), data = data)
#plot(resid(density_lm) ~ fitted(density_lm))
anova(density_lm)

(plot_density <- ggplot(data = data, aes(x = Elevation, y = mean_distance)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) +
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1500)) +
    scale_y_continuous(name = "Distance Between Neighbors (m)", limits = c(0, 6)) + 
    guides(color = "none"))

(plots_allometry <- plot_canopy + plot_diam + plot_density + 
    plot_layout(guides = 'keep') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## foliar isotopes
### d13C
d13C_lm <- lm(as.formula(paste(dep_variables[5],
                               paste(ind_variables, collapse = "*"),
                               sep = "~")), data = data)
#plot(resid(d13C_lm) ~ fitted(d13C_lm))
anova(d13C_lm)

(plot_d13C_elev <- ggplot(data = data, aes(x = Elevation, y = d13C)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression(delta^{"13"}*"C (‰)")) +
    guides(color = guide_legend("Fire History")))

d13C_slope_fire <- summary(emtrends(d13C_lm, ~ Fire, var = "Slope"))[1, 2] # slope = -0.1200
d13C_intercept_fire <- summary(emmeans(d13C_lm, ~ Fire, at = list(Slope = 0)))[1, 2] # intercept = -25.5
d13C_seq_fire <- seq(min(data$Slope, na.rm = T), max(data$Slope, na.rm = T), 0.01)
d13C_trend_fire <- d13C_intercept_fire + d13C_seq_fire * d13C_slope_fire
d13C_trend_fire <- as.data.frame(cbind(d13C_seq_fire, d13C_trend_fire))

d13C_slope_nofire <- summary(emtrends(d13C_lm, ~ Fire, var = "Slope"))[2, 2] # slope = 0.0778
d13C_intercept_nofire <- summary(emmeans(d13C_lm, ~ Fire, at = list(Slope = 0)))[2, 2] # intercept = -29.3
d13C_seq_nofire <- seq(min(data$Slope, na.rm = T), max(data$Slope, na.rm = T), 0.01)
d13C_trend_nofire <- d13C_intercept_nofire + d13C_seq_nofire * d13C_slope_nofire
d13C_trend_nofire <- as.data.frame(cbind(d13C_seq_nofire, d13C_trend_nofire))

(plot_d13C_slope <- ggplot(data = data, aes(x = Slope, y = d13C)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line(data = d13C_trend_fire, aes(x = d13C_seq_fire, y = d13C_trend_fire), 
              col = 'red', lwd = 2, alpha = 0.5) +
    geom_line(data = d13C_trend_nofire, aes(x = d13C_seq_nofire, y = d13C_trend_nofire), 
              col = 'blue', lwd = 2, alpha = 0.5) +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Slope (°)") +
    ylab(expression(delta^{"13"}*"C (‰)")) +
    guides(color = guide_legend("Fire History")))

### d15N
d15N_lm <- lm(as.formula(paste(dep_variables[6],
                               paste(ind_variables, collapse = "*"),
                               sep = "~")), data = data)
#plot(resid(d15N_lm) ~ fitted(d15N_lm))
anova(d15N_lm)

(plots_isotopes <- plot_d13C_elev + plot_d13C_slope +
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## foliar organics
### C_foliar
C_foliar_lm <- lm(as.formula(paste(dep_variables[7],
                                  paste(ind_variables, collapse = "*"),
                                  sep = "~")), data = data)
#plot(resid(C_foliar_lm) ~ fitted(C_foliar_lm))
anova(C_foliar_lm)

(plot_C_foliar <- ggplot(data = data, aes(x = Fire, y = C_foliar)) +
    geom_jitter(height = 0, aes(color = Fire), size = 2) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    scale_color_manual(values = c('red', 'blue')) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Fire History") +
    ylab(expression("Foliar Carbon (g g"^{-1}*")")) +
    guides(color = "none"))

### N_foliar
N_foliar_lm <- lm(as.formula(paste(dep_variables[8],
                                  paste(ind_variables, collapse = "*"),
                                  sep = "~")), data = data)
#plot(resid(N_foliar_lm) ~ fitted(N_foliar_lm))
anova(N_foliar_lm)

### CN_foliar
CN_foliar_lm <- lm(as.formula(paste(dep_variables[9],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(CN_foliar_lm) ~ fitted(CN_foliar_lm))
anova(CN_foliar_lm)

## foliar inorganics
### Ca_foliar
Ca_foliar_lm <- lm(as.formula(paste(dep_variables[10],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(Ca_foliar_lm) ~ fitted(Ca_foliar_lm))
anova(Ca_foliar_lm)

(plot_Ca_foliar <- ggplot(data = data, aes(x = Elevation, y = Ca_foliar)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Foliar Calcium (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### P_foliar
P_foliar_lm <- lm(as.formula(paste(dep_variables[11],
                                  paste(ind_variables, collapse = "*"),
                                  sep = "~")), data = data)
#plot(resid(P_foliar_lm) ~ fitted(P_foliar_lm))
anova(P_foliar_lm)

(plot_P_foliar <- ggplot(data = data, aes(x = Fire, y = P_foliar)) +
    geom_jitter(height = 0, aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Fire History") +
    ylab(expression("Foliar Phosphorus (g g"^{-1}*")")) +
    guides(color = "none"))

### K_foliar
K_foliar_lm <- lm(as.formula(paste(dep_variables[12],
                                  paste(ind_variables, collapse = "*"),
                                  sep = "~")), data = data)
#plot(resid(K_foliar_lm) ~ fitted(K_foliar_lm))
anova(K_foliar_lm)

ggplot(data, aes(x = Elevation, y = K_foliar)) + geom_jitter() + geom_smooth(method = lm)

(plot_K_foliar <- ggplot(data = data, aes(x = Fire, y = K_foliar)) +
    geom_jitter(height = 0, aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    theme_few(base_size = 16) +
    scale_x_discrete(name = "Fire History") +
    ylab(expression("Foliar Potassium (mg g"^{-1}*")")) +
    guides(color = "none"))

### Mg_foliar
Mg_foliar_lm <- lm(as.formula(paste(dep_variables[13],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(Mg_foliar_lm) ~ fitted(Mg_foliar_lm))
anova(Mg_foliar_lm)

### Al_foliar
Al_foliar_lm <- lm(as.formula(paste(dep_variables[14],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(Al_foliar_lm) ~ fitted(Al_foliar_lm))
anova(Al_foliar_lm)

### Zn_foliar
Zn_foliar_lm <- lm(as.formula(paste(dep_variables[15],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(Zn_foliar_lm) ~ fitted(Zn_foliar_lm))
anova(Zn_foliar_lm)

(plot_Zn_foliar <- ggplot(data = data, aes(x = Elevation, y = Zn_foliar)) +
    geom_jitter(height = 0, aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Foliar Zinc (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

(plots_foliar <- (plot_C_foliar + plot_P_foliar + plot_K_foliar) /
    (plot_Ca_foliar + plot_Zn_foliar) + 
    plot_layout(guides = 'collect') +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## soil organics
### C_soil
C_soil_lm <- lm(as.formula(paste(dep_variables[24],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
#plot(resid(C_soil_lm) ~ fitted(C_soil_lm))
anova(C_soil_lm)

(plot_C_soil <- ggplot(data = data, aes(x = Elevation, y = C_soil)) +
    geom_jitter(height = 0, aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Carbon (g g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

### N_soil
N_soil_lm <- lm(as.formula(paste(dep_variables[25],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
#plot(resid(N_soil_lm) ~ fitted(N_soil_lm))
anova(N_soil_lm)

### CN_soil
CN_soil_lm <- lm(as.formula(paste(dep_variables[26],
                                 paste(ind_variables, collapse = "*"),
                                 sep = "~")), data = data)
#plot(resid(CN_soil_lm) ~ fitted(CN_soil_lm))
anova(CN_soil_lm)

(plot_CN_soil_elev <- ggplot(data = data, aes(x = Elevation, y = CN_soil)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    scale_y_continuous(name = "Soil Carbon/Nitrogen") +
    guides(color = guide_legend("Fire History")))

(plot_CN_soil_slope <- ggplot(data = data, aes(x = Slope, y = CN_soil)) +
    geom_jitter(size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Slope ()") +
    scale_y_continuous(name = "Soil Carbon/Nitrogen"))

## soil inorganics
### Ca_soil
Ca_soil_lm <- lm(as.formula(paste(dep_variables[16],
                                 paste(ind_variables, collapse = "*"),
                                 sep = "~")), data = data)
#plot(resid(Ca_soil_lm) ~ fitted(Ca_soil_lm))
anova(Ca_soil_lm)

(plot_Ca_soil_elev <- ggplot(data = data, aes(x = Elevation, y = Ca_soil)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Calcium (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

# Ca_slope_fire <- summary(emtrends(Ca_soil_lm, ~ Fire, var = "Slope"))[1, 2]
# Ca_intercept_fire <- summary(emmeans(Ca_soil_lm, ~ Fire, at = list(Slope = 0)))[1, 2]
# Ca_seq_fire <- seq(min(data$Slope, na.rm = T), max(data$Slope, na.rm = T), 0.01)
# Ca_trend_fire <- Ca_intercept_fire + Ca_seq_fire * Ca_slope_fire
# Ca_trend_fire <- as.data.frame(cbind(Ca_seq_fire, Ca_trend_fire))
# 
# Ca_slope_nofire <- summary(emtrends(Ca_soil_lm, ~ Fire, var = "Slope"))[2, 2]
# Ca_intercept_nofire <- summary(emmeans(Ca_soil_lm, ~ Fire, at = list(Slope = 0)))[2, 2]
# Ca_seq_nofire <- seq(min(data$Slope, na.rm = T), max(data$Slope, na.rm = T), 0.01)
# Ca_trend_nofire <- Ca_intercept_nofire + Ca_seq_nofire * Ca_slope_nofire
# Ca_trend_nofire <- as.data.frame(cbind(Ca_seq_nofire, Ca_trend_nofire))
# 
# (plot_Ca_soil_slope <- ggplot(data = data, aes(x = Slope, y = Ca_soil)) +
#     geom_jitter(aes(color = Fire), size = 2) +
#     scale_color_manual(values = c('red', 'blue')) +
#     geom_line(data = Ca_trend_fire, aes(x = Ca_seq_fire, y = Ca_trend_fire), 
#               col = 'red', lwd = 2, alpha = 0.5) +
#     geom_line(data = Ca_trend_nofire, aes(x = Ca_seq_nofire, y = Ca_trend_nofire), 
#               col = 'blue', lwd = 2, alpha = 0.5) +
#     theme_few(base_size = 16) + 
#     scale_x_continuous(name = "Slope ()") +
#     ylab(expression("Soil Calcium (mg g"^{-1}*")")) +
#     guides(color = guide_legend("Fire History")))

### P_soil
P_soil_lm <- lm(as.formula(paste(dep_variables[17],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
#plot(resid(P_soil_lm) ~ fitted(P_soil_lm))
anova(P_soil_lm)

(plot_P_soil_elev <- ggplot(data = data, aes(x = Elevation, y = P_soil)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Phosphorus (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

# P_slope_fire <- summary(emtrends(P_soil_lm, ~ Fire, var = "Slope"))[1, 2] # slope = 0.2066
# P_intercept_fire <- summary(emmeans(P_soil_lm, ~ Fire, at = list(Slope = 0)))[1, 2] # intercept = -10.4
# P_seq_fire <- seq(min(data$Slope, na.rm = T), max(data$Slope, na.rm = T), 0.01)
# P_trend_fire <- P_intercept_fire + P_seq_fire * P_slope_fire
# P_trend_fire <- as.data.frame(cbind(P_seq_fire, P_trend_fire))
# 
# P_slope_nofire <- summary(emtrends(P_soil_lm, ~ Fire, var = "Slope"))[2, 2] # slope = -0.0287
# P_intercept_nofire <- summary(emmeans(P_soil_lm, ~ Fire, at = list(Slope = 0)))[2, 2] # intercept = -5.59
# P_seq_nofire <- seq(min(data$Slope, na.rm = T), max(data$Slope, na.rm = T), 0.01)
# P_trend_nofire <- P_intercept_nofire + P_seq_nofire * P_slope_nofire
# P_trend_nofire <- as.data.frame(cbind(P_seq_nofire, P_trend_nofire))
# 
# (plot_P_soil_slope <- ggplot(data = data, aes(x = Slope, y = P_soil)) +
#     geom_jitter(aes(color = Fire), size = 2) +
#     scale_color_manual(values = c('red', 'blue')) +
#     geom_abline(slope = 0.2066, intercept = -10.4) +
#     geom_abline(slope = -0.0287, intercept = -5.59) +
#     # geom_line(data = P_trend_fire, aes(x = P_seq_fire, y = P_trend_fire), 
#     #           col = 'red', lwd = 2, alpha = 0.5) +
#     # geom_line(data = P_trend_nofire, aes(x = P_seq_nofire, y = P_trend_nofire), 
#     #           col = 'blue', lwd = 2, alpha = 0.5) +
#     theme_few(base_size = 16) + 
#     scale_x_continuous(name = "Slope ()") +
#     ylab(expression("Soil Phosphorus (mg g"^{-1}*")")) +
#     guides(color = guide_legend("Fire History")))

### K_soil
K_soil_lm <- lm(as.formula(paste(dep_variables[18],
                                paste(ind_variables, collapse = "*"),
                                sep = "~")), data = data)
#plot(resid(K_soil_lm) ~ fitted(K_soil_lm))
anova(K_soil_lm)

(plot_K_soil_elev <- ggplot(data = data, aes(x = Elevation, y = K_soil)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Potassium (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

(plot_K_soil_fire <- ggplot(data = data, aes(x = Fire, y = K_soil)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    theme_few(base_size = 16) + 
    scale_x_discrete(name = "Fire History") +
    ylab(expression("Soil Potassium (mg g"^{-1}*")")) +
    guides(color = "none"))

ggplot(data, aes(x = Fire, y = K_soil)) + geom_boxplot()

### Mg_soil
Mg_soil_lm <- lm(as.formula(paste(dep_variables[19],
                                 paste(ind_variables, collapse = "*"),
                                 sep = "~")), data = data)
#plot(resid(Mg_soil_lm) ~ fitted(Mg_soil_lm))
anova(Mg_soil_lm)

(plot_Mg_soil_elev <- ggplot(data = data, aes(x = Elevation, y = Mg_soil)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Magnesium (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

# Mg_slope_fire <- summary(emtrends(Mg_soil_lm, ~ Fire, var = "Slope"))[1, 2] # slope = 0.03235
# Mg_intercept_fire <- summary(emmeans(Mg_soil_lm, ~ Fire, at = list(Slope = 0)))[1, 2] # intercept = -0.517
# Mg_seq_fire <- seq(min(data$Slope, na.rm = T), max(data$Slope, na.rm = T), 0.01)
# Mg_trend_fire <- Mg_intercept_fire + Mg_seq_fire * Mg_slope_fire
# Mg_trend_fire <- as.data.frame(cbind(Mg_seq_fire, Mg_trend_fire))
# 
# Mg_slope_nofire <- summary(emtrends(Mg_soil_lm, ~ Fire, var = "Slope"))[2, 2] # slope = -0.00979
# Mg_intercept_nofire <- summary(emmeans(Mg_soil_lm, ~ Fire, at = list(Slope = 0)))[2, 2] # intercept = 0.243
# Mg_seq_nofire <- seq(min(data$Slope, na.rm = T), max(data$Slope, na.rm = T), 0.01)
# Mg_trend_nofire <- Mg_intercept_nofire + Mg_seq_nofire * Mg_slope_nofire
# Mg_trend_nofire <- as.data.frame(cbind(Mg_seq_nofire, Mg_trend_nofire))
# 
# (plot_Mg_soil_slope <- ggplot(data = data, aes(x = Slope, y = Mg_soil)) +
#     geom_jitter(aes(color = Fire), size = 2) +
#     scale_color_manual(values = c('red', 'blue')) +
#     geom_line(data = Mg_trend_fire, aes(x = Mg_seq_fire, y = Mg_trend_fire),
#               col = 'red', lwd = 2, alpha = 0.5) +
#     geom_line(data = Mg_trend_nofire, aes(x = Mg_seq_nofire, y = Mg_trend_nofire),
#               col = 'blue', lwd = 2, alpha = 0.5) +
#     theme_few(base_size = 16) + 
#     scale_x_continuous(name = "Slope ()") +
#     ylab(expression("Soil Magnesium (mg g"^{-1}*")")) +
#     guides(color = guide_legend("Fire History")))

### Al_soil
Al_soil_lm <- lm(as.formula(paste(dep_variables[20],
                                 paste(ind_variables, collapse = "*"),
                                 sep = "~")), data = data)
#plot(resid(Al_soil_lm) ~ fitted(Al_soil_lm))
anova(Al_soil_lm)

(plot_Al_soil_elev <- ggplot(data = data, aes(x = Elevation, y = Al_soil)) +
    geom_jitter(aes(color = Fire), size = 2) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Elevation (m)", limits = c(0, 1000)) +
    ylab(expression("Soil Aluminum (mg g"^{-1}*")")) +
    guides(color = guide_legend("Fire History")))

# Al_slope_fire <- summary(emtrends(Al_soil_lm, ~ Fire, var = "Slope"))[1, 2] # slope = -0.1056
# Al_intercept_fire <- summary(emmeans(Al_soil_lm, ~ Fire, at = list(Slope = 0)))[1, 2] # intercept = -0.0278
# Al_seq_fire <- seq(min(data$Slope, na.rm = T), max(data$Slope, na.rm = T), 0.01)
# Al_trend_fire <- Al_intercept_fire + Al_seq_fire * Al_slope_fire
# Al_trend_fire <- as.data.frame(cbind(Al_seq_fire, Al_trend_fire))
# 
# Al_slope_nofire <- summary(emtrends(Al_soil_lm, ~ Fire, var = "Slope"))[2, 2] # slope = -0.0247
# Al_intercept_nofire <- summary(emmeans(Al_soil_lm, ~ Fire, at = list(Slope = 0)))[2, 2] # intercept = -2.0219
# Al_seq_nofire <- seq(min(data$Slope, na.rm = T), max(data$Slope, na.rm = T), 0.01)
# Al_trend_nofire <- Al_intercept_nofire + Al_seq_nofire * Al_slope_nofire
# Al_trend_nofire <- as.data.frame(cbind(Al_seq_nofire, Al_trend_nofire))
# 
# (plot_Al_soil_slope <- ggplot(data = data, aes(x = Slope, y = Al_soil)) +
#     geom_jitter(aes(color = Fire), size = 2) +
#     scale_color_manual(values = c('red', 'blue')) +
#     geom_line(data = Al_trend_fire, aes(x = Al_seq_fire, y = Al_trend_fire),
#               col = 'red', lwd = 2, alpha = 0.5) +
#     geom_line(data = Al_trend_nofire, aes(x = Al_seq_nofire, y = Al_trend_nofire),
#               col = 'blue', lwd = 2, alpha = 0.5) +
#     theme_few(base_size = 16) + 
#     scale_x_continuous(name = "Slope ()") +
#     ylab(expression("Soil Aluminum (mg g"^{-1}*")")) +
#     guides(color = guide_legend("Fire History")))

### Zn_soil
Zn_soil_lm <- lm(as.formula(paste(dep_variables[21],
                                 paste(ind_variables, collapse = "*"),
                                 sep = "~")), data = data)
#plot(resid(Zn_soil_lm) ~ fitted(Zn_soil_lm))
anova(Zn_soil_lm)

(plots_soil_elev <- plot_C_soil + plot_Ca_soil_elev + 
    plot_P_soil_elev + plot_K_soil_elev + plot_Mg_soil_elev + plot_Al_soil_elev +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 16)))

## soil characteristics
### pH
pH_lm <- lm(as.formula(paste(dep_variables[22],
                            paste(ind_variables, collapse = "*"),
                            sep = "~")), data = data)
#plot(resid(pH_lm) ~ fitted(pH_lm))
anova(pH_lm)

### CEC
CEC_lm <- lm(as.formula(paste(dep_variables[23],
                             paste(ind_variables, collapse = "*"),
                             sep = "~")), data = data)
#plot(resid(CEC_lm) ~ fitted(CEC_lm))
anova(CEC_lm)

## soil characteristics
### retention
retention_lm <- lm(as.formula(paste(dep_variables[27],
                                   paste(ind_variables, collapse = "*"),
                                   sep = "~")), data = data)
#plot(resid(retention_lm) ~ fitted(retention_lm))
anova(retention_lm)

(plot_retention <- ggplot(data = data, aes(x = Slope, y = Retention)) +
    geom_jitter(size = 2) +
    geom_smooth(method = lm, color = "black") +
    theme_few(base_size = 16) + 
    scale_x_continuous(name = "Slope (°)") +
    scale_y_continuous(name = "Soil Water Retention (%)"))

#### tables and posthoc ####

### topography
#### create table with mean latitude, longitude, elevation, slope, and aspect for each site
topography <- data %>% group_by(Site) %>% summarise_at(vars(Latitude, Longitude, Elevation, Slope, Aspect), mean, na.rm = TRUE)
write.csv(topography, "tables/topography.csv")

### allometry
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(anova(height_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(canopy_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(diam_lm)[, c(1, 4, 5)]),
                as.matrix(anova(density_lm)[, c(1, 4, 5)])),
          'tables/allometry.csv')

### foliar organics
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(anova(C_foliar_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(N_foliar_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(CN_foliar_lm)[, c(1, 4, 5)])),
          'tables/foliar_cn.csv')

### foliar inorganics
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(anova(Ca_foliar_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(P_foliar_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(K_foliar_lm)[, c(1, 4, 5)]),
                as.matrix(anova(Mg_foliar_lm)[, c(1, 4, 5)]),
                as.matrix(anova(Al_foliar_lm)[, c(1, 4, 5)]),
                as.matrix(anova(Zn_foliar_lm)[, c(1, 4, 5)])),
          'tables/foliar_inorganics.csv')

### foliar isotopes
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(anova(d13C_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(d15N_lm)[, c(4, 5)])),
          'tables/foliar_isotope.csv')

### soil organics
#### create table with degrees of f reedom, f-value, p-value results from linear models
write.csv(cbind(as.matrix(anova(C_soil_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(N_soil_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(CN_soil_lm)[, c(1, 4, 5)])),
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
write.csv(cbind(as.matrix(anova(Ca_soil_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(P_soil_lm)[, c(4, 5)]), 
                as.matrix(anova(K_soil_lm)[, c(4, 5)]),
                as.matrix(anova(Mg_soil_lm)[, c(4, 5)]),
                as.matrix(anova(Al_soil_lm)[, c(4, 5)]),
                as.matrix(anova(Zn_soil_lm)[, c(4, 5)])),
          'tables/soil_inorganics.csv')

(summary(emmeans(K_soil_lm, ~fire))[1,2] - summary(emmeans(K_soil_lm, ~fire))[2,2])/ summary(emmeans(K_soil_lm, ~fire))[2,2]
(summary(emmeans(Ca_soil_lm, ~elevation_fac))[1,2] - summary(emmeans(Ca_soil_lm, ~elevation_fac))[2,2])/ summary(emmeans(Ca_soil_lm, ~elevation_fac))[2,2]

## soil characteristics
write.csv(cbind(as.matrix(anova(retention_lm)[, c(1, 4, 5)]), 
                as.matrix(anova(pH_lm)[, c(4, 5)]), 
                as.matrix(anova(CEC_lm)[, c(4, 5)])),
          'tables/soil_char.csv')

#### save graphs ####

## allometry
ggsave("plots/plots_allometry.jpeg", plot = plots_allometry,
       width = 28, height = 18, units = "cm", dpi = 600) # 3 panels

## foliar isotopes
ggsave("plots/plots_isotopes.jpeg", plot = plots_isotopes,
       width = 42, height = 25, units = "cm", dpi = 600) # 2 panels

## foliar nutrients
ggsave("plots/plots_foliar.jpeg", plot = plots_foliar,
       width = 45, height = 25, units = "cm", dpi = 600) # 5 panels

## soil nutrients
ggsave("plots/plots_soil.jpeg", plot = plots_soil_elev,
       width = 42, height = 25, units = "cm", dpi = 600) # 6 panels

## soil retention
ggsave("plots/plots_retention.jpeg", plot = plot_retention,
       width = 42, height = 25, units = "cm", dpi = 600) # 1 panel



