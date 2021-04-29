# script to clean the mdi pitch pine data (and combine as possible)

library(dplyr)

## read in data
foliar = read.csv('raw/pitchpine_mdi_foliar_final.csv')[, 1:8]
colnames(foliar) = c('ID', 'Name', 'Ca_foliar', 'P_foliar', 
                     'K_foliar', 'Mg_foliar', 'Al_foliar', 'Zn_foliar')
foliar$Ca_foliar <- foliar$Ca_foliar/1000
foliar$P_foliar <- foliar$P_foliar/1000
foliar$K_foliar <- foliar$K_foliar/1000
foliar$Mg_foliar <- foliar$Mg_foliar/1000
foliar$Al_foliar <- foliar$Al_foliar/1000
foliar$Zn_foliar <- foliar$Zn_foliar/1000
biophysical = read.csv('raw/pitchpine_mdi_biophysical.csv', skip = 1)[, c(1:2, 4:6)]
colnames(biophysical) = c('ID', 'Name', 'Height', 'Canopy', 'Diam')
biophysical$Height <- biophysical$Height/100
biophysical$Canopy <- biophysical$Canopy/100
isotope = read.csv('raw/pitchpine_mdi_isotope.csv')[, 1:6]
colnames(isotope) = c('ID', 'Name', 'd13C', 'd15N', 'C_foliar', 'N_foliar')
isotope[isotope$N_foliar > 15,]$N_foliar = NA
retention = read.csv('raw/pitchpine_mdi_moistretent.csv', skip = 1)[, c(1, 2, 7)]
colnames(retention) = c('ID', 'Name', 'Retention')
retention$Retention = as.numeric(gsub("[\\%,]", "", retention$Retention))
soil = read.csv('raw/pitchpine_mdi_soil.csv')[, 1:12]
colnames(soil) = c('ID', 'Name', 'Ca_soil', 'P_soil', 'K_soil', 'Mg_soil', 'Al_soil', 'Zn_soil',
                   'pH', 'CEC', 'C_soil', 'N_soil')
soil$Ca_soil <- soil$Ca_soil/1000
soil$P_soil <- soil$P_soil/1000
soil$K_soil <- soil$K_soil/1000
soil$Mg_soil <- soil$Mg_soil/1000
soil$Al_soil <- soil$Al_soil/1000
soil$Zn_soil <- soil$Zn_soil/1000
geo = read.csv('raw/MDITreeGeoRasterData_2021.csv')
colnames(geo) = c('ID1', 'Name', 'Longitude', 'Latitude', 'Label', 'Elevation', 'Slope', 'Aspect')
geo$ID = NA
geo$ID[geo$Name == 'CAD'] = paste(geo$ID1[geo$Name == 'CAD'], '-', 'HIGHELEV-DIST', 
                                  sep = '')
geo$ID[geo$Name == 'STSAUV'] = paste(geo$ID1[geo$Name == 'STSAUV'], '-', 'HIGHELEV-DIST', 
                                  sep = '')
geo$ID[geo$Name == 'CADCLIFFS'] = paste(geo$ID1[geo$Name == 'CADCLIFFS'], '-', 'LOWELEV-DIST', 
                                     sep = '')
geo$ID[geo$Name == 'WOND'] = paste(geo$ID1[geo$Name == 'WOND'], '-', 'LOWELEV-DIST', 
                                        sep = '')
geo$ID1 <- NULL
geo$Label <- NULL

## combine data
geo_biophysical <- full_join(geo, biophysical, by = c("ID", "Name"))
geo_biophysical_isotope <- full_join(geo_biophysical, isotope, by = c("ID", "Name"))
geo_biophysical_isotope_foliar <- full_join(geo_biophysical_isotope, foliar, by = c("ID", "Name"))
geo_biophysical_isotope_foliar_soil <- full_join(geo_biophysical_isotope_foliar, soil, by = c("ID", "Name"))
mdi_all = full_join(geo_biophysical_isotope_foliar_soil, retention, by = c('ID', 'Name'))
mdi_all <- mdi_all %>% select(ID, everything())

write.csv(mdi_all, 'mdi_all_clean.csv', row.names = F)



