# script to clean the mdi pitch pine data (and combine as possible)

library(dplyr)

## read in data
foliar = read.csv('raw/pitchpine_mdi_foliar_final.csv')[, 1:8]
colnames(foliar) = c('ID', 'Name', 'Ca_foliar', 'P_foliar', 
                     'K_foliar', 'Mg_foliar', 'Al_foliar', 'Zn_foliar')
biophysical = read.csv('raw/pitchpine_mdi_biophysical.csv', skip = 1)[, c(1:2, 4:6)]
colnames(biophysical) = c('ID', 'Name', 'Height', 'Canopy', 'Diam')
isotope = read.csv('raw/pitchpine_mdi_isotope.csv')[, 1:6]
colnames(isotope) = c('ID', 'Name', 'd13C', 'd15N', 'C_foliar', 'N_foliar')
isotope[isotope$N_foliar > 15,]$N_foliar = NA
retention = read.csv('raw/pitchpine_mdi_moistretent.csv', skip = 1)[, c(1, 2, 7)]
colnames(retention) = c('ID', 'Name', 'Retention')
retention$Retention = as.numeric(gsub("[\\%,]", "", retention$Retention))
soil = read.csv('raw/pitchpine_mdi_soil.csv')[, 1:12]
colnames(soil) = c('ID', 'Name', 'Ca_soil', 'P_soil', 'K_soil', 'Mg_soil', 'Al_soil', 'Zn_soil',
                   'pH', 'CEC', 'C_soil', 'N_soil')
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



