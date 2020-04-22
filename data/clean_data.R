# script to clean the mdi pitch pine data (and combine as possible)

library(dplyr)

## read in data
foliar = read.csv('raw/pitchpine_mdi_foliar_final.csv')[1:45, 1:8]
colnames(foliar) = c('ID', 'Name', 'Ca_foliar', 'P_foliar', 
                     'K_foliar', 'Mg_foliar', 'Al_foliar', 'Zn_foliar')
biophysical = read.csv('raw/pitchpine_mdi_biophysical.csv', skip = 1)[1:39, c(1:2, 4:6)]
isotope = read.csv('raw/pitchpine_mdi_isotope.csv')[1:41, 1:6]
colnames(isotope) = c('ID', 'Name', 'd13C', 'd15N', 'C_foliar', 'N_foliar')
retention = read.csv('raw/pitchpine_mdi_moistretent.csv', skip = 1)[1:41, c(1, 2, 7)]
colnames(retention) = c('ID', 'Name', 'retention')
retention$retention = as.numeric(gsub("[\\%,]", "", retention$retention))
soil = read.csv('raw/pitchpine_mdi_soil.csv')[1:32, 1:12]
colnames(soil) = c('ID', 'Name', 'Ca_soil', 'P_soil', 'K_soil', 'Mg_soil', 'Al_soil', 'Zn_soil',
                   'pH', 'CEC', 'C_soil', 'N_soil')

## combine data
biophysical_isotope = full_join(biophysical, isotope, by = c('ID', 'Name'))
biophysical_isotope_foliar = full_join(biophysical_isotope, foliar, by = c('ID', 'Name'))
biophysical_isotope_foliar_soil = full_join(biophysical_isotope_foliar, soil, by = c('ID', 'Name'))
mdi_all = full_join(biophysical_isotope_foliar_soil, retention, by = c('ID', 'Name'))

mdi_all_sorted = arrange(mdi_all, Name, ID)

write.csv(mdi_all_sorted, 'mdi_all_clean.csv')



