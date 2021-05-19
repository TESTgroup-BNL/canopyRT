####################################################################################################
#
#    --- Last updated: 05.12.2021 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Load libraries
devtools::install_github(repo = "TESTgroup-BNL/spectratrait", dependencies=TRUE)

list.of.packages <- c("devtools","readr","RCurl","httr","dplyr","spectratrait")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=c("Depends", "Imports",
                                                                       "LinkingTo"))
invisible(lapply(list.of.packages, library, character.only = TRUE))


# not in
`%notin%` <- Negate(`%in%`)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### output location
output_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/canopyRT/")
if (! file.exists(output_dir)) dir.create(output_dir,recursive=TRUE)
  outdir <- file.path(path.expand(output_dir))
setwd(outdir) # set working directory
getwd()  # check wd
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### PLSR Coefficients
plsr_coeff_dir <- file.path(output_dir,"plsr_coefficients")
leaf_age_plsr_coeff <- read.csv(file = file.path(plsr_coeff_dir, 
                                                 "Wu_etal_spectra_leaf_age_plsr_coefficients.csv"),
                                header = T)
head(leaf_age_plsr_coeff)[1:10]

mean_LeafAge.plsr.intercept <- mean(leaf_age_plsr_coeff[,2])
mean_LeafAge.plsr.coeffs.vec <- colMeans(leaf_age_plsr_coeff[,which(names(leaf_age_plsr_coeff) %in% 
                                                                      paste0("Wave_",seq(400,2400,1)))])
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
###### Panama datasets
# 2016
ecosis_id <- "22dc6b53-5d4a-45c6-9d02-000d0f0ec5a0"
dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:45]

# clean up
dat_raw[dat_raw==-9999]=NA
head(dat_raw)

panama_2016_leaf_refl <- dat_raw %>%
  select(Site=Location,Sample_ID=BNL_Barcode,Sample_Date=Sample_Date,Species_Code=`Species_Code`,
         Genus=`Latin Genus`,Species=`Latin Species`,Canopy_position,Instrument=`Instrument Model`,`350`:`2500`)
head(panama_2016_leaf_refl)[,1:10]
spectra <- panama_2016_leaf_refl %>% 
  select(`350`:`2500`) %>%
  setNames(paste0('Wave_', names(.)))
temp_sub_spec <- droplevels(spectra[,which(names(spectra) %in% 
                                          paste0("Wave_",seq(400,2400,1)))])*0.01
plsr_leaf_age_days <- ((as.matrix(temp_sub_spec) %*% mean_LeafAge.plsr.coeffs.vec)+mean_LeafAge.plsr.intercept)^2
plsr_leaf_age_days[plsr_leaf_age_days<0] <- NA
hist(plsr_leaf_age_days)
spec_info <- panama_2016_leaf_refl[,names(panama_2016_leaf_refl) %notin% seq(350,2500,1)]
spec_info$plsr_leaf_age_days <- plsr_leaf_age_days
spec_info <- spec_info %>%
  select(Site, Sample_ID, Sample_Date, Species_Code, Genus, Species, Canopy_position, 
         PLSR_LeafAge_days=plsr_leaf_age_days, Instrument)
panama_2016_leaf_refl_out <- data.frame(spec_info,spectra)
write.csv(panama_2016_leaf_refl_out, 
          file = file.path(output_dir,"data","compiled_data",
                           "NGEE-Tropics_Panama_2016_Leaf_Spectral_Reflectance.csv"), row.names = F)
rm(ecosis_id,dat_raw,panama_2016_leaf_refl,panama_2016_leaf_refl_out,spectra,spec_info,temp_sub_spec)


# 2017
ecosis_id <- "c4ce128e-7984-4325-8e8d-e8d54b3f783e"
dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:55]

# clean up
dat_raw[dat_raw==-9999]=NA
head(dat_raw)

panama_2017_leaf_refl <- dat_raw %>%
  select(Site=`Location Name`,Sample_ID=SampleID,Sample_Date=`Measurement Date`,Species_Code=Species_Code,
         Genus=`Latin Genus`,Species=`Latin Species`,Canopy_position,Branch_Number,GasEx_Leaf_Age,
         Relative_leaf_age,GasEx_Leaf_Age_Revised,Leaf_age_revised,
         Sample_type,Instrument=`Instrument Model`,`350`:`2500`)
head(panama_2017_leaf_refl)[,1:35]
spectra <- panama_2017_leaf_refl %>% 
  select(`350`:`2500`) %>%
  setNames(paste0('Wave_', names(.)))
temp_sub_spec <- droplevels(spectra[,which(names(spectra) %in% 
                                             paste0("Wave_",seq(400,2400,1)))])*0.01
plsr_leaf_age_days <- ((as.matrix(temp_sub_spec) %*% mean_LeafAge.plsr.coeffs.vec)+mean_LeafAge.plsr.intercept)^2
plsr_leaf_age_days[plsr_leaf_age_days<0] <- NA
hist(plsr_leaf_age_days)
spec_info <- panama_2017_leaf_refl[,names(panama_2017_leaf_refl) %notin% seq(350,2500,1)]
spec_info$plsr_leaf_age_days <- plsr_leaf_age_days
spec_info <- spec_info %>%
  select(Site, Sample_ID, Sample_Date, Species_Code, Genus, Species, Canopy_position,
         Branch_Number,GasEx_Leaf_Age,Relative_leaf_age,GasEx_Leaf_Age_Revised,Leaf_age_revised,
         Sample_type,PLSR_LeafAge_days=plsr_leaf_age_days, Instrument)
panama_2017_leaf_refl_out <- data.frame(spec_info,spectra)
write.csv(panama_2017_leaf_refl_out, 
          file = file.path(output_dir,"data","compiled_data",
                           "NGEE-Tropics_Panama_2017_Leaf_Spectral_Reflectance.csv"), row.names = F)
rm(ecosis_id,dat_raw,panama_2017_leaf_refl,panama_2017_leaf_refl_out,spectra,spec_info,temp_sub_spec)

###### Puerto Rico datasets
# 2017 - reflectance
ecosis_id <- "7b5fd8dc-fa19-4fc0-8435-58acff826e9d"
dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:55]

# clean up
dat_raw[dat_raw==-9999]=NA
head(dat_raw)

pr_2017_leaf_refl <- dat_raw %>% 
  filter(Measurement=="Reflectance") %>%
  select(Site=`Location Name`,Sample_ID=SampleID,Sample_Date=`Measurement Date`,Species_Code=Species_ID,
         Genus=`Latin Genus`,Species=`Latin Species`,Canopy_position=Canopy_Position,Light_Exposure,
         Relative_Leaf_Age=Leaf_Age,Sample_Type,Latitude,Longitude,Instrument=`Instrument Model`,
         `350`:`2500`)
head(pr_2017_leaf_refl)[,1:35]
spectra <- pr_2017_leaf_refl %>% 
  select(`350`:`2500`) %>%
  setNames(paste0('Wave_', names(.)))
temp_sub_spec <- droplevels(spectra[,which(names(spectra) %in% 
                                             paste0("Wave_",seq(400,2400,1)))])*0.01
plsr_leaf_age_days <- ((as.matrix(temp_sub_spec) %*% mean_LeafAge.plsr.coeffs.vec)+mean_LeafAge.plsr.intercept)^2
plsr_leaf_age_days[plsr_leaf_age_days<0] <- NA
hist(plsr_leaf_age_days)    
spec_info <- pr_2017_leaf_refl[,names(pr_2017_leaf_refl) %notin% seq(350,2500,1)]
spec_info$plsr_leaf_age_days <- plsr_leaf_age_days
spec_info <- spec_info %>%
  select(Site, Sample_ID, Sample_Date, Species_Code, Genus, Species, Canopy_position,
         Light_Exposure,Relative_Leaf_Age,Sample_Type,PLSR_LeafAge_days=plsr_leaf_age_days,
         Latitude,Longitude,Instrument)
pr_2017_leaf_refl_out <- data.frame(spec_info,spectra)
write.csv(pr_2017_leaf_refl_out, 
          file = file.path(output_dir,"data","compiled_data",
                           "NGEE-Tropics_PuertoRico_2017_Leaf_Spectral_Reflectance.csv"), row.names = F)
rm(ecosis_id,dat_raw,pr_2017_leaf_refl,pr_2017_leaf_refl_out,spectra,spec_info,temp_sub_spec)


# 2017 - transmittance
ecosis_id <- "7b5fd8dc-fa19-4fc0-8435-58acff826e9d"
dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:55]

# clean up
dat_raw[dat_raw==-9999]=NA
head(dat_raw)

pr_2017_leaf_trans <- dat_raw %>% 
  filter(Measurement=="Transmittance") %>%
  select(Site=`Location Name`,Sample_ID=SampleID,Sample_Date=`Measurement Date`,Species_Code=Species_ID,
         Genus=`Latin Genus`,Species=`Latin Species`,Canopy_position=Canopy_Position,Light_Exposure,
         Relative_Leaf_Age=Leaf_Age,Sample_Type,Latitude,Longitude,Instrument=`Instrument Model`,
         `350`:`2500`)
head(pr_2017_leaf_trans)[,1:35]

spectra <- pr_2017_leaf_trans %>% 
  select(`350`:`2500`) %>%
  setNames(paste0('Wave_', names(.)))

write.csv(pr_2017_leaf_trans, 
          file = file.path(output_dir,"data","compiled_data",
                           "NGEE-Tropics_PuertoRico_2017_Leaf_Spectral_Transmittance.csv"), row.names = F)
rm(ecosis_id,dat_raw,pr_2017_leaf_trans)


# 2017 - trait data
ecosis_id <- "7b5fd8dc-fa19-4fc0-8435-58acff826e9d"
dat_raw <- get_ecosis_data(ecosis_id = ecosis_id)
head(dat_raw)
names(dat_raw)[1:55]

# clean up
dat_raw[dat_raw==-9999]=NA
head(dat_raw)

pr_2017_leaf_traits <- dat_raw %>% 
  select(Site=`Location Name`,Sample_ID=SampleID,Sample_Date=`Measurement Date`,Species_Code=Species_ID,
         Genus=`Latin Genus`,Species=`Latin Species`,Canopy_position=Canopy_Position,Light_Exposure,
         Sample_Type,Latitude,Longitude,Cmass_mg_g=Cmass,Nmass_mg_g=Nmass,CNratio=CNratio,
         Carea_g_m2=Carea,Narea_g_m2=Narea,LMA_gDW_m2=LMA)
head(pr_2017_leaf_traits)

write.csv(pr_2017_leaf_traits, 
          file = file.path(output_dir,"data","compiled_data",
                           "NGEE-Tropics_PuertoRico_2017_Leaf_Traits.csv"), row.names = F)
rm(ecosis_id,dat_raw,pr_2017_leaf_traits)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

### Combine datasets - here we should standardize the columns/names
input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/canopyRT/data/compiled_data")
panama_2016_leaf_refl <- read.csv(file = file.path(input_dir,"NGEE-Tropics_Panama_2016_Leaf_Spectral_Reflectance.csv"))
panama_2017_leaf_refl <- read.csv(file = file.path(input_dir,"NGEE-Tropics_Panama_2017_Leaf_Spectral_Reflectance.csv"))
pr_2017_leaf_refl <- read.csv(file = file.path(input_dir,"NGEE-Tropics_PuertoRico_2017_Leaf_Spectral_Reflectance.csv"))
pr_2017_leaf_trans <- read.csv(file = file.path(input_dir,"NGEE-Tropics_PuertoRico_2017_Leaf_Spectral_Transmittance.csv"))

NGEETropics_leaf_reflectance <- list(Panama2016_leaf_refl = panama_2016_leaf_refl,
                                    Panama2017_leaf_refl = panama_2017_leaf_refl,
                                    PuertoRico_2017_leaf_refl = pr_2017_leaf_refl)
  
save(NGEETropics_leaf_reflectance, file = file.path(input_dir,"NGEETropics_Leaf_Reflectance.RData"))

NGEETropics_leaf_transmittance <- list(PuertoRico_2017_leaf_trans = pr_2017_leaf_trans)

save(NGEETropics_leaf_transmittance, file = file.path(input_dir,"NGEETropics_Leaf_Transmittance.RData"))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF