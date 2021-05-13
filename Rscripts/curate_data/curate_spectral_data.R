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
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files

### Combine datasets - here we should standardize the columns/names
input_dir <- file.path("~/Data/Dropbox/MANUSCRIPTS/BNL_TEST/canopyRT/data/compiled_data")
panama_2016_leaf_refl <- read.csv(file = file.path(input_dir,"NGEE-Tropics_Panama_2016_Leaf_Spectral_Reflectance.csv"))
panama_2017_leaf_refl <- read.csv(file = file.path(input_dir,"NGEE-Tropics_Panama_2017_Leaf_Spectral_Reflectance.csv"))


NGEETropics_leaf_reflectance <- list(Panama2016_leaf_refl = panama_2016_leaf_refl,
                                    Panama2017_leaf_refl = panama_2017_leaf_refl)
  
save(NGEETropics_leaf_reflectance, file = file.path(input_dir,"NGEETropics_Leaf_Reflectance.RData"))
#--------------------------------------------------------------------------------------------------#



#--------------------------------------------------------------------------------------------------#
### EOF
