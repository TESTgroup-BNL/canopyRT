
####################################################################################################
#
#  Apply Wu et al 2017 to 2020 Panama samples
#
####################################################################################################



#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Load libraries
list.of.packages <- c("pls","dplyr","reshape2","here","plotrix","ggplot2","gridExtra",
                      "spectratrait")
invisible(lapply(list.of.packages, library, character.only = TRUE))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup other functions and options
# not in
`%notin%` <- Negate(`%in%`)

# Default par options
opar <- par(no.readonly = T)

# Specify output directory, output_dir 
# Options: 
# tempdir - use a OS-specified temporary directory 
# user defined PATH - e.g. "~/scratch/PLSR"
output_dir <- file.path("~/Data/ngee_tropics/")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Set working directory
if (output_dir=="tempdir") {
  outdir <- tempdir()
} else {
  if (! file.exists(output_dir)) dir.create(output_dir,recursive=TRUE)
  outdir <- file.path(path.expand(output_dir))
}
setwd(outdir) # set working directory
getwd()  # check wd
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### PLSR Coefficients
plsr_coeff_dir <- file.path("/Volumes/Test/Projects/NGEE-Tropics/Data/Panama/2020/R_Scripts/PLSR/plsr_coefficients/")
leaf_age_plsr_coeff <- read.csv(file = file.path(plsr_coeff_dir, 
                                                 "Wu_etal_spectra_leaf_age_plsr_coefficients.csv"),
                                header = T)
head(leaf_age_plsr_coeff)[1:10]
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
spectra_dir <- file.path("~/Data/ngee_tropics/")
load(file.path(spectra_dir,"Leaf_Spec_PA_2020.Rdata"))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Prepare new data for estimation
Start.wave <- 400
End.wave <- 2400
wv <- seq(Start.wave,End.wave,1)
Spectra <- as.matrix(Leaf_Spec_PA_2020[,names(Leaf_Spec_PA_2020) %in% paste0("Wave_",wv)])
colnames(Spectra) <- c(paste0("Wave_",wv))
head(Spectra)[1:6,1:10]
sample_info <- Leaf_Spec_PA_2020[,names(Leaf_Spec_PA_2020) %notin% paste0("Wave_",seq(350,2500,1))]
head(sample_info)

sample_info2 <- sample_info %>%
  select(Spectra_Name=Spectra,SampleID, Long.identifier, Leaf_Age=Age, Species)
head(sample_info2)

plsr_data <- data.frame(sample_info2,Spectra)
rm(sample_info,sample_info2,Spectra)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
print("**** Applying PLSR model to estimate LMA from spectral observations ****")
# setup model
dims <- dim(leaf_age_plsr_coeff)
mean_LeafAge.plsr.intercept <- mean(leaf_age_plsr_coeff[,2])
mean_LeafAge.plsr.coeffs.vec <- colMeans(leaf_age_plsr_coeff[,which(names(leaf_age_plsr_coeff) %in% 
                                                                      paste0("Wave_",seq(Start.wave,End.wave,1)))])

sub_spec <- droplevels(plsr_data[,which(names(plsr_data) %in% 
                                          paste0("Wave_",seq(Start.wave,End.wave,1)))])*0.01

temp <- as.matrix(sub_spec) %*% mean_LeafAge.plsr.coeffs.vec
LeafAge <- data.frame(rowSums(temp))+mean_LeafAge.plsr.intercept
LeafAge <- (LeafAge[,1])^2
LeafAge[LeafAge<0] <- NA
hist(LeafAge)

# organize output
sample_info <- plsr_data[,which(names(plsr_data) %notin% 
                                  paste0("Wave_",seq(Start.wave,End.wave,1)))]
LeafAge.PLSR.dataset <- data.frame(sample_info, FS_PLSR_LeafAge_days=LeafAge)
head(LeafAge.PLSR.dataset)


# older code to do this
LeafAge.plsr.jk.coeffs <- leaf_age_plsr_coeff[,which(names(leaf_age_plsr_coeff) %in% 
                                                       paste0("Wave_",seq(Start.wave,End.wave,1)))]
head(LeafAge.plsr.jk.coeffs)[1:10]
intercepts <- leaf_age_plsr_coeff[,2]
print("**** Deriving uncertainty estimates ****")
dims <- dim(LeafAge.plsr.jk.coeffs)
jk.leaf.age.est <- array(data=NA,dim=c(dim(sub_spec)[1],dims[1]))
for (i in 1:length(intercepts)){
  coefs <- as.vector(unlist(LeafAge.plsr.jk.coeffs[i,]))
  temp <- as.matrix(sub_spec) %*% coefs
  values <- data.frame(rowSums(temp))+intercepts[i]
  jk.leaf.age.est[,i] <- (values[,1]^2)
  rm(temp, values, coefs)
}

jk.leaf.age.est.quant <- apply(jk.leaf.age.est,1,quantile,probs=c(0.025,0.975))
jk.leaf.age.est.quant2 <- data.frame(t(jk.leaf.age.est.quant))
names(jk.leaf.age.est.quant2) <- c("FS_PLSR_LeafAge_L5","FS_PLSR_LeafAge_U95")
jk.leaf.age.est.sd <- apply(jk.leaf.age.est,1,sd)

## Combine into final dataset
stats <- data.frame(jk.leaf.age.est.sd,jk.leaf.age.est.quant2)
names(stats) <- c("FS_PLSR_LeafAge_days_Sdev","FS_PLSR_LeafAge_L5","FS_PLSR_LeafAge_U95")
LeafAge.PLSR.dataset.out <- data.frame(LeafAge.PLSR.dataset,stats)
head(LeafAge.PLSR.dataset.out)


# output results
write.csv(x = LeafAge.PLSR.dataset.out, file = file.path(output_dir,"PLSR_estimated_LeafAge_data.csv"),
          row.names = F)
#--------------------------------------------------------------------------------------------------#

