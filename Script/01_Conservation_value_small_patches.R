# sessinInfo()
# R version 4.1.0 (2021-05-18)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252    LC_MONETARY=English_Canada.1252 LC_NUMERIC=C                   
# [5] LC_TIME=English_Canada.1252    
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] RColorBrewer_1.1-2 AICcmodavg_2.3-1   glmmTMB_1.0.2.1    ggpubr_0.4.0       ggthemes_4.2.4     ggplot2_3.3.5      vegan_2.5-7       
# [8] lattice_0.20-44    permute_0.9-5      tidyr_1.1.3        data.table_1.14.2 
# 
# loaded via a namespace (and not attached):
#   [1] VGAM_1.1-5       splines_4.1.0    carData_3.0-4    assertthat_0.2.1 sp_1.4-5         stats4_4.1.0     cellranger_1.1.0
# [8] pillar_1.6.4     backports_1.2.1  glue_1.4.2       ggsignif_0.6.1   minqa_1.2.4      colorspace_2.0-2 sandwich_3.0-1  
# [15] plyr_1.8.6       Matrix_1.3-3     pkgconfig_2.0.3  raster_3.4-10    broom_0.7.6      haven_2.4.1      purrr_0.3.4     
# [22] xtable_1.8-4     mvtnorm_1.1-2    scales_1.1.1     openxlsx_4.2.3   rio_0.5.26       lme4_1.1-27      emmeans_1.7.0   
# [29] tibble_3.1.2     mgcv_1.8-35      generics_0.1.1   car_3.0-10       ellipsis_0.3.2   TH.data_1.0-10   withr_2.4.2     
# [36] TMB_1.7.20       survival_3.2-11  magrittr_2.0.1   crayon_1.4.2     readxl_1.3.1     estimability_1.3 fansi_0.5.0     
# [43] nlme_3.1-152     MASS_7.3-54      rstatix_0.7.0    forcats_0.5.1    foreign_0.8-81   tools_4.1.0      hms_1.1.1       
# [50] lifecycle_1.0.1  multcomp_1.4-17  stringr_1.4.0    munsell_0.5.0    cluster_2.1.2    zip_2.1.1        compiler_4.1.0  
# [57] tinytex_0.35     rlang_0.4.11     grid_4.1.0       nloptr_1.2.2.2   boot_1.3-28      gtable_0.3.0     codetools_0.2-18
# [64] abind_1.4-5      DBI_1.1.1        curl_4.3.1       R6_2.5.1         zoo_1.8-9        dplyr_1.0.7      utf8_1.2.2      
# [71] stringi_1.7.6    parallel_4.1.0   unmarked_1.1.1   Rcpp_1.0.7       vctrs_0.3.8      tidyselect_1.1.1 xfun_0.23       
# [78] coda_0.19-4     


#### packages used
library(data.table)
library(tidyr)
library(vegan)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(glmmTMB)
library(AICcmodavg)
library(RColorBrewer)

################################################
#### analysis preparation
################################################

#### open functions used in the analysis
source("Script//00_functions.R")

#### import data
## open FragSAD dataset, available at https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2861
data_fragsad = fread("Data\\FragSAD_data.csv", header = TRUE)
metadata_fragsad = fread("Data\\FragSAD_metadata.csv", header = TRUE)

#### open IUCN declining species; Search on 2021-06-09 at 13:30:36 at https://www.iucnredlist.org/
data_iucn = fread("Data\\Red_list_declining.csv", header = TRUE)

## remove studies with equal effort in patches of different sizes, and patches for which patch size was not available
metadata_fragsad2 <- subset(metadata_fragsad, metadata_fragsad$site_size != 'continuous')
metadata_fragsad2 <- subset(metadata_fragsad2, metadata_fragsad2$sampling_effort != 1)
metadata_fragsad2$refshort <- as.factor(metadata_fragsad2$refshort)
metadata_fragsad2$site_size <- as.numeric(metadata_fragsad2$site_size)

# to repeat the analysis assessing declining species, run the following two lines of code (removes species that are not declining based on IUCN Red List)
#data_fragsad_declining <- data_fragsad[data_fragsad$scientific_name %in% data_iucn$scientificName,]
#data_fragsad <- data_fragsad_declining

################################################
#### creation of table for analysis
################################################

#  Baldi_1999 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Baldi_1999") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Baldi_1999") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

# evaluate which patches should be removed
# hist(data_sampling$site_size)
# mean(data_sampling$site_size) + 3* sd(data_sampling$site_size)
# mean(data_sampling$site_size) - 3* sd(data_sampling$site_size)

data_sampling <- subset(data_sampling, data_sampling$site_size < 26)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## evenness
H <- diversity(data_sampling$site_size) ## Shannon's diversity of patches
J <- H/log(nrow(data_sampling)) ## divide Shannon's diversity by the natural logarithm of the number of patches

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 26)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Baldi_1999 <-   data.frame("study" = c( data_meta[1,2 ]),
                           "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                           "gamma" = max(result$SL$species),
                           "number of patches" =  nrow(sloss_data), 
                           "minimum patch area" =   min(sloss_data$site_size), 
                           "maximum patch area" =   max(sloss_data$site_size), 
                           "sd patch area" = sd(sloss_data$site_size),
                           "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                           "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)

# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Bell_2006_a ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Bell_2006_a") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Bell_2006_a") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
# repeat after removing the 1000 ha patch
data_sampling <- subset(data_sampling, data_sampling$site_size < 1000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 1000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Bell_2006_a <-   data.frame("study" = c( data_meta[1,2 ]),
                            "comparison" = c("inconclusive"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                            "gamma" = max(result$SL$species),
                            "number of patches" =  nrow(sloss_data), 
                            "minimum patch area" =   min(sloss_data$site_size), 
                            "maximum patch area" =   max(sloss_data$site_size), 
                            "sd patch area" = sd(sloss_data$site_size),
                            "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                            "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
                            
)


# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Bell_2006_b ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Bell_2006_b") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Bell_2006_b") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 1000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 1000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


Bell_2006_b <-   data.frame("study" = c( data_meta[1,2 ]),
                            "comparison" = c("inconclusive"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                            "gamma" = max(result$SL$species),
                            "number of patches" =  nrow(sloss_data), 
                            "minimum patch area" =   min(sloss_data$site_size), 
                            "maximum patch area" =   max(sloss_data$site_size), 
                            "sd patch area" = sd(sloss_data$site_size),
                            "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                            "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)


# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Berg_1997 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Berg_1997") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Berg_1997") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 999)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 999)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Berg_1997 <-   data.frame("study" = c( data_meta[1,2 ]),
                          "comparison" = c("inconclusive"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                          "gamma" = max(result$SL$species),
                          "number of patches" =  nrow(sloss_data), 
                          "minimum patch area" =   min(sloss_data$site_size), 
                          "maximum patch area" =   max(sloss_data$site_size), 
                          "sd patch area" = sd(sloss_data$site_size),
                          "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                          "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
                          
)


# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Bernard_2007 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Bernard_2007") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Bernard_2007") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 300)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 300)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


Bernard_2007 <-   data.frame("study" = c( data_meta[1,2 ]),
                             "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                             "gamma" = max(result$SL$species),
                             "number of patches" =  nrow(sloss_data), 
                             "minimum patch area" =   min(sloss_data$site_size), 
                             "maximum patch area" =   max(sloss_data$site_size), 
                             "sd patch area" = sd(sloss_data$site_size),
                             "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                             "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
                             
)



# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Bolger_1997 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Bolger_1997") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Bolger_1997") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 84)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

H <- diversity(data_sampling$site_size)
J <- H/log(nrow(data_sampling))
## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 84)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


Bolger_1997 <-   data.frame("study" = c( data_meta[1,2 ]),
                            "comparison" = c("inconclusive"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                            "gamma" = max(result$SL$species),
                            "number of patches" =  nrow(sloss_data), 
                            "minimum patch area" =   min(sloss_data$site_size), 
                            "maximum patch area" =   max(sloss_data$site_size), 
                            "sd patch area" = sd(sloss_data$site_size),
                            "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                            "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Bossart_2006 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Bossart_2006") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Bossart_2006") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 6000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 6000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


Bossart_2006 <-   data.frame("study" = c( data_meta[1,2 ]),
                             "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                             "gamma" = max(result$SL$species),
                             "number of patches" =  nrow(sloss_data), 
                             "minimum patch area" =   min(sloss_data$site_size), 
                             "maximum patch area" =   max(sloss_data$site_size), 
                             "sd patch area" = sd(sloss_data$site_size),
                             "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                             "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)


# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Bossart_2016 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Bossart_2016") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Bossart_2016") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 6000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 6000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Bossart_2016 <-   data.frame("study" = c( data_meta[1,2 ]),
                             "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                             "gamma" = max(result$SL$species),
                             "number of patches" =  nrow(sloss_data), 
                             "minimum patch area" =   min(sloss_data$site_size), 
                             "maximum patch area" =   max(sloss_data$site_size), 
                             "sd patch area" = sd(sloss_data$site_size),
                             "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                             "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Chiarello_1999 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Chiarello_1999") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Chiarello_1999") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Chiarello_1999 <-   data.frame("study" = c( data_meta[1,2 ]),
                               "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                               "gamma" = max(result$SL$species),
                               "number of patches" =  nrow(sloss_data), 
                               "minimum patch area" =   min(sloss_data$site_size), 
                               "maximum patch area" =   max(sloss_data$site_size), 
                               "sd patch area" = sd(sloss_data$site_size),
                               "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                               "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)

# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Dami_2012 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Dami_2012") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Dami_2012") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 34)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 34)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Dami_2012 <-   data.frame("study" = c( data_meta[1,2 ]),
                          "comparison" = c("inconclusive"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                          "gamma" = max(result$SL$species),
                          "number of patches" =  nrow(sloss_data), 
                          "minimum patch area" =   min(sloss_data$site_size), 
                          "maximum patch area" =   max(sloss_data$site_size), 
                          "sd patch area" = sd(sloss_data$site_size),
                          "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                          "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)

# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Dauber_2006 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Dauber_2006") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Dauber_2006") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(13,31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling

# must aggregate at the site level
sampling_effort <- aggregate(. ~site_id, data=data_sampling, sum, na.rm=TRUE)
patch_area <- aggregate(. ~site_id, data=data_sampling, min, na.rm=TRUE)
sampling_and_patch_area <- cbind(sampling_effort, patch_area)
sampling_and_patch_area <- sampling_and_patch_area[,c(1, 2,6)]
data_sampling <- sampling_and_patch_area

data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

sloss_data_abundance <- aggregate(. ~site_id, data=sloss_data, sum, na.rm=TRUE)
sloss_data_abundance[,2] <- data_sampling[,3]
sloss_data <- sloss_data_abundance
# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

sloss_data$site_id<- as.factor(sloss_data$site_id)
row.names(sloss_data) <- sloss_data$site_id


Dauber_2006 <-   data.frame("study" = c( data_meta[1,2 ]),
                            "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                            "gamma" = max(result$SL$species),
                            "number of patches" =  nrow(sloss_data), 
                            "minimum patch area" =   min(sloss_data$site_size), 
                            "maximum patch area" =   max(sloss_data$site_size), 
                            "sd patch area" = sd(sloss_data$site_size),
                            "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                            "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
                            
)


# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Ding_2013 ---- 
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Ding_2013") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Ding_2013") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 100)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 100)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)
plot(result, show.legend=TRUE, show.index = FALSE, main = "Empirical SLOSS comparison", ylab = "Species richness", xlab = "Area (ha)")

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


Ding_2013 <-   data.frame("study" = c( data_meta[1,2 ]),
                          "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                          "gamma" = max(result$SL$species),
                          "number of patches" =  nrow(sloss_data), 
                          "minimum patch area" =   min(sloss_data$site_size), 
                          "maximum patch area" =   max(sloss_data$site_size), 
                          "sd patch area" = sd(sloss_data$site_size),
                          "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                          "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)


# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Gavish_2012_a ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Gavish_2012_a") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Gavish_2012_a") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 100)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 100)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Gavish_2012_a <-   data.frame("study" = c( data_meta[1,2 ]),
                              "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                              "gamma" = max(result$SL$species),
                              "number of patches" =  nrow(sloss_data), 
                              "minimum patch area" =   min(sloss_data$site_size), 
                              "maximum patch area" =   max(sloss_data$site_size), 
                              "sd patch area" = sd(sloss_data$site_size),
                              "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                              "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)


## plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Gavish_2012_b ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Gavish_2012_b") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Gavish_2012_b") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 100)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 100)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Gavish_2012_b <-   data.frame("study" = c( data_meta[1,2 ]),
                              "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                              "gamma" = max(result$SL$species),
                              "number of patches" =  nrow(sloss_data), 
                              "minimum patch area" =   min(sloss_data$site_size), 
                              "maximum patch area" =   max(sloss_data$site_size), 
                              "sd patch area" = sd(sloss_data$site_size),
                              "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                              "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)

## plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)


####################################################################################################
####################################################################################################

#  Gavish_2012_c ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Gavish_2012_c") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Gavish_2012_c") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 300)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 300)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Gavish_2012_c <-   data.frame("study" = c( data_meta[1,2 ]),
                              "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                              "gamma" = max(result$SL$species),
                              "number of patches" =  nrow(sloss_data), 
                              "minimum patch area" =   min(sloss_data$site_size), 
                              "maximum patch area" =   max(sloss_data$site_size), 
                              "sd patch area" = sd(sloss_data$site_size),
                              "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                              "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



## plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Guadagnin_2005 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Guadagnin_2005") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Guadagnin_2005") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Guadagnin_2005 <-   data.frame("study" = c( data_meta[1,2 ]),
                               "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                               "gamma" = max(result$SL$species),
                               "number of patches" =  nrow(sloss_data), 
                               "minimum patch area" =   min(sloss_data$site_size), 
                               "maximum patch area" =   max(sloss_data$site_size), 
                               "sd patch area" = sd(sloss_data$site_size),
                               "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                               "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



# # plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Henry_2007 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Henry_2007") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Henry_2007") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


Henry_2007 <-   data.frame("study" = c( data_meta[1,2 ]),
                           "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                           "gamma" = max(result$SL$species),
                           "number of patches" =  nrow(sloss_data), 
                           "minimum patch area" =   min(sloss_data$site_size), 
                           "maximum patch area" =   max(sloss_data$site_size), 
                           "sd patch area" = sd(sloss_data$site_size),
                           "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                           "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



# # plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Jauker_2019 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Jauker_2019") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Jauker_2019") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


Jauker_2019 <-   data.frame("study" = c( data_meta[1,2 ]),
                            "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                            "gamma" = max(result$SL$species),
                            "number of patches" =  nrow(sloss_data), 
                            "minimum patch area" =   min(sloss_data$site_size), 
                            "maximum patch area" =   max(sloss_data$site_size), 
                            "sd patch area" = sd(sloss_data$site_size),
                            "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                            "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)

# # plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Kapoor_2008 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Kapoor_2008") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Kapoor_2008") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Kapoor_2008 <-   data.frame("study" = c( data_meta[1,2 ]),
                            "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                            "gamma" = max(result$SL$species),
                            "number of patches" =  nrow(sloss_data), 
                            "minimum patch area" =   min(sloss_data$site_size), 
                            "maximum patch area" =   max(sloss_data$site_size), 
                            "sd patch area" = sd(sloss_data$site_size),
                            "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                            "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)

# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  MacNally_2001 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "MacNally_2001") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "MacNally_2001") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


MacNally_2001 <-   data.frame("study" = c( data_meta[1,2 ]),
                              "comparison" = c("inconclusive"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                              "gamma" = max(result$SL$species),
                              "number of patches" =  nrow(sloss_data), 
                              "minimum patch area" =   min(sloss_data$site_size), 
                              "maximum patch area" =   max(sloss_data$site_size), 
                              "sd patch area" = sd(sloss_data$site_size),
                              "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                              "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)

# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Martensen_2012 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Martensen_2012") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Martensen_2012") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Martensen_2012 <-   data.frame("study" = c( data_meta[1,2 ]),
                               "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                               "gamma" = max(result$SL$species),
                               "number of patches" =  nrow(sloss_data), 
                               "minimum patch area" =   min(sloss_data$site_size), 
                               "maximum patch area" =   max(sloss_data$site_size), 
                               "sd patch area" = sd(sloss_data$site_size),
                               "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                               "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



## plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)


####################################################################################################
####################################################################################################

#  Manu_2007 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Manu_2007") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Manu_2007") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(13,31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling

# must aggregate at the site level
sampling_effort <- aggregate(. ~site_id, data=data_sampling, sum, na.rm=TRUE)
patch_area <- aggregate(. ~site_id, data=data_sampling, min, na.rm=TRUE)
sampling_and_patch_area <- cbind(sampling_effort, patch_area)
sampling_and_patch_area <- sampling_and_patch_area[,c(1, 2,6)]
data_sampling <- sampling_and_patch_area

data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

sloss_data_abundance <- aggregate(. ~site_id, data=sloss_data, sum, na.rm=TRUE)
sloss_data_abundance[,2] <- data_sampling[,3]
sloss_data <- sloss_data_abundance
# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


Manu_2007 <-   data.frame("study" = c( data_meta[1,2 ]),
                          "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                          "gamma" = max(result$SL$species),
                          "number of patches" =  nrow(sloss_data), 
                          "minimum patch area" =   min(sloss_data$site_size), 
                          "maximum patch area" =   max(sloss_data$site_size), 
                          "sd patch area" = sd(sloss_data$site_size),
                          "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                          "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)

# # plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  McCollin_1993 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "McCollin_1993") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "McCollin_1993") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

McCollin_1993 <-   data.frame("study" = c( data_meta[1,2 ]),
                              "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                              "gamma" = max(result$SL$species),
                              "number of patches" =  nrow(sloss_data), 
                              "minimum patch area" =   min(sloss_data$site_size), 
                              "maximum patch area" =   max(sloss_data$site_size), 
                              "sd patch area" = sd(sloss_data$site_size),
                              "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                              "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



# # plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Montgomery_2014 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Montgomery_2014") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Montgomery_2014") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 3000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 3000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


Montgomery_2014 <-   data.frame("study" = c( data_meta[1,2 ]),
                                "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                                "gamma" = max(result$SL$species),
                                "number of patches" =  nrow(sloss_data), 
                                "minimum patch area" =   min(sloss_data$site_size), 
                                "maximum patch area" =   max(sloss_data$site_size), 
                                "sd patch area" = sd(sloss_data$site_size),
                                "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                                "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)


# # plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)


####################################################################################################
####################################################################################################

#  Rocha_2016 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Rocha_2016") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Rocha_2016") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


Rocha_2016 <-   data.frame("study" = c( data_meta[1,2 ]),
                           "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                           "gamma" = max(result$SL$species),
                           "number of patches" =  nrow(sloss_data), 
                           "minimum patch area" =   min(sloss_data$site_size), 
                           "maximum patch area" =   max(sloss_data$site_size), 
                           "sd patch area" = sd(sloss_data$site_size),
                           "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                           "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)


# # plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)


####################################################################################################
####################################################################################################

#  Sridhar_2008 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Sridhar_2008") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Sridhar_2008") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Sridhar_2008 <-   data.frame("study" = c( data_meta[1,2 ]),
                             "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                             "gamma" = max(result$SL$species),
                             "number of patches" =  nrow(sloss_data), 
                             "minimum patch area" =   min(sloss_data$site_size), 
                             "maximum patch area" =   max(sloss_data$site_size), 
                             "sd patch area" = sd(sloss_data$site_size),
                             "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                             "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Stireman_2014 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Stireman_2014") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Stireman_2014") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(13,31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling

# must aggregate at the site level
sampling_effort <- aggregate(. ~site_id, data=data_sampling, sum, na.rm=TRUE)
patch_area <- aggregate(. ~site_id, data=data_sampling, min, na.rm=TRUE)
sampling_and_patch_area <- cbind(sampling_effort, patch_area)
sampling_and_patch_area <- sampling_and_patch_area[,c(1, 2,6)]
data_sampling <- sampling_and_patch_area

data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

sloss_data_abundance <- aggregate(. ~site_id, data=sloss_data, sum, na.rm=TRUE)
sloss_data_abundance[,2] <- data_sampling[,3]
sloss_data <- sloss_data_abundance
# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

sloss_data$site_id<- as.factor(sloss_data$site_id)
row.names(sloss_data) <- sloss_data$site_id

Stireman_2014 <-   data.frame("study" = c( data_meta[1,2 ]),
                              "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                              "gamma" = max(result$SL$species),
                              "number of patches" =  nrow(sloss_data), 
                              "minimum patch area" =   min(sloss_data$site_size), 
                              "maximum patch area" =   max(sloss_data$site_size), 
                              "sd patch area" = sd(sloss_data$site_size),
                              "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                              "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)


# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Ulrich_2016 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Ulrich_2016") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Ulrich_2016") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Ulrich_2016 <-   data.frame("study" = c( data_meta[1,2 ]),
                            "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                            "gamma" = max(result$SL$species),
                            "number of patches" =  nrow(sloss_data), 
                            "minimum patch area" =   min(sloss_data$site_size), 
                            "maximum patch area" =   max(sloss_data$site_size), 
                            "sd patch area" = sd(sloss_data$site_size),
                            "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                            "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



# # plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Vallan_2000 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Vallan_2000") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Vallan_2000") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Vallan_2000 <-   data.frame("study" = c( data_meta[1,2 ]),
                            "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                            "gamma" = max(result$SL$species),
                            "number of patches" =  nrow(sloss_data), 
                            "minimum patch area" =   min(sloss_data$site_size), 
                            "maximum patch area" =   max(sloss_data$site_size), 
                            "sd patch area" = sd(sloss_data$site_size),
                            "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                            "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)


# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Wang_2012 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Wang_2012") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Wang_2012") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 1000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 1000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Wang_2012 <-   data.frame("study" = c( data_meta[1,2 ]),
                          "comparison" = c("inconclusive"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                          "gamma" = max(result$SL$species),
                          "number of patches" =  nrow(sloss_data), 
                          "minimum patch area" =   min(sloss_data$site_size), 
                          "maximum patch area" =   max(sloss_data$site_size), 
                          "sd patch area" = sd(sloss_data$site_size),
                          "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                          "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

#  Williams_2011 ----
## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Williams_2011") 
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Williams_2011") 
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts 
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset 
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao


Williams_2011 <-   data.frame("study" = c( data_meta[1,2 ]),
                              "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive 
                              "gamma" = max(result$SL$species),
                              "number of patches" =  nrow(sloss_data), 
                              "minimum patch area" =   min(sloss_data$site_size), 
                              "maximum patch area" =   max(sloss_data$site_size), 
                              "sd patch area" = sd(sloss_data$site_size),
                              "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                              "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

###Cayuela_2006

## subset the data
data <- subset(data_fragsad, data_fragsad$refshort == "Cayuela_2006")
data_meta <- subset(metadata_fragsad2, metadata_fragsad2$refshort == "Cayuela_2006")
data_meta2 <- data_meta[, c(13,15)]
data_sampling <- data_meta[, c(31,15)]

## check if larger sites are sample at much lower efforts
data_sampling
data_sampling <- subset(data_sampling, data_sampling$site_size < 60000)
data_sampling$effort_per_ha <- as.numeric(data_sampling$sampling_effort)/as.numeric(data_sampling$site_size)
data_sampling$effort_norm_per_ha <- (data_sampling$effort_per_ha / max(data_sampling$effort_per_ha)) #/ data_sampling$site_size

effort <- min(data_sampling$effort_norm_per_ha)/max(data_sampling$effort_norm_per_ha)
effort_sd <- sd(data_sampling$effort_norm_per_ha)

## reshape data to create a table of species per site
data <- data[,c(3,5,6)]
data2 <- reshape(data, idvar = "site_id", timevar = "scientific_name", direction = "wide")
data2[is.na(data2)] = 0
sloss_data <- merge.data.frame( data_meta2, data2, by = 'site_id')
sloss_data$site_size <-as.numeric(sloss_data$site_size)
sloss_data <- subset(sloss_data, sloss_data$site_size < 60000)

# create occurrence dataset
sloss_data_occu <- sloss_data
sloss_data_occu[,3:ncol(sloss_data_occu)][sloss_data_occu[,3:ncol(sloss_data_occu)]>0] <- 1

# richness per site
richness_site <- data.frame(richness = rowSums(sloss_data_occu[,3:(ncol(sloss_data_occu))]))

# empirical sloss comparison
result <- sloss(sloss_data[,3:(ncol(sloss_data))], area = sloss_data$site_size)

# number of species in the study
diversity_study <- specpool(sloss_data[,3:(ncol(sloss_data))])
study_area_extrapolated_richness <- diversity_study$chao

Cayuela_2006 <-   data.frame("study" = c( data_meta[1,2 ]),
                             "comparison" = c("SS > SL"), # whether the sloss comparison was SS > SL, SL > SS, or inconclusive
                             "gamma" = max(result$SL$species),
                             "number of patches" =  nrow(sloss_data),
                             "minimum patch area" =   min(sloss_data$site_size),
                             "maximum patch area" =   max(sloss_data$site_size),
                             "sd patch area" = sd(sloss_data$site_size),
                             "effort cor" = cor((rowSums(sloss_data[,3:(ncol(sloss_data))])/sloss_data$site_size), sloss_data$site_size),
                             "patch evenness" = diversity(data_sampling$site_size)/log(nrow(data_sampling))
)



# plots
plot.SLOSS(result, show.legend=TRUE, show.index = FALSE, main = "", ylab = "Species richness", xlab = "Area (ha)")
hist(sloss_data$site_size, main = "", ylab = "Number of patches", xlab = "Patch area (ha)")

# correlation patch size and density
abundance <- data.frame(size= sloss_data$site_size, abundance=rowSums(sloss_data[,3:(ncol(sloss_data))]))
sampl_effort <- data.frame(size= data_sampling$site_size, effort=data_sampling$sampling_effort)
correlation_data <- merge(abundance, sampl_effort)
correlation_data$density <- (correlation_data$abundance/correlation_data$effort)
cor.test(correlation_data$size, correlation_data$density)

####################################################################################################
####################################################################################################

################################################
#### analysis
################################################

# assemble dataset
database_analysis <- rbind(Baldi_1999, Bell_2006_a, Bell_2006_b, Berg_1997, Bernard_2007, Bolger_1997,
                           Bossart_2016, Bossart_2006, Chiarello_1999, Dami_2012, Dauber_2006, Ding_2013,
                           Gavish_2012_a, Gavish_2012_b, Gavish_2012_c, Guadagnin_2005, Henry_2007,
                           Jauker_2019, Kapoor_2008, MacNally_2001, Manu_2007, Martensen_2012, McCollin_1993,
                           Montgomery_2014, Rocha_2016, Sridhar_2008, Stireman_2014, Ulrich_2016,
                           Vallan_2000, Wang_2012, Williams_2011, Cayuela_2006)

# preparing the table for analysis
database_analysis$comparison <- as.factor(database_analysis$comparison)

study_characteristics <- metadata_fragsad2[, c(2, 5:12)]
study_characteristics <- aggregate(. ~refshort, data=study_characteristics, min, na.rm=TRUE)
study_characteristics <- subset(study_characteristics, refshort %in% database_analysis$refshort)

database_analysis <- merge(database_analysis, study_characteristics, by ="refshort")
database_analysis$log.minimum.patch.area <- log(database_analysis$minimum.patch.area)
database_analysis$log.maximum.patch.area <- log(database_analysis$maximum.patch.area)

database_analysis$comparison_n <- database_analysis$comparison
levels(database_analysis$comparison_n)[levels(database_analysis$comparison_n)=="inconclusive"] <- "0"
levels(database_analysis$comparison_n)[levels(database_analysis$comparison_n)=="SS > SL"] <- "1"
database_analysis$comparison_n <- as.numeric(database_analysis$comparison_n)
database_analysis$comparison_n[database_analysis$comparison_n == 1] <- 0
database_analysis$comparison_n[database_analysis$comparison_n == 2] <- 1

#### model fitting and plotting

# model including only sampling bias covariate
fit <- glmmTMB(comparison_n ~  effort.cor, 
               data=database_analysis,
               family="binomial")
summary(fit)
AICc(fit)

# model including both sampling bias and patch size evenness
fit <- glmmTMB(comparison_n ~  effort.cor + patch.evenness, 
               data=database_analysis,
               family="binomial")
summary(fit)
AICc(fit)


# plot model predictions
database_analysis$order = findInterval(database_analysis$effort.cor, sort(database_analysis$effort.cor))
pal = colorRampPalette(c("grey90", "grey0"))
# Specify font
windowsFonts(A = windowsFont("Times New Roman"))  


plotdat <- data.frame(patch.evenness=seq(0,1,0.01),
                      effort.cor = mean(database_analysis$effort.cor)) 

preddat <- predict(fit, newdata=plotdat, se.fit=TRUE)
with(database_analysis, plot(patch.evenness, comparison_n, type="n", 
                             ylim=c(0, 1), ylab="Probability of SS > SL", xlab="Patch size evenness", family = "A", cex.lab =1.5, cex.axis = 1.5))
with(preddat, lines(seq(0,1,0.01), exp(fit)/(1+exp(fit)), col="red", lwd = 2))
with(preddat, lines(seq(0,1,0.01), exp(fit+1.96*se.fit)/(1+exp(fit+1.96*se.fit)), col = "darkred", lty=2))
with(preddat, lines(seq(0,1,0.01), exp(fit-1.96*se.fit)/(1+exp(fit-1.96*se.fit)), col = "darkred", lty=2))
points(database_analysis$patch.evenness, 
       database_analysis$comparison_n,
       pch = 19,
       col = pal(nrow(database_analysis))[database_analysis$order]) # color based on sampling bias covariate

# predictions 
mean(exp(preddat$fit)/(1+exp(preddat$fit)))
mean(exp(preddat$fit-1.96*preddat$se.fit)/(1+exp(preddat$fit-1.96*preddat$se.fit)))
mean(exp(preddat$fit+1.96*preddat$se.fit)/(1+exp(preddat$fit+1.96*preddat$se.fit)))

min(exp(preddat$fit)/(1+exp(preddat$fit)))
min(exp(preddat$fit-1.96*preddat$se.fit)/(1+exp(preddat$fit-1.96*preddat$se.fit)))
min(exp(preddat$fit+1.96*preddat$se.fit)/(1+exp(preddat$fit+1.96*preddat$se.fit)))

max(exp(preddat$fit)/(1+exp(preddat$fit)))
max(exp(preddat$fit-1.96*preddat$se.fit)/(1+exp(preddat$fit-1.96*preddat$se.fit)))
max(exp(preddat$fit+1.96*preddat$se.fit)/(1+exp(preddat$fit+1.96*preddat$se.fit)))

################################################
#### plots and models for supplementary material
################################################
hist(database_analysis$number.of.patches, main = "", xlab = "Number of patches", breaks = 20)
hist(database_analysis$gamma, main = "", xlab = "Gamma diversity", breaks = 20)
hist(log10(database_analysis$minimum.patch.area), main = "", xlab = "log10 of minimum patch area", breaks = 20)
hist(log10(database_analysis$maximum.patch.area), main = "", xlab = "log10 of maximum patch area", breaks = 20)

plot(database_analysis$number.of.patches ~ database_analysis$patch.evenness, 
     main = "", 
     ylab = "Number of patches",
     xlab = "Evenness")

plot(database_analysis$gamma ~ database_analysis$patch.evenness, 
     main = "", 
     ylab = "Gamma diversity",
     xlab = "Evenness")

plot(log10(database_analysis$minimum.patch.area) ~ database_analysis$patch.evenness, 
     main = "", 
     ylab = "log10 minimum patch area",
     xlab = "Evenness")

plot(log10(database_analysis$maximum.patch.area) ~ database_analysis$patch.evenness, 
     main = "", 
     ylab = "log10 maximum patch area",
     xlab = "Evenness")


# Plots for covariate of interest; replace maximum patch area with the other covariates of interest
plot(log10(database_analysis$maximum.patch.area) ~ database_analysis$patch.evenness, xlab = "Patch size evenness", ylab = "Log10 maximum patch area")
hist(database_analysis$minimum.patch.area, main = "Minimum patch area across datasets", xlab = "Logarithm of patch area", breaks = 20)

# correlation between patch size evenness and other covariates of interest
cor(database_analysis$patch.evenness, database_analysis$gamma, method=c("pearson"))
cor(database_analysis$patch.evenness, database_analysis$number.of.patches, method=c("pearson"))
cor(database_analysis$patch.evenness, log10(database_analysis$minimum.patch.area), method=c("pearson"))
cor(database_analysis$patch.evenness, log10(database_analysis$maximum.patch.area), method=c("pearson"))

# models including other covariates of interest 
fit <- glmmTMB(comparison_n ~  effort.cor + patch.evenness + number.of.patches, 
               data=database_analysis,
               family="binomial")
summary(fit)
AICc(fit)

fit <- glmmTMB(comparison_n ~  effort.cor + patch.evenness + gamma, 
               data=database_analysis,
               family="binomial")
summary(fit)
AICc(fit)


fit <- glmmTMB(comparison_n ~  effort.cor + patch.evenness + log(minimum.patch.area), 
               data=database_analysis,
               family="binomial")
summary(fit)
AICc(fit)

fit <- glmmTMB(comparison_n ~  effort.cor + patch.evenness + log(maximum.patch.area), 
               data=database_analysis,
               family="binomial")
summary(fit)
AICc(fit)

fit <- glmmTMB(comparison_n ~  effort.cor +  biome, 
               data=database_analysis,
               family="binomial")
summary(fit)
AICc(fit)

