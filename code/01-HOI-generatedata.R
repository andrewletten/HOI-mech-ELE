#######################################################################################################
## 01-HOI-generatedata: Generate 'observed' data from simulations of consumer-resource models #########
## Letten and Stouffer ################################################################################
## The mechanistic basis for higher-order interactions and non-additivity in competitive communities ##
## Ecology Letters ####################################################################################
#######################################################################################################

library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Cairo)
library(cowplot)
library(gtools)
library(DescTools)

source("code/HOI-funcs-generic.R") 
source("code/HOI-gendat-functions.R") 

percap.save.all = gendat_save(folder = "data/", 
                              essential = TRUE, 
                              chemo = FALSE, 
                              linear = TRUE, 
                              spnum = 2, 
                              resnum = 2, 
                              time.total = 1500, 
                              save = TRUE, 
                              perturb = 10,
                              ccfrac = 1)
  