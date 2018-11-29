##################################################################################################################
## 02-HOI-fitmodels: Fit (higher-order) Lotka-Volterra models to observed data and resim to compare fits #########
## Letten and Stouffer ###########################################################################################
## The mechanistic basis for higher-order interactions and non-additivity in competitive communities #############
## Ecology Letters ###############################################################################################
##################################################################################################################

library(deSolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Cairo)
library(cowplot)
library(gtools)
library(DescTools)
library(foreach)
library(parallel)
library(tictoc)

source("code/HOI-funcs-generic.R") 
source("code/HOI-fitmod-functions.R") 

modelfit.dat = run_code(folder = "data/", 
                      N = 2, 
                      linear = TRUE, 
                      parallel = FALSE, 
                      show = TRUE, 
                      test = FALSE, 
                      essential = TRUE, 
                      chemo = FALSE,
                      save = TRUE)




