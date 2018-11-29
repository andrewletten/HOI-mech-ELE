#######################################################################################################
## HOI-gendat-functions: misc functions for generating 'observed' data ################################
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


source("code/HOI-funcs-generic.R") # model functions

make_par_list = function(spnum = 2, resnum = 2, linear = TRUE, essential = FALSE){
  pars = list()
  pars$mu = matrix(rep(0.1, times = spnum*resnum), nrow = spnum, byrow = TRUE) # updated in loop
  pars$Ks = matrix(rep(1, times = spnum*resnum), nrow = spnum, byrow = TRUE) # updated in loop
  pars$Qs = pars$mu # updated in loop
  if(linear == TRUE){
    pars$phi= matrix(rep(0, times = spnum*resnum), nrow = spnum, byrow = TRUE)
  } else {
    pars$phi= matrix(rep(1, times = spnum*resnum), nrow = spnum, byrow = TRUE)
  }
  pars$type3 = matrix(rep(1/2, times = spnum*resnum), nrow = spnum, byrow = TRUE)
  if(essential == TRUE){
    pars$all_d = 0.015 # 0.03 for substitutable, 0.015 for essential
  } else{
    pars$all_d = 0.03
  }
  pars$nconsumers = spnum
  pars$res.speed = 50 # resource supply rate (chemo) # in original submission this was set to 10 (sped up for consistent pop sizes compared to logistic growth)
  pars$S = rep(1, times = resnum) # resource concentration (chemo) / carrying capacity (logis)
  pars$r.res = rep(100, times = resnum) # resource intrinsic rate of increase (logis)
  pars$cc = pars$S # logis
  pars$spnum = spnum
  pars$resnum = resnum
  return(pars)
  }

get_perms = function(essential = FALSE, chemo = FALSE, pars, times, state.ini, perturb = 10, ccfrac = 1){
  stateN = state.ini[1:pars$spnum]
  ccN = NULL
  for(l in 1:length(stateN)){
    cc.state.ini = state.ini
    cc.state.ini[1:pars$spnum][-l]  = 0
    out.df = as.data.frame(ode(func=mech.multi.ode,
                               y=cc.state.ini,
                               parms=pars,
                               times=times,
                               essential = essential,
                               chemo = chemo,
                               method = "lsoda"))
    cc =  tail(out.df)[1,l+1]
    ccN[l] = cc
  }
  print(paste("CCs", ccN))
  minseq = c(seq(1, min(ccN*ccfrac), length.out = perturb))
  maxseq = seq(1, max(ccN*ccfrac), length.out = perturb) 
  percapfrac.seq = c(minseq, maxseq[maxseq > min(ccN*ccfrac)])
  percapfrac.loop = data.frame(permutations(n = length(percapfrac.seq), r = pars$spnum, v = percapfrac.seq))
  names(percapfrac.loop) = stateN
  for (i in 1:pars$spnum){
    percapfrac.loop = percapfrac.loop[percapfrac.loop[i] < (((ccN[i])*ccfrac)+1),] 
  }
  percapfrac.loop = as.data.frame(percapfrac.loop)
  print(nrow(percapfrac.loop))
  return(percapfrac.loop)
}

test_sim = function(essential = FALSE, chemo = FALSE, pars, times, state.ini, loopnum){
  testrun.df=as.data.frame(ode(func=mech.multi.ode,
                             y=state.ini,
                             parms=pars,
                             times=times,
                             essential = essential,
                             chemo = chemo,
                             method = "lsoda"))
  names(testrun.df)[-1] =  c(paste0("N", 1:pars$nconsumers), paste0("R", letters[c(1:((ncol(testrun.df))-(pars$nconsumers+1)))]))
  # loop.sim[[k]] = tail(testrun.df)[1,]
  print(paste("Random param simulation:", loopnum))
  return(testrun.df)
  }

perturb_fun = function(essential = FALSE, chemo = FALSE, percapfrac.loop, pars, times, state.ini, loopnum){
  percap.temp = NULL
  denslevels = NULL
  percap.save = data.frame(matrix(nrow = nrow(percapfrac.loop), ncol = (pars$spnum*2 + pars$resnum)))
  for (i in 1:nrow(percapfrac.loop)){
  perturb.frac = unlist(percapfrac.loop[i,])
  pars$perturb.frac = perturb.frac
  # Run numerical simulations
  pulse = 1000 
  out.df=as.data.frame(ode(func=mech.multi.ode,
                           y=state.ini,
                           parms=pars,
                           times=times,
                           essential = essential,
                           chemo = chemo,
                           method = "lsoda", events = list(func = eventfun.mortpulse.generic, time = pulse)))
  percapest = out.df[(pulse*10+-10):(pulse*10+300),] #change back to 300
  names(percapest)[-1] =  c(paste0("N", 1:pars$nconsumers), paste0("R", letters[c(1:((ncol(out.df))-(pars$nconsumers+1)))]))
  baseplot.ode.HOI(percapest, spnum = pars$spnum, resnum = pars$resnum, Cxlim = NULL)
  percapest.maxres = percapest[21:30,] # change back to 30 for main sims
  
  for (j in 1:(pars$spnum)){
    fitsN = smooth.spline(x = percapest.maxres[,1], y = log(percapest.maxres[,j+1]))
    percap.temp[j] = predict(fitsN,median(percapest.maxres$time),deriv=1)$y # (log(percapest.maxres[100,j+1]) - log(percapest.maxres[1,j+1]))/(percapest.maxres$time[100] - percapest.maxres$time[1]) # 
    denslevels[j] = percapest.maxres[6,j+1] # percapest.maxres[21,j+1] #  500 for midpoint
  }
  percap.save[i,] = c(percap.temp, denslevels, reslevels = percapest.maxres[6,(pars$spnum+2):(pars$spnum+pars$resnum+1)]) # c(percap.temp, denslevels, reslevels = percapest.maxres[6,(pars$spnum+2):(pars$spnum+pars$resnum+1)])
  print(paste("Random param simulation:", loopnum, "perturb", i))
  }
  return(percap.save)
}


main_sim_func = function(essential = FALSE, chemo = FALSE, linear = TRUE, spnum = 2, resnum = 2, time.total = 1500, perturb = 10, ccfrac = 1){
  percap.save.all = list()
  #################################################################
  ### Make parameters
  pars = make_par_list(essential = essential, spnum = spnum, resnum = resnum, linear = linear)
  if (linear == TRUE){
    paramloop = matrix(runif(spnum*resnum*100, 0.015,0.1), ncol = spnum*resnum, byrow = TRUE) #matrix(runif(spnum*resnum*100, pars$all_d,0.1), ncol = spnum*resnum, byrow = TRUE)
  } else {
    paramloop = matrix(runif(spnum*resnum*100, 0.01,5), ncol = spnum*resnum, byrow = TRUE)
  }
  for (k in c(1:nrow(paramloop))){
  tryCatch({
    loopnum = k
    state.ini = c(rep(10, pars$spnum), pars$S)
    if(linear == TRUE){
      pars$mu = matrix(paramloop[k,], nrow = pars$spnum, byrow = TRUE)
      if(spnum > 2){
        pars$Qs = 100*t(apply(pars$mu, 1, sample))
      } else {
        if(essential == TRUE){
          pars$Qs = 100*Rev(pars$mu, margin = 2)
        } else{
          pars$Qs = 100*pars$mu
        }
      }
    } else {
      pars$Ks = matrix(paramloop[k,], nrow = spnum, byrow = TRUE) 
      if(spnum > 2){
        pars$Qs = t(apply(pars$Ks, 1, sample))
      } else {
        if(essential == TRUE){
          pars$Qs = 1*pars$Ks
        } else{
          pars$Qs = Rev(pars$Ks, margin = 2)
        }
      }
    }
    #################################################################
    ### Run full simulation for each paramterisation and print figure 
    times=round(seq(0.1,time.total,by=0.1),1) # total number of time-steps
    testrun.df = test_sim(essential = essential, chemo = chemo, pars, times, state.ini, loopnum)
    source("code/HOI-plots-generic.R", echo = FALSE, print = TRUE, local = TRUE) #
    #################################################################
    ### Get monoculture carrying capacities and generate table of density permutations 
    percapfrac.loop = get_perms(essential = essential, chemo = chemo, pars, times, state.ini, perturb, ccfrac)
    #################################################################
    percap.save = perturb_fun(essential = essential, chemo = chemo, percapfrac.loop, pars, times, state.ini, loopnum)
    names(percap.save) = c(paste0("rN", 1:pars$spnum), paste0("N", 1:pars$spnum, "dens"), paste0("R", letters[c(1:pars$resnum)]))
    #################################################################
    datasave = list()
    datasave$fulldynamics = testrun.df
    datasave$df = percap.save
    datasave$meta = list()
    datasave$meta$state.ini = state.ini
    datasave$spnum = pars$spnum
    spparams = pars 
    datasave$meta$params = pars
    datasave$pulsedf = percapfrac.loop
    percap.save.all[[k]] = datasave
  }, error=function(e){"say blah"})
}
  return(percap.save.all)
}

gendat_save = function(folder, essential = FALSE, chemo = FALSE, linear = TRUE, spnum = 1, resnum = 1, time.total = 1500, perturb = 10, save = TRUE, ccfrac = 1){
  percap.save.all = main_sim_func(essential = essential, 
                                  chemo = chemo, 
                                  linear = linear,
                                  spnum = spnum, 
                                  resnum = resnum,
                                  time.total = time.total,
                                  perturb = perturb,
                                  ccfrac = ccfrac)
  if (linear){
    filename <- paste(rep("lin", times = spnum), collapse = "")
  }else{
    filename <- paste(rep("nonlin", times = spnum), collapse = "")
    }
  filename <- paste(folder, filename, "-", spnum, "sp", sep="")
  if (chemo){
    filename <- paste(filename, "chemo", sep="-")
    }else{
      filename <- paste(filename, "logis", sep="-")
  } 
  if (essential){
    filename <- paste(filename, "essent.Rdata", sep="-") # check file ending!!!!!
    }else{
      filename <- paste(filename, "subs.Rdata", sep="-")
    }
  if (save == TRUE){save(percap.save.all, file = filename)}
  return(percap.save.all)
}
