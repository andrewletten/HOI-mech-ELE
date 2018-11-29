#######################################################################################################
## HOI-fitmod-functions: misc functions for fitting models and resimulating ###########################
## Letten and Stouffer ################################################################################
## The mechanistic basis for higher-order interactions and non-additivity in competitive communities ##
## Ecology Letters ####################################################################################
#######################################################################################################

Rss = function(datasave){
  pars = datasave$meta$params
  LVpars.alphas = LVpars.inters = LVpars.intras = LVpars.allbetas = NULL
  Rsq.alpha = Rsq.inter = Rsq.intra = Rsq.allbetas = NULL
  
  tryCatch({
    Rsq.alpha = LV.rsquares(datasave=datasave, LVpars=allfoc.fit.stat.model(datasave, resp = c(1:datasave$spnum), fit.intrabetas = F, fit.interbetas = F, dynamic.coefs = FALSE))
  }, error=function(e) print("unable to fit alphas"))
  tryCatch({
    Rsq.inter = LV.rsquares(datasave=datasave, LVpars=allfoc.fit.stat.model(datasave, resp = c(1:datasave$spnum), fit.intrabetas = F, fit.interbetas = T, dynamic.coefs = FALSE))
  }, error=function(e) print("unable to fit inters"))
  tryCatch({
    Rsq.intra = LV.rsquares(datasave=datasave, LVpars=allfoc.fit.stat.model(datasave, resp = c(1:datasave$spnum), fit.intrabetas = T, fit.interbetas = F, dynamic.coefs = FALSE))
  }, error=function(e) print("unable to fit intras"))
  tryCatch({
    Rsq.allbetas = LV.rsquares(datasave=datasave, LVpars=allfoc.fit.stat.model(datasave, resp = c(1:datasave$spnum), fit.intrabetas = T, fit.interbetas = T, dynamic.coefs = FALSE))
  }, error=function(e) print("unable to fit allbetas"))
  
  allrsquares = c(alpha = mean(Rsq.alpha), alphainter = mean(Rsq.inter, na.rm = TRUE), alphaintra = mean(Rsq.intra), betas = mean(Rsq.allbetas, na.rm = TRUE))
  return(allrsquares) 
}

gen_params = function(datasave){
  LVpars.alphas = LVpars.inters = LVpars.intras = LVpars.allbetas = NULL
  
  tryCatch({
    LVpars.alphas = allfoc.fit.stat.model(datasave, resp = c(1:datasave$spnum), fit.intrabetas = F, fit.interbetas = F, dynamic.coefs = TRUE)
  }, error=function(e) print("unable to fit alphas"))
  
  tryCatch({
    LVpars.inters = allfoc.fit.stat.model(datasave, resp = c(1:datasave$spnum), fit.intrabetas = F, fit.interbetas = T, dynamic.coefs = TRUE)
  }, error=function(e) print("unable to fit inters"))
  
  tryCatch({
    LVpars.intras = allfoc.fit.stat.model(datasave, resp = c(1:datasave$spnum), fit.intrabetas = T, fit.interbetas = F, dynamic.coefs = TRUE)
  }, error=function(e) print("unable to fit intras"))
  
  tryCatch({
    LVpars.allbetas = allfoc.fit.stat.model(datasave, resp = c(1:datasave$spnum), fit.intrabetas = T, fit.interbetas = T, dynamic.coefs = TRUE)
  }, error=function(e) print("unable to fit allbetas"))
  
  LVpars.loop = list(alphas = LVpars.alphas, inters = LVpars.inters, intras = LVpars.intras, allbetas = LVpars.allbetas)
  return(LVpars.loop)
}

testR = function(datasave, params){
  #######################################
  # Test run
  time.total = 1500 # 2000 # real length of simulation
  times=round(seq(0.1,time.total,by=0.1),1) # total number of time-steps
  stateLV.ini = c(rep(10, datasave$spnum))
  out.LV=as.data.frame(ode(LV.HOI,
                           y=stateLV.ini,
                           parms=params$alphas,
                           times=times,
                           method = "lsoda"))
  baseplot.ode.LV(out.LV)
  #######################################
}



gen_test = function(datasave, LVpars.loop, chemo, essential, indices, times, show){
  out.LVstore = NULL
  meanerrsN.iter = NULL
  meanerrsN = NULL #data.frame(matrix(nrow = nrow(ini.loop.modelcomp), ncol = (datasave$spnum)))
  meanerrsN.eachLVmod = NULL
  meanerrsN.eachsp = NULL
  meanerrsN.eachini = NULL
  mechvar.loop = NULL
  mechnull.loop = NULL

  ini.loop.modelcomp = datasave$df[,c((datasave$spnum+1):(datasave$spnum+datasave$spnum))]
  ini.loop.modelcomp = as.data.frame(ini.loop.modelcomp)

  names(ini.loop.modelcomp) = paste0("N", 1:datasave$spnum)

  meanerrsN.alphas = data.frame(matrix(nrow = nrow(ini.loop.modelcomp), ncol = (datasave$spnum)))
  meanerrsN.inters = data.frame(matrix(nrow = nrow(ini.loop.modelcomp), ncol = (datasave$spnum)))
  meanerrsN.intras = data.frame(matrix(nrow = nrow(ini.loop.modelcomp), ncol = (datasave$spnum)))
  meanerrsN.allbetas = data.frame(matrix(nrow = nrow(ini.loop.modelcomp), ncol = (datasave$spnum)))
  mechvar = data.frame(matrix(nrow = nrow(ini.loop.modelcomp), ncol = (datasave$spnum)))
  mechnull = data.frame(matrix(nrow = nrow(ini.loop.modelcomp), ncol = (datasave$spnum)))

  for (i in 1:nrow(ini.loop.modelcomp)){
    state.ini = c(unlist(ini.loop.modelcomp[i,]), datasave$meta$params$S)
    out.df=as.data.frame(ode(func=mech.multi.ode,
                             y=state.ini,
                             parms=datasave$meta$params,
                             times=times,
                             essential = essential,
                             chemo = chemo,
                             method = "lsoda",
                             events = NULL))
    names(out.df)[-1] =  c(paste0("N", 1:datasave$spnum), paste0("R", letters[c(1:((ncol(out.df))-(datasave$spnum+1)))]))
    stateLV.ini = unlist(ini.loop.modelcomp[i,])
    for (j in indices){

      LVpars = LVpars.loop[[j]]
      out.LV = as.data.frame(ode(func=LV.HOI,
                                 y=stateLV.ini,
                                 parms=LVpars,
                                 times=times,
                                 method = "lsoda"))
      names(out.LV)[-1] =  c(paste0("N", 1:datasave$spnum))
      out.LVstore[[j]] = out.LV
      for (k in 1:datasave$spnum){
        meanerrsN.eachsp[k] = mean(abs(out.df[,k+1] - out.LV[,k+1]))
        #mechvar.loop[k] = mean(abs(out.df[,k+1] - mean(out.df[,k+1])))
        mechnull.loop[k] = mean(abs(out.df[,k+1] - state.ini[k]))
      }
      meanerrsN.eachLVmod[[j]] = meanerrsN.eachsp
    }
    #meanerrsN.eachini[i,] = unlist(meanerrsN.eachLVmod)
    #mechvar[i,] = mechvar.loop
    mechnull[i,] = mechnull.loop
    if(datasave$spnum > 1){
      meanerrsN.alphas[i,] = meanerrsN.eachLVmod[[1]]
      meanerrsN.inters[i,] = meanerrsN.eachLVmod[[2]]
      meanerrsN.intras[i,] = meanerrsN.eachLVmod[[3]]
      meanerrsN.allbetas[i,] = meanerrsN.eachLVmod[[4]]
    } else {
      meanerrsN.alphas[i,] = meanerrsN.eachLVmod[[1]]
      meanerrsN.intras[i,] = meanerrsN.eachLVmod[[3]]}
    if(show){
      tryCatch({
        baseplot.ode.multiLV.mech(datasave, out.LVstore, out.df, Cylim = c(min(ini.loop.modelcomp), max(ini.loop.modelcomp)*8))
      }, error=function(e){"say blah"})
    }
  }

  meanerrs.df = bind_cols(meanerrsN.alphas, meanerrsN.inters, meanerrsN.intras, meanerrsN.allbetas, mechnull)
  meanerrs.df = meanerrs.df[apply(meanerrs.df, MARGIN = 1, function(x) all(x < 1e5, na.rm = TRUE)), ]
  if(datasave$spnum > 1){
    meanerrs.df = meanerrs.df[complete.cases(meanerrs.df),]
  } else {meanerrs.df = meanerrs.df[complete.cases(meanerrs.df[,c(1,3)]),]
  }
  return(meanerrs.df)
}



main_f = function(datasave, N=2, essential=TRUE, chemo=TRUE, test=TRUE, show=TRUE){
  ######### Define Variables
  if (essential){
    coff= -0.0149
  }else{
    coff= -0.0299
  }
  time.total = 20 # real length of simulation
  times=round(seq(0.1,time.total,by=0.1),1) # total number of time-steps

  bestmod = as.data.frame(matrix(nrow = 4, ncol = 4))
  names(bestmod) = c("AIC", "loglikelihood", "pseudo.Rsquare", "numparam")
  rownames(bestmod) = c("alphas only", "alphas and inters", "alphas and intras", "alphas, intras and inters")
  indices = 1:4
  if (datasave$spnum==1){
    datasave$df = datasave$df[(datasave$df$rN1 > coff),]
    indices = c(1,3)
  }else if(datasave$spnum==2){
    indices = 1:4
    datasave$df = datasave$df[(datasave$df$rN1 > coff) & (datasave$df$rN2 > coff),]
  }else if(datasave$spnum==3){
    indices = 1:4
    datasave$df = datasave$df[(datasave$df$rN1 > coff) & (datasave$df$rN2 > coff) & (datasave$df$rN3 > coff),]
  }else{
    indices = 1:4
    datasave$df = datasave$df[(datasave$df$rN1 > coff) & (datasave$df$rN2 > coff) & (datasave$df$rN3 > coff) & (datasave$df$rN4 > coff),]
  }

  ######### training data

  allrsquares = Rss(datasave=datasave)
  LVpars.loop = gen_params(datasave)

  if(show){
    testR(datasave=datasave, params=LVpars.loop)
  }


  ######### test data

  meanerrs.df = gen_test(datasave=datasave,  LVpars.loop=LVpars.loop, chemo=chemo, essential=essential, indices=indices, times=times, show=show)

  meanerrsN.list = list(meanerrsN.alphas = meanerrs.df[,1:(datasave$spnum)],
                        meanerrsN.inters = meanerrs.df[,(1+(datasave$spnum)):((datasave$spnum)*2)],
                        meanerrsN.intras = meanerrs.df[,(1+(datasave$spnum)*2):((datasave$spnum)*3)],
                        meanerrsN.allbetas = meanerrs.df[,(1+(datasave$spnum)*3):((datasave$spnum)*4)])


  ######## Reformat data

  for (i in 1:4){
    meanerr.df = meanerrsN.list[[i]]
    if (N == 1){
      AICwithoutcorrection = (length(meanerr.df))*(log((sum((meanerr.df), na.rm = TRUE)/length(meanerr.df))))
      Rsquare.pseudo = 1 - ((sum((meanerr.df))))/sum((meanerrs.df[,(1+(datasave$spnum)*4):((datasave$spnum)*5)]))
    } else {
      AICwithoutcorrection = (nrow(meanerr.df))*(log((sum((rowMeans(meanerr.df)), na.rm = TRUE)/nrow(meanerr.df))))
      Rsquare.pseudo = 1 - ((sum((rowMeans(meanerr.df)), na.rm = TRUE)))/sum(rowMeans(as.data.frame(meanerrs.df[,(1+(datasave$spnum)*4):((datasave$spnum)*5)])))
      }
    #AICwithoutcorrection = nrow(meanerr.df)*(log(sum((rowMeans(meanerr.df))^2, na.rm = TRUE)/nrow(meanerr.df)))
    #Rsquare = 1 - ((sum((rowMeans(meanerr.df)), na.rm = TRUE))/nrow(meanerr.df))/var(rowMeans(ini.loop.modelcomp))
    numberLVparams = length(unlist(LVpars.loop[i])[unlist(LVpars.loop[i]) != 0])
    AIC = AICwithoutcorrection + 2*numberLVparams
    bestmod[i,] = c(AIC, AICwithoutcorrection, Rsquare.pseudo, numberLVparams)
  }

  bestmod$traindata.Rsq = allrsquares

  return(list(savebestmod=bestmod, ndat=nrow(meanerrs.df)))
}




run_code = function(folder, N, linear=TRUE, essential=TRUE, chemo=TRUE, test=TRUE, parallel=FALSE, show=FALSE, save = TRUE){
  if (linear){
    filename = paste(rep("lin", times = N), collapse = "")
  }else{
    filename = paste(rep("nonlin", times = N), collapse = "")
  }
  filename = paste(folder, filename, "-", N, "sp", sep="")
  if (chemo){
    filename = paste(filename, "chemo", sep="-")
  }else{
    filename = paste(filename, "logis", sep="-")
  }
  if (essential){
    filename = paste(filename, "essent.Rdata", sep="-")
  }else{
    filename = paste(filename, "subs.Rdata", sep="-")
  }

  if(!file.exists(filename)){
    print(filename)
    stop(paste("The file named ", filename, " doesn't exist, mate", sep="'"))
  }
  print(filename)
  load(file = filename)
  savebestmod = NULL
  ndat = NULL
  tic(filename)
  if (parallel){
    no_cores = detectCores() - 1
    cl = makeCluster(6, type="FORK")
    res = parLapply(cl, 1:length(percap.save.all),
              function(l){
                datasave = percap.save.all[[l]]
                if (is.null(datasave)){
                  return(list(savebestmod=NULL, ndat=NULL))
                }
                pres = main_f(datasave, N=N, essential=essential, chemo=chemo, show=FALSE)
                return(pres)
            })
    stopCluster(cl)
    savebestmod = lapply(res, function(x) x$savebestmod)
    ndat = lapply(res, function(x) x$ndat)
  }else{
    pb = txtProgressBar(min = 1, max = length(percap.save.all)+1, style = 3)
    for (l in 1:length(percap.save.all)){
      datasave = percap.save.all[[l]]
      if (is.null(datasave)){
        next
      }
      setTxtProgressBar(pb, l)
      res = main_f(datasave, N=N, essential=essential, chemo=chemo, show=show)

      ndat[l] = res$ndat
      savebestmod[[l]] = res$savebestmod
    }
    setTxtProgressBar(pb, length(percap.save.all)+1)
    close(pb)
  }
  if (save == TRUE){
    modelfit.dat = list(savebestmod=savebestmod, ndat=ndat)
    if (linear){
      savefilename = paste(rep("lin", times = N), collapse = "")
    }else{
      savefilename = paste(rep("nonlin", times = N), collapse = "")
    }
    savefilename = paste(folder, savefilename, "-", N, "sp", sep="")
    if (chemo){
      savefilename = paste(savefilename, "chemo", sep="-")
    }else{
      savefilename = paste(savefilename, "logis", sep="-")
    }
    if (essential){
      savefilename = paste(savefilename, "essent-modelfits.Rdata", sep="-")
    }else{
      savefilename = paste(savefilename, "subs-modelfits.Rdata", sep="-")
    }
    save(modelfit.dat, file = savefilename)
  }
  return(list(savebestmod=savebestmod, ndat=ndat))
  print(toc())
}