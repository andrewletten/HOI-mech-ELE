#######################################################################################################
## HOI-funcs-generic: misc functions for simulating, fitting and plotting #############################
## Letten and Stouffer ################################################################################
## The mechanistic basis for higher-order interactions and non-additivity in competitive communities ##
## Ecology Letters ####################################################################################
#######################################################################################################


### mech.multi.ode
mech.multi.ode = function(Time, State, Pars, essential = FALSE, chemo = TRUE)
{
  with(as.list(c(State, Pars)),
       {
         nconsumers = Pars$nconsumers
         N = State[1:nconsumers]
         R = State[(1+nconsumers):length(State)]
         
         # Consumer dynamics
         dN.perR = mu 
         dN = N
         for(i in 1:length(N)){
           for(j in 1:length(R)){
             dN.perR[i,j] = (mu[i,j] * N[i] * (R[j])^(2*type3[i,j]))/((Ks[i,j])^(2*type3[i,j]) + phi[i,j]*(R[j])^(2*type3[i,j])) 
           }
           if(essential == TRUE){
             dN[i] = min(dN.perR[i,]) - (all_d * N[i])
           } else {
             dN[i] = sum(dN.perR[i,]) - (all_d * N[i])
           }
         }
         
         # Resource dynamics
         dR.perN = mu
         dR = R
         if(essential == TRUE){
           for(j in 1:length(R)){
             for(i in 1:length(N)){
               dR.perN[i,] = (min(dN.perR[i,]))*Qs[i,]
             }
             if (chemo == TRUE){
               dR[j] = res.speed*(S[j]-R[j]) - sum(dR.perN[,j])
             } else{
               dR[j] = (r.res[j]*R[j]*(1-(R[j]/cc[j]))) - sum(dR.perN[,j]) 
             }
           }
         }
         else {
           for(j in 1:length(R)){
             for(i in 1:length(N)){
               dR.perN[i,j] = dN.perR[i,j]*Qs[i,j]
             }
             if (chemo == TRUE){
               dR[j] = res.speed*(S[j]-R[j]) - sum(dR.perN[,j])
             } else{
               dR[j] = (r.res[j]*R[j]*(1-(R[j]/cc[j]))) - sum(dR.perN[,j]) 
             }
           }
         }
         return(list(c(dN, dR)))
       })
}

### Lotka-Volterra funcs

LV.HOI = function(Time, State, Pars)
{
  N = State
  r = Pars$r
  alphas = Pars$alphas
  betas = Pars$betas
  NN = matrix(kronecker(N,N),length(N),length(N))
  dN = N
  for(i in 1:length(State)){
    dN[i] = r[i] * N[i] * (1 + alphas[i,] %*% N + sum((betas[,,i] * NN)[lower.tri(NN,diag=TRUE)]))
  }
  return(list(c(dN)))
}

#########################################################################
# ODE plotiing funcs


baseplot.ode.HOI = function(out.df, Cxlim = NULL, spnum = 2, resnum = 2){
  par(mfrow=c(1,2))
  # Convert pop sizes less than 0.0001 to zero.
  allspnum = c(1:spnum)
  allspids = paste0("N", allspnum)
  for (i in 1:length(allspids)){
    out.df[allspids[i]][out.df[allspids[i]] < 0.01 & out.df[allspids[i]] > 0] = 0
  }
  allreslets = letters[c(1:resnum)]
  allresids = paste0("R", allreslets)
  matplot(out.df[,"time"], out.df[,allspids], log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),lty=1,las=1,bty='l', xlim = Cxlim)
  matplot(out.df[,"time"], out.df[,allresids],type="l",xlab="Time",ylab="Resource abundance",col=c("lightgreen", "black", "yellowgreen"),lty=1,las=1,bty='l', xlim = Cxlim)
  par(mfrow = c(1,1))
}

baseplot.ode.LV = function(out.LV, Cxlim = NULL){
  par(mfrow=c(1,1))
  # Convert pop sizes less than 0.0001 to zero.
  matplot(out.LV[,1],(out.LV[,2:ncol(out.LV)]), log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),lty=1,las=1,bty='l', xlim = Cxlim)
}



baseplot.ode.LV.mech = function(out.LV, out.df, Cxlim = NULL){
  par(mfrow=c(1,2))
  # Convert pop sizes less than 0.0001 to zero.
  allspnum = c(1:spnum)
  allspids = paste0("N", allspnum)
  for (i in 1:length(allspids)){
    out.df[allspids[i]][out.df[allspids[i]] < 0.01 & out.df[allspids[i]] > 0] = 0
  }
  # allreslets = letters[c(1:resnum)]
  # allresids = paste0("R", allreslets)
  matplot(out.LV[,"time"], out.LV[,allspids], log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",
          col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),
          lty=1,las=1,bty='l', xlim = Cxlim)
  matplot(out.df[,"time"], out.df[,allspids], log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",
          col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),
          lty=1,las=1,bty='l', xlim = Cxlim)
  #matplot(out.df[,"time"], out.df[,allresids],type="l",xlab="Time",ylab="Resource abundance",col=c("lightgreen", "black", "yellowgreen"),lty=1,las=1,bty='l', xlim = Cxlim)
  par(mfrow = c(1,1))
}

baseplot.ode.multiLV.mech = function(datasave=datasave, out.LV, out.df, Cxlim = NULL, Cylim = NULL){
  if(datasave$spnum > 1){
  par(mfrow=c(1,5))
  # Convert pop sizes less than 0.0001 to zero.
  allspnum = c(1:datasave$spnum)
  allspids = paste0("N", allspnum)
  # for (i in 1:length(allspids)){
  #   out.df[allspids[i]][out.df[allspids[i]] < 0.01 & out.df[allspids[i]] > 0] = 0
  # }
  # allreslets = letters[c(1:resnum)]
  # allresids = paste0("R", allreslets)
  matplot(out.df[,"time"], out.df[,allspids], log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",
          col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),
          lty=1,las=1,bty='l', xlim = Cxlim, ylim = Cylim, main = "mechanistic model")
  matplot(out.LV[[1]][,"time"], out.LV[[1]][,allspids], log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",
          col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),
          lty=1,las=1,bty='l', xlim = Cxlim, ylim = Cylim, main = "alphas only")
  matplot(out.LV[[2]][,"time"], out.LV[[2]][,allspids], log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",
          col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),
          lty=1,las=1,bty='l', xlim = Cxlim, ylim = Cylim, main = "alphas + betas (inters only)")
  matplot(out.LV[[3]][,"time"], out.LV[[3]][,allspids], log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",
          col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),
          lty=1,las=1,bty='l', xlim = Cxlim, ylim = Cylim, main = "alphas + betas (intras only)")
  matplot(out.LV[[4]][,"time"], out.LV[[4]][,allspids], log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",
          col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),
          lty=1,las=1,bty='l', xlim = Cxlim, ylim = Cylim, main = "alphas + betas (intras & inters)")
  #matplot(out.df[,"time"], out.df[,allresids],type="l",xlab="Time",ylab="Resource abundance",col=c("lightgreen", "black", "yellowgreen"),lty=1,las=1,bty='l', xlim = Cxlim)
  } else {
    par(mfrow=c(1,3))
    # Convert pop sizes less than 0.0001 to zero.
    allspnum = c(1:datasave$spnum)
    allspids = paste0("N", allspnum)
    for (i in 1:length(allspids)){
      out.df[allspids[i]][out.df[allspids[i]] < 0.01 & out.df[allspids[i]] > 0] = 0
    }
    # allreslets = letters[c(1:resnum)]
    # allresids = paste0("R", allreslets)
    matplot(out.df[,"time"], out.df[,allspids], log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",
            col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),
            lty=1,las=1,bty='l', xlim = Cxlim, ylim = Cylim, main = "mechanistic model")
    matplot(out.LV[[1]][,"time"], out.LV[[1]][,allspids], log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",
            col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),
            lty=1,las=1,bty='l', xlim = Cxlim, ylim = Cylim, main = "alphas only")
    matplot(out.LV[[3]][,"time"], out.LV[[3]][,allspids], log = "y", type="l",xlab="Time",ylab="Consumer abundance (log scale)",
            col=c("blue","red", "yellow", "grey", "black", "orange", "pink", "green"),
            lty=1,las=1,bty='l', xlim = Cxlim, ylim = Cylim, main = "alphas + betas (intras only)")
  }
  par(mfrow = c(1,1))
}

############################################################################

####
# function to fit interaction coefficients within a regression framework 
####


fit.stat.model = function(
  data,
  resp = "1",
  fit.alphas=TRUE,
  fit.intrabetas=FALSE,
  fit.interbetas=FALSE,
  sample.alphas=NULL,
  sample.betas=NULL,
  dynamic.coefs = TRUE,
  prime = FALSE,
  predictors = 1:data$spnum,
  ...
){
  if(dynamic.coefs == TRUE & prime == TRUE){cat("Alphas Saavedra style")}
  # make sure the user has provided an eligible model
  #match.arg(type)
  
  # lets figure out who the observed competitors are
  competitors = paste0("N", predictors, "dens")
  
  
  # start with no competition model
  base.model.formula = paste(paste0("rN", resp), "~ ", sep = " ") 
  
  # in case it's the null
  model.formula = base.model.formula
  
  # add pairwise coefficients for all competitor species
  all.alphas = c()
  if((fit.alphas)){
    # all competitors are allowed an alpha
    all.alphas = competitors
    # if there are alphas to add, add them to the model formula
    if(length(all.alphas)>0){
      model.formula = paste0(model.formula, paste0(all.alphas, collapse=" + "))
    }
  }
  
  # determine potential higher order interactions if so desired and add them to the model specification
  if(fit.intrabetas){
    # betas that have alphas being estimated (principle of marginality)
    possible.intrabetas = all.alphas
    possible.intrabetas = unlist(lapply(possible.intrabetas, function(x){paste0("I(",x,"^2)")}))
    if(length(possible.intrabetas)>0){
      model.formula = paste0(model.formula, " + ", paste0(possible.intrabetas, collapse=" + "))
    }
  }
  
  if(fit.interbetas){
    # betas between heterospecifics
    possible.interbetas <- combn(all.alphas,2)
    possible.interbetas <- apply(possible.interbetas,2,paste,collapse=":")
    if(length(possible.interbetas)>0){
      model.formula = paste0(model.formula, " + ", paste0(possible.interbetas, collapse=" + "))
    }
  }# combine all betas together into a single variable
  
  ######
  # fit the model
  ######
  m = lm(as.formula(model.formula),
         data = data$df,
         ...
  )

  if(dynamic.coefs){
    modcoefs = (m$coefficients)
    LVrs = modcoefs["(Intercept)"]
    names(LVrs) = paste0("r", resp)
    LV.alphas.prime = modcoefs[-c(grep("(Intercept)", names(modcoefs)), grep("\\:", names(modcoefs)), grep("\\^2", names(modcoefs)))]
    alphas.ID = gsub("N(\\d)dens", "\\1", names(LV.alphas.prime), perl = TRUE)
    names(LV.alphas.prime) = paste0("a", resp, alphas.ID)
    LV.alphas = LV.alphas.prime/LVrs
    LV.intra.betas = rep(0, times = length(LV.alphas))
    LV.inter.betas = rep(0, times = length(LV.alphas))
    LV.intra.betas.prime = NA
    LV.inter.betas.prime = NA
    
    if(fit.intrabetas){
      LV.intra.betas.prime = modcoefs[grep("\\^2", names(modcoefs))]
      intra.betas.ID = gsub("I\\(N(\\d)dens\\^2\\)", "\\1", names(LV.intra.betas.prime), perl = TRUE)
      names(LV.intra.betas.prime) = paste0("a", resp, paste0(intra.betas.ID, intra.betas.ID))
      LV.intra.betas = LV.intra.betas.prime/LVrs
    } 
    
    if(fit.interbetas){
      LV.inter.betas.prime = modcoefs[grep("\\:", names(modcoefs))]
      inter.betas.ID = gsub("N(\\d)dens\\:N(\\d)dens", "\\1\\2", names(LV.inter.betas.prime), perl = TRUE)
      names(LV.inter.betas.prime) = paste0("a", resp, inter.betas.ID)
      LV.inter.betas = LV.inter.betas.prime/LVrs
    }
    ifelse(prime == FALSE, return(list(LVrs, LV.alphas, LV.intra.betas, LV.inter.betas)), 
           return(list(LVrs, LV.alphas.prime, LV.intra.betas.prime, LV.inter.betas.prime)))
  } else{
    return(m)
  }
}

allfoc.fit.stat.model= function(
  data,
  resp = c("1","2","3"),
  dynamic.coefs = TRUE,
  ...)
{
  allfits = NULL
  r = rep(NA,length(resp))
  alphas = matrix(NA,length(resp),length(resp))
  betas = array(NA,c(length(resp),length(resp),length(resp)))
  if(dynamic.coefs){
    for (i in 1:length(resp)){
    resp.iter = resp[i]
    allfits.perresp = unlist(fit.stat.model(data,
                                            resp = resp.iter,
                                            ...))
    r[i] = allfits.perresp[paste0("r",i)]
    alphas[i,] = allfits.perresp[paste0("a",i,1:length(resp))]
    if(length(resp) > 1){
    diag(betas[,,i]) = allfits.perresp[paste0("a",i,paste0(1:length(resp),1:length(resp)))]
    betas[,,i][lower.tri(betas[,,i])] = allfits.perresp[paste0("a",i,apply(combn(c(1:data$spnum),2),2,paste,collapse=""))] # have a look at this to check if 1,2,3 needs to be genralised using resp
    } else {
    betas[,,i][i] = allfits.perresp[paste0("a",i,paste0(1:length(resp),1:length(resp)))]
    }
    betas[,,i][is.na(betas[,,i])] = 0
  }
  return(list(r = r, alphas = alphas, betas = betas))}
  else{
    allfits.perresp = NULL
    for (i in 1:length(resp)){
      resp.iter = resp[i]
      allfits.perresp[[i]] = fit.stat.model(data,resp = resp.iter, dynamic.coefs = FALSE, ...)
    }
    return(allfits.perresp)}
}


LV.rsquares = function(datasave = datasave, LVpars){
  LV.Rsq = rep(NA,length = datasave$spnum)
  for (i in 1:datasave$spnum){
    LV.Rsq[i] = summary(LVpars[[i]])$r.squared
  }
  return(LV.Rsq)
}

LV.summary = function(datasave = datasave, LVpars){
  LV.sums = list()
  for (i in 1:datasave$spnum){
    LV.sums[[i]] = summary(LVpars[[i]])
  }
  return(LV.sums)
}


############################################################################

#############################################
##### Pulses, forcings  and perturbations ###
#############################################

# Pulse abundances only
eventfun.mortpulse.generic = function(Time, State, Pars, ...){
  with(as.list(State),{
    nconsumers = Pars$nconsumers
    N = Pars$perturb.frac
    R = State[(1+nconsumers):length(State)]
    return(c(N, R))
  })
}

# Pulse resources only
eventfun.respulse = function(Time, State, Pars){
  with(as.list(State),{
    N1 = N1 
    N2 = N2 
    N3 = N3
    R = (a)*R + resource*(1-a)
    return(c(N1, N2, N3, R))
  })
}

################################################
## Plotting functional forms

funcform = function(R = R, mu = mu, K = K, phi = phi, type3 = type3){
    per.cap = mu * (R)^(2*type3)/((K)^(2*type3) + phi*(R)^(2*type3))
    return(per.cap)
}

monod.sim = function(x = R, r = r, k = k){
  per.cap = r*x/(k+x)
  return(per.cap)
}

holling.sim.type3 = function(x = R, r = r, k = k){
  per.cap = r*(x)^2/((k)^2+(x)^2)
  return(per.cap)
}

