#######################################################################################################
## 03-HOI-summarizefits: Summarize static and dynamic fits (Figs 3-8 in ms) ###########################
## Letten and Stouffer ################################################################################
## The mechanistic basis for higher-order interactions and non-additivity in competitive communities ##
## Ecology Letters ####################################################################################
#######################################################################################################

library(ggplot2)

cbPalette <- c( "#56B4E9", "#D55E00", "#009E73", "#E69F00",  "#F0E442",  "#CC79A7")
dodge = position_dodge(.35)

# change folder below to local model fits destination
allmodelfits = grep(".Rdata", list.files("data/modelfits"), value = TRUE)

allcoefs = list()
for (j in 1:length(allmodelfits)){
  load(file = paste0("data/modelfits/",allmodelfits[j]))
  savebestmod = modelfit.dat[[1]]
  pseudo.df = data.frame(matrix(ncol = 8, nrow = length(savebestmod)))
  for (i in 1:length(savebestmod)){
    if (is.null(savebestmod[[i]])){
      next
    }
    pseudo.df[i,] = c(savebestmod[[i]]$pseudo.Rsquare[c(1:4)], savebestmod[[i]]$traindata.Rsq[c(1:4)])
  }
  names(pseudo.df) = c("test.alpha", "test.alphainter", "test.alphaintra", "test.betas", "train.alpha", "train.alphainter", "train.alphaintra", "train.betas")
  quant95Rsq = sapply(pseudo.df, quantile, c(0.025,0.975), na.rm = TRUE)
  seRsq = sapply(pseudo.df,function(x)sd(x, na.rm = TRUE)/sqrt(length(x)))
  meanRsq = colMeans(pseudo.df, na.rm = TRUE)
  coef.df = data.frame(LV = names(meanRsq), meanRsq, lowerSE = meanRsq - 1.96*seRsq, upperSE = meanRsq + 1.96*seRsq)
  coef.df$traintest = c("test", "test", "test", "test", "train", "train", "train", "train")
  coef.df$coef = c("alpha", "alpha.inter","alpha.intra", "betas", "alpha", "alpha.inter","alpha.intra", "betas")
  coef.df$islinear = strsplit(allmodelfits[[j]], "-")[[1]][1]
  coef.df$spnum = strsplit(allmodelfits[[j]], "-")[[1]][2]
  coef.df$ischemo = strsplit(allmodelfits[[j]], "-")[[1]][3]
  coef.df$isessent = strsplit(allmodelfits[[j]], "-")[[1]][4]
  coef.df$islinear[coef.df$islinear == "lin"] = "Linear"
  coef.df$islinear[coef.df$islinear == "nonlin"] = "Nonlinear"
  coef.df$islinear[coef.df$islinear == "linlin"] = "Linear"
  coef.df$islinear[coef.df$islinear == "nonlinnonlin"] = "Nonlinear"
  coef.df$islinear[coef.df$islinear == "linlinlin"] = "Linear"
  coef.df$islinear[coef.df$islinear == "nonlinnonlinnonlin"] = "Nonlinear"
  coef.df$coeflab = coef.df$coef
  coef.df$coeflab[coef.df$coef == "alpha.inter"] = "...+ beta (inter)"
  coef.df$coeflab[coef.df$coef == "alpha.intra"] = "...+ beta (intra)"
  coef.df$coeflab[coef.df$coef == "betas"] = "...+ beta (inter + intra)"
  alphastd.v0 = coef.df$meanRsq[coef.df$coef == "alpha"]
  coef.df$alphastd = rep(alphastd.v0, each = 4)
  coef.df$meanRsq.std = coef.df$meanRsq -  coef.df$alphastd
  coef.df$lowerSE.std = coef.df$meanRsq.std - (coef.df$meanRsq -  coef.df$lowerSE)
  coef.df$upperSE.std = coef.df$meanRsq.std + (coef.df$upperSE - coef.df$meanRsq)
  coef.df$statdyn = coef.df$traintest
  coef.df$statdyn[coef.df$statdyn == "train"] = "Static"
  coef.df$statdyn[coef.df$statdyn == "test"] = "Dynamic"
  allcoefs[[j]] = coef.df
}

allcoefs.df = do.call("rbind", allcoefs)
allcoefs.df$coeflab = factor(allcoefs.df$coeflab)
allcoefs.df$coeflab = factor(allcoefs.df$coeflab, levels = levels(allcoefs.df$coeflab)[c(4,2,3,1)])
allcoefs.df$ischemo = factor(allcoefs.df$ischemo)
allcoefs.df$ischemo = factor(allcoefs.df$ischemo, levels = levels(allcoefs.df$ischemo)[c(2,1)])
allcoefs.df$statdyn = factor(allcoefs.df$statdyn, levels = unique(allcoefs.df$statdyn)[c(2,1)])
levels(allcoefs.df$ischemo) = c("Logistic", "Constant")
allcoefs.df$spnum[allcoefs.df$spnum == "1sp"] = "n = 1"
allcoefs.df$spnum[allcoefs.df$spnum == "2sp"] = "n = 2"
allcoefs.df$spnum[allcoefs.df$spnum == "3sp"] = "n = 3"

####################################
## Multi-species (figs 4-7 in ms) ##
####################################


coef.df.sub = allcoefs.df[(allcoefs.df$isessent == "essent" & allcoefs.df$ischemo == "Logistic") & allcoefs.df$spnum != "n = 1",]
# For other scenarios `$isessent == "subs"' and/or "$ischemo == "Constant"'
alphatestR2 = coef.df.sub[coef.df.sub$LV == "train.alpha" & coef.df.sub$traintest == "train",]

comboplot = ggplot(coef.df.sub, aes(y = meanRsq.std, x = statdyn)) +
  geom_hline(yintercept=0, lty=2, lwd=0.2, colour="black", alpha = 1) +
  geom_errorbar(aes(ymin=lowerSE.std, ymax=upperSE.std, colour = coef), lwd=1, width=0, position = dodge) +
  geom_point(aes(shape = traintest, fill = coef), size=3, position = dodge) +
  #coord_cartesian(ylim= c(0.5,1)) + 
  facet_grid(islinear ~ spnum) +
  theme_bw() +
  theme(legend.title=element_blank(),
        legend.position = c(0.15, 0.9),
        legend.text=element_text(size=10),
        legend.background = element_rect(fill=NA),
        text=element_text(size=16, family="Helvetica"),
        axis.title.x=element_blank(),
        strip.background = element_blank(),
        axis.text.x=element_text(size = 10, hjust = 0.5, vjust = 0.5, colour = "black"),
        axis.text.y=element_text(colour = "black"),
        plot.background = element_rect(fill = "white"),
        #axis.line.x = element_line(colour = "white"),
        #axis.line.y = element_line(colour = "white"),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black")) +
  #guides(colour=FALSE) +
  labs(y=expression(paste(Delta, italic(R^{2})))) +
  scale_fill_manual(values=cbPalette, guide = FALSE) +
  scale_colour_manual(values=cbPalette, labels=c(bquote(alpha[ij]), 
                                                 bquote(alpha[ij] + beta[ijk]), 
                                                 bquote(alpha[ij] + beta[ijj]), 
                                                 bquote(alpha[ij] + beta[ijk] + beta[ijj])), guide = guide_legend(reverse = TRUE)) + 
  scale_shape_manual(values = c(23,21), guide = FALSE) + scale_x_discrete(breaks = NULL) +
  geom_text(data = alphatestR2, 
            aes(label = round(meanRsq, 2), y = -0.01, x = 0.65), size = 3, colour = "#666666")

comboplot

# pdf(file = "figures/ms-rsq-logis-essent.pdf", width = 5, height = 7)
# comboplot
# dev.off()

#####################################
## Single species (figs 4-7 in ms) ##
#####################################

coef.df.sub = allcoefs.df[(allcoefs.df$coef == "alpha" | allcoefs.df$coef == "alpha.intra") & allcoefs.df$spnum == "n = 1",]
alphatestR2 = coef.df.sub[coef.df.sub$LV == "train.alpha" & coef.df.sub$traintest == "train",]

comboplot = ggplot(coef.df.sub, aes(y = meanRsq.std, x = statdyn)) +
  geom_hline(yintercept=0, lty=2, lwd=0.2, colour="black", alpha = 1) +
  geom_errorbar(aes(ymin=lowerSE.std, ymax=upperSE.std, colour = coef), lwd=1, width=0, position = dodge) +
  geom_point(aes(shape = traintest, fill = coef), size=3, position = dodge) +
  #coord_cartesian(ylim= c(0.5,1)) + 
  facet_grid(islinear ~ ischemo) +
  theme_bw() +
  theme(legend.title=element_blank(),
        legend.position = c(0.12, 0.93),
        legend.text=element_text(size=10),
        legend.background = element_rect(fill=NA),
        text=element_text(size=16, family="Helvetica"),
        axis.title.x=element_blank(),
        strip.background = element_blank(),
        axis.text.x=element_text(size = 10, hjust = 0.5, vjust = 0.5, colour = "black"),
        axis.text.y=element_text(colour = "black"),
        plot.background = element_rect(fill = "white"),
        #axis.line.x = element_line(colour = "white"),
        #axis.line.y = element_line(colour = "white"),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black")) +
  #guides(colour=FALSE) +
  labs(y=expression(paste(Delta, italic(R^{2})))) +
  scale_fill_manual(values=cbPalette[c(1,3)], guide = FALSE) +
  scale_colour_manual(values=cbPalette[c(1,3)], labels=c(bquote(alpha[ij]), 
                                                         bquote(alpha[ij] + beta[ijj])), guide = guide_legend(reverse = TRUE)) + 
  scale_shape_manual(values = c(23,21), guide = FALSE) + scale_x_discrete(breaks = NULL) + 
  geom_text(data = alphatestR2, aes(label = round(meanRsq, 2), y = -0.015, x = 0.65), size = 2, colour = "#666666")

comboplot

# pdf(file = "figures/ms-rsq-1sp.pdf", width = 5, height = 5)
# comboplot
# dev.off()
