######################################################################
## Plotting mechanistic dynamics and mechanistic functional forms ####
######################################################################
cbbPalette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

#########################################

plot.df = testrun.df
## Dynamics
comp.gg = gather(plot.df, state.var, count, -time)
comp.gg$state.var.type = "Consumers"
comp.gg$state.var.type[grep("R", comp.gg$state.var)] = "Resources"
fulldynamics.plot = ggplot(comp.gg, aes(y = count, x= time)) +
  geom_line(aes(group = state.var, col = state.var), size = 1) + scale_y_continuous(trans="log1p") +
  theme(legend.position="none") + 
  ylab("Density") + xlab("Time (hours)") +
  panel_border(colour = "black") + #xlim(0,576) + 
  theme(axis.text = element_text(size = 10),
        axis.title= element_text(size = 10)) + 
  facet_wrap(~ state.var.type, scales = "free") + scale_colour_manual(values=cbbPalette) 

#fulldynamics.plot

#########################################

# pars = datasave$meta$params

## Functional forms
resource.levels = seq(0,max(pars$S),length.out = 1000) # (max(k1.time,k2.time))
resp.iter = data.frame(matrix(nrow = length(resource.levels), ncol = length(pars$S)))
resp.iter$resource.levels = resource.levels
resp.list = list()

for (i in 1:pars$nconsumers){
  for (j in 1:length(pars$S)){
    resp.iter[,j] = funcform(R = resource.levels, 
                            mu = pars$mu[i,j], 
                            K = pars$Ks[i,j], 
                            phi = pars$phi[i,j],
                            type3 = pars$type3[i,j])
    names(resp.iter)[j] = paste0("R", letters[c(1:length(pars$S))[j]])
  }
  resp.iter$sp = paste0("N", c(1:pars$nconsumers)[i])
  resp.list[[i]] = resp.iter 
}

resp.all = do.call("rbind", resp.list)
resp.gg = gather(resp.all, key = resource, value = growth, -resource.levels, -sp)
func.plot = ggplot(resp.gg, aes(y = growth, x = resource.levels)) + 
  geom_line(aes(col = sp, linetype = sp), size = 1) + 
  theme(legend.position="none") + 
  xlab("Resource concentration") + ylab("Per capita growth") +
  panel_border(colour = "black") + geom_hline(yintercept = pars$all_d, linetype = "dashed") + 
  theme(axis.text = element_text(size = 10),
        axis.title= element_text(size = 10)) + 
  coord_cartesian(expand = FALSE) + facet_wrap(~resource) + scale_colour_manual(values=cbbPalette)

#func.plot
p = plot_grid(fulldynamics.plot, func.plot)
p

