###########################################################################################
### Cleaning the data to examine the nonlinearity of predator functional responses at 
### realistic prey densities
###########################################################################################

### load libraries

library(dplyr); library(ggplot2); library(cowplot); library(brms);

### load bayesian regressions for mass-abundance relationships

load('MassAbundance_Fits.RData')

### need to connect these mass-abundance regressions to FoRAGE to get
### abundance estimates of the prey

### load FoRAGE

forage <- read.csv('FoRAGE_db_V2_3_Oct_26_2022_mod.csv', stringsAsFactors = TRUE)

### drop rows to not include

forage <- forage %>% filter(Include_Biotic_Prey == 1)

### drop levels of Prey Major groups

forage$Prey_Major_Group_1 <- droplevels(forage$Prey_Major_Group_1)

### drop studies without prey masses

forage <- forage %>% filter(Prey_mass_mg > 0)

### add abundance group column

forage <- forage %>% mutate(PreyAbundanceGroup = ifelse(Prey_Vert.Invert== 'Protozoan', 'Protist',
                                                        ifelse(Prey_Major_Group_1 == 'Mammal', 'Mammal',
                                                               ifelse(Prey_Major_Group_1 == 'Bird', 'Bird',
                                                                      ifelse(Prey_Vert.Invert == 'Prokaryote' | Prey_Vert.Invert == 'Algae', 'Prokaryote', 'Ectotherm')))))  

### add columns for masses of predators and prey in grams

forage$Mass_g <- forage$Prey_mass_mg/1000

forage$PredMass_g <- forage$Pred_mass_mg/1000

Abundance <- matrix(nrow = 2164, ncol = 9)

for(i in 1:2164) {
  
  if(forage$PreyAbundanceGroup[i] == 'Protist') {
    newdat <-  forage[i,]
    Abundance[i,] <-  posterior_predict(object = protist_fit, newdata = newdat, draws = 1000) %>% 
      quantile(probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  } else if(forage$PreyAbundanceGroup[i] == 'Mammal') {
    newdat <-  forage[i,]
    Abundance[i,] <- posterior_predict(object = mammal_fit, newdata = newdat, draws = 1000) %>% 
      quantile(probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  } else if(forage$PreyAbundanceGroup[i] == 'Bird') {
    newdat <- forage[i,]
    Abundance[i,] <- posterior_predict(object = bird_fit, newdata = newdat, draws = 1000) %>% 
      quantile(probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  } else if(forage$PreyAbundanceGroup[i] == 'Ectotherm'){
    newdat <- forage[i,]
    Abundance[i,] <- posterior_predict(object = ectotherm_fit, newdata = newdat, draws = 1000) %>% 
      quantile(probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  } else {
    newdat <- forage[i,]
    Abundance[i,] <- posterior_predict(object = prokaryote_fit, newdata = newdat, draws = 1000) %>% 
      quantile(probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
  }
  
  print(i)
  
}


Abundance <- exp(Abundance)

Abundance <- as.data.frame(Abundance)

colnames(Abundance) <- c('PreyAbundance_10', 'PreyAbundance_20', 'PreyAbundance_30',
                         'PreyAbundance_40', 'PreyAbundance_50', 'PreyAbundance_60',
                         'PreyAbundance_70', 'PreyAbundance_80', 'PreyAbundance_90')

### now need to add Abundance to forage

forage <- forage %>% bind_cols(Abundance)

### now we want to calculate the fraction feeding at each estimated prey density

### write a quick function to calculate the feeding rate saturation index

fracfeed <- function(a, h, Robs){
  
  frac <- a*Robs*h/(1 + a*Robs*h)
  
  frac
  
}

### calculate feeding rate saturation for each predator-prey pair

frac_10 <- vector(length = length(forage[,1]))

for(i in 1:length(forage[,1])) {
  
  frac_10[i] <- fracfeed(a = forage$Fitted_a[i], h = forage$Fitted_h[i], Robs = forage$PreyAbundance_10[i])
  
}

frac_20 <- vector(length = length(forage[,1]))

for(i in 1:length(forage[,1])) {
  
  frac_20[i] <- fracfeed(a = forage$Fitted_a[i], h = forage$Fitted_h[i], Robs = forage$PreyAbundance_20[i])
  
}

frac_30 <- vector(length = length(forage[,1]))

for(i in 1:length(forage[,1])) {
  
  frac_30[i] <- fracfeed(a = forage$Fitted_a[i], h = forage$Fitted_h[i], Robs = forage$PreyAbundance_30[i])
  
}

frac_40 <- vector(length = length(forage[,1]))

for(i in 1:length(forage[,1])) {
  
  frac_40[i] <- fracfeed(a = forage$Fitted_a[i], h = forage$Fitted_h[i], Robs = forage$PreyAbundance_40[i])
  
}

frac_50 <- vector(length = length(forage[,1]))

for(i in 1:length(forage[,1])) {
  
  frac_50[i] <- fracfeed(a = forage$Fitted_a[i], h = forage$Fitted_h[i], Robs = forage$PreyAbundance_50[i])
  
}

frac_60 <- vector(length = length(forage[,1]))

for(i in 1:length(forage[,1])) {
  
  frac_60[i] <- fracfeed(a = forage$Fitted_a[i], h = forage$Fitted_h[i], Robs = forage$PreyAbundance_60[i])
  
}

frac_70 <- vector(length = length(forage[,1]))

for(i in 1:length(forage[,1])) {
  
  frac_70[i] <- fracfeed(a = forage$Fitted_a[i], h = forage$Fitted_h[i], Robs = forage$PreyAbundance_70[i])
  
}

frac_80 <- vector(length = length(forage[,1]))

for(i in 1:length(forage[,1])) {
  
  frac_80[i] <- fracfeed(a = forage$Fitted_a[i], h = forage$Fitted_h[i], Robs = forage$PreyAbundance_80[i])
  
}

frac_90 <- vector(length = length(forage[,1]))

for(i in 1:length(forage[,1])) {
  
  frac_90[i] <- fracfeed(a = forage$Fitted_a[i], h = forage$Fitted_h[i], Robs = forage$PreyAbundance_90[i])
  
}

### combine fraction feeding and forage

forage <- forage %>% bind_cols(frac_10, frac_20, frac_30,
                               frac_40, frac_50, frac_60,
                               frac_70, frac_80, frac_90)

colnames(forage)[59:67] <- c('frac_10', 'frac_20', 'frac_30',
                             'frac_40', 'frac_50', 'frac_60',
                             'frac_70', 'frac_80', 'frac_90')


# ### save this data set as a .csv

write.csv(forage, file = 'forage_fracfeed.csv')

