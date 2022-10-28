################################################################################
### Formal Analysis -- Model with birds and mammals
################################################################################

### load libraries

library(tidyverse); library(cowplot); library(RColorBrewer); library(rstan); library(brms);library(grid);library(gridExtra); library(loo)

### options for stan/brms

ncores = parallel::detectCores()
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

### define expit function

expit <- function(x){
  out <- exp(x)/(1 + exp(x))
  return(out)
}

### load data 

forage <- read.csv('forage_fracfeed.csv', stringsAsFactors = TRUE)

### modify data to get to only the studies we want to consider

# first drop the studies with unit conversion mistakes for prey densities

################################################################################
### set up for the formal analysis
################################################################################


# drop everything with NA temperatures, prey masses, and predator masses

forage_formal <- forage %>% filter(is.na(Temp_C) == FALSE) %>% filter(is.na(Mass_g) == FALSE) %>% 
  filter(is.na(PredMass_g) == FALSE)

# drop data with very low handling time estimates

forage_formal <- forage_formal %>% filter(log(Fitted_h_day) > log(1e-6))

# find which predator and prey major groups have 15 or more observations

# start with the predators

pred_taxa_n <- forage_formal %>% group_by(Pred_Major_Grouping_1) %>% summarise(Number = n())

pred_drop <- pred_taxa_n$Pred_Major_Grouping_1[which(pred_taxa_n$Number < 15)]

# drop these taxa

forage_formal <- forage_formal[-which(forage_formal$Pred_Major_Grouping_1 %in% pred_drop),]

# now the prey taxa

prey_taxa_n <- forage_formal %>% group_by(Prey_Major_Group_1) %>% summarise(Number = n())

prey_drop <- prey_taxa_n$Prey_Major_Group_1[which(prey_taxa_n$Number < 15)]

# drop these taxa

forage_formal <- forage_formal[-which(forage_formal$Prey_Major_Group_1 %in% prey_drop),]

# drop unnecessary levels from prey and predator major groups

forage_formal$Pred_Major_Grouping_1 <- droplevels(forage_formal$Pred_Major_Grouping_1)

forage_formal$Prey_Major_Group_1 <- droplevels(forage_formal$Prey_Major_Group_1)

# down to 1528 observations

### prepare some variables for the analysis

# create a version of dimension that is a factor

forage_formal <- forage_formal %>% mutate(Dim_factor = as.factor(paste(Dim)))

forage_formal$Dim_factor <- relevel(forage_formal$Dim_factor, ref = '3')

# create a combined habitat variable 

forage_formal <- forage_formal %>% mutate(Habitat_Aquatic = as.factor(paste(Habitat, Fresh_Marine)))

################################################################################
### Now perform the formal analysis with the saturation index at the 
### median predicted prey density
################################################################################

fit <- brm(frac_50 ~ Prey_Major_Group_1 + Pred_Major_Grouping_1 + Dim_factor + Habitat_Aquatic +
             log(PredMass_g) + log(Mass_g) + Temp_C + I(Temp_C^2) + (1|Prey_Major_Group_2) + (1|Pred_Major_Grouping_2),
           data = forage_formal, family = Beta(link = 'logit'), inits = '0', cores = ncores)

plot(fit)

summary(fit, prob = 0.9)

fixef(fit, probs = c(0.05, 0.95))


################################################################################
### Plots for the formal analysis
################################################################################

################################################################################
### plots to show relationships between covariates and saturation
################################################################################

### plots for quantitative variables

# get residuals

residuals_fit <- as.data.frame(residuals(fit))

residuals_fit <- residuals_fit %>% mutate(partial_predmass = Estimate + expit(0.04*log(forage_formal$PredMass_g)),
                                          partial_preymass = Estimate + expit(-0.25*log(forage_formal$Mass_g)),
                                          partial_temp = Estimate + expit(0.053*forage_formal$Temp_C + (-0.0018*forage_formal$Temp_C^2)),
                                          predmass = forage_formal$PredMass_g,
                                          preymass = forage_formal$Mass_g,
                                          temp = forage_formal$Temp_C,
                                          prey_major_taxa = forage_formal$Prey_Major_Group_1,
                                          pred_major_taxa = forage_formal$Pred_Major_Grouping_1) 



Pred_Mass_PreyTaxa <- ggplot(data = residuals_fit, aes(x = log(predmass), y = partial_predmass, color = prey_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Prey Taxa') + geom_function(fun = function(x) expit(0.04*x), color = 'black') +
  xlab('log Predator Mass (g)') + ylab('Partial Residual of Saturation Index') + theme(axis.title.y = element_blank())

Pred_Mass_PredTaxa <- ggplot(data = residuals_fit, aes(x = log(predmass), y = partial_predmass, color = pred_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Predator Taxa') + geom_function(fun = function(x) expit(0.04*x), color = 'black') + 
  xlab('log Predator Mass (g)') + ylab('Partial Residual of Saturation Index') + theme(axis.title.y = element_blank())


Prey_Mass_PreyTaxa <- ggplot(data = residuals_fit, aes(x = log(preymass), y = partial_preymass, color = prey_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Prey Taxa') + geom_function(fun = function(x) expit(-0.25*x), color = 'black') +
  xlab('log Prey Mass (g)') + ylab('') + theme(legend.position = 'none')+ theme(axis.title.y = element_blank())

Prey_Mass_PredTaxa <- ggplot(data = residuals_fit, aes(x = log(preymass), y = partial_preymass, color = pred_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Predator Taxa') + geom_function(fun = function(x) expit(-0.25*x), color = 'black') +
  xlab('log Prey Mass (g)') + ylab('') + theme(legend.position = 'none')+ theme(axis.title.y = element_blank())

Temp_PreyTaxa <- ggplot(data = residuals_fit, aes(x = temp, y = partial_temp, color = prey_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Prey Taxa') + geom_function(fun = function(x) expit(0.053*x - 0.0018*x^2), color = 'black') +
  xlab('Temperature (C)') + ylab('') + theme(legend.position = 'none') + theme(axis.title.y = element_blank())

Temp_PredTaxa <- ggplot(data = residuals_fit, aes(x = temp, y = partial_temp, color = pred_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Predator Taxa') + geom_function(fun = function(x) expit(0.05*x - 0.0018*x^2), color = 'black')  + 
  xlab('Temperature (C)') + ylab('') + theme(legend.position = 'none') + theme(axis.title.y = element_blank())


### want to put these together with predator taxa in one column and prey taxa in one column with
### only the legend with the middle figure

PreyQuant <- plot_grid(Prey_Mass_PreyTaxa, Pred_Mass_PreyTaxa, Temp_PreyTaxa, ncol = 1, nrow = 3, align = 'hv', axis = 'tblr', labels = c('E', 'F', 'G'), label_x = 0.65)

Quant_ylabel <- textGrob('Partial Residual of Saturation Index', gp = gpar(fontsize = 15) ,rot = 90)

PreyQuant <- grid.arrange(arrangeGrob(PreyQuant, left = Quant_ylabel))

PredQuant <- plot_grid(Prey_Mass_PredTaxa, Pred_Mass_PredTaxa, Temp_PredTaxa, ncol = 1, nrow = 3, align = 'hv', axis = 'tblr', labels = c('H', 'I', 'J'), label_x = 0.65)

PredQuant <- grid.arrange(arrangeGrob(PredQuant, left = Quant_ylabel))

QuantPlots <- plot_grid(PreyQuant, PredQuant, ncol = 2, nrow = 1)

### Now for the discrete factors

# prey taxa effects

fit_prey_effect <- data.frame(Taxa = as.factor(c('Intercept', levels(forage_formal$Prey_Major_Group_1)[-1])),
                              Effect = fixef(fit)[1:10,1],
                              Upper = fixef(fit, probs = c(0.05, 0.95))[1:10,4],
                              Lower = fixef(fit, probs = c(0.05, 0.95))[1:10,3])

fit_prey_effect$Taxa <- factor(fit_prey_effect$Taxa, levels = c('Intercept', 'Amphibian',
                                                                'Arachnid', 'Crustacean',
                                                                'Fish', 'Insect', 'Mammal',
                                                                'Mollusk', 'Protozoan',
                                                                'Rotifer')) 

PreyTaxaEffect <- ggplot(data = fit_prey_effect, aes(x = Taxa, y = Effect, color = Taxa)) + geom_point() + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), color = 'black') + geom_point(size = 3) + theme_cowplot() + 
  scale_color_brewer(palette = 'Paired') + 
  geom_hline(yintercept = 0, color = 'gray') + geom_vline(xintercept = 1.5, color = 'gray', linetype = 'dashed') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') +
  xlab('Prey Taxa') + ylab(' Partial Effect (logit scale)') + theme(axis.title.y = element_blank())

# pred taxa effects

fit_pred_effect <- data.frame(Taxa = as.factor(levels(forage_formal$Pred_Major_Grouping_1)),
                              Effect = c(0,fixef(fit)[11:19,1]),
                              Upper = c(0,fixef(fit, probs = c(0.05, 0.95))[11:19,4]),
                              Lower = c(0,fixef(fit, probs = c(0.05, 0.95))[11:19,3]))

PredTaxaEffect <- ggplot(data = fit_pred_effect, aes(x = Taxa, y = Effect, color = Taxa)) + geom_point() + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), color = 'black') + geom_point(size = 3) + theme_cowplot() + 
  scale_color_brewer(palette = 'Paired') + 
  geom_hline(yintercept = 0, color = 'gray') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') +
  xlab('Predator Taxa') + ylab(' Partial Effect (logit scale)') + theme(axis.title.y = element_blank())

# habitat effect

fit_habitat_effect <- data.frame(Habitat = as.factor(levels(forage_formal$Habitat_Aquatic)),
                                 Effect = c(0,fixef(fit)[22:23,1]),
                                 Upper = c(0,fixef(fit, probs = c(0.05, 0.95))[22:23,4]),
                                 Lower = c(0,fixef(fit, probs = c(0.05, 0.95))[22:23,3]))

HabitatEffect <- ggplot(data = fit_habitat_effect, aes(x = Habitat, y = Effect)) + geom_point(size = 3) +
  theme_cowplot() + geom_hline(yintercept = 0, color = 'gray') + xlab('Habitat') + ylab('Effect (logit scale)') + geom_errorbar(aes(ymin = Lower, ymax = Upper)) +
  theme(axis.title.y = element_blank())

# dimension effect

fit_dimension_effect <- data.frame(Dimension = as.factor(levels(forage_formal$Dim_factor)),
                                   Effect = c(0,fixef(fit)[20:21,1]),
                                   Upper = c(0,fixef(fit, probs = c(0.05, 0.95))[20:21,4]),
                                   Lower = c(0,fixef(fit, probs = c(0.05, 0.95))[20:21,3]))

DimensionEffect <- ggplot(data = fit_dimension_effect, aes(x = Dimension, y = Effect)) + geom_point(size = 3) +
  theme_cowplot() + geom_hline(yintercept = 0, color = 'gray') + xlab('Dimension') + ylab('Effect (logit scale)') + geom_errorbar(aes(ymin = Lower, ymax = Upper)) + 
  theme(axis.title.y = element_blank())

# put plots together

QualPlots <- plot_grid(PreyTaxaEffect, PredTaxaEffect, HabitatEffect, DimensionEffect, ncol = 1, nrow = 4, axis = 'tblr', align = 'v', rel_heights = c(1.6,1.6,1,1), labels = c('A', 'B', 'C', 'D'), label_x = 0.9)

QualPlot_yaxis <- textGrob('Partial Effect (logit scale)', gp = gpar(fontsize = 15) ,rot = 90)

QualPlots <- grid.arrange(arrangeGrob(QualPlots, left = QualPlot_yaxis))

### put qualitative and quantitative plots together

TogetherPlots <- plot_grid(QualPlots, QuantPlots, nrow = 1, rel_widths = c(1,2.2))


save_plot(filename = 'CovariatesPlot_BirdsandMammals.png', plot = TogetherPlots, ncol =2.4, nrow = 2.4)










