################################################################################
### Formal Analysis of FoRAGE and figures
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

drop_units <- c(237, 565:568, 592, 640:645, 675, 691,825, 1226:1243, 1322:1341, 1470,
                1494, 1586:1588, 1593:1596, 1638, 1644, 1694, 1791:1794, 1908, 2277, 2308, 2309,
                2368, 2584:2593)

forage <- forage %>% filter(!(Data.set %in% drop_units))


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

### need to create a column for arena size

forage_arena <- forage_formal %>% filter(!is.na(X3D_arena_size_cm3) | !is.na(X2D_Arena_Size_cm2))

forage_arena <- forage_arena %>% filter(Field != 'Y')

forage_arena <- forage_arena %>% mutate(Arena_size = ifelse(Dim == 3, X3D_arena_size_cm3,
                                                            ifelse(Dim == 2, X2D_Arena_Size_cm2,
                                                                   ifelse(Dim == 2.5 & !is.na(X3D_arena_size_cm3), X3D_arena_size_cm3, X2D_Arena_Size_cm2))))

forage_arena <- forage_arena %>% filter(Pred_Major_Grouping_1 != 'Mammal')

forage_arena$Pred_Major_Grouping_1 <- droplevels(forage_arena$Pred_Major_Grouping_1)

forage_arena$Prey_Major_Group_1 <- droplevels(forage_arena$Prey_Major_Group_1)

forage_arena <- forage_arena %>% filter(!(is.na(Arena_size)))

################################################################################
# want to make a plot of the empirical cumulative distribution functions 
# at each of the quantiles of prey densities
################################################################################

### just the studies used in the formal analysis

arena_plot_data <- data.frame(Percentile = as.factor(rep(c('10th', '20th', '30th',
                                                         '40th', '50th', '60th',
                                                         '70th', '80th', '90th'), each = length(forage_arena[,1]))),
                              Saturation = c(forage_arena$frac_10, forage_arena$frac_20,
                                             forage_arena$frac_30, forage_arena$frac_40,
                                             forage_arena$frac_50, forage_arena$frac_60,
                                             forage_arena$frac_70, forage_arena$frac_80,
                                             forage_arena$frac_90))

Reduced_Saturation_Plot <- ggplot(data = arena_plot_data, aes(x = Saturation, color = Percentile)) +
  geom_hline(yintercept = 0.5, color = 'gray') + stat_ecdf(pad = FALSE, size = 1) +
  theme_cowplot() + scale_color_brewer(palette = 'Paired', name = 'Prey\nDensity\nPercentile') + 
  xlab('Saturation Index') + ylab('Cumulative Proportion of\nFunctional Responses') + 
  ggtitle('Reduced Dataset')  + theme(plot.title = element_text(hjust = 0.5)) 

### plots with all of the data 

full_plot_data <- data.frame(Percentile = as.factor(rep(c('10th', '20th', '30th',
                                                          '40th', '50th', '60th',
                                                          '70th', '80th', '90th'), each = length(forage[,1]))),
                             Saturation = c(forage$frac_10, forage$frac_20,
                                            forage$frac_30, forage$frac_40,
                                            forage$frac_50, forage$frac_60,
                                            forage$frac_70, forage$frac_80,
                                            forage$frac_90))


Full_Saturation_Plot <- ggplot(data = full_plot_data, aes(x = Saturation, color = Percentile)) + 
  geom_hline(yintercept = 0.5, color = 'gray') + stat_ecdf(pad = FALSE, size = 1) +
  theme_cowplot() + scale_color_brewer(palette = 'Paired', name = 'Prey\nDensity\nPercentile') + 
  xlab('Saturation Index') + ylab('Cumulative Proportion of\nFunctional Responses') + 
  ggtitle('Full Dataset') + theme(plot.title = element_text(hjust = 0.5), legend.position = 'none') 

### example histograms from the full data set

# 10th percentile

full_10_hist <- ggplot(data = filter(arena_plot_data, Percentile == '10th'), aes(x = Saturation)) + 
  geom_histogram(bins = 30, fill = '#a6cee3') + theme_cowplot() + 
  xlab('Saturation Index') + ylab('Count') + ggtitle('10th Percentile\nReduced Dataset') +
  theme(plot.title = element_text(hjust = 0.5)) 

# 50th percentile

full_50_hist <- ggplot(data = filter(arena_plot_data, Percentile == '50th'), aes(x = Saturation)) + 
  geom_histogram(bins = 30, fill = '#fb9a99') + theme_cowplot() + 
  xlab('Saturation Index') + ylab('Count') + ggtitle('50th Percentile\nReduced Dataset') +
  theme(plot.title = element_text(hjust = 0.5)) 

### 90th percentile

full_90_hist <- ggplot(data = filter(arena_plot_data, Percentile == '90th'), aes(x = Saturation)) + 
  geom_histogram(bins = 30, fill = '#cab2d6') + theme_cowplot() + 
  xlab('Saturation Index') + ylab('Count') + ggtitle('90th Percentile\nReduced Dataset') +
  theme(plot.title = element_text(hjust = 0.5)) 

### put saturation plots together

Saturation_Upper <- plot_grid(Full_Saturation_Plot, Reduced_Saturation_Plot,
                             nrow = 1, ncol = 2, labels = 'AUTO', axis = 'tblr', align = 'hv')

Saturation_Lower <- plot_grid(full_10_hist, full_50_hist, full_90_hist, 
                              nrow = 1, ncol = 3, labels = c('C', 'D', 'E'))

Saturation_Plot <- plot_grid(Saturation_Upper, Saturation_Lower, ncol = 1, nrow = 2)

save_plot(filename = 'Saturation_Plot.png', plot = Saturation_Plot, ncol = 2, nrow = 2)

################################################################################
### Now perform the formal analysis with the saturation index at the 
### median predicted prey density
################################################################################

fit <- brm(frac_50 ~ Prey_Major_Group_1 + Pred_Major_Grouping_1 + Dim_factor + Habitat_Aquatic + log(Arena_size) +
                           log(PredMass_g) + log(Mass_g) + Temp_C + I(Temp_C^2) + (1|Prey_Major_Group_2) + (1|Pred_Major_Grouping_2),
                         data = forage_arena, family = Beta(link = 'logit'), inits = '0', cores = ncores)

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

residuals_fit <- as.data.frame(residuals(fit, type = 'ordinary'))

residuals_fit <- residuals_fit %>% mutate(partial_predmass = Estimate + expit(0.04*log(forage_arena$PredMass_g)),
                                          partial_preymass = Estimate + expit(-0.26*log(forage_arena$Mass_g)),
                                          partial_temp = Estimate + expit(0.0534*forage_arena$Temp_C + (-0.0018*forage_arena$Temp_C^2)),
                                          partial_arena = Estimate + expit(0.0816*log(forage_arena$Arena_size)),
                                          predmass = forage_arena$PredMass_g,
                                          preymass = forage_arena$Mass_g,
                                          temp = forage_arena$Temp_C,
                                          arena_size = forage_arena$Arena_size,
                                          prey_major_taxa = forage_arena$Prey_Major_Group_1,
                                          pred_major_taxa = forage_arena$Pred_Major_Grouping_1) 



Pred_Mass_PreyTaxa <- ggplot(data = residuals_fit, aes(x = log(predmass), y = partial_predmass, color = prey_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Prey Taxa') + geom_function(fun = function(x) expit(0.04*x), color = 'black') +
  xlab('log Predator Mass (g)') + ylab('Partial Residual of Saturation Index') + theme(axis.title.y = element_blank())

Pred_Mass_PredTaxa <- ggplot(data = residuals_fit, aes(x = log(predmass), y = partial_predmass, color = pred_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Predator Taxa') + geom_function(fun = function(x) expit(0.04*x), color = 'black') + 
  xlab('log Predator Mass (g)') + ylab('Partial Residual of Saturation Index') + theme(axis.title.y = element_blank())


Prey_Mass_PreyTaxa <- ggplot(data = residuals_fit, aes(x = log(preymass), y = partial_preymass, color = prey_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Prey Taxa') + geom_function(fun = function(x) expit(-0.26*x), color = 'black') +
  xlab('log Prey Mass (g)') + ylab('') + theme(legend.position = 'none')+ theme(axis.title.y = element_blank())

Prey_Mass_PredTaxa <- ggplot(data = residuals_fit, aes(x = log(preymass), y = partial_preymass, color = pred_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Predator Taxa') + geom_function(fun = function(x) expit(-0.245*x), color = 'black') +
  xlab('log Prey Mass (g)') + ylab('') + theme(legend.position = 'none')+ theme(axis.title.y = element_blank())

Temp_PreyTaxa <- ggplot(data = residuals_fit, aes(x = temp, y = partial_temp, color = prey_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Prey Taxa') + geom_function(fun = function(x) expit(0.054*x - 0.0018*x^2), color = 'black') +
  xlab(expression(Temperature~(degree*C))) + ylab('') + theme(legend.position = 'none') + theme(axis.title.y = element_blank())

Temp_PredTaxa <- ggplot(data = residuals_fit, aes(x = temp, y = partial_temp, color = pred_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Predator Taxa') + geom_function(fun = function(x) expit(0.054*x - 0.0018*x^2), color = 'black')  + 
  xlab(expression(Temperature~(degree*C))) + ylab('') + theme(legend.position = 'none') + theme(axis.title.y = element_blank())

Arena_PreyTaxa <- ggplot(data = residuals_fit, aes(x = log(arena_size), y = partial_arena, color = prey_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Prey Taxa') + geom_function(fun = function(x) expit(0.0817*x), color = 'black') +
  xlab(expression(log~Arena~Size~(cm^2~or~cm^3))) + ylab('') + theme(legend.position = 'none') + theme(axis.title.y = element_blank())

Arena_PredTaxa <- ggplot(data = residuals_fit, aes(x = log(arena_size), y = partial_arena, color = pred_major_taxa)) + geom_point() + 
  theme_cowplot() + scale_colour_brewer(palette = 'Paired', name = 'Prey Taxa') + geom_function(fun = function(x) expit(0.0817*x), color = 'black') +
  xlab(expression(log~Arena~Size~(cm^2~or~cm^3))) + ylab('') + theme(legend.position = 'none') + theme(axis.title.y = element_blank())
  
  
### want to put these together with predator taxa in one column and prey taxa in one column with
### only the legend with the middle figure

PreyQuant <- plot_grid(Prey_Mass_PreyTaxa, Pred_Mass_PreyTaxa, Temp_PreyTaxa, Arena_PreyTaxa, ncol = 1, nrow = 4, align = 'hv', axis = 'tblr', labels = c('E', 'F', 'G', 'H'), label_x = 0.65)

Quant_ylabel <- textGrob('Partial Residual of Saturation Index', gp = gpar(fontsize = 15) ,rot = 90)

PreyQuant <- grid.arrange(arrangeGrob(PreyQuant, left = Quant_ylabel))

PredQuant <- plot_grid(Prey_Mass_PredTaxa, Pred_Mass_PredTaxa, Temp_PredTaxa, Arena_PredTaxa, ncol = 1, nrow = 4, align = 'hv', axis = 'tblr', labels = c('I', 'J', 'K', 'L'), label_x = 0.65)

PredQuant <- grid.arrange(arrangeGrob(PredQuant, left = Quant_ylabel))

QuantPlots <- plot_grid(PreyQuant, PredQuant, ncol = 2, nrow = 1)

### Now for the discrete factors

# prey taxa effects

fit_prey_effect <- data.frame(Taxa = as.factor(levels(forage_arena$Prey_Major_Group_1)),
                              Effect = c(0, fixef(fit)[2:9,1]),
                              Upper = c(0,fixef(fit, probs = c(0.05, 0.95))[2:9,4]),
                              Lower = c(0,fixef(fit, probs = c(0.05, 0.95))[2:9,3]))

fit_prey_effect$Taxa <- factor(fit_prey_effect$Taxa, levels = c('Algae', 'Amphibian',
                                                                'Arachnid', 'Crustacean',
                                                                'Fish', 'Insect',
                                                                'Mollusk', 'Protozoan',
                                                                'Rotifer')) 

PreyTaxaEffect <- ggplot(data = fit_prey_effect, aes(x = Taxa, y = Effect, color = Taxa)) + geom_point() + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), color = 'black') + geom_point(size = 3) + theme_cowplot() + 
  scale_color_brewer(palette = 'Paired') + 
  geom_hline(yintercept = 0, color = 'gray') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') +
  xlab('Prey Taxa') + ylab(' Partial Effect (logit scale)') + theme(axis.title.y = element_blank())

# pred taxa effects

fit_pred_effect <- data.frame(Taxa = as.factor(levels(forage_arena$Pred_Major_Grouping_1)),
                              Effect = c(0,fixef(fit)[10:16,1]),
                              Upper = c(0,fixef(fit, probs = c(0.05, 0.95))[10:16,4]),
                              Lower = c(0,fixef(fit, probs = c(0.05, 0.95))[10:16,3]))

PredTaxaEffect <- ggplot(data = fit_pred_effect, aes(x = Taxa, y = Effect, color = Taxa)) + geom_point() + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), color = 'black') + geom_point(size = 3) + theme_cowplot() + 
  scale_color_brewer(palette = 'Paired') + 
  geom_hline(yintercept = 0, color = 'gray') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') +
  xlab('Predator Taxa') + ylab(' Partial Effect (logit scale)') + theme(axis.title.y = element_blank())

# habitat effect

fit_habitat_effect <- data.frame(Habitat = as.factor(levels(forage_arena$Habitat_Aquatic)),
                                 Effect = c(0,fixef(fit)[19:20,1]),
                                 Upper = c(0,fixef(fit, probs = c(0.05, 0.95))[19:20,4]),
                                 Lower = c(0,fixef(fit, probs = c(0.05, 0.95))[19:20,3]))

HabitatEffect <- ggplot(data = fit_habitat_effect, aes(x = Habitat, y = Effect)) + geom_point(size = 3) +
  theme_cowplot() + geom_hline(yintercept = 0, color = 'gray') + xlab('Habitat') + ylab('Effect (logit scale)') + geom_errorbar(aes(ymin = Lower, ymax = Upper)) +
  theme(axis.title.y = element_blank())

# dimension effect

fit_dimension_effect <- data.frame(Dimension = as.factor(levels(forage_formal$Dim_factor)),
                                   Effect = c(0,fixef(fit)[17:18,1]),
                                   Upper = c(0,fixef(fit, probs = c(0.05, 0.95))[17:18,4]),
                                   Lower = c(0,fixef(fit, probs = c(0.05, 0.95))[17:18,3]))

DimensionEffect <- ggplot(data = fit_dimension_effect, aes(x = Dimension, y = Effect)) + geom_point(size = 3) +
  theme_cowplot() + geom_hline(yintercept = 0, color = 'gray') + xlab('Dimension') + ylab('Effect (logit scale)') + geom_errorbar(aes(ymin = Lower, ymax = Upper)) + 
  theme(axis.title.y = element_blank())

# put plots together

QualPlots <- plot_grid(PreyTaxaEffect, PredTaxaEffect, HabitatEffect, DimensionEffect, ncol = 1, nrow = 4, axis = 'tblr', align = 'v', rel_heights = c(1.6,1.6,1,1), labels = c('A', 'B', 'C', 'D'), label_x = 0.9)

QualPlot_yaxis <- textGrob('Partial Effect (logit scale)', gp = gpar(fontsize = 15) ,rot = 90)

QualPlots <- grid.arrange(arrangeGrob(QualPlots, left = QualPlot_yaxis))

### put qualitative and quantitative plots together

TogetherPlots <- plot_grid(QualPlots, QuantPlots, nrow = 1, rel_widths = c(1,2.2))


save_plot(filename = 'CovariatesPlot.png', plot = TogetherPlots, ncol =2.4, nrow = 2.4)







