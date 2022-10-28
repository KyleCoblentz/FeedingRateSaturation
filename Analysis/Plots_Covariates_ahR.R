###############################################################################
### Effects of covariates through a and h
###############################################################################

### load libraries

library(tidyverse); library(cowplot); library(brms); library(grid); library(gridExtra)

### load and manipulate data

### define expit function

expit <- function(x){
  out <- exp(x)/(1 + exp(x))
  return(out)
}

### load data 

forage <- read.csv('forage_fracfeed.csv', stringsAsFactors = TRUE)

### modify data to get to only the studies we want to consider

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

forage_arena <- forage_arena %>% mutate(Habitat_Dimension = paste0(Habitat_Aquatic, Dim_factor))

##################################################################################################
#### go through covariates and make plots for the relationships between the covariate a, h, and R
##################################################################################################

library(brms)

a_fit <-  brm(log(Fitted_a) ~ Prey_Major_Group_1 + Pred_Major_Grouping_1 + Habitat_Aquatic + Dim_factor + log(Arena_size) +
          log(PredMass_g) + log(Mass_g) + Temp_C + I(Temp_C^2) + (1|Prey_Major_Group_2) + (1|Pred_Major_Grouping_2),
        data = forage_arena)

summary(a_fit, prob = 0.9)

h_fit <- brm(log(Fitted_h_day) ~ Prey_Major_Group_1 + Pred_Major_Grouping_1 + Habitat_Aquatic + Dim_factor + log(Arena_size) +
                   log(PredMass_g) + log(Mass_g) + Temp_C + I(Temp_C^2) + (1|Prey_Major_Group_2) + (1|Pred_Major_Grouping_2),
                 data = forage_arena)

summary(h_fit, prob = 0.9)

### put together plots that are similar to the plots that we made for the overall 
### covariate analyses but for each particular functional response parameter

##################################################################################################################
### plots for a 
##################################################################################################################

residuals_a_fit <- as.data.frame(residuals(a_fit))

residuals_a_fit <- residuals_a_fit %>% mutate(partial_predmass = Estimate + 0.35*log(forage_arena$PredMass_g),
                                          partial_preymass = Estimate + 0.04*log(forage_arena$Mass_g),
                                          partial_temp = Estimate + 0.24*forage_arena$Temp_C + -0.01*forage_arena$Temp_C^2,
                                          partial_arena = Estimate + 0.02*log(forage_arena$Arena_size),
                                          predmass = forage_arena$PredMass_g,
                                          preymass = forage_arena$Mass_g,
                                          temp = forage_arena$Temp_C,
                                          arena_size = forage_arena$Arena_size,
                                          prey_major_taxa = forage_arena$Prey_Major_Group_1,
                                          pred_major_taxa = forage_arena$Pred_Major_Grouping_1) 

### make plots 

### quantitative variables

# effect of predator mass

a_predmass_plot <- ggplot(data = residuals_a_fit, aes(x = log(predmass), y = partial_predmass)) + geom_point() + 
  theme_cowplot() + geom_function(fun = function(x) 0.35*x, color = 'black', size = 1) +
  xlab('log Predator Mass (g)') + ylab('Partial Residual of\nSpace Clearance Rate') +   theme(axis.title.y = element_blank())

# effect of prey mass

a_preymass_plot <- ggplot(data = residuals_a_fit, aes(x = log(preymass), y = partial_preymass)) + geom_point() + 
  theme_cowplot() + geom_function(fun = function(x) 0.04*x, color = 'black', size = 1) +
  xlab('log Prey Mass (g)') + ylab('Partial Residual of\nSpace Clearance Rate') +   theme(axis.title.y = element_blank()) 

# effect of temperature

a_temp_plot <- ggplot(data = residuals_a_fit, aes(x = temp, y = partial_temp)) + geom_point() + 
  theme_cowplot() + geom_function(fun = function(x) 0.24*x + -0.01*x^2, color = 'black', size = 1) +
  xlab(expression(Temperature~(degree*C))) + ylab('Partial Residual of\nSpace Clearance Rate') +   theme(axis.title.y = element_blank())

# effect of arena size

a_arena_size_plot <- ggplot(data = residuals_a_fit, aes(x = log(arena_size), y = partial_arena)) + geom_point() + 
  theme_cowplot() + geom_function(fun = function(x) 0.02*x, color = 'black', size = 1) +
  xlab(expression(log~Arena~Size~(cm^2~or~cm^3))) + ylab('Partial Residual of\nSpace Clearance Rate') +   theme(axis.title.y = element_blank())

### discrete variables

# prey taxa effects

a_fit_prey_effect <- data.frame(Taxa = as.factor(levels(forage_arena$Prey_Major_Group_1)),
                              Effect = c(0,fixef(a_fit)[2:9,1]),
                              Upper = c(0, fixef(a_fit, probs = c(0.05, 0.95))[2:9,4]),
                              Lower = c(0,fixef(a_fit, probs = c(0.05, 0.95))[2:9,3]))

a_fit_prey_effect$Taxa <- factor(a_fit_prey_effect$Taxa, levels = c('Algae', 'Amphibian',
                                                                'Arachnid', 'Crustacean',
                                                                'Fish', 'Insect',
                                                                'Mollusk', 'Protozoan',
                                                                'Rotifer')) 

a_PreyTaxaEffect <- ggplot(data = a_fit_prey_effect, aes(x = Taxa, y = Effect, color = Taxa)) + geom_point() + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), color = 'black') + geom_point(size = 3) + theme_cowplot() + 
  scale_color_brewer(palette = 'Paired') + 
  geom_hline(yintercept = 0, color = 'gray') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') +
  xlab('Prey Taxa') + ylab(' Partial Effect on\nSpace Clearance Rate (log scale)') +   theme(axis.title.y = element_blank())

# pred taxa effects

a_fit_pred_effect <- data.frame(Taxa = as.factor(levels(forage_arena$Pred_Major_Grouping_1)),
                              Effect = c(0,fixef(a_fit)[10:16,1]),
                              Upper = c(0,fixef(a_fit, probs = c(0.05, 0.95))[10:16,4]),
                              Lower = c(0,fixef(a_fit, probs = c(0.05, 0.95))[10:16,3]))

a_PredTaxaEffect <- ggplot(data = a_fit_pred_effect, aes(x = Taxa, y = Effect, color = Taxa)) + geom_point() + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), color = 'black') + geom_point(size = 3) + theme_cowplot() + 
  scale_color_brewer(palette = 'Paired') + 
  geom_hline(yintercept = 0, color = 'gray') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') +
  xlab('Predator Taxa') + ylab(' Partial Effect on\nSpace Clearance Rate(log scale)')  +   theme(axis.title.y = element_blank())

# habitat effect

a_fit_habitat_effect <- data.frame(Habitat = as.factor(levels(forage_arena$Habitat_Aquatic)),
                                 Effect = c(0,fixef(a_fit)[17:18,1]),
                                 Upper = c(0,fixef(a_fit, probs = c(0.05, 0.95))[17:18,4]),
                                 Lower = c(0,fixef(a_fit, probs = c(0.05, 0.95))[17:18,3]))

a_HabitatEffect <- ggplot(data = a_fit_habitat_effect, aes(x = Habitat, y = Effect)) + geom_point(size = 3) +
  theme_cowplot() + geom_hline(yintercept = 0, color = 'gray') + xlab('Habitat') + ylab('Partial Effect on\nSpace Clearance Rate\n(log scale)') + geom_errorbar(aes(ymin = Lower, ymax = Upper)) + 
  theme(axis.title.y = element_blank())

# dimension effect

a_fit_dimension_effect <- data.frame(Dimension = as.factor(levels(forage_formal$Dim_factor)),
                                   Effect = c(0,fixef(a_fit)[19:20,1]),
                                   Upper = c(0,fixef(a_fit, probs = c(0.05, 0.95))[19:20,4]),
                                   Lower = c(0,fixef(a_fit, probs = c(0.05, 0.95))[19:20,3]))

a_DimensionEffect <- ggplot(data = a_fit_dimension_effect, aes(x = Dimension, y = Effect)) + geom_point(size = 3) +
  theme_cowplot() + geom_hline(yintercept = 0, color = 'gray') + xlab('Dimension') + ylab('Partial Effect on\nSpace Clearance Rate\n(log scale)') + geom_errorbar(aes(ymin = Lower, ymax = Upper)) + 
  theme(axis.title.y = element_blank())

###############################################################################################
### Handling times 
###############################################################################################

residuals_h_fit <- as.data.frame(residuals(h_fit))

residuals_h_fit <- residuals_h_fit %>% mutate(partial_predmass = Estimate + -0.24*log(forage_arena$PredMass_g),
                                              partial_preymass = Estimate + 0.14*log(forage_arena$Mass_g),
                                              partial_temp = Estimate - 0.1*forage_arena$Temp_C + 0.0011*forage_arena$Temp_C^2,
                                              partial_arena = Estimate + 0.096*log(forage_arena$Arena_size),
                                              predmass = forage_arena$PredMass_g,
                                              preymass = forage_arena$Mass_g,
                                              temp = forage_arena$Temp_C,
                                              arena_size = forage_arena$Arena_size,
                                              prey_major_taxa = forage_arena$Prey_Major_Group_1,
                                              pred_major_taxa = forage_arena$Pred_Major_Grouping_1) 

### make plots 

### quantitative variables

# effect of predator mass

h_predmass_plot <- ggplot(data = residuals_h_fit, aes(x = log(predmass), y = partial_predmass)) + geom_point() + 
  theme_cowplot() + geom_function(fun = function(x) -0.24*x, color = 'black', size = 1) +
  xlab('log Predator Mass (g)') + ylab('Partial Residual of\nHandling Time') +   theme(axis.title.y = element_blank())

# effect of prey mass

h_preymass_plot <- ggplot(data = residuals_h_fit, aes(x = log(preymass), y = partial_preymass)) + geom_point() + 
  theme_cowplot() + geom_function(fun = function(x) 0.14*x, color = 'black', size = 1) +
  xlab('log Prey Mass (g)') + ylab('Partial Residual of\nHandling Time') +   theme(axis.title.y = element_blank())

# effect of temperature

h_temp_plot <- ggplot(data = residuals_h_fit, aes(x = temp, y = partial_temp)) + geom_point() + 
  theme_cowplot() + geom_function(fun = function(x) -0.1*x + 0.0011*x^2, color = 'black', size = 1) +
  xlab(expression(Temperature~(degree*C))) + ylab('Partial Residual of\nHandling Time') +   theme(axis.title.y = element_blank())

# effect of arena size

h_arena_size_plot <- ggplot(data = residuals_h_fit, aes(x = log(arena_size), y = partial_arena)) + geom_point() + 
  theme_cowplot() + geom_function(fun = function(x) 0.096*x, color = 'black', size = 1) +
  xlab(expression(log~Arena~Size~(cm^2~or~cm^3))) + ylab('Partial Residual of\nHandling Time') +   theme(axis.title.y = element_blank())

### discrete variables

# prey taxa effects

h_fit_prey_effect <- data.frame(Taxa = as.factor(levels(forage_arena$Prey_Major_Group_1)),
                                Effect = c(0, fixef(h_fit)[2:9,1]),
                                Upper = c(0, fixef(h_fit, probs = c(0.05, 0.95))[2:9,4]),
                                Lower = c(0, fixef(h_fit, probs = c(0.05, 0.95))[2:9,3]))

h_fit_prey_effect$Taxa <- factor(h_fit_prey_effect$Taxa, levels = c('Algae', 'Amphibian',
                                                                    'Arachnid', 'Crustacean',
                                                                    'Fish', 'Insect',
                                                                    'Mollusk', 'Protozoan',
                                                                    'Rotifer')) 

h_PreyTaxaEffect <- ggplot(data = h_fit_prey_effect, aes(x = Taxa, y = Effect, color = Taxa)) + geom_point() + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), color = 'black') + geom_point(size = 3) + theme_cowplot() + 
  scale_color_brewer(palette = 'Paired') + 
  geom_hline(yintercept = 0, color = 'gray')+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') +
  xlab('Prey Taxa') + ylab(' Partial Effect on\nHandling Time (log scale)') +   theme(axis.title.y = element_blank())

# pred taxa effects

h_fit_pred_effect <- data.frame(Taxa = as.factor(levels(forage_arena$Pred_Major_Grouping_1)),
                                Effect = c(0,fixef(h_fit)[10:16,1]),
                                Upper = c(0,fixef(h_fit, probs = c(0.05, 0.95))[10:16,4]),
                                Lower = c(0,fixef(h_fit, probs = c(0.05, 0.95))[10:16,3]))

h_PredTaxaEffect <- ggplot(data = h_fit_pred_effect, aes(x = Taxa, y = Effect, color = Taxa)) + geom_point() + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), color = 'black') + geom_point(size = 3) + theme_cowplot() + 
  scale_color_brewer(palette = 'Paired') + 
  geom_hline(yintercept = 0, color = 'gray') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = 'none') +
  xlab('Predator Taxa') + ylab(' Partial Effect on\nHandling Time (log scale)') +   theme(axis.title.y = element_blank())

# habitat effect

h_fit_habitat_effect <- data.frame(Habitat = as.factor(levels(forage_arena$Habitat_Aquatic)),
                                   Effect = c(0,fixef(h_fit)[17:18,1]),
                                   Upper = c(0,fixef(h_fit, probs = c(0.05, 0.95))[17:18,4]),
                                   Lower = c(0,fixef(h_fit, probs = c(0.05, 0.95))[17:18,3]))

h_HabitatEffect <- ggplot(data = h_fit_habitat_effect, aes(x = Habitat, y = Effect)) + geom_point(size = 3) +
  theme_cowplot() + geom_hline(yintercept = 0, color = 'gray') + xlab('Habitat') + ylab('Partial Effect on Handling Time\n(log scale)') + geom_errorbar(aes(ymin = Lower, ymax = Upper)) +
  theme(axis.title.y = element_blank())

# dimension effect

h_fit_dimension_effect <- data.frame(Dimension = as.factor(levels(forage_formal$Dim_factor)),
                                     Effect = c(0,fixef(h_fit)[19:20,1]),
                                     Upper = c(0,fixef(h_fit, probs = c(0.05, 0.95))[19:20,4]),
                                     Lower = c(0,fixef(h_fit, probs = c(0.05, 0.95))[19:20,3]))

h_DimensionEffect <- ggplot(data = h_fit_dimension_effect, aes(x = Dimension, y = Effect)) + geom_point(size = 3) +
  theme_cowplot() + geom_hline(yintercept = 0, color = 'gray') + xlab('Dimension') + ylab('Partial Effect on Handling Time\n(log scale)') + geom_errorbar(aes(ymin = Lower, ymax = Upper)) +
  theme(axis.title.y = element_blank())




### add create plots for manuscript
### make plots similar to the plots for saturation index but for space clearance 
### rates and handling times 

### space clearance rate plots

a_QuantPlots <- plot_grid(a_preymass_plot, a_predmass_plot, a_temp_plot, a_arena_size_plot,
                          ncol = 1, nrow = 4, labels = c('E', 'F', 'G', 'H'), label_x = 0.95)

a_Quant_ylabel <- textGrob('Partial Residual of\nlog Space Clearance Rate', gp = gpar(fontsize = 15) ,rot = 90)

a_Quant <- grid.arrange(arrangeGrob(a_QuantPlots, left = a_Quant_ylabel))

a_QualPlots <- plot_grid(a_PreyTaxaEffect, a_PredTaxaEffect, a_HabitatEffect, a_DimensionEffect, ncol = 1, nrow = 4,
                         labels = c('A', 'B', 'C', 'D'), label_x = 0.95)

a_Qual_ylabel <- textGrob('Partial Effect (log scale)', gp = gpar(fontsize = 15) ,rot = 90)

a_Qual <- grid.arrange(arrangeGrob(a_QualPlots, left = a_Qual_ylabel))

# put plots together

a_plots <- plot_grid(a_Qual, a_Quant, nrow = 1, ncol = 2)

# save plot

save_plot(filename = 'SpaceClearanceRatePlot.png', plot = a_plots, ncol = 2, nrow = 2.85)


### handling time plots

h_QuantPlots <- plot_grid(h_preymass_plot, h_predmass_plot, h_temp_plot, h_arena_size_plot,
                          ncol = 1, nrow = 4, labels = c('E', 'F', 'G', 'H'), label_x = 0.95)

h_Quant_ylabel <- textGrob('Partial Residual of\nlog Handling Time', gp = gpar(fontsize = 15) ,rot = 90)

h_Quant <- grid.arrange(arrangeGrob(h_QuantPlots, left = h_Quant_ylabel))

h_QualPlots <- plot_grid(h_PreyTaxaEffect, h_PredTaxaEffect, h_HabitatEffect, h_DimensionEffect, ncol = 1, nrow = 4,
                         labels = c('A', 'B', 'C', 'D'), label_x = 0.95)

h_Qual_ylabel <- textGrob('Partial Effect (log scale)', gp = gpar(fontsize = 15) ,rot = 90)

h_Qual <- grid.arrange(arrangeGrob(h_QualPlots, left = h_Qual_ylabel))

# put plots together

h_plots <- plot_grid(h_Qual, h_Quant, nrow = 1, ncol = 2)

save_plot(filename = 'HandlingTimePlot.png', plot = h_plots, ncol = 2, nrow = 2.85)





