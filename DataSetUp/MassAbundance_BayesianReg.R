#####################################################################################
### Bayesian reanalysis of Hatton et al. 2019. PNAS
#####################################################################################

### load packages

library(dplyr); library(ggplot2); library(cowplot); library(brms);

### load data

scale_data <- read.csv('MassAbundanceScalingData.csv')

### want to estimate the mass-abundance models done by Hatton et al. (but only 
### the ones that we will need for FoRAGE)

### will just be log-log regressions

### Mammals

mammal <- scale_data %>% filter(Major_taxa == 'Mammal')

mammal_fit <- brm(formula = log(Density_Nm2) ~ log(Mass_g), data = mammal)

### summary

summary(mammal_fit, prob = 0.9)

bayes_R2(mammal_fit)

### plots

# generate prediction interval dataset

mammal_pred_interval <- data.frame(logMass_g = seq(min(log(mammal$Mass_g)), max(log(mammal$Mass_g)), length.out = 1000))

mammal_pred_interval$Mass_g <- exp(mammal_pred_interval$logMass_g)

mammal_predict <- posterior_predict(mammal_fit, newdata = mammal_pred_interval)

mammal_pred_interval <- mammal_pred_interval %>% mutate(Median = apply(mammal_predict, 2, median), 
                                                        Upper = apply(mammal_predict, 2, function(x) quantile(x, probs = c(0.90))),
                                                        Lower = apply(mammal_predict, 2, function(x) quantile(x, probs = c(0.10))))
 
Mammal_plot <- ggplot(data = mammal, aes(x = log(Mass_g), y = log(Density_Nm2))) + geom_point() + 
  ggtitle('Mammals') + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line(data = mammal_pred_interval, aes(x = log(Mass_g), y = Median), color = 'red') +
  geom_ribbon(data = mammal_pred_interval, aes(ymin = Lower, ymax = Upper, x = log(Mass_g)), inherit.aes = FALSE, alpha = 0.5) + 
  xlab('log Mass (g)') + ylab('log Density (Number m^-2)')

### Protists

protist <- scale_data %>% filter(Major_taxa == 'Protist')

protist_fit <- brm(formula = log(Density_Nm2) ~ log(Mass_g), data = protist)

### summary

summary(protist_fit, prob = 0.9)

bayes_R2(protist_fit)

### plots 

protist_pred_interval <- data.frame(logMass_g = seq(min(log(protist$Mass_g)), max(log(protist$Mass_g)), length.out = 1000))

protist_pred_interval$Mass_g <- exp(protist_pred_interval$logMass_g)

protist_predict <- posterior_predict(protist_fit, newdata = protist_pred_interval)

protist_pred_interval <- protist_pred_interval %>% mutate(Median = apply(protist_predict, 2, median), 
                                                        Upper = apply(protist_predict, 2, function(x) quantile(x, probs = c(0.90))),
                                                        Lower = apply(protist_predict, 2, function(x) quantile(x, probs = c(0.10))))

Protist_plot <- ggplot(data = protist, aes(x = log(Mass_g), y = log(Density_Nm2))) + geom_point() + 
  ggtitle('Protists') + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line(data = protist_pred_interval, aes(x = log(Mass_g), y = Median), color = 'red') +
  geom_ribbon(data = protist_pred_interval, aes(ymin = Lower, ymax = Upper, x = log(Mass_g)), inherit.aes = FALSE, alpha = 0.5) + 
  xlab('log Mass (g)') + ylab('log Density (Number m^-2')

### Ectotherms

ectotherm <- scale_data %>% filter(Major_taxa == 'Invertebrate'| Major_taxa == 'Reptilia' | Major_taxa == 'Amphibia' | Major_taxa == 'Fish')

ectotherm_fit <- brm(formula = log(Density_Nm2) ~ log(Mass_g), data = ectotherm)

### summary

summary(ectotherm_fit, prob = 0.9)

bayes_R2(ectotherm_fit)

### plot

ectotherm_pred_interval <- data.frame(logMass_g = seq(min(log(ectotherm$Mass_g)), max(log(ectotherm$Mass_g)), length.out = 1000))

ectotherm_pred_interval$Mass_g <- exp(ectotherm_pred_interval$logMass_g)

ectotherm_predict <- posterior_predict(ectotherm_fit, newdata = ectotherm_pred_interval)

ectotherm_pred_interval <- ectotherm_pred_interval %>% mutate(Median = apply(ectotherm_predict, 2, median), 
                                                          Upper = apply(ectotherm_predict, 2, function(x) quantile(x, probs = c(0.90))),
                                                          Lower = apply(ectotherm_predict, 2, function(x) quantile(x, probs = c(0.10))))

Ectotherm_plot <- ggplot(data = ectotherm, aes(x = log(Mass_g), y = log(Density_Nm2))) + geom_point() + 
  ggtitle('Ectotherms') + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line(data = ectotherm_pred_interval, aes(x = log(Mass_g), y = Median), color = 'red') +
  geom_ribbon(data = ectotherm_pred_interval, aes(ymin = Lower, ymax = Upper, x = log(Mass_g)), inherit.aes = FALSE, alpha = 0.5) + 
  xlab('log Mass (g)') + ylab('log Density (Number m^-2)')


### birds

bird <- scale_data %>% filter(Major_taxa == 'Bird')   

bird_fit <- brm(formula = log(Density_Nm2) ~ log(Mass_g), data = bird)

summary(bird_fit, prob = 0.9)

bayes_R2(bird_fit)

### plot

bird_pred_interval <- data.frame(logMass_g = seq(min(log(bird$Mass_g)), max(log(bird$Mass_g)), length.out = 1000))

bird_pred_interval$Mass_g <- exp(bird_pred_interval$logMass_g)

bird_predict <- posterior_predict(bird_fit, newdata = bird_pred_interval)

bird_pred_interval <- bird_pred_interval %>% mutate(Median = apply(bird_predict, 2, median), 
                                                              Upper = apply(bird_predict, 2, function(x) quantile(x, probs = c(0.90))),
                                                              Lower = apply(bird_predict, 2, function(x) quantile(x, probs = c(0.10))))

Bird_plot <- ggplot(data = bird, aes(x = log(Mass_g), y = log(Density_Nm2))) + geom_point() + 
  ggtitle('Birds') + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line(data = bird_pred_interval, aes(x = log(Mass_g), y = Median), color = 'red') +
  geom_ribbon(data = bird_pred_interval, aes(ymin = Lower, ymax = Upper, x = log(Mass_g)), inherit.aes = FALSE, alpha = 0.5) + 
  xlab('log Mass (g)') + ylab('log Density (Number m^-2)')

### Note: Hatton et al. give some skepticism about the bird density data because 
### estimating bird densities is apparently really hard.

### bacteria -- note that this is really mostly marine phytoplankton <20um from one study
### also note that this is the aerial density

prokaryote <- scale_data %>% filter(Major_taxa == 'Prokaryote')

prokaryote_fit <- brm(formula = log(Density_Nm2) ~ log(Mass_g), data = prokaryote)

summary(prokaryote_fit, prob = 0.9)

bayes_R2(prokaryote_fit)

#### plots

prokaryote_pred_interval <- data.frame(logMass_g = seq(min(log(prokaryote$Mass_g)), max(log(prokaryote$Mass_g)), length.out = 1000))

prokaryote_pred_interval$Mass_g <- exp(prokaryote_pred_interval$logMass_g)

prokaryote_predict <- posterior_predict(prokaryote_fit, newdata = prokaryote_pred_interval)

prokaryote_pred_interval <- prokaryote_pred_interval %>% mutate(Median = apply(prokaryote_predict, 2, median), 
                                                    Upper = apply(prokaryote_predict, 2, function(x) quantile(x, probs = c(0.90))),
                                                    Lower = apply(prokaryote_predict, 2, function(x) quantile(x, probs = c(0.10))))

Prokaryote_plot <- ggplot(data = prokaryote, aes(x = log(Mass_g), y = log(Density_Nm2))) + geom_point() +
  ggtitle('Prokaryotes/Algae') + theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_line(data = prokaryote_pred_interval, aes(x = log(Mass_g), y = Median), color = 'red') +
  geom_ribbon(data = prokaryote_pred_interval, aes(ymin = Lower, ymax = Upper, x = log(Mass_g)), inherit.aes = FALSE, alpha = 0.5) + 
  xlab('log Mass (g)') + ylab('log Density (Number m^-2)')

### All taxa plot

### change major_taxa to ectotherm for some stuff

scale_data$Major_taxa[which(scale_data$Major_taxa == 'Invertebrate' | scale_data$Major_taxa == 'Reptilia' | scale_data$Major_taxa == 'Amphibia' | scale_data$Major_taxa == 'Fish')] <- 'Ectotherm'

### change prokaryote to prokaryote/algae

scale_data$Major_taxa[which(scale_data$Major_taxa == 'Prokaryote')] <- 'Prokaryote/Algae'

### drop plants

scale_data_mplants <- scale_data %>% filter(Major_taxa != 'Plant')

### Make plot

All_taxa_plot <- ggplot(data = scale_data_mplants, aes(x = log(Mass_g), y = log(Density_Nm2), color = Major_taxa)) + geom_point() + 
  geom_smooth(method = 'lm', se = FALSE) + theme_cowplot() + xlab('log Mass (g)') + ylab('log Density (Number m^-2)') +
  scale_color_discrete(name = 'Major Taxa')

### put plots together

Mass_Abundance_plot <- plot_grid(Mammal_plot, Bird_plot, Ectotherm_plot, Protist_plot, Prokaryote_plot, All_taxa_plot, nrow = 3,
                                 ncol = 2, labels = 'AUTO')

save_plot(filename = 'Mass_Abundance_plot.png', plot = Mass_Abundance_plot, ncol = 2, nrow = 3)






