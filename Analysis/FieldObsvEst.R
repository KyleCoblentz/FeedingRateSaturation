################################################################################
### Supplemental Analysis -- Comparing estimated and observed saturation 
### index values for field studies
################################################################################

### load libraries

library(tidyverse); library(cowplot);

### load data

# saturation data

forage <- read.csv('forage_fracfeed.csv')

# original data 

orig_curves <- read.csv('FoRAGE_db_V2_3_Oct_26_2022_original_curves.csv')

### Want to compare the estimated saturation index values for the field studies to 
### the observed values for each prey density in the study 

### get just the field studies from the dataset

forage_field <- forage %>% filter(Field == 'Y') 

### now get just the field studies from the original curves

orig_curves_field <- orig_curves %>% filter(X %in% unique(forage_field$Data.set))

### make single prey density field

orig_curves_field$PreyDensity <- as.numeric(ifelse(is.na(orig_curves_field$X2D.Prey.density) == FALSE,
                                        orig_curves_field$X2D.Prey.density, 
                                        orig_curves_field$X3D.Prey.density))

### from both data sets drop everything that we don't need 

### for forage, we need estimates of a and h, the study, and the saturation index values at the 
### 10th, 50th, and 90th percentiles of estimated prey densities

forage_field <- forage_field %>% select(Data.set, Source, Fitted_a, Fitted_h_day, frac_10, frac_50, frac_90,
                                        PreyAbundance_10, PreyAbundance_90)

### for original curves, need prey density and the study

colnames(orig_curves_field)[1] = 'Data.set'

orig_curves_field <- orig_curves_field %>% select(Data.set, Source, PreyDensity)

### left join these two data sets into a single data set

combined_data <- left_join(orig_curves_field, forage_field, by = 'Data.set')

### now create a column that is the observed amount of saturation in the field for that data point

combined_data <- combined_data %>% mutate(ObsSat = Fitted_a*Fitted_h_day*PreyDensity/(1 + Fitted_a*Fitted_h_day*PreyDensity))

combined_data$Study_Dataset <- paste(combined_data$Source.x, combined_data$Data.set, sep = '_')

combined_data$in_limit <- combined_data$ObsSat >= combined_data$frac_10 & combined_data$ObsSat <= combined_data$frac_90

combined_data_group <- combined_data %>% group_by(Study_Dataset) %>% summarize(min_sat = min(frac_10), max_sat = max(frac_90))

ObsvEstFieldPlot <- ggplot(data = combined_data, aes(x = Study_Dataset, y = ObsSat, color = in_limit)) + geom_point() + 
  geom_segment(data = combined_data_group, aes(x = Study_Dataset, xend = Study_Dataset, y = min_sat, yend = max_sat), inherit.aes = FALSE) + theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(breaks = c(TRUE, FALSE), values = c('black', 'red'), name = 'In Interval?') + 
  xlab('Data Source') + ylab('Saturation Index Value')

save_plot(filename = 'ObsvEstFieldPlot.png', plot = ObsvEstFieldPlot, base_height = 8.25, base_width = 8.75)

