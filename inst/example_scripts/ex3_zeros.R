# Exploring alternative assumptions for zero biomass observations.
# Example: other rockfish (all species except shortspine thornyhead) in the Bering Sea

# Example R scripts and data files:
# inst/example_scripts
# inst/example_data

# set up ----
library(rema)
library(dplyr)
library(cowplot) # provides helpful plotting utilities

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# create directory for analysis
# out_path <- "test_rema"
if(!exists("out_path")) out_path = getwd()
if(!dir.exists(out_path)) dir.create(out_path)

# copy all data files to working directory
rema_path <- find.package('rema')
example_data_files <- list.files(path = file.path(rema_path, "example_data"))
example_data_files
file.copy(from = file.path(path = file.path(rema_path, "example_data"),
                           example_data_files),
          to = file.path(file.path(out_path), example_data_files),
          overwrite = TRUE)

setwd(out_path)

