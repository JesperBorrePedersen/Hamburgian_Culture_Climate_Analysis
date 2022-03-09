# Hamburgian Culture Climate Analysis - First data Exploraiton 
# Jakob Assmann j.assmann@bios.au.dk 17 August 2020

## 1) Housekeeping ----
# |_ Global matters ----
# Dependencies
library(sf)
library(raster)
library(tidyverse)
library(cowplot)
library(rasterVis)
library(colorspace)
# library(rnaturalearth)
# library(rnaturalearthdata)
library(dismo)
library(magick)

# Raster options
rasterOptions(progress = "text")

# Set sf without s2 (the original analysis was developed before the move to s2 in sf)
sf_use_s2(FALSE)

# Set area of interest
area_of_interest <- extent(c(-11,32,49,60)) 

# |_ Set colour + shape ----
# Colour
col_pulse1 <- "black"
col_pulse1_transparent <- "#1C87E6FF"
col_pulse2 <- "black"
col_pulse2_transparent <- "white" # "#9317FCFF"
col_grey_light <- "#D3D3D3FF"
col_grey_light_transparent <- "#D3D3D344"
col_grey_dark <- "#696969FF"
col_grey_dark_transparent <- "#69696944"

# shape
shape_pulse_1 <- 21
shape_pulse_2 <- 21
shape_uncertain <- 63

# |_ Load data ----
# Load in CHELSA LGM files
temp <- list.files(
  "O:/Nat_Ecoinformatics/B_Read/LegacyData/World/Climate/PaleoClimate/Original/CHELSA_TRACE/BIO_01/",
  pattern = "*.tif",
  full.names = T)
precip <- list.files(
  "O:/Nat_Ecoinformatics/B_Read/LegacyData/World/Climate/PaleoClimate/Original/CHELSA_TRACE/BIO_12/",
  pattern = "*.tif",
  full.names = T)

# Load rasters as stack
temp <- stack(temp)
precip <- stack(precip)

# Rename layers to something a bit more sensible
names_temp <- gsub(".*bio01_(.*)_.*", "temp_\\1", names(temp))
names_temp <- gsub("temp_\\.(.*)", "temp_\\1_BCE", names_temp)
names_temp <- gsub("temp_([0-9]*)$", "temp_\\1_CE", names_temp)
names(temp) <- names_temp

names_precip <- gsub(".*bio12_(.*)_.*", "precip_\\1", names(precip))
names_precip <- gsub("precip_\\.(.*)", "precip_\\1_BCE", names_precip)
names_precip <- gsub("precip_([0-9]*)$", "precip_\\1_CE", names_precip)
names(precip) <- names_precip

# Load mean rasters
hamburgian_mean_temp <- raster("data/mean_rasters/hamburgian_mean_temp.tif")
hamburgian_mean_temp_pulse_1 <- raster("data/mean_rasters/hamburgian_mean_temp_pulse_1.tif")
hamburgian_mean_temp_pulse_2 <- raster("data/mean_rasters/hamburgian_mean_temp_pulse_2.tif")

hamburgian_mean_precip <- raster("data/mean_rasters/hamburgian_mean_precip.tif")
hamburgian_mean_precip_pulse_1 <- raster("data/mean_rasters/hamburgian_mean_precip_pulse_1.tif")
hamburgian_mean_precip_pulse_2 <- raster("data/mean_rasters/hamburgian_mean_precip_pulse_2.tif")

# Generate land from Bølling/Meiendorf Map
land <- read_sf("data/shp/EPHA_Mei_v110.shp") %>%
  filter(typ == "land area") %>%
  st_union()
ocean <- read_sf("data/shp/EPHA_Mei_v110.shp") %>%
  filter(typ == "ocean" | typ == "lake") %>%
  st_union()
ice <- read_sf("data/shp/EPHA_Mei_v110.shp") %>%
  filter(typ == "ice sheet" | typ == "ice lake") %>%
  st_union()

land_for_maps <- land %>%
  st_transform(st_crs(hamburgian_mean_temp)) %>%
  as_Spatial() %>%
  crop(hamburgian_mean_temp)

ocean_for_maps <- ocean %>%
  st_transform(st_crs(hamburgian_mean_temp)) %>%
  as_Spatial() %>%
  crop(hamburgian_mean_temp)

ice_for_maps <- ice %>%
  st_transform(st_crs(hamburgian_mean_temp)) %>%
  as_Spatial() %>%
  crop(hamburgian_mean_temp)

# Create look up tables for years using BP 
temp_look_up <- data.frame(
  names_temp = names_temp,
  year_CE = as.numeric(gsub(".*_([0-9]*)_.*", "\\1", names_temp)) * 100,
  epoch = gsub(".*_([BCE]*)$", "\\1", names_temp),
  year_BP = NA)

temp_look_up$year_BP <- sapply(
  names_temp, 
  function(x){
    year_BP <- temp_look_up[temp_look_up$names_temp == x,]$year_CE
    if (temp_look_up[temp_look_up$names_temp == x,]$epoch == "BCE"){
      year_BP <- year_BP + 2000
    } else {
      year_BP <- 2000 - year_BP 
    }
    return(year_BP)
  }
)

temp_look_up <- arrange(temp_look_up, year_BP)

precip_look_up <- data.frame(
  names_precip = names_precip,
  year_CE = as.numeric(gsub(".*_([0-9]*)_.*", "\\1", names_precip)) * 100,
  epoch = gsub(".*_([BCE]*)$", "\\1", names_precip),
  year_BP = NA)

precip_look_up$year_BP <- sapply(
  names_precip, 
  function(x){
    year_BP <- precip_look_up[precip_look_up$names_precip == x,]$year_CE
    if (precip_look_up[precip_look_up$names_precip == x,]$epoch == "BCE"){
      year_BP <- year_BP + 2000
    } else {
      year_BP <- 2000 - year_BP 
    }
    return(year_BP)
  }
)

precip_look_up <- arrange(precip_look_up, year_BP)

# Load site locations 
hamburgian_sites <- read.csv("data/ham_data.csv") %>%
  mutate(pulse_1_start = 14500,
         pulse_1_end = 14300,
         pulse_2_start = 14200,
         pulse_2_end = 14100)


# Covnert to sf object
hamburgian_sites <- st_as_sf(hamburgian_sites, 
                             coords = c("long", "lat"),
                             crs = 4326,
                             remove = F)

# Load climate data extractions 
load("data/temp_df.Rda")
load("data/precip_df.Rda")

# Transform hamburgian sites
hamburgian_sites <- st_transform(hamburgian_sites, st_crs(hamburgian_mean_temp))
hamburgian_sites_sp <- as_Spatial(hamburgian_sites)

## 2) Prep Climate Variables ----

# Mean temp / precip in the two Hamburgian Culture periods
mean_temp <- temp_df %>% 
  group_by(site) %>%
  filter(year_BP >= 14100 & year_BP < 14500) %>%
  summarise(mean_temp = mean(temp))
mean_temp_pulse_1 <- temp_df %>% 
  group_by(site) %>%
  filter(year_BP >= 14300 & year_BP < 14500) %>%
  summarise(mean_temp_pulse_1 = mean(temp))
mean_temp_pulse_2 <- temp_df %>% 
  group_by(site) %>%
  filter(year_BP >= 14100 & year_BP < 14200) %>%
  summarise(mean_temp_pulse_2 = mean(temp))

mean_precip <- precip_df %>% 
  group_by(site) %>%
  filter(year_BP >= 14100 & year_BP < 14500) %>%
  summarise(mean_precip = mean(precip))
mean_precip_pulse_1 <- precip_df %>% 
  group_by(site) %>%
  filter(year_BP >= 14300 & year_BP < 14500) %>%
  summarise(mean_precip_pulse_1 = mean(precip))
mean_precip_pulse_2 <- precip_df %>% 
  group_by(site) %>%
  filter(year_BP >= 14100 & year_BP < 14200) %>%
  summarise(mean_precip_pulse_2 = mean(precip))

# merge with hamburgian_sites df
hamburgian_sites <- hamburgian_sites %>%
  full_join(mean_temp) %>%
  full_join(mean_temp_pulse_1) %>%
  full_join(mean_temp_pulse_2) %>%
  full_join(mean_precip) %>%
  full_join(mean_precip_pulse_1) %>%
  full_join(mean_precip_pulse_2) 

# tidy up
rm(mean_temp)
rm(mean_temp_pulse_1)
rm(mean_temp_pulse_2)
rm(mean_precip)
rm(mean_precip_pulse_1)
rm(mean_precip_pulse_2)

# set mean temp to NA where no evidence is available
hamburgian_sites$mean_temp_pulse_1[hamburgian_sites$chron_association != "pulse_1"] <- NA
hamburgian_sites$mean_temp_pulse_2[hamburgian_sites$chron_association != "pulse_2"] <- NA
hamburgian_sites$mean_precip_pulse_1[hamburgian_sites$chron_association != "pulse_1"] <- NA
hamburgian_sites$mean_precip_pulse_2[hamburgian_sites$chron_association != "pulse_2"] <- NA

hamburgian_sites$pulse_1_start[hamburgian_sites$chron_association != "pulse_1"] <- NA
hamburgian_sites$pulse_1_end[hamburgian_sites$chron_association != "pulse_1"] <- NA
hamburgian_sites$pulse_2_start[hamburgian_sites$chron_association != "pulse_2"] <- NA
hamburgian_sites$pulse_2_end[hamburgian_sites$chron_association != "pulse_2"] <- NA

# Remove "creswellian" sites
hamburgian_sites <- hamburgian_sites %>% filter(arch_association != "creswellian")

# Order sites by latitude
hamburgian_sites <- hamburgian_sites %>%
  mutate(site = ordered(site, levels = levels(fct_reorder(hamburgian_sites$site, hamburgian_sites$lat))))

# 3) Fit glm Models ----

# Add cell_id to main geometry to allow for exclusion of multiple samples from
# the a single cell
hamburgian_sites <- raster::extract(hamburgian_mean_temp, 
                as_Spatial(hamburgian_sites),
                cellnumbers = T,
                df = T) %>% 
  mutate(site = hamburgian_sites$site,
         cell_id = cells) %>%
  select(site, cell_id) %>%
  full_join(hamburgian_sites, .)

# Gather temperature and precip columns into training data
bioclim_training_all <- hamburgian_sites %>%
  st_drop_geometry() %>%
  filter(chron_association != "uncertain") %>%
  select(site, 
         lat,
         long,
         cell_id,
         mean_temp_pulse_1,
         mean_temp_pulse_2) %>%
  gather("temp_var", "mean_temp", mean_temp_pulse_1, mean_temp_pulse_2) %>%
  na.omit() %>%
  distinct(cell_id, temp_var, mean_temp) %>%
  mutate(period = gsub("mean_temp_(.*)", "\\1", temp_var))
bioclim_training_all <- hamburgian_sites %>%
  st_drop_geometry() %>%
  filter(chron_association != "uncertain") %>%
  select(site, 
         lat,
         long,
         cell_id,
         mean_precip_pulse_1,
         mean_precip_pulse_2) %>%
  gather("precip_var", "mean_precip", mean_precip_pulse_1, mean_precip_pulse_2) %>%
  na.omit() %>%
  distinct(cell_id, precip_var, mean_precip) %>%
  mutate(period = gsub("mean_precip_(.*)", "\\1", precip_var)) %>%
  full_join(bioclim_training_all, .)
bioclim_training_all <- hamburgian_sites %>% filter(!(cell_id %in% bioclim_training_all$cell_id)) %>% 
  st_drop_geometry() %>%
  select(site,
         lat,
         long,
         cell_id,
         mean_temp,
         mean_precip,
         chron_association) %>%
  mutate(temp_var = "mean_temp", 
         precip_var = "mean_precip",
         period = "uncertain") %>%
  distinct(cell_id, temp_var, mean_temp, precip_var, mean_precip, period) %>%
  full_join(bioclim_training_all)

# Quality control
view(bioclim_training_all)
# bioclim_training_all <- filter(bioclim_training_all, period == "pulse_2")

# Retain climate columns only
bioclim_training <- select(bioclim_training_all, "mean_temp", "mean_precip")

# Prepare raster stacks for predictions and background data
climate_mean <-  stack(hamburgian_mean_temp,
                       hamburgian_mean_precip)
climate_mean <- setNames(climate_mean, 
                            c("mean_temp", "mean_precip"))
climate_pulse_1 <- stack(hamburgian_mean_temp_pulse_1,
                             hamburgian_mean_precip_pulse_1)
climate_pulse_1 <- setNames(climate_pulse_1, 
                                c("mean_temp", "mean_precip"))
climate_pulse_2 <- stack(hamburgian_mean_temp_pulse_2,
                             hamburgian_mean_precip_pulse_2)
climate_pulse_2 <- setNames(climate_pulse_2,
                                c("mean_temp", "mean_precip"))

# Mask out everything but land
climate_mean_masked <- mask(climate_mean,
                            land_for_maps)
climate_pulse_1_masked <- mask(climate_pulse_1,
                                   land_for_maps)
climate_pulse_2_masked <- mask(climate_pulse_2,
                                   land_for_maps)

# Generate random background points for evaluation (absence data)
set.seed(50)
absence_data <- randomPoints(climate_pulse_1_masked,
                             900)
colnames(absence_data) <- c('lon', 'lat')
absence_data <- as.data.frame(absence_data)
absence_data$period <- c(rep("uncertain", nrow(absence_data)/3),
                         rep("pulse_1", nrow(absence_data)/3),
                         rep("pulse_2", nrow(absence_data)/3))
# Add cell_ids
absence_data$cell_id <- st_as_sf(absence_data,
                                 coords = c("lon", "lat"),
                                 crs = st_crs(climate_pulse_1)) %>%
  raster::extract(climate_pulse_1,
                  .,
                  df = T,
                  cellnumbers = T) %>% 
  pull(cells)

# filter cell_ids to avoid coccurence with presence data
absence_data <- absence_data[!(absence_data$cell_id %in% bioclim_training_all$cell_id),]

# fill in absence data, 1/3 from each climate_data
absence_data$mean_temp <- c(raster::extract(climate_mean_masked[["mean_temp"]], 
                                            absence_data$cell_id[absence_data$period == "uncertain"]),
                            raster::extract(climate_pulse_1_masked[["mean_temp"]], 
                                            absence_data$cell_id[absence_data$period == "pulse_1"]),
                            raster::extract(climate_pulse_2_masked[["mean_temp"]], 
                                            absence_data$cell_id[absence_data$period == "pulse_2"]))
absence_data$mean_precip <- c(raster::extract(climate_mean_masked[["mean_precip"]], 
                                            absence_data$cell_id[absence_data$period == "uncertain"]),
                            raster::extract(climate_pulse_1_masked[["mean_precip"]], 
                                            absence_data$cell_id[absence_data$period == "pulse_1"]),
                            raster::extract(climate_pulse_2_masked[["mean_precip"]], 
                                            absence_data$cell_id[absence_data$period == "pulse_2"]))
## Combine into one dataframe
absence_data$presence <- 0
bioclim_training_all$presence <- 1

# Define function to carry out model fitting 
fit_glm <- function(training_data, absence_data, k){
    # Status
    cat("Fitting glm model with", k, "-fold cross validation.\n" )
  
    # Add grouping variable
    training_data$group <- kfold(nrow(training_data), k)
    absence_data$group <- kfold(nrow(absence_data), k)
    
    # Set glm formula
    glm_formula <- formula(presence ~ mean_temp + mean_precip + 
                             I(mean_temp^2) + I(mean_precip^2) + 
                             mean_temp:mean_precip)
    # Prepare output lists
    glm_models <- list()
    glm_evals <- list()
    
    # Carry out cross validation
    for(i in 1:5){
      # Status
      cat("CV round", i, "\n")
      
      # Subset training data
      training_group <- filter(training_data, group != i) %>%
        bind_rows(filter(absence_data, group != i))

      # Fit Bioclim envlope model to training data
      glm_models[[i]] <- glm(glm_formula,
                             data = training_group,
                             family = "binomial")
            # Evaluate model
      glm_evals[[i]] <- evaluate(p = training_group %>% 
                                    filter(presence  == 1),
                                  a = training_group %>%  
                                    filter(presence == 0), 
                                   model = glm_models[[i]],
                                 type = "response")
    }
    
    # Status 
    cat("Calculating summary stats...\n")
    
    # Calc mean auc
    mean_auc <- mean(unlist(lapply(glm_evals, 
                                   function(x) x@auc))) 
    
    # Determine optimum threshold 
    # (by maximizing specific sensitivity / TSS)
    # and extract mean kappa at this threshold for all models
    mean_kappa <- mean(unlist(lapply(glm_evals, 
                                     function(x) {
                                       thresh <- threshold(x)["spec_sens"][1,1]
                                       x@kappa[which(x@t == thresh)]
                                     })))
    
    # Determine optimum threshold 
    # (by maximizing specific sensitivity / TSS)
    # and extract mean TSS at this threshold for all models
    mean_tss <- mean(unlist(lapply(glm_evals, 
                                   function(x) {
                                     thresh <- threshold(x)["spec_sens"][1,1,]
                                     x@TPR[which(x@t == thresh)] + x@TNR[which(x@t == thresh)] - 1
                                   }))) 
    
    # Status
    cat("Fitting final model...\n")
    
    # fit single model 
    glm_model <- glm(glm_formula,
                     data = bind_rows(training_data,
                                      absence_data),
                     family = "binomial")
    glm_eval <- evaluate(p = bind_rows(training_data,
                                          absence_data) %>%
                              filter(presence == 1),
                            a = bind_rows(training_data,
                                          absence_data) %>%
                              filter(presence == 0),
                            glm_model,
                         type = "response")
    
    # Determine threshold with max specific sensitivity / TSS
    glm_thresh <-  threshold(glm_eval, 'spec_sens') 
    
    # Extract AUC, Kappa, TSS
    auc <- glm_eval@auc
    kappa <-  glm_eval@kappa[which(glm_eval@t == glm_thresh)]
    tss <-  glm_eval@TPR[which(glm_eval@t == glm_thresh)] + glm_eval@TNR[which(glm_eval@t == glm_thresh)] - 1
    
    # Status
    cat("Done.\n")
    
    # return all outputs
    list(glm_models = glm_models,
         glm_evals = glm_evals,
         mean_auc = mean_auc,
         mean_kappa = mean_kappa,
         mean_tss = mean_tss,
         glm_model = glm_model,
         glm_eval = glm_eval,
         glm_thresh = glm_thresh,
         auc = auc,
         kappa = kappa,
         tss = tss)
}


# Fit model for all three time_periods
set.seed(7)

# Global model
model_global <- fit_glm(bioclim_training_all,
                           absence_data,
                           5)
# Pulse 1
model_pulse_1 <- fit_glm(bioclim_training_all %>% 
                          filter(period == "pulse_1"),
                        absence_data %>% 
                          filter(period == "pulse_1"),
                        5)

# Pulse 2
model_pulse_2 <- fit_glm(bioclim_training_all %>% 
                          filter(period == "pulse_2"),
                          absence_data %>% 
                          filter(period == "pulse_2"),
                        5)

# Evaluation summary
model_list <- list(model_global, model_pulse_1, model_pulse_2)
eval_summary <- data.frame(
  model_name = c("model_global", "model_pulse_1", "model_pulse_2"),
  n_cross = rep(5,3),
  mean_auc = sapply(model_list, function(x) round(x$mean_auc,2)),
  mean_kappa = sapply(model_list, function(x) round(x$mean_kappa,2)),
  mean_tss = sapply(model_list, function(x) round(x$mean_tss,2)),
  auc = sapply(model_list, function(x) round(x$auc,2)),
  kappa = sapply(model_list, function(x) round(x$kappa,2)),
  tss = sapply(model_list, function(x) round(x$tss,2))) 

write_csv(eval_summary, "tables/supplementary/model_evaluation_glm.csv")

# Write out p-values for model parameters
mapply(function(x, y){
  x <- summary(x$glm_model)$coefficients %>%
    as.data.frame()
  x <- rownames_to_column(x, var = "predictor")
  x[,5] <- round(x[,5], 8)
  x[,2] <- round(x[,2], 5)
  select(x, c(1,2,5)) %>%
    setNames(c("Predictor", 
               paste0("Estimate (", y, ")"), 
               paste0("p-value (", y, ")"))) %>%
    return()
  },
  model_list,
  c("global", "pulse_1", "pulse_2"),
  SIMPLIFY = F) %>% 
  reduce(full_join) %>%
  mutate(Predictor = c("Intercept",
                       "Temperature",
                       "Precipitation",
                       "Temperature²",
                       "Precipitation²",
                       "Temperature :  Precipitation")) %>%
  write_csv("tables/supplementary/model_coefs_glm.csv")

# Predict for core time periods
predictions_global <- predict(climate_mean,
                              model_global[["glm_model"]],
                              type="response")
predictions_pulse_1 <- predict(climate_pulse_1,
                                   model_pulse_1[["glm_model"]],
                               type="response") 
predictions_pulse_2 <- predict(climate_pulse_2, 
                                   model_pulse_2[["glm_model"]],
                               type="response")

# Apply threshold
predictions_global_thresh <- predictions_global > model_global[["glm_thresh"]]
predictions_pulse_1_thresh <- predictions_pulse_1 > model_pulse_1[["glm_thresh"]]
predictions_pulse_2_thresh <- predictions_pulse_2 > model_pulse_2[["glm_thresh"]]

predictions_global <- mask(predictions_global, reclassify(predictions_global_thresh, c(-0.1, 0.1, NA)))
predictions_pulse_1 <- mask(predictions_pulse_1, reclassify(predictions_pulse_1_thresh, c(-0.1, 0.1, NA)))
predictions_pulse_2 <- mask(predictions_pulse_2, reclassify(predictions_pulse_2_thresh, c(-0.1, 0.1, NA)))

# Find out minimum and maximum values
min_value <- list(predictions_global,
                  predictions_pulse_1,
                  predictions_pulse_2) %>%
  sapply(function(x) x@data@min) %>% min()
max_value <- list(predictions_global,
                  predictions_pulse_1,
                  predictions_pulse_2) %>%
  sapply(function(x) x@data@max) %>% max()

# Set breaks for key (21)
key_breaks <- seq(floor(min_value * 0.9 * 100) / 100 ,
                  ceiling(max_value * 1.1 * 100) / 100, length.out = 21) 

# Visualise predictions

# Define function to visualise predictions
plot_preds_map <- function(base_raster,
                           main_title,
                           file_name,
                           palette = "default",
                           colour_key = F,
                           key_at = key_breaks,
                           to_file = F) {
  # set pallette if needed
  if(palette[1] == "default"){
    palette <- sequential_hcl(
      103, 
      "Inferno")
    palette <- palette[c(-1,-2,-3)]
  }
  # create plot
  map_plot <- levelplot(base_raster, 
                        margin = F,
                        main = main_title,
                        colorkey = colour_key,
                        at = key_at,
                        scales = list(draw = F),
                        xlab = NULL,
                        ylab = NULL,
                        par.settings = rasterTheme(
                          region = palette,
                          layout.widths = list(left.padding=-1,
                                               right.padding=-1),
                          layout.heights = list(top.padding=-1, 
                                                bottom.padding=-1)),
                        panel = function(...) {
                          panel.fill(col = "grey30")
                          panel.levelplot(...)
                        }) +
    latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 1)) +
    latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
    latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
    latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "uncertain",],
                                  col = col_grey_light, pch = shape_uncertain, cex = 0.8)) +
    latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_1",],
                                  col = col_pulse1, fill = col_pulse1_transparent, pch = shape_pulse_1)) +
    latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_2",],
                                  col = col_pulse1, fill = col_pulse2_transparent, pch = shape_pulse_2)) # +


  # Write out file if needed
  if(to_file == T){
    png(file_name, 
        width = 6,
        height = 3,
        units = "in",
        res = 300,
        type = "cairo",
        antialias = "none")
    print(map_plot)
    dev.off()
  }
  return(map_plot)
}

# |_ Mean ----
preds_map_global <- plot_preds_map(base_raster = predictions_global,
                                   main = NULL)

# |_ Pulse 1----
preds_map_pulse_1 <- plot_preds_map(base_raster = predictions_pulse_1,
               main = NULL)

# |_ Pulse 2 ----
preds_map_pulse_2 <- plot_preds_map(base_raster = predictions_pulse_2,
               main = NULL)

# |_ Colourkey ----
palette <- sequential_hcl(
  103, 
  "Inferno")
palette <- palette[c(-1,-2,-3)]
predictions_scale <-  levelplot(predictions_global, 
                         margin = F,
                         main = NULL,
                         colorkey = list(space="bottom", labels = list(cex = 0.5)),
                         at = key_breaks,
                         scales = list(draw = F),
                         xlab = NULL,
                         ylab = NULL,
                         par.settings = rasterTheme(
                           region = palette,
                           layout.widths = list(left.padding=0.5,
                                                right.padding=0.5),
                           layout.heights = list(top.padding=-1, 
                                                 bottom.padding=-1)),
                         panel = function(...) {
                           panel.fill(col = "grey30")
                           panel.levelplot(...)
                         }) 

# Export to file
png("figures/helper_figures/glm_predictions_scale.png", 
    width = 1.5,
    height = 3,
    units = "in",
    res = 300,
    type = "cairo",
    antialias = "none")
print(predictions_scale)
dev.off()

# Cut out scale
predictions_scale <- image_read("figures/helper_figures/glm_predictions_scale.png")
predictions_scale <- image_crop(predictions_scale, "450x180+0+480")
image_write(predictions_scale, "figures/helper_figures/glm_predictions_scale.png") 
predictions_scale <- ggdraw() + draw_image(predictions_scale,
                                    x = -0.105,
                                    y = 0.5,
                                    hjust = 0,
                                    vjust = 0.5,
                                    scale = 1)

# |_ Legend ----
predictions_legend <- levelplot(predictions_global, 
                         margin = F,
                         main = NULL,
                         colorkey = F,
                         at = seq(0,1,0.05),
                         scales = list(draw = F),
                         xlab = NULL,
                         ylab = NULL,
                         par.settings = rasterTheme(
                           region = palette,
                           layout.widths = list(left.padding=-1,
                                                right.padding=-1),
                           layout.heights = list(top.padding=-1, 
                                                 bottom.padding=-1)),
                         panel = function(...) {
                           panel.fill(col = "grey30")
                           panel.levelplot(...)
                         }) +
  
  latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "uncertain",],
                                col = col_grey_light, pch = shape_uncertain, cex = 0.8)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_1",],
                                col = col_pulse1, fill = col_pulse1_transparent, pch = shape_pulse_1)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_2",],
                                col = col_pulse1, fill = col_pulse2_transparent, pch = shape_pulse_2))  +

latticeExtra::layer({
  centre_x <- 27.5
  centre_y <- 58.3
  grid.rect(x=centre_x, y=centre_y,
            width=7.5, height=2.5,
            gp=gpar(fill="white",
                    col = "black",
                    alpha = 1),
            default.units='native')


  grid.points(#label = "Pulse 1",
    x=centre_x-2.75, y=centre_y+2.5/4,
    pch = shape_pulse_1,
    gp=gpar(cex = 0.45,
            col = col_pulse1,
            fill = col_pulse1_transparent),
    default.units='native')
  grid.points(#label = "Pulse 2",
    x=centre_x-2.75, y=centre_y,
    pch = shape_pulse_2,
    gp=gpar(cex = 0.45,
            col = col_pulse2,
            fill = col_pulse2_transparent),
    default.units='native')
  grid.points(#label = "Uncertain",
    x=centre_x-2.75, y=centre_y-2.5/4,
    pch = shape_uncertain,
    gp=gpar(cex = 0.65,
            col = col_grey_dark),
    default.units='native')

  grid.text(label = "Pulse 1",
            x=centre_x-1.75, y=centre_y+2.5/4,
            #just = "left",
            hjust = 0,
            vjust = 0.4,
            gp=gpar(cex=0.5,
                    col = "black"),
            default.units='native')
  grid.text(label = "Pulse 2",
            x=centre_x-1.75, y=centre_y,
            #just = "left",
            hjust = 0,
            vjust = 0.4,
            gp=gpar(cex=0.5,
                    col = "black"),
            default.units='native')
  grid.text(label = "Uncertain",
            x=centre_x-1.75, y=centre_y-2.5/4,
            #just = "left",
            hjust = 0,
            vjust = 0.4,
            gp=gpar(cex=0.5,
                    col = "black"),
            default.units='native')
})
# Export to file
png("figures/helper_figures/glm_predictions_legend.png", 
    width = 5,
    height = 3/6*5,
    units = "in",
    res = 300,
    type = "cairo",
    antialias = "none")
print(predictions_legend)
dev.off()

# Cut out legend
predictions_legend <- image_read("figures/helper_figures/glm_predictions_legend.png")
predictions_legend <- image_crop(predictions_legend, "263x152+1211+70")
image_write(predictions_legend, "figures/helper_figures/glm_predictions_legend.png") 
predictions_legend <- ggdraw() + draw_image(predictions_legend,
                                     x = 0.44,
                                     y = 0.6,
                                     hjust = 0.5,
                                     vjust = 0.5,
                                     scale = 0.75)

predictions_legend

# 4) Time-Series for Predictions ----

# Set time range
start_cent <- 150-20
end_cent <- 130-20
time_range <- seq(end_cent, start_cent, 1)

# Crop and stack rasters
time_series_climate <- lapply(time_range, function(cent_bce){
  cat(paste0(cent_bce, ":\n"))
  temp_raster_name <- paste0("temp_", cent_bce, "_BCE")
  precip_raster_name <- paste0("precip_", cent_bce, "_BCE")
  climate_bce <- stack(crop(temp[[temp_raster_name]],area_of_interest),
                       crop(precip[[precip_raster_name]], area_of_interest))
}) %>% setNames(time_range)

# Define function to calculate predictions for a given model
predict_time_series <- function(glm_model_list, model_name, file_name){
  # Status
  cat(paste0("Starting ", model_name, " Model ...\n"))
  
  # create folder
  out_folder <- paste0("figures/supplementary/glm_preds_ts_", file_name)
  dir.create(out_folder)
  
  # Calculate predictions
  cat(paste0("Calculating predictions ... \n"))
  glm_model <- glm_model_list[["glm_model"]]
  glm_thresh <- glm_model_list[["glm_thresh"]]
  
  preds_time_series <- lapply(time_series_climate, function(climate_raster){
    cat(paste0(names(climate_raster), " \n"))
    
    # update climate raster names to match model
    climate_raster <- setNames(climate_raster, c("mean_temp", "mean_precip"))
    # Predict for both time periods
    predictions <- predict(climate_raster,
                           glm_model,
                           type = "response")
    # Apply threshold
    predictions_thresh <- predictions > glm_thresh
    predictions <- mask(predictions, reclassify(predictions_thresh, c(-0.1, 0.1, NA) ))
    
    return(predictions)
  }) %>% setNames(time_range)
  
  # Generate Plots
  cat(paste0("\nGenerating maps ... \n"))
  lapply(seq_along(time_range), function(index){
    k_year_bp_start <- formatC(time_range[index] / 10 + 2.1, 1, format = "f")
    k_year_bp_end <- formatC(time_range[index] / 10 + 2, 1, format = "f")
    cat(paste0( k_year_bp_start, "-", k_year_bp_end, "k BP\n"))
    plot_preds_map(base_raster = preds_time_series[[index]],
                   main = NULL, 
                   file_name = paste0(out_folder, "/predictions_", 
                                      k_year_bp_start, "-", k_year_bp_end, "k_BP.png"),
                   to_file = T
                 #palette = c("grey30", "#f38f32")
                 )
  })
  
  cat(paste0("Done. \n"))
  return(preds_time_series)
}

# Generate time-series
ts_global <- predict_time_series(model_global, "Global", "global")
ts_pulse_1 <- predict_time_series(model_pulse_1, "Pulse 1", "pulse_1")
ts_pulse_2 <- predict_time_series(model_pulse_2, "Pulse 2", "pulse_2")

# 5) Final Figures ----

# |_ Figure S3 - Suitabillity Predictions ----
# List of plots 
plot_list <- c("preds_map_global",
               "preds_map_pulse_1",
               "preds_map_pulse_2")

# Check all object exist
lapply(plot_list, function(x) exists(x))

# Helper function to add standardised title to plot
add_title <- function(plot, title, face = "plain", font_size = 8, 
                      rel_height = 0.075, top_margin = 4, left_margin = 15.5){
  title_panel <- ggdraw() + 
    draw_label(
      title,
      fontface = face,
      x = 0,
      hjust = 0,
      size = font_size
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(top_margin, 0, 0, left_margin)
    )
  plot_with_title <- plot_grid(
    title_panel, plot,
    ncol = 1,
    rel_heights = c(rel_height, 1)
  )
  return(plot_with_title)
}
# Combine legend and colour scale
predictions_legend_panel <- plot_grid(add_title(predictions_scale, 
                                         "Predicted Propability (glm)",
                                         "plain",
                                         8,
                                         rel_height = 0.2,
                                         top_margin = 2,
                                         left_margin = 15.5),
                               predictions_legend,
                               rel_widths = c(0.66,0.33))

# Add titles and combine into one figure
figure_s3 <- plot_grid(
  plotlist = list(add_title(preds_map_global, 
                           "Global Period Model Predictions 14.5k-14.1k BP",
                           top_margin = 3.5),
                 add_title(preds_map_pulse_1,
                           "Pulse 1 Model Predictions 14.5k-14.3k BP",
                           top_margin = 3.5),
                 add_title(preds_map_pulse_2,
                           "Pulse 2 Model Predictions 14.2k-14.1k BP",
                           top_margin = 3.5),
                 predictions_legend_panel),
  labels = "AUTO",
  ncol = 1,
  label_size = 8,
  rel_heights = c(1,1,1,0.5))
save_plot("figures/supplementary/figure_S3_suitability_predictions_glm.png", 
          figure_s3,
          base_height = 2,
          ncol = 1,
          nrow = 3,
          base_asp = 1.85,
          bg = "white")
save_plot("figures/supplementary/figure_S3_suitability_predictions_glm.eps", 
          figure_s3,
          base_height = 2,
          ncol = 1,
          nrow = 3,
          base_asp = 1.85,
          bg = "white")

# |_ Figure S4 - Suitability Time-Series ----
# read in time-series plots (centuries 14.7-14.1)
ts_maps_global <- lapply(rev(list.files("figures/supplementary/glm_preds_ts_global/", pattern = ".png", full.names = T)[11:17]),
                         function(x) ggdraw() + draw_image(image_read(x), scale = 0.95))
ts_maps_pulse_1 <- lapply(rev(list.files("figures/supplementary/glm_preds_ts_pulse_1/", pattern = ".png", full.names = T)[11:17]),
                          function(x) ggdraw() + draw_image(image_read(x), scale = 0.95))
ts_maps_pulse_2 <- lapply(rev(list.files("figures/supplementary/glm_preds_ts_pulse_2/", pattern = ".png", full.names = T)[11:17]),
                          function(x) ggdraw() + draw_image(image_read(x), scale = 0.95))
# Arrange as panels in column
ts_maps_global <- plot_grid(plotlist = ts_maps_global,
                            ncol = 1)
ts_maps_pulse_1 <- plot_grid(plotlist = ts_maps_pulse_1,
                            ncol = 1)
ts_maps_pulse_2 <- plot_grid(plotlist = ts_maps_pulse_2,
                            ncol = 1)
# Generate a timeline plot
ts_time_scale <- ggplot() + 
  scale_y_continuous(limits = c(14.0, 14.7),
                     breaks = seq(14.0, 14.7, 0.1),
                     labels = paste0(formatC(seq(14.0, 14.7, 0.1), digits = 1, format = "f"), " kyr BP"),
                     expand = c(0,0),
                     ) +
  geom_blank() +
  theme_cowplot() +
  theme(axis.line.x = element_line(colour = NA),
        axis.text.x = element_text(colour = NA),
        axis.ticks.x = element_line(colour = NA),
        axis.text.y = element_text(size = 12))

# Update the legend
predictions_scale <- image_read("figures/helper_figures/glm_predictions_scale.png")
predictions_scale <- ggdraw() + draw_image(predictions_scale,
                                           x = -0.231,
                                           y = 0.5,
                                           hjust = 0,
                                           vjust = 0.5,
                                           scale = 1)
predictions_legend_panel <- plot_grid(add_title(predictions_scale, 
                                                "Predicted propability (glm)",
                                                "plain",
                                                12,
                                                rel_height = 0.2,
                                                top_margin = 2,
                                                left_margin = 7),
                                      plot_grid(predictions_legend),
                                      rel_widths = c(0.66,0.3))
# Assemble plot
figure_s4_body <- plot_grid(ts_time_scale,
                      ts_maps_global,
                      ts_maps_pulse_1,
                      ts_maps_pulse_2,
                      ncol = 4,
                      rel_widths = c(0.3,1,1,1))
# add labels and legend
figure_s4_labels <- lapply(
  c("Global Model", "Pulse 1 Model", "Pulse 2 Model"),
  function(title){
    ggdraw() + 
      draw_label(
        title,
        fontface = "bold",
        x = 0,
        hjust = 0,
        size = 12) + 
      theme(plot.margin = margin(4, 0, 0, 7))
    
  }
)
figure_s4_labels <- plot_grid(
  NULL, figure_s4_labels[[1]],
  figure_s4_labels[[2]], 
  figure_s4_labels[[3]],
  rel_widths =c(0.3,1,1,1),
  nrow= 1)
figure_s4_legend <- plot_grid(
  NULL, predictions_legend_panel, NULL,
  rel_widths = c(0.3,1.8,1.2),
  nrow = 1
)
figure_s4 <- plot_grid(figure_s4_labels,
                      figure_s4_body,
                      figure_s4_legend,
                      ncol = 1,
                      rel_heights = c(0.015,1,0.1))
save_plot("figures/supplementary/figure_S4-suitability_time_series_glm.png",
          figure_s4,
          base_height = 1.9,
          ncol = 3,
          nrow = 7,
          base_asp = 2,
          bg = "white")
save_plot("figures/supplementary/figure_S4-suitability_time_series_glm.eps",
          figure_s4,
          base_height = 1.9,
          ncol = 3,
          nrow = 7,
          base_asp = 2,
          bg = "white")

# |_ Animations ----

# Generate Animation
animate_time_series <- function(ts_name){
  # Status
  cat(paste0("Starting ", ts_name, " Model ...\n"))
  
  # set folder
  out_folder <- paste0("figures/supplementary/glm_preds_ts_", ts_name)
  
  cat(paste0("\nAnimating time-series ... \n"))
  # List files and loade as images
  imgs <- list.files(out_folder, 
                     pattern = ".png", full.names = TRUE)
  
  img_list <- lapply(imgs, image_read)
  
  # Add legend and Time Stamp
  prediction_legend <- image_read("figures/helper_figures/glm_predictions_legend.png")
  max_index <- max(seq_along(time_range))
  img_list <- lapply(seq_along(time_range), function(index){
    cat(paste0(index,"."))
    k_year_bp_start <- formatC(time_range[index] / 10 + 2.1, 1, format = "f")
    k_year_bp_end <- formatC(time_range[index] / 10 + 2, 1, format = "f")
    image1 <- img_list[[index]] %>% image_composite(prediction_legend, offset = "+1490+100") %>%
      image_annotate(paste0(k_year_bp_start, " kyr BP"), size = 70, color = "white", location = "+60+90") %>%
      image_crop("1800x796+0+52")
    return(image1)
    # if(index == 10) image_morphed <- image1
    # else {
    #   image2 <- img_list[[index + 1]] %>% image_composite(prediction_legend, offset = "+1490+100") %>%
    #   image_annotate(paste0(k_year_bp_end, " kyr BP"), size = 70, color = "white", location = "+60+90")
    # 
    #   image_morphed <- image_morph(c(image2,image1), frames = 5)
    # }
    # return(image_morphed)
    })
  
  # Status
  cat("\nJoining images...")
  # join images 
  img_joined <- image_join(rev(img_list))
  
  # Status
  cat("\nExporting GIF...")
  # Export as gif
  image_write_gif(image = img_joined,
            path = paste0(out_folder, "/glm_preds_time_series_", ts_name, ".gif"),
            delay = 0.5)
  cat("\nDone.")
  return(NULL)
}
animate_time_series("global")
animate_time_series("pulse_1")
animate_time_series("pulse_2")
