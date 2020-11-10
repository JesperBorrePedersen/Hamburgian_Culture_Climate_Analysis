# Hamburgian Culture Climate Analysis - First data Exploraiton 
# Jakob Assmann j.assmann@bios.au.dk 17 August 2020

## 1) Housekeeping ----
# Dependencies
library(sf)
library(raster)
library(tidyverse)
library(cowplot)
library(rasterVis)
library(colorspace)
library(rnaturalearth)
library(rnaturalearthdata)
library(dismo)

rasterOptions(progress = "text")

# Load in CHELSA LGM files
temp <- list.files(
  "O:/Nat_Ecoinformatics/B_Read/World/Climate/PaleoClimate/Original/CHELSA_TRACE/BIO_01/",
  pattern = "*.tif",
  full.names = T)
precip <- list.files(
  "O:/Nat_Ecoinformatics/B_Read/World/Climate/PaleoClimate/Original/CHELSA_TRACE/BIO_12/",
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
# # Make time ranges interpretable
# hamburgian_sites$CalBP.Pulses_begin <- as.numeric(
#   gsub("^([0-9]*\\.[0-9]*) -.*", 
#        "\\1", 
#        
#        hamburgian_sites$CalBP.Pulses)) * 1000
# hamburgian_sites$CalBP.Pulses_end <- as.numeric(
#   gsub(".* - ([0-9]*\\.[0-9]*$)",
#        "\\1", 
#        
#        hamburgian_sites$CalBP.Pulses))* 1000
# hamburgian_sites$CalBP.Broad_begin <- as.numeric(
#   gsub("^([0-9]*\\.[0-9]*) -.*", 
#        "\\1", 
#        hamburgian_sites$CalBP.Broad))* 1000
# hamburgian_sites$CalBP.Broad_end <- as.numeric(
#   gsub(".* - ([0-9]*\\.[0-9]*$)", 
#        "\\1", 
#        hamburgian_sites$CalBP.Broad)) * 1000

# Extract data from layers
temp_df <- as.data.frame(raster::extract(temp, 
                                         as_Spatial(hamburgian_sites)))
precip_df <-  as.data.frame(raster::extract(precip, 
                                            as_Spatial(hamburgian_sites)))

# Add site column
temp_df$site <- hamburgian_sites$site
precip_df$site <- hamburgian_sites$site

# Pivot longer
temp_df <- pivot_longer(temp_df, 1:(ncol(temp_df)-1),
                        names_to = "names_temp", 
                        values_to = "temp")

precip_df <- pivot_longer(precip_df, 1:(ncol(precip_df)-1), 
                          names_to = "names_precip", 
                          values_to = "precip")

# Add years BP from look up table
temp_df <- full_join(temp_df, temp_look_up)
precip_df <- full_join(precip_df, precip_look_up)

# Save dfs
save(temp_df, file = "data/temp_df.Rda")
save(precip_df, file = "data/precip_df.Rda")

# Load dfs
load("data/temp_df.Rda")
load("data/precip_df.Rda")

## 2) Time-series Plots ----
# Calculate key stats

# Mean temp / precip in the two Hamburgian Culture periods
mean_temp <- temp_df %>% 
  group_by(site) %>%
  filter(year_BP >= 13800 & year_BP <= 15000) %>%
  summarise(mean_temp = mean(temp))
mean_temp_pulse_1 <- temp_df %>% 
  group_by(site) %>%
  filter(year_BP >= 14300 & year_BP <= 14500) %>%
  summarise(mean_temp_pulse_1 = mean(temp))
mean_temp_pulse_2 <- temp_df %>% 
  group_by(site) %>%
  filter(year_BP >= 14100 & year_BP <= 14200) %>%
  summarise(mean_temp_pulse_2 = mean(temp))


mean_precip <- precip_df %>% 
  group_by(site) %>%
  filter(year_BP >= 13800 & year_BP <= 15000) %>%
  summarise(mean_precip = mean(precip))
mean_precip_pulse_1 <- precip_df %>% 
  group_by(site) %>%
  filter(year_BP >= 14300 & year_BP <= 14500) %>%
  summarise(mean_precip_pulse_1 = mean(precip))
mean_precip_pulse_2 <- precip_df %>% 
  group_by(site) %>%
  filter(year_BP >= 14100 & year_BP <= 14200) %>%
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


# set mean temp to NA where no evidence is avialable
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

# |_ Temp Time-Series Plot ----
temp_plot <- ggplot(temp_df %>% filter(year_BP >= 13000 & year_BP <= 15500) %>%
                      filter(site %in% unique(hamburgian_sites$site)) %>%
                      mutate(site = ordered(site, levels = levels(fct_reorder(hamburgian_sites$site, hamburgian_sites$lat)))),
       aes(x = year_BP,
           y = temp,
           colour = site)) +
  geom_line() +
  geom_segment(data = hamburgian_sites,
             mapping = aes(x = pulse_1_end,
                           xend = pulse_1_start,
                           y = mean_temp_pulse_1, 
                           yend = mean_temp_pulse_1),
             colour = "black") +
  geom_segment(data = hamburgian_sites,
    mapping = aes(x = pulse_2_end,
                  xend = pulse_2_start,
                  y = mean_temp_pulse_2, 
                  yend = mean_temp_pulse_2),
    colour = "black") +
  geom_text(data = hamburgian_sites,
            mapping = aes(x = pulse_1_start + 100,
                          y = mean_temp_pulse_1,
                          label = paste0("mean:\n", round(mean_temp_pulse_1,1), "°C")),
            vjust = 0.5,
            hjust = 1,
            size = 2,
            colour = "black") +
  geom_text(data = hamburgian_sites,
            mapping = aes(x = pulse_2_start - 200,
                          y = mean_temp_pulse_2,
                          label = paste0("mean:\n", round(mean_temp_pulse_2,1), "°C")),
            colour = "black",
            vjust = 0.5,
            hjust = 0,
            size = 2) +
  labs(x = "Year BP", y = "Annual Mean Temp °C (Bio01)") +
  scale_x_reverse(limits = c(15500, 13000),
                  breaks = seq(15500, 13000, -500),
                  labels = rev(c("13.0k", "13.5k", "14.0k", "14.5k", "15.0k", "15.5k"))) +
  scale_y_continuous(limits = c(-7, 10),
                     breaks = c(-5,0,5,10)) +
  scale_colour_discrete_qualitative(palette = "Dark2") +
  theme_cowplot(8) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

save_plot("figures/figure1_temp.png",
          temp_plot,
          base_aspect_ratio = 1.6,
          base_height = 4)

# Plot per site
temp_plot <- temp_plot + facet_wrap(vars(site),
                                    scales = "free") +
  theme(axis.title = element_text(face = "bold"))
save_plot("figures/figure1_temp_by_site.png",
          temp_plot,
          base_asp = 1.6,
          base_height = 20)

# |_ Precip Time-Series Plot ----
precip_plot <- ggplot(precip_df %>% filter(year_BP >= 13000 & year_BP <= 15500) %>%
                      mutate(site = ordered(site, levels = levels(fct_reorder(hamburgian_sites$site, hamburgian_sites$lat)))),
                    aes(x = year_BP,
                        y = precip,
                        colour = site)) +
  geom_line() +
  geom_segment(data = hamburgian_sites,
               mapping = aes(x = pulse_1_end,
                             xend = pulse_1_start,
                             y = mean_precip_pulse_1, 
                             yend = mean_precip_pulse_1),
               colour = "black") +
  geom_segment(data = hamburgian_sites,
               mapping = aes(x = pulse_2_end,
                             xend = pulse_2_start,
                             y = mean_precip_pulse_2, 
                             yend = mean_precip_pulse_2),
               colour = "black") +
  geom_text(data = hamburgian_sites,
            mapping = aes(x = pulse_1_start + 100,
                          y = mean_precip_pulse_1,
                          label = paste0("mean:\n", round(mean_precip_pulse_1,1), "mm")),
            vjust = 0.5,
            hjust = 1,
            size = 2,
            colour = "black") +
  geom_text(data = hamburgian_sites,
            mapping = aes(x = pulse_2_start - 200,
                          y = mean_precip_pulse_2,
                          label = paste0("mean:\n", round(mean_precip_pulse_2,1), "mm")),
            colour = "black",
            vjust = 0.5,
            hjust = 0,
            size = 2,) +
  labs(x = "Year BP", y = "Annual Precipitation [mm] (Bio12)") +
  scale_x_reverse(limits = c(15500, 13000),
                  breaks = seq(15500, 13000, -500),
                  labels = rev(c("13.0k", "13.5k", "14.0k", "14.5k", "15.0k", "15.5k"))) +
  scale_y_continuous(limits = c(400, 1100),
                     breaks = seq(400,1200,200)) +
  scale_colour_discrete_qualitative(palette = "Dark2") +
  theme_cowplot(8) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

save_plot("figures/figure1_precip.png",
          precip_plot,
          base_aspect_ratio = 1.6,
          base_height = 4)

# Plot per site
precip_plot <- precip_plot + facet_wrap(vars(site),
                                    scales = "free") +
  theme(axis.title = element_text(face = "bold"))
save_plot("figures/figure1_precip_by_site.png",
          precip_plot,
          base_asp = 1.6,
          base_height = 20)

## 3) Raster Maps ----
area_of_interest <- extent(c(-11,32,49,60)) 
time_range_bp <- data.frame(min_bp = 15000-2000,
                         max_bp = 13800-2000,
                         min_cent = 150-20,
                         max_cent = 138-20,
                         pulse_1_start = 14500-2000,
                         pulse_1_end = 14300-2000,
                         pulse_2_start = 14200-2000,
                         pulse_2_end = 14100-2000,
                         pulse_1_min_cent = 145-20,
                         pulse_1_max_cent = 143-20,
                         pulse_2_min_cent = 142-20,
                         pulse_2_max_cent = 141-20)



## calculate average values for the broad time-range
# Temperature
temp_time_range_bp_all_cents <- paste0(
  "temp_",
  seq(time_range_bp$max_cent, time_range_bp$min_cent, 1),
  "_BCE")
hamburgian_mean_temp <- mean(crop(temp[[temp_time_range_bp_all_cents]], area_of_interest))

temp_time_range_bp_pulse_1 <- paste0(
  "temp_",
  seq(time_range_bp$pulse_1_max_cent, time_range_bp$pulse_1_min_cent, 1),
  "_BCE")
hamburgian_mean_temp_pulse_1 <- mean(crop(temp[[temp_time_range_bp_pulse_1]], area_of_interest))

temp_time_range_bp_pulse_2 <- paste0(
  "temp_",
  seq(time_range_bp$pulse_2_max_cent, time_range_bp$pulse_2_min_cent, 1),
  "_BCE")
hamburgian_mean_temp_pulse_2 <- mean(crop(temp[[temp_time_range_bp_pulse_2]], area_of_interest))

# Precip
precip_time_range_bp_all_cents <- paste0(
  "precip_",
  seq(time_range_bp$max_cent, time_range_bp$min_cent, 1),
  "_BCE")
hamburgian_mean_precip <- mean(crop(precip[[precip_time_range_bp_all_cents]], area_of_interest))

precip_time_range_bp_pulse_1 <- paste0(
  "precip_",
  seq(time_range_bp$pulse_1_max_cent, time_range_bp$pulse_1_min_cent, 1),
  "_BCE")
hamburgian_mean_precip_pulse_1 <- mean(crop(precip[[precip_time_range_bp_pulse_1]], area_of_interest))

precip_time_range_bp_pulse_2 <- paste0(
  "precip_",
  seq(time_range_bp$pulse_2_max_cent, time_range_bp$pulse_2_min_cent, 1),
  "_BCE")
hamburgian_mean_precip_pulse_2 <- mean(crop(precip[[precip_time_range_bp_pulse_2]], area_of_interest))

## Save rasters for future use
writeRaster(hamburgian_mean_temp, 
            "data/mean_rasters/hamburgian_mean_temp.tif",
            overwrite = T)
writeRaster(hamburgian_mean_temp_pulse_1, 
            "data/mean_rasters/hamburgian_mean_temp_pulse_1.tif",
            overwrite = T)
writeRaster(hamburgian_mean_temp_pulse_2, 
            "data/mean_rasters/hamburgian_mean_temp_pulse_2.tif",
            overwrite = T)

writeRaster(hamburgian_mean_precip, 
            "data/mean_rasters/hamburgian_mean_precip.tif",
            overwrite = T)
writeRaster(hamburgian_mean_precip_pulse_1, 
            "data/mean_rasters/hamburgian_mean_precip_pulse_1.tif",
            overwrite = T)
writeRaster(hamburgian_mean_precip_pulse_2, 
            "data/mean_rasters/hamburgian_mean_precip_pulse_2.tif",
            overwrite = T)

# Load rasters
hamburgian_mean_temp <- raster("data/mean_rasters/hamburgian_mean_temp.tif")
hamburgian_mean_temp_pulse_1 <- raster("data/mean_rasters/hamburgian_mean_temp_pulse_1.tif")
hamburgian_mean_temp_pulse_2 <- raster("data/mean_rasters/hamburgian_mean_temp_pulse_2.tif")

hamburgian_mean_precip <- raster("data/mean_rasters/hamburgian_mean_precip.tif")
hamburgian_mean_precip_pulse_1 <- raster("data/mean_rasters/hamburgian_mean_precip_pulse_1.tif")
hamburgian_mean_precip_pulse_2 <- raster("data/mean_rasters/hamburgian_mean_precip_pulse_2.tif")

## Plot mean maps

# Generate political boundaries for maps
# countries <- ne_countries(scale = "medium", returnclass = "sp")
# countries <- spTransform(countries, crs(hamburgian_mean_temp))
# countries <- crop(countries, hamburgian_mean_temp)

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
  
# Transform hamburgian sites
hamburgian_sites <- st_transform(hamburgian_sites, st_crs(hamburgian_mean_temp))
hamburgian_sites_sp <- as_Spatial(hamburgian_sites)


# |_ Temp Maps ----
# Define Function to plot temp map
plot_ham_temp_map <- function(base_raster,
                              main_title,
                              file_name){
  cat(paste0("Plotting ", file_name, "...\n"))
  
  png(file_name, 
      width = 6,
      height = 3,
      units = "in",
      res = 300)
  print(levelplot(get(base_raster, envir = .GlobalEnv), 
                  margin = F,
                  main = main_title,
                  at = seq(-12, 12, 1)) +
          latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 1)) +
          latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
          latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
          latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "uncertain",],
                                        col = "gray90", pch = 63, cex = 0.8)) +
          latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_1",],
                                        col = "dodgerblue2", pch = 1, cex = 1.05)) +
          latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_2",],
                                        col = "black", pch = 3)) +
          
          latticeExtra::layer({
            centre_x <- 27.5
            centre_y <- 58.3
            grid.rect(x=centre_x, y=centre_y,
                      width=7.5, height=2.5,
                      gp=gpar(fill="white",
                              col = "black",
                              alpha = 0.9),
                      default.units='native')
            
            
            grid.points(#label = "Pulse 1",
              x=centre_x-2.75, y=centre_y+2.5/4,
              pch = 1,
              gp=gpar(cex = 0.65,
                      col = "dodgerblue2"),
              default.units='native')
            grid.points(#label = "Pulse 2",
              x=centre_x-2.75, y=centre_y,
              pch = 3,
              gp=gpar(cex = 0.4,
                      col = "black"),
              default.units='native')
            grid.points(#label = "Uncertain",
              x=centre_x-2.75, y=centre_y-2.5/4,
              pch = 63,
              gp=gpar(cex = 0.65,
                      col = "grey30"),
              default.units='native')
            
            grid.text(label = "Pulse 1",
                      x=centre_x-1.75, y=centre_y+2.5/4,
                      #just = "left",
                      hjust = 0,
                      vjust = 0.4,
                      gp=gpar(cex=0.5,
                              col = "black"),
                      default.units='native')
            grid.text(label = "Pusle 2",
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
          }))
  dev.off()
  return("Done.")
}

# Mean preciptation map 15k-13.5k BP
plot_ham_temp_map(base_raster = "hamburgian_mean_temp",
                  main_title = "Mean Temperature 15k-13.5k BP (°C)",
                  file_name = "figures/mean_temp.png")
# Mean temperature map Pulse 1 14.5k-14.3k BP
plot_ham_temp_map(base_raster = "hamburgian_mean_temp_pulse_1",
                  main_title = "Mean Temperature 14.5k-14.3k BP (°C)",
                  file_name = "figures/mean_temp_pulse_1.png")

# Mean temperature map Pulse 2 14.2k-14.1k BP
plot_ham_temp_map(base_raster = "hamburgian_mean_temp_pulse_2",
                  main_title = "Mean Temperature 14.2k-14.1k BP (°C)",
                  file_name = "figures/mean_temp_pulse_2.png")

# |_ Precip Maps ----

# Define Function to plot precip map
plot_ham_precip_map <- function(base_raster,
                              main_title,
                              file_name){
  cat(paste0("Plotting ", file_name, "...\n"))

  png(file_name, 
      width = 6,
      height = 3,
      units = "in",
      res = 300)
  print(levelplot(get(base_raster, envir = .GlobalEnv), 
                  margin = F,
                  main = main_title,
                  at = seq(0, 2500, 2500/20),
                  par.settings = rasterTheme(
                    region = sequential_hcl(
                      100, 
                      "Blues",
                      rev = T))) + 
          latticeExtra::layer(sp.polygons(land_for_maps, col = "white", alpha = 1)) +
          latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
          latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
          latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "uncertain",],
                                        col = "grey30", pch = 63, cex = 0.8)) +
          latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_1",],
                                        col = "dodgerblue2", pch = 1, cex = 1.05)) +
          latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_2",],
                                        col = "black", pch = 3)) +
          latticeExtra::layer({
            centre_x <- 27.5
            centre_y <- 58.3
            grid.rect(x=centre_x, y=centre_y,
                      width=7.5, height=2.5,
                      gp=gpar(fill="white", 
                              col = "black",
                              alpha = 0.9),
                      default.units='native')
            
            grid.points(#label = "Hamburgian",
              x=centre_x-2.75, y=centre_y+2.5/4,
              pch = 1,
              gp=gpar(cex = 0.65, 
                      col = "dodgerblue2"),
              default.units='native')
            grid.points(#label = "Havelte Group",
              x=centre_x-2.75, y=centre_y,
              pch = 3,
              gp=gpar(cex = 0.4, 
                      col = "black"),
              default.units='native')
            grid.points(#label = "Possibly Ham.",
              x=centre_x-2.75, y=centre_y-2.5/4,
              pch = 63,
              gp=gpar(cex = 0.65, 
                      col = "grey30"),
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
          }))
        dev.off()
}

# Mean preciptation map 15k-13.5k BP
plot_ham_precip_map(base_raster = "hamburgian_mean_precip",
                  main_title = "Mean Precipitation 15k-13.5k BP (mm)",
                  file_name = "figures/mean_precip.png")

# Mean preciptation map 14.5k-14.3k BP
plot_ham_precip_map(base_raster = "hamburgian_mean_precip_pulse_1",
                    main_title = "Mean Precipitation 14.5k-14.3k BP (mm)",
                    file_name = "figures/mean_precip_pulse_1.png")


# Mean preciptation map 14.7k-14.5k BP
plot_ham_precip_map(base_raster = "hamburgian_mean_precip_pulse_2",
                    main_title = "Mean Precipitation 14.2k-14.1k BP (mm)",
                    file_name = "figures/mean_precip_pulse_2.png")

# 3) Temp vs. Precip Plot ---- 

# Cross check with extraction from rasters
mean_both_ham_extract <- data.frame(
  site = unique(hamburgian_sites$site),
  mean_temp = raster::extract(hamburgian_mean_temp, as_Spatial(hamburgian_sites)),
  mean_precip = raster::extract(hamburgian_mean_precip, as_Spatial(hamburgian_sites)),
  mean_temp_pulse_1 = raster::extract(hamburgian_mean_temp_pulse_1, as_Spatial(hamburgian_sites)),
  mean_precip_pulse_1 = raster::extract(hamburgian_mean_precip_pulse_1, as_Spatial(hamburgian_sites)),
  mean_temp_pulse_2 = raster::extract(hamburgian_mean_temp_pulse_2, as_Spatial(hamburgian_sites)),
  mean_precip_pulse_2 = raster::extract(hamburgian_mean_precip_pulse_2, as_Spatial(hamburgian_sites))
) %>% 
  as_tibble() %>%
  mutate(site = ordered(site, levels = levels(fct_reorder(hamburgian_sites$site, hamburgian_sites$lat))))
all_equal(hamburgian_sites %>%
            st_drop_geometry() %>% 
            select(site,
                 mean_temp,
                 mean_precip) %>% tibble(),
          mean_both_ham_extract %>% 
            select(site,
                   mean_temp,
                   mean_precip))

set.seed(6)

# Mask rasters
hamburgian_mean_temp_masked <- mask(
  hamburgian_mean_temp,
  land_for_maps)
hamburgian_mean_precip_masked <- mask(
  hamburgian_mean_precip,
  land_for_maps)

hamburgian_mean_temp_pulse_1_masked <- mask(
  hamburgian_mean_temp_pulse_1,
  land_for_maps)
hamburgian_mean_precip_pulse_1_masked <- mask(
  hamburgian_mean_precip_pulse_1,
  land_for_maps)

hamburgian_mean_temp_pulse_2_masked <- mask(
  hamburgian_mean_temp_pulse_2,
  land_for_maps)
hamburgian_mean_precip_pulse_2_masked <- mask(
  hamburgian_mean_precip_pulse_2,
  land_for_maps)

# Optain random sample for land areas from masked rasters
set.seed(56)
random_sample <- data.frame(
  site = paste0("Random_", 1:100000),
  mean_temp = sampleRandom(hamburgian_mean_temp_masked, 100000),
  mean_precip = sampleRandom(hamburgian_mean_precip_masked, 100000),
  mean_temp_pulse_1 = sampleRandom(hamburgian_mean_temp_pulse_1_masked, 100000),
  mean_precip_pulse_1 = sampleRandom(hamburgian_mean_precip_pulse_1_masked, 100000),
  mean_temp_pulse_2 = sampleRandom(hamburgian_mean_temp_pulse_2_masked, 100000),
  mean_precip_pulse_2 = sampleRandom(hamburgian_mean_precip_pulse_2_masked, 100000))

# Plot mean
temp_vs_precip_plot <- ggplot() +  
  geom_point(data = random_sample,
             aes(x = mean_temp,
                 y = mean_precip),
             colour = "grey",
             fill = "grey",
             alpha = 0.25,
             shape = 21
  ) +
  geom_point(data = hamburgian_sites %>% mutate(chron_association = ordered(chron_association,
                                                                            levels = c("uncertain",
                                                                                       "pulse_2",
                                                                                       "pulse_1"))) %>%
               arrange(chron_association),
             aes(x = mean_temp,
                 y = mean_precip,
                 colour = chron_association,
                 fill = chron_association),
             shape = 21,
             alpha = 0.75) +
  labs(x = "Mean Temperature (°C)",
       y = "Mean Precipitation (mm)") +
  scale_x_continuous(limits = c(-7, 10),
                     breaks = seq(-6, 10, 2)) +
  scale_y_continuous(limits = c(400, 2500)) +
  scale_colour_manual(values = c("grey30", "black", "dodgerblue2")) +
  scale_fill_manual(values = c("grey30", "black", "dodgerblue2")) +
  annotate("text", x = 8, y = 2500, 
           colour = "dodgerblue2", hjust = 0, vjust = 0.4,
           label = "Pulse 1") +
  annotate("text", x = 8, y = 2375,
           colour = "black", hjust = 0, vjust = 0.4,
           label = "Pulse 2") +
  annotate("text", x = 8, y = 2250,
           colour = "grey30", hjust = 0, vjust = 0.4,
           label = "Uncertain") +
  annotate("point", x = 7.7, y = 2500, 
           colour = "dodgerblue2") +
  annotate("point", x = 7.7, y = 2375,
           colour = "black") +
  annotate("point", x = 7.7, y = 2250,
           colour = "grey30") +
  theme_cowplot(14) +
  theme(legend.position = "none")
save_plot("figures/temp_vs_precip.png",
          temp_vs_precip_plot,
          base_aspect_ratio = 1.6,
          base_height = 4)

# plot the two pusles
temp_vs_precip_plot_by_pulse <- ggplot() +  
  geom_point(data = random_sample,
             aes(x = mean_temp_pulse_1,
                 y = mean_precip_pulse_1),
             colour = "lightskyblue2",
             fill = "lightskyblue2",
             alpha = 0.25,
             shape = 21
             ) +
  geom_point(data = random_sample,
             aes(x = mean_temp_pulse_2,
                 y = mean_precip_pulse_2),
             colour = "grey",
             fill = "grey",
             alpha = 0.25,
             shape = 21
             ) +
  geom_point(data = hamburgian_sites,
             aes(x = mean_temp_pulse_1,
                 y = mean_precip_pulse_1),
             colour = "dodgerblue2",
             shape = 16
             ) +
  geom_point(data = hamburgian_sites,
             aes(x = mean_temp_pulse_2,
                 y = mean_precip_pulse_2),
             colour = "black",
             shape = 16
             ) +
  labs(x = "Mean Temperature (°C)",
       y = "Mean Precipitation (mm)") +
  scale_x_continuous(limits = c(-7, 14),
                     breaks = seq(-6, 14, 2)) +
  scale_y_continuous(limits = c(400, 2500)) +
  annotate("text", x = 12, y = 2500, 
           colour = "dodgerblue2", hjust = 0, vjust = 0.4,
           label = "Pulse 1") +
  annotate("text", x = 12, y = 2400,
           colour = "black", hjust = 0, vjust = 0.4,
           label = "Pulse 2") +
  annotate("point", x = 11.7, y = 2500, 
           colour = "dodgerblue2") +
  annotate("point", x = 11.7, y = 2400,
           colour = "black") +
  theme_cowplot(14) 
save_plot("figures/temp_vs_precip_by_pulse.png",
          temp_vs_precip_plot_by_pulse,
          base_aspect_ratio = 1.6,
          base_height = 4)

# 4) BIOCLIM Model ----

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

# Quality control
view(bioclim_training_all)
# bioclim_training_all <- filter(bioclim_training_all, period == "pulse_2")

# Retain climate columns only
bioclim_training <- select(bioclim_training_all, "mean_temp", "mean_precip")

# repeat for mean values
bioclim_training_mean  <- hamburgian_sites %>%
  st_drop_geometry() %>%
  filter(chron_association != "uncertain") %>%
  select(site, 
         lat,
         long,
         cell_id,
         mean_temp,
         mean_precip) %>%
  distinct(cell_id, mean_temp, mean_precip) %>%
  select(-cell_id)

# Prepare raster stacks for predictions
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
climate_pulse_1_masked <- mask(climate_pulse_1,
                                   land_for_maps)
climate_pulse_2_masked <- mask(climate_pulse_2,
                                   land_for_maps)

# Generate random background points for evaluation (absence data)
set.seed(50)
absence_data <- randomPoints(climate_pulse_1_masked,
                             200)
colnames(absence_data) <- c('lon', 'lat')
absence_data <- as.data.frame(absence_data)

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

# fill in absence data, half from classic, half from halvelte
absence_data$mean_temp[1:floor(nrow(absence_data)/2)] <- raster::extract(
  climate_pulse_1[["mean_temp"]], absence_data$cell_id[1:floor(nrow(absence_data)/2)])
absence_data$mean_precip[1:floor(nrow(absence_data)/2)] <- raster::extract(
  climate_pulse_1[["mean_precip"]], absence_data$cell_id[1:floor(nrow(absence_data)/2)])

absence_data$mean_temp[(floor(nrow(absence_data)/2)+1):nrow(absence_data)] <- raster::extract(
  climate_pulse_2[["mean_temp"]], absence_data$cell_id[(floor(nrow(absence_data)/2)+1):nrow(absence_data)])
absence_data$mean_precip[(floor(nrow(absence_data)/2)+1):nrow(absence_data)] <- raster::extract(
  climate_pulse_2[["mean_precip"]], absence_data$cell_id[(floor(nrow(absence_data)/2)+1):nrow(absence_data)])

# Ten-fold cross validation
set.seed(7)
# Group training data
bioclim_training_mean$group <- kfold(nrow(bioclim_training_mean), 5)
# Prepare output lists
bioclim_models <- list()
bioclim_evals <- list()
bioclim_threshs <- list() 

# Carry out cross validation
for(i in 1:5){
  # Fit Bioclim envlope model to training data
  bioclim_models[i] <- bioclim(bioclim_training_mean[bioclim_training_mean$group != i, 1:2])

  # Evaluate model
  bioclim_evals[i] <- evaluate(p = bioclim_training_mean[bioclim_training_mean$group == i, 1:2],
                           a = absence_data[,4:5],
                           model = bioclim_models[[i]])
  # Determine threshold
  bioclim_threshs[i] <- threshold(bioclim_evals[[i]], 'spec_sens')
}
mean(unlist(lapply(bioclim_evals, function(x) x@auc))) 

# fit single model 
bioclim_model <- bioclim(bioclim_training_mean[,1:2])
bioclim_eval <- evaluate(bioclim_training_mean[,1:2],
                         absence_data[,4:5],
                         bioclim_model)
bioclim_thresh <-  threshold(bioclim_eval, 'spec_sens')

# Predict for both time periods
predictions_mean <- predict(climate_mean,
                            bioclim_model)
predictions_pulse_1 <- predict(climate_pulse_1,
                                   bioclim_model) 
predictions_pulse_2 <- predict(climate_pulse_2, 
                                   bioclim_model)

# Apply threshold
predictions_mean_thresh <- predictions_mean > bioclim_thresh
predictions_pulse_1_thresh <- predictions_pulse_1 > bioclim_thresh
predictions_pulse_2_thresh <- predictions_pulse_2 > bioclim_thresh

# Visualise predictions

# Define function to visualise predictions
plot_preds_map <- function(base_raster,
                           main_title,
                           file_name,
                           palette = "default" ) {
  if(palette[1] == "default"){
    palette <- sequential_hcl(
      100, 
      "Inferno")
  }
  png(file_name, 
      width = 6,
      height = 3,
      units = "in",
      res = 300)
  print(levelplot(base_raster, 
                  margin = F,
                  main = main_title,
                  colorkey = F,
                  par.settings = rasterTheme(
                    region = palette)) +
          latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 1)) +
          latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
          latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
          latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "uncertain",],
                                        col = "grey90", pch = 63, cex = 0.8)) +
          latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_1",],
                                        col = "dodgerblue2", pch = 1, cex = 1.05)) +
          latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_2",],
                                        col = "black", pch = 3)) +
          
          latticeExtra::layer({
            centre_x <- 27.5
            centre_y <- 58.3
            grid.rect(x=centre_x, y=centre_y,
                      width=7.5, height=2.5,
                      gp=gpar(fill="white",
                              col = "black",
                              alpha = 0.9),
                      default.units='native')
            
            
            grid.points(#label = "Pulse 1",
              x=centre_x-2.75, y=centre_y+2.5/4,
              pch = 1,
              gp=gpar(cex = 0.65,
                      col = "dodgerblue2"),
              default.units='native')
            grid.points(#label = "Pulse 2",
              x=centre_x-2.75, y=centre_y,
              pch = 3,
              gp=gpar(cex = 0.4,
                      col = "black"),
              default.units='native')
            grid.points(#label = "Uncertain",
              x=centre_x-2.75, y=centre_y-2.5/4,
              pch = 63,
              gp=gpar(cex = 0.65,
                      col = "grey30"),
              default.units='native')
            
            grid.text(label = "Pulse 1",
                      x=centre_x-1.75, y=centre_y+2.5/4,
                      #just = "left",
                      hjust = 0,
                      vjust = 0.4,
                      gp=gpar(cex=0.5,
                              col = "black"),
                      default.units='native')
            grid.text(label = "Pusle 2",
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
          }))
  dev.off()
  return("Done")
}

# |_ Mean ----
plot_preds_map(base_raster = predictions_mean,
               main = "Raw Predictions 15.0k-13.8k BP ", 
               file_name = "figures/predictions_mean1.png")
plot_preds_map(base_raster = predictions_mean_thresh,
               main = "Presence/Absence Predictions 15.0k-13.8k BP ", 
               file_name = "figures/predictions_mean_thresh.png",
               palette = c("grey30", "#f38f32"))

# |_ Pulse 1----
plot_preds_map(base_raster = predictions_pulse_1,
               main = "Raw Predictions 14.5k-14.3k BP ", 
               file_name = "figures/predictions_pulse_1.png")
plot_preds_map(base_raster = predictions_pulse_1_thresh,
               main = "Presence/Absence  Predictions 14.5k-14.3k BP ", 
               file_name = "figures/predictions_pulse_1_thresh.png",
               palette = c("grey30", "#f38f32"))

# |_ Pulse 2 ----
plot_preds_map(base_raster = predictions_pulse_2,
               main = "Raw Predictions 14.2k-14.1k BP ", 
               file_name = "figures/predictions_pulse_2.png")
plot_preds_map(base_raster = predictions_pulse_2_thresh,
               main = "Presence/Absence  Predictions 14.2k-14.1k BP ", 
               file_name = "figures/predictions_pulse_2_thresh.png",
               palette = c("grey30", "#f38f32"))

# |_ Predictions Time-Series ----

# Set range
time_range <- seq(time_range_bp$max_cent, time_range_bp$min_cent, 1)

# Crop and stack rasters
time_series_climate <- lapply(time_range, function(cent_bce){
  cat(paste0(cent_bce, " ... \n"))
  temp_raster_name <- paste0("temp_", cent_bce, "_BCE")
  precip_raster_name <- paste0("precip_", cent_bce, "_BCE")
  climate_bce <- stack(crop(temp[[temp_raster_name]],area_of_interest),
                       crop(precip[[precip_raster_name]], area_of_interest))
}) %>% setNames(time_range)

# Calculate predictions
preds_time_series <- lapply(time_series_climate, function(climate_raster){
  cat(paste0(names(climate_raster), " ... \n"))
  
  # update climate raster names to match model
  climate_raster <- setNames(climate_raster, c("mean_temp", "mean_precip"))
  # Predict for both time periods
  predictions <- predict(climate_raster,
                                 bioclim_model)
  # Apply threshold
  predictions <- predictions > bioclim_thresh
  
  return(predictions)
}) %>% setNames(time_range)

# Generate Plots
lapply(seq_along(time_range), function(index){
  k_year_bp <- time_range[index] / 10 + 2
  cat(paste0(k_year_bp, "k BP...\n"))
  plot_preds_map(base_raster = preds_time_series[[index]],
                 main = paste0("Predictions ", formatC(k_year_bp, 1, format = "f"), "k BP"), 
                 file_name = paste0("figures/preds_time-series/predictions_", 
                                    formatC(k_year_bp, 1, format = "f"), "k_BP.png"),
                 palette = c("grey30", "#f38f32"))
  return("Done")
})

# |_ Generate Animation ----
library("magick")
# List files and loade as images
imgs <- list.files("figures/preds_time-series/", 
                   pattern = ".png", full.names = TRUE)
img_list <- lapply(rev(imgs), image_read)

# join images 
img_joined <- image_join(img_list)

# animate at 1 frame per second
img_animated <- image_animate(img_joined, fps = 0.5)

# Export as gif
image_write(image = img_animated,
            path = "figures/preds_time-series/preds_time_series.gif")



# GEOFACETS playground (experimental) ----

# Create Grid
x_grid <- seq(floor(min(hamburgian_sites$Longitude)),
              ceiling(max(hamburgian_sites$Longitude)),
              1)
y_grid <- seq(floor(min(hamburgian_sites$Latitude)),
              ceiling(max(hamburgian_sites$Latitude)),
              1)
hamburgian_sites <- hamburgian_sites %>%
  arrange(Longitude, Latitude)
hamburgian_sites$col <- ceiling(hamburgian_sites$Longitude*10) - 4 *10
hamburgian_sites$row <- ceiling(hamburgian_sites$Latitude*10) - 55 *10
ggplot(hamburgian_sites, aes(x = row)) + geom_histogram(binwidth =  1)
ggplot(hamburgian_sites, aes(x = col)) + geom_histogram(binwidth =  1)
hamburgian_sites %>% group_by(col) %>% summarise(n = n()) %>% filter(n> 1)
hamburgian_sites2 <- hamburgian_sites %>%
  st_drop_geometry() %>%
  mutate(col2 = col, row2 = row) %>%
  group_by(col) %>%
  group_modify(function(x, y){
    print(x$row)
    print(y)
    print(x$row2)
    if(nrow(x) > 1) {
      x$col2 <- x$col2 + rep(seq(0, nrow(x)/10, nrow(x)/10/nrow(x)), ceiling(nrow(x)/2))[1:nrow(x)] 
      x$row2 <- x$row2 + c(rep(0, ceiling(nrow(x)/2)),
                           rep(1, ceiling(nrow(x)/2)))[1:nrow(x)] 
      return(x)
    } else return(x)
  }) %>% 
  ungroup() %>%
  arrange(col2,row2) %>%
  mutate(col = as.numeric(factor(col2)),
         row = as.numeric(factor(row2)))

library(geofacet)
grid_design(data = dplyr::select(hamburgian_sites, Site, Latitude, Longitude, col, row),
            img = "http://bit.ly/us-grid")

write.csv(select(hamburgian_sites2,
                 Site,
                 Latitude,
                 Longitude,
                 col,
                 row), 
          "test_grid.csv", row.names = F)
