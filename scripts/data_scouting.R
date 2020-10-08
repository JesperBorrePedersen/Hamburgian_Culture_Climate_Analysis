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
hamburgian_sites <- read.csv("data/HamburgianData.csv") %>%
  select(1:20) %>%
  mutate(classic_ham = X14.520.Ð.14.100.BP,
         havelte_grp = X14.750.Ð.14.470.BP,
         classic_ham_start = 14520,
         classic_ham_end = 14100,
         havelte_grp_start = 14750,
         havelte_grp_end = 14470) %>%
  select(-X14.520.Ð.14.100.BP,
         -X14.750.Ð.14.470.BP)


# Covnert to sf object
hamburgian_sites <- st_as_sf(hamburgian_sites, 
                             coords = c("Longitude", "Latitude"),
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
temp_df$Site <- hamburgian_sites$Site
precip_df$Site <- hamburgian_sites$Site

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
  group_by(Site) %>%
  filter(year_BP >= 13800 & year_BP <= 15000) %>%
  summarise(mean_temp = mean(temp))
mean_temp_classic <- temp_df %>% 
  group_by(Site) %>%
  filter(year_BP >= 14100 & year_BP <= 14520) %>%
  summarise(mean_temp_classic = mean(temp))
mean_temp_havelte <- temp_df %>% 
  group_by(Site) %>%
  filter(year_BP >= 14470 & year_BP <= 14750) %>%
  summarise(mean_temp_havelte = mean(temp))

mean_precip <- precip_df %>% 
  group_by(Site) %>%
  filter(year_BP >= 13800 & year_BP <= 15000) %>%
  summarise(mean_precip = mean(precip))
mean_precip_classic <- precip_df %>% 
  group_by(Site) %>%
  filter(year_BP >= 14100 & year_BP <= 14520) %>%
  summarise(mean_precip_classic = mean(precip))
mean_precip_havelte <- precip_df %>% 
  group_by(Site) %>%
  filter(year_BP >= 14470 & year_BP <= 14750) %>%
  summarise(mean_precip_havelte = mean(precip))

# merge with hamburgian_sites df
hamburgian_sites <- hamburgian_sites %>%
  full_join(mean_temp) %>%
  full_join(mean_temp_classic) %>%
  full_join(mean_temp_havelte) %>%
  full_join(mean_precip) %>%
  full_join(mean_precip_classic) %>%
  full_join(mean_precip_havelte) 

# tidy up
rm(mean_temp)
rm(mean_temp_classic)
rm(mean_temp_havelte)
rm(mean_precip)
rm(mean_precip_classic)
rm(mean_precip_havelte)


# set mean temp to NA where no evidence is avialable
hamburgian_sites$mean_temp_classic[hamburgian_sites$classic_ham == 0] <- NA
hamburgian_sites$mean_precip_classic[hamburgian_sites$classic_ham == 0] <- NA
hamburgian_sites$mean_temp_halvelte[hamburgian_sites$havelte_grp == 0] <- NA
hamburgian_sites$mean_precip_halvelte[hamburgian_sites$havelte_grp == 0] <- NA

hamburgian_sites$classic_ham_start[hamburgian_sites$classic_ham == 0] <- NA
hamburgian_sites$classic_ham_end[hamburgian_sites$classic_ham == 0] <- NA
hamburgian_sites$havelte_grp_start[hamburgian_sites$havelte_grp == 0] <- NA
hamburgian_sites$havelte_grp_end[hamburgian_sites$havelte_grp == 0] <- NA

# Order sites by latitude
hamburgian_sites <- hamburgian_sites %>%
  mutate(Site = ordered(Site, levels = levels(fct_reorder(hamburgian_sites$Site, hamburgian_sites$Latitude))))

# |_ Temp Time-Series Plot ----
temp_plot <- ggplot(temp_df %>% filter(year_BP >= 13000 & year_BP <= 15500) %>%
                      mutate(Site = ordered(Site, levels = levels(fct_reorder(hamburgian_sites$Site, hamburgian_sites$Latitude)))),
       aes(x = year_BP,
           y = temp,
           colour = Site)) +
  geom_line() +
  geom_segment(data = hamburgian_sites,
             mapping = aes(x = classic_ham_end,
                           xend = classic_ham_start,
                           y = mean_temp_classic, 
                           yend = mean_temp_classic),
             colour = "black") +
  geom_segment(data = hamburgian_sites,
    mapping = aes(x = havelte_grp_end,
                  xend = havelte_grp_start,
                  y = mean_temp_havelte, 
                  yend = mean_temp_havelte),
    colour = "black") +
  geom_text(data = hamburgian_sites,
            mapping = aes(x = havelte_grp_start + 100,
                             y = mean_temp_havelte, 
                             label = paste0("mean:\n", round(mean_temp_havelte,1), "°C")),
               colour = "black",
            vjust = 0.5,
            hjust = 1) +
  geom_text(data = hamburgian_sites,
            mapping = aes(x = classic_ham_end - 100, 
                          y = mean_temp_classic,
                          label = paste0("mean:\n", round(mean_temp_classic,1), "°C")),
            vjust = 0.5,
            hjust = 0,
            colour = "black") +
  labs(x = "Year BP", y = "Annual Mean Temp °C (Bio01)") +
  scale_x_reverse(limits = c(15500, 13000),
                  breaks = seq(15500, 13000, -500),
                  labels = rev(c("13.0k", "13.5k", "14.0k", "14.5k", "15.0k", "15.5k"))) +
  scale_y_continuous(limits = c(-7, 10),
                     breaks = seq(-6,9,3)) +
  scale_colour_discrete_qualitative(palette = "Dark2") +
  theme_cowplot(15) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

save_plot("figures/figure1_temp.png",
          temp_plot,
          base_aspect_ratio = 1.6,
          base_height = 4)

# Plot per site
temp_plot <- temp_plot + facet_wrap(vars(Site),
                                    scales = "free") +
  theme(axis.title = element_text(face = "bold"))
save_plot("figures/figure1_temp_by_site.png",
          temp_plot,
          base_asp = 1.6,
          base_height = 20)

# |_ Precip Time-Series Plot ----
precip_plot <- ggplot(precip_df %>% filter(year_BP >= 13000 & year_BP <= 15500) %>%
                      mutate(Site = ordered(Site, levels = levels(fct_reorder(hamburgian_sites$Site, hamburgian_sites$Latitude)))),
                    aes(x = year_BP,
                        y = precip,
                        colour = Site)) +
  geom_line() +
  geom_segment(data = hamburgian_sites,
               mapping = aes(x = classic_ham_end,
                             xend = classic_ham_start,
                             y = mean_precip_classic, 
                             yend = mean_precip_classic),
               colour = "black") +
  geom_segment(data = hamburgian_sites,
               mapping = aes(x = havelte_grp_end,
                             xend = havelte_grp_start,
                             y = mean_precip_havelte, 
                             yend = mean_precip_havelte),
               colour = "black") +
  geom_text(data = hamburgian_sites,
            mapping = aes(x = havelte_grp_start + 100,
                          y = mean_precip_havelte, 
                          label = paste0("mean:\n", round(mean_precip_havelte,1), "°C")),
            colour = "black",
            vjust = 0.5,
            hjust = 1) +
  geom_text(data = hamburgian_sites,
            mapping = aes(x = classic_ham_end - 100, 
                          y = mean_precip_classic,
                          label = paste0("mean:\n", round(mean_precip_classic,1), "°C")),
            vjust = 0.5,
            hjust = 0,
            colour = "black") +
  labs(x = "Year BP", y = "Annual Precipitation mm (Bio12)") +
  scale_x_reverse(limits = c(15500, 13000),
                  breaks = seq(15500, 13000, -500),
                  labels = rev(c("13.0k", "13.5k", "14.0k", "14.5k", "15.0k", "15.5k"))) +
  scale_y_continuous(limits = c(400, 1100),
                     breaks = seq(400,1200,200)) +
  scale_colour_discrete_qualitative(palette = "Dark2") +
  theme_cowplot(15) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = NA),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

save_plot("figures/figure1_precip.png",
          precip_plot,
          base_aspect_ratio = 1.6,
          base_height = 4)

# Plot per site
precip_plot <- precip_plot + facet_wrap(vars(Site),
                                    scales = "free") +
  theme(axis.title = element_text(face = "bold"))
save_plot("figures/figure1_precip_by_site.png",
          precip_plot,
          base_asp = 1.6,
          base_height = 20)

## 3) Raster Maps ----
area_of_interest <- extent(c(-11,30,50,60)) 
time_range_bp <- data.frame(min_bp = 15000-2000,
                         max_bp = 13800-2000,
                         min_cent = 150-20,
                         max_cent = 138-20,
                         classic_ham_start = 14520-2000,
                         classic_ham_end = 14100-2000,
                         havelte_grp_start = 14750-2000,
                         havelte_grp_end = 14470-2000,
                         classic_ham_min_cent = 145-20,
                         classic_ham_max_cent = 141-20,
                         havelte_grp_min_cent = 147-20,
                         havelte_grp_max_cent = 145-20)

## calculate average values for the broad time-range
# Temperature
temp_time_range_bp_all_cents <- paste0(
  "temp_",
  seq(time_range_bp$max_cent, time_range_bp$min_cent, 1),
  "_BCE")
hamburgian_mean_temp <- mean(crop(temp[[temp_time_range_bp_all_cents]], area_of_interest))

temp_time_range_bp_classic_ham <- paste0(
  "temp_",
  seq(time_range_bp$classic_ham_max_cent, time_range_bp$classic_ham_min_cent, 1),
  "_BCE")
hamburgian_mean_temp_classic_ham <- mean(crop(temp[[temp_time_range_bp_classic_ham]], area_of_interest))

temp_time_range_bp_havelte_grp <- paste0(
  "temp_",
  seq(time_range_bp$havelte_grp_max_cent, time_range_bp$havelte_grp_min_cent, 1),
  "_BCE")
hamburgian_mean_temp_havelte_grp <- mean(crop(temp[[temp_time_range_bp_havelte_grp]], area_of_interest))

# Precip
precip_time_range_bp_all_cents <- paste0(
  "precip_",
  seq(time_range_bp$max_cent, time_range_bp$min_cent, 1),
  "_BCE")
hamburgian_mean_precip <- mean(crop(precip[[precip_time_range_bp_all_cents]], area_of_interest))

precip_time_range_bp_classic_ham <- paste0(
  "precip_",
  seq(time_range_bp$classic_ham_max_cent, time_range_bp$classic_ham_min_cent, 1),
  "_BCE")
hamburgian_mean_precip_classic_ham <- mean(crop(precip[[precip_time_range_bp_classic_ham]], area_of_interest))

precip_time_range_bp_havelte_grp <- paste0(
  "precip_",
  seq(time_range_bp$havelte_grp_max_cent, time_range_bp$havelte_grp_min_cent, 1),
  "_BCE")
hamburgian_mean_precip_havelte_grp <- mean(crop(precip[[precip_time_range_bp_havelte_grp]], area_of_interest))

## Save rasters for future use
writeRaster(hamburgian_mean_temp, 
            "data/mean_rasters/hamburgian_mean_temp.tif",
            overwrite = T)
writeRaster(hamburgian_mean_temp_classic_ham, 
            "data/mean_rasters/hamburgian_mean_temp_classic_ham.tif",
            overwrite = T)
writeRaster(hamburgian_mean_temp_havelte_grp, 
            "data/mean_rasters/hamburgian_mean_temp_havelte_grp.tif",
            overwrite = T)

writeRaster(hamburgian_mean_precip, 
            "data/mean_rasters/hamburgian_mean_precip.tif",
            overwrite = T)
writeRaster(hamburgian_mean_precip_classic_ham, 
            "data/mean_rasters/hamburgian_mean_precip_classic_ham.tif",
            overwrite = T)
writeRaster(hamburgian_mean_precip_havelte_grp, 
            "data/mean_rasters/hamburgian_mean_precip_havelte_grp.tif",
            overwrite = T)


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
# Mean preciptation map 15k-13.5k BP
png("figures/mean_temp.png", 
    width = 6,
    height = 3,
    units = "in",
    res = 300)
levelplot(hamburgian_mean_temp, 
          margin = F,
          main = "Mean Temperature 15k-13.5k BP (°C)",
          at = seq(-12, 12, 1)) +
  latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 0.75)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Classic Hamburgian",],
                                col = "dodgerblue2", pch = 1, cex = 1.05)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Havelte Group",],
                                col = "black", pch = 3)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Possibly Hamburgian",],
                                col = "gray90", pch = 63, cex = 1.05)) +
  latticeExtra::layer({
    centre_x <- 25.5
    centre_y <- 51.5
    grid.rect(x=centre_x, y=centre_y,
              width=7.5, height=2.5,
              gp=gpar(fill="white", 
                      col = "black",
                      alpha = 0.75),
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
    
    grid.text(label = "Hamburgian",
              x=centre_x-1.75, y=centre_y+2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Havelte Group",
              x=centre_x-1.75, y=centre_y,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Possibly Ham.",
              x=centre_x-1.75, y=centre_y-2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
  })
dev.off()


# Mean temperature map Classic Hamburgian 14.5k-14.1k BP
png("figures/mean_temp_classic_ham.png", 
    width = 6,
    height = 3,
    units = "in",
    res = 300)
levelplot(hamburgian_mean_temp_classic_ham, 
          margin = F,
          main = "Mean Temperature 14.5k-14.1k BP (°C)",
          at = seq(-12, 12, 1)) +
  latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 0.75)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Classic Hamburgian",],
                                col = "dodgerblue2", pch = 1, cex = 1.05)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Havelte Group",],
                                col = "black", pch = 3)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Possibly Hamburgian",],
                                col = "gray90", pch = 63, cex = 1.05)) +
  latticeExtra::layer({
    centre_x <- 25.5
    centre_y <- 51.5
    grid.rect(x=centre_x, y=centre_y,
              width=7.5, height=2.5,
              gp=gpar(fill="white", 
                      col = "black",
                      alpha = 0.75),
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
    
    grid.text(label = "Hamburgian",
              x=centre_x-1.75, y=centre_y+2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Havelte Group",
              x=centre_x-1.75, y=centre_y,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Possibly Ham.",
              x=centre_x-1.75, y=centre_y-2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
  })
dev.off()


# Mean temperature map Havelte 14.7k-14.5k BP
png("figures/mean_temp_havelte_grp.png", 
    width = 6,
    height = 3,
    units = "in",
    res = 300)
levelplot(hamburgian_mean_temp_havelte_grp, 
          margin = F,
          main = "Mean Temperature 14.7k-14.5k BP (°C)",
          at = seq(-12, 12, 1)) +
  latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 0.75)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Classic Hamburgian",],
                                col = "dodgerblue2", pch = 1, cex = 1.05)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Havelte Group",],
                                col = "black", pch = 3)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Possibly Hamburgian",],
                                col = "gray90", pch = 63, cex = 1.05)) +
  latticeExtra::layer({
    centre_x <- 25.5
    centre_y <- 51.5
    grid.rect(x=centre_x, y=centre_y,
              width=7.5, height=2.5,
              gp=gpar(fill="white", 
                      col = "black",
                      alpha = 0.75),
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
    
    grid.text(label = "Hamburgian",
              x=centre_x-1.75, y=centre_y+2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Havelte Group",
              x=centre_x-1.75, y=centre_y,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Possibly Ham.",
              x=centre_x-1.75, y=centre_y-2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
  })
dev.off()

# |_ Precip Maps ----

# Mean preciptation map 15k-13.5k BP
png("figures/mean_precip.png", 
    width = 6,
    height = 3,
    units = "in",
    res = 300)
levelplot(hamburgian_mean_precip, 
          margin = F,
          main = "Mean Precipitation 15k-13.5k BP (mm)",
          at = seq(0, 2500, 2500/20),
          par.settings = rasterTheme(
            region = sequential_hcl(
              100, 
              "Blues",
              rev = T))) + 
  latticeExtra::layer(sp.polygons(land_for_maps, col = "white", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 0.75)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Classic Hamburgian",],
                                col = "dodgerblue2", pch = 1, cex = 1.05)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Havelte Group",],
                                col = "black", pch = 3)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Possibly Hamburgian",],
                                col = "grey30", pch = 63, cex = 1.05)) +
  latticeExtra::layer({
    centre_x <- 25.5
    centre_y <- 51.5
    grid.rect(x=centre_x, y=centre_y,
              width=7.5, height=2.5,
              gp=gpar(fill="white", 
                      col = "black",
                      alpha = 0.75),
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
    
    grid.text(label = "Hamburgian",
              x=centre_x-1.75, y=centre_y+2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Havelte Group",
              x=centre_x-1.75, y=centre_y,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Possibly Ham.",
              x=centre_x-1.75, y=centre_y-2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
  })
dev.off()

# Mean preciptation map 14.5k-14.1k BP
png("figures/mean_precip_classic_ham.png", 
    width = 6,
    height = 3,
    units = "in",
    res = 300)
levelplot(hamburgian_mean_precip_classic_ham, 
          margin = F,
          main = "Mean Precipitation 14.5k-14.1k BP (mm)",
          at = seq(0, 2500, 2500/20),
          par.settings = rasterTheme(
            region = sequential_hcl(
              100, 
              "Blues",
              rev = T))) + 
  latticeExtra::layer(sp.polygons(land_for_maps, col = "white", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 0.75)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Classic Hamburgian",],
                                col = "dodgerblue2", pch = 1, cex = 1.05)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Havelte Group",],
                                col = "black", pch = 3)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Possibly Hamburgian",],
                                col = "grey30", pch = 63, cex = 1.05)) +
  latticeExtra::layer({
    centre_x <- 25.5
    centre_y <- 51.5
    grid.rect(x=centre_x, y=centre_y,
              width=7.5, height=2.5,
              gp=gpar(fill="white", 
                      col = "black",
                      alpha = 0.75),
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
    
    grid.text(label = "Hamburgian",
              x=centre_x-1.75, y=centre_y+2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Havelte Group",
              x=centre_x-1.75, y=centre_y,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Possibly Ham.",
              x=centre_x-1.75, y=centre_y-2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
  })
dev.off()


# Mean preciptation map 14.7k-14.5k BP
png("figures/mean_precip_havelte_grp.png", 
    width = 6,
    height = 3,
    units = "in",
    res = 300)
levelplot(hamburgian_mean_precip_havelte_grp, 
          margin = F,
          main = "Mean Precipitation 14.7k-14.5k BP (mm)",
          at = seq(0, 2500, 2500/20),
          par.settings = rasterTheme(
            region = sequential_hcl(
              100, 
              "Blues",
              rev = T))) + 
  latticeExtra::layer(sp.polygons(land_for_maps, col = "white", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 0.5)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 0.75)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Classic Hamburgian",],
                                col = "dodgerblue2", pch = 1, cex = 1.05)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Havelte Group",],
                                col = "black", pch = 3)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Possibly Hamburgian",],
                                col = "grey30", pch = 63, cex = 1.05)) +
  latticeExtra::layer({
    centre_x <- 25.5
    centre_y <- 51.5
    grid.rect(x=centre_x, y=centre_y,
              width=7.5, height=2.5,
              gp=gpar(fill="white", 
                      col = "black",
                      alpha = 0.75),
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
    
    grid.text(label = "Hamburgian",
              x=centre_x-1.75, y=centre_y+2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Havelte Group",
              x=centre_x-1.75, y=centre_y,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Possibly Ham.",
              x=centre_x-1.75, y=centre_y-2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
  })
dev.off()

# 3) Temp vs. Precip Plot ---- 

# Cross check with extraction from rasters
mean_both_ham_extract <- data.frame(
  Site = unique(hamburgian_sites$Site),
  mean_temp = raster::extract(hamburgian_mean_temp, as_Spatial(hamburgian_sites)),
  mean_precip = raster::extract(hamburgian_mean_precip, as_Spatial(hamburgian_sites)),
  mean_temp_classic = raster::extract(hamburgian_mean_temp_classic_ham, as_Spatial(hamburgian_sites)),
  mean_precip_classic = raster::extract(hamburgian_mean_precip_classic_ham, as_Spatial(hamburgian_sites)),
  mean_temp_havelte = raster::extract(hamburgian_mean_temp_havelte_grp, as_Spatial(hamburgian_sites)),
  mean_precip_havelte = raster::extract(hamburgian_mean_precip_havelte_grp, as_Spatial(hamburgian_sites))
) %>% 
  as_tibble() %>%
  mutate(Site = ordered(Site, levels = levels(fct_reorder(hamburgian_sites$Site, hamburgian_sites$Latitude))))
all_equal(hamburgian_sites %>%
            st_drop_geometry() %>% 
            select(Site,
                 mean_temp,
                 mean_precip,
                 mean_temp_classic,
                 mean_temp_havelte,
                 mean_precip_classic,
                 mean_precip_havelte),
          mean_both_ham_extract)

set.seed(6)

# Mask rasters
hamburgian_mean_temp_classic_ham_masked <- mask(
  hamburgian_mean_temp_classic_ham,
  land_for_maps)
hamburgian_mean_precip_classic_ham_masked <- mask(
  hamburgian_mean_precip_classic_ham,
  land_for_maps)

hamburgian_mean_temp_havelte_grp_masked <- mask(
  hamburgian_mean_temp_havelte_grp,
  land_for_maps)
hamburgian_mean_precip_havelte_grp_masked <- mask(
  hamburgian_mean_precip_havelte_grp,
  land_for_maps)

# Optain random sample for land areas from masked rasters
random_sample <- data.frame(
  site = paste0("Random_", 1:100000),
  mean_temp_classic = sampleRandom(hamburgian_mean_temp_classic_ham_masked, 100000),
  mean_precip_classic = sampleRandom(hamburgian_mean_precip_classic_ham_masked, 100000),
  mean_temp_havelte = sampleRandom(hamburgian_mean_temp_havelte_grp_masked, 100000),
  mean_precip_havelte = sampleRandom(hamburgian_mean_precip_havelte_grp_masked, 100000))

temp_vs_precip_plot <- ggplot() +  
  geom_point(data = random_sample,
             aes(x = mean_temp_classic,
                 y = mean_precip_classic),
             colour = "lightskyblue2",
             fill = "lightskyblue2",
             alpha = 0.25,
             shape = 21
             ) +
  geom_point(data = random_sample,
             aes(x = mean_temp_havelte,
                 y = mean_precip_havelte),
             colour = "grey",
             fill = "grey",
             alpha = 0.25,
             shape = 21
             ) +
  geom_point(data = hamburgian_sites,
             aes(x = mean_temp_classic,
                 y = mean_precip_classic),
             colour = "dodgerblue2",
             shape = 16
             ) +
  geom_point(data = hamburgian_sites,
             aes(x = mean_temp_havelte,
                 y = mean_precip_havelte),
             colour = "black",
             shape = 16
             ) +
  labs(x = "Mean Temperature (°C)",
       y = "Mean Precipitation (mm)") +
  scale_x_continuous(limits = c(-7, 12),
                     breaks = seq(-6, 12, 2)) +
  scale_y_continuous(limits = c(400, 2500)) +
  annotate("text", x = 7, y = 2500, 
           colour = "dodgerblue2", hjust = 0, vjust = 0.4,
           label = "Classic Hamburgian") +
  annotate("text", x = 7, y = 2400,
           colour = "black", hjust = 0, vjust = 0.4,
           label = "Havelte Group") +
  annotate("point", x = 6.7, y = 2500, 
           colour = "dodgerblue2") +
  annotate("point", x = 6.7, y = 2400,
           colour = "black") +
  theme_cowplot(14) 
save_plot("figures/temp_vs_precip.png",
          temp_vs_precip_plot,
          base_aspect_ratio = 1.6,
          base_height = 4)

# 4) BIOCLIM Model ----

# Add cell_id to main geometry to allow for exclusion of multiple samples from
# the a single cell
hamburgian_sites <- raster::extract(hamburgian_mean_temp, 
                as_Spatial(hamburgian_sites),
                cellnumbers = T,
                df = T) %>% 
  mutate(Site = hamburgian_sites$Site,
         cell_id = cells) %>%
  select(Site, cell_id) %>%
  full_join(hamburgian_sites, .)

# Gather temperature and precip columns into training data
bioclim_training_all <- hamburgian_sites %>%
  st_drop_geometry() %>%
  filter(Classic.Hamburgian.Havelte.Group != "Possibly Hamburgian") %>%
  select(Site, 
         Latitude,
         Longitude,
         cell_id,
         mean_temp_classic,
         mean_temp_havelte) %>%
  gather("temp_var", "mean_temp", mean_temp_classic, mean_temp_havelte) %>%
  na.omit() %>%
  distinct(cell_id, temp_var, mean_temp) %>%
  mutate(period = gsub("mean_temp_(.*)", "\\1", temp_var))
bioclim_training_all <- hamburgian_sites %>%
  st_drop_geometry() %>%
  filter(Classic.Hamburgian.Havelte.Group != "Possibly Hamburgian") %>%
  select(Site, 
         Latitude,
         Longitude,
         cell_id,
         mean_precip_classic,
         mean_precip_havelte) %>%
  gather("precip_var", "mean_precip", mean_precip_classic, mean_precip_havelte) %>%
  na.omit() %>%
  distinct(cell_id, precip_var, mean_precip) %>%
  mutate(period = gsub("mean_precip_(.*)", "\\1", precip_var)) %>%
  full_join(bioclim_training_all, .)

# Quality control
view(bioclim_training_all)

# Retain climate columns only
bioclim_training <- select(bioclim_training_all, "mean_temp", "mean_precip")

# Split into training and validation data
set.seed(0)
sample_rows <- sample(1:nrow(bioclim_training),
                      floor(nrow(bioclim_training) / 3)) %>%
  sort()
bioclim_validation <- bioclim_training[sample_rows,]
bioclim_training <- bioclim_training[-sample_rows,]

# Fit Bioclim envlope model to training data
bioclim_model <- bioclim(bioclim_training)

# Prepare raster stacks for predictions
climate_classic_ham <- stack(hamburgian_mean_temp_classic_ham,
                             hamburgian_mean_precip_classic_ham)
climate_classic_ham <- setNames(climate_classic_ham, 
                                c("mean_temp", "mean_precip"))
climate_havelte_grp <- stack(hamburgian_mean_temp_havelte_grp,
                             hamburgian_mean_precip_havelte_grp)
climate_havelte_grp <- setNames(climate_havelte_grp,
                                c("mean_temp", "mean_precip"))

# Mask out everything but land
climate_classic_ham_masked <- mask(climate_classic_ham,
                                   land_for_maps)
climate_havelte_grp_masked <- mask(climate_havelte_grp,
                                   land_for_maps)

# Generate random background points for evaluation (absence data)
set.seed(50)
absence_data <- randomPoints(climate_classic_ham_masked,
                             200)
colnames(absence_data) <- c('lon', 'lat')
absence_data <- as.data.frame(absence_data)

# Add cell_ids
absence_data$cell_id <- st_as_sf(absence_data,
                                 coords = c("lon", "lat"),
                                 crs = st_crs(climate_classic_ham)) %>%
  raster::extract(climate_classic_ham,
                  .,
                  df = T,
                  cellnumbers = T) %>% 
  pull(cells)

# filter cell_ids to avoid coccurence with presence data
absence_data <- absence_data[!(absence_data$cell_id %in% bioclim_training_all$cell_id),]

# fill in absence data, half from classic, half from halvelte
absence_data$mean_temp[1:floor(nrow(absence_data)/2)] <- raster::extract(
  climate_classic_ham[["mean_temp"]], absence_data$cell_id[1:floor(nrow(absence_data)/2)])
absence_data$mean_precip[1:floor(nrow(absence_data)/2)] <- raster::extract(
  climate_classic_ham[["mean_precip"]], absence_data$cell_id[1:floor(nrow(absence_data)/2)])

absence_data$mean_temp[(floor(nrow(absence_data)/2)+1):nrow(absence_data)] <- raster::extract(
  climate_havelte_grp[["mean_temp"]], absence_data$cell_id[(floor(nrow(absence_data)/2)+1):nrow(absence_data)])
absence_data$mean_precip[(floor(nrow(absence_data)/2)+1):nrow(absence_data)] <- raster::extract(
  climate_havelte_grp[["mean_precip"]], absence_data$cell_id[(floor(nrow(absence_data)/2)+1):nrow(absence_data)])


# Evaluate model
bioclim_eval <- evaluate(bioclim_training,
                         absence_data[,c("mean_temp", "mean_precip")],
                         bioclim_model)
# Determine threshold
bioclim_thresh <- threshold(bioclim_eval, 'spec_sens') 

# Predict for both time periods
predictions_classic_ham <- predict(climate_classic_ham,
                                   bioclim_model) 
predictions_havelte_grp <- predict(climate_havelte_grp, 
                                   bioclim_model)

# Apply threshold
predictions_classic_ham_thresh <- predictions_classic_ham > bioclim_thresh
predictions_havelte_grp_thresh <- predictions_havelte_grp > bioclim_thresh

# Visualise predictions

# |_ Maps Classic Ham ----

png("figures/predictions_classic_ham.png", 
    width = 6,
    height = 3,
    units = "in",
    res = 300)
levelplot(predictions_classic_ham, 
          margin = F,
          main = "Raw Predictions 14.5k-14.1k BP ",
          colorkey = F,
          par.settings = rasterTheme(
            region = sequential_hcl(
              100, 
              "Inferno"))) +
  latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Classic Hamburgian",],
                                col = "dodgerblue2", pch = 1, cex = 1.05)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Havelte Group",],
                                col = "white", pch = 3)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Possibly Hamburgian",],
                                col = "gray90", pch = 63, cex = 1.05)) +
  latticeExtra::layer({
    centre_x <- 25.5
    centre_y <- 51.5
    grid.rect(x=centre_x, y=centre_y,
              width=7.5, height=2.5,
              gp=gpar(fill="white", 
                      col = "black",
                      alpha = 0.75),
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
              col = "white"),
      default.units='native')
    grid.points(#label = "Possibly Ham.",
      x=centre_x-2.75, y=centre_y-2.5/4,
      pch = 63,
      gp=gpar(cex = 0.65, 
              col = "grey30"),
      default.units='native')
    
    grid.text(label = "Hamburgian",
              x=centre_x-1.75, y=centre_y+2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Havelte Group",
              x=centre_x-1.75, y=centre_y,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Possibly Ham.",
              x=centre_x-1.75, y=centre_y-2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
  })
dev.off()


png("figures/predictions_classic_ham_thresh.png", 
    width = 6,
    height = 3,
    units = "in",
    res = 300)
levelplot(predictions_classic_ham_thresh, 
          margin = F,
          main = "Presence/Absence Predictions 14.5k-14.1k BP ",
          colorkey = F,
          par.settings = rasterTheme(
            region = sequential_hcl(
              5, 
              "Inferno")[c(1,3)])) +
  latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Classic Hamburgian",],
                                col = "dodgerblue2", pch = 1, cex = 1.05)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Havelte Group",],
                                col = "white", pch = 3)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Possibly Hamburgian",],
                                col = "gray90", pch = 63, cex = 1.05)) +
  latticeExtra::layer({
    centre_x <- 25.5
    centre_y <- 51.5
    grid.rect(x=centre_x, y=centre_y,
              width=7.5, height=2.5,
              gp=gpar(fill="white", 
                      col = "black",
                      alpha = 0.75),
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
              col = "white"),
      default.units='native')
    grid.points(#label = "Possibly Ham.",
      x=centre_x-2.75, y=centre_y-2.5/4,
      pch = 63,
      gp=gpar(cex = 0.65, 
              col = "grey30"),
      default.units='native')
    
    grid.text(label = "Hamburgian",
              x=centre_x-1.75, y=centre_y+2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Havelte Group",
              x=centre_x-1.75, y=centre_y,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Possibly Ham.",
              x=centre_x-1.75, y=centre_y-2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
  })
dev.off()

# |_ Maps Havelte Group ----

png("figures/predictions_havelte_grp.png", 
    width = 6,
    height = 3,
    units = "in",
    res = 300)
levelplot(predictions_havelte_grp, 
          margin = F,
          main = "Raw Predictions 14.7k-14.5k BP",
          colorkey = F,
          par.settings = rasterTheme(
            region = sequential_hcl(
              100, 
              "Inferno"))) +
  latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Classic Hamburgian",],
                                col = "dodgerblue2", pch = 1, cex = 1.05)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Havelte Group",],
                                col = "white", pch = 3)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Possibly Hamburgian",],
                                col = "gray90", pch = 63, cex = 1.05)) +
  latticeExtra::layer({
    centre_x <- 25.5
    centre_y <- 51.5
    grid.rect(x=centre_x, y=centre_y,
              width=7.5, height=2.5,
              gp=gpar(fill="white", 
                      col = "black",
                      alpha = 0.75),
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
              col = "white"),
      default.units='native')
    grid.points(#label = "Possibly Ham.",
      x=centre_x-2.75, y=centre_y-2.5/4,
      pch = 63,
      gp=gpar(cex = 0.65, 
              col = "grey30"),
      default.units='native')
    
    grid.text(label = "Hamburgian",
              x=centre_x-1.75, y=centre_y+2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Havelte Group",
              x=centre_x-1.75, y=centre_y,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Possibly Ham.",
              x=centre_x-1.75, y=centre_y-2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
  })
dev.off()

png("figures/predictions_havelte_grp_thresh.png", 
    width = 6,
    height = 3,
    units = "in",
    res = 300)
levelplot(predictions_havelte_grp_thresh, 
          margin = F,
          main = "Presence/Absence Predictions 14.7k-14.5k BP",
          colorkey = F,
          par.settings = rasterTheme(
            region = sequential_hcl(
              5, 
              "Inferno")[c(1,3)])) +
  latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Classic Hamburgian",],
                                col = "dodgerblue2", pch = 1, cex = 1.05)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Havelte Group",],
                                col = "white", pch = 3)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$Classic.Hamburgian.Havelte.Group == "Possibly Hamburgian",],
                                col = "gray90", pch = 63, cex = 1.05)) +
  latticeExtra::layer({
    centre_x <- 25.5
    centre_y <- 51.5
    grid.rect(x=centre_x, y=centre_y,
              width=7.5, height=2.5,
              gp=gpar(fill="white", 
                      col = "black",
                      alpha = 0.75),
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
              col = "white"),
      default.units='native')
    grid.points(#label = "Possibly Ham.",
      x=centre_x-2.75, y=centre_y-2.5/4,
      pch = 63,
      gp=gpar(cex = 0.65, 
              col = "grey30"),
      default.units='native')
    
    grid.text(label = "Hamburgian",
              x=centre_x-1.75, y=centre_y+2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Havelte Group",
              x=centre_x-1.75, y=centre_y,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
    grid.text(label = "Possibly Ham.",
              x=centre_x-1.75, y=centre_y-2.5/4,
              #just = "left",
              hjust = 0,
              vjust = 0.4,
              gp=gpar(cex=0.5, 
                      col = "black"),
              default.units='native')
  })
dev.off()

levelplot(predictions_classic_ham_thresh, margins = F)
levelplot(predictions_havelte_grp, margins = F)
levelplot(predictions_havelte_grp_thresh, margins = F)
plot(climate_classic_ham[[1]])
points(absence_data)
points(rbind(bioclim_training, bioclim_validation))
preds <- predict(stack(hamburgian_mean_temp, hamburgian_mean_precip),
            bioclim_envelope_model)
predictions_classic_ham+ 
  layer(sp.points(as_Spatial(hamburgian_sites), col = "white")) +
  layer(sp.polygons(countries, col = "darkgrey", ))
png("figures/bioclim_env_mean_raw_preds.png", 
    width = 6,
    height = 3,
    units = "in",
    res = 300)
levelplot(preds, 
          margin = F,
          main = "BIOCLIM Envelope Raw Predictions") +
  layer(sp.points(as_Spatial(hamburgian_sites), col = "white")) +
  layer(sp.polygons(countries, col = "darkgrey", ))
dev.off()

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
