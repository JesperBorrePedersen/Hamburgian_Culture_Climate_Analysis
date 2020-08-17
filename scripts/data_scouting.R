# Hamburgian Culture Climate Analysis - First data Exploraiton 
# Jakob Assmann j.assmann@bios.au.dk 17 August 2020

## 1) Housekeeping ----
# Dependencies
library(sf)
library(raster)
library(tidyverse)
library(cowplot)

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
hamburgian_sites <- read.csv("data/Hamburgian Data.csv")
hamburgian_sites <- hamburgian_sites[!is.na(hamburgian_sites$Longitude),]
# Covnert to sf object
hamburgian_sites <- st_as_sf(hamburgian_sites, 
                             coords = c("Longitude", "Latitude"),
                             crs = 4326)
# Make time ranges interpretable
hamburgian_sites$CalBP.Pulses_begin <- as.numeric(
  gsub("^([0-9]*\\.[0-9]*) -.*", 
       "\\1", 
       
       hamburgian_sites$CalBP.Pulses)) * 1000
hamburgian_sites$CalBP.Pulses_end <- as.numeric(
  gsub(".* - ([0-9]*\\.[0-9]*$)",
       "\\1", 
       
       hamburgian_sites$CalBP.Pulses))* 1000
hamburgian_sites$CalBP.Broad_begin <- as.numeric(
  gsub("^([0-9]*\\.[0-9]*) -.*", 
       "\\1", 
       hamburgian_sites$CalBP.Broad))* 1000
hamburgian_sites$CalBP.Broad_end <- as.numeric(
  gsub(".* - ([0-9]*\\.[0-9]*$)", 
       "\\1", 
       hamburgian_sites$CalBP.Broad)) * 1000

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

## 2) Figure 1 - Time-series of LGM climate per site ----
# Calculate key stats

# Mean temp / precip in Hamburgian Culture
mean_temp_ham <- temp_df %>% 
  group_by(Site) %>%
  filter(year_BP >= 13800 & year_BP <= 15000) %>%
  summarise(mean_temp_ham = mean(temp))
mean_precip_ham <- precip_df %>% 
  group_by(Site) %>%
  filter(year_BP >= 13800 & year_BP <= 15000) %>%
  summarise(mean_precip_ham = mean(precip))

# Plots of Temp and Precip
temp_plot <- ggplot(temp_df,
       aes(x = year_BP,
           y = temp,
           colour = Site)) +
  geom_line() +
  geom_segment(data = mean_temp_ham,
             mapping = aes(x = 13500,
                           xend = 15000,
                           y = mean_temp_ham, 
                           yend = mean_temp_ham,
                           colour = Site)) +
  geom_vline(data = hamburgian_sites,
             mapping = aes(xintercept = CalBP.Pulses_begin, 
                           colour = Site))+
  geom_vline(data = hamburgian_sites,
             mapping = aes(xintercept = CalBP.Pulses_end)) +
  labs(x = "Year BP", y = "Temp °C (Bio01)") +
  scale_x_reverse() +
  theme_cowplot(15) +
  theme(legend.position = "none")

save_plot("figures/figure1_temp.png",
          temp_plot,
          base_aspect_ratio = 1.6,
          base_height = 4)

# Plot per site
temp_plot <- temp_plot + facet_wrap(vars(Site))
save_plot("figures/figure1_temp_by_site.png",
          temp_plot,
          base_aspect_ratio = 1.6,
          base_height = 20)

# Precip
precip_plot <- ggplot(precip_df,
                    aes(x = year_BP,
                        y = precip,
                        colour = Site)) +
  geom_line() +
  geom_segment(data = mean_precip_ham,
               mapping = aes(x = 13500,
                             xend = 15000,
                             y = mean_precip_ham, 
                             yend = mean_precip_ham,
                             colour = Site)) +
  geom_vline(data = hamburgian_sites,
             mapping = aes(xintercept = CalBP.Pulses_begin, 
                           colour = Site))+
  geom_vline(data = hamburgian_sites,
             mapping = aes(xintercept = CalBP.Pulses_end)) +
  labs(x = "Year BP", y = "Precipitation mm (Bio12)") +
  scale_x_reverse() +
  theme_cowplot(15) +
  theme(legend.position = "none") 

save_plot("figures/figure1_precip.png",
          precip_plot,
          base_aspect_ratio = 1.6,
          base_height = 4)

# Plot per site
precip_plot <- precip_plot + facet_wrap(vars(Site))
save_plot("figures/figure1_precip.png_by_site.png",
          precip_plot,
          base_aspect_ratio = 1.3,
          base_height = 20)

