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
library(landscapemetrics)

# Raster options
rasterOptions(progress = "text")

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

# |_ Extract climate variables ----
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

## 2) Climate Variables ----
# |_ Prep Data ----

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

# |_ Histograms ----
variables_to_plot <- c("mean_temp", 
                       "mean_temp_pulse_1",
                       "mean_temp_pulse_2",
                       "mean_precip",
                       "mean_precip_pulse_1",
                       "mean_precip_pulse_2")
x_titles <-  c("All Sites\nMean Temp 14.5k-14.1k BP (°C)", 
               "Pulse 1\nMean Temp 14.5k-14.3k BP (°C)",
               "Pulse 2\nMean Temp 14.2k-14.1k BP (°C)",
               "All Sites\nMean Precip 14.5k-14.1k BP (mm)",
               "Pulse 1\nMean Precip 14.5k-14.3k BP (mm)",
               "Pulse 2\nMean Precip 14.2k-14.1k BP (mm)")
colours <- c(rep(c("black",
                col_pulse1,
                col_pulse2),
              2))
fills <- c(rep(c(col_grey_light,
                   col_pulse1_transparent,
                   col_pulse2_transparent),
                 2))
histogram_plots <- lapply(
  seq_along(variables_to_plot),
  function(index){
    is_temp <- str_detect(variables_to_plot[index], "temp")
    if(is_temp){ 
      x_limits <- c(-3,6)
      x_step <- c(0.5)
    } else{
      x_limits <- c(500,1000)
      x_step <- c(50)
    }
    ggplot(hamburgian_sites) +
      geom_histogram(aes_string(variables_to_plot[index]), 
                     breaks = seq(x_limits[1], x_limits[2], x_step),
                     colour = colours[index],
                     fill = fills[index]) +
      scale_y_continuous(limits = c(0,85),
                         sec.axis=sec_axis(~., breaks = NULL),
                         expand = c(0,0)) +
      scale_x_continuous(sec.axis=sec_axis(~., breaks = NULL)) +
      #(sec.axis=sec_axis(~., breaks = NULL))
      labs(x = x_titles[index], y = "Number of Sites") +
      theme_cowplot(12) 
  })

# Save plots (depricated, see end of script)
# save_plot("figures/climate_histograms.png", 
#           plot_grid(plotlist = histogram_plots,
#                     labels = "AUTO"),
#           base_height = 6)
# save_plot("figures/climate_histograms.eps", 
#           plot_grid(plotlist = histogram_plots,
#                     labels = "AUTO"),
#           base_height = 6)

# |_ Summary Statistics ----

summary_stats <- hamburgian_sites %>%
  st_drop_geometry() %>%
  summarise(n_sites = sum(!is.na(mean_temp)),
            mean = mean(mean_temp, na.rm = T),
            min = min(mean_temp, na.rm = T),
            max = max(mean_temp, na.rm = T),
            sd = sd(mean_temp, na.rm = T)) %>% 
  t()
summary_stats <- hamburgian_sites %>%
  st_drop_geometry() %>%
  summarise(n_sites = sum(!is.na(mean_temp_pulse_1)),
            mean = mean(mean_temp_pulse_1, na.rm = T),
            min = min(mean_temp_pulse_1, na.rm = T),
            max = max(mean_temp_pulse_1, na.rm = T),
            sd = sd(mean_temp_pulse_1, na.rm = T)) %>% 
  t() %>%
  cbind(summary_stats, .)
summary_stats <- hamburgian_sites %>%
  st_drop_geometry() %>%
  summarise(n_sites = sum(!is.na(mean_temp_pulse_2)),
            mean = mean(mean_temp_pulse_2, na.rm = T),
            min = min(mean_temp_pulse_2, na.rm = T),
            max = max(mean_temp_pulse_2, na.rm = T),
            sd = sd(mean_temp_pulse_2, na.rm = T)) %>% 
  t() %>%
  cbind(summary_stats, .)
summary_stats <- hamburgian_sites %>%
  st_drop_geometry() %>%
  summarise(n_sites = sum(!is.na(mean_precip)),
            mean = mean(mean_precip, na.rm = T),
            min = min(mean_precip, na.rm = T),
            max = max(mean_precip, na.rm = T),
            sd = sd(mean_precip, na.rm = T)) %>% 
  t() %>%
  cbind(summary_stats, .)
summary_stats <- hamburgian_sites %>%
  st_drop_geometry() %>%
  summarise(n_sites = sum(!is.na(mean_precip_pulse_1)),
            mean = mean(mean_precip_pulse_1, na.rm = T),
            min = min(mean_precip_pulse_1, na.rm = T),
            max = max(mean_precip_pulse_1, na.rm = T),
            sd = sd(mean_precip_pulse_1, na.rm = T)) %>% 
  t() %>%
  cbind(summary_stats, .)
summary_stats <- hamburgian_sites %>%
  st_drop_geometry() %>%
  summarise(n_sites = sum(!is.na(mean_precip_pulse_2)),
            mean = mean(mean_precip_pulse_2, na.rm = T),
            min = min(mean_precip_pulse_2, na.rm = T),
            max = max(mean_precip_pulse_2, na.rm = T),
            sd = sd(mean_precip_pulse_2, na.rm = T)) %>% 
  t() %>%
  cbind(summary_stats, .) %>% 
  as.data.frame() %>%
  setNames(., c("temp_all", "temp_pulse_1", "temp_pulse_2", 
                "precip_all", "precip_pulse_1", "precip_pulse_2"))
summary_stats[1,] <- apply(summary_stats[1,],
                           2,
                           function(x) formatC(x, digits = 0, format = "f"))
summary_stats[2:4,] <- apply(summary_stats[2:4,],
                             2,
                             function(x) formatC(as.numeric(x), digits = 1, format = "f"))
write_csv(summary_stats, "tables/summary_stats.csv")

# |_ Difference tests ----
# Prepare temperature data
temp_summary <- hamburgian_sites %>%
  st_drop_geometry() %>%
  pivot_longer(c("mean_temp", "mean_temp_pulse_1", "mean_temp_pulse_2"),
               names_to = "variable",
               values_to = "temp") %>%
  select(variable, temp) %>% 
  na.omit()

# Kruskal wallis test for temperature
kruskal <- kruskal.test(temp ~ variable, temp_summary)
kruskal
# Kruskal-Wallis rank sum test
# 
# data:  temp by variable
# Kruskal-Wallis chi-squared = 106.37, df = 2, p-value < 2.2e-16
# => The means are significantly different between the periods

# Run post-hoc wilcox test with adjustment for multiple testing
pw_wilcox <-  pairwise.wilcox.test(temp_summary$temp, temp_summary$variable,
                                   p.adjust.method = "bonferroni")
pw_wilcox
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  temp_summary$temp and temp_summary$variable 
# 
# mean_temp mean_temp_pulse_1
# mean_temp_pulse_1 4.1e-15   -                
#   mean_temp_pulse_2 2.0e-12   2.3e-08          
# 
# P value adjustment method: bonferroni 

# Prepare precipitation data
precip_summary <- hamburgian_sites %>%
  st_drop_geometry() %>%
  pivot_longer(c("mean_precip", "mean_precip_pulse_1", "mean_precip_pulse_2"),
               names_to = "variable",
               values_to = "precip") %>%
  select(variable, precip) %>% 
  na.omit()

# Kruskal wallis test for temperature
kruskal <- kruskal.test(precip ~ variable, precip_summary)
kruskal
# Kruskal-Wallis rank sum test
# 
# data:  precip by variable
# Kruskal-Wallis chi-squared = 8.9768, df = 2, p-value = 0.01124
# => The means are not significantly different between the periods

# Run post-hoc wilcox test with adjustment for multiple testing
pw_wilcox <-  pairwise.wilcox.test(precip_summary$precip, precip_summary$variable,
                                   p.adjust.method = "bonferroni")
pw_wilcox
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  precip_summary$precip and precip_summary$variable 
# 
# mean_precip mean_precip_pulse_1
# mean_precip_pulse_1 0.0095      -                  
#   mean_precip_pulse_2 1.0000      0.1193             
# 
# P value adjustment method: bonferroni 

# |_ Temp Time-Series Plot ----
temp_plot <- ggplot(temp_df %>% filter(year_BP >= 13000 & year_BP < 15500) %>%
                      filter(site %in% unique(hamburgian_sites$site)) %>%
                      full_join(select(mutate(hamburgian_sites, 
                                              site = as.character(site),
                                              chron_association = ordered(chron_association,
                                                                         levels = c("uncertain", "pulse_1", "pulse_2"))), site, chron_association)),
       aes(x = year_BP + 50,
           y = temp,
           group = site)) +
  annotate("rect", xmin = 14100,
           xmax = 14500,
           ymin = -10,
           ymax = 10,
           fill = col_grey_light,
           colour = "black",
           alpha = 1) +
  annotate("rect", xmin = 14300,
           xmax = 14500,
           ymin = -9.5,
           ymax = 9.5,
           fill = col_pulse1_transparent,
           colour = col_pulse1,
           alpha = 1) +
  annotate("rect", xmin = 14100,
           xmax = 14200,
           ymin = -9.5,
           ymax = 9.5,
           fill = col_pulse2_transparent,
           colour = col_pulse2,
           alpha = 1) +
  geom_line(colour = "red3") +
  labs(x = "kyr BP", y = "")  +
  scale_x_reverse(limits = c(15500, 13000),
                  breaks = seq(15500, 13000, -100),
                  labels = rev(c("13.0", "","","","", 
                                 "13.5",  "","","","", 
                                 "14.0",  "","","","", 
                                 "14.5",  "","","","", 
                                 "15.0",  "","","","", 
                                 "15.5")),
                  sec.axis=sec_axis(~., breaks = NULL)) +
  scale_y_continuous(limits = c(-10, 10),
                     breaks = c(-10,-5,0,5,10),
                     #labels = c("-10 °C", "-5 °C", "0 °C", "5 °C", "10 °C"),
                     sec.axis=sec_axis(~., breaks = NULL)) +
  theme_cowplot(10) +
  theme(legend.position = "none",
        plot.margin = margin(5,16,5,4),
        strip.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  )

# Save plot (depricated, see end of script)
# save_plot("figures/figure1_temp.png",
#           temp_plot,
#           base_aspect_ratio = 1.6,
#           base_height = 4)

## Plot per site (depricated)
# temp_plot_by_site <- ggplot(temp_df %>% filter(year_BP >= 13000 & year_BP < 15500) %>%
#                                  filter(site %in% unique(hamburgian_sites$site)) %>%
#                                  mutate(site = ordered(site, levels = levels(fct_reorder(hamburgian_sites$site, hamburgian_sites$lat)))),
#                                aes(x = year_BP + 50,
#                                    y = temp,
#                                    colour = site)) +
#   geom_line() +
#   geom_segment(data = hamburgian_sites,
#                mapping = aes(x = pulse_1_end,
#                              xend = pulse_1_start,
#                              y = mean_temp_pulse_1, 
#                              yend = mean_temp_pulse_1),
#                colour = "black") +
#   geom_segment(data = hamburgian_sites,
#                mapping = aes(x = pulse_2_end,
#                              xend = pulse_2_start,
#                              y = mean_temp_pulse_2, 
#                              yend = mean_temp_pulse_2),
#                colour = "black") +
#   labs(x = "Year BP", y = "Annual Mean Temp (°C)") +
#   scale_x_reverse(limits = c(15500, 13000),
#                   breaks = seq(15500, 13000, -500),
#                   labels = rev(c("13.0k", "13.5k", "14.0k", "14.5k", "15.0k", "15.5k"))) +
#   scale_y_continuous(limits = c(-7, 10),
#                      breaks = c(-5,0,5,10)) +
#   scale_colour_discrete_qualitative(palette = "Dark2") + 
#   geom_text(data = hamburgian_sites,
#             mapping = aes(x = pulse_1_start + 100,
#                           y = mean_temp_pulse_1,
#                           label = paste0("mean:\n", round(mean_temp_pulse_1,1), "°C")),
#             vjust = 0.5,
#             hjust = 1,
#             size = 2,
#             colour = "black") +
#   geom_text(data = hamburgian_sites,
#             mapping = aes(x = pulse_2_start - 200,
#                           y = mean_temp_pulse_2,
#                           label = paste0("mean:\n", round(mean_temp_pulse_2,1), "°C")),
#             colour = "black",
#             vjust = 0.5,
#             hjust = 0,
#             size = 2) +
#   facet_wrap(vars(site),
#              scales = "free") +  
#   theme_cowplot(8) +
#   theme(legend.position = "none",
#         strip.background = element_rect(fill = NA),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#   theme(axis.title = element_text(face = "bold"))
# save_plot("figures/figure1_temp_by_site.png",
#           temp_plot_by_site,
#           base_asp = 1.6,
#           base_height = 20)

# |_ Precip Time-Series Plot ----
precip_plot <- ggplot(precip_df %>% filter(year_BP >= 13000 & year_BP < 15500) %>%
                      mutate(site = ordered(site, levels = levels(fct_reorder(hamburgian_sites$site, hamburgian_sites$lat)))),
                    aes(x = year_BP + 50,
                        y = precip,
                        group = site)) +
  annotate("rect", xmin = 14100,
           xmax = 14500,
           ymin = 200,
           ymax = 1200,
           fill = col_grey_light,
           colour = "black",
           alpha = 1) +
  annotate("rect", xmin = 14300,
           xmax = 14500,
           ymin = 225,
           ymax = 1175,
           fill = col_pulse1_transparent,
           colour = col_pulse1,
           alpha = 1) +
  annotate("rect", xmin = 14100,
           xmax = 14200,
           ymin = 225,
           ymax = 1175,
           fill = col_pulse2_transparent,
           colour = col_pulse2,
           alpha = 1) +
  geom_line(colour = "blue") +
  labs(x = "kyr BP", y = "") +
  scale_x_reverse(limits = c(15500, 13000),
                  breaks = seq(15500, 13000, -100),
                  labels = rev(c("13.0", "","","","", 
                                 "13.5",  "","","","", 
                                 "14.0",  "","","","", 
                                 "14.5",  "","","","", 
                                 "15.0",  "","","","", 
                                 "15.5")),
                  sec.axis=sec_axis(~., breaks = NULL)) +
  scale_y_continuous(limits = c(200, 1200),
                     breaks = seq(200,1200,200),
                     #labels = paste0(seq(200,1200,200), " mm"),
                     sec.axis=sec_axis(~., breaks = NULL)) +
  scale_colour_discrete_qualitative(palette = "Dark2") +
  theme_cowplot(10) +
  theme(legend.position = "none",
        plot.margin = margin(5,16,5,4),
        strip.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Save plot (depricated, see end of script)
# save_plot("figures/figure1_precip.png",
#           precip_plot,
#           base_aspect_ratio = 1.6,
#           base_height = 4)

# Temperature and precipitaton combined (depriacted)
# save_plot("figures/temp_precip_time_series.png",
#           plot_grid(temp_plot, precip_plot, 
#                     labels = "AUTO", 
#                     align = "v", 
#                     axis = "tb"),
#           base_aspect_ratio = 3.2)
# save_plot("figures/temp_precip_time_series.eps",
#           plot_grid(temp_plot, precip_plot, 
#                     labels = "AUTO", 
#                     align = "v", 
#                     axis = "tb"),
#           base_aspect_ratio = 3.2)

# Plot per site (depricated)
# precip_plot_by_site <- ggplot(precip_df %>% filter(year_BP >= 13000 & year_BP < 15500) %>%
#                                 mutate(site = ordered(site, levels = levels(fct_reorder(hamburgian_sites$site, hamburgian_sites$lat)))),
#                               aes(x = year_BP + 50,
#                                   y = precip,
#                                   colour = site)) +
#   geom_line() +
#   geom_segment(data = hamburgian_sites,
#                mapping = aes(x = 14100,
#                              xend = 14500,
#                              y = mean_precip, 
#                              yend = mean_precip),
#                colour = "grey30") +
#   geom_segment(data = hamburgian_sites,
#                mapping = aes(x = pulse_1_end,
#                              xend = pulse_1_start,
#                              y = mean_precip_pulse_1, 
#                              yend = mean_precip_pulse_1),
#                colour = "dodgerblue2") +
#   geom_segment(data = hamburgian_sites,
#                mapping = aes(x = pulse_2_end,
#                              xend = pulse_2_start,
#                              y = mean_precip_pulse_2, 
#                              yend = mean_precip_pulse_2),
#                colour = "black") +
#   labs(x = "Year BP", y = "Annual Precipitation (mm)") +
#   scale_x_reverse(limits = c(15500, 13000),
#                   breaks = seq(15500, 13000, -500),
#                   labels = rev(c("13.0k", "13.5k", "14.0k", "14.5k", "15.0k", "15.5k"))) +
#   scale_y_continuous(limits = c(400, 1100),
#                      breaks = seq(400,1200,200)) +
#   scale_colour_discrete_qualitative(palette = "Dark2") + 
#   geom_text(data = hamburgian_sites,
#             mapping = aes(x = pulse_1_start + 100,
#                           y = mean_precip_pulse_1,
#                           label = paste0("mean:\n", round(mean_precip_pulse_1,1), "mm")),
#             vjust = 0.5,
#             hjust = 1,
#             size = 2,
#             colour = "black") +
#   geom_text(data = hamburgian_sites,
#             mapping = aes(x = pulse_2_start - 200,
#                           y = mean_precip_pulse_2,
#                           label = paste0("mean:\n", round(mean_precip_pulse_2,1), "mm")),
#             colour = "black",
#             vjust = 0.5,
#             hjust = 0,
#             size = 2,) +
#   facet_wrap(vars(site),
#                                     scales = "free") +
#   theme_cowplot(8) +
#   theme(legend.position = "none",
#         strip.background = element_rect(fill = NA),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#   theme(axis.title = element_text(face = "bold"))
# save_plot("figures/figure1_precip_by_site.png",
#           precip_plot_by_site,
#           base_asp = 1.6,
#           base_height = 20)

## 3) Raster Maps ----
area_of_interest <- extent(c(-11,32,49,60)) 
time_range_bp <- data.frame(min_cent = 144-20,
                            max_cent = 141-20,
                            pulse_1_min_cent = 144-20,
                            pulse_1_max_cent = 143-20,
                            pulse_2_min_cent = 141-20,
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
  cat(paste0("Generating ", file_name, "...\n"))
  
  # Create plot object
  map_plot <- levelplot(get(base_raster, envir = .GlobalEnv), 
                        margin = F,
                        main = main_title,
                        at = seq(-12, 12, 1),
                        colorkey = F,
                        xlab = NULL,
                        ylab = NULL,
                        scales =list(draw = F),
                        par.settings = rasterTheme(
                          layout.widths = list(left.padding=-1,
                                               right.padding=-1),
                          layout.heights = list(top.padding=-1, 
                                                bottom.padding=-1))) +
    latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 1)) +
    latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
    latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
    latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "uncertain",],
                                  col = "white", pch = shape_uncertain, cex = 0.8)) +
    latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_1",],
                                  col = col_pulse1, fill = col_pulse1_transparent, pch = shape_pulse_1)) +
    latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_2",],
                                  col = col_pulse2, fill = col_pulse2_transparent, pch = shape_pulse_2)) # +
    # 
    # latticeExtra::layer({
    #   centre_x <- 27.5
    #   centre_y <- 58.3
    #   grid.rect(x=centre_x, y=centre_y,
    #             width=7.5, height=2.5,
    #             gp=gpar(fill="white",
    #                     col = "black",
    #                     alpha = 0.9),
    #             default.units='native')
    #   
    #   
    #   grid.points(#label = "Pulse 1",
    #     x=centre_x-2.75, y=centre_y+2.5/4,
    #     pch = shape_pulse_1,
    #     gp=gpar(cex = 0.45,
    #             col = col_pulse1,
    #             fill = col_pulse1_transparent),
    #     default.units='native')
    #   grid.points(#label = "Pulse 2",
    #     x=centre_x-2.75, y=centre_y,
    #     pch = shape_pulse_2,
    #     gp=gpar(cex = 0.45,
    #             col = col_pulse2,
    #             fill = col_pulse2_transparent),
    #     default.units='native')
    #   grid.points(#label = "Uncertain",
    #     x=centre_x-2.75, y=centre_y-2.5/4,
    #     pch = shape_uncertain,
    #     gp=gpar(cex = 0.65,
    #             col = col_grey_dark),
    #     default.units='native')
    #   
    #   grid.text(label = "Pulse 1",
    #             x=centre_x-1.75, y=centre_y+2.5/4,
    #             #just = "left",
    #             hjust = 0,
    #             vjust = 0.4,
    #             gp=gpar(cex=0.5,
    #                     col = "black"),
    #             default.units='native')
    #   grid.text(label = "Pusle 2",
    #             x=centre_x-1.75, y=centre_y,
    #             #just = "left",
    #             hjust = 0,
    #             vjust = 0.4,
    #             gp=gpar(cex=0.5,
    #                     col = "black"),
    #             default.units='native')
    #   grid.text(label = "Uncertain",
    #             x=centre_x-1.75, y=centre_y-2.5/4,
    #             #just = "left",
    #             hjust = 0,
    #             vjust = 0.4,
    #             gp=gpar(cex=0.5,
    #                     col = "black"),
    #             default.units='native')
    # })
  # Export to file (depricated)
  # png(file_name, 
  #     width = 6,
  #     height = 3,
  #     units = "in",
  #     res = 300)
  # print(map_plot)
  # dev.off()
  return(map_plot)
}

# |__ Maps ----

# Mean temperature map 14.5k-14.1k BP
temp_map_global <- plot_ham_temp_map(base_raster = "hamburgian_mean_temp",
                  main_title = NULL, #"Mean Temperature 14.5k-14.1k BP (°C)",
                  file_name = "figures/mean_temp.png")
# Mean temperature map Pulse 1 14.5k-14.3k BP
temp_map_pulse1 <- plot_ham_temp_map(base_raster = "hamburgian_mean_temp_pulse_1",
                  main_title = NULL, #"Mean Temperature 14.5k-14.3k BP (°C)",
                  file_name = "figures/mean_temp_pulse_1.png")

# Mean temperature map Pulse 2 14.2k-14.1k BP
temp_map_pulse2 <- plot_ham_temp_map(base_raster = "hamburgian_mean_temp_pulse_2",
                  main_title = NULL, #"Mean Temperature 14.2k-14.1k BP (°C)",
                  file_name = "figures/mean_temp_pulse_2.png")

# Plot legend components for alter use
# |__ Colourkey ----
temp_scale <- levelplot(hamburgian_mean_temp, 
                        margin = F,
                        main = NULL,
                        at = seq(-12, 12, 1),
                        colorkey = list(space="bottom", labels = list(cex = 0.5)),
                        xlab = NULL,
                        ylab = NULL,
                        scales =list(draw = F),
                        par.settings = rasterTheme(
                          layout.widths = list(left.padding=0.5,
                                               right.padding=0.5),
                          layout.heights = list(top.padding=-1, 
                                                bottom.padding=-1))) 

# Export to file
png("figures/helper_figures/temp_scale.png", 
    width = 1.5,
    height = 3,
    units = "in",
    res = 300,
    type = "cairo",
    antialias = "none")
print(temp_scale)
dev.off()

# Cut out scale
temp_scale <- image_read("figures/helper_figures/temp_scale.png")
temp_scale <- image_crop(temp_scale, "450x180+0+480")
image_write(temp_scale, "figures/helper_figures/temp_scale.png") 
temp_scale <- ggdraw() + draw_image(temp_scale,
                                    x = -0.105,
                                    y = 0.5,
                                    hjust = 0,
                                    vjust = 0.5,
                                    scale = 1)

# |__ Legend ----
temp_legend <- levelplot(hamburgian_mean_temp, 
                         margin = F,
                         main = NULL,
                         at = seq(-12, 12, 1),
                         colorkey = F,
                         xlab = NULL,
                         ylab = NULL,
                         scales =list(draw = F),
                         par.settings = rasterTheme(
                           layout.widths = list(left.padding=-1,
                                                right.padding=-1),
                           layout.heights = list(top.padding=-1, 
                                                 bottom.padding=-1))) +
  latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "uncertain",],
                                col = "white", pch = shape_uncertain, cex = 0.8)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_1",],
                                col = col_pulse1, fill = col_pulse1_transparent, pch = shape_pulse_1)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_2",],
                                col = col_pulse2, fill = col_pulse2_transparent, pch = shape_pulse_2)) +
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
png("figures/helper_figures/temp_legend.png", 
    width = 5,
    height = 3/6*5,
    units = "in",
    res = 300,
    type = "cairo",
    antialias = "none")
print(temp_legend)
dev.off()

# Cut out legend
temp_legend <- image_read("figures/helper_figures/temp_legend.png")
temp_legend <- image_crop(temp_legend, "263x152+1211+70")
image_write(temp_legend, "figures/helper_figures/temp_legend.png") 
temp_legend <- ggdraw() + draw_image(temp_legend,
                                     x = 0.44,
                                     y = 0.6,
                                     hjust = 0.5,
                                     vjust = 0.5,
                                     scale = 0.75)

# |_ Precip Maps ----

# Define Function to plot precip map
plot_ham_precip_map <- function(base_raster,
                              main_title,
                              file_name){
  cat(paste0("Plotting ", file_name, "...\n"))

  # Create plot object
  map_plot <- levelplot(get(base_raster, envir = .GlobalEnv), 
                        margin = F,
                        main = main_title,
                        at = seq(0, 2500, 2500/20),
                        colorkey = F,
                        xlab = NULL,
                        ylab = NULL,
                        scales = list(draw = F),
                        par.settings = rasterTheme(
                          region = sequential_hcl(
                            100, 
                            "Blues",
                            rev = T),
                          
                          layout.widths = list(left.padding=-1,
                                               right.padding=-1),
                          layout.heights = list(top.padding=-1, 
                                                bottom.padding=-1))) + 
    latticeExtra::layer(sp.polygons(land_for_maps, col = "white", alpha = 1)) +
    latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
    latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
    latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "uncertain",],
                                  col = col_grey_dark, pch = shape_uncertain, cex = 0.8)) +
    latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_1",],
                                  col = col_pulse1, fill = col_pulse1_transparent, pch = shape_pulse_1)) +
    latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_2",],
                                  col = col_pulse2, fill = col_pulse2_transparent, pch = shape_pulse_2)) # +
    # latticeExtra::layer({
    #   centre_x <- 27.5
    #   centre_y <- 58.3
    #   grid.rect(x=centre_x, y=centre_y,
    #             width=7.5, height=2.5,
    #             gp=gpar(fill="white", 
    #                     col = "black",
    #                     alpha = 0.9),
    #             default.units='native')
    #   
    #   grid.points(#label = "Hamburgian",
    #     x=centre_x-2.75, y=centre_y+2.5/4,
    #     pch = shape_pulse_1,
    #     gp=gpar(cex = 0.45, 
    #             col = col_pulse1,
    #             fill = col_pulse1_transparent),
    #     default.units='native')
    #   grid.points(#label = "Havelte Group",
    #     x=centre_x-2.75, y=centre_y,
    #     pch = shape_pulse_2,
    #     gp=gpar(cex = 0.44, 
    #             col = col_pulse2,
    #             fill = col_pulse2_transparent),
    #     default.units='native')
    #   grid.points(#label = "Possibly Ham.",
    #     x=centre_x-2.75, y=centre_y-2.5/4,
    #     pch = 63,
    #     gp=gpar(cex = 0.65, 
    #             col = col_grey_dark),
    #     default.units='native')
    #   
    #   grid.text(label = "Pulse 1",
    #             x=centre_x-1.75, y=centre_y+2.5/4,
    #             #just = "left",
    #             hjust = 0,
    #             vjust = 0.4,
    #             gp=gpar(cex=0.5, 
    #                     col = "black"),
    #             default.units='native')
    #   grid.text(label = "Pulse 2",
    #             x=centre_x-1.75, y=centre_y,
    #             #just = "left",
    #             hjust = 0,
    #             vjust = 0.4,
    #             gp=gpar(cex=0.5, 
    #                     col = "black"),
    #             default.units='native')
    #   grid.text(label = "Uncertain",
    #             x=centre_x-1.75, y=centre_y-2.5/4,
    #             #just = "left",
    #             hjust = 0,
    #             vjust = 0.4,
    #             gp=gpar(cex=0.5, 
    #                     col = "black"),
    #             default.units='native')
    # })
  # Save file (depricated)
  # png(file_name, 
  #     width = 6,
  #     height = 3,
  #     units = "in",
  #     res = 300)
  # print(map_plot)
  # dev.off()
  return(map_plot)
}

# |__ Maps ----
# Mean preciptation map 14.5k-14.1k BP
precip_map_global <- plot_ham_precip_map(base_raster = "hamburgian_mean_precip",
                  main_title = NULL, # "Mean Precipitation 14.5k-14.1k BP (mm)",
                  file_name = "figures/mean_precip.png")

# Mean preciptation map 14.5k-14.3k BP
precip_map_pulse1 <- plot_ham_precip_map(base_raster = "hamburgian_mean_precip_pulse_1",
                    main_title = NULL, #  "Mean Precipitation 14.5k-14.3k BP (mm)",
                    file_name = "figures/mean_precip_pulse_1.png")

# Mean preciptation map 14.2k-14.1k BP
precip_map_pulse2 <- plot_ham_precip_map(base_raster = "hamburgian_mean_precip_pulse_2",
                    main_title = NULL, #  "Mean Precipitation 14.2k-14.1k BP (mm)",
                    file_name = "figures/mean_precip_pulse_2.png")

# Plot legend components for alter use
# |__ Colourkey ----
precip_scale <-  levelplot(hamburgian_mean_precip, 
                           margin = F,
                           main = NULL,
                           at = seq(0, 2500, 2500/20),
                           colorkey = list(space="bottom", labels = list(cex = 0.5)),
                           xlab = NULL,
                           ylab = NULL,
                           scales = list(draw = F),
                           par.settings = rasterTheme(
                             region = sequential_hcl(
                               100, 
                               "Blues",
                               rev = T),
                             
                             layout.widths = list(left.padding=0.5,
                                                  right.padding=0.5),
                             layout.heights = list(top.padding=-1, 
                                                   bottom.padding=-1)))

# Export to file
png("figures/helper_figures/precip_scale.png", 
    width = 1.5,
    height = 3,
    units = "in",
    res = 300,
    type = "cairo",
    antialias = "none",
    )
print(precip_scale)
dev.off()

# Cut out scale
precip_scale <- image_read("figures/helper_figures/precip_scale.png")
precip_scale <- image_crop(precip_scale, "450x180+0+480")
image_write(precip_scale, "figures/helper_figures/precip_scale.png") 
precip_scale <- ggdraw() + draw_image(precip_scale,
                                      x = -0.105,
                                      y = 0.5,
                                      hjust = 0,
                                      vjust = 0.5,
                                      scale = 1)


# |__ Legend ----
precip_legend <- levelplot(hamburgian_mean_precip,
                           margin = F,
                           main = NULL,
                           at = seq(-12, 12, 1),
                           colorkey = F,
                           xlab = NULL,
                           ylab = NULL,
                           scales =list(draw = F),
                           par.settings = rasterTheme(
                             layout.widths = list(left.padding=-1,
                                                  right.padding=-1),
                             layout.heights = list(top.padding=-1, 
                                                   bottom.padding=-1))) +
  latticeExtra::layer(sp.polygons(land_for_maps, col = "black", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ocean_for_maps, col = "NA", fill = "darkblue", alpha = 1)) +
  latticeExtra::layer(sp.polygons(ice_for_maps, col = "darkgrey", fill = "white", alpha = 1)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "uncertain",],
                                col = "white", pch = shape_uncertain, cex = 0.8)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_1",],
                                col = col_pulse1, fill = col_pulse1_transparent, pch = shape_pulse_1)) +
  latticeExtra::layer(sp.points(hamburgian_sites_sp[hamburgian_sites_sp$chron_association == "pulse_2",],
                                col = col_pulse2, fill = col_pulse2_transparent, pch = shape_pulse_2)) +
  
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
png("figures/helper_figures/precip_legend.png", 
    width = 5,
    height = 3/6*5,
    units = "in",
    res = 300,
    type = "cairo",
    antialias = "none")
print(precip_legend)
dev.off()

# Cut out legend
precip_legend <- image_read("figures/helper_figures/precip_legend.png")
precip_legend <- image_crop(precip_legend, "263x152+1211+70")
image_write(precip_legend, "figures/helper_figures/precip_legend.png") 
precip_legend <- ggdraw() + draw_image(precip_legend,
                                       x = 0.44,
                                       y = 0.6,
                                       hjust = 0.5,
                                       vjust = 0.5,
                                       scale = 0.75)

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
             colour = col_grey_dark,
             fill = col_grey_dark,
             alpha = 0.25,
             size = 2,
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
                 fill = chron_association,
                 shape = chron_association,
                 size = chron_association),
             alpha = 1,
             stroke = 1) +
  labs(x = "Mean Temperature 14.5-14.1k BP (°C)",
       y = "\nMean Precipitation 14.5-14.1k BP (mm)") +
  scale_x_continuous(limits = c(-4, 12),
                     breaks = seq(-4, 12, 2),
                     sec.axis=sec_axis(~., breaks = NULL)) +
  scale_y_continuous(limits = c(400, 2600),
                     sec.axis=sec_axis(~., breaks = NULL)) +
  scale_colour_manual(values = c(col_grey_light, col_pulse2, col_pulse1)) +
  scale_size_manual(values = c(2.5,2,2)) +
  scale_fill_manual(values = c(col_grey_light_transparent, col_pulse2_transparent, col_pulse1_transparent)) +
  scale_shape_manual(values = c(shape_uncertain, shape_pulse_2, shape_pulse_1)) +
  annotate("text", x = 9, y = 2600, 
           colour = "black", hjust = 0, vjust = 0.4,
           label = "Pulse 1") +
  annotate("text", x = 9, y = 2475,
           colour = "black", hjust = 0, vjust = 0.4,
           label = "Pulse 2") +
  annotate("text", x = 9, y = 2350,
           colour = "black", hjust = 0, vjust = 0.4,
           label = "Uncertain") +
  annotate("text", x = 9, y = 2225,
           colour = "black", hjust = 0, vjust = 0.4,
           label = "Backg. Sample") +
  annotate("point", x = 8.7, y = 2600, shape = shape_pulse_1, size = 2, stroke = 1,
           colour = col_pulse1, fill = col_pulse1_transparent) +
  annotate("point", x = 8.7, y = 2475, shape = shape_pulse_2, size = 2, stroke = 1,
           colour = col_pulse2, fill = col_pulse2_transparent) +
  annotate("point", x = 8.7, y = 2350, shape = shape_uncertain, size = 3, stroke = 1,
           colour = col_grey_light, fill = col_grey_light) +
  annotate("point", x = 8.7, y = 2225, shape = 21, size = 2,
           colour = col_grey_dark, fill = col_grey_dark, alpha = 0.25) +
  theme_cowplot(14) +
  theme(legend.position = "none")
# save_plot("figures/temp_vs_precip.png",
#           temp_vs_precip_plot,
#           base_aspect_ratio = 1.6,
#           base_height = 4)

# plot the two pusles
temp_vs_precip_plot_by_pulse <- ggplot() +  
  geom_point(data = random_sample,
             aes(x = mean_temp_pulse_1,
                 y = mean_precip_pulse_1),
             colour = "grey60",
             fill = "grey60",
             alpha = 0.25,
             size = 2,
             shape = 21
             ) +
  geom_point(data = random_sample,
             aes(x = mean_temp_pulse_2,
                 y = mean_precip_pulse_2),
             colour = col_grey_dark,
             fill = col_grey_dark,
             alpha = 0.25,
             shape = 21,
             size = 2
             ) +
  geom_point(data = hamburgian_sites %>% filter(chron_association == "uncertain"),
             aes(x = mean_temp,
                 y = mean_precip),
             colour = col_grey_light,
             fill = col_grey_light,
             alpha = 1,
             shape = shape_uncertain,
             size = 2.5,
             stroke = 1
  ) +
  geom_point(data = hamburgian_sites,
             aes(x = mean_temp_pulse_1,
                 y = mean_precip_pulse_1),
             colour = col_pulse1,
             fill = col_pulse1_transparent,
             alpha = 1,
             shape = shape_pulse_1,
             size = 2,
             stroke = 1
             ) +
  geom_point(data = hamburgian_sites,
             aes(x = mean_temp_pulse_2,
                 y = mean_precip_pulse_2),
             colour = col_pulse2,
             fill = col_pulse2_transparent,
             alpha = 1,
             shape = shape_pulse_2,
             size = 2,
             stroke = 1
             ) +
  labs(x = "Mean Temperature of Pulse (°C)",
       y = "\nMean Precipitation of Pulse (mm)") +
  scale_x_continuous(limits = c(-7, 14),
                     breaks = seq(-6, 14, 2),
                     sec.axis=sec_axis(~., breaks = NULL)) +
  scale_y_continuous(limits = c(400, 3100),
                     breaks = seq(0,3000,500),
                     sec.axis=sec_axis(~., breaks = NULL)) +
  annotate("text", x = 5.5, y = 3100, 
           colour = "black", hjust = 0, vjust = 0.4, 
           label = "Pulse 1 (14.5-14.3k BP)") +
  annotate("text", x = 5.5, y = 2950,
           colour = "black", hjust = 0, vjust = 0.4,
           label = "Pulse 2 (14.2-14.1k BP)") +
  annotate("text", x = 5.5, y = 2800,
           colour = "black", hjust = 0, vjust = 0.4,
           label = "Uncertain (14.5-14.1k BP)") +
  annotate("text", x = 5.5, y = 2650,
           colour = "black", hjust = 0, vjust = 0.4,
           label = "Backg. Sample (14.5-14.3k BP)") +
  annotate("text", x = 5.5, y = 2500,
           colour = "black", hjust = 0, vjust = 0.4,
           label = "Backg. Sample (14.2-14.1k BP)") +
  annotate("point", x = 5.1, y = 3100, shape = shape_pulse_1, size = 2, stroke = 1,
           colour = col_pulse1, fill = col_pulse1_transparent) +
  annotate("point", x = 5.1, y = 2950, shape = shape_pulse_2, size = 2, stroke = 1,
           colour = col_pulse2, fill = col_pulse2_transparent) +
  annotate("point", x = 5.1, y = 2800, shape = shape_uncertain, size = 3,
           colour = col_grey_light, fill = col_grey_light) +
  annotate("point", x = 5.1, y = 2650, shape = 21, size = 2,
           colour = "grey60", fill = "grey60", alpha = 0.25) +
  annotate("point", x = 5.1, y = 2500, shape = 21, size = 2,
           colour = col_grey_dark, fill = col_grey_dark, alpha = 0.25) +
  theme_cowplot(14) 
# save_plot("figures/temp_vs_precip_by_pulse.png",
#           temp_vs_precip_plot_by_pulse,
#           base_aspect_ratio = 1.6,
#           base_height = 4)

# 4) BIOCLIM Models ----

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
                             300)
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

# Define function to carry out model fitting and cross validation
fit_bioclim <- function(training_data, absence_data, k){
  training_data$group <- kfold(nrow(training_data), k)
  
  # Prepare output lists
  bioclim_models <- list()
  bioclim_evals <- list()
  
  # Carry out cross validation
  for(i in 1:5){
    # Fit Bioclim envlope model to training data
    bioclim_models[i] <- bioclim(training_data[training_data$group != i, 1:2])
    
    # Evaluate model
    bioclim_evals[i] <- evaluate(p = training_data[training_data$group == i, 1:2],
                                 a = absence_data,
                                 model = bioclim_models[[i]])
  }
  
  # Calc mean auc
  mean_auc <- mean(unlist(lapply(bioclim_evals, 
                                 function(x) x@auc))) 
  
  # Determine optimum thershold 
  # (by maximising specific sensitivit / TSS)
  # and extract mean kappa at this thershold for all models
  mean_kappa <- mean(unlist(lapply(bioclim_evals, 
                                 function(x) {
                                   thresh <- threshold(x)["spec_sens"][1,1]
                                   x@kappa[which(x@t == thresh)]
                                   })))
  
  # Determine optimum thershold 
  # (by maximising specific sensitivit / TSS)
  # and extract mean TSS at this thershold for all models
  mean_tss <- mean(unlist(lapply(bioclim_evals, 
                                 function(x) {
                                   thresh <- threshold(x)["spec_sens"][1,1,]
                                   x@TPR[which(x@t == thresh)] + x@TNR[which(x@t == thresh)] - 1
                                 }))) 
  # fit single model 
  bioclim_model <- bioclim(training_data[,1:2])
  bioclim_eval <- evaluate(training_data[,1:2],
                           absence_data,
                           bioclim_model)
  # Determine threshold with max specific sensitivity / TSS
  bioclim_thresh <-  threshold(bioclim_eval, 'spec_sens') 
  
  # Extract AUC, Kappa, TSS
  auc <- bioclim_eval@auc
  kappa <-  bioclim_eval@kappa[which(bioclim_eval@t == bioclim_thresh)]
  tss <-  bioclim_eval@TPR[which(bioclim_eval@t == bioclim_thresh)] + bioclim_eval@TNR[which(bioclim_eval@t == bioclim_thresh)] - 1
  
  # return all outputs
  list(bioclim_models = bioclim_models,
       bioclim_evals = bioclim_evals,
       mean_auc = mean_auc,
       mean_kappa = mean_kappa,
       mean_tss = mean_tss,
       bioclim_model = bioclim_model,
       bioclim_eval = bioclim_eval,
       bioclim_thresh = bioclim_thresh,
       auc = auc,
       kappa = kappa,
       tss = tss)
}

# Fit model for all three time_periods
set.seed(7)

# Global model
model_global <- fit_bioclim(select(bioclim_training_all, mean_temp, mean_precip),
                            select(absence_data, mean_temp, mean_precip),
                            5)
# Pulse 1
model_pulse_1 <- fit_bioclim(bioclim_training_all %>% 
                               filter(period == "pulse_1") %>%
                               select(mean_temp, mean_precip),
                             absence_data %>% 
                               filter(period == "pulse_1") %>%
                               select(mean_temp, mean_precip),
                             5)

# Pulse 2
model_pulse_2 <- fit_bioclim(bioclim_training_all %>% 
                               filter(period == "pulse_2") %>%
                               select(mean_temp, mean_precip),
                             absence_data %>% 
                               filter(period == "pulse_2") %>%
                               select(mean_temp, mean_precip),
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

write_csv(eval_summary, "tables/model_evaluation.csv")

# Predict for core time periods
predictions_global <- predict(climate_mean,
                              model_global[["bioclim_model"]])
predictions_pulse_1 <- predict(climate_pulse_1,
                                   model_pulse_1[["bioclim_model"]]) 
predictions_pulse_2 <- predict(climate_pulse_2, 
                                   model_pulse_2[["bioclim_model"]])

# Apply threshold
predictions_global_thresh <- predictions_global > model_global[["bioclim_thresh"]]
predictions_pulse_1_thresh <- predictions_pulse_1 > model_pulse_1[["bioclim_thresh"]]
predictions_pulse_2_thresh <- predictions_pulse_2 > model_pulse_2[["bioclim_thresh"]]

predictions_global <- mask(predictions_global, reclassify(predictions_global_thresh, c(-0.1, 0.1, NA)))
predictions_pulse_1 <- mask(predictions_pulse_1, reclassify(predictions_pulse_1_thresh, c(-0.1, 0.1, NA)))
predictions_pulse_2 <- mask(predictions_pulse_2, reclassify(predictions_pulse_2_thresh, c(-0.1, 0.1, NA)))

# Visualise predictions

# Define function to visualise predictions
plot_preds_map <- function(base_raster,
                           main_title,
                           file_name,
                           palette = "default",
                           colour_key = F,
                           key_at = seq(0,1,0.05),
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
    
    # latticeExtra::layer({
    #   centre_x <- 27.5
    #   centre_y <- 58.3
    #   grid.rect(x=centre_x, y=centre_y,
    #             width=7.5, height=2.5,
    #             gp=gpar(fill="white",
    #                     col = "black",
    #                     alpha = 0.9),
    #             default.units='native')
    #   
    #   
    #   grid.points(#label = "Pulse 1",
    #     x=centre_x-2.75, y=centre_y+2.5/4,
    #     pch = shape_pulse_1,
    #     gp=gpar(cex = 0.45,
    #             col = col_pulse1,
    #             fill = col_pulse1_transparent),
    #     default.units='native')
    #   grid.points(#label = "Pulse 2",
    #     x=centre_x-2.75, y=centre_y,
    #     pch = shape_pulse_2,
    #     gp=gpar(cex = 0.45,
    #             col = col_pulse2,
    #             fill = col_pulse2_transparent),
    #     default.units='native')
    #   grid.points(#label = "Uncertain",
    #     x=centre_x-2.75, y=centre_y-2.5/4,
    #     pch = shape_uncertain,
    #     gp=gpar(cex = 0.65,
    #             col = col_grey_dark),
    #     default.units='native')
    #   
    #   grid.text(label = "Pulse 1",
    #             x=centre_x-1.75, y=centre_y+2.5/4,
    #             #just = "left",
    #             hjust = 0,
    #             vjust = 0.4,
    #             gp=gpar(cex=0.5,
    #                     col = "black"),
    #             default.units='native')
    #   grid.text(label = "Pusle 2",
    #             x=centre_x-1.75, y=centre_y,
    #             #just = "left",
    #             hjust = 0,
    #             vjust = 0.4,
    #             gp=gpar(cex=0.5,
    #                     col = "black"),
    #             default.units='native')
    #   grid.text(label = "Uncertain",
    #             x=centre_x-1.75, y=centre_y-2.5/4,
    #             #just = "left",
    #             hjust = 0,
    #             vjust = 0.4,
    #             gp=gpar(cex=0.5,
    #                     col = "black"),
    #             default.units='native')
    # })
  
  # Write out file if needed
  if(to_file == T){
    png(file_name, 
        width = 6,
        height = 3,
        units = "in",
        res = 300
        type = "cairo",
        antialias = "none")
    print(map_plot)
    dev.off()
  }
  return(map_plot)
}

# |_ Mean ----
preds_map_global <- plot_preds_map(base_raster = predictions_global,
                                   main = NULL, #"Global Model Suitability 14.5k-14.1k BP ", 
                                   file_name = "figures/predictions_global.png")
# plot_preds_map(base_raster = predictions_global_thresh,
#                main = "Global Model Pres./Abs. 14.5k-14.1k BP ", 
#                file_name = "figures/predictions_global_thresh.png",
#                palette = c("grey30", "#f38f32"),
#                colour_key = F)

# |_ Pulse 1----
preds_map_pulse_1 <- plot_preds_map(base_raster = predictions_pulse_1,
               main = NULL, #"Pulse 1 Model Suitability 14.5k-14.3k BP ", 
               file_name = "figures/predictions_pulse_1.png")
# plot_preds_map(base_raster = predictions_pulse_1_thresh,
#                main = "Pulse 1 Model Pres./Abs. 14.5k-14.3k BP ", 
#                file_name = "figures/predictions_pulse_1_thresh.png",
#                palette = c("grey30", "#f38f32"),
#                colour_key = F)

# |_ Pulse 2 ----
preds_map_pulse_2 <- plot_preds_map(base_raster = predictions_pulse_2,
               main = NULL, #"Pulse 2 Model Suitability 14.2k-14.1k BP ", 
               file_name = "figures/predictions_pulse_2.png")
# plot_preds_map(base_raster = predictions_pulse_2_thresh,
#                main = "Pulse 2 Model Pres./Abs. 14.2k-14.1k BP ", 
#                file_name = "figures/predictions_pulse_2_thresh.png",
#                palette = c("grey30", "#f38f32"),
#                colour_key = F)

# |_ Colourkey ----
palette <- sequential_hcl(
  103, 
  "Inferno")
palette <- palette[c(-1,-2,-3)]
predictions_scale <-  levelplot(predictions_global, 
                         margin = F,
                         main = NULL,
                         colorkey = list(space="bottom", labels = list(cex = 0.5)),
                         at = seq(0,1,0.05),
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
png("figures/helper_figures/predictions_scale.png", 
    width = 1.5,
    height = 3,
    units = "in",
    res = 300,
    type = "cairo",
    antialias = "none")
print(predictions_scale)
dev.off()

# Cut out scale
predictions_scale <- image_read("figures/helper_figures/predictions_scale.png")
predictions_scale <- image_crop(predictions_scale, "450x180+0+480")
image_write(predictions_scale, "figures/helper_figures/predictions_scale.png") 
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
png("figures/helper_figures/predictions_legend.png", 
    width = 5,
    height = 3/6*5,
    units = "in",
    res = 300,
    type = "cairo",
    antialias = "none")
print(predictions_legend)
dev.off()

# Cut out legend
predictions_legend <- image_read("figures/helper_figures/predictions_legend.png")
predictions_legend <- image_crop(predictions_legend, "263x152+1211+70")
image_write(predictions_legend, "figures/helper_figures/predictions_legend.png") 
predictions_legend <- ggdraw() + draw_image(predictions_legend,
                                     x = 0.44,
                                     y = 0.6,
                                     hjust = 0.5,
                                     vjust = 0.5,
                                     scale = 0.75)

predictions_legend

# 5) Time-Series for Predictions ----


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
predict_time_series <- function(bioclim_model_list, model_name, file_name){
  # Status
  cat(paste0("Starting ", model_name, " Model ...\n"))
  
  # create folder
  out_folder <- paste0("figures/preds_ts_", file_name)
  dir.create(out_folder)
  
  # Calculate predictions
  cat(paste0("Calculating predictions ... \n"))
  bioclim_model <- bioclim_model_list[["bioclim_model"]]
  bioclim_thresh <- bioclim_model_list[["bioclim_thresh"]]
  
  preds_time_series <- lapply(time_series_climate, function(climate_raster){
    cat(paste0(names(climate_raster), " \n"))
    
    # update climate raster names to match model
    climate_raster <- setNames(climate_raster, c("mean_temp", "mean_precip"))
    # Predict for both time periods
    predictions <- predict(climate_raster,
                           bioclim_model)
    # Apply threshold
    predictions_thresh <- predictions > bioclim_thresh
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


# 6) Fragstats ----
# Prepare time series of thresholded maps converted to an arbitary 
# reference system with 
threshold_ts <- function(raw_preds){
  cat(".")
  # Reclassify raster
  reclassed <- reclassify(raw_preds, c(NA, NA, 0, 0,1,1))
  # Mask land
  reclassed <-  mask(reclassed, land_for_maps) 
  # Copy to raster with arbitary reference system in m
  thresholded <- raster(nrows = nrow(reclassed),
                       ncols = ncol(reclassed),
                       xmn = 0,
                       ymn = 0,
                       xmx = 1000 * ncol(reclassed),
                       ymx = 1000 * nrow(reclassed),
                       vals = getValues(reclassed),
                       crs = "+proj=longlat +datum=WGS84 +no_defs +units=m"
                       )
  return(thresholded)
}
ts_global_thresh <- lapply(ts_global, threshold_ts)
ts_pulse_1_thresh <- lapply(ts_pulse_1, threshold_ts)
ts_pulse_2_thresh <- lapply(ts_pulse_2, threshold_ts)

# Specify functions to calculate landscape metrics
calc_lsm <- function(preds_thresholded){
  cat(".")
  # Proportion of cells suitable
  prop_cells <- sum(getValues(preds_thresholded) == 1, na.rm = T) / sum(!is.na(getValues(preds_thresholded)))
  
  # return NA if prop = 0
  if(prop_cells == 0) {
    return(data.frame(prop_cells = NA,
                      n_patches = NA,
                      patch_area_mean = NA,
                      patch_area_sd = NA,
                      mean_nearest_neighbour = NA,
                      patch_cohesion = NA))
  }
  # Number of patches
  n_patches <- lsm_l_np(preds_thresholded) %>% pull(value)
  # Mean patch area
  patch_area_mean <- lsm_c_area_mn(preds_thresholded) %>% 
    filter(class == 1) %>%
    pull(value) * 0.01
  # Path area sd
  patch_area_sd <- lsm_c_area_sd(preds_thresholded) %>% 
    filter(class == 1) %>%
    pull(value) * 0.01
  # Nearest neighbour distance
  mean_nearest_neighbour <- lsm_p_enn(preds_thresholded) %>% 
    filter(class == 1) %>%
    pull(value) %>% mean(.) / 1000 
  # Patch cohesion
  patch_cohesion <- lsm_c_cohesion(preds_thresholded) %>% 
    filter(class == 1) %>% pull(value)
  
  return(data.frame(prop_cells,
                    n_patches,
                    patch_area_mean,
                    patch_area_sd,
                    mean_nearest_neighbour,
                    patch_cohesion))
}

# Calculate metrics
# Main predictions
global_stats <- calc_lsm(threshold_ts(predictions_global))
pulse_1_stats <- calc_lsm(threshold_ts(predictions_pulse_1))
pulse_2_stats <- calc_lsm(threshold_ts(predictions_pulse_2))
stats_all <- bind_rows(global_stats,
                       pulse_1_stats,
                       pulse_2_stats) %>%
  mutate(model = c("Global", "Pulse 1", "Pulse 2"),
         time = c("14.5k-14.1k BP",
                  "14.5k-14.3k BP",
                  "14.2k-14.1k BP")) %>%
  select(model, time, prop_cells:patch_cohesion)
names(stats_all) <- c("Model",
                      "Time-Window",
                      "Prop. cells suitable",
                      "Number of Patches",
                      "Mean Patch Area (km2)",
                      "Std. Dev. Patch Area (km2)",
                      "Mean Dist to neares neighbour (km)",
                      "Patch Cohesion")
write_csv(stats_all, "tables/landscape_stats.csv")

# Time-series
ts_global_patch_stats <- bind_rows(lapply(ts_global_thresh, calc_lsm))
ts_pulse_1_patch_stats <- bind_rows(lapply(ts_pulse_1_thresh, calc_lsm))
ts_pulse_2_patch_stats <- bind_rows(lapply(ts_pulse_2_thresh, calc_lsm))

ts_global_patch_stats$start_bp <- formatC(time_range / 10 + 2.1, 1, format = "f")
ts_global_patch_stats$end_bp <- formatC(time_range / 10 + 2, 1, format = "f")
ts_pulse_1_patch_stats$start_bp <- formatC(time_range / 10 + 2.1, 1, format = "f")
ts_pulse_1_patch_stats$end_bp <- formatC(time_range / 10 + 2, 1, format = "f")
ts_pulse_2_patch_stats$start_bp <- formatC(time_range / 10 + 2.1, 1, format = "f")
ts_pulse_2_patch_stats$end_bp <- formatC(time_range / 10 + 2, 1, format = "f")
# Save
write_csv(ts_global_patch_stats, "tables/landscape_stats_ts_global.csv")
write_csv(ts_pulse_1_patch_stats, "tables/landscape_stats_ts_pulse_1.csv")
write_csv(ts_pulse_2_patch_stats, "tables/landscape_stats_ts_pulse_2.csv")

# 7) Final Figures ----

# |_ Figure 1 - Histograms & Points ----
# Check whether histograms exists
exists("histogram_plots")

# Add point drawings to panels b,c,e,f
point_drawings <- image_read("figures/helper_figures/ham_points.jpeg")
# Crop Pulse 1 point ("classic hamburgian )
point_pulse_1 <- image_crop(point_drawings, "410x1165+1280+650")
# Crop Pulse 2 point ("Havelte Group")
point_pulse_2 <- image_crop(point_drawings, "500x1814+0+10")

# Add to the respective histogram plots
histogram_plots[[2]] <- ggdraw() + 
  draw_plot(histogram_plots[[2]]) + 
  draw_image(point_pulse_1,
             scale = 0.25,
             x = 0.675,
             y = 0.2,
             width = 410 / 1000,
             height = 1164 / 1000)
histogram_plots[[3]] <- ggdraw() + 
  draw_plot(histogram_plots[[3]]) + 
  draw_image(point_pulse_2,
             scale = 0.25,
             x = 0.625,
             y = -0.2,
             width = 500 / 1000,
             height = 1804 / 1000)
histogram_plots[[5]] <- ggdraw() + 
  draw_plot(histogram_plots[[5]]) + 
  draw_image(point_pulse_1,
             scale = 0.25,
             x = 0.675,
             y = 0.2,
             width = 410 / 1000,
             height = 1164 / 1000)
histogram_plots[[6]] <- ggdraw() + 
  draw_plot(histogram_plots[[6]]) + 
  draw_image(point_pulse_2,
             scale = 0.25,
             x = 0.625,
             y = -0.2,
             width = 500 / 1000,
             height = 1804 / 1000)

# Save figure
save_plot("figures/figure_1_climate_hist_and_points.png", 
          plot_grid(plotlist = histogram_plots,
                    labels = "AUTO"),
          base_height = 6)
save_plot("figures/figure_1_climate_hist_and_points.eps", 
          plot_grid(plotlist = histogram_plots,
                    labels = "AUTO"),
          base_height = 6)

# |_ Figure 2 - Climate Context ----
# List of plots by row
plot_list <- c("temp_plot", "precip_plot",
               "temp_map_global", "precip_map_global",
               "temp_map_pulse1", "precip_map_pulse1",
               "temp_map_pulse2", "precip_map_pulse2",
               "temp_scale", "precip_scale",
               "temp_legend", "precip_legend")

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

# Prepare legend panels for the bottom
temp_legend_panel <- plot_grid(add_title(temp_scale, 
                                         "Mean Annual Temperature (°C)",
                                         "plain",
                                         8,
                                         rel_height = 0.2,
                                         top_margin = 0,
                                         left_margin = 15.5),
                               temp_legend,
                               rel_widths = c(0.66,0.33))
precip_legend_panel <- plot_grid(add_title(precip_scale, 
                                         "Mean Annual Precipitation (mm)",
                                         "plain",
                                         8,
                                         rel_height = 0.2,
                                         top_margin = 0,
                                         left_margin = 15.5),
                               precip_legend,
                               rel_widths = c(0.66,0.33))
figure_2 <- plot_grid(
  plotlist = list(add_title(temp_plot, "Mean Annual Temperature (°C)", "bold", rel_height = 0.1*0.6, top_margin = 3),
                  add_title(precip_plot, "Mean Annual Precipitation (mm)", "bold", rel_height = 0.1*0.6, top_margin = 3),
                  add_title(temp_map_global, "Global Period (14.5 - 14.1 kyr BP)", "plain"),
                  add_title(precip_map_global, "Global Period (14.5 - 14.1 kyr BP)","plain"),
                  add_title(temp_map_pulse1, "Pulse 1 (14.5 - 14.3 kyr BP)", "plain"),
                  add_title(precip_map_pulse1, "Pulse 1 (14.5 - 14.3 kyr BP)", "plain"),
                  add_title(temp_map_pulse2, "Pulse 2 (14.2 - 14.1 kyr BP)", "plain"),
                  add_title(precip_map_pulse2, "Pulse 2 (14.2 - 14.1 kyr BP)", "plain"),
                  temp_legend_panel,
                  precip_legend_panel
                  ),
  ncol = 2,
  labels = c("A", "B", "C", "D", rep("",6)),
  rel_heights = c(1.3,
                  rep(1,3),
                  0.5),
  label_size = 8)
save_plot("figures/figure_2-climate_context.png", 
          figure_2, 
          base_height = 2,
          ncol = 2,
          nrow = 4,
          base_asp = 1.8)
save_plot("figures/figure_2-climate_context.eps", 
          figure_2, 
          base_height = 2,
          ncol = 2,
          nrow = 4,
          base_asp = 1.8)

# |_ Figure 3 - Climate Space ----
# temperature vs precipitation plots
# Panel both figures up into one
save_plot("figures/figure_3-climate_space.png",
          plot_grid(temp_vs_precip_plot, 
                    temp_vs_precip_plot_by_pulse, 
                    labels = "AUTO", 
                    align = "v", 
                    axis = "tb"),
          base_aspect_ratio = 3.2,
          base_height = 4)
save_plot("figures/figure_3-climate_space.eps",
          plot_grid(temp_vs_precip_plot, 
                    temp_vs_precip_plot_by_pulse, 
                    labels = "AUTO", 
                    align = "v", 
                    axis = "tb"),
          base_aspect_ratio = 3.2,
          base_height = 4)

# |_ Figure 4 - Suitabillity Predictions ----
# List of plots 
plot_list <- c("preds_map_global",
               "preds_map_pulse_1",
               "preds_map_pulse_2")

# Check all object exist
lapply(plot_list, function(x) exists(x))

# Combine legend and colour scale
predictions_legend_panel <- plot_grid(add_title(predictions_scale, 
                                         "Predicted Suitabillity (BIOCLIM)",
                                         "plain",
                                         8,
                                         rel_height = 0.2,
                                         top_margin = 2,
                                         left_margin = 15.5),
                               predictions_legend,
                               rel_widths = c(0.66,0.33))

# Add titles and combine into one figure
figure_4 <- plot_grid(
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
save_plot("figures/figure_4_suitability_predictions.png", 
          figure_4,
          base_height = 2,
          ncol = 1,
          nrow = 3,
          base_asp = 1.85)
save_plot("figures/figure_4_suitability_predictions.eps", 
          figure_4,
          base_height = 2,
          ncol = 1,
          nrow = 3,
          base_asp = 1.85)

# |_ Figure 5 - Suitabillity Time-Series ----
# read in time-series plots (centuries 14.7-14.1)
ts_maps_global <- lapply(rev(list.files("figures/preds_ts_global/", pattern = ".png", full.names = T)[11:17]),
                         function(x) ggdraw() + draw_image(image_read(x), scale = 0.95))
ts_maps_pulse_1 <- lapply(rev(list.files("figures/preds_ts_pulse_1/", pattern = ".png", full.names = T)[11:17]),
                          function(x) ggdraw() + draw_image(image_read(x), scale = 0.95))
ts_maps_pulse_2 <- lapply(rev(list.files("figures/preds_ts_pulse_2/", pattern = ".png", full.names = T)[11:17]),
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
predictions_scale <- image_read("figures/helper_figures/predictions_scale.png")
predictions_scale <- ggdraw() + draw_image(predictions_scale,
                                           x = -0.231,
                                           y = 0.5,
                                           hjust = 0,
                                           vjust = 0.5,
                                           scale = 1)
predictions_legend_panel <- plot_grid(add_title(predictions_scale, 
                                                "Predicted Suitabillity (BIOCLIM)",
                                                "plain",
                                                12,
                                                rel_height = 0.2,
                                                top_margin = 2,
                                                left_margin = 7),
                                      plot_grid(predictions_legend),
                                      rel_widths = c(0.66,0.3))
# Assemble plot
figure_5_body <- plot_grid(ts_time_scale,
                      ts_maps_global,
                      ts_maps_pulse_1,
                      ts_maps_pulse_2,
                      ncol = 4,
                      rel_widths = c(0.3,1,1,1))
# add labels and legend
figure_5_labels <- lapply(
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
figure_5_labels <- plot_grid(
  NULL, figure_5_labels[[1]],
  figure_5_labels[[2]], 
  figure_5_labels[[3]],
  rel_widths =c(0.3,1,1,1),
  nrow= 1)
figure_5_legend <- plot_grid(
  NULL, predictions_legend_panel, NULL,
  rel_widths = c(0.3,1.8,1.2),
  nrow = 1
)
figure_5 <- plot_grid(figure_5_labels,
                      figure_5_body,
                      figure_5_legend,
                      ncol = 1,
                      rel_heights = c(0.015,1,0.1))
save_plot("figures/figure_5-suitability_time_series.png",
          figure_5,
          base_height = 1.9,
          ncol = 3,
          nrow = 7,
          base_asp = 2)
save_plot("figures/figure_5-suitability_time_series.eps",
          figure_5,
          base_height = 1.9,
          ncol = 3,
          nrow = 7,
          base_asp = 2)

# |_ Animations ----

# Generate Animation
animate_time_series <- function(ts_name){
  # Status
  cat(paste0("Starting ", ts_name, " Model ...\n"))
  
  # set folder
  out_folder <- paste0("figures/preds_ts_", ts_name)
  
  cat(paste0("\nAnimating time-series ... \n"))
  # List files and loade as images
  imgs <- list.files(out_folder, 
                     pattern = ".png", full.names = TRUE)
  
  img_list <- lapply(imgs, image_read)
  
  # Add legend and Time Stamp
  prediction_legend <- image_read("figures/helper_figures/predictions_legend.png")
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
            path = paste0(out_folder, "/preds_time_series_", ts_name, ".gif"),
            delay = 0.5)
  cat("\nDone.")
  return(NULL)
}
animate_time_series("global")
animate_time_series("pulse_1")
animate_time_series("pulse_2")
