
library(tidyverse)
library(cowplot)
library(magick)

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

# Set time range
start_cent <- 150-20
end_cent <- 130-20
time_range <- seq(end_cent, start_cent, 1)

# # Prepare Legend
# predictions_legend <- image_read("figures/helper_figures/predictions_legend.png")
# predictions_legend <- ggdraw() + draw_image(predictions_legend,
#                                             x = 0.44,
#                                             y = 0.6,
#                                             hjust = 0.5,
#                                             vjust = 0.5,
#                                             scale = 0.75)
# 
# # Prepare Scale
# predictions_scale <- image_read("figures/helper_figures/predictions_scale.png")
# predictions_scale <- ggdraw() + draw_image(predictions_scale,
#                                            x = -0.231,
#                                            y = 0.5,
#                                            hjust = 0,
#                                            vjust = 0.5,
#                                            scale = 1)
# 
# # Prepare Predictons Legend Panel
# predictions_legend_panel <- plot_grid(add_title(predictions_scale, 
#                                                 "Predicted Suitabillity (BIOCLIM)",
#                                                 "plain",
#                                                 12,
#                                                 rel_height = 0.2,
#                                                 top_margin = 2,
#                                                 left_margin = 7),
#                                       plot_grid(predictions_legend),
#                                       rel_widths = c(0.66,0.3))

# Animate time series function 
animate_time_series <- function(ts_name){
  # Status
  cat(paste0("Starting ", ts_name, " Model ...\n"))
  
  # set folder
  out_folder <- paste0("figures/preds_ts_", ts_name)
  
  cat(paste0("\nAnimating time-series ... \n"))
  # List files and loade as images
  imgs <- list.files(out_folder, 
                     pattern = ".png", full.names = TRUE)
  
  img_list <- rev(lapply(imgs, image_read))[1:12]
  time_range <- rev(time_range)[1:12]
  # Add legend and Time Stamp
  prediction_legend <- image_read("figures/helper_figures/predictions_legend.png")
  max_index <- max(seq_along(time_range))
  img_list <- lapply(seq_along(time_range), function(index){
    cat(paste0(index,"."))
    k_year_bp_start <- formatC(time_range[index] / 10 + 2.1, 1, format = "f")
    k_year_bp_end <- formatC(time_range[index] / 10 + 2.0, 1, format = "f")
    image1 <- img_list[[index]] %>% image_composite(prediction_legend, offset = "+1490+100") %>%
      image_annotate(paste0(k_year_bp_start, " kyr BP"), size = 70, color = "white", location = "+60+90") %>%
      image_crop("1800x796+0+52")
    # return(image1)
    if(index == max_index) image_morphed <- image_join(rep(image1,10))
    else {
      image2 <- img_list[[index + 1]] %>% image_composite(prediction_legend, offset = "+1490+100") %>%
      image_annotate(paste0(k_year_bp_end, " kyr BP"), size = 70, color = "white", location = "+60+90") %>%
        image_crop("1800x796+0+52")
      image_morphed <- image_morph(c(image1,image2), frames = 10)
      image_morphed <- image_join(c(rep(image1,10), image_morphed))
    }
    return(image_morphed)
  })
  
  # Status
  cat("\nJoining images...")
  # join images 
  img_joined <- image_join(img_list)
  
  # Status
  cat("\nExporting GIF...")
  # Export as gif
  image_write_gif(image = img_joined,
                  path = paste0(out_folder, "/preds_time_series_", ts_name, "_pretty.gif"),
                  delay = 0.025)
  cat("\nDone.")
  return(NULL)
}

animate_time_series("pulse_1")
