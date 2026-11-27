library(tidyverse)

cbbPalette <- c("#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
tol.cbb <- c("#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255", "#332288", "#D55E00", "#000000")


Optimized_test <- read.csv("KEC13,12_specificity_finalized_2024-06-03-Quantification Amplification Results.csv")

head (Optimized_test)

# Define the sample names in the correct order
sample_names <- c("DH44", "DH35", "DH45", "JJ591", "JJ1831", 
                  "JJ64", "JJ1883", "JJ190", "DH46", "NTC")

# Pivot the data and create new variables
Optimized_test_long <- Optimized_test %>%
  pivot_longer(
    cols = -Cycle, 
    names_to = "PlatePosition", 
    values_to = "Fluorescence"
  ) %>%
  # Create a Primer column based on the first letter of the plate position
  mutate(
    Primer = case_when(
      str_starts(PlatePosition, "A") ~ "KEC12",
      str_starts(PlatePosition, "C") ~ "KEC13",
      TRUE ~ NA_character_
    ),
    # Extract the numeric part from the plate position to determine sample order
    Sample_order = as.numeric(str_extract(PlatePosition, "\\d+")),
    #Assign Sample names based on sample order
    Sample = sample_names[Sample_order]
  ) %>%
  # Keep only wells 1 through 10 for each primer set
  filter(Sample_order <= 10) %>%
  # Convert Sample to a factor with levels = sample_names for consistent color mapping
  mutate(Sample = factor(Sample, levels = sample_names))


Optimized_test_plot <- ggplot(Optimized_test_long, aes(x = Cycle, y = Fluorescence, color = Sample)) +
  # Use linewidth for lines in newer ggplot2
  geom_line(linewidth = 0.7, key_glyph = "path") +
  
  scale_color_manual(values = tol.cbb) +
  facet_wrap(~ Primer, scales = "fixed", ncol = 2, nrow = 5) +
  labs(x = "Cycle", y = "RFU") +
  # Dashed line at RFU=100 (optional)
  geom_hline(yintercept = 100, linetype = "dashed", color = "black",linewidth = 0.7 ) +
  # Use a black-and-white base theme
  theme_bw(base_size = 14) +
  # Optionally specify major/minor breaks to match your example grid spacing
  scale_x_continuous(breaks = seq(0, 40, 10), minor_breaks = seq(0, 40, 5)) +
  scale_y_continuous(breaks = seq(0, 1200, 200), minor_breaks = seq(0, 1200, 100)) +
  # Customize theme elements
  theme(
    panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.2),
    panel.grid.minor = element_line(color = "black", linetype = "dotted", size = 0.2),
    
    # Draw a black border around the plotting area
    panel.border = element_rect(color = "black", fill = NA, size = 1),
  
    axis.title       = element_text(face = "bold", size = 10),
    axis.text        = element_text(face = "bold", size = 10),
    legend.title     = element_text(face = "bold", size = 10),
    legend.text      = element_text(face = "bold", size = 10),
    strip.text       = element_text(face = "bold", size = 10),
    legend.key.size  = unit(1.2, "lines")   # overall legend key area
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 2,    # Thicker lines in the legend
        shape = NA,       # No points, show lines
        linetype = 1      # Solid line
      )
    )
  )


# Display the plot
Optimized_test_plot


# Save at high resolution
ggsave("Optimized_test_plot_KEC12,13.png", Optimized_test_plot, width = 10, height = 10, dpi = 400)




#optimized_melting _plot

optimized_melting_peak <- read.csv("KEC13,12_specificity_2024-06-03 - Finalized_Melt Curve Derivative Results.csv")
head (optimized_melting_peak)

# Define the sample names in the correct order
sample_names <- c("DH44", "DH35", "DH45", "JJ591", "JJ1831", 
                  "JJ64", "JJ1883", "JJ190", "DH46", "NTC")

# Pivot the data and create new variables
optimized_melting_peak_long <- optimized_melting_peak %>%
  pivot_longer(
    cols = -Temperature, 
    names_to = "PlatePosition", 
    values_to = "neg_dRFU_dT"
  ) %>%
  # Create a Primer column based on the first letter of the plate position
  mutate(
    Primer = case_when(
      str_starts(PlatePosition, "A") ~ "KEC12",
      str_starts(PlatePosition, "C") ~ "KEC13",
      TRUE ~ NA_character_
    ),
    # Extract the numeric part from the plate position to determine sample order
    Sample_order = as.numeric(str_extract(PlatePosition, "\\d+")),
    #Assign Sample names based on sample order
    Sample = sample_names[Sample_order]
  ) %>%
  # Keep only wells 1 through 10 for each primer set
  filter(Sample_order <= 10) %>%
  #Convert Sample to a factor with levels = sample_names for consistent color mapping
  mutate(Sample = factor(Sample, levels = sample_names))


optimized_melting_peak_plot <- ggplot(optimized_melting_peak_long, aes(x = Temperature, y = neg_dRFU_dT , color = Sample)) +
  geom_line(linewidth = 0.7, key_glyph = "path") +
  
  scale_color_manual(values = tol.cbb) +
  facet_wrap(~ Primer, scales = "fixed", ncol = 2, nrow = 5) +
  labs(
    x = "Temperature (Celsius)",
    y = "-d(RFU)/dT"
  ) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black", size = 0.7) +
  theme_bw(base_size = 14) +  # start with a clean base
  # Make the grid lines dotted and black
  theme(
    
    panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.2),
    panel.grid.minor = element_line(color = "black", linetype = "dotted", size = 0.2),
    
    # Draw a black border around the plotting area
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    
    # Make axis and legend text bold
    axis.title = element_text(face = "bold", size = 10),
    axis.text  = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold", size = 10),
    legend.text  = element_text(face = "bold", size = 10),
    
    # Make facet strip labels bold
    strip.text = element_text(face = "bold", size = 10),
    
    # Increase legend key size
    legend.key.size  = unit(1.2, "lines")   # overall legend key area
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 2,    # Thicker lines in the legend
        shape = NA,       # No points, show lines
        linetype = 1      # Solid line
      )
    )
  )
optimized_melting_peak_plot

# Assuming your plot is stored in the variable `plot`
ggsave("optimized_melting_peak_plot2.png", plot = optimized_melting_peak_plot, width = 10, height = 10, dpi = 400)


library(ggpubr)
combo <- ggarrange(Optimized_test_plot, optimized_melting_peak_plot,
                   ncol = 1, nrow = 2,
                   labels = "AUTO",
                   align = "v",
                   common.legend = T, legend        = "right")
combo
ggsave(
  filename = "Final_optimized_KEC12_13.png",  # path/filename
  plot     = combo,              # the plot object
  width    = 10,                # width in inches
  height   = 10,                # height in inches
  dpi      = 600               # dots per inch (300–600 is standard for publication)
)
