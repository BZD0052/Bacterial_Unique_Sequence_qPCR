library(tidyverse)


cbbPalette <- c("#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
tol.cbb <- c("#117733", "#44AA99", "#88CCEE", "#000000")


# 2) Read raw data
qPCR_LOD <- read.csv("Final_KEC12&13_primer_LOD_2024-06-18-QuantificationAmplificationResults.csv", check.names = FALSE)

# 3) Pivot to long
qPCR_LOD_long <- qPCR_LOD %>%
  pivot_longer(
    cols = -Cycle,
    names_to = "Well",
    values_to = "Fluorescence"
  )

# 4) Create plate layout info (example)
qPCR_LOD_plate_layout <- expand.grid(
  Row = LETTERS[1:8],
  Col = 1:11
) %>%
  mutate(
    Well = paste0(Row, Col),
    Primer = case_when(
      Row %in% c("A","B","C") ~ "Primer1",
      Row %in% c("E","F","G") ~ "Primer2",
      TRUE ~ NA_character_
    ),
    Dilution = case_when(
      Row == "A" ~ "0.1pg",
      Row == "B" ~ "0.01pg",
      Row == "C" ~ "0.001pg (1fg)",
      Row == "E" ~ "0.1pg",
      Row == "F" ~ "0.01pg",
      Row == "G" ~ "0.001pg (1fg)",
      TRUE ~ "Unknown"
    ),
    IsNTC = if_else(Col == 11, TRUE, FALSE)
  )

# 5) Merge
qPCR_LOD_long_merged <- qPCR_LOD_long %>%
  left_join(qPCR_LOD_plate_layout, by = "Well") %>%
  mutate(Sample = if_else(IsNTC, "NTC", Dilution))

qPCR_LOD_plot <- ggplot(qPCR_LOD_long_merged, aes(x = Cycle, y = Fluorescence,
                      color = Sample, group = Well)) +
  geom_line(linewidth = 0.4, key_glyph = "path") +
  coord_cartesian(xlim = c(20, 40)) +
  facet_wrap(~ Primer, scales = "fixed",
             labeller = as_labeller(c("Primer1" = "KEC 12", "Primer2" = "KEC 13"))) +
  scale_color_manual(
    values = c(name ="Samples",
               "0.1pg" = "#117733",
               "0.01pg" = "#D55E00",
               "0.001pg (1fg)" = "#88CCEE",
               "NTC"   = "black" # or another color for the NTC
    )  # or another color for the NTC
    ) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.6) +
  theme_bw(base_size = 14) +
   scale_x_continuous(breaks = seq(0, 40, 10), minor_breaks = seq(0, 40, 5), expand = expansion(mult = c(0.1, 0.3))) +
  scale_y_continuous(name ="RFU",breaks = seq(0, 1200, 200), minor_breaks = seq(0, 1200, 100)) +
  theme(
    panel.grid.major = element_line(color = "black", linetype = "dotted", linewidth = 0.2),
    panel.grid.minor = element_line(color = "black", linetype = "dotted", linewidth = 0.2),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title       = element_text(face = "bold", size = 9),
    axis.text        = element_text(face = "bold", size = 10),
    legend.position = "top",
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
      ))
    )
      
qPCR_LOD_plot



# Assuming your plot is stored in the variable `plot`
ggsave("qPCR_LOD_plot2.png", plot = qPCR_LOD_plot, width = 10, height = 6, dpi = 400)


# Set the threshold
threshold <- 100

# Calculate the Cq value for each individual well (each replicate)
qPCR_LOD_Cq_value <- qPCR_LOD_long_merged %>%
  group_by(Primer, Sample, Well) %>%  # each well is its own group
  arrange(Cycle) %>%
  summarize(
    crossing_cycle = if (max(Fluorescence, na.rm = TRUE) < threshold) {
      NA_real_
    } else {
      approx(x = Fluorescence, y = Cycle, xout = threshold)$y
    },
    .groups = "drop"
  )

# Calculate the average Cq for each Primer-Sample combination
qPCR_LOD_Cq_joined <- qPCR_LOD_long_merged %>%
  group_by(Primer, Sample, Well) %>%  # each well is its own group
  arrange(Cycle) %>%
  summarize(
    crossing_cycle = if (max(Fluorescence, na.rm = TRUE) < threshold) {
      NA_real_
    } else {
      approx(x = Fluorescence, y = Cycle, xout = threshold)$y
    },
    .groups = "drop_last"  # Drop the last grouping (Well), keeping Primer and Sample groups
  ) %>%
  mutate(avg_Cq = mean(crossing_cycle, na.rm = TRUE)) %>%  # Average across replicates
  ungroup()

# View the joined result (with average Cq for each Primer-Sample combination)
qPCR_LOD_Cq_joined
qPCR_LOD_Cq_value
#save the computed table (e.g., the cq_values data frame) in CSV format
write.csv(qPCR_LOD_Cq_joined, "qPCR_LOD_Cq_joined.csv", row.names = FALSE)
write.csv(qPCR_LOD_Cq_value, "qPCR_LOD_Cq_value.csv", row.names = FALSE)

#####################standard curve melting peak
library(tidyverse)

# ----- 1. Read the melt-curve derivative data -----
LOD_melting_peak <- read.csv(
  "Final_KEC12&13_primer_LOD_2024-06-18-MeltCurveDerivativeResults.csv",
  check.names = FALSE,
  na.strings = "na"
)

# ----- 2. Pivot the data to long format -----
LOD_melt_long <- LOD_melting_peak %>%
  pivot_longer(
    cols = -Temperature,      # keep Temperature as is
    names_to = "Well",
    values_to = "Derivative"
  )

# ----- 3. Create plate layout info (if not already created) -----
qPCR_LOD_plate_layout <- expand.grid(
  Row = LETTERS[1:8],
  Col = 1:11
) %>%
  mutate(
    Well = paste0(Row, Col),
    Primer = case_when(
      Row %in% c("A", "B", "C") ~ "Primer1",
      Row %in% c("E", "F", "G") ~ "Primer2",
      TRUE ~ NA_character_
    ),
    Dilution = case_when(
      Row == "A" ~ "0.1pg",
      Row == "B" ~ "0.01pg",
      Row == "C" ~ "0.001pg (1fg)",
      Row == "E" ~ "0.1pg",
      Row == "F" ~ "0.01pg",
      Row == "G" ~ "0.001pg (1fg)",
      TRUE ~ "Unknown"
    ),
    IsNTC = if_else(Col == 11, TRUE, FALSE)
  )

# ----- 4. Join melt data with plate layout and create the 'Sample' column -----
LOD_melt_long_merged <- LOD_melt_long %>%
  left_join(qPCR_LOD_plate_layout, by = "Well") %>%
  mutate(Sample = if_else(IsNTC, "NTC", Dilution))

# ----- 5. Plot the melt curves using the same theme and style as your qPCR LOD plot -----
LOD_melt_curve_plot <- ggplot(LOD_melt_long_merged, aes(x = Temperature, y = Derivative, color = Sample, group = Well)) +
  geom_line(linewidth = 0.4, key_glyph = "path") +
  facet_wrap(~ Primer, scales = "fixed",
             labeller = as_labeller(c("Primer1" = "KEC 12", "Primer2" = "KEC 13"))) +
  labs(
    x = "Temperature (Celsius)",
    y = "-d(RFU)/dT"
  ) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black", size = 0.4) +
  theme_bw(base_size = 14) +
  scale_color_manual(
    name = "Samples",
    values = c(
      "0.1pg" = "#117733",
      "0.01pg" = "#D55E00",
      "0.001pg (1fg)" = "#88CCEE",
      "NTC"   = "black"
    )
  ) + # Make the grid lines dotted and black
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


print(LOD_melt_curve_plot)
# Assuming your plot is stored in the variable `plot`
ggsave("LOD_melt_curve_plot.png", plot = LOD_melt_curve_plot, width = 10, height = 6, dpi = 400)


library(ggpubr)

#Set aspect ratio in each individual plot:
  qPCR_LOD_plot <- qPCR_LOD_plot +
  theme(aspect.ratio = 1.2)

LOD_melt_curve_plot <- LOD_melt_curve_plot +
  theme(aspect.ratio = 1.2)
combined_plot <- ggarrange(
  qPCR_LOD_plot, LOD_melt_curve_plot,
  ncol = 1, nrow = 2,
  common.legend = TRUE, legend = "bottom"
)

print(combined_plot)
ggsave("Combo_LOD_melt_curve_plot.png", plot = combined_plot, width = 10, height = 6, dpi = 400)
