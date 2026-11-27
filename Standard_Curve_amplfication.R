# Load libraries
library("tidyverse")
library("ggplot2")


#Read the CSV file
qPCR_efficiency <- read.csv("bbk_efficiency_12&13_2024-06-04 -Finalized_Quantification Amplification Results.csv")
head (qPCR_efficiency)

# Define the sample names in the correct order
sample_names <- c(
  "10 ng",
  "1 ng",   # 1 × 10^-9 g
  "0.1 ng", # 1 × 10^-10 g
  "0.01 ng",# 1 × 10^-11 g
  "1 pg",   # 1 × 10^-12 g
  "0.1 pg", # 1 × 10^-13 g
  "0.01 pg",# 1 × 10^-14 g
  "1 fg",   # 1 × 10^-15 g
  "0.1 fg", # 1 × 10^-16 g
  "0.01 fg",# 1 × 10^-17 g
  "NTC"
)

# Pivot the data and create new variables
qPCR_efficiency_long <- qPCR_efficiency %>%
  pivot_longer(
    cols = -Cycle, 
    names_to = "PlatePosition", 
    values_to = "Fluorescence"
  ) %>%
  # Create a Primer column based on the first letter of the plate position
  mutate(
    Primer = case_when(
      str_starts(PlatePosition, "A") ~ "KEC12",
      str_starts(PlatePosition, "B") ~ "KEC12",
      str_starts(PlatePosition, "C") ~ "KEC12",
      str_starts(PlatePosition, "E") ~ "KEC13",
      str_starts(PlatePosition, "F") ~ "KEC13",
      str_starts(PlatePosition, "G") ~ "KEC13",
     TRUE ~ NA_character_
    ),
    # Extract the numeric part from the plate position to determine sample order
    Sample_order = as.numeric(str_extract(PlatePosition, "\\d+")),
    #Assign Sample names based on sample order
    Sample = sample_names[Sample_order]
  ) %>%
   #Convert Sample to a factor with levels = sample_names for consistent color mapping
  mutate(Sample = factor(Sample, levels = sample_names))




cbbPalette <- c("#0072B2", "#D55E00", "#CC79A7", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
tol.cbb <- c("#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255", "#332288", "#D55E00","#0072B2", "#000000")


##plot

qPCR_efficiency_plot <- ggplot(
  data = qPCR_efficiency_long,
  aes(
    x = Cycle,
    y = Fluorescence,
    color = Sample,
    group = PlatePosition  # ensures each well is its own line
  )
) +
  geom_line(linewidth = 0.7, key_glyph = "path") +
  scale_color_manual(values = tol.cbb) +
  facet_wrap(~ Primer, scales = "fixed") +
  labs(x = "Cycle", y = "RFU") +
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", linewidth = 0.7) +
  
  theme_bw(base_size = 16) +
  scale_x_continuous(breaks = seq(0, 40, 10), minor_breaks = seq(0, 40, 5)) +
  scale_y_continuous(breaks = seq(0, 1200, 200), minor_breaks = seq(0, 1200, 100)) +
  theme(
    panel.grid.major = element_line(color = "black", linetype = "dotted", linewidth = 0.2),
    panel.grid.minor = element_line(color = "black", linetype = "dotted", linewidth = 0.2),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 14),
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 14),
    strip.text       = element_text(size = 14),
    legend.key.size  = unit(1.2, "lines")   # overall legend key area
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        linewidth = 2,    # Thicker lines in the legend
        shape = NA,       # No points, show lines
        linetype = 1      # Solid line
      )
    ))
    

qPCR_efficiency_plot


# Assuming your plot is stored in the variable `plot`
ggsave("qPCR_efficiency_plot6_no_bold.png", plot = qPCR_efficiency_plot, width = 10, height = 6, dpi = 400)


# We define our threshold
threshold <- 100


# For each (Primer, Sample, well position) combination, find the crossing cycle at threshold=100

qPCR_efficiency_Cq_value <- qPCR_efficiency_long %>%
  group_by(Primer, Sample, PlatePosition) %>%  #now each well is its own group
  arrange(Cycle) %>%
  summarize(
    crossing_cycle = if (max(Fluorescence, na.rm = TRUE) < threshold) {
      NA_real_
    } else {
      approx(x = Fluorescence, y = Cycle, xout = threshold)$y
    },
    .groups = "drop"
  )

threshold <- 100

qPCR_efficiency_Cq_joined <- qPCR_efficiency_long %>%
  group_by(Primer, Sample, PlatePosition) %>%
  arrange(Cycle) %>%
  summarize(
    crossing_cycle = if (max(Fluorescence, na.rm = TRUE) < threshold) {
      NA_real_
    } else {
      approx(x = Fluorescence, y = Cycle, xout = threshold)$y
    },
    .groups = "drop_last"   # Keep grouping by Primer and Sample for the next step
  ) %>%
  mutate(avg_Cq = mean(crossing_cycle, na.rm = TRUE)) %>%
  ungroup() 
  
qPCR_efficiency_Cq_joined #used this method
#try1 # not used
threshold <- 100

qPCR_efficiency_Cq_merged <- qPCR_efficiency_long %>%
  group_by(Primer, Sample, PlatePosition) %>%
  arrange(Cycle) %>%
  summarize(
    crossing_cycle = if (max(Fluorescence, na.rm = TRUE) < threshold) {
      NA_real_
    } else {
      approx(x = Fluorescence, y = Cycle, xout = threshold)$y
    },
    .groups = "drop"  # now we have one row per replicate well
  ) %>%
  group_by(Primer, Sample) %>%
  summarize(
    # Combine replicate crossing cycles into a single comma-separated string
    crossing_cycles = paste(round(crossing_cycle, 2), collapse = ", "),
    avg_Cq          = mean(crossing_cycle, na.rm = TRUE),
    .groups         = "drop"
  )


#save the computed table (e.g., the cq_values data frame) in CSV format
write.csv(qPCR_efficiency_Cq_joined, "qPCR_efficiency_Cq_joined.csv", row.names = FALSE)


#####################standard curve melting peak
#Read the CSV file
Standard_melting_peak <- read.csv("bbk_efficiency_12&13_2024-06-04 - Finalized_ Melt Curve Derivative Results.csv", na.strings = "na")
head (Standard_melting_peak)

# Define the sample names in the correct order
sample_names <- c(
  "10 ng",
  "1 ng",   # 1 × 10^-9 g
  "0.1 ng", # 1 × 10^-10 g
  "0.01 ng",# 1 × 10^-11 g
  "1 pg",   # 1 × 10^-12 g
  "0.1 pg", # 1 × 10^-13 g
  "0.01 pg",# 1 × 10^-14 g
  "1 fg",   # 1 × 10^-15 g
  "0.1 fg", # 1 × 10^-16 g
  "0.01 fg",# 1 × 10^-17 g
  "NTC"
)

# Pivot the data and create new variables
Standard_melting_peak_long <- Standard_melting_peak %>%
  pivot_longer(
    cols = -Temperature, 
    names_to = "PlatePosition", 
    values_to = "neg_dRFU_dT"
  ) %>%
  # Create a Primer column based on the first letter of the plate position
  mutate(
    Primer = case_when(
      str_starts(PlatePosition, "A") ~ "KEC12",
      str_starts(PlatePosition, "B") ~ "KEC12",
      str_starts(PlatePosition, "C") ~ "KEC12",
      str_starts(PlatePosition, "E") ~ "KEC13",
      str_starts(PlatePosition, "F") ~ "KEC13",
      str_starts(PlatePosition, "G") ~ "KEC13",
      TRUE ~ NA_character_
    ),
    # Extract the numeric part from the plate position to determine sample order
    Sample_order = as.numeric(str_extract(PlatePosition, "\\d+")),
    #Assign Sample names based on sample order
    Sample = sample_names[Sample_order]
  ) %>%
  #Convert Sample to a factor with levels = sample_names for consistent color mapping
  mutate(Sample = factor(Sample, levels = sample_names))


Standard_melting_peak_plot <- ggplot(Standard_melting_peak_long, aes(x = Temperature, y = neg_dRFU_dT , color = Sample)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(tol.cbb)) +
  facet_wrap(~ Primer, scales = "fixed") +
  labs(
    x = "Temperature, Celsius)",
    y = "-d(RFU)/dT"
  ) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "black", size = 0.8) +
  theme_bw(base_size = 14) +  # start with a clean base
  # Make the grid lines dotted and black
  theme(
    panel.grid.major = element_line(color = "black", linetype = "dotted", size = 0.5),
    panel.grid.minor = element_line(color = "black", linetype = "dotted", size = 0.3),
    
    # Draw a black border around the plotting area
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    
    # Make axis and legend text bold
    axis.title = element_text(face = "bold", size = 14),
    axis.text  = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.text  = element_text(face = "bold", size = 10),
    
    # Make facet strip labels bold
    strip.text = element_text(face = "bold", size = 14),
    
    # Increase legend key size
    legend.key.size = unit(2, "lines")
  )

Standard_melting_peak_plot


 ##try2

Standard_melting_peak_plot <- ggplot(Standard_melting_peak_long, aes(x = Temperature, y = neg_dRFU_dT , color = Sample)) +
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

Standard_melting_peak_plot
# Assuming your plot is stored in the variable `plot`
ggsave("Standard_melting_peak_plot2.png", plot = Standard_melting_peak_plot, width = 10, height = 6, dpi = 400)





###########STANDARD EQUATION CURVE################
# 1) Load libraries
library(tidyverse)
library(ggplot2)
library(ggpmisc)
library(ggpp)


# 2) Read the CSV
#    Adjust the file path as needed
Standard_eq <- read.csv("Standard_eq_excel.csv", na.strings = "na")

# Optional: Peek at the data
head(Standard_eq)
str(Standard_eq)
Standard_eq$KEC13_Average_Cq <- as.numeric(Standard_eq$KEC13_Average_Cq)
colnames(Standard_eq)


Standard_eq <- Standard_eq %>%
  rename(
    LOG10_Genomic_Equivalence = "LOG10.Genomic_Equivalence." ,
    KEC12 = "KEC12_Average_Cq" ,
    KEC13 = "KEC13_Average_Cq" 
    # etc. if needed
  )

Standard_eq_long <- Standard_eq %>%
  pivot_longer(
    cols = c("KEC12", "KEC13"),
    names_to = "Primer",
    values_to = "Cq"
  ) %>%
  filter(!is.na(Cq))  # drop rows where the Cq is missing

# recode names first:
#df_long <- df_long %>%
 # mutate(Primer = recode(Primer,
                       #  "KEC12_Average_Cq" = "KEC12",
                      #   "KEC13_Average_Cq" = "KEC13"
 # ))

# 2) Plot
Standard_curve_plot<- ggplot(Standard_eq_long, aes(x = LOG10_Genomic_Equivalence, y = Cq, color = Primer)) +
  geom_point(size = 4) +
  
  # Draw regression lines
  stat_poly_line(method = "lm") +
  scale_color_manual(values = c("KEC12" = "#E76F51","KEC13" = "#2A9D8F")) +
  
  # Basic axis labels
  labs(
    title = NULL,
    x = "Log10 (Genomic Equivalence)",
    y = "Average Cq",
    color = "Primer"    # legend title
  ) +
  
  theme_grey(base_size = 16) +
  
  # Facet by primer
  facet_wrap(~ Primer) +
  theme(
    legend.position = "top",
    legend.title     = element_text(size = 14),
    legend.text      = element_text(size = 14),
    strip.text       = element_text(size = 14),
    axis.text.x = element_text(size = 14),  # Adjust the x-axis text size
    axis.text.y = element_text(size = 16),   # Adjust the y-axis text size
    axis.title = element_text(size = 16)
  )

Standard_curve_plot

# Assuming plot is stored in the variable `Standard_curve_plot`
ggsave("Standard_curve_plot5_no_bold.png", plot = Standard_curve_plot, width = 10, height = 6.5, dpi = 400)

###For table:::
library(dplyr)

Standard_Table <- Standard_eq_long %>%
  group_by(Primer) %>%
  do({
    mod <- lm(Cq ~ LOG10_Genomic_Equivalence, data = .)
    slope <- coef(mod)[2]
    intercept <- coef(mod)[1]
    r_squared <- summary(mod)$r.squared
    
    # Calculate efficiency
    efficiency <- (10^(-1/slope) - 1) * 100
    
    # Return a small data frame with the results
    data.frame(
      Primer = unique(.$Primer),
      slope = slope,
      intercept = intercept,
      r2 = r_squared,
      efficiency = efficiency
    )
  }) %>%
  ungroup()

Standard_Table


#save the computed table in CSV format
write.csv(Standard_Table, "Standard_Table.csv", row.names = FALSE)


