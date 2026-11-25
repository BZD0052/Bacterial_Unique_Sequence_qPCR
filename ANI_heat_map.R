library(tidyverse)

# Read the ANI matrix from a CSV file.
# Assumes the first column contains row names (genome names).
ani <- read.csv("ANIm_percentage_identity_Payni.csv", header = TRUE, row.names = 1)


# Convert the data frame (matrix) into long format.
ani_long <- as.data.frame(as.table(as.matrix(ani)))
colnames(ani_long) <- c("Genome1", "Genome2", "PercentIdentity")

# Create the heatmap.
ANI_heatmap_plot <- ggplot(ani_long, aes(x = Genome1, y = Genome2, fill = PercentIdentity)) +
  geom_tile(color = "white") +                       # Add white borders between tiles
  scale_fill_gradient(low = "white", high = "Blue") +  # Set color gradient
  theme_bw(base_size = 14, base_family = "Arial") +    # Use Arial as base font
  theme(aspect.ratio = .75,
    axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),  # Increase x-axis text size and rotate
    axis.text.y = element_text(size = 12, color = "black"),                         # Increase y-axis text size
    axis.title = element_text(size = 18, face = "bold"),           # Increase axis title size
    legend.title = element_text(size = 14, face = "bold"),         # Increase legend title size
    legend.text  = element_text(size = 12)                         # Increase legend text size
  ) +
  labs(
    x = "",
    y = "",
    fill = "Percent Similarity"
  )

# Display the plot
print(ANI_heatmap_plot)



# Save the plot to a file at high resolution.
ggsave("ANI_heatmap_plot3.png", ANI_heatmap_plot, width = 10, height = 8, dpi = 400)



library(tidyverse)

# --- your code up to ani_long already ran ---

# 1) Order of codes (left→right on x, bottom→top on y)
code_order <- c("DH35","DH44_assembly","DH45","DH46","JJ1831","JJ1883","JJ190","JJ591","JJ64")

# 2) Pretty labels matching that order
label_map <- c(
  "DH35"          = "Bacillus cereus AG DH35",
  "DH44_assembly" = "P. riograndensis DH44",   # assuming DH44 (you wrote DH4$)
  "DH45"          = "Paenibacillus spp. DH45",
  "DH46"          = "P. riograndensis DH46",
  "JJ1831"        = "Paenibacillus spp. JJ1831",
  "JJ1883"        = "P. silagei JJ1883",
  "JJ190"         = "P. riograndensis JJ190",
  "JJ591"         = "P. graminis JJ591",
  "JJ64"          = "P. rhizoplanae JJ64"      # normalized spacing/punctuation
)


label_map <- c(
  "DH35"          = "italic('Bacillus A cereus')~AG~DH35",
  "DH44_assembly" = "italic('P. riograndensis')~DH44",
  "DH45"          = "italic('Paenibacillus spp.')~DH45",
  "DH46"          = "italic('P. riograndensis')~DH46",
  "JJ1831"        = "italic('Paenibacillus spp.')~JJ1831",
  "JJ1883"        = "italic('P. silagei')~JJ1883",
  "JJ190"         = "italic('P. riograndensis')~JJ190",
  "JJ591"         = "italic('P. graminis')~JJ591",
  "JJ64"          = "italic('P. rhizoplanae')~JJ64"
)

# 3) Lock factor order to codes
ani_long <- ani_long %>%
  mutate(
    Genome1 = factor(Genome1, levels = code_order),
    Genome2 = factor(Genome2, levels = code_order)
  )

# 4) Plot with axis labels remapped at the very end
ANI_heatmap_plot <- ggplot(ani_long, aes(x = Genome1, y = Genome2, fill = PercentIdentity)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "Blue") +
  scale_x_discrete(labels = function(x) parse(text = label_map[x])) +
  scale_y_discrete(labels = function(x) parse(text = label_map[x])) +
  theme_bw(base_size = 14, base_family = "Arial") +
  theme(
    aspect.ratio = .75,
    axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title  = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 14)
  ) +
  labs(x = "", y = "", fill = "Percent Similarity")

print(ANI_heatmap_plot)
ggsave("ANI_heatmap_plot4.png", ANI_heatmap_plot, width = 10, height = 8, dpi = 400)
ANI_heatmap_plot4.png



