library(ggplot2)

kmer_unique <- read.csv("KEC_table_final.csv")
head(kmer_unique)


###plot; Kmer size and the unique sequence output

kmer_vs_uniq <-ggplot(kmer_unique, aes(x =  K.mers_size, y = No..of.unique.regions.)) +
  geom_col(fill = "#1f77b4", width = 0.7, alpha = 0.8) +  # Use a professional blue color
  geom_text(aes(label = No..of.unique.regions.),          # Add labels above bars
            vjust = -0.5, size = 4, color = "black", fontface = "bold") +
  labs(
    x = "K-mer Size",
    y = "Number of Unique Regions",
  ) +
  theme_grey(base_size = 16) +                        # Increase base font size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),      # Rotate x-axis labels
    axis.text.y = element_text(size = 12),                            # Adjust y-axis text size
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),  # Bold x-axis title
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),  # Bold y-axis title
    panel.grid.major.x = element_blank(),                             # Remove vertical grid lines
    panel.grid.minor = element_blank(),                               # Remove minor grid lines
    plot.caption = element_text(size = 10, color = "gray50", hjust = 1),  # Add caption
    plot.margin = unit(c(1, 1, 1, 1), "cm")                           # Adjust plot margins
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +          # Add space above bars
  scale_x_discrete(labels = function(x) gsub("kec-", "", x))          # Clean up x-axis labels
# Print the plot
kmer_vs_uniq


colors <- c(
  "#c6dbef", "#9ecae1", "#3182bd", "#08519c", "#6baed6",
  "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#005a32",
  "#fdd0a2", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#8c2d04",
  "#dadaeb", "#bcbddc", "#9e9ac8", "#807dba", "#6a51a3", "#4a1486",
  "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d",
  "#d9d9d9", "#bdbdbd", "#969696", "#737373", "#525252", "#252525"
)


# Create the plot
plot_kec_vs_uniqueseq <- ggplot(kmer_unique, aes(x = K.mers_size, y = No..of.unique.regions.)) +
  geom_col(aes(fill = K.mers_size), width = 0.7, alpha = 0.8) +  # Use custom palette for fill
  geom_text(aes(label = No..of.unique.regions.),  
            vjust = -0.5, size = 3, fontface = "bold") +  # Smaller labels for more KEC sizes
  scale_fill_manual(values = colors) +  # Use custom color palette
  scale_color_identity() +  # Use the colors directly from the data for labels
  labs(
    x = "K-mer Size",
    y = "Number of Unique Sequences"
  ) +
  theme_gray(base_size = 14) +  # Clean and minimal theme
  theme(
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),  # Rotate x-axis labels 90 degrees
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12, face= "bold"),
    axis.title.y = element_text(size = 12, face= "bold"),
    panel.grid.major.x = element_blank(),  # Remove vertical grid lines
    panel.grid.minor = element_blank(),    # Remove minor grid lines
    legend.position = "none"              # Remove legend
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  # Add space above bars
plot_kec_vs_uniqueseq


# Assuming your plot is stored in the variable `plot`
ggsave("Figures/unique_sequences_plot.png", plot = plot_kec_vs_uniqueseq, width = 10, height = 6, dpi = 400)




##Figure : BLAST hits; DH44 coding regions  vs KEC regions 
#Load the data
Blast_results <- read.csv("Final_BLAST.csv")

str(results)
library("tidyverse")
filtered.results <- Blast_results %>%
  subset(percent_id == 100 & alignment_length >= 150) %>%
  group_by(Unique_genes, alignment_length) %>%
  count() 

# these are unique single copy regions in DH44 with a kmer size of 13
filtered.results[filtered.results$n == 1,] %>%
  print(n = 100)

single_blast_hits <-ggplot(filtered.results, aes(x = reorder(Unique_genes, -n), y = n)) + 
  geom_point() +
  theme_classic() +
  ylab("Number of exact blast hits") +
  xlab("KEC identified regions (kmer = 16)") +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust=1, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"), 
        axis.title.x = element_text(size = 14, face= "bold"),
        axis.title.y = element_text(size = 14, face= "bold"))+
  geom_vline(xintercept = 7.5, color = "red", linetype = "dashed",size = 0.8)

single_blast_hits
# Assuming your plot is stored in the variable `plot`
ggsave("Figures/single_blast_hits.png", plot = single_blast_hits, width = 10, height = 6, dpi = 400)


