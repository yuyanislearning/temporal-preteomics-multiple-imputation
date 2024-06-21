# Load necessary libraries
install.packages("viridis")
library(viridis)
install.packages("reshape2")
library(reshape2)
library(tidyverse)
library(ggplot2)
library(dplyr)
install.packages("gridExtra")
library(gridExtra)
install.packages("patchwork")
library(patchwork)
library(jsonlite)

average_data = read_csv('/content/iso_complexome_turnovers.csv')
average_data_ctrl = read_csv("/content/ctrl_turnovers_for_heterocomplexes_shared_with_iso.csv")

# filter out Complexes in average_data  that are no in average_data_ctrl
average_data_filt <- average_data %>%
  filter(Complex %in% unique(average_data_ctrl$Complex))

#check to ensure all complexes are the same in ctrl and iso
print(unique(average_data_filt$Complex) == unique(average_data_ctrl$Complex))


# Ensure that "Interactors' is a factor and its levels are in the same order as the dataframe
average_data_filt$Interactors <- factor(average_data_filt$Interactors, levels = unique(average_data_filt$Interactors))

# Now, plot with ggplot2
iso_complex_plot <- ggplot(average_data_filt, aes(x = AvgTurnover, y = Interactors, color = imp)) +
  geom_point(alpha=0.5, size = 1.5) +
  scale_x_continuous(limits = c(-6, 3), breaks = seq(-7, 4, 0.5)) +
  geom_errorbar(aes(xmin = AvgTurnover - SEM, xmax = AvgTurnover + SEM), width=0.2) +
  scale_color_manual(values = c('before' = 'blue', 'after' = 'orange')) +
  labs(x = "Average Turnover Rate (log2k)", y = "Protein Interactors", color = "Imputation Status") +
  ggtitle("ISO") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size =4.5, angle = 0, hjust = 1)  # Adjust text angle and horizontal justification as needed
  )

iso_complex_plot

# Ensure that "Interactors' is a factor and its levels are in the same order as the dataframe
average_data_ctrl$Interactors <- factor(average_data_ctrl$Interactors, levels = unique(average_data_ctrl$Interactors))

# Now, plot with ggplot2
ctrl_complex_plot <- ggplot(average_data_ctrl , aes(x = AvgTurnover, y = Interactors, color = imp)) +
  geom_point(alpha=0.5, size = 1.5) +
  scale_x_continuous(limits = c(-6, 3), breaks = seq(-7, 4, 0.5)) +
  geom_errorbar(aes(xmin = AvgTurnover - SEM, xmax = AvgTurnover + SEM), width=0.2) +
  scale_color_manual(values = c('before' = 'blue', 'after' = 'orange')) +
  labs(x = "Average Turnover Rate (log2k)", y = "Protein Interactors") +
  ggtitle("CTRL") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(size =4.5, angle = 0, hjust = 1)  # Adjust text angle and horizontal justification as needed
  )

ctrl_complex_plot

colors = c('#629da3',
 '#5c272b',
 '#cf2467',
 '#20d62c',
 '#2e7edc',
 '#cb9da1',
 '#cddb28',
 '#2a2375',
 '#d973da',
 '#2229cf',
 '#a02d22',
 '#72d722',
 '#18c08c',
 '#257a1d',
 '#7a811e',
 '#8f8de2',
 '#dd7368',
 '#d8d2d9',
 '#25638f',
 '#e1dc7b',
 '#751d8b',
 '#9be190',
 '#398f55',
 '#7130de',
 '#29d8d5',
 '#c823e2',
 '#865f5e',
 '#e03322',
 '#54e47e',
 '#d48f20',
 '#99a959',
 '#da2fa7',
 '#81dada',
 '#1c2f2e',
 '#8c5da7')

length(colors)


# Create the annotation data frame
annotation_df <- average_data_ctrl %>%
  select(Interactors, Complex) %>%
  distinct()

# Create a base ggplot object for the annotation
annotation_plot <- ggplot(annotation_df, aes(x = factor(1), y = Interactors, fill = Complex)) +
  geom_tile() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(fill = "Annotation")

# Combine the plots using patchwork
combined_plot <- iso_complex_plot + ctrl_complex_plot + annotation_plot

# Use plot_layout to specify the widths of the plots
combined_plot <- combined_plot +
                 plot_layout(widths = c(10,10, 1))

# Draw the combined plot
combined_plot
#ggsave("combined_complex_plot.pdf", combined_plot)

annotation_plot <- ggplot(annotation_df, aes(x = factor(1), y = Interactors, fill = Complex)) +
  geom_tile() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(fill = "Annotation")

annotation_plot
#ggsave("annot.pdf", annotation_plot)


# filter average_data for just complexes of interest
average_data_filtered <- average_data %>%
  filter(Complex %in% c("CPX-3041", "CPX-5823", "CPX-3009", "CPX-5152"))

# Ensure that "Interactors' is a factor and its levels are in the same order as the dataframe
average_data_filtered$Interactors <- factor(average_data_filtered$Interactors, levels = unique(average_data_filtered$Interactors))

# Now, plot with ggplot2
iso_plot <- ggplot(average_data_filtered, aes(x = AvgTurnover, y = Interactors, color = imp)) +
  geom_point(alpha=0.5) +
  # scale_x_continuous(limits = c(-5.5, 0), breaks = seq(-5.5, 0, 0.5)) +
  geom_errorbar(aes(xmin = AvgTurnover - SEM, xmax = AvgTurnover + SEM), width=0.2) +
  scale_color_manual(values = c('before' = 'blue', 'after' = 'orange')) +
  labs(x = "Average Turnover Rate (log2k)", y = "Protein Interactors") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(angle = 0, hjust = 1)  # Adjust text angle and horizontal justification as needed
  )

###### CTRL

# filter for complexes of interest
average_data_filtered_ctrl <- average_data_ctrl %>%
  filter(Complex %in% c("CPX-3041", "CPX-5823", "CPX-3009", "CPX-5152"))

# Ensure that "Interactors' is a factor and its levels are in the same order as the dataframe
average_data_filtered_ctrl$Interactors <- factor(average_data_filtered_ctrl$Interactors, levels = unique(average_data_filtered_ctrl$Interactors))

# Now, plot with ggplot2
ctrl_plot <- ggplot(average_data_filtered_ctrl, aes(x = AvgTurnover, y = Interactors, color = imp)) +
  geom_point(alpha=0.5) +
  # scale_x_continuous(limits = c(-5.5, 0), breaks = seq(-5.5, 0, 0.5)) +
  geom_errorbar(aes(xmin = AvgTurnover - SEM, xmax = AvgTurnover + SEM), width=0.2) +
  scale_color_manual(values = c('before' = "blue", 'after' = 'orange')) +
  labs(x = "Average Turnover Rate (log2k)", y = "Protein Interactors", color = "Imputation") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(angle = 0, hjust = 1)  # Adjust text angle and horizontal justification as needed
  )

  colors = c("#e41a1c", "#377eb8", "#4daf4a",  '#7a811e')

# Create the annotation data frame
annotation_df <- average_data_filtered_ctrl %>%
  select(Interactors, Complex) %>%
  distinct()

# Create a base ggplot object for the annotation
annotation_plot <- ggplot(annotation_df, aes(x = factor(1), y = Interactors, fill = Complex)) +
  geom_tile() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(fill = "Complex")

# Combine the plots using patchwork
combined_plot <- ctrl_plot + iso_plot + annotation_plot

# Use plot_layout to specify the widths of the plots
combined_plot <- combined_plot +
                 plot_layout(widths = c(50, 50, 1))

# Draw the combined plot
combined_plot

ggsave("zoom_in_V2.pdf", combined_plot)

json_data <- read_json('/content/coherence_scores.json')

# Convert the JSON data to a data frame
ctrl <- data.frame(json_data$coherence_ctrl)
iso <- data.frame(json_data$coherence_iso)

ctrl <- data.frame(t(ctrl))
iso <- data.frame(t(iso))
colnames(ctrl) <- "pval_Ctrl"
colnames(iso) <- "pval_Iso"
ctrl_iso_coherence <- cbind(ctrl, iso)


# First, reset the row names to be a column in the dataframe
ctrl_iso_coherence$ID <- rownames(ctrl_iso_coherence)

# Now, melt the data frame using the new ID column
long_data <- melt(ctrl_iso_coherence, id.vars = 'ID')

# Rename the 'variable' and 'value' columns generated by melt for clarity
names(long_data) <- c('ID', 'Condition', 'Score')
long_data$Score <- log10(long_data$Score)

# Create the heatmap
coherence_plt <- ggplot(long_data, aes(x = Condition, y = ID, fill = Score)) +
  geom_tile(color = "white") +  # Add a border to make the labels easier to read
  geom_text(aes(label = round(Score, 8)), size = 3) +
  scale_fill_gradient(low = "lightyellow", high = "darkgreen") +  # Use a red color gradient
  labs(x = "", y = "", fill = "log10(p-value) w/ ANOVA") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(angle = 0, hjust = 1)  # Adjust text angle and horizontal justification as needed
  )

coherence_plt

ggsave("coherence_scores.pdf", coherence_plt )