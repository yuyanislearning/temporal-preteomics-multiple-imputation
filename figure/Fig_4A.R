library(tidyverse)
library(ggplot2)
library(dplyr)
install.packages("gridExtra")
library(gridExtra)
install.packages("patchwork")
library(patchwork)

strains = c('dba','fvb','cej','c57', 'balbc','aj')
conditions = c('iso', 'ctrl')
imputations = c('before', 'after')
pathway = '/content/'

all_data <- list()
for(strain in strains){
  for(cond in conditions){
    for(imp in imputations){
      dat = read_csv(paste0(pathway, paste(cond, strain, imp, sep='_'), '.csv'))
      dat$cond = cond
      dat$strain = strain
      dat$imp = imp
      all_data[[length(all_data) + 1]] <- dat
    }
  }
}
combined_data <- do.call(rbind, all_data)

iso_data <- combined_data %>%
  filter(cond == 'iso')

ctrl_data <- combined_data %>%
  filter(cond == 'ctrl')

########## ISO ####################

iso_data$log2_median.K <- log2(iso_data$median.K)
# Summarize the data while handling proteins with a single or no observation
average_data <- iso_data %>%
  group_by(UniProt, imp) %>%
  summarise(
    AvgTurnover = mean(log2_median.K, na.rm = TRUE),  # Calculate average turnover rate
    SEM = sd(log2_median.K, na.rm = FALSE) / sqrt(n()),  # Calculate SEM, will be NA if n < 2
    Count = n(),  # Count the number of observations
    .groups = 'drop'  # Drop the grouping
  ) %>%
  # Replace NA SEM values with 0
  mutate(SEM = replace_na(SEM, 0))

# Pivot the data to a wide format
df_wide <- average_data %>%
  pivot_wider(
    names_from = imp,
    values_from = c(AvgTurnover, SEM, Count),
    names_sep = "_"  # This will create column names like AvgTurnover_before, AvgTurnover_after, etc.
  )

# Replace NA with 0 and calculate the before_after_diff
df_wide <- df_wide %>%
  mutate(
    AvgTurnover_before = replace_na(AvgTurnover_before, 0),
    AvgTurnover_after = replace_na(AvgTurnover_after, 0),
    before_after_diff = AvgTurnover_after - AvgTurnover_before
  )

# Prepare a lookup dataframe with just UniProt and before_after_diff
lookup_df <- df_wide %>%
  select(UniProt, before_after_diff)


# Join the before_after_diff column to average_data
average_data <- average_data %>%
  left_join(lookup_df, by = "UniProt")

# Now we create the scatterplot with error bars only when SEM is available
iso_plot <- ggplot(average_data, aes(x = AvgTurnover, y = reorder(UniProt, before_after_diff), color = imp)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(limits = c(-7, 9), breaks = seq(-7, 9, 1)) +
  geom_errorbar(aes(xmin = AvgTurnover - SEM, xmax = AvgTurnover + SEM), width = 0.2, na.rm = FALSE) +
  scale_color_manual(values = c('before' = 'blue', 'after' = 'orange')) +
  labs(x = "Average Turnover Rate", y = "Unique Proteins") +
  #theme_minimal() +
  ggtitle("ISO") +
  theme(
    legend.position = "bottom",
    # axis.text.y = element_text(size =1, angle = 0, hjust = 1)
    axis.text.y = element_blank()
  )

iso_plot
# Save the plot if desired

ggsave("iso.pdf", iso_plot)

########## CTRL ####################

ctrl_data$log2_median.K <- log2(ctrl_data$median.K)
# Summarize the data while handling proteins with a single or no observation
average_data <- ctrl_data %>%
  group_by(UniProt, imp) %>%
  summarise(
    AvgTurnover = mean(log2_median.K, na.rm = TRUE),  # Calculate average turnover rate
    SEM = sd(log2_median.K, na.rm = FALSE) / sqrt(n()),  # Calculate SEM, will be NA if n() < 2
    Count = n(),  # Count the number of observations
    .groups = 'drop'  # Drop the grouping
  ) %>%
  # Replace NA SEM values with 0
  mutate(SEM = replace_na(SEM, 0))

# Pivot the data to a wide format
df_wide <- average_data %>%
  pivot_wider(
    names_from = imp,
    values_from = c(AvgTurnover, SEM, Count),
    names_sep = "_"  # This will create column names like AvgTurnover_before, AvgTurnover_after, etc.
  )

# Replace NA with 0 and calculate the before_after_diff
df_wide <- df_wide %>%
  mutate(
    AvgTurnover_before = replace_na(AvgTurnover_before, 0),
    AvgTurnover_after = replace_na(AvgTurnover_after, 0),
    before_after_diff = AvgTurnover_after - AvgTurnover_before
  )

# Prepare a lookup dataframe with just UniProt and before_after_diff
lookup_df <- df_wide %>%
  select(UniProt, before_after_diff)


# Join the before_after_diff column to average_data
average_data <- average_data %>%
  left_join(lookup_df, by = "UniProt")

# Now we create the scatterplot with error bars only when SEM is available
ctrl_plot <- ggplot(average_data, aes(x = AvgTurnover, y = reorder(UniProt, before_after_diff), color = imp)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(limits = c(-7, 9), breaks = seq(-7, 9, 1)) +
  geom_errorbar(aes(xmin = AvgTurnover - SEM, xmax = AvgTurnover + SEM), width = 0.2, na.rm = FALSE) +
  scale_color_manual(values = c('before' = 'blue', 'after' = 'orange')) +
  labs(x = "Average Turnover Rate") +
  ggtitle("CTRL") +
  # theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank()
  )
ctrl_plot
# Save the plot if desired
#ggsave("protein_turnover_scatterplot_ctrl.png", width = 10, height = 12, dpi = 300)


ggsave("ctrl.pdf", ctrl_plot)

# Combine the plots using patchwork
combined_plot <- iso_plot + ctrl_plot

# Use plot_layout to specify the widths of the plots
combined_plot <- combined_plot +
                 plot_layout(widths = c(20, 20))

# Draw the combined plot
combined_plot

# prompt: set range and ticks in g plot with scale x contiuous

# Set the x-axis limits
ggplot(data, aes(x = x, y = y)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 10))


ggsave("protein_turnover_scatterplot.pdf", width = 10, height = 12, dpi = 300)

