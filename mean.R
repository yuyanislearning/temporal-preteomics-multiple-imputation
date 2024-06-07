#!/usr/bin/env Rscript

# Load necessary libraries
library(tidyr)
library(glue)
library(data.table)
library(missMethods)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set default values if arguments are not provided
if (length(args) != 2) {
  stop("Usage: Rscript mean.R <n_masked> <missing_data>", call. = FALSE)
}

# Parse arguments
n_masked <- as.numeric(args[1])
dat1_file <- args[2]

# Print the inputs for verification
cat("n_masked:", n_masked, "\n")
cat("dat1_file:", dat1_file, "\n")

# Load in data 
dat1 <- fread(dat1_file)

# Filter for columns ID, t, A0
dat1 <- dat1[, .(ID, t, A0)]

# Aggregate A0 by t and ID, calculating the mean
unique_A0 <- aggregate(dat1$A0, by = list(dat1$t, dat1$ID), mean)
colnames(unique_A0) <- c("t", "ID", "A0")

# Convert long format to wide format
data <- spread(unique_A0, t, A0)

# Remove negative halflife
data[data < 0] <- NA

# This vector will contain the number of missing values
cnt_na <- apply(data, 1, function(z) sum(is.na(z)))

data_filt <- data
names(data_filt) <- c("ID", "d0", "d1", "d3", "d5", "d7", "d10", "d14") #nolint 

# Store IDs as a vector
id <- data_filt$ID

# Remove ID column for imputation
data_filt1 <- subset(data_filt, select = -c(ID))

# Convert all columns to numeric
data_filt1[] <- lapply(data_filt1, function(x) as.numeric(as.character(x)))

# Calculate the percentage of missing values
mis_per <- sum(is.na(data_filt1)) / prod(dim(data_filt1))
print(glue("{mis_per*100}% of the data is missing"))

#### Perform mean imputation
mean_impute <- impute_mean(data_filt1)

#file output directory 
directory <- unlist(strsplit(dat1_file, "/"))[1]

write.csv(mean_impute, glue("{directory}/si_mean_{n_masked}.csv"))
