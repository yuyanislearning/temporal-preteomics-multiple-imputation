#!/usr/bin/env Rscript

# Load necessary libraries
library(tidyr)
library(glue)
library(mice)
library(data.table)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set default values if arguments are not provided
if (length(args) != 3) {
  stop("Usage: Rscript mice.R <n_masked> <hl.out> <hl-data>", call. = FALSE)
}

# Parse arguments
n_masked <- as.numeric(args[1])
dat_file <- args[2]
dat1_file <- args[3]

# Print the inputs for verification
cat("n_masked:", n_masked, "\n")
cat("dat_file:", dat_file, "\n")
cat("dat1_file:", dat1_file, "\n")

# Load in data
dat <- fread(dat_file)
dat1 <- fread(dat1_file)

# Filter for columns ID, t, A0 
dat1 <- dat1[, .(ID, t, A0)]

# Aggregate A0 by t and ID, calculating the mean
unique_A0 <- aggregate(dat1$A0, by = list(dat1$t, dat1$ID), mean)
colnames(unique_A0) <- c("t", "ID", "A0")

# Convert long format to wide format
data <- spread(unique_A0, t, A0)

# Merge peptide A0 data with additional metadata stored in data (can merge additional metadata here if desired)
data_merge <- merge(data, dat, by = "ID")
data_merge <- data_merge[, -c(9:21)]
#data_merge <- data_merge[, -c(10,11,13,14,15)]
data <- data_merge
data[data < 0] <- NA # Remove negative halflife

# This vector will contain the number of missing values
cnt_na <- apply(data, 1, function(z) sum(is.na(z)))

data_filt <- data
names(data_filt) <- c("ID", "d0", "d1", "d3", "d5", "d7", "d10", "d14") #, "a", "N") #nolint 

# Store IDs as a vector
id <- data_filt$ID

# Remove ID column for imputation
data_filt1 <- subset(data_filt, select = -c(ID))

# Convert all columns to numeric
data_filt1[] <- lapply(data_filt1, function(x) as.numeric(as.character(x)))

# Calculate the percentage of missing values
mis_per <- sum(is.na(data_filt1)) / prod(dim(data_filt1))
print(glue("{mis_per*100}% of the data is missing"))

# Number of columns in data_filt1
num_cols <- ncol(data_filt1)

# Create a method vector with 'pmm' for each column,
# method_vector <- c(rep("norm.predict", num_cols - 2), "", "") #this is if N and a are included
method_vector <- c(rep("norm.predict"))

# Create a custom predictor matrix
predictor_matrix <- make.predictorMatrix(data_filt1)

# Exclude N and a from being imputed if included 
# predictor_matrix["N", ] <- 0
# predictor_matrix["a", ] <- 0

#### Perform multiple imputation using MICE.
# The arguments used are -->
# data_filt: The data frame to be imputed.
# m: The number of imputations to be performed
# method: imputation method for each column.
# The pmm (predictive mean matching) method 
# maxiter: The maximum number of iterations
# for the imputation algorithm (10 in this case).
mice_impute <- mice(data_filt1, m = 10, method = method_vector, maxiter = 50, predictorMatrix = predictor_matrix) # nolint

# Summary of the imputation process
summary(mice_impute)

# Output the imputed data matrix with the 10 iterations
data_all <- complete(mice_impute, "long")

#file output directory 
directory <- unlist(strsplit(dat_file, "/"))[1]

write.csv(data_all, glue("{directory}/test_{n_masked}.csv"))
#write.csv(id, glue("{directory}/IDs_{n_masked}.csv"))
