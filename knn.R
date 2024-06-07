#!/usr/bin/env Rscript

# Load necessary libraries
library(tidyr)
library(glue)
library(bnstruct)
library(data.table)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set default values if arguments are not provided
if (length(args) != 2) {
  stop("Usage: Rscript mice.R <n_masked> <missing_data>", call. = FALSE)
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
data[data < 0] <- NA # Remove negative halflife

# This vector will contain the number of missing values
cnt_na <- apply(data, 1, function(z) sum(is.na(z)))

data_filt <- data
names(data_filt) <- c("ID", "d0", "d1", "d3", "d5", "d7", "d10", "d14") # nolint 

# Store IDs as a vector
id <- data_filt$ID

# Remove ID column for imputation
data_filt1 <- subset(data_filt, select = -c(ID))

# Convert all columns to numeric
data_filt1[] <- lapply(data_filt1, function(x) as.numeric(as.character(x)))

#transpose the df so we can impute the time points since knn is doing by rows
data_filt1 <- t(data_filt1)
print(colnames(data_filt1))
print(rownames(data_filt1))

# Ensure data_filt1 is a matrix
data_matrix <- as.matrix(data_filt1)

# Calculate the percentage of missing values
mis_per <- sum(is.na(data_matrix)) / prod(dim(data_matrix))
print(glue("{mis_per*100}% of the data is missing"))

# Number of columns in data_matrix
num_cols <- ncol(data_matrix)


# Perform KNN imputation 
knn_impute <- knn.impute(data_matrix, k = 10, cat.var = 1:num_cols, to.impute = 1:nrow(data_matrix)) # nolint 

# Transpose the imputed matrix
knn_impute <- t(knn_impute)
print(colnames(knn_impute))

# File output directory
directory <- unlist(strsplit(dat1_file, "/"))[1]

write.csv(knn_impute, glue("{directory}/knn_{n_masked}.csv"))
