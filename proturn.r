library(tidyverse)
library(dplyr)
library(glue)
rm(list = ls(all = T))  # Remove all objects from the current workspace

# Define the directory and file locations
dir <- "proturn_output_aj_ctrl"
hl.out_location <- glue("{dir}/hl.out")
hl.data.out_location <- glue("{dir}/dmi_np_1.csv")
ids <- glue("{dir}/IDs_1.csv")

# Read in IDs from the specified CSV file
id_df <- read.csv(ids, header = FALSE, stringsAsFactors = FALSE, skip = 1)

output_file <- "ProTurn_Output_aj_ctrl_imp1.txt"

Refit <- TRUE  # Flag to determine whether to refit the kinetic curve
# Optimization method used, Brent or Nelder-Mead
optim_method <- "Brent"

# Define the header for the output file and write it
oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
write(oput, file = output_file, append = T) 

# Read the input data files
hl.out <- read.table(hl.out_location, header = TRUE, fill = TRUE, skip = 0)  # Read hl.out data
dp_init <- read.csv(hl.data.out_location)  # Read the initial data

# Filter dp_init to include only rows where .imp column is 1
dp1 <- dp_init[dp_init[".imp"] == "1", ]

# Add a new column "ID" to dp1 using the second column from id_df
dp1 <- cbind(dp1, ID = id_df$V2)

# Rename columns for better readability
colnames(dp1)[colnames(dp1) == "d0"] <- "0"
colnames(dp1)[colnames(dp1) == "d1"] <- "1"
colnames(dp1)[colnames(dp1) == "d3"] <- "3"
colnames(dp1)[colnames(dp1) == "d5"] <- "5"
colnames(dp1)[colnames(dp1) == "d7"] <- "7"
colnames(dp1)[colnames(dp1) == "d10"] <- "10"
colnames(dp1)[colnames(dp1) == "d14"] <- "14"

# Remove unnecessary columns
dp1 <- dp1[ , !(names(dp1) %in% c("X", ".imp", ".id"))]

# Convert dp1 to long format
dp <- dp1 %>%
  pivot_longer(cols = c(`0`, `1`, `3`, `5`, `7`, `10`, `14`), 
               names_to = "t", 
               values_to = "A0") %>%
  select(ID, t, A0)

dp$t <- as.numeric(as.character(dp$t))  # Convert t to numeric

# Remove specific columns from hl.out and save the truncated data in dt
dt <- hl.out[,1:11];dt <- dt[,-c(7)]
# Save the original hl.out data
hl.out_truth <- hl.out

print(head(dp))
print(head(dt))

# Function to calculate Fractional Synthesis (FS)
Calculate_FS <- function(x) {
  A0. <- dt$a[c]
  Ainf. <- dt$a[c] * (1 - dt$pss[c]) ^ dt$N[c]
  FS. <- (x - A0.) / (Ainf. - A0.)
  return(FS.)
}

# Function for refitting the kinetic curve
Refitting_Function <- function(ki) {
  current_k <<- ki  # Set the current k globally for use in the core function
  Refitting_Predicted <- sapply(ds$t, Refitting_Function_Core)
  Refitting_R2 <- 1 - (sum((ds$FS - Refitting_Predicted)^2)) / (sum((ds$FS - mean(ds$FS))^2))
  return(1 - Refitting_R2)
}

# Core function for reoptimization with KL
Refitting_Function_Core <- function(x) {
  z <- 0
  N <- dt$N[c]
  kp <- dt$kp[c]
  pss <- dt$pss[c]
  a <- dt$a[c]
  for (n in 0:N) {
    b <- factorial(N) / (factorial(n) * factorial(N - n)) * (1 - pss)^(N - n) * (pss)^n
    bp <- (current_k / (current_k - n * kp)) * b
    y <- (bp * exp(-n * kp * x) + exp(-current_k * x) * (1 / (N + 1) - bp))
    z <- y + z
  }
  final <- (z - 1) / ((1 - pss)^N - 1)
  return(final)
}

# Function to model the kinetic curve using optimized k
Refitted_Model <- function(x) {
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c]) / (factorial(n) * factorial(dt$N[c] - n)) * (1 - dt$pss[c])^(dt$N[c] - n) * (dt$pss[c])^n
    bp <- (Optimize$par / (Optimize$par - n * dt$kp[c])) * b
    y <- dt$a[c] * (bp * exp(-n * dt$kp[c] * x) + exp(-Optimize$par * x) * (1 / (dt$N[c] + 1) - bp))
    z <- y + z
  }
  return(z)
}

# Function to model the kinetic curve without optimization
Model <- function(t) {
  z <- 0
  for (n in 0:dt$N[c]) {
    k <- dt$k[c]
    N <- dt$N[c]
    kp <- dt$kp[c]
    pss <- dt$pss[c]
    a <- dt$a[c]
    
    b <- factorial(N) / (factorial(n) * factorial(N - n)) * (1 - pss)^(N - n) * (pss)^n
    bp <- (k / (k - n * kp)) * b
    y <- a * (bp * exp(-n * kp * t) + exp(-k * t) * (1 / (N + 1) - bp))
    z <- y + z
  }
  return(z)
}

# Function to calculate dk
Calculate_dk <- function(x) {
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c]) / (factorial(n) * factorial(dt$N[c] - n)) * (1 - dt$pss[c])^(dt$N[c] - n) * (dt$pss[c])^n
    bp <- (Optimize$par / (Optimize$par - n * dt$kp[c])) * b
    y <- dt$a[c] * ((n * dt$kp[c]) / (Optimize$par * (Optimize$par - n * dt$kp[c])) * bp * (exp(-Optimize$par * x) - exp(-n * dt$kp[c] * x)) - x * (1 / (dt$N[c] + 1) - bp) * exp(-Optimize$par * x))
    z <- y + z
  }
  return(se_new / z)
}

# Refitting process if Refit is TRUE
if (Refit == TRUE) {
  for (c in 1:nrow(dt)) {
    print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100, 2), "% done."))
    id <- dt$ID[c] 
    ds <- dp[ which(dp$ID == id),]  # Subsetting for the particular ID
    
    if (nrow(ds) != 0) {  # Ensure the subset is not empty
      # Calculate Fractional Synthesis
      FS <- sapply(ds$A0, Calculate_FS)
      ds <- cbind(ds, FS)
      
      # Optimization function to find the best k value
      Optimize <- optim(0.3, Refitting_Function)
      
      # Calculate the refitted model and R2
      Refitted_Predicted <- sapply(ds$t, Refitted_Model)
      Refitted_R2 <- 1 - (sum((ds$A0 - Refitted_Predicted)^2)) / (sum((ds$A0 - mean(ds$A0))^2))
      
      # Calculate the standard error
      se_new <- (sum((ds$A0 - Refitted_Predicted)^2) / 6)^0.5  # Assuming 7 time points
      
      # Calculate dk
      dkTable <- sapply(ds$t, Calculate_dk)
      dk <- min(abs(dkTable))
      
      oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c],dt$DP[c],dt$z[c],"0",round(se_new,4),round(dt$a[c],3), round(dt$pss[c],4),dt$kp[c], dt$N[c],round(Optimize$par,3), round(dk,3),  round(Refitted_R2,3),      sep = "\t")
      write(oput, file=output_file, append=T)
        }
  }
  
}
