setwd('/Users/simhasankar/Documents/human_plasma_proteome_imputation_turnover/human') #nolint 
rm(list=ls())
library(tidyverse)
library(glue)

imp_files <- c('HS2_RE2/ProTurn_Output_hs2_imp1.txt', 'HS3RE_sub1/ProTurn_Output_hs3_imp1.txt', 'HS4/ProTurn_Output_hs4_imp1.txt', #nolint
              'HS5/ProTurn_Output_hs5_imp1.txt', 'HS6/ProTurn_Output_hs6_imp1.txt', 'HS7/ProTurn_Output_hs7_imp1.txt', 'HS8/ProTurn_Output_hs8_imp1.txt', #nolint
              'HS9/ProTurn_Output_hs9_imp1.txt', 'HS10/ProTurn_Output_hs10_imp1.txt', 'HS11/ProTurn_Output_hs11_imp1.txt') #nolint 

raw_files <- c('HS2_RE2/hl.out', 'HS3RE_sub1/hl.out', 'HS4/hl.out', #nolint 
              'HS5/hl.out', 'HS6/hl.out', 'HS7/hl.out', 'HS8/hl.out', #nolint 
              'HS9/hl.out', 'HS10/hl.out', 'HS11/hl.out') #nolint 

raw_pro <- list()
out_pro <- list()

for (i in 1:10) {
  raw <- read_delim(raw_files[i], delim = "\t")
  # Apply filters
  raw <- raw %>%
    filter(R2 >= 0.8 | SS < 0.05) %>%
    filter(DP > 4) %>%
    group_by(Uniprot)
  raw_pro[[i]] <- as.vector(unique(raw$Uniprot))
  out <- read_delim(imp_files[i], delim = "\t")
  out <- out %>%
    filter(R2 >= 0.8 | SS < 0.05) %>%
    filter(DP >= 2) %>%
    group_by(Uniprot)
  out_pro[[i]] <- as.vector(unique(out$Uniprot))
}

sum <- 0
for (i in 1:10) {
  print(paste(i, length(raw_pro[[i]])))
  sum <- sum + length(raw_pro[[i]])
}
print(sum)

for (i in 1:10) {
  print(paste(i, length(out_pro[[i]])))
}

#make a directory called Pre_Post_Imputation_Proteins
dir.create("Pre_Post_Imputation_Proteins")

#write outpro[[i]] to a file titled out_pro[[i]].txt
#write raw_pro[[i]] to a file titled raw_pro[[i]].txt
for (i in 1:10) {
  write.csv(raw_pro[[i]], file=glue("Pre_Post_Imputation_Proteins/raw_pro_{i}.csv"), row.names=FALSE) #nolint 
  write.csv(out_pro[[i]], file=glue("Pre_Post_Imputation_Proteins/out_pro_{i}.csv"), row.names=FALSE) #nolint 
}