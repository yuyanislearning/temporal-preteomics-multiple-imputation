rm(list = ls(all = T)) 

#R2_threshold <- 0.8
#SE_threshold <- 0.05 

#home_directory <- "~/Documents"
hl.out_location <- "hl.out"
#hl.data.out_location <- "hl-data_verified.out" 
hl.data.out_location <- "hs3_imp.1.out"   ############################3

output_file <- "ProTurn_Output_hs3_imp1.txt"  ######################################
#output_file <- "TestSS_SE.txt"

Refit <- TRUE           # Should the Grapher attempt to gather the data points for each peptide and refit the kinetic curve in R.
#optim_method = "Brent"  # Optimization method used, refer to ?optim for details
optim_method = "Nelder-Mead"


#oput <- paste("ID","Uniprot","Peptide","z","R2","SS","SE","k","dk","DP","a","pss","kp","N", sep = "\t")
oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
write(oput, file=output_file, append=T)

hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
#dp <- dp[,1:3]

dt <- hl.out[,1:11];dt <- dt[,-c(7)]
hl.out_truth<-hl.out
hl.out<-0


Calculate_FS <- function(x){
  A0. <- dt$a[c]
  Ainf. <- dt$a[c]*(1-dt$pss[c])^dt$N[c]
  FS. <- (x-A0.)/(Ainf.-A0.)
  return(FS.)
}



Refitting_Function <- function(ki){
  current_k <<- ki                # May have to set the current k from within the function as a global parameter for use in the Core function.
  Refitting_Predicted <- sapply(ds$t,Refitting_Function_Core)
  Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
  return(1-Refitting_R2)
}

## This is the core function for reoptimization with KL - given a k and and t, plot out the predicted A0
Refitting_Function_Core <-function(x){
  z <- 0
  N <- dt$N[c]
  kp <- dt$kp[c]
  pss <- dt$pss[c]
  a <- dt$a[c]
  for (n in 0:N) {
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (current_k/(current_k-n*kp))*b
    #y<- a*(bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))                                                                  
    y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
    z <- y + z}
  final <- (z-1)/((1-pss)^N-1)
  return(final)
}


# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Optimize$par*x)*(1/(dt$N[c]+1)-bp))                                                                  
    z <- y + z}
  return(z)
}




Model <- function(t){
  z <- 0
  for (n in 0:dt$N[c]) {
    k <- dt$k[c]
    N <- dt$N[c]
    kp <- dt$kp[c]
    pss <- dt$pss[c]
    a <- dt$a[c]
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (k/(k-n*kp))*b
    y<- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))                                                                  
    z <- y + z}
  return(z)
}


Calculate_dk <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*((n*dt$kp[c])/(Optimize$par*(Optimize$par-n*dt$kp[c]))*bp*(exp(-Optimize$par*x)-exp(-n*dt$kp[c]*x))-x*(1/(dt$N[c]+1)-bp)*exp(-Optimize$par*x))                                                              
    z <- y + z}
  return(se_new/z)}




if (Refit == TRUE){
  for (c in 1:nrow(dt)) {
    
    print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
    id <- dt$ID[c] 
    ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
    
    #dkTable <- sapply(ds$t,Calculate_dk)
    #dk <- min(abs(dkTable))
    if (nrow(ds)!=0){ ################################################
      #if (nrow(ds)>1){
      ### Calculate Fractional Synthesis
      FS <- sapply(ds$A0, Calculate_FS)
      ds <- cbind(ds,FS)
      
      #Upper <- dt$k[c]+dk
      #Lower <- dt$k[c]^2/(dt$k[c]+dk)
      #Model_Predicted <- sapply(ds$t,Model)
      #Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
      
      # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
      #Optimize <- optim(0.3,Refitting_Function,method=optim_method,lower=0, upper = 1) #optim(start value, fxn) # Use optim() for Nelder-Mead
      Optimize <- optim(0.3,Refitting_Function)
      #if(Optimize$par > 10){   ##############################3
      #Optimize$par <- 10}    #########################################
      
      # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
      Refitted_Predicted <- sapply(ds$t,Refitted_Model)
      Refitted_R2 <- 1- (sum((ds$A0-Refitted_Predicted)^2))/(sum((ds$A0-mean(ds$A0))^2))
      
      #print(dt$k[c])
      print(Optimize$par)
      print(Refitted_R2)
      #ss_new<-sum((ds$A0-mean(ds$A0))^2)
      #se_new<-sum((ds$A0-Refitted_Predicted)^2)
      #se_new<-(sum((ds$A0-Refitted_Predicted)^2)/(dt$DP[c]-1))^0.5  ##################################33
      se_new<-(sum((ds$A0-Refitted_Predicted)^2)/9)^0.5  #6 for seven tps....5 for 6 tps  (female n male)
      print(se_new)
      
      
      
      dkTable <- sapply(ds$t,Calculate_dk)
      dk <- min(abs(dkTable))
      
      #oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")
      oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c],dt$DP[c],dt$z[c],"0",round(se_new,4),round(dt$a[c],3), round(dt$pss[c],4),dt$kp[c], dt$N[c],round(Optimize$par,3), round(dk,3),  round(Refitted_R2,3),      sep = "\t")
      #oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
      write(oput, file=output_file, append=T)
    }
  }
  
}







rm(list = ls(all = T)) 

#R2_threshold <- 0.8
#SE_threshold <- 0.05 

#home_directory <- "~/Documents"
hl.out_location <- "hl.out"
#hl.data.out_location <- "hl-data_verified.out" 
hl.data.out_location <- "hs3_imp.2.out"   ############################3

output_file <- "ProTurn_Output_hs3_imp2.txt"  ######################################
#output_file <- "TestSS_SE.txt"

Refit <- TRUE           # Should the Grapher attempt to gather the data points for each peptide and refit the kinetic curve in R.
#optim_method = "Brent"  # Optimization method used, refer to ?optim for details
optim_method = "Nelder-Mead"


#oput <- paste("ID","Uniprot","Peptide","z","R2","SS","SE","k","dk","DP","a","pss","kp","N", sep = "\t")
oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
write(oput, file=output_file, append=T)

hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
#dp <- dp[,1:3]

dt <- hl.out[,1:11];dt <- dt[,-c(7)]
hl.out_truth<-hl.out
hl.out<-0


Calculate_FS <- function(x){
  A0. <- dt$a[c]
  Ainf. <- dt$a[c]*(1-dt$pss[c])^dt$N[c]
  FS. <- (x-A0.)/(Ainf.-A0.)
  return(FS.)
}



Refitting_Function <- function(ki){
  current_k <<- ki                # May have to set the current k from within the function as a global parameter for use in the Core function.
  Refitting_Predicted <- sapply(ds$t,Refitting_Function_Core)
  Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
  return(1-Refitting_R2)
}

## This is the core function for reoptimization with KL - given a k and and t, plot out the predicted A0
Refitting_Function_Core <-function(x){
  z <- 0
  N <- dt$N[c]
  kp <- dt$kp[c]
  pss <- dt$pss[c]
  a <- dt$a[c]
  for (n in 0:N) {
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (current_k/(current_k-n*kp))*b
    #y<- a*(bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))                                                                  
    y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
    z <- y + z}
  final <- (z-1)/((1-pss)^N-1)
  return(final)
}


# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Optimize$par*x)*(1/(dt$N[c]+1)-bp))                                                                  
    z <- y + z}
  return(z)
}




Model <- function(t){
  z <- 0
  for (n in 0:dt$N[c]) {
    k <- dt$k[c]
    N <- dt$N[c]
    kp <- dt$kp[c]
    pss <- dt$pss[c]
    a <- dt$a[c]
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (k/(k-n*kp))*b
    y<- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))                                                                  
    z <- y + z}
  return(z)
}


Calculate_dk <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*((n*dt$kp[c])/(Optimize$par*(Optimize$par-n*dt$kp[c]))*bp*(exp(-Optimize$par*x)-exp(-n*dt$kp[c]*x))-x*(1/(dt$N[c]+1)-bp)*exp(-Optimize$par*x))                                                              
    z <- y + z}
  return(se_new/z)}




if (Refit == TRUE){
  for (c in 1:nrow(dt)) {
    
    print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
    id <- dt$ID[c] 
    ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
    
    #dkTable <- sapply(ds$t,Calculate_dk)
    #dk <- min(abs(dkTable))
    if (nrow(ds)!=0){ ################################################
      #if (nrow(ds)>1){
      ### Calculate Fractional Synthesis
      FS <- sapply(ds$A0, Calculate_FS)
      ds <- cbind(ds,FS)
      
      #Upper <- dt$k[c]+dk
      #Lower <- dt$k[c]^2/(dt$k[c]+dk)
      #Model_Predicted <- sapply(ds$t,Model)
      #Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
      
      # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
      #Optimize <- optim(0.3,Refitting_Function,method=optim_method,lower=0, upper = 1) #optim(start value, fxn) # Use optim() for Nelder-Mead
      Optimize <- optim(0.3,Refitting_Function)
      #if(Optimize$par > 10){   ##############################3
      #Optimize$par <- 10}    #########################################
      
      # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
      Refitted_Predicted <- sapply(ds$t,Refitted_Model)
      Refitted_R2 <- 1- (sum((ds$A0-Refitted_Predicted)^2))/(sum((ds$A0-mean(ds$A0))^2))
      
      #print(dt$k[c])
      print(Optimize$par)
      print(Refitted_R2)
      #ss_new<-sum((ds$A0-mean(ds$A0))^2)
      #se_new<-sum((ds$A0-Refitted_Predicted)^2)
      #se_new<-(sum((ds$A0-Refitted_Predicted)^2)/(dt$DP[c]-1))^0.5  ##################################33
      se_new<-(sum((ds$A0-Refitted_Predicted)^2)/9)^0.5  #6 for seven tps....5 for 6 tps  (female n male)
      print(se_new)
      
      
      
      dkTable <- sapply(ds$t,Calculate_dk)
      dk <- min(abs(dkTable))
      
      #oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")
      oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c],dt$DP[c],dt$z[c],"0",round(se_new,4),round(dt$a[c],3), round(dt$pss[c],4),dt$kp[c], dt$N[c],round(Optimize$par,3), round(dk,3),  round(Refitted_R2,3),      sep = "\t")
      #oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
      write(oput, file=output_file, append=T)
    }
  }
  
}








rm(list = ls(all = T)) 

#R2_threshold <- 0.8
#SE_threshold <- 0.05 

#home_directory <- "~/Documents"
hl.out_location <- "hl.out"
#hl.data.out_location <- "hl-data_verified.out" 
hl.data.out_location <- "hs3_imp.3.out"   ############################3

output_file <- "ProTurn_Output_hs3_imp3.txt"  ######################################
#output_file <- "TestSS_SE.txt"

Refit <- TRUE           # Should the Grapher attempt to gather the data points for each peptide and refit the kinetic curve in R.
#optim_method = "Brent"  # Optimization method used, refer to ?optim for details
optim_method = "Nelder-Mead"


#oput <- paste("ID","Uniprot","Peptide","z","R2","SS","SE","k","dk","DP","a","pss","kp","N", sep = "\t")
oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
write(oput, file=output_file, append=T)

hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
#dp <- dp[,1:3]

dt <- hl.out[,1:11];dt <- dt[,-c(7)]
hl.out_truth<-hl.out
hl.out<-0


Calculate_FS <- function(x){
  A0. <- dt$a[c]
  Ainf. <- dt$a[c]*(1-dt$pss[c])^dt$N[c]
  FS. <- (x-A0.)/(Ainf.-A0.)
  return(FS.)
}



Refitting_Function <- function(ki){
  current_k <<- ki                # May have to set the current k from within the function as a global parameter for use in the Core function.
  Refitting_Predicted <- sapply(ds$t,Refitting_Function_Core)
  Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
  return(1-Refitting_R2)
}

## This is the core function for reoptimization with KL - given a k and and t, plot out the predicted A0
Refitting_Function_Core <-function(x){
  z <- 0
  N <- dt$N[c]
  kp <- dt$kp[c]
  pss <- dt$pss[c]
  a <- dt$a[c]
  for (n in 0:N) {
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (current_k/(current_k-n*kp))*b
    #y<- a*(bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))                                                                  
    y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
    z <- y + z}
  final <- (z-1)/((1-pss)^N-1)
  return(final)
}


# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Optimize$par*x)*(1/(dt$N[c]+1)-bp))                                                                  
    z <- y + z}
  return(z)
}




Model <- function(t){
  z <- 0
  for (n in 0:dt$N[c]) {
    k <- dt$k[c]
    N <- dt$N[c]
    kp <- dt$kp[c]
    pss <- dt$pss[c]
    a <- dt$a[c]
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (k/(k-n*kp))*b
    y<- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))                                                                  
    z <- y + z}
  return(z)
}


Calculate_dk <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*((n*dt$kp[c])/(Optimize$par*(Optimize$par-n*dt$kp[c]))*bp*(exp(-Optimize$par*x)-exp(-n*dt$kp[c]*x))-x*(1/(dt$N[c]+1)-bp)*exp(-Optimize$par*x))                                                              
    z <- y + z}
  return(se_new/z)}




if (Refit == TRUE){
  for (c in 1:nrow(dt)) {
    
    print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
    id <- dt$ID[c] 
    ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
    
    #dkTable <- sapply(ds$t,Calculate_dk)
    #dk <- min(abs(dkTable))
    if (nrow(ds)!=0){ ################################################
      #if (nrow(ds)>1){
      ### Calculate Fractional Synthesis
      FS <- sapply(ds$A0, Calculate_FS)
      ds <- cbind(ds,FS)
      
      #Upper <- dt$k[c]+dk
      #Lower <- dt$k[c]^2/(dt$k[c]+dk)
      #Model_Predicted <- sapply(ds$t,Model)
      #Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
      
      # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
      #Optimize <- optim(0.3,Refitting_Function,method=optim_method,lower=0, upper = 1) #optim(start value, fxn) # Use optim() for Nelder-Mead
      Optimize <- optim(0.3,Refitting_Function)
      #if(Optimize$par > 10){   ##############################3
      #Optimize$par <- 10}    #########################################
      
      # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
      Refitted_Predicted <- sapply(ds$t,Refitted_Model)
      Refitted_R2 <- 1- (sum((ds$A0-Refitted_Predicted)^2))/(sum((ds$A0-mean(ds$A0))^2))
      
      #print(dt$k[c])
      print(Optimize$par)
      print(Refitted_R2)
      #ss_new<-sum((ds$A0-mean(ds$A0))^2)
      #se_new<-sum((ds$A0-Refitted_Predicted)^2)
      #se_new<-(sum((ds$A0-Refitted_Predicted)^2)/(dt$DP[c]-1))^0.5  ##################################33
      se_new<-(sum((ds$A0-Refitted_Predicted)^2)/9)^0.5  #6 for seven tps....5 for 6 tps  (female n male)
      print(se_new)
      
      
      
      dkTable <- sapply(ds$t,Calculate_dk)
      dk <- min(abs(dkTable))
      
      #oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")
      oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c],dt$DP[c],dt$z[c],"0",round(se_new,4),round(dt$a[c],3), round(dt$pss[c],4),dt$kp[c], dt$N[c],round(Optimize$par,3), round(dk,3),  round(Refitted_R2,3),      sep = "\t")
      #oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
      write(oput, file=output_file, append=T)
    }
  }
  
}





rm(list = ls(all = T)) 

#R2_threshold <- 0.8
#SE_threshold <- 0.05 

#home_directory <- "~/Documents"
hl.out_location <- "hl.out"
#hl.data.out_location <- "hl-data_verified.out" 
hl.data.out_location <- "hs3_imp.4.out"   ############################3

output_file <- "ProTurn_Output_hs3_imp4.txt"  ######################################
#output_file <- "TestSS_SE.txt"

Refit <- TRUE           # Should the Grapher attempt to gather the data points for each peptide and refit the kinetic curve in R.
#optim_method = "Brent"  # Optimization method used, refer to ?optim for details
optim_method = "Nelder-Mead"


#oput <- paste("ID","Uniprot","Peptide","z","R2","SS","SE","k","dk","DP","a","pss","kp","N", sep = "\t")
oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
write(oput, file=output_file, append=T)

hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
#dp <- dp[,1:3]

dt <- hl.out[,1:11];dt <- dt[,-c(7)]
hl.out_truth<-hl.out
hl.out<-0


Calculate_FS <- function(x){
  A0. <- dt$a[c]
  Ainf. <- dt$a[c]*(1-dt$pss[c])^dt$N[c]
  FS. <- (x-A0.)/(Ainf.-A0.)
  return(FS.)
}



Refitting_Function <- function(ki){
  current_k <<- ki                # May have to set the current k from within the function as a global parameter for use in the Core function.
  Refitting_Predicted <- sapply(ds$t,Refitting_Function_Core)
  Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
  return(1-Refitting_R2)
}

## This is the core function for reoptimization with KL - given a k and and t, plot out the predicted A0
Refitting_Function_Core <-function(x){
  z <- 0
  N <- dt$N[c]
  kp <- dt$kp[c]
  pss <- dt$pss[c]
  a <- dt$a[c]
  for (n in 0:N) {
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (current_k/(current_k-n*kp))*b
    #y<- a*(bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))                                                                  
    y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
    z <- y + z}
  final <- (z-1)/((1-pss)^N-1)
  return(final)
}


# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Optimize$par*x)*(1/(dt$N[c]+1)-bp))                                                                  
    z <- y + z}
  return(z)
}




Model <- function(t){
  z <- 0
  for (n in 0:dt$N[c]) {
    k <- dt$k[c]
    N <- dt$N[c]
    kp <- dt$kp[c]
    pss <- dt$pss[c]
    a <- dt$a[c]
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (k/(k-n*kp))*b
    y<- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))                                                                  
    z <- y + z}
  return(z)
}


Calculate_dk <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*((n*dt$kp[c])/(Optimize$par*(Optimize$par-n*dt$kp[c]))*bp*(exp(-Optimize$par*x)-exp(-n*dt$kp[c]*x))-x*(1/(dt$N[c]+1)-bp)*exp(-Optimize$par*x))                                                              
    z <- y + z}
  return(se_new/z)}




if (Refit == TRUE){
  for (c in 1:nrow(dt)) {
    
    print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
    id <- dt$ID[c] 
    ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
    
    #dkTable <- sapply(ds$t,Calculate_dk)
    #dk <- min(abs(dkTable))
    if (nrow(ds)!=0){ ################################################
      #if (nrow(ds)>1){
      ### Calculate Fractional Synthesis
      FS <- sapply(ds$A0, Calculate_FS)
      ds <- cbind(ds,FS)
      
      #Upper <- dt$k[c]+dk
      #Lower <- dt$k[c]^2/(dt$k[c]+dk)
      #Model_Predicted <- sapply(ds$t,Model)
      #Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
      
      # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
      #Optimize <- optim(0.3,Refitting_Function,method=optim_method,lower=0, upper = 1) #optim(start value, fxn) # Use optim() for Nelder-Mead
      Optimize <- optim(0.3,Refitting_Function)
      #if(Optimize$par > 10){   ##############################3
      #Optimize$par <- 10}    #########################################
      
      # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
      Refitted_Predicted <- sapply(ds$t,Refitted_Model)
      Refitted_R2 <- 1- (sum((ds$A0-Refitted_Predicted)^2))/(sum((ds$A0-mean(ds$A0))^2))
      
      #print(dt$k[c])
      print(Optimize$par)
      print(Refitted_R2)
      #ss_new<-sum((ds$A0-mean(ds$A0))^2)
      #se_new<-sum((ds$A0-Refitted_Predicted)^2)
      #se_new<-(sum((ds$A0-Refitted_Predicted)^2)/(dt$DP[c]-1))^0.5  ##################################33
      se_new<-(sum((ds$A0-Refitted_Predicted)^2)/9)^0.5  #6 for seven tps....5 for 6 tps  (female n male)
      print(se_new)
      
      
      
      dkTable <- sapply(ds$t,Calculate_dk)
      dk <- min(abs(dkTable))
      
      #oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")
      oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c],dt$DP[c],dt$z[c],"0",round(se_new,4),round(dt$a[c],3), round(dt$pss[c],4),dt$kp[c], dt$N[c],round(Optimize$par,3), round(dk,3),  round(Refitted_R2,3),      sep = "\t")
      #oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
      write(oput, file=output_file, append=T)
    }
  }
  
}




rm(list = ls(all = T)) 

#R2_threshold <- 0.8
#SE_threshold <- 0.05 

#home_directory <- "~/Documents"
hl.out_location <- "hl.out"
#hl.data.out_location <- "hl-data_verified.out" 
hl.data.out_location <- "hs3_imp.5.out"   ############################3

output_file <- "ProTurn_Output_hs3_imp5.txt"  ######################################
#output_file <- "TestSS_SE.txt"

Refit <- TRUE           # Should the Grapher attempt to gather the data points for each peptide and refit the kinetic curve in R.
#optim_method = "Brent"  # Optimization method used, refer to ?optim for details
optim_method = "Nelder-Mead"


#oput <- paste("ID","Uniprot","Peptide","z","R2","SS","SE","k","dk","DP","a","pss","kp","N", sep = "\t")
oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
write(oput, file=output_file, append=T)

hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
#dp <- dp[,1:3]

dt <- hl.out[,1:11];dt <- dt[,-c(7)]
hl.out_truth<-hl.out
hl.out<-0


Calculate_FS <- function(x){
  A0. <- dt$a[c]
  Ainf. <- dt$a[c]*(1-dt$pss[c])^dt$N[c]
  FS. <- (x-A0.)/(Ainf.-A0.)
  return(FS.)
}



Refitting_Function <- function(ki){
  current_k <<- ki                # May have to set the current k from within the function as a global parameter for use in the Core function.
  Refitting_Predicted <- sapply(ds$t,Refitting_Function_Core)
  Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
  return(1-Refitting_R2)
}

## This is the core function for reoptimization with KL - given a k and and t, plot out the predicted A0
Refitting_Function_Core <-function(x){
  z <- 0
  N <- dt$N[c]
  kp <- dt$kp[c]
  pss <- dt$pss[c]
  a <- dt$a[c]
  for (n in 0:N) {
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (current_k/(current_k-n*kp))*b
    #y<- a*(bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))                                                                  
    y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
    z <- y + z}
  final <- (z-1)/((1-pss)^N-1)
  return(final)
}


# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Optimize$par*x)*(1/(dt$N[c]+1)-bp))                                                                  
    z <- y + z}
  return(z)
}




Model <- function(t){
  z <- 0
  for (n in 0:dt$N[c]) {
    k <- dt$k[c]
    N <- dt$N[c]
    kp <- dt$kp[c]
    pss <- dt$pss[c]
    a <- dt$a[c]
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (k/(k-n*kp))*b
    y<- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))                                                                  
    z <- y + z}
  return(z)
}


Calculate_dk <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*((n*dt$kp[c])/(Optimize$par*(Optimize$par-n*dt$kp[c]))*bp*(exp(-Optimize$par*x)-exp(-n*dt$kp[c]*x))-x*(1/(dt$N[c]+1)-bp)*exp(-Optimize$par*x))                                                              
    z <- y + z}
  return(se_new/z)}




if (Refit == TRUE){
  for (c in 1:nrow(dt)) {
    
    print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
    id <- dt$ID[c] 
    ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
    
    #dkTable <- sapply(ds$t,Calculate_dk)
    #dk <- min(abs(dkTable))
    if (nrow(ds)!=0){ ################################################
      #if (nrow(ds)>1){
      ### Calculate Fractional Synthesis
      FS <- sapply(ds$A0, Calculate_FS)
      ds <- cbind(ds,FS)
      
      #Upper <- dt$k[c]+dk
      #Lower <- dt$k[c]^2/(dt$k[c]+dk)
      #Model_Predicted <- sapply(ds$t,Model)
      #Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
      
      # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
      #Optimize <- optim(0.3,Refitting_Function,method=optim_method,lower=0, upper = 1) #optim(start value, fxn) # Use optim() for Nelder-Mead
      Optimize <- optim(0.3,Refitting_Function)
      #if(Optimize$par > 10){   ##############################3
      #Optimize$par <- 10}    #########################################
      
      # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
      Refitted_Predicted <- sapply(ds$t,Refitted_Model)
      Refitted_R2 <- 1- (sum((ds$A0-Refitted_Predicted)^2))/(sum((ds$A0-mean(ds$A0))^2))
      
      #print(dt$k[c])
      print(Optimize$par)
      print(Refitted_R2)
      #ss_new<-sum((ds$A0-mean(ds$A0))^2)
      #se_new<-sum((ds$A0-Refitted_Predicted)^2)
      #se_new<-(sum((ds$A0-Refitted_Predicted)^2)/(dt$DP[c]-1))^0.5  ##################################33
      se_new<-(sum((ds$A0-Refitted_Predicted)^2)/9)^0.5  #6 for seven tps....5 for 6 tps  (female n male)
      print(se_new)
      
      
      
      dkTable <- sapply(ds$t,Calculate_dk)
      dk <- min(abs(dkTable))
      
      #oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")
      oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c],dt$DP[c],dt$z[c],"0",round(se_new,4),round(dt$a[c],3), round(dt$pss[c],4),dt$kp[c], dt$N[c],round(Optimize$par,3), round(dk,3),  round(Refitted_R2,3),      sep = "\t")
      #oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
      write(oput, file=output_file, append=T)
    }
  }
  
}



rm(list = ls(all = T)) 

#R2_threshold <- 0.8
#SE_threshold <- 0.05 

#home_directory <- "~/Documents"
hl.out_location <- "hl.out"
#hl.data.out_location <- "hl-data_verified.out" 
hl.data.out_location <- "hs3_imp.6.out"   ############################3

output_file <- "ProTurn_Output_hs3_imp6.txt"  ######################################
#output_file <- "TestSS_SE.txt"

Refit <- TRUE           # Should the Grapher attempt to gather the data points for each peptide and refit the kinetic curve in R.
#optim_method = "Brent"  # Optimization method used, refer to ?optim for details
optim_method = "Nelder-Mead"


#oput <- paste("ID","Uniprot","Peptide","z","R2","SS","SE","k","dk","DP","a","pss","kp","N", sep = "\t")
oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
write(oput, file=output_file, append=T)

hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
#dp <- dp[,1:3]

dt <- hl.out[,1:11];dt <- dt[,-c(7)]
hl.out_truth<-hl.out
hl.out<-0


Calculate_FS <- function(x){
  A0. <- dt$a[c]
  Ainf. <- dt$a[c]*(1-dt$pss[c])^dt$N[c]
  FS. <- (x-A0.)/(Ainf.-A0.)
  return(FS.)
}



Refitting_Function <- function(ki){
  current_k <<- ki                # May have to set the current k from within the function as a global parameter for use in the Core function.
  Refitting_Predicted <- sapply(ds$t,Refitting_Function_Core)
  Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
  return(1-Refitting_R2)
}

## This is the core function for reoptimization with KL - given a k and and t, plot out the predicted A0
Refitting_Function_Core <-function(x){
  z <- 0
  N <- dt$N[c]
  kp <- dt$kp[c]
  pss <- dt$pss[c]
  a <- dt$a[c]
  for (n in 0:N) {
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (current_k/(current_k-n*kp))*b
    #y<- a*(bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))                                                                  
    y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
    z <- y + z}
  final <- (z-1)/((1-pss)^N-1)
  return(final)
}


# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Optimize$par*x)*(1/(dt$N[c]+1)-bp))                                                                  
    z <- y + z}
  return(z)
}




Model <- function(t){
  z <- 0
  for (n in 0:dt$N[c]) {
    k <- dt$k[c]
    N <- dt$N[c]
    kp <- dt$kp[c]
    pss <- dt$pss[c]
    a <- dt$a[c]
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (k/(k-n*kp))*b
    y<- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))                                                                  
    z <- y + z}
  return(z)
}


Calculate_dk <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*((n*dt$kp[c])/(Optimize$par*(Optimize$par-n*dt$kp[c]))*bp*(exp(-Optimize$par*x)-exp(-n*dt$kp[c]*x))-x*(1/(dt$N[c]+1)-bp)*exp(-Optimize$par*x))                                                              
    z <- y + z}
  return(se_new/z)}




if (Refit == TRUE){
  for (c in 1:nrow(dt)) {
    
    print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
    id <- dt$ID[c] 
    ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
    
    #dkTable <- sapply(ds$t,Calculate_dk)
    #dk <- min(abs(dkTable))
    if (nrow(ds)!=0){ ################################################
      #if (nrow(ds)>1){
      ### Calculate Fractional Synthesis
      FS <- sapply(ds$A0, Calculate_FS)
      ds <- cbind(ds,FS)
      
      #Upper <- dt$k[c]+dk
      #Lower <- dt$k[c]^2/(dt$k[c]+dk)
      #Model_Predicted <- sapply(ds$t,Model)
      #Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
      
      # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
      #Optimize <- optim(0.3,Refitting_Function,method=optim_method,lower=0, upper = 1) #optim(start value, fxn) # Use optim() for Nelder-Mead
      Optimize <- optim(0.3,Refitting_Function)
      #if(Optimize$par > 10){   ##############################3
      #Optimize$par <- 10}    #########################################
      
      # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
      Refitted_Predicted <- sapply(ds$t,Refitted_Model)
      Refitted_R2 <- 1- (sum((ds$A0-Refitted_Predicted)^2))/(sum((ds$A0-mean(ds$A0))^2))
      
      #print(dt$k[c])
      print(Optimize$par)
      print(Refitted_R2)
      #ss_new<-sum((ds$A0-mean(ds$A0))^2)
      #se_new<-sum((ds$A0-Refitted_Predicted)^2)
      #se_new<-(sum((ds$A0-Refitted_Predicted)^2)/(dt$DP[c]-1))^0.5  ##################################33
      se_new<-(sum((ds$A0-Refitted_Predicted)^2)/9)^0.5  #6 for seven tps....5 for 6 tps  (female n male)
      print(se_new)
      
      
      
      dkTable <- sapply(ds$t,Calculate_dk)
      dk <- min(abs(dkTable))
      
      #oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")
      oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c],dt$DP[c],dt$z[c],"0",round(se_new,4),round(dt$a[c],3), round(dt$pss[c],4),dt$kp[c], dt$N[c],round(Optimize$par,3), round(dk,3),  round(Refitted_R2,3),      sep = "\t")
      #oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
      write(oput, file=output_file, append=T)
    }
  }
  
}



rm(list = ls(all = T)) 

#R2_threshold <- 0.8
#SE_threshold <- 0.05 

#home_directory <- "~/Documents"
hl.out_location <- "hl.out"
#hl.data.out_location <- "hl-data_verified.out" 
hl.data.out_location <- "hs3_imp.7.out"   ############################3

output_file <- "ProTurn_Output_hs3_imp7.txt"  ######################################
#output_file <- "TestSS_SE.txt"

Refit <- TRUE           # Should the Grapher attempt to gather the data points for each peptide and refit the kinetic curve in R.
#optim_method = "Brent"  # Optimization method used, refer to ?optim for details
optim_method = "Nelder-Mead"


#oput <- paste("ID","Uniprot","Peptide","z","R2","SS","SE","k","dk","DP","a","pss","kp","N", sep = "\t")
oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
write(oput, file=output_file, append=T)

hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
#dp <- dp[,1:3]

dt <- hl.out[,1:11];dt <- dt[,-c(7)]
hl.out_truth<-hl.out
hl.out<-0


Calculate_FS <- function(x){
  A0. <- dt$a[c]
  Ainf. <- dt$a[c]*(1-dt$pss[c])^dt$N[c]
  FS. <- (x-A0.)/(Ainf.-A0.)
  return(FS.)
}



Refitting_Function <- function(ki){
  current_k <<- ki                # May have to set the current k from within the function as a global parameter for use in the Core function.
  Refitting_Predicted <- sapply(ds$t,Refitting_Function_Core)
  Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
  return(1-Refitting_R2)
}

## This is the core function for reoptimization with KL - given a k and and t, plot out the predicted A0
Refitting_Function_Core <-function(x){
  z <- 0
  N <- dt$N[c]
  kp <- dt$kp[c]
  pss <- dt$pss[c]
  a <- dt$a[c]
  for (n in 0:N) {
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (current_k/(current_k-n*kp))*b
    #y<- a*(bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))                                                                  
    y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
    z <- y + z}
  final <- (z-1)/((1-pss)^N-1)
  return(final)
}


# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Optimize$par*x)*(1/(dt$N[c]+1)-bp))                                                                  
    z <- y + z}
  return(z)
}




Model <- function(t){
  z <- 0
  for (n in 0:dt$N[c]) {
    k <- dt$k[c]
    N <- dt$N[c]
    kp <- dt$kp[c]
    pss <- dt$pss[c]
    a <- dt$a[c]
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (k/(k-n*kp))*b
    y<- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))                                                                  
    z <- y + z}
  return(z)
}


Calculate_dk <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*((n*dt$kp[c])/(Optimize$par*(Optimize$par-n*dt$kp[c]))*bp*(exp(-Optimize$par*x)-exp(-n*dt$kp[c]*x))-x*(1/(dt$N[c]+1)-bp)*exp(-Optimize$par*x))                                                              
    z <- y + z}
  return(se_new/z)}




if (Refit == TRUE){
  for (c in 1:nrow(dt)) {
    
    print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
    id <- dt$ID[c] 
    ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
    
    #dkTable <- sapply(ds$t,Calculate_dk)
    #dk <- min(abs(dkTable))
    if (nrow(ds)!=0){ ################################################
      #if (nrow(ds)>1){
      ### Calculate Fractional Synthesis
      FS <- sapply(ds$A0, Calculate_FS)
      ds <- cbind(ds,FS)
      
      #Upper <- dt$k[c]+dk
      #Lower <- dt$k[c]^2/(dt$k[c]+dk)
      #Model_Predicted <- sapply(ds$t,Model)
      #Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
      
      # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
      #Optimize <- optim(0.3,Refitting_Function,method=optim_method,lower=0, upper = 1) #optim(start value, fxn) # Use optim() for Nelder-Mead
      Optimize <- optim(0.3,Refitting_Function)
      #if(Optimize$par > 10){   ##############################3
      #Optimize$par <- 10}    #########################################
      
      # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
      Refitted_Predicted <- sapply(ds$t,Refitted_Model)
      Refitted_R2 <- 1- (sum((ds$A0-Refitted_Predicted)^2))/(sum((ds$A0-mean(ds$A0))^2))
      
      #print(dt$k[c])
      print(Optimize$par)
      print(Refitted_R2)
      #ss_new<-sum((ds$A0-mean(ds$A0))^2)
      #se_new<-sum((ds$A0-Refitted_Predicted)^2)
      #se_new<-(sum((ds$A0-Refitted_Predicted)^2)/(dt$DP[c]-1))^0.5  ##################################33
      se_new<-(sum((ds$A0-Refitted_Predicted)^2)/9)^0.5  #6 for seven tps....5 for 6 tps  (female n male)
      print(se_new)
      
      
      
      dkTable <- sapply(ds$t,Calculate_dk)
      dk <- min(abs(dkTable))
      
      #oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")
      oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c],dt$DP[c],dt$z[c],"0",round(se_new,4),round(dt$a[c],3), round(dt$pss[c],4),dt$kp[c], dt$N[c],round(Optimize$par,3), round(dk,3),  round(Refitted_R2,3),      sep = "\t")
      #oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
      write(oput, file=output_file, append=T)
    }
  }
  
}


rm(list = ls(all = T)) 

#R2_threshold <- 0.8
#SE_threshold <- 0.05 

#home_directory <- "~/Documents"
hl.out_location <- "hl.out"
#hl.data.out_location <- "hl-data_verified.out" 
hl.data.out_location <- "hs3_imp.8.out"   ############################3

output_file <- "ProTurn_Output_hs3_imp8.txt"  ######################################
#output_file <- "TestSS_SE.txt"

Refit <- TRUE           # Should the Grapher attempt to gather the data points for each peptide and refit the kinetic curve in R.
#optim_method = "Brent"  # Optimization method used, refer to ?optim for details
optim_method = "Nelder-Mead"


#oput <- paste("ID","Uniprot","Peptide","z","R2","SS","SE","k","dk","DP","a","pss","kp","N", sep = "\t")
oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
write(oput, file=output_file, append=T)

hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
#dp <- dp[,1:3]

dt <- hl.out[,1:11];dt <- dt[,-c(7)]
hl.out_truth<-hl.out
hl.out<-0


Calculate_FS <- function(x){
  A0. <- dt$a[c]
  Ainf. <- dt$a[c]*(1-dt$pss[c])^dt$N[c]
  FS. <- (x-A0.)/(Ainf.-A0.)
  return(FS.)
}



Refitting_Function <- function(ki){
  current_k <<- ki                # May have to set the current k from within the function as a global parameter for use in the Core function.
  Refitting_Predicted <- sapply(ds$t,Refitting_Function_Core)
  Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
  return(1-Refitting_R2)
}

## This is the core function for reoptimization with KL - given a k and and t, plot out the predicted A0
Refitting_Function_Core <-function(x){
  z <- 0
  N <- dt$N[c]
  kp <- dt$kp[c]
  pss <- dt$pss[c]
  a <- dt$a[c]
  for (n in 0:N) {
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (current_k/(current_k-n*kp))*b
    #y<- a*(bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))                                                                  
    y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
    z <- y + z}
  final <- (z-1)/((1-pss)^N-1)
  return(final)
}


# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Optimize$par*x)*(1/(dt$N[c]+1)-bp))                                                                  
    z <- y + z}
  return(z)
}




Model <- function(t){
  z <- 0
  for (n in 0:dt$N[c]) {
    k <- dt$k[c]
    N <- dt$N[c]
    kp <- dt$kp[c]
    pss <- dt$pss[c]
    a <- dt$a[c]
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (k/(k-n*kp))*b
    y<- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))                                                                  
    z <- y + z}
  return(z)
}


Calculate_dk <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*((n*dt$kp[c])/(Optimize$par*(Optimize$par-n*dt$kp[c]))*bp*(exp(-Optimize$par*x)-exp(-n*dt$kp[c]*x))-x*(1/(dt$N[c]+1)-bp)*exp(-Optimize$par*x))                                                              
    z <- y + z}
  return(se_new/z)}




if (Refit == TRUE){
  for (c in 1:nrow(dt)) {
    
    print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
    id <- dt$ID[c] 
    ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
    
    #dkTable <- sapply(ds$t,Calculate_dk)
    #dk <- min(abs(dkTable))
    if (nrow(ds)!=0){ ################################################
      #if (nrow(ds)>1){
      ### Calculate Fractional Synthesis
      FS <- sapply(ds$A0, Calculate_FS)
      ds <- cbind(ds,FS)
      
      #Upper <- dt$k[c]+dk
      #Lower <- dt$k[c]^2/(dt$k[c]+dk)
      #Model_Predicted <- sapply(ds$t,Model)
      #Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
      
      # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
      #Optimize <- optim(0.3,Refitting_Function,method=optim_method,lower=0, upper = 1) #optim(start value, fxn) # Use optim() for Nelder-Mead
      Optimize <- optim(0.3,Refitting_Function)
      #if(Optimize$par > 10){   ##############################3
      #Optimize$par <- 10}    #########################################
      
      # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
      Refitted_Predicted <- sapply(ds$t,Refitted_Model)
      Refitted_R2 <- 1- (sum((ds$A0-Refitted_Predicted)^2))/(sum((ds$A0-mean(ds$A0))^2))
      
      #print(dt$k[c])
      print(Optimize$par)
      print(Refitted_R2)
      #ss_new<-sum((ds$A0-mean(ds$A0))^2)
      #se_new<-sum((ds$A0-Refitted_Predicted)^2)
      #se_new<-(sum((ds$A0-Refitted_Predicted)^2)/(dt$DP[c]-1))^0.5  ##################################33
      se_new<-(sum((ds$A0-Refitted_Predicted)^2)/9)^0.5  #6 for seven tps....5 for 6 tps  (female n male)
      print(se_new)
      
      
      
      dkTable <- sapply(ds$t,Calculate_dk)
      dk <- min(abs(dkTable))
      
      #oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")
      oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c],dt$DP[c],dt$z[c],"0",round(se_new,4),round(dt$a[c],3), round(dt$pss[c],4),dt$kp[c], dt$N[c],round(Optimize$par,3), round(dk,3),  round(Refitted_R2,3),      sep = "\t")
      #oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
      write(oput, file=output_file, append=T)
    }
  }
  
}


rm(list = ls(all = T)) 

#R2_threshold <- 0.8
#SE_threshold <- 0.05 

#home_directory <- "~/Documents"
hl.out_location <- "hl.out"
#hl.data.out_location <- "hl-data_verified.out" 
hl.data.out_location <- "hs3_imp.9.out"   ############################3

output_file <- "ProTurn_Output_hs3_imp9.txt"  ######################################
#output_file <- "TestSS_SE.txt"

Refit <- TRUE           # Should the Grapher attempt to gather the data points for each peptide and refit the kinetic curve in R.
#optim_method = "Brent"  # Optimization method used, refer to ?optim for details
optim_method = "Nelder-Mead"


#oput <- paste("ID","Uniprot","Peptide","z","R2","SS","SE","k","dk","DP","a","pss","kp","N", sep = "\t")
oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
write(oput, file=output_file, append=T)

hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
#dp <- dp[,1:3]

dt <- hl.out[,1:11];dt <- dt[,-c(7)]
hl.out_truth<-hl.out
hl.out<-0


Calculate_FS <- function(x){
  A0. <- dt$a[c]
  Ainf. <- dt$a[c]*(1-dt$pss[c])^dt$N[c]
  FS. <- (x-A0.)/(Ainf.-A0.)
  return(FS.)
}



Refitting_Function <- function(ki){
  current_k <<- ki                # May have to set the current k from within the function as a global parameter for use in the Core function.
  Refitting_Predicted <- sapply(ds$t,Refitting_Function_Core)
  Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
  return(1-Refitting_R2)
}

## This is the core function for reoptimization with KL - given a k and and t, plot out the predicted A0
Refitting_Function_Core <-function(x){
  z <- 0
  N <- dt$N[c]
  kp <- dt$kp[c]
  pss <- dt$pss[c]
  a <- dt$a[c]
  for (n in 0:N) {
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (current_k/(current_k-n*kp))*b
    #y<- a*(bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))                                                                  
    y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
    z <- y + z}
  final <- (z-1)/((1-pss)^N-1)
  return(final)
}


# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Optimize$par*x)*(1/(dt$N[c]+1)-bp))                                                                  
    z <- y + z}
  return(z)
}




Model <- function(t){
  z <- 0
  for (n in 0:dt$N[c]) {
    k <- dt$k[c]
    N <- dt$N[c]
    kp <- dt$kp[c]
    pss <- dt$pss[c]
    a <- dt$a[c]
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (k/(k-n*kp))*b
    y<- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))                                                                  
    z <- y + z}
  return(z)
}


Calculate_dk <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*((n*dt$kp[c])/(Optimize$par*(Optimize$par-n*dt$kp[c]))*bp*(exp(-Optimize$par*x)-exp(-n*dt$kp[c]*x))-x*(1/(dt$N[c]+1)-bp)*exp(-Optimize$par*x))                                                              
    z <- y + z}
  return(se_new/z)}




if (Refit == TRUE){
  for (c in 1:nrow(dt)) {
    
    print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
    id <- dt$ID[c] 
    ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
    
    #dkTable <- sapply(ds$t,Calculate_dk)
    #dk <- min(abs(dkTable))
    if (nrow(ds)!=0){ ################################################
      #if (nrow(ds)>1){
      ### Calculate Fractional Synthesis
      FS <- sapply(ds$A0, Calculate_FS)
      ds <- cbind(ds,FS)
      
      #Upper <- dt$k[c]+dk
      #Lower <- dt$k[c]^2/(dt$k[c]+dk)
      #Model_Predicted <- sapply(ds$t,Model)
      #Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
      
      # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
      #Optimize <- optim(0.3,Refitting_Function,method=optim_method,lower=0, upper = 1) #optim(start value, fxn) # Use optim() for Nelder-Mead
      Optimize <- optim(0.3,Refitting_Function)
      #if(Optimize$par > 10){   ##############################3
      #Optimize$par <- 10}    #########################################
      
      # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
      Refitted_Predicted <- sapply(ds$t,Refitted_Model)
      Refitted_R2 <- 1- (sum((ds$A0-Refitted_Predicted)^2))/(sum((ds$A0-mean(ds$A0))^2))
      
      #print(dt$k[c])
      print(Optimize$par)
      print(Refitted_R2)
      #ss_new<-sum((ds$A0-mean(ds$A0))^2)
      #se_new<-sum((ds$A0-Refitted_Predicted)^2)
      #se_new<-(sum((ds$A0-Refitted_Predicted)^2)/(dt$DP[c]-1))^0.5  ##################################33
      se_new<-(sum((ds$A0-Refitted_Predicted)^2)/9)^0.5  #6 for seven tps....5 for 6 tps  (female n male)
      print(se_new)
      
      
      
      dkTable <- sapply(ds$t,Calculate_dk)
      dk <- min(abs(dkTable))
      
      #oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")
      oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c],dt$DP[c],dt$z[c],"0",round(se_new,4),round(dt$a[c],3), round(dt$pss[c],4),dt$kp[c], dt$N[c],round(Optimize$par,3), round(dk,3),  round(Refitted_R2,3),      sep = "\t")
      #oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
      write(oput, file=output_file, append=T)
    }
  }
  
}


rm(list = ls(all = T)) 

#R2_threshold <- 0.8
#SE_threshold <- 0.05 

#home_directory <- "~/Documents"
hl.out_location <- "hl.out"
#hl.data.out_location <- "hl-data_verified.out" 
hl.data.out_location <- "hs3_imp.10.out"   ############################3

output_file <- "ProTurn_Output_hs3_imp10.txt"  ######################################
#output_file <- "TestSS_SE.txt"

Refit <- TRUE           # Should the Grapher attempt to gather the data points for each peptide and refit the kinetic curve in R.
#optim_method = "Brent"  # Optimization method used, refer to ?optim for details
optim_method = "Nelder-Mead"


#oput <- paste("ID","Uniprot","Peptide","z","R2","SS","SE","k","dk","DP","a","pss","kp","N", sep = "\t")
oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
write(oput, file=output_file, append=T)

hl.out <- read.table(hl.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
dp <- read.table(hl.data.out_location, header=TRUE, fill=TRUE, skip = 0) #skip 1
#dp <- dp[,1:3]

dt <- hl.out[,1:11];dt <- dt[,-c(7)]
hl.out_truth<-hl.out
hl.out<-0


Calculate_FS <- function(x){
  A0. <- dt$a[c]
  Ainf. <- dt$a[c]*(1-dt$pss[c])^dt$N[c]
  FS. <- (x-A0.)/(Ainf.-A0.)
  return(FS.)
}



Refitting_Function <- function(ki){
  current_k <<- ki                # May have to set the current k from within the function as a global parameter for use in the Core function.
  Refitting_Predicted <- sapply(ds$t,Refitting_Function_Core)
  Refitting_R2 <- 1- (sum((ds$FS-Refitting_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
  return(1-Refitting_R2)
}

## This is the core function for reoptimization with KL - given a k and and t, plot out the predicted A0
Refitting_Function_Core <-function(x){
  z <- 0
  N <- dt$N[c]
  kp <- dt$kp[c]
  pss <- dt$pss[c]
  a <- dt$a[c]
  for (n in 0:N) {
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (current_k/(current_k-n*kp))*b
    #y<- a*(bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))                                                                  
    y<- (bp*exp(-n*kp*x)+exp(-current_k*x)*(1/(N+1)-bp))
    z <- y + z}
  final <- (z-1)/((1-pss)^N-1)
  return(final)
}


# This function takes in time information and returns the KL equation y value using the now optimized k
Refitted_Model <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*(bp*exp(-n*dt$kp[c]*x)+exp(-Optimize$par*x)*(1/(dt$N[c]+1)-bp))                                                                  
    z <- y + z}
  return(z)
}




Model <- function(t){
  z <- 0
  for (n in 0:dt$N[c]) {
    k <- dt$k[c]
    N <- dt$N[c]
    kp <- dt$kp[c]
    pss <- dt$pss[c]
    a <- dt$a[c]
    
    b <- factorial(N)/(factorial(n)*factorial(N-n))*(1-pss)^(N-n)*(pss)^n
    bp <- (k/(k-n*kp))*b
    y<- a*(bp*exp(-n*kp*t)+exp(-k*t)*(1/(N+1)-bp))                                                                  
    z <- y + z}
  return(z)
}


Calculate_dk <-function(x){
  z <- 0
  for (n in 0:dt$N[c]) {
    b <- factorial(dt$N[c])/(factorial(n)*factorial(dt$N[c]-n))*(1-dt$pss[c])^(dt$N[c]-n)*(dt$pss[c])^n
    bp <- (Optimize$par/(Optimize$par-n*dt$kp[c]))*b
    y<- dt$a[c]*((n*dt$kp[c])/(Optimize$par*(Optimize$par-n*dt$kp[c]))*bp*(exp(-Optimize$par*x)-exp(-n*dt$kp[c]*x))-x*(1/(dt$N[c]+1)-bp)*exp(-Optimize$par*x))                                                              
    z <- y + z}
  return(se_new/z)}




if (Refit == TRUE){
  for (c in 1:nrow(dt)) {
    
    print(paste("Now refitting peptide ", c, " of ", nrow(dt), ". ", round(c/nrow(dt)*100,2), "% done."))
    id <- dt$ID[c] 
    ds <- dp[ which(dp$ID == id),]  # Subsetting that particular ID
    
    #dkTable <- sapply(ds$t,Calculate_dk)
    #dk <- min(abs(dkTable))
    if (nrow(ds)!=0){ ################################################
      #if (nrow(ds)>1){
      ### Calculate Fractional Synthesis
      FS <- sapply(ds$A0, Calculate_FS)
      ds <- cbind(ds,FS)
      
      #Upper <- dt$k[c]+dk
      #Lower <- dt$k[c]^2/(dt$k[c]+dk)
      #Model_Predicted <- sapply(ds$t,Model)
      #Model_R2 <- 1- (sum((ds$FS-Model_Predicted)^2))/(sum((ds$FS-mean(ds$FS))^2))
      
      # This is the optimization function that tries out different k values and plug into the Refitting function. Optimzation$par will give you the now optimized K. 
      #Optimize <- optim(0.3,Refitting_Function,method=optim_method,lower=0, upper = 1) #optim(start value, fxn) # Use optim() for Nelder-Mead
      Optimize <- optim(0.3,Refitting_Function)
      #if(Optimize$par > 10){   ##############################3
      #Optimize$par <- 10}    #########################################
      
      # This is the KL function that takes in different values of k and return 1 minus R2 value, which the optim function will try to minimize
      Refitted_Predicted <- sapply(ds$t,Refitted_Model)
      Refitted_R2 <- 1- (sum((ds$A0-Refitted_Predicted)^2))/(sum((ds$A0-mean(ds$A0))^2))
      
      #print(dt$k[c])
      print(Optimize$par)
      print(Refitted_R2)
      #ss_new<-sum((ds$A0-mean(ds$A0))^2)
      #se_new<-sum((ds$A0-Refitted_Predicted)^2)
      #se_new<-(sum((ds$A0-Refitted_Predicted)^2)/(dt$DP[c]-1))^0.5  ##################################33
      se_new<-(sum((ds$A0-Refitted_Predicted)^2)/9)^0.5  #6 for seven tps....5 for 6 tps  (female n male)
      print(se_new)
      
      
      
      dkTable <- sapply(ds$t,Calculate_dk)
      dk <- min(abs(dkTable))
      
      #oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c], dt$z[c], round(dt$R2[c],3), round(dt$SS[c],4), round(dt$SE[c],4), round(log2(dt$k[c]),3), round(log2(dt$k[c]+dk)-log2(dt$k[c]),3), dt$DP[c], round(dt$a[c],3), round(dt$pss[c],4), dt$kp[c], dt$N[c], sep = "\t")
      oput <- paste(dt$ID[c],dt$Uniprot[c], dt$Peptide[c],dt$DP[c],dt$z[c],"0",round(se_new,4),round(dt$a[c],3), round(dt$pss[c],4),dt$kp[c], dt$N[c],round(Optimize$par,3), round(dk,3),  round(Refitted_R2,3),      sep = "\t")
      #oput <- paste("ID","Uniprot","Peptide","DP","z","mi","SS", "a","pss","kp","N", "k", "dk","R2", sep = "\t")
      write(oput, file=output_file, append=T)
    }
  }
  
}