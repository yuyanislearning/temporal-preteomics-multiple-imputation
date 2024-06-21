rm(list = ls(all = T)) 


library(data.table)
dat <- fread("hl.out")
dat1 <- fread("hl-data.out")

dat1<-dat1[,c(1,2,3)]


#dat_fil <- subset(dat,   DP>=4 )
#median_hl_sage <- aggregate(dat_fil$k,by = list(dat_fil$Uniprot),median) #489

#######howard

library(tidyverse)

tmp_dat<-read_csv('HU4.csv', col_names =F)

tmp_dat=tmp_dat %>% filter(!is.na(X2))

sum(!is.na(match(tmp_dat$X1, (unique(dat_fil$Uniprot)) )))
257/377
####################verify previous results
dat_sage <- fread("hl.out")
dat_sage <- subset(dat_sage,   DP>=4 & (R2>0.8 | SS<0.05 ))
#dat_sage <- subset(dat_sage,   DP>=4 & R2>0.8 )
median_hl_sage <- aggregate(dat_sage$k,by = list(dat_sage$Uniprot),median) #489

####################verify

unique_A0 <- aggregate(dat1$A0,by = list(dat1$t,dat1$ID),mean)
colnames(unique_A0)<-c("t","ID","A0")
library(tidyr)
data <- spread(unique_A0, t, A0) 


data[data < 0]<-NA #REMOVE - halflife
minv<-min(data, na.rm = TRUE)

cnt_na <- apply(data, 1, function(z) sum(is.na(z)))

data_filt<-data[cnt_na < 9,]  # 14 FOR 15 tp...9 for 10 tps
#data_filt<-data[cnt_na < 5,]  # male is one less 
data_filt1<-subset(data_filt,select=-c(ID))
#data_filt1[data_filt1 > 1]<-1

mis_per<-sum(is.na(data_filt1))/prod(dim(data_filt1)) #37% missing aj ctrl
#originially >56% data is missing (at least 1 tp)
#at least 1 tp..56% missing  (35945 observaions total)
#at least 2 tp..40.4% missing (23198)
#at least 3 tp..30%   missing (17372)
#at least 4 tp..21% missing (13101)


#x<-abs(data_filt1)
#x<-(data_filt1)
#minv<-min(x, na.rm = TRUE)
#maxv<-max(x, na.rm = TRUE)


library(mice)
mice_impute <- mice(data_filt1,m=10,method=c('pmm','pmm','pmm','pmm','pmm',
                                             'pmm','pmm','pmm','pmm','pmm'), maxiter=10)     #female
#mice_impute <- mice(data_filt1,m=10,method=c('pmm','pmm','pmm','pmm','pmm',
#                                            'pmm'), maxiter=10)     #male

summary(mice_impute)
data_all <- complete(mice_impute, "long")
write.csv(data_all,"pom_hl_m_mice5.csv")


for(i in 1:10) {
  
  mice_imp<- complete(mice_impute,i)
  mice_imp<-cbind(data_filt[,c("ID")],mice_imp)
  #colnames(mice_imp)[1] <- "ID"
  colnames(mice_imp) <- c("ID","0","1","2","4","5","8","9","11","12","14")
  #colnames(mice_imp) <- c("ID","d00","d03","d05","d07","d10","d14") #male
  long_data<-gather(mice_imp, key = t, value = A0,
                    "0","1","2","4","5","8","9","11","12","14")
  long_data <- long_data[order(long_data$ID),] 
  
  name <- paste("hs4_imp", i,"out", sep = ".")
  write.table(long_data, name, sep="\t",row.names = F)
  #write.csv(mice_imp,name)
  
}




long_data<-gather(mice_imp, key = t, value = A0,
                  d00,d03,d05,d07,d10,d14)
long_data <- long_data[order(long_data$ID),] 

write.table(long_data, "mydata.out", sep="\t",row.names = F)

for(i in 1:5) {
  
  mice_imp<- complete(mice_impute,i)
  mice_imp<-cbind(data_filt[,c("ID")],mice_imp)
  #colnames(mice_imp)[1] <- "ID"
  #colnames(mice_imp) <- c("ID","d00","d01","d03","d05","d07","d10","d14")
  colnames(mice_imp) <- c("ID","d00","d03","d05","d07","d10","d14") #male
  name <- paste("cm_m_hl_imp", i,"csv", sep = ".")
  write.csv(mice_imp,name)
  
}


library(Amelia)
#amelia_impute <- amelia(data_filt1, m = 5)
#all<-as.data.frame(amelia_impute$imputations)
bds <- matrix(c(1, 2, 3, 4,5,6,7,0,0,0,0,0,0,0,1 ,1,1,1,1,1,1), nrow = 7, ncol = 3) #getting 0 to 1 bound for each varaible
amelia_impute <- amelia(data_filt1, m = 5,  bounds = bds, max.resample = 1000)

j=1
for(i in 1:5) {
  amel_imp<- all[,c(j,j+1,j+2,j+3,j+4,j+5,j+6)]
  
  #amel_imp[amel_imp < 0] <- minv
  #amel_imp[amel_imp > 1] <- 1
  amel_imp<-cbind(data_filt[,c("ID")], amel_imp)
  colnames(amel_imp) <- c("ID","X0","X1", "X3","X5","X7","X10","X14")
  name <- paste("amel_imp", i,"csv", sep = ".")
  write.csv(amel_imp,name)
  j=j+7
  
}


library(Amelia)
amelia_impute <- amelia(data_filt1, m = 5)
all<-as.data.frame(amelia_impute$imputations)

j=1
for(i in 1:5) {
  amel_imp<- all[,c(j,j+1,j+2,j+3,j+4,j+5,j+6)]
  
  #amel_imp[amel_imp < 0] <- minv
  #amel_imp[amel_imp > 1] <- 1
  amel_imp<-cbind(data_filt[,c("ID")], amel_imp)
  colnames(amel_imp) <- c("ID","X0","X1", "X3","X5","X7","X10","X14")
  name <- paste("amel_imp", i,"csv", sep = ".")
  write.csv(amel_imp,name)
  j=j+7
  
}

































long_data<-gather(mice_imp, key = t, value = A0,
                  d00,d03,d05,d07,d10,d14)
long_data <- long_data[order(long_data$ID),] 

write.table(long_data, "mydata.out", sep="\t",row.names = F)


library(data.table)
dat_noimp <- fread("hl.out")
dat_imp <- fread("ProTurn_Output_imp1_fcs.txt")


#newdata <- subset(dat, R2 >= 0.8 & DP >= 4)
newdata_noimp <- subset(dat_noimp, R2 >= 0.8 & DP >= 4)
median_hl <- aggregate(newdata_noimp$k,by = list(newdata_noimp$Uniprot),median) #1302
colnames(median_hl) <- c("UniProt","median.K")

newdata_imp <- subset(dat_imp, R2 >= 0.8 )
median_hl_imp <- aggregate(newdata_imp$k,by = list(newdata_imp$Uniprot),median) #1428
















