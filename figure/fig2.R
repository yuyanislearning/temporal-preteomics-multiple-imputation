setwd('~/Projects/peipei/Turnover/Data/')
rm(list = ls())
library(tidyverse)
library(ggplot2)
library(data.table)

tt_d = read_delim('./all_hl-data.txt', delim='\t')
tt = read_delim('./all_hl.txt', delim = '\t')
tt_d = tt_d %>% group_by(ID) %>% summarise(count=n()) 
tt_d$uniprot = tt$Uniprot[match(tt_d$ID, tt$ID)]
tt_d = tt_d %>% group_by(uniprot) %>% summarise(count=max(count)) 

all_dat = data.frame(uniprot=NA,  count=NA, strain=NA, cond=NA)
for(strain in c('aj','balb','c57', 'cej','dba','fvb')){
  for(cond in c('iso', 'ctr')){
    tt_d = read_delim(paste0(strain,cond,'/hl-data.out'), delim='\t')
    tt = read_delim(paste0(strain,cond,'/hl.out'), delim = '\t')

    tt_d = tt_d %>% group_by(ID) %>% summarise(count=n()) 
    
    tt_d$uniprot = tt$Uniprot[match(tt_d$ID, tt$ID)]
    tt_d$R2 = tt$R2[match(tt_d$ID, tt$ID)]
    tt_d = tt_d %>% filter(R2>0.9) %>% group_by(uniprot) %>% summarise(count=max(count)) %>% 
      mutate(strain=strain, cond=cond)
    
    all_dat = rbind(all_dat, tt_d)
    
}}
all_dat = as.tibble(all_dat)
all_dat = all_dat[2:dim(all_dat)[1],]

# all proteins that will be imputated
temp3 = all_dat %>% 
  filter(count>1) %>%
  group_by(uniprot) %>%
  summarise(before=sum(count))

# all proteins that will not be imputated
temp2 = all_dat %>%
  filter(count<2) %>%
  group_by(uniprot) %>%
  summarise(before=sum(count)) %>%
  filter(!uniprot %in% unique(temp3$uniprot))

# impute the data :D
temp4 = all_dat %>%
  filter(count>1) %>%
  mutate(count=7 ) %>%
  group_by(uniprot) %>%
  summarise(after=sum(count)) %>%
  left_join(temp3, by = c('uniprot'='uniprot'))

temp4 = temp4 %>%
  arrange(desc(before))

# plot(1:dim(temp4)[1], temp4$before, col = 'red', type = 'l',
#      xlab = 'Proteins', ylab = '# Samples')
# lines(1:dim(temp4)[1], temp4$after, col = 'blue', type = 'l', lty='dashed')
# legend('bottomleft', legend=c('Before Imp', 'After Imp'), col = c('red','blue'), lty=c(1,2))
# polygon(c(1:dim(temp4)[1], rev(1:dim(temp4)[1])), c(temp4$before, rev(temp4$after)), col = "#6BD7AF")
# 
# blue <- rgb(0, 0, 1, alpha=0.5)
# red <- rgb(1, 0, 0, alpha=0.5)
# 
# barplot(temp4$before,1:dim(temp4)[1], space=0, col=red, border=NA)
# barplot(temp4$after,1:dim(temp4)[1], space=0,  col=blue,  border=NA, add=TRUE)


temp4$id = 1:dim(temp4)[1]
temp4$diff = temp4$after-temp4$before
temp4 = temp4 %>% select(id, diff,before)
names(temp4) = c('Proteins', 'After Imp', 'Before Imp')
df <- melt(temp4, id.vars = c('Proteins'),  value.name='Samples')
temp2 = temp2 %>% 
  arrange(desc(before)) 

temp2$variable = 'Before Imp'
names(temp2) = c('Proteins','Samples','variable')
temp2$Proteins = (length(unique(temp4$Proteins))+1):(length(unique(temp4$Proteins))+dim(temp2)[1])

temp2 = temp2 %>% select(Proteins, variable, Samples)

df = rbind(df, temp2)

ggplot(df, aes(x=Proteins, y=Samples, fill = variable)) + 
  geom_bar(stat='identity') +
  geom_vline(xintercept = 3214) +
   geom_vline(xintercept = 4469) 
  # geom_vline(xintercept = 4082)

print(dim(all_dat %>% filter(count>3) %>% group_by(uniprot) %>% summarise(count=n()))[1])
print(unique(length(temp3$uniprot)))
print(length(unique(df$Proteins)))





total_samples = length(unique(df$Proteins)) * 84
before_imp_samples = sum(df$Samples[df$variable=='Before Imp'])
missingness = 1 - before_imp_samples / total_samples
imputed_samples = sum(df$Samples[df$variable=='After Imp'])

imputed_samples / before_imp_samples



