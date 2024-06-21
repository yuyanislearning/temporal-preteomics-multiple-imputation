rm(list = ls())
setwd('~/Projects/peipei/Turnover/')
library(xlsx)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)



dat = read_csv('NSAF fractionate_DW_Processed.csv')

tmp = dat %>% select(day, strain, group, NSAF, `Uniprot ID`) %>%
  group_by(strain, group, `Uniprot ID`) %>%
  summarise(mean_nsaf = mean(NSAF))
tmp = tmp %>%
  mutate_at(c('strain', 'group'), .funs = toupper)



abun = read_csv('Strain_specific_turnover_Feb26.csv')

abun = abun %>%
  gather(key = 'con', value = 'prot') %>%
  separate(col = con, into = c('strain', 'group', 'impute'),sep = '_') %>%
  filter(!is.na(prot))

keep_prot = abun %>%
  group_by(strain, group, prot) %>%
  summarise(count = n())
abun = abun %>%
  left_join(keep_prot, by = c('strain', 'group', 'prot')) %>%
  filter(impute=='Before' | (impute=='After' & count==1))
  #filter(impute=='Before' | (impute=='After' & prot %in% keep_prot$prot))

abun = abun %>%
  mutate_at(c('strain'), .funs = toupper)

names(tmp)[3] = 'prot'
abun = abun %>%
  left_join(tmp , by = c('strain', 'group', 'prot'))

abun$impute = as.factor(abun$impute)
abun$impute = ordered(abun$impute, levels=c('Before','After'))


for(s in names(table(abun$strain))){
  pdf(paste0('figure/', s, '_violin.pdf'))
  print(ggplot(abun[abun$strain==s,], aes(x = group, y = log10(mean_nsaf), fill = impute)) +
    geom_violin(position=position_dodge(1), scale='count')+ 
      ylim(c(-5.5,-1.5)) +
      geom_signif(y_position=c(-1.5,-1.5), xmin = c(0.7,1.7), xmax=c(1.3,2.3),
                  annotation = c('***','***'), tip_length = 0))
  print(s)
  print('CTRL')
  a = wilcox.test(abun[abun$strain==s & abun$group=='CTRL' & abun$impute=='After',]$mean_nsaf,
              abun[abun$strain==s & abun$group=='CTRL' & abun$impute=='Before',]$mean_nsaf,)$p.value
  print(a)
  print('ISO')
  a=wilcox.test(abun[abun$strain==s & abun$group=='ISO' & abun$impute=='After',]$mean_nsaf,
              abun[abun$strain==s & abun$group=='ISO' & abun$impute=='Before',]$mean_nsaf,)$p.value
  print(a)
  dev.off()
}


ttt = abun %>% 
  group_by( prot, strain, group) %>%
  summarise(count = n()) %>%
  filter(count>1)

