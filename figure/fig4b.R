rm(list=ls())
setwd('~/Projects/peipei/Turnover/complex/')
library(tidyverse)
library(ggplot2)
library(ggsignif)


set.seed(123)
dat_iso = read_csv('iso_complexome_turnovers.csv')
dat_ctrl = read_csv('ctrl_turnovers_for_heterocomplexes_shared_with_iso.csv')

strains = c('dba','fvb','cej','c57', 'balbc','aj')
conditions = c('ctrl', 'iso')
imputations = c('before', 'after')

imp_interest = 'after'

all_data <- list()
for(strain in strains){
  for(cond in conditions){
    for(imp in imputations){
      dat = read_csv(paste0('/Users/yuyan/Library/Mobile Documents/com~apple~CloudDocs/Projects/peipei/Turnover/Data/processed/res_/',
                            paste(cond, strain, imp,sep='_'), '.csv'))
      dat$cond = cond
      dat$strain = strain
      dat$imp = imp
      all_data[[length(all_data) + 1]] <- dat
    }
  }
}
combined_data <- do.call(rbind, all_data)

quick_func = function(dat_iso, combined_data, cond_interest){
  iso_after_summary = dat_iso %>% filter(imp==imp_interest) %>%
    group_by(Complex) %>% summarise(std=sd(AvgTurnover))
  complex_count = dat_iso %>% filter(imp==imp_interest) %>% group_by(Complex) %>%
    summarise(count=n()) 
  complex_freq = table(complex_count$count)/sum(complex_count$count)
  iso_boot = rep(0, 100)
  
  for(i in 1:100){
    comb_iso = combined_data %>% filter(cond==cond_interest & imp==imp_interest) 
    n_sample = sample(as.integer(names(complex_freq)), 1, 
                      prob = complex_freq)
    comb_iso = comb_iso %>% sample_n(n_sample) %>%
      group_by(cond) %>% summarise(std = sd(log2(median.K)))
    iso_boot[i] = comb_iso$std
  }
  
  iso_after_summary$group = 'Protein Complex'
  iso_after_summary = iso_after_summary %>% select(std, group)
  
  iso_boot = data.frame(group='Proteome', std=iso_boot)
  return(list('iso'=iso_after_summary, 'boot'=iso_boot))
}

quick_res = quick_func(dat_iso, combined_data, 'iso')
iso_after_summary = quick_res$iso
iso_boot = quick_res$boot
iso_after_summary$cond = 'ISO'
iso_boot$cond = 'ISO'
df = rbind(iso_after_summary, iso_boot)

quick_res = quick_func(dat_ctrl, combined_data, 'ctrl')
iso_after_summary = quick_res$iso
iso_boot = quick_res$boot
iso_after_summary$cond = 'CTRL'
iso_boot$cond = 'CTRL'
df = rbind(df, iso_after_summary)
df = rbind(df, iso_boot)

df_ctrl_com = df %>% filter(cond=='CTRL' & group == 'Protein Complex')
df_ctrl_prot = df %>% filter(cond=='CTRL' & group == 'Proteome')
wilcox.test(df_ctrl_com$std, df_ctrl_prot$std)

df_ctrl_com = df %>% filter(cond=='ISO' & group == 'Protein Complex')
df_ctrl_prot = df %>% filter(cond=='ISO' & group == 'Proteome')
wilcox.test(df_ctrl_com$std, df_ctrl_prot$std)

df = df %>% filter(cond=='CTRL')

g = ggplot(df, aes(x=group, y=std, fill=group)) +
  geom_violin(position=position_dodge(1)) +
  facet_wrap(~cond, scales = "free_x") +
  geom_signif(y_position=c(4), xmin = c(1), xmax=c(2),
              annotation = c('***'), tip_length = 0) +
  labs(title = "Stability of Protein Complex vs Proteome", 
       x = "Group", y = "Protein Complex stability (Standard deviation of turnover)")
g
















