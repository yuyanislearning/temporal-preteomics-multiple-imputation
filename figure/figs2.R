rm(list=ls())
library(tidyverse)
library(ggplot2)
source('~/Projects/Tools/GO_analysis/GeneOntology.R')
setwd('~/Projects/peipei/Turnover/reactome/')

ctrlwimp = read_tsv('ctrl_after_bp.txt', skip=11)
names(ctrlwimp) = c('GO', 'ref_list', 'input_list', 
               'expected', 'over_under', 'fold_enrich',
               'raw_p', 'fdr')
ctrlwimp$fold_enrich = as.numeric(ctrlwimp$fold_enrich)

ctrlwoimp = read_tsv('ctrl_before_bp.txt',skip=11)
names(ctrlwoimp) = c('GO', 'ref_list', 'input_list', 
                'expected', 'over_under', 'fold_enrich',
                'raw_p', 'fdr')
ctrlwoimp$fold_enrich = as.numeric(ctrlwoimp$fold_enrich)

# show BP that only enrich after 
imp_only = setdiff(ctrlwimp$GO, ctrlwoimp$GO)
ctrl_imp_only = ctrlwimp %>%
  filter(GO %in% imp_only)
ctrl_imp_only = ctrl_imp_only %>% 
  filter(fdr<0.05 & fold_enrich>1) %>%
  arrange(fdr)


ctrl_imp_only$GO[1:20]
select_index = c(1, 3, 4:6, 9, 11, 14, 15, 20)
#assign_large_bp('0010506', BP_df$id)
#trace_ancesstry('0010506')
#select_size = 10
ctrl_imp_only = ctrl_imp_only[select_index,]
ctrl_imp_only = ctrl_imp_only[order(ctrl_imp_only$fold_enrich),]
ctrl_imp_only$GO = factor(ctrl_imp_only$GO, levels=ctrl_imp_only$GO)
p = ggplot(ctrl_imp_only,
       aes(y=GO, x=log2(fold_enrich), size=input_list,
           color=-log10(fdr))) +#dim(ctrl_imp_only)[1]
  geom_point(alpha=0.5) +
  scale_fill_continuous(limits=c(1,5)) +
  scale_size(range = c(5, 24), name='# proteins')

pdf('ctrl_imp.pdf',width = 10, height = 6)
print(p)
dev.off()

isowimp = read_tsv('iso_after_bp.txt', skip=11)
names(isowimp) = c('GO', 'ref_list', 'input_list', 
                    'expected', 'over_under', 'fold_enrich',
                    'raw_p', 'fdr')

isowoimp = read_tsv('iso_before_bp.txt',skip=11)
names(isowoimp) = c('GO', 'ref_list', 'input_list', 
                     'expected', 'over_under', 'fold_enrich',
                     'raw_p', 'fdr')

imp_only = setdiff(isowimp$GO, isowoimp$GO)
iso_imp_only = isowimp %>%
  filter(GO %in% imp_only)
iso_imp_only = iso_imp_only %>% 
  filter(fdr<0.05 & fold_enrich>1) %>%
  arrange(fdr)


iso_imp_only$GO[1:20]
select_index = c(1:5, 7, 8, 14, 15, 19)
#assign_large_bp('0010506', BP_df$id)
#trace_ancesstry('0010506')
#select_size = 10
iso_imp_only = iso_imp_only[select_index,]
iso_imp_only = iso_imp_only[order(iso_imp_only$fold_enrich),]
iso_imp_only$GO = factor(iso_imp_only$GO, levels=iso_imp_only$GO)
p = ggplot(iso_imp_only,
           aes(y=GO, x=log2(fold_enrich), size=input_list,
               color=-log10(fdr))) +#dim(iso_imp_only)[1]
  geom_point(alpha=0.5) +
  scale_fill_continuous(limits=c(1,5)) +
  scale_size(range = c(5, 24), name='# proteins')

pdf('iso_imp.pdf',width = 7, height = 6)
print(p)
dev.off()

ctrl_imp_only$group = 'control group'
iso_imp_only$group = 'ISO group'
all = rbind(ctrl_imp_only, iso_imp_only)

p = ggplot(all,
           aes(y=GO, x=log2(fold_enrich), size=input_list,
               color=-log10(fdr))) +#dim(iso_imp_only)[1]
  geom_point(alpha=0.5) +
  scale_fill_continuous(limits=c(1,5)) +
  scale_size(range = c(5, 24), name='# proteins') +
  facet_wrap(vars(group), scales = 'free_y')

