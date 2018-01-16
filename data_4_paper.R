library(tidyverse)
library(broom)
library(stringr)
library(data.table)

knitr::opts_chunk$set(echo = FALSE, warning=FALSE, fig.align = 'center')

# load raw data (with preprocessing summaries)
exp_data <- readRDS('data/exp_data.rds')
exp_summary <- readRDS('data/exp_summary.rds')

# original data are separated for individual experiment
# here we combine them together
# 1. combine subject-wise summary table
sdata <- tibble(exp = 1:3, sdata = map(exp_summary,'sdata')) %>% 
  unnest() %>% # combine 3 exps together
  # make 3 experiments consistent with xfreq, target
  mutate(xfreq = ifelse(is.na(BlkType), 2, as.numeric(BlkType)),
         target = ifelse(is.na(target), as.character(dimension), as.character(target)),
         # change block frequency for 'absent' and 'orientation'
         BlkFreq = ifelse(target %in% c('Absent','Orientation') & exp < 3, 1-xfreq/4, xfreq/4),
         exp = factor(exp, levels = 1:3, labels = c('Exp. 1', 'Exp. 2', 'Exp. 3'))
         )

errors <- tibble(exp = 1:3, errors = map(exp_summary,'errors')) %>% unnest()  %>% 
  mutate(xfreq = ifelse(is.na(BlkType), 2, as.numeric(BlkType)),
         target = ifelse(is.na(target), as.character(dimension), as.character(target)),
         # change block frequency for 'absent' and 'orientation'
         BlkFreq = ifelse(target %in% c('Absent','Orientation') & exp < 3, 1-xfreq/4, xfreq/4),
         exp = factor(exp, levels = 1:3, labels = c('Exp. 1', 'Exp. 2', 'Exp. 3'))
  ) 

# 2. grand mean RTs and errors
ssdata <-  sdata %>% group_by(exp, target, BlkFreq) %>% 
  summarise(mmRT = mean(mRT), 
            n= n(), 
            seRT = sd(mRT)/sqrt(n-1)) 

serrors <- errors %>% group_by(exp, target, BlkFreq) %>% 
  summarise(merr = mean(err), 
            n= n(), 
            se_err = sd(err)/sqrt(n-1)) 

# 3. restructure the inter-trial effects
inttrial <-  tibble(exp = 1:3, data = map(exp_summary, 'inttrial')) %>% 
  unnest() %>%  #combine 3 exps together
  mutate(BlkType = ifelse(is.na(BlkType),2,BlkType), # adding missing in Exp. 3
         slope = 1/recinorm_sigma, # LATER slope
         intercept = recinorm_mu/recinorm_sigma, # LATER intercept
         resp = ifelse(d1==dn,'Repeat','Switch'), 
         dim = resp, 
         # split exp/ 3 response vs. dimension
         resp = ifelse(exp <3, resp, #remain the same for Exp. 1 & 2
                       ifelse((d1=='Absent' | dn == 'Absent') & d1 != dn, 'Switch','Repeat')),
         dim = ifelse(exp < 3, dim, 
                      ifelse((d1 != 'Absent' & dn!='Absent' & d1 == dn), 'Repeat','Switch')))
inttrial$exp <- factor(inttrial$exp, labels = c('Exp. 1', 'Exp. 2', 'Exp. 3'))
inttrial$BlkType <- factor(inttrial$BlkType, labels = c('1:3','1:1','3:1'))

# 4. grand mean inter-trial effects
mint <- inttrial %>% group_by(exp, dn, resp) %>% 
  summarise(mslope = mean(slope),mintercept = mean(intercept), n=n(),
            seslope = sd(slope)/sqrt(n-1), seintercept = sd(intercept)/sqrt(n-1)) 

# 5. combine inter-trial predictions based on adaptive model

# 5.1 model data and target names
pred_data = list(fname = c("data/Bayes_pred_e1.txt",
                           "data/Bayes_pred_e2.txt",
                           "data/Bayes_pred_e3_pa_att.txt"),
                 tar_name = c('target','dimension','dimension'))

# 5.2. a function to calculate inter-trial effects from prediction
iti_predict <- function(fname, tar_name){
  tsall <- read.table(fname) # read predicted data
  expno <- str_sub(fname,18,18) # get experiment no
  if (expno=='3'){
    tsall$BlkType <-  '1:1' # missing in Exp. 3
    blkTrls <- 56 # 56 trials per block in Exp. 3 
  } else {
    blkTrls <- 40 # 40 trials per block in Exp. 1 and 2 
  }
  
  tsall <-  tsall %>% 
    mutate_(dn = tar_name) %>% #current target, unify name to dn
    mutate(blksep = row_number()%%blkTrls !=1, # not first trial
           d1 = lag(dn)) # add n-1 trial condition
  tsall %>%  dplyr::filter(blksep) %>% # remove the first trial of each block
    group_by(sub, BlkType, dn, d1) %>% 
    summarise(pmRT = mean(rt, na.rm = TRUE)) %>%
    mutate(exp = paste('Exp.',expno)) # add Exp. #
}

# 5.3. Now combine 3 experiment data together with predictions
iti_preds <- map2_df(pred_data$fname,pred_data$tar_name,iti_predict) 
# due to use factors and characters in different codes, so convert all to character
iti_preds <- inttrial %>% map_if(is.factor, as.character) %>% # convert to chacter
  as.tibble() %>% left_join(., iti_preds, by = c('exp','sub','BlkType','d1','dn'))

# figure theme
figTheme <- function(pos) {
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     legend.title = element_blank(),legend.position = pos) }
