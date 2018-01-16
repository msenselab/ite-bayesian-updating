library(tidyverse)
library(data.table)
library(cowplot) 

source('codes/bayesian_updates.R')

# first subject with 54 models, run for 347 seconds (~ 6 minutes)
exp1_pars <- readRDS('data/exp1_fit_pars_v4.rds')
exp2_pars <- readRDS('data/exp2_fit_pars_v4.rds')
exp3_pars <- readRDS('data/exp3_fit_pars_v4.rds')

subs_e1 <- c('aleg', 'alis', 'elis', 'jarw', 'julm', 'katm', 'lens', 'marb', 'nadf', 'olir', 'ricg', 'thep')
row.names(exp1_pars) <- apply(expand.grid(subs_e1, 1:144), 1, paste, collapse=" ")
subs_e2 <- c('adrl', 'aleg', 'anas', 'elis', 'fran', 'jons', 'julk', 'julm', 'lens', 'marm', 'olir', 'ricg')
row.names(exp2_pars) <- apply(expand.grid(subs_e2, 1:144), 1, paste, collapse=" ")
subs_e3 <- c('aleg', 'elis', 'fra', 'grel', 'irij', 'jul', 'julm', 'lens', 'maxw', 'olivr', 'ricg', 'rosr')
row.names(exp3_pars) <- apply(expand.grid(subs_e3, 1:144), 1, paste, collapse=" ")


names(exp1_pars) <- c("mem_resp", "beta_resp", "mem_dim", "beta_dim", "scale_resp", "u_mem_resp", "u_delta_resp", "scale_dim", "u_mem_dim", 
                      "u_delta_dim", "ter", "later_diffusion", "resp_update", "dim_update", "theta1", "mu1", "sig1", "theta2", "mu2", "sig2", 
                      "theta3","mu3", "sig3", "nll", "AIC", "BIC", "sub", "model")
names(exp2_pars) <- c("mem_resp", "beta_resp", "mem_dim", "beta_dim", "scale_resp", "u_mem_resp", "u_delta_resp", "scale_dim", "u_mem_dim", 
                      "u_delta_dim", "ter", "later_diffusion", "resp_update", "dim_update", "theta1", "mu1", "sig1", "theta2", "mu2", "sig2", 
                      "nll", "AIC", "BIC", "sub", "model")
names(exp3_pars) <- c("mem_resp", "beta_resp", "mem_dim", "beta_dim", "scale_resp", "u_mem_resp", "u_delta_resp", "scale_dim", "u_mem_dim", 
                      "u_delta_dim", "ter", "later_diffusion", "resp_update", "dim_update", "theta1", "mu1", "sig1", "theta2", "mu2", "sig2", 
                      "theta3","mu3", "sig3", "nll", "AIC", "BIC", "sub", "model")

exp1_pars <- dplyr::filter(exp1_pars, dim_update != 1) # Remove dimension updating based on S0 without forgetting since this was 
exp1_pars <- dplyr::filter(exp1_pars, later_diffusion<=1) # causing problems. 
exp2_pars <- dplyr::filter(exp2_pars, dim_update != 1) # Note: Some of the filtering criteria look different from simply removing one
exp2_pars <- dplyr::filter(exp2_pars, as.numeric(nll)<1e9) # form of dimension updating. This is necessary because things have ended
exp3_pars <- dplyr::filter(exp3_pars, dim_update != 1) # up in wrong columns...
exp3_pars <- dplyr::filter(exp3_pars, !is.na(theta1)) 

#exp1_pars$s0_update <- factor(exp1_pars$s0_update, levels=1:3, labels=c("No update", "Simulation", "Analytical"))
#exp2_pars$s0_update <- factor(exp2_pars$s0_update, levels=1:3, labels=c("No update", "Simulation", "Analytical"))
#exp3_pars$s0_update <- factor(exp3_pars$s0_update, levels=1:3, labels=c("No update", "Simulation", "Analytical"))
exp1_pars$later_diffusion <- factor(exp1_pars$later_diffusion, levels=0:1, labels=c("DDM", "LATER"))
exp2_pars$later_diffusion <- factor(exp2_pars$later_diffusion, levels=0:1, labels=c("DDM", "LATER"))
exp3_pars$later_diffusion <- factor(exp3_pars$later_diffusion, levels=0:1, labels=c("DDM", "LATER"))
exp1_pars$resp_update <- factor(exp1_pars$resp_update, levels=0:5, labels=c("No update", "S0 with full memory", "S0 with decay", "Binary rate", 
                                                                            "Rate with decay", "Weighted rate"))
exp2_pars$resp_update <- factor(exp2_pars$resp_update, levels=0:5, labels=c("No update", "S0 with full memory", "S0 with decay", "Binary rate", 
                                                                            "Rate with decay", "Weighted rate"))
exp3_pars$resp_update <- factor(exp3_pars$resp_update, levels=0:5, labels=c("No update", "S0 with full memory", "S0 with decay", "Binary rate", 
                                                                            "Rate with decay", "Weighted rate"))
exp1_pars$dim_update <- factor(exp1_pars$dim_update, levels=0:5, labels=c("No update", "S0 with full memory", "S0 with decay", "Binary rate", 
                                                                          "Rate with decay", "Weighted rate"))
exp2_pars$dim_update <- factor(exp2_pars$dim_update, levels=0:5, labels=c("No update", "S0 with full memory", "S0 with decay", "Binary rate", 
                                                                          "Rate with decay", "Weighted rate"))
exp3_pars$dim_update <- factor(exp3_pars$dim_update, levels=0:5, labels=c("No update", "S0 with full memory", "S0 with decay", "Binary rate", 
                                                                          "Rate with decay", "Weighted rate"))

# check update on simulation/analytical and no update

# s0 update with simulation or analytical did not make any difference. 
# Both are significant better than no-update 
# rate update with drift-like works best
# exp1_pars %>% mutate(beta_fix = beta_a == 1) %>% group_by(s0_update, u_update, beta_fix) %>% 
#   summarise(AIC = mean(AIC)) %>% 
#   ggplot(aes(x = s0_update, y = AIC, group = as.factor(u_update), color = as.factor(u_update))) + 
#   geom_line() + facet_wrap(~ beta_fix)
# 
# exp1_pars %>% mutate(mem_fix = w_m == 1) %>% group_by(s0_update, u_update, isLater) %>% 
#   summarise(AIC = mean(AIC)) %>% 
#   ggplot(aes(x = s0_update, y = AIC, group = as.factor(u_update), color = as.factor(u_update))) + 
#   geom_line() + geom_point(size=3) + facet_wrap(~ isLater) + theme_bw() + 
#   labs(x = "S0 updating", y="AIC", color="Rate updating") + theme(text = element_text(size=20))
# 
# exp2_pars %>% mutate(mem_fix = w_m == 1) %>% group_by(s0_update, u_update, isLater) %>% 
#   summarise(AIC = mean(AIC)) %>% 
#   ggplot(aes(x = s0_update, y = AIC, group = as.factor(u_update), color = as.factor(u_update))) + 
#   geom_line() + geom_point(size=3) + facet_wrap(~ isLater) + theme_bw() + 
#   labs(x = "S0 updating", y="AIC", color="Rate updating") + theme(text = element_text(size=20))
# 
# exp3_pars %>% mutate(mem_fix = w_m == 1) %>% filter(nll<1e8) %>% group_by(s0_update, u_update, isLater) %>% 
#   summarise(AIC = mean(AIC)) %>% 
#   ggplot(aes(x = s0_update, y = AIC, group = as.factor(u_update), color = as.factor(u_update))) + 
#   geom_line() + geom_point(size=3) + facet_wrap(~ isLater) + theme_bw() + 
#   labs(x = "S0 updating", y="AIC", color="Rate updating") + theme(text = element_text(size=20))

fig_aic1 <- exp1_pars %>% dplyr::filter(substr(model,3,3)=='V') %>% group_by(resp_update, dim_update, later_diffusion) %>% 
  summarize(mAIC = mean(as.numeric(AIC))) %>% 
  ggplot(aes(x = resp_update, y=mAIC, group=dim_update, color=dim_update)) + theme_bw() + facet_wrap(~ later_diffusion) + 
  geom_line() + geom_point(size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x="RDF based updating", y="average AIC", color="TDD based updating")
fig_aic2 <- exp2_pars %>% dplyr::filter(substr(model,3,3)=='V') %>% group_by(resp_update, dim_update, later_diffusion) %>% 
  summarize(mAIC = mean(as.numeric(AIC))) %>% 
  ggplot(aes(x = resp_update, y=mAIC, group=dim_update, color=dim_update)) + theme_bw() + facet_wrap(~ later_diffusion) + 
  geom_line() + geom_point(size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x="RDF based updating", y="average AIC", color="TDD based updating")
fig_aic3 <- exp3_pars %>% dplyr::filter(substr(model,3,3)=='V') %>% group_by(resp_update, dim_update, later_diffusion) %>% 
  summarize(mAIC = mean(as.numeric(AIC))) %>% 
  ggplot(aes(x = resp_update, y=mAIC, group=dim_update, color=dim_update)) + theme_bw() + facet_wrap(~ later_diffusion) + 
  geom_line() + geom_point(size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x="RDF based updating", y="average AIC", color="TDD based updating")

fig_noNDT_aic1 <- exp1_pars %>% dplyr::filter(substr(model,3,3)=='0') %>% group_by(resp_update, dim_update, later_diffusion) %>% 
  summarize(mAIC = mean(as.numeric(AIC))) %>% 
  ggplot(aes(x = resp_update, y=mAIC, group=dim_update, color=dim_update)) + theme_bw() + facet_wrap(~ later_diffusion) + 
  geom_line() + geom_point(size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x="RDF based updating", y="average AIC", color="TDD based updating")
fig_noNDT_aic2 <- exp2_pars %>% dplyr::filter(substr(model,3,3)=='0') %>% group_by(resp_update, dim_update, later_diffusion) %>% 
  summarize(mAIC = mean(as.numeric(AIC))) %>% 
  ggplot(aes(x = resp_update, y=mAIC, group=dim_update, color=dim_update)) + theme_bw() + facet_wrap(~ later_diffusion) + 
  geom_line() + geom_point(size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x="RDF based updating", y="average AIC", color="TDD based updating")
fig_noNDT_aic3 <- exp3_pars %>% dplyr::filter(substr(model,3,3)=='0') %>% group_by(resp_update, dim_update, later_diffusion) %>% 
  summarize(mAIC = mean(as.numeric(AIC))) %>% 
  ggplot(aes(x = resp_update, y=mAIC, group=dim_update, color=dim_update)) + theme_bw() + facet_wrap(~ later_diffusion) + 
  geom_line() + geom_point(size=3) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(x="RDF based updating", y="average AIC", color="TDD based updating")

# 
# # check beta_a, some benefits (D_AIC = 3)
# exp1_pars %>% filter(s0_update > 1) %>% group_by( beta_a == 1) %>% 
#   summarise(AIC = mean(AIC))
# 
# # check u_update within s0_update
# # no update is a bit better than step update (d_AIC = 2),
# # but significant better with drift update (d_AIC = 20)
# exp1_pars %>% filter(s0_update > 1) %>% group_by(u_update) %>% 
#   summarise(AIC = mean(AIC))
# 
# exp1_pars %>% group_by(sub) %>% filter(AIC == min(AIC)) -> exp1_best
# # best model(s) for Exp.1: LS_MVT0R1
# # best model(s): LS_M1BVT0R3, LS_MVB1T0R3

# Create lists of models, sorted by total AIC
group_by(exp1_pars, model) %>% summarize(totAIC=sum(as.numeric(AIC)), totBIC=sum(as.numeric(BIC))) %>% arrange(totAIC) -> model_list_e1
group_by(exp2_pars, model) %>% summarize(totAIC=sum(as.numeric(AIC)), totBIC=sum(as.numeric(BIC))) %>% arrange(totAIC) -> model_list_e2
group_by(exp3_pars, model) %>% summarize(totAIC=sum(as.numeric(AIC)), totBIC=sum(as.numeric(BIC))) %>% arrange(totAIC) -> model_list_e3

# Create lists of best and worst models based on total AIC
model_list_e1[1:5,] -> top5_models_e1
model_list_e1[104:108,] -> bottom5_models_e1
model_list_e2[1:5,] -> top5_models_e2
model_list_e2[104:108,] -> bottom5_models_e2
model_list_e3[1:5,] -> top5_models_e3
model_list_e3[104:108,] -> bottom5_models_e3

# Top 5 lists are all based on the LATER model rather than the DDM and have updating of the prior with forgetting (MV not M1).
# For exp 1 all top models have beta functions parameter as a free parameter, but for exp 2 and 3 there is no clear 
# preference either way. Neither the analytical nor the simulation based updating is clearly preferred over the other. 
# Including a non-decision time parameter in the LATER model seems to have been beneficial in exps 2 and 3 but there is not
# such a clear preference in exp 1. For the rate updating, the memory model seems to work best for exps 1 and 2 while the weight
# model is better for exp 3.

# Top lists based on number of participants for which a model is the best fitting one in terms of AIC
topmod_e1 = list()
topmod_e2 = list()
topmod_e3 = list()

for (i in 1:12) { dplyr::filter(exp1_pars, sub==i) %>% arrange(AIC) %>% dplyr::filter(AIC==min(AIC)) %>% dplyr::select(model) -> topmod_e1[i] }
for (i in 1:12) { dplyr::filter(exp2_pars, sub==i) %>% arrange(AIC) %>% dplyr::filter(AIC==min(AIC)) %>% dplyr::select(model) -> topmod_e2[i] }
for (i in 1:12) { dplyr::filter(exp3_pars, sub==i) %>% arrange(AIC) %>% dplyr::filter(AIC==min(AIC)) %>% dplyr::select(model) -> topmod_e3[i] }


# Here we select the best model for each experiment by for each factor picking the level of that factor that provides the best fit for the most
# subjects. This leads to the same outcome as picking the model with the lowest average (or equivalently total) AIC in all respects except that
# the DR step model is favored over DR weight for dimension based updating in experiment 3. 

data_e1 = readRDS('data/exp1.rds') %>% mutate(targets = target == 'Present') %>% 
  mutate(inttrial_resp=ifelse((target)==lag(target), 1, 0), 
         targets = ifelse(targets, as.numeric(dimension), 0), 
         inttrial_dim=ifelse((dimension)==lag(dimension), 1, 0)) 
data_e1$inttrial_dim[data_e1$targets==FALSE]<-NA
data_e1$inttrial_dim[lag(data_e1$targets)==FALSE]<-NA
data_e1 <- group_by(data_e1, sub) %>% 
  mutate(outlier = rt < (quantile(rt)[2] - 1.5*(quantile(rt)[4] - quantile(rt)[2])) | rt > 2 | c(1,diff(blkNo)))

data_e1$dimension <- as.numeric(data_e1$dimension)
data_e1$dimension[data_e1$target=="Absent"] <- 0
top_model_e1 <- top5_models_e1$model[1]
model_params <- dplyr::filter(exp1_pars, model==top_model_e1)
model_params$resp_update <- as.numeric(model_params$resp_update)-1
model_params$dim_update <- as.numeric(model_params$dim_update)-1
subs <- c('aleg', 'alis', 'elis', 'jarw', 'julm', 'katm', 'lens', 'marb', 'nadf', 'olir', 'ricg', 'thep')
model_params$sub <- subs
data_pred_e1 <- data.frame(sub=character(0))
for(s in subs) {
  dp <- genSeq(dplyr::filter(model_params, sub==s), dplyr::filter(data_e1, sub==s))
  dp$tno <- 1:nrow(dp)
  data_pred_e1 <- full_join(data_pred_e1, dp)
}

group_by(data_pred_e1, blkNo, sub) %>% mutate(pt=lag(target)) %>% dplyr::filter(error==FALSE, outlier==FALSE, !is.na(pt)) %>% 
  group_by(pt, target, BlkType, sub) %>% summarize(mRT=mean(rt), mpRT=mean(predRT)) -> data_pred_summary_e1

data_e2 = readRDS('data/exp2.rds') %>% mutate(targets = dimension == 'Color') %>% 
  mutate(inttrial_dim=ifelse((dimension)==lag(dimension), 1, 0),
         inttrial_resp=inttrial_dim) 
data_e2 <- group_by(data_e2, sub) %>% 
  mutate(outlier = rt < (quantile(rt)[2] - 1.5*(quantile(rt)[4] - quantile(rt)[2])) | rt > 2 | c(1,diff(blkNo))) 
top_model_e2 <- top5_models_e2$model[1]
model_params <- dplyr::filter(exp2_pars, model==top_model_e2)
model_params$resp_update <- as.numeric(model_params$resp_update)-1
model_params$dim_update <- as.numeric(model_params$dim_update)-1
subs <- c('adrl', 'aleg', 'anas', 'elis', 'fran', 'jons', 'julk', 'julm', 'lens', 'marm', 'olir', 'ricg')
model_params$sub <- subs
data_pred_e2 <- data.frame(sub=character(0))
for(s in subs) {
  dp <- genSeq(dplyr::filter(model_params, sub==s), dplyr::filter(data_e2, sub==s))
  dp$tno <- 1:nrow(dp)
  data_pred_e2 <- full_join(data_pred_e2, dp)
}

group_by(data_pred_e2, blkNo, sub) %>% mutate(pd=lag(dimension)) %>% dplyr::filter(error==FALSE, outlier==FALSE, !is.na(pd)) %>% 
  group_by(pd, dimension, BlkType, sub) %>% summarize(mRT=mean(rt), mpRT=mean(predRT)) -> data_pred_summary_e2

data_e3 = readRDS('data/exp3.rds') %>% mutate(targets = as.numeric(dimension)-1, 
                                              resp = ifelse(targets>0, 1, 0),
                                              inttrial_resp = ifelse(resp == lag(resp) ,1, 0),
                                              inttrial_dim=ifelse((dimension)==lag(dimension), 1, 0))
data_e3 <- group_by(data_e3, sub) %>% 
  mutate(outlier = rt < (quantile(rt)[2] - 1.5*(quantile(rt)[4] - quantile(rt)[2])) | rt > 2 | c(1,diff(blkNo))) 
data_e3$inttrial_dim[data_e3$targets==FALSE]<-NA
data_e3$inttrial_dim[lag(data_e3$targets)==FALSE]<-NA
data_e3$dimension <- as.numeric(data_e3$dimension) - 1

top_model_e3 <- top5_models_e3$model[1]
model_params <- dplyr::filter(exp3_pars, model==top_model_e3)
model_params$resp_update <- as.numeric(model_params$resp_update)-1
model_params$dim_update <- as.numeric(model_params$dim_update)-1
subs <- c('aleg', 'elis', 'fra', 'grel', 'irij', 'jul', 'julm', 'lens', 'maxw', 'olivr', 'ricg', 'rosr')
model_params$sub <- subs
data_pred_e3 <- data.frame(sub=character(0))
for(s in subs) {
  dp <- genSeq(dplyr::filter(model_params, sub==s), dplyr::filter(data_e3, sub==s))
  dp$tno <- 1:nrow(dp)
  data_pred_e3 <- full_join(data_pred_e3, dp)
}

group_by(data_pred_e3, blkNo, sub) %>% mutate(pd=lag(dimension)) %>% dplyr::filter(error==FALSE, outlier==FALSE, !is.na(pd)) %>% 
  group_by(pd, dimension, sub) %>% summarize(mRT=mean(rt), mpRT=mean(predRT)) -> data_pred_summary_e3

data_pred_summary_e1$exp <- 1
data_pred_summary_e2$exp <- 2
data_pred_summary_e3$exp <- 3
data_pred_summary_e1 <- ungroup(data_pred_summary_e1)
data_pred_summary_e2 <- ungroup(data_pred_summary_e2)
data_pred_summary_e3 <- ungroup(data_pred_summary_e3)
data_pred_summary <- rbind(dplyr::select(data_pred_summary_e1, exp, mRT, mpRT), dplyr::select(data_pred_summary_e2, exp, mRT, mpRT), 
                           dplyr::select(data_pred_summary_e3, exp, mRT, mpRT))
data_pred_summary$exp <- factor(data_pred_summary$exp, levels=1:3, labels=c("Exp. 1", "Exp. 2", "Exp. 3"))

data_vs_pred_plot <- ggplot(data_pred_summary, aes(x = mRT, y = mpRT, color=exp)) +
  geom_point() + 
  geom_smooth(method = 'lm', fullrange = TRUE, se = FALSE) +  
  labs(x="measured mean RT", y="predicted mean RT") +
  coord_fixed() +
  scale_x_log10(breaks = seq(0.4,1, by=0.1)) + 
  scale_y_log10(breaks = seq(0.4,1, by=0.1)) + theme_bw()

sample_sub_e1 <- 1
sample_sub_e2 <- 3
sample_sub_e3 <- 2

e1_sample <- dplyr::filter(data_pred_e1, sub==subs_e1[sample_sub_e1], BlkType=="3:1")[11:110,] %>% mutate(s0=s0 * ((targets>0) * 2 - 1)) %>% 
  mutate(tar_start = tno, tar_end = tno+1)

p_e1_s0 <- e1_sample %>% ggplot(aes(x=tno, y=s0))  + geom_line() + theme_bw() +
  geom_rect( aes(xmin = tar_start, xmax = tar_end,x = NULL, y = NULL, ymin = -0.3, ymax = 0.3, fill = target), alpha=0.5) + 
  xlab('Trial no.') + theme(legend.position="bottom")

rates_e1 <- c(dplyr::filter(exp1_pars, sub==sample_sub_e1, model==top_model_e1)$mu1, 
              dplyr::filter(exp1_pars, sub==sample_sub_e1, model==top_model_e1)$mu2, 
              dplyr::filter(exp1_pars, sub==sample_sub_e1, model==top_model_e1)$mu3) 

p_e1_rate <- e1_sample %>% ggplot(aes(x=tno, y=rate)) + geom_line() + theme_bw()  +
  geom_rect( aes(xmin = tar_start, xmax = tar_end,x =NULL, y = NULL, ymin = 3, ymax = 5, fill = dimension), alpha=0.5) + 
  xlab('Trial no.') + geom_hline(yintercept = rates_e1[1], linetype="dashed", color="red") +
  geom_hline(yintercept = rates_e1[2], linetype="dashed", color="darkgreen") +
  geom_hline(yintercept = rates_e1[3], linetype="dashed", color="blue" ) + theme(legend.position="bottom")

e2_sample <- dplyr::filter(data_pred_e2, sub==subs_e2[sample_sub_e2], BlkType=="1:3")[11:110,] %>% mutate(s0=s0 * ((targets>0) * 2 - 1)) %>%
  mutate(tar_start = tno, tar_end = tno+1)

p_e2_s0 <- e2_sample %>% ggplot(aes(x=tno, y=s0)) + geom_line() + theme_bw() +
  geom_rect( aes(xmin = tar_start, xmax = tar_end,x =NULL, y = NULL, ymin = -1.2, ymax = 1.2, fill = dimension), alpha=0.5) + 
  xlab('Trial no.') + theme(legend.position="bottom")

rates_e2 <- c(dplyr::filter(exp2_pars, sub==sample_sub_e2, model==top_model_e2)$mu1, 
              dplyr::filter(exp2_pars, sub==sample_sub_e2, model==top_model_e2)$mu2, 
              dplyr::filter(exp2_pars, sub==sample_sub_e2, model==top_model_e2)$mu3) 

p_e2_rate <- e2_sample %>% ggplot(aes(x=tno, y=rate)) + geom_line() + theme_bw()  +
  geom_rect( aes(xmin = tar_start, xmax = tar_end,x =NULL, y = NULL, ymin = 11, ymax = 15, fill = dimension), alpha=0.5) + 
  xlab('Trial no.') + geom_hline(yintercept = rates_e2[1], linetype="dashed", color="red") +
  geom_hline(yintercept = rates_e2[2], linetype="dashed", color="blue")  + theme(legend.position="bottom")

e3_sample <- dplyr::filter(data_pred_e3, sub==subs_e3[sample_sub_e3])[11:110,] %>% mutate(s0=s0 * ((targets>0) * 2 - 1)) %>%
  mutate(tar_start = tno, tar_end = tno+1)

p_e3_s0 <- e3_sample %>% ggplot(aes(x=tno, y=s0)) + geom_line() + theme_bw() +
  geom_rect( aes(xmin = tar_start, xmax = tar_end,x =NULL, y = NULL, ymin = -0.3, ymax = 0.3, fill = target), alpha=0.5) + 
  xlab('Trial no.') + theme(legend.position="bottom")

rates_e3 <- c(dplyr::filter(exp3_pars, sub==sample_sub_e3, model==top_model_e3)$mu1, 
              dplyr::filter(exp3_pars, sub==sample_sub_e3, model==top_model_e3)$mu2, 
              dplyr::filter(exp3_pars, sub==sample_sub_e3, model==top_model_e3)$mu3) 

p_e3_rate <- e3_sample %>% ggplot(aes(x=tno, y=rate)) + geom_line() + theme_bw()  +
  geom_rect( aes(xmin = tar_start, xmax = tar_end,x =NULL, y = NULL, ymin = 2, ymax = 9, fill = dimension), alpha=0.5) + 
  xlab('Trial no.') + geom_hline(yintercept = rates_e3[1], linetype="dashed", color="red") +
  geom_hline(yintercept = rates_e3[2], linetype="dashed", color="darkgreen") +
  geom_hline(yintercept = rates_e3[3], linetype="dashed", color="blue") + theme(legend.position="bottom")

updating_example_plot <- plot_grid(p_e1_s0, p_e1_rate, p_e2_s0, p_e2_rate, p_e3_s0, p_e3_rate, ncol=2, nrow=3,
                                   labels=c("A", "B", "C", "D", "E", "F"))

set.seed(92)
t <- seq(0,0.53,0.01)
y_ddm <- c(0, cumsum(rnorm(53))+20*t[1:53])
y_l <- seq(0, 10.32, length.out=length(t)) 
x_poly <- c(0, max(t)-0.1, max(t)+0.17)
y_poly <- c(0, max(y_l), max(y_l))
d <- data.frame(t=t, y_ddm=y_ddm, y_l=y_l)
poly <- data.frame(x=x_poly, y=y_poly)
DDM_plot <- ggplot(d, aes(x=t, y=y_ddm)) + geom_line(color="blue") + geom_line(color="red", aes(x=t, y=y_l)) + 
  geom_polygon(data=poly, aes(x=x, y=y), linetype=0, fill="red", alpha=0.15) + theme_void() + geom_hline(yintercept=10.32) + geom_hline(yintercept=-10.32) + 
  geom_hline(yintercept=0, linetype=2) + coord_cartesian(xlim=c(0,0.8))
