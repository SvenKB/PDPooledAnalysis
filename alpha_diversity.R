######################################################################
#### Microbiome Meta-Analysis - Bayesian alpha diversity analyses ####
######################################################################

#### load packages ####
source('packages.R')

#### load functions ####
source("functions.R")

library(brms)
library(tidybayes)
library(broom)
library(ggdist)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(glue)


##########################
#### Data preparation ####
##########################

dat <- readRDS("Data/processed_data/phyl_list.RDS")

alpha_dat <- estimate_alpha(dat)

author.names <- c("Pietrucci et al., 2019",
                  "Qian et al., 2018",
                  "Keshavarzian et al., 2015",
                  "Aho et al., 2019",
                  "Heintz-Buschart et al., 2018",
                  "Weis et al., 2019",
                  "Barichella et al., 2019",
                  "Wallen et al., 2020a",
                  "Wallen et al., 2020b")

names <- cbind(study=paste0("Study",1:9),author.names)


###########################
#### Inspect diversity ####
###########################

#### Same region ####

ggplot(aes(x=as.factor(author),y=(richness)),data=alpha_dat) +
  geom_boxplot(aes(fill=status)) +
  scale_fill_jama() +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Richness")

ggplot(aes(x=as.factor(author),y=(shannon)),data=alpha_dat) +
  geom_boxplot(aes(fill=status)) +
  scale_fill_jama() +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Shannon")

ggplot(aes(x=as.factor(author),y=inv.simpson),data=alpha_dat) +
  geom_boxplot(aes(fill=status)) +
  scale_fill_jama() +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Inv. Simpson")


#### Combined plots ####
plot_alpha <- alpha_dat %>%
  pivot_longer(cols=c("richness","shannon","inv.simpson"),names_to="measure") %>%
  mutate(measure = factor(measure,levels=c("richness","shannon","inv.simpson")))


ggplot(aes(x=author,y=value,fill=status),data=plot_alpha) +
  geom_jitter(aes(group=status,color=status),position = position_jitterdodge(),alpha=.2) +
  geom_boxplot() +
  coord_flip() +
  facet_wrap(~measure,nrow=1,scales="free_x") +
  theme_pubclean() +
  ylab("Effective number of taxa") +
  scale_fill_jama() +
  scale_color_jama() +
  theme(axis.text.x = element_text(angle=90))


#### Correlation of measures with seq depth ####
alpha_dat %>%
  group_by(author) %>%
  dplyr::summarize("Richness" = round(cor(richness,N),2),
                   "Shannon" = round(cor(shannon,N),2),
                   "Simpson" = round(cor(inv.simpson,N),2)) %>%
  data.frame(row.names=NULL) %>%
  column_to_rownames("author") %>%
  t %>% kable(caption = "Correlation of alpha diversity indices with library size - Same region sequences") %>%
  kableExtra::kable_styling("striped")


ggplot(aes(x=N,y=value),data=plot_alpha) +
  geom_point(alpha=.5) +
  geom_smooth(colour="red") +
  facet_grid(author~measure,scales="free") +
  theme_hc()


######################################
#### Bayesian random effect model ####
######################################

#### Same region sequences ####

alpha_dat <- alpha_dat %>%
  mutate(author = case_when(author==names[[1,2]] ~ names[[1,1]],
                            author==names[[2,2]] ~ names[[2,1]],
                            author==names[[3,2]] ~ names[[3,1]],
                            author==names[[4,2]] ~ names[[4,1]],
                            author==names[[5,2]] ~ names[[5,1]],
                            author==names[[6,2]] ~ names[[6,1]],
                            author==names[[7,2]] ~ names[[7,1]],
                            author==names[[8,2]] ~ names[[8,1]],
                            author==names[[9,2]] ~ names[[9,1]]))


# Richness
fit_SR_RF_R <- brms::brm(richness~1+status+scale(N)+(1+status|author),data=alpha_dat,family=negbinomial,cores=4)

study.draws <- fit_SR_RF_R %>% spread_draws(b_Intercept,
                                            b_statusPD,
                                            r_author[author,effect]) %>%
  filter(effect=="statusPD") %>%
  mutate(b_statusPD = r_author + b_statusPD) %>%
  mutate(b_statusPD = exp(b_statusPD)) %>%
  ungroup %>%
  dplyr::select(.chain,.iteration,.draw,b_statusPD,author) %>% convert(fct(author))

pooled.draws <- fit_SR_RF_R %>% spread_draws(b_statusPD) %>%
  mutate(b_statusPD = exp(b_statusPD),
         author = "Pooled Effect")


forest.data.richness <- bind_rows(study.draws, pooled.draws) %>% 
  ungroup() %>%
  mutate(author = case_when(author==names[[1,1]] ~ names[[1,2]],
                            author==names[[2,1]] ~ names[[2,2]],
                            author==names[[3,1]] ~ names[[3,2]],
                            author==names[[4,1]] ~ names[[4,2]],
                            author==names[[5,1]] ~ names[[5,2]],
                            author==names[[6,1]] ~ names[[6,2]],
                            author==names[[7,1]] ~ names[[7,2]],
                            author==names[[8,1]] ~ names[[8,2]],
                            author==names[[9,1]] ~ names[[9,2]],
                            author=="Pooled Effect" ~ "Pooled Effect")) %>%
  mutate(author = reorder(author,b_statusPD)) %>%
  mutate(measure="Richness")

# Shannon
fit_SR_RF_SH <- brms::brm(shannon~1+status+scale(N)+(1+status|author),data=alpha_dat,family=negbinomial,cores=4)

study.draws <- fit_SR_RF_SH %>% spread_draws(b_Intercept,
                                             b_statusPD,
                                             r_author[author,effect]) %>%
  filter(effect=="statusPD") %>%
  mutate(b_statusPD = r_author + b_statusPD) %>%
  mutate(b_statusPD = exp(b_statusPD)) %>%
  ungroup %>%
  dplyr::select(.chain,.iteration,.draw,b_statusPD,author) %>% convert(fct(author))

pooled.draws <- fit_SR_RF_SH %>% spread_draws(b_statusPD) %>%
  mutate(b_statusPD = exp(b_statusPD),
         author = "Pooled Effect")


forest.data.shannon <- bind_rows(study.draws, pooled.draws) %>% 
  ungroup() %>%
  mutate(author = case_when(author==names[[1,1]] ~ names[[1,2]],
                            author==names[[2,1]] ~ names[[2,2]],
                            author==names[[3,1]] ~ names[[3,2]],
                            author==names[[4,1]] ~ names[[4,2]],
                            author==names[[5,1]] ~ names[[5,2]],
                            author==names[[6,1]] ~ names[[6,2]],
                            author==names[[7,1]] ~ names[[7,2]],
                            author==names[[8,1]] ~ names[[8,2]],
                            author==names[[9,1]] ~ names[[9,2]],
                            author=="Pooled Effect" ~ "Pooled Effect")) %>%
  mutate(author = reorder(author,b_statusPD)) %>%
  mutate(measure="Shannon")

# Shannon
fit_SR_RF_IS <- brms::brm(inv.simpson~1+status+scale(N)+(1+status|author),data=alpha_dat,family=negbinomial,cores=4)

study.draws <- fit_SR_RF_IS %>% spread_draws(b_Intercept,
                                             b_statusPD,
                                             r_author[author,effect]) %>%
  filter(effect=="statusPD") %>%
  mutate(b_statusPD = r_author + b_statusPD) %>%
  mutate(b_statusPD = exp(b_statusPD)) %>%
  ungroup %>%
  dplyr::select(.chain,.iteration,.draw,b_statusPD,author) %>% convert(fct(author))

pooled.draws <- fit_SR_RF_IS %>% spread_draws(b_statusPD) %>%
  mutate(b_statusPD = exp(b_statusPD),
         author = "Pooled Effect")


forest.data.simpson <- bind_rows(study.draws, pooled.draws) %>% 
  ungroup() %>%
  mutate(author = case_when(author==names[[1,1]] ~ names[[1,2]],
                            author==names[[2,1]] ~ names[[2,2]],
                            author==names[[3,1]] ~ names[[3,2]],
                            author==names[[4,1]] ~ names[[4,2]],
                            author==names[[5,1]] ~ names[[5,2]],
                            author==names[[6,1]] ~ names[[6,2]],
                            author==names[[7,1]] ~ names[[7,2]],
                            author==names[[8,1]] ~ names[[8,2]],
                            author==names[[9,1]] ~ names[[9,2]],
                            author=="Pooled Effect" ~ "Pooled Effect")) %>%
  mutate(author = reorder(author,b_statusPD)) %>%
  mutate(measure="Inv.Simpson")

res_all_pooled <- bind_rows(forest.data.richness,
                            forest.data.shannon,
                            forest.data.simpson) %>% mutate(measure = factor(measure,levels=c("Richness","Shannon","Inv.Simpson")))

saveRDS(res_all_pooled,"results_pooled_random_effect.RDS")

res_all_pooled <- readRDS("results_pooled_random_effect.RDS")

res_all_pooled %>%
  mutate(measure = recode_factor(measure,"Inv.Simpson" = "Simpson")) %>%
  mutate(measure = factor(measure,levels=c("Richness","Shannon","Simpson"))) -> res_all_pooled


ggplot(aes(x=b_statusPD,relevel(author,"Pooled Effect")),data=res_all_pooled) +
  geom_vline(xintercept = 1, color = "darkgrey", size = 1) +
  stat_gradientinterval(position ="dodge") +
  #stat_halfeye(position=position_dodge(width=1),alpha=1,color="black") +
  labs(x = "Rate Ratio",
       y = element_blank()) +
  theme_minimal() +
  scale_x_log10() +
  scale_fill_jama() +
  #scale_fill_manual(values=pal) +
  facet_wrap(~measure,ncol=1,strip.position = "right",scales="free") +
  theme(strip.text = element_text(size=16))

summary_data <- res_all_pooled %>%
  mutate(b_statusPD = log(b_statusPD)) %>%
  group_by(author,measure) %>%
  mean_qi(b_statusPD,.width = .95) %>%
  mutate(b_statusPD = round(exp(b_statusPD),2),
         .lower = round(exp(.lower),2),
         .upper = round(exp(.upper),2)) %>%
  mutate(label = glue("{format(b_statusPD,digits=3)} [{format(.lower,digits=3)}, {format(.upper,digits=3)}]")) 


levels(res_all_pooled$author)

res_all_pooled %>% mutate(author = factor(author,levels = levels(author)[c(6,5,9,7,10,1,8,3,4,2)])) -> res_all_pooled

data.frame(author=NA)

p <- ggplot(aes(x=b_statusPD,author),data=res_all_pooled) +
  geom_vline(xintercept = 1, color = "darkgrey", size = 1) +
  stat_gradientinterval(position ="dodge",aes(group=author,slab_alpha=stat(-pmax(abs(1-2*cdf),.95)),fill = stat(x > 0))) +
  geom_vline(aes(xintercept=b_statusPD),color="black",linetype=1,data= summary_data %>% filter(author=="Pooled Effect")) +
  geom_vline(aes(xintercept=.lower),color="black",linetype=2,data= summary_data %>% filter(author=="Pooled Effect")) +
  geom_vline(aes(xintercept=.upper),color="black",linetype=2,data= summary_data %>% filter(author=="Pooled Effect")) +
  #stat_halfeye(aes(fill = stat(exp(x)>1)),position=position_dodge(width=1),alpha=1,color="black") +
  #stat_halfeye(aes(fill = stat(cut_cdf_qi(cdf,.width=c(.95,1)))),position=position_dodge(width=1),alpha=1,color="black") +
  labs(x = "Rate Ratio",
       y = element_blank()) +
  #theme_pubr() +
  theme_minimal() +
  scale_x_log10(limits=c(0.5,4),breaks=c(0.5,0.75,1,2)) +
  scale_fill_jama() +
  geom_text(data = summary_data,
            aes(label = label, x = 2.2),hjust=0,colour="black") +
  #scale_fill_manual(values=pal) +
  facet_wrap(~measure,ncol=1,strip.position = "right",scales="free_x") +
  theme(strip.text = element_text(size=16),
        legend.position = "none") 


ggsave(p,filename = "Plots/alpha_diversity_intervalslab.PDF",device = "pdf",width = 20,height = 8,units = "in",dpi = 600)


res_all_pooled$author <-  factor(res_all_pooled$author, levels=c(levels(res_all_pooled$author)[1]," ",levels(res_all_pooled$author)[-1]))

res_all_pooled$shaper <- as.factor(ifelse(res_all_pooled$author=="Pooled Effect",1,0))
res_all_pooled$sizer <- ifelse(res_all_pooled$author=="Pooled Effect",16,10)
  
p <- ggplot(aes(x=b_statusPD,author),data=res_all_pooled) +
 # geom_vline(xintercept = 1, color = "darkgrey", size = 1) +
  #stat_gradientinterval(position ="dodge",aes(group=author,slab_alpha=stat(-pmax(abs(1-2*cdf),.95)),fill = stat(x > 0))) +
  #stat_gradientinterval(position ="dodge",aes(group=author)) +
  stat_halfeye(aes(fill = stat(exp(x)>1)),position=position_dodge(width=1),alpha=1,color="black",scale=1.4) +
  #stat_density_ridges(aes(fill = stat(exp(x)>1)),
  #                    calc_ecdf = TRUE,geom="density_ridges_gradient",
  #                    quantiles = c(0.05, 0.95),
  #                    scale=.9) +
  #stat_pointinterval(.width = c(0.66, 0.95,.99)) +
  geom_vline(aes(xintercept=b_statusPD),color="black",linetype=1,data= summary_data %>% filter(author=="Pooled Effect")) +
  geom_vline(aes(xintercept=.lower),color="black",linetype=2,data= summary_data %>% filter(author=="Pooled Effect")) +
  geom_vline(aes(xintercept=.upper),color="black",linetype=2,data= summary_data %>% filter(author=="Pooled Effect")) +
  labs(x = "Rate Ratio",
       y = element_blank()) +
  theme_minimal() +
  #theme_minimal() +
  scale_x_log10(limits=c(0.68,4),breaks=c(0.75,1,1.5)) +
  scale_fill_jama(drop=FALSE) +
  scale_y_discrete(drop=F) +
  geom_text(data = summary_data,
            aes(label = label, x = 2),hjust=0,colour="black",size=7) +
  #scale_fill_manual(values=pal) +
  facet_wrap(~measure,nrow=1,strip.position = "top",scales="free_x") +
  theme(strip.text = element_text(size=24),
        legend.position = "none",
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(hjust = 0, size=20,face = "bold"),
        axis.title.x = element_text(size=24))

ggsave(p,filename = "Plots/alpha_diversity_intervalslab1.PDF",device = "pdf",width = 20,height = 10,units = "in",dpi = 600)

