##################################################
#### Bayesian differential abundance analyses ####
##################################################


#####################################
#### Set parameters for analyses ####
#####################################


#### Taxonomy ####

level <- "Genus"

#### Filter ####

#### Minimum overall prevalence - ASV prevalent in more than X percent of all samples
prev <- .2



###############################
#### Load data and scripts ####
###############################

#### load packages ####
source('packages.R')

library(tidybayes)
library(ggridges)
library(glue)
library(hrbrthemes)
library(ggdist)
library(ggsci)
library(brms)
library(parallel)



#### load functions ####
source("functions.R")

#### Load and prepare data ####

dat <- readRDS("Data/processed_data/phyl_list.RDS")


#### Select datasets with data available
dat <- dat[c(3,4,6,7)]



#### Collect vairables names
HY_stage <- c("HY.stage..","hys","Hoehn.Yahr.stage")
duration <- c("Disease.duration..yrs.","Disease.duration..months.","duration")
smoking <- c("smoking","Smoker")
constipation <- c("Constipation","constipation")
comt <- c("meds_COMT_inhibitor","comt")
ldopa <- c("ldopa","meds_dopa","Average.L.dopa.dose..last.2.years.")


meta <- c("status","author","BMI","gender","age",
          "meda_MAO_inhibitor","meda_dopamine_agonist",
          "meds_statin","meda_Warfarin","meda_ca_antagonist","Wexner_total","MMSE_total",
          "UPDRS_I_total","UPDRS_II_total","UPDRS_III_total_ON","UPDRS_III_total_OFF","UPDRS_IV_total",
          "UPDRS_V_ON","UPDRS_V_OFF","hypert",
          "diabetes","led","nmse","nmsq","updrs.1", "updrs.2","updrs.3",
          "updrs.4") %>% c(HY_stage,duration,smoking,constipation,comt,ldopa)


full.dat <- prepare_DA_data(dat,tax=level,meta=meta,preval = prev,relab=T)


### Wrangle data

full.dat %>%
  # Disease duration
  mutate(disease_duration = Disease.duration..months.,.keep="unused") %>%
  mutate(disease_duration = ifelse(is.na(disease_duration), 12*Disease.duration..yrs.,disease_duration),.keep="unused") %>%
  mutate(disease_duration = ifelse(is.na(disease_duration), 12*duration,disease_duration),.keep="unused") %>%
  # Smoking
  mutate(smoking = ifelse(smoking == 1,"yes","no")) %>%
  mutate(Smoker = ifelse(is.na(Smoker),smoking,Smoker),.keep = "unused") %>%
  # HY stage
  convert(num(all_of(HY_stage))) %>%
  mutate(hys = dplyr::coalesce(HY.stage..,hys,Hoehn.Yahr.stage),.keep = "unused") %>%
  # Constipation
  mutate(constipation = ifelse(constipation==1,"yes","no")) %>%
  mutate(Constipation = ifelse(is.na(Constipation),constipation,Constipation),.keep = "unused") %>%
  ## COMT inhibitor
  mutate(comt = coalesce(meds_COMT_inhibitor,comt),.keep = "unused") %>%
  #ldopa
  mutate(ldopa_bin = ifelse(ldopa > 0,1,0)) %>%
  mutate(avg_ldopa_bin = ifelse(Average.L.dopa.dose..last.2.years. > 0,1,0)) %>%
  mutate(ldopa_bin = coalesce(ldopa_bin,avg_ldopa_bin,meds_dopa),.keep="unused") -> full.dat



#### If necessary, select only specific taxa

taxa <- c("Synergistetes","Verrucomicrobia","Proteobacteria","Lactobacillaceae","Synergistaceae","Verrucomicrobiaceae","Peptococcaceae.2","Clostridiales_Incertae.Sedis.XII",
          "Bifidobacteriaceae","Porphyromonadaceae","Enterobacteriaceae","Lachnospiraceae","Ruminococcaceae","Akkermansia","Lactobacillus","Porphyromonas","Eisenbergiella","Desulfurispora","Acidaminobacter","Bifidobacterium",
          "Anaerotruncus","Blautia","Prevotella","Faecalibacterium","Butyricicoccus","Fusicatenibacter","Roseburia","Catabacter",
          "Hungatella")

phyla <- c("Synergistetes","Verrucomicrobia","Proteobacteria")


families <- c("Lactobacillaceae","Synergistaceae","Verrucomicrobiaceae","Peptococcaceae.2","Clostridiales_Incertae.Sedis.XII",
          "Bifidobacteriaceae","Porphyromonadaceae","Enterobacteriaceae","Lachnospiraceae","Ruminococcaceae")


genera <- c("Akkermansia","Lactobacillus","Eisenbergiella","Desulfurispora","Acidaminobacter","Bifidobacterium",
          "Anaerotruncus","Blautia","Prevotella","Faecalibacterium","Butyricicoccus","Fusicatenibacter","Roseburia","Catabacter",
          "Hungatella")

##################
#### HY stage ####
##################

full.dat %>%
  filter(!is.na(hys)) %>%
  pivot_longer(cols = all_of(genera),names_to = "taxa") %>%
  ggplot(aes(x=as.factor(hys),y = value)) +
  geom_boxplot() +
  facet_grid(author~taxa,scales = "free_y")

hy <- c("Akkermansia","Lactobacillus","Bifidobacterium","Faecalibacterium","Blautia","Roseburia")


full.dat %>%
  filter(!is.na(hys)) %>%
  mutate(hys = as.factor(hys)) %>%
  filter(author == "Barichella et al., 2019") %>%
  pivot_longer(cols = all_of(hy),names_to = "taxa") %>%
  dplyr::select(author,value,taxa,hys) %>%
  ggplot(aes(x=hys,y = value)) +
  geom_boxplot(fill="#DF8F44FF",alpha=.8,outlier.shape = NA) + 
  geom_jitter(position=position_jitter(0.2),alpha=.4) +
  theme_minimal() +
  facet_wrap(~taxa,scales = "free_y",nrow = 1) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Relative abundance") -> p1

full.dat %>%
  filter(!is.na(hys)) %>%
  mutate(hys = as.factor(hys)) %>%
  filter(author == "Keshavarzian et al., 2015") %>%
  pivot_longer(cols = all_of(hy),names_to = "taxa") %>%
  dplyr::select(author,value,taxa,hys) %>%
  ggplot(aes(x=hys,y = value)) +
  geom_boxplot(fill="#DF8F44FF") +
  geom_jitter(position=position_jitter(0.2),alpha=.4) +
  theme_minimal() +
  facet_wrap(~taxa,scales = "free_y",nrow = 1) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank()) +
  ylab("Relative abundance") +
  scale_x_discrete(breaks=factor(c(1,1.5,2,2.5,3,4,5)), drop=FALSE) -> p2

full.dat %>%
  filter(!is.na(hys)) %>%
  mutate(hys = as.factor(hys)) %>%
  filter(author == "Weis et al., 2019") %>%
  pivot_longer(cols = all_of(hy),names_to = "taxa") %>%
  dplyr::select(author,value,taxa,hys) %>%
  ggplot(aes(x=as.factor(hys),y = value)) +
  geom_boxplot(fill="#DF8F44FF") +
  geom_jitter(position=position_jitter(0.2),alpha=.4) +
  theme_minimal() +
  facet_wrap(~taxa,scales = "free_y",nrow = 1) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.x = element_text("Hoehn and Yahr Stage"),
    axis.text.x = element_text(angle=45)) +
  ylab("Relative abundance") +
  xlab("Hoehn and Yahr Stage") +
  scale_x_discrete(breaks=factor(c(1.0,1.5,2.0,2.5,3.0,4.0,5.0)), drop=FALSE) -> p3

p_hy <- ggarrange(p1 +scale_y_log10(breaks=scales::breaks_pretty()),
                  p2 +scale_y_log10(breaks=scales::breaks_pretty()),
                  p3 +scale_y_log10(breaks=scales::breaks_pretty()),ncol=1,labels = c("A","B","C"))


ggsave(p_hy,filename = "Plots/hy_stage.PDF",device = cairo_pdf,width = 10,height = 8.3,units = "in",dpi = 600)


##########################
#### Disease duration ####
##########################

genera <- c("Akkermansia","Lactobacillus","Bifidobacterium",
            "Blautia","Prevotella","Faecalibacterium","Roseburia")

genera <- c("Ruminococcus")

full.dat %>%
  filter(!is.na(disease_duration)) %>%
  filter(author == "Barichella et al., 2019") %>%
  pivot_longer(cols = all_of(genera),names_to = "taxa") %>%
  filter(value > 0) %>%
  #mutate(disease_duration = (as.numeric(disease_duration)/12) ) %>%
  ggplot(aes(x=as.numeric(disease_duration),y = value)) +
  geom_point(colour="#374E55FF") +
  geom_smooth(method="lm",color="#DF8F44FF") +
  facet_wrap(~taxa,scales = "free",nrow=1) +
  theme_minimal() +
  scale_y_log10(breaks=scales::pretty_breaks()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle=45)) +
  ylab("Relative abundance") -> p1


full.dat %>%
  filter(!is.na(disease_duration)) %>%
  filter(author == "Keshavarzian et al., 2015") %>%
  pivot_longer(cols = all_of(genera),names_to = "taxa") %>%
  filter(value > 0) %>%
  #mutate(disease_duration = (as.numeric(disease_duration)/12) ) %>%
  ggplot(aes(x=as.numeric(disease_duration),y = value)) +
  geom_point(colour="#374E55FF") +
  geom_smooth(method="lm",color="#DF8F44FF") +
  facet_wrap(~taxa,scales = "free",nrow=1) +
  theme_minimal() +
  scale_y_log10(breaks=scales::pretty_breaks()) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45)) +
  ylab("Relative abundance") -> p2


full.dat %>%
  filter(!is.na(disease_duration)) %>%
  filter(author == "Weis et al., 2019") %>%
  pivot_longer(cols = all_of(genera),names_to = "taxa") %>%
  filter(value > 0) %>%
  #mutate(disease_duration = (as.numeric(disease_duration)/12) ) %>%
  ggplot(aes(x=as.numeric(disease_duration),y = value)) +
  geom_point(colour="#374E55FF") +
  geom_smooth(method="lm",color="#DF8F44FF") +
  facet_wrap(~taxa,scales = "free",nrow=1) +
  theme_minimal() +
  scale_y_log10(breaks=scales::pretty_breaks()) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle=45)) +
  ylab("Relative abundance") +
  xlab("Disease duration in months")-> p3


p_dur <- ggarrange(p1,p2,p3,ncol=1,labels = c("A","B","C"))


ggsave(p_dur,filename = "Plots/disease_dur.PDF",device = cairo_pdf,width = 10,height = 10,units = "in",dpi = 600)


####################
#### medication ####
####################

genera <- c("Lactobacillus")

full.dat %>%
  filter(!is.na(ldopa_bin)) %>%
  pivot_longer(cols = all_of(genera),names_to = "taxa") %>%
  ggplot(aes(x=as.factor(ldopa_bin),y = value)) +
  geom_boxplot() +
  facet_wrap(author~taxa,scales = "free")

full.dat %>%
  filter(!is.na(comt)) %>%
  pivot_longer(cols = all_of(genera),names_to = "taxa") %>%
  ggplot(aes(x=as.factor(comt),y = value)) +
  geom_boxplot() +
  facet_grid(author~taxa,scales = "free")

full.dat %>%
  filter(!is.na(ldopa_bin)) %>%
  dplyr::select(author,ldopa_bin,Lactobacillus) %>%
  group_by(author) %>%
  summarize(co = cor(ldopa_bin,Lactobacillus,method = "spearman"))
  

full.dat %>%
  filter(!is.na(comt)) %>%
  dplyr::select(author,comt,Lactobacillus) %>%
  group_by(author) %>%
  summarize(co = cor(comt,Lactobacillus,method = "spearman"))


### COMT AND LDOPA

full.dat %>%
  filter(!is.na(ldopa_bin)) %>%
  mutate(ldopa_bin = ifelse(ldopa_bin == 1,"Yes","No")) %>%
  pivot_longer(cols = all_of(genera),names_to = "taxa") %>%
  ggplot(aes(x=as.factor(ldopa_bin),y = value,fill = as.factor(ldopa_bin))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = .6) +
  facet_grid(taxa~author,scales = "free") +
  scale_y_log10(breaks = scales::pretty_breaks()) +
  scale_fill_jama() +
  ylab("Relative abundance") +
  xlab("Treatment with LDOPA") +
  theme(legend.position = "none") -> p1


full.dat %>%
  filter(!is.na(comt)) %>%
  pivot_longer(cols = all_of(genera),names_to = "taxa") %>%
  mutate(comt = ifelse(comt == 1,"Yes","No")) %>%
  ggplot(aes(x=as.factor(comt),y = value,fill = as.factor(comt))) +
  facet_grid(taxa~author,scales = "free") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = .6) +
  scale_y_log10(breaks = scales::pretty_breaks()) +
  scale_fill_jama() +
  xlab("Treatment with COMT inhibitors") +
  ylab("Relative abundance") +
  theme(legend.position = "none") -> p2


p_meds <- ggarrange(p1,p2,ncol=2,labels = c("A","B"))


ggsave(p_meds,filename = "Plots/meds.PDF",device = cairo_pdf,width = 10,height = 6,units = "in",dpi = 600)