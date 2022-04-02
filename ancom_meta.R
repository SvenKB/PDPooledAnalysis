#############################
#### ANCOM Meta-Analysis ####
#############################


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
library(ANCOMBC)
library(metafor)


#### load functions ####
source("functions.R")

phyl <- readRDS("Data/processed_data/phyl_list.RDS")
names <- unique(unlist(lapply(phyl,function(x) sample_data(x)[,"author"])))
names(phyl) <- names


fitANCOMBC <- function(x) {
  fit <- ancombc(phyl[[x]],
                 formula = "status",
                 p_adj_method = "BH",
                 zero_cut = .9,
                 group = "status",struc_zero = T,
                 neg_lb = T)
  res <- data.frame("author" = x,
                    "Taxon" = rownames(fit$res$beta),
                    "Y" = fit$res$beta[[1]],
                    "V" = fit$res$se[[1]],
                    "Sig" = fit$res$diff_abn[[1]])
  
  return(res)
  
}

## Genus
phyl <- lapply(phyl,function(x) aggregate_taxa(x,level="Genus"))
res_genus <- reduce(lapply(names,fitANCOMBC),bind_rows)

## Family
phyl <- lapply(phyl,function(x) aggregate_taxa(x,level="Family"))
res_family <- reduce(lapply(names,fitANCOMBC),bind_rows)

## Phylum
phyl <- lapply(phyl,function(x) aggregate_taxa(x,level="Phylum"))
res_phylum <- reduce(lapply(names,fitANCOMBC),bind_rows)


## Combine study-wise ANCOM results
res <- bind_rows(res_phylum,
                 res_family,
                 res_genus)


## Apply random-effect meta-analysis to all taxa
res %>%
  group_by(Taxon) %>%
  nest() %>%
  mutate(re_model = map(data, function(x) rma(yi = Y, vi = V, data=x,method = "REML"))) -> results



## Summarize models in a tidy way
tidyParameters <- function(x) {
  out <- tibble(coef(summary(x)),
                "Studies"=x$k,
                "I2" = x$I2,
                "tau2" = x$tau2,
                )
  return(out)
}

results %>%
  mutate(tidied = map(re_model, tidyParameters)) %>%
  unnest(cols = c(Taxon,tidied)) %>%
  dplyr::select(-c(data,re_model)) %>%
  ungroup %>%
  mutate(qval = p.adjust(pval,"BH"),.after="pval") -> res_summary


## Prepare forst plots

library(ggforestplot)


#### Look at primary results

taxa <- c("Synergistetes","Verrucomicrobia","Proteobacteria","Lactobacillaceae","Synergistaceae","Verrucomicrobiaceae","Peptococcaceae.2","Clostridiales_Incertae.Sedis.XII",
          "Bifidobacteriaceae","Porphyromonadaceae","Enterobacteriaceae","Lachnospiraceae","Ruminococcaceae","Akkermansia","Lactobacillus","Porphyromonas","Eisenbergiella","Desulfurispora","Acidaminobacter","Bifidobacterium",
          "Anaerotruncus","Blautia","Prevotella","Faecalibacterium","Butyricicoccus","Fusicatenibacter","Roseburia","Catabacter",
          "Hungatella")

phyla <- c("Synergistetes","Verrucomicrobia","Proteobacteria","Actinobacteria")


families <- c("Lactobacillaceae","Synergistaceae","Verrucomicrobiaceae","Peptococcaceae.2","Clostridiales_Incertae.Sedis.XII",
              "Bifidobacteriaceae","Porphyromonadaceae","Enterobacteriaceae","Lachnospiraceae","Ruminococcaceae")


genera <- c("Akkermansia","Lactobacillus","Porphyromonas","Eisenbergiella","Desulfurispora","Acidaminobacter","Bifidobacterium",
            "Anaerotruncus","Blautia","Prevotella","Faecalibacterium","Butyricicoccus","Fusicatenibacter","Roseburia","Catabacter",
            "Hungatella")

### Function for single taxa forest plots
plotForest <- function(x,taxon="Proteobacteria") {
  x %>%
    filter(Taxon==taxon) %>%
    pull(re_model) -> dat
  dat <- dat[[1]]
  forest(dat,annotate = T)
  text(-0, -1.5, pos=2, cex=0.75, bquote(paste("(Q = ",
                                              .(formatC(dat$QE, digits=2, format="f")), ", df = ", .(dat$k - dat$p),
                                              ", p = ", .(formatC(dat$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                              .(formatC(dat$I2, digits=1, format="f")), "%, ", tau^2, " = ",
                                              .(formatC(dat$tau2, digits=1, format="f")),")")))
}

plotForest(results,"Akkermansia")



### Summarized forest plot

#### ANCOM 2-stage results

taxa <- c(taxa,"Actinobacteria")


res_summary %>%
  filter(Taxon %in% taxa) %>%
  mutate(across(is.numeric,.fns=function(x) round(x,4))) %>%
  mutate(Rank = case_when(Taxon %in% phyla ~ "Phylum",
                          Taxon %in% families ~ "Family",
                          Taxon %in% genera ~ "Genus"),
         Rank = factor(Rank,levels = c("Phylum","Family","Genus")),
         Sig = ifelse(qval < 0.05, TRUE,FALSE),
         Label = glue("{round(estimate,2)} ({round(ci.lb,2)},{round(ci.ub,2)}) ;"),
         I_label = glue("I2 = {format(round(I2,2),digits=3)}%, {Studies} Studies")) -> prim_res


forestplot(
  df = prim_res,
  name = Taxon,
  estimate = estimate,
  se = se,
  pval = qval,
  colour = Sig,
  fatten=4,
  xlab = "Beta"
) + ggforce::facet_col(
  facets = ~Rank,
  space = "free",
  scales = "free_y"
) +
  geom_text(aes(x= 2.0,label=Label),hjust=1,size=3) +
  geom_text(aes(x= 2.01,label=I_label),hjust=0,size=3) +
  coord_cartesian(xlim = c(-1.5, 1.3), # This focuses the x-axis on the range of interest
                  clip = 'off') + theme(legend.position = "none") +
  theme(plot.margin = unit(c(.5,12,.5,.5), "lines")) -> p
  

ggsave(p,filename = "Plots/ANCOM_meta.PDF",device = "pdf",width = 10.5,height = 6,units = "in",dpi = 600)


###

taxa <- c("Roseburia",
          "Fusicatenibacter",
          "Eisenbergiella",
          "Blautia",
          "Faecalibacterium",
          "Butyricicoccus",
          "Anaerotruncus",
          "Desulfurispora",
          "Acidaminobacter",
          "Lactobacillus",
          "Prevotella",
          "Porphyromonas",
          "Akkermansia",
          "Synergistaceae",
          "Enterobacteriaceae",
          "Bifidobacterium"
          )



forestplot(
  df = prim_res %>% filter(Taxon %in% taxa) %>% mutate(Taxon = ordered(Taxon, levels=taxa)) %>% arrange(Taxon),
  name = Taxon,
  estimate = estimate,
  se = se,
  pval = qval,
  colour = Sig,
  fatten=4,
)  + geom_vline(aes(xintercept=0),size=1) + theme(legend.position = "none",
           axis.text.y = element_blank(),
           axis.title.x = element_blank(),
           plot.margin =  unit(c(1.8,1,1,0), "lines")) +
  scale_color_jama() -> p_ancom


  
ggsave(p,filename = "Plots/ANCOM_effects.PDF",device = "pdf",width = 3,height = 6,units = "in",dpi = 600)


### Pooled analysis

bayes_dat <- readRDS("Data/other/summary_bayes_pooled.RDS")


bayes_dat %>%
  mutate("Taxon" = Outcome,
         "estimate" = b_statusPD,
         "se" = (.upper-.lower) / 5.15,
         .keep="unused",
         Rank = case_when(Taxon %in% phyla ~ "Phylum",
                          Taxon %in% families ~ "Family",
                          Taxon %in% genera ~ "Genus"),
         Rank = factor(Rank,levels = c("Phylum","Family","Genus"))) %>% dplyr::select(-c(author,.width,.point,.interval)) -> bayes_dat


forestplot(
  df = bayes_dat,
  name = Taxon,
  estimate = estimate,
  se = se,
  #colour = Sig,
  title = "2-stage meta analysis based on ANCOM-BC",
  xlab = "Beta"
) + ggforce::facet_col(
  facets = ~Rank,
  space = "free",
  scales = "free_y"
) +
  scale_x_log10() +
  geom_text(aes(x= 2.7,label=label),hjust=1,size=3) +
  #coord_cartesian(xlim = c(-1.5, 1.5), # This focuses the x-axis on the range of interest
  #                clip = 'off') + theme(legend.position = "none") +
  theme(plot.margin = unit(c(.5,12,.5,.5), "lines"))

res %>%
  mutate(estimate = Y,
         se = V,.keep = "unused",
         study = author) %>%
  dplyr::select(-Sig) %>%
  filter(Taxon %in% taxa)-> study_res

prim_res %>%
  mutate(study = "Pooled") %>%
  bind_rows(study_res) %>%
  mutate(Rank = case_when(Taxon %in% phyla ~ "Phylum",
                                                 Taxon %in% families ~ "Family",
                                                 Taxon %in% genera ~ "Genus"),
                                Rank = factor(Rank,levels = c("Phylum","Family","Genus"))) -> forest_data

forestplot(
  df = prim_res ,
  name = Taxon,
  estimate = estimate,
  se = se,
  pval = qval,
  colour = study
) + ggforce::facet_col(
  facets = ~ Rank,
  space = "free",
  scales = "free_y"
) + theme(legend.position = "bottom") -> p


ggsave(p,filename = "Plots/ANCOM_meta.PDF",device = "pdf",width = 8.2,height = 11.6,units = "in",dpi = 600)
