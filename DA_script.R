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

full.dat <- prepare_DA_data(dat,tax=level,meta=c("status","author"),preval = prev)


#### Change author names to avoid errors using stan

name <- unique(full.dat$author)

names <- cbind(paste0("Study",1:9),name)

full.dat <- full.dat %>%
  mutate(author = case_when(author==names[[1,2]] ~ names[[1,1]],
                            author==names[[2,2]] ~ names[[2,1]],
                            author==names[[3,2]] ~ names[[3,1]],
                            author==names[[4,2]] ~ names[[4,1]],
                            author==names[[5,2]] ~ names[[5,1]],
                            author==names[[6,2]] ~ names[[6,1]],
                            author==names[[7,2]] ~ names[[7,1]],
                            author==names[[8,2]] ~ names[[8,1]],
                            author==names[[9,2]] ~ names[[9,1]]))

#### Prepare model function ####
taxa <- setdiff(colnames(full.dat),c("ID","status","author","N"))

indicator_df <-  data.frame(full.dat[,c("ID","status","author","N")],ifelse(full.dat[,taxa] > 0,1,0))

filter_studies <- function(x) {
  indicator_df %>%
    dplyr::select(author,all_of(x)) %>%
    group_by(author) %>%
    summarise(n = n(),
              prev = sum(.data[[x]])) %>%
    filter( prev/n > .2) %>% pull(author) -> study_filter
  return(study_filter)
}

study_filter <- lapply(taxa, filter_studies)
names(study_filter) <- taxa

myrename <- function(x) {for(i in 1:nrow(names)) {
  x[x==names[i,1]] <- names[i,2]
}
  return(x)
}

used_studies <- lapply(study_filter, myrename)

fit_brms <- function(x) {
  
  flt <- study_filter[[x]]
  
  df <- full.dat[full.dat$author %in% flt,]
  
  
  fit <- brms::brm(data = df, family = zero_inflated_negbinomial,
                   as.formula(paste0(x,"~status+offset(log(N))+(status|author)")),
                   prior = c(prior(gamma(0.01, 0.01), class = shape)),  
                   iter = 8000, warmup = 1000, cores = 4, chains = 4,
                   seed = 11,control = list(adapt_delta = .99))
  return(fit)
}


######################
#### Run analyses ####
######################


# Prepare parallel computing
cores <- detectCores()/2
cl <- makeCluster(cores)

clusterExport(cl,c("zero_inflated_negbinomial","full.dat","prior","study_filter","%>%"))
b_models <- parLapply(cl,taxa,fit_brms)

names(b_models) <- taxa

stopCluster(cl)

output <- paste0("Data/models/",level,"_DA_results_studywise_filter_selected.RDS")

saveRDS(b_models,output)




