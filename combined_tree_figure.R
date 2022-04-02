library(mia)
library(miaViz)
library(vegan)
library(scater)
library(tidyverse)
library(reshape2)
library(ggtree)
library(treeio)
library(phyloseq)
library(ggsci)
library(ggnewscale)
library(ggtreeExtra)
library(hablar)
library(ggsci)
library(ggpubr)

## Function

`%notin%` <- Negate(`%in%`)

author.names <- c("Pietrucci et al., 2019",
                  "Qian et al., 2018",
                  "Keshavarzian et al., 2015",
                  "Aho et al., 2019",
                  "Heintz-Buschart et al., 2018",
                  "Weis et al., 2019",
                  "Barichella et al., 2019",
                  "Wallen et al., 2020a",
                  "Wallen et al., 2020b")

## Load data
dat <- readRDS("Data/processed_data/merged_tse.RDS")
res_help <- readRDS("Data/processed_data/result_quick.RDS")

## Prepare identifier
phylum <- res_help %>% filter(Rank == "Phylum") %>% pull(Outcome)
family <- res_help %>% filter(Rank == "Family") %>% pull(Outcome)
genus <- res_help %>% filter(Rank == "Genus") %>% pull(Outcome)


seq1 <- rownames(rowData(dat)[rowData(dat)[,"Phylum"] %in% phylum,])
seq2 <- rownames(rowData(dat)[rowData(dat)[,"Family"] %in% family,])
seq3 <- rownames(rowData(dat)[rowData(dat)[,"Genus"] %in% genus,])
seqs <- unique(c(seq1,seq2,seq3))

#### Taxonomy as tree data

phyl <- mia::makePhyloseqFromTreeSummarizedExperiment(dat)
custom_tree <- read_delim("Data/other/custom_result_tree.csv",delim = ";",col_types = list("d","d","d","c","f","c","l","f"))

result_tax <- readxl::read_xlsx("Data/other/selected_tax.xlsx")


new_phyl <- subset_taxa(phyl, Genus %in% result_tax$Genus)

taxa <- tax_table(new_phyl)[rownames(tax_table(new_phyl)) %in% seqs,]

new_tree <- MicrobiotaProcess::as.treedata(taxa)

grp <- taxa %>%
  data.frame() %>%
  dplyr::select(Phylum,Genus) %>%
  pivot_wider(names_from = Phylum,values_from = Genus) %>% as.list() %>%
  lapply(unlist) %>%
  lapply(function(x) paste0("g__",x))

grp <- lapply(grp,unique)

new_tree <- groupOTU(new_tree,grp)

selected <- c(paste0("p__",phylum),
              paste0("f__",family),
              paste0("g__",genus))


selected <- c(paste0("p__",c("Synergistetes","Verrucomicrobia"," ")),
              paste0("f__",c("Synergistaceae","Verrucomicrobiaceae","Peptococcaceae.2",
                             "Porphyromonadaceae","Lachnospiraceae")),
              paste0("g__",c("Eisenbergiella","Akkermansia","Desulfurispora",
                             "Acidaminobacter","Faecalibacterium","Roseburia")))




meta <- data.frame("label" = unique(c("r_root",paste0("k__",taxa[,"Kingdom"]),
                                      paste0("p__",taxa[,"Phylum"]),
                                      paste0("c__",taxa[,"Class"]),
                                      paste0("o__",taxa[,"Order"]),
                                      paste0("f__",taxa[,"Family"]),
                                      paste0("g__",taxa[,"Genus"])))) %>% mutate(selected = label %in% selected,
                                                                                 selected = ifelse(selected==T,"Yes","No"))

new_tree %>%
  as_tibble() %>%
  left_join(meta) %>%
  mutate(selected = ifelse(label=="r__root","No",selected)) %>%
  mutate(label = case_when(label=="f__Enterobacteriaceae" ~ "ff__Enterobacteriaceae",
                           label=="f__Synergistaceae" ~ "ff__Synergistaceae",
                           label=="g__Cloacibacillus" ~ "f__Synergistaceae",
                           label=="g__Enterobacter" ~ "f__Enterobacteriaceae",
                           TRUE ~ label)) %>%
  mutate(show = case_when(label=="f__Enterobacteriaceae" ~ FALSE,
                          label=="f__Synergistaceae" ~ FALSE,
                          TRUE ~ TRUE)) %>% 
  as.treedata() -> final_tree


## Circular tree

## Prepare data to indicate whether study was included in analysis

dat_gen <- agglomerateByRank(dat,"Genus")
dat_fam <- agglomerateByRank(dat,"Family")
met <- colData(dat) %>% data.frame() %>% dplyr::select(author,status) %>% rownames_to_column("ID")

genera <- data.frame(t(assays(dat_gen)$counts)) %>% dplyr::select(all_of(genus)) %>% rename_all(function(x) paste0("g__",x))
genera <- data.frame(ifelse(genera > 0,1,0)) 

families <- data.frame(t(assays(dat_fam)$counts)) %>%
  dplyr::select(all_of(c("Synergistaceae","Enterobacteriaceae"))) %>%
  rename_all(function(x) paste0("f__",x))
families <- data.frame(ifelse(families > 0,1,0))

taxa <- genera %>%
  bind_cols(families) %>% rownames_to_column("ID")

met %>%
  left_join(taxa) %>%
  dplyr::select(-c(ID,status)) %>%
  group_by(author) %>%
  summarize(across(everything(),sum)) %>%
  mutate(across(where(is.numeric),function(x) ifelse(x>0,1,0))) -> prevalent

prevalent %>%
  mutate(author_num = case_when(author==author.names[[1]] ~ 1,
                                author==author.names[[2]] ~ 2,
                                author==author.names[[3]] ~ 3,
                                author==author.names[[4]] ~ 4,
                                author==author.names[[5]] ~ 5,
                                author==author.names[[6]] ~ 6,
                                author==author.names[[7]] ~ 7,
                                author==author.names[[8]] ~ 8,
                                author==author.names[[9]] ~ 9,)) %>%
  pivot_longer(!c("author","author_num"),names_to = "label") %>% dplyr::select(-value) -> prevalent

used_studies1 <- readRDS("Data/other/studywise_filter_studies.RDS")
used_studies2 <- readRDS("Data/other/studywise_filter_studies_fam.RDS")


names(used_studies1) <- paste0("g__",names(used_studies1))
names(used_studies2) <- paste0("f__",names(used_studies2))


used_studies <- c(used_studies1,used_studies2) 


maxl <- max(sapply(used_studies, length))

used_studies <- sapply(used_studies, FUN = function(x, ml) {
  difference <- ml - length(x)
  c(x, rep(NA, difference))
}, ml = maxl, simplify = FALSE)

studies_df <- data.frame(used_studies)

used_studies %>%
  lapply(function(x) author.names %in% x) %>%
  data.frame() %>%
  mutate(author = author.names,.before=1,
         f__Enterobacteriaceae = rep(TRUE,9)) %>%
  melt("author") %>%
  mutate(Included = ifelse(value==T,"Yes","No"),
         label = variable,.keep="unused") %>%
  convert(chr(author,label,Included)) %>%
  as_tibble() -> studies_df


prevalent %>%
  left_join(studies_df,by=c("label","author")) -> prevalent



#prevalent$random <- sample(1:9,nrow(prevalent),replace = T)


# Prepare abundance data  

genera <- data.frame(t(assays(dat_gen)$counts)) %>%
  dplyr::select(all_of(genus)) %>%
  rename_all(function(x) paste0("g__",x))

families <- data.frame(t(assays(dat_fam)$counts)) %>%
  dplyr::select(all_of(c("Synergistaceae","Enterobacteriaceae"))) %>%
  rename_all(function(x) paste0("f__",x))

taxa <- genera %>%
  bind_cols(families) %>%
  apply(1,function(x) log((x+1)/sum(x+1)) ) %>%
  t() %>%
  data.frame() %>% 
  rownames_to_column("ID")


met %>%
  left_join(taxa) %>%
  dplyr::select(-ID) %>%
  melt() %>% mutate(label = as.character(variable),.keep="unused")-> abundance


#### Prepare draws for ggdist ####

draws <- readRDS("Data/other/results_draws.RDS")


taxa <- c("f__Synergistaceae",
          "g__Porphyromonadaceae","g__Lachnospiraceae","g__Akkermansia","g__Eisenbergiella","g__Desulfurispora","g__Acidaminobacter",
          "g__Faecalibacterium","g__Roseburia")


library(ggdist)

draws %>%
  mutate(label = Outcome,
         sig = ifelse(Outcome %in% taxa, "Yes","No")) %>%
  dplyr::select(.chain,.iteration,.draw,b_statusPD,author,label,sig) -> draws

draws %>%
  filter(label %in% prevalent$label) %>%
  filter(author == "Pooled Effect") -> draws

draws$label <- as.factor(draws$label)

new_order <- levels(draws$label)[c(6,1,2,4,14,15,13,3,9,5,8,11,7,10,12,16)]


draws$label <- factor(draws$label, levels = new_order)


tab_res <- readRDS(file = "Data/other/tab_dat.RDS") %>%
  dplyr::mutate(label = case_when(Outcome %in% c("Enterobacteriaceae","Synergistaceae") ~ paste0("f__",Outcome),
                                  Outcome %notin% c("Enterobacteriaceae","Synergistaceae") ~ paste0("g__",Outcome))) %>%
  dplyr::filter(label %in% draws$label) %>%
  dplyr::select(label,RR, `99% CR`, Reported, Increased, Decreased)


ggplot(aes(x=b_statusPD,y=label),data=draws ) +
  geom_vline(xintercept = 1,size=1,colour="black") +
  #stat_gradientinterval(aes(slab_alpha=stat(-pmax(abs(1-2*cdf),.95)),fill = stat(x > 0))) +
  #stat_halfeye(aes(slab_alpha=stat(-pmax(abs(1-2*cdf),.95)),fill = stat(x > 0))) +
  #stat_halfeye(aes(fill = stat(exp(x)>1)),position=position_dodge(width=1),alpha=1,color="black",scale=1.2) +
  stat_pointinterval(aes(colour=sig),.width = .99,size=1) +
  scale_fill_jama() +
  theme_forest() +
  coord_cartesian(xlim=c(.4,5)) +
  scale_x_log10() +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.x = element_blank(),
        legend.position = "none",
        plot.margin =  unit(c(1.8,0,1,1), "lines")) +
  geom_stripes() +
  scale_color_jama() -> p1


tab_res$label <- as.factor(tab_res$label)
tab_res$label <- factor(tab_res$label, levels =new_order)

tab_res %>%
  mutate(new_label = glue("{RR} {`99% CR`}")) -> tab_res


ggplot(aes(x=0,y=label,label=new_label),data=tab_res) +
  geom_text(hjust=0,size=5) +
  coord_cartesian(xlim=c(0,0.02)) +
  theme_void() + 
  theme(plot.margin =  unit(c(1.8,1,1.8,1), "lines"))-> p3

tab_res %>%
  dplyr::select(label,Increased,Decreased) %>%
  pivot_longer(!label,names_to = "Reported",values_to = "number")  -> reporting_data


# Circular
p <- 
  ggtree(final_tree,size=1,layout="roundrect",aes(alpha=show)) +
  scale_color_jama() +
  ggnewscale::new_scale_color() +
  geom_point(aes(colour=selected,alpha=show),size=3) +
  scale_colour_jama(name="Differentially abundant") +
  scale_alpha_discrete(range=c(0,1))

#
p <- p +
  new_scale_fill() + 
  geom_fruit(
    data=prevalent,
    geom=geom_tile,
    colour="white",
    mapping=aes(y=label, x=author,fill=Included),
    lwd=1,
    linetype=1,
    offset=.4,   # The distance between external layers, default is 0.03 times of x range of tree.
    pwidth=.6 # width of the external layer, default is 0.2 times of x range of tree.
  ) +
  geom_tiplab(size=4,alpha=1,offset = .5,fontface =2,hjust = 0) +
  scale_fill_jama() +
  theme(plot.margin =  unit(c(0,0,0,0), "lines"))


reporting_data %>%
  mutate(number = if_else(Reported == "Decreased", -number, number)) %>%
  ggplot(aes(x=label,y=number,fill=Reported)) +
  geom_col(width=.5) +
  coord_flip() +
  geom_hline(yintercept = 0,size=2,colour = "black") +
  scale_y_continuous(breaks = seq(-5, 5, 1), 
                     labels = abs(seq(-5, 5, 1))) +
  theme_forest() +   
  theme(panel.spacing.x = unit(0, "mm"),
        legend.position = "top",
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin =  unit(c(1.8,0,1,0), "lines")) +
  scale_fill_jama() -> p_reported

p2 <- ggarrange(p,p_reported,p1,p_ancom,widths=c(1,.25,.5,.5),common.legend = T,nrow=1)


ggsave("Plots/results_summary.PDF",p2,device=cairo_pdf(),width = 17,height = 10,units = "in",dpi = 1200)
