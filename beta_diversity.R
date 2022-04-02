library(mia)
library(miaViz)
library(vegan)
library(scater)
library(tidyverse)
library(reshape2)


dat <- readRDS("Data/processed_data/merged_phyl.RDS")

otu <- phyloseq::otu_table(dat)
meta <- phyloseq::sample_data(dat)
taxa <- phyloseq::tax_table(dat)
tree <- phyloseq::phy_tree(dat)


tse <- TreeSummarizedExperiment(assays = list(counts = otu),
                                colData = meta,
                                rowData = taxa,
                                rowTree = tree)

colData(tse)$type <- "trimmed"
colData(tse)$depth <- colSums(assays(tse)$counts)


tse_FS <- readRDS("Data/processed_data/merged_tse_FS.RDS")
colData(tse_FS)$type <- "Full"


#### Alpha diversity 
colData(tse)$Richness <- hillR::hill_taxa(comm=assays(tse)$counts,MARGIN=2,q=0)
colData(tse)$Shannon <-  hillR::hill_taxa(comm=assays(tse)$counts,MARGIN=2,q=1)
colData(tse)$Simpson <-  hillR::hill_taxa(comm=assays(tse)$counts,MARGIN=2,q=2)

colData(tse_FS)$Richness <- hillR::hill_taxa(comm=assays(tse_FS)$counts,MARGIN=2,q=0)
colData(tse_FS)$Shannon<-  hillR::hill_taxa(comm=assays(tse_FS)$counts,MARGIN=2,q=1)
colData(tse_FS)$Simpson<-  hillR::hill_taxa(comm=assays(tse_FS)$counts,MARGIN=2,q=2)

as.data.frame(colData(tse)) %>%
  dplyr::select(author,status,region,Richness,Shannon,Simpson) %>%
  melt() %>%
  ggplot(aes(x=author,y=value,fill=status)) +
  geom_boxplot() +
  facet_wrap(~variable,scales="free",nrow=1) +
  theme(axis.text.x = element_text(angle=90))

as.data.frame(colData(tse_FS)) %>%
  dplyr::select(author,status,region,Richness,Shannon,Simpson) %>%
  melt() %>%
  ggplot(aes(x=author,y=value,fill=status)) +
  geom_boxplot() +
  facet_wrap(~variable,scales="free",nrow=1) +
  theme(axis.text.x = element_text(angle=90))


as.data.frame(colData(tse)) %>%
  dplyr::select(author,status,region,Richness,Shannon,Simpson) %>%
  melt() %>%
  ggplot(aes(x=region,y=value,fill=author)) +
  geom_boxplot() +
  facet_wrap(~variable,scales="free") +
  theme(axis.text.x = element_text(angle=90))


trim <- as.data.frame(colData(tse)) %>%
  dplyr::select(author,status,region,type,Richness,Shannon,Simpson) %>%
  melt()

full <- as.data.frame(colData(tse_FS)) %>%
  dplyr::select(author,status,region,type,Richness,Shannon,Simpson) %>%
  melt()


trim %>%
  bind_rows(full) %>%
  ggplot(aes(x=region,y=value,fill=type)) +
  geom_boxplot() +
  facet_wrap(~variable,scales="free") +
  theme(axis.text.x = element_text(angle=90))



#calculateUniFrac(tse,weighted=T)


uf_dist <- phyloseq::distance(mia::makePhyloseqFromTreeSummarizedExperiment(tse), method = "unifrac")
wuf_dist <- phyloseq::distance(mia::makePhyloseqFromTreeSummarizedExperiment(tse), method = "wunifrac")
bc_dist <- phyloseq::distance(mia::makePhyloseqFromTreeSummarizedExperiment(tse), method = "bray")

otu <- t(assays(tse)$counts)

otu <- data.frame(t(apply(otu,1, function(x) x/sum(x))))
otu <- compositions::clr(otu)

clr_dist <- dist(otu)


uf_pcoa <- ecodist::pco(uf_dist)
uf_explained_var <- round((uf_pcoa$values/sum(uf_pcoa$values))*100,2)

wuf_pcoa <- ecodist::pco(wuf_dist)
wuf_explained_var <- round((wuf_pcoa$values/sum(wuf_pcoa$values))*100,2)

bc_pcoa <- ecodist::pco(bc_dist)
bc_explained_var <- round((bc_pcoa$values/sum(bc_pcoa$values))*100,2)

clr_pcoa <- ecodist::pco(clr_dist)
clr_explained_var <- round((clr_pcoa$values/sum(clr_pcoa$values))*100,2)



uf_pcoa_df <- data.frame("PCo1" = uf_pcoa$vectors[,1],
                         "PCo2" = uf_pcoa$vectors[,2],
                         "status" = colData(tse)$status,
                         "region" = colData(tse)$region,
                         "author" = colData(tse)$author)

wuf_pcoa_df <- data.frame("PCo1" = wuf_pcoa$vectors[,1],
                         "PCo2" = wuf_pcoa$vectors[,2],
                         "status" = colData(tse)$status,
                         "region" = colData(tse)$region,
                         "author" = colData(tse)$author)

bc_pcoa_df <- data.frame("PCo1" = bc_pcoa$vectors[,1],
                          "PCo2" = bc_pcoa$vectors[,2],
                          "status" = colData(tse)$status,
                          "region" = colData(tse)$region,
                          "author" = colData(tse)$author)

clr_pcoa_df <- data.frame("PCo1" = clr_pcoa$vectors[,1],
                         "PCo2" = clr_pcoa$vectors[,2],
                         "status" = colData(tse)$status,
                         "region" = colData(tse)$region,
                         "author" = colData(tse)$author)

## Status

p1 <- ggplot(aes(x=PCo1,y=PCo2,color=status),data=bc_pcoa_df) +
  geom_point() +
  theme_pubr() +
  xlab(paste0("PC1 ",bc_explained_var[1],"%")) +
  ylab(paste0("PC2 ",bc_explained_var[2],"%")) + 
  ggtitle("Bray-Curtis dissimilarity") +
  scale_color_discrete(name="Disease status")


p2 <- ggplot(aes(x=PCo1,y=PCo2,color=status),data=wuf_pcoa_df) +
  geom_point() +
  theme_pubr() +
  xlab(paste0("PC1 ",wuf_explained_var[1],"%")) +
  ylab(paste0("PC2 ",wuf_explained_var[2],"%")) + 
  ggtitle("Weighted UniFrac") +
  scale_color_discrete(name="Disease status")

p3 <- ggplot(aes(x=PCo1,y=PCo2,color=status),data=clr_pcoa_df) +
  geom_point() +
  theme_pubr() +
  xlab(paste0("PC1 ",clr_explained_var[1],"%")) +
  ylab(paste0("PC2 ",clr_explained_var[2],"%")) + 
  ggtitle("CLR Euclidean") +
  scale_color_discrete(name="Disease status")



p11 <- ggarrange(p1,p2,p3,common.legend = T,legend="right",nrow=1)


## Author

p1 <- ggplot(aes(x=PCo1,y=PCo2,color=author),data=bc_pcoa_df) +
  geom_point() +
  theme_pubr() +
  xlab(paste0("PC1 ",bc_explained_var[1],"%")) +
  ylab(paste0("PC2 ",bc_explained_var[2],"%")) + 
  #ggtitle("Bray-Curtis dissimilarity") +
  scale_color_discrete(name="Study") 


p2 <- ggplot(aes(x=PCo1,y=PCo2,color=author),data=wuf_pcoa_df) +
  geom_point() +
  theme_pubr() +
  xlab(paste0("PC1 ",wuf_explained_var[1],"%")) +
  ylab(paste0("PC2 ",wuf_explained_var[2],"%")) + 
  #ggtitle("Weighted UniFrac") +
  scale_color_discrete(name="Study")


p3 <- ggplot(aes(x=PCo1,y=PCo2,color=author),data=clr_pcoa_df) +
  geom_point() +
  theme_pubr() +
  xlab(paste0("PC1 ",clr_explained_var[1],"%")) +
  ylab(paste0("PC2 ",clr_explained_var[2],"%")) + 
  #ggtitle("Weighted UniFrac") +
  scale_color_discrete(name="Study")


p22 <- ggarrange(p1,p2,p3,common.legend = T,legend="right",nrow=1)


## Region

p1 <- ggplot(aes(x=PCo1,y=PCo2,color=region),data=bc_pcoa_df) +
  geom_point() +
  theme_pubr() +
  xlab(paste0("PC1 ",bc_explained_var[1],"%")) +
  ylab(paste0("PC2 ",bc_explained_var[2],"%")) + 
  #ggtitle("Bray-Curtis dissimilarity") +
  scale_color_discrete(name="Variable region")


p2 <- ggplot(aes(x=PCo1,y=PCo2,color=region),data=wuf_pcoa_df) +
  geom_point() +
  theme_pubr() +
  xlab(paste0("PC1 ",wuf_explained_var[1],"%")) +
  ylab(paste0("PC2 ",wuf_explained_var[2],"%")) + 
  #ggtitle("Weighted UniFrac") +
  scale_color_discrete(name="Variable region")

p3 <- ggplot(aes(x=PCo1,y=PCo2,color=region),data=clr_pcoa_df) +
  geom_point() +
  theme_pubr() +
  xlab(paste0("PC1 ",clr_explained_var[1],"%")) +
  ylab(paste0("PC2 ",clr_explained_var[2],"%")) + 
  #ggtitle("Weighted UniFrac") +
  scale_color_discrete(name="Variable region")

p33 <- ggarrange(p1,p2,p3,common.legend = T,legend="right",nrow=1)


p <- ggarrange(p11,p22,p33,ncol=1)

ggsave(p,filename = "Plots/beta_div_noleg.PDF",device = "pdf",width = 12,height = 8.2,units = "in",dpi = 600)

### Appendix

### wuf

ggplot(aes(x=PCo1,y=PCo2,color=status),data=wuf_pcoa_df) +
  geom_point() +
  facet_wrap(~author,scales="free") +
  theme_pubr() +
  xlab("PC1") +
  ylab("PC2") + 
  ggtitle("Weighted UniFrac") +
  scale_colour_jama() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())  +
  guides(colour=guide_legend(title="Disease status")) -> p1

ggplot(aes(x=PCo1,y=PCo2,color=status),data=bc_pcoa_df) +
  geom_point() +
  facet_wrap(~author,scales="free") +
  theme_pubr() +
  xlab("PC1") +
  ylab("PC2") + 
  ggtitle("Bray-Curtis dissimilarity") +
  scale_colour_jama() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  guides(colour=guide_legend(title="Disease status")) -> p2


ggplot(aes(x=PCo1,y=PCo2,color=status),data=clr_pcoa_df) +
  geom_point() +
  facet_wrap(~author,scales="free") +
  theme_pubr() +
  xlab("PC1") +
  ylab("PC2") + 
  ggtitle("CLR Euclidean") +
  scale_colour_jama() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  guides(colour=guide_legend(title="Disease status")) -> p3

ggsave(p1,filename = "Plots/beta_authors_wuf.PDF",device = "pdf",width = 8.3,height = 8.3,units = "in",dpi = 600)

ggsave(p2,filename = "Plots/beta_authors_bc.PDF",device = "pdf",width = 8.3,height = 8.3,units = "in",dpi = 600)

ggsave(p3,filename = "Plots/beta_authors_clr.PDF",device = "pdf",width = 8.3,height = 8.3,units = "in",dpi = 600)


####

com <- assays(tse)$counts+1
env <- as.data.frame(colData(tse)[c("region")])

#convert com to a matrix
m_com = as.matrix(t(com))

#nmds code
set.seed(123)

nmds = metaMDS(uf_dist,k=2)

en = envfit(nmds, env, permutations = 999, na.rm = TRUE)

plot(nmds)
plot(en)
