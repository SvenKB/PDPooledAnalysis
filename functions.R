#### Specify custom functions ####

##########################
#### Helper functions ####
##########################

`%notin%` <- Negate(`%in%`)
`not.na` <- Negate(is.na)


author.names <- c("Pietrucci et al., 2019",
                  "Qian et al., 2018",
                  "Keshavarzian et al., 2015",
                  "Aho et al., 2019",
                  "Heintz-Buschart et al., 2018",
                  "Weis et al., 2019",
                  "Barichella et al., 2019",
                  "Wallen et al., 2020a",
                  "Wallen et al., 2020b")


####################################
#### Data preparation functions ####
####################################

tax <- c('Kingdom','Phylum','Class','Order','Family','Genus',"sequence")

prepare_taxonomy <- function(data,indicator="A") {
  d <- data %>% dplyr::select(ID,all_of(tax)) %>% dplyr::rename("Sequence" = "sequence") %>% mutate(ID = paste0("sp",ID))
  return(d)
}


prepare_taxonomy_sameregion <- function(x) {
  x %>%
    dplyr::select("ID"=X1,"Kingdom"=X3,"Phylum"=X6,"Class"=X9,"Order"=X12,"Family"=X15,"Genus"=X18) %>%
    mutate(ID = paste0("sp",sub("^[^_]*_", "", ID)))
  }

## Transform the data into a tidy format
tidy_up <- function(data,trim=FALSE,indicator="A") {
  d <- data %>%
    dplyr::select(-all_of(tax)) %>%
    rownames_to_column %>% 
    tidyr::gather(var, value, -rowname) %>% 
    tidyr::spread(rowname, value) %>%
    filter(var!="ID") %>%
    setNames(c("ID",paste0(indicator,names(.[,-1]))))
  
  if(trim==TRUE) {d$ID <- unlist(strsplit(d$ID, split='-', fixed=TRUE))[seq(2, nrow(d)*2, 2)]}
  
  return(d)
}


## Load JSON from ENA 
loadENAMeta <- function(path) {
  d <- rjson::fromJSON(file=path,simplify = F)
  col_name <- names(d[[1]])
  
  meta <- as_tibble(data.frame(matrix(unlist(d), nrow=length(d), byrow=T),stringsAsFactors=FALSE))
  colnames(meta) <- col_name
  return(meta)
}
  

##  Data inspection

inspect_depth <- function(d) {
  names <- d %>% dplyr::select(ID) %>% pull
  d <- d %>% dplyr::select(-ID)
  depth <- apply(d,1,sum)
  names(depth) <- names
  return(depth)
}

inspect_reads <- function(d) {
  d <- d %>% dplyr::select(-ID)
  count <- apply(d,2,sum)
  return(count)
}

inspect_sparsity <- function(d,sp=c("asv")) {
  names <- d %>% dplyr::select(ID) %>% pull
  d <- d %>% dplyr::select(-ID)
  
  if(sp=="obs")   {spar <- apply(d,1,function(x) (sum(x==0)/length(x))*100)}
  if(sp=='asv')   {spar <- apply(d,2,function(x) (sum(x==0)/length(x))*100)}
  if(sp=='total') {spar <- sum(d==0) / (nrow(d)*ncol(d))*100}
 
  return(spar)
 }
 
prepare_DA_data <- function(data,tax="Genus",preval =.2 ,meta=c("status","author"),relab=F) {
   
   d <- lapply(data, function(x) as.data.frame(t(phyloseq::otu_table(aggregate_taxa(x,tax)))) %>% rownames_to_column("ID"))
   full_d <- reduce(d,full_join)
   full_d[is.na(full_d)] <- 0
   
   full_ID <- full_d[,1] # extraxct patient IDs
   rownames(full_d) <- full_ID # set rownames to patient IDs
   d <- full_d[,-1] # filter ID variable
   d <- d[,apply(d,2,sum)!=0] # Filter empty taxa
   
   #### Filter minimum prevalence ####
   d1 <- d
   d1[d>0] <- 1
   selected <- colSums(d1)/nrow(d1) > preval
   d <- d[,selected]
   
   
   filt <- apply(d,1,sum) != 0
   d <- d[filt,] # Filter observations without any OTU count
   d <- d %>% rownames_to_column("ID") %>% mutate(N=rowSums(.[,-1])) 
   
   ## Relative abundance
   if(relab ==TRUE) {
     N <- d$N
     ID <- d$ID
     d$N <- NULL
     d$ID <- NULL
     
     d <- data.frame(t(apply(d,1,function(x) x / sum(x))))
     d$ID <- ID
     d$N <- N
   }
   
   ## Extract metadata
   met <- lapply(data,function(x) sample_data(x) %>% data.frame() %>% dplyr::select(any_of(meta)) %>% rownames_to_column("ID"))
   met1 <- reduce(met,full_join) %>% filter(filt)
   
   out.dat <- reduce(list(met1,d),full_join)
   names(out.dat) <- make.names(names(out.dat)) # make tidy variable names
   
   return(out.dat)
 }

 
prepare_beta_data <- function(dat,meta=c("status","author"),tax="Genus",method="bray",filter_id = c()) {
   
   ## Combine data to a single dataset
   d <- lapply(dat, function(x) as.data.frame(t(otu_table(aggregate_taxa(x,tax)))) %>% rownames_to_column("ID"))
   full_d <- reduce(d,full_join)
   full_d[is.na(full_d)] <- 0 # Set NA OTU counts to 0, these are taxa which where not detected in other datasets
   full_d <- full_d %>% filter(ID %notin% filter_id) # Filter IDs
   
   
   full_ID <- full_d[,1] # extraxct patient IDs
   rownames(full_d) <- full_ID # set rownames to patient IDs
   d <- full_d[,-1] # filter ID variable
   filt <- apply(d,1,sum) != 0
   d <- d[filt,] # Filter observations without any OTU count
   
   ## Extract metadata
   met <- lapply(dat,function(x) sample_data(x) %>% data.frame() %>% dplyr::select(meta) %>% rownames_to_column("ID") %>% filter(ID %notin% filter_id))
   met_N <- lapply(dat,function(x) data.frame("N"=sample_sums(x))  %>% rownames_to_column("ID") %>% filter(ID %notin% filter_id))
   met11 <- reduce(met,full_join) %>% filter(filt)
   met_N1 <- reduce(met_N,full_join) %>% filter(filt) 
   met1 <- met11 %>% left_join(met_N1,by="ID")
   
   ## Calculate distance matrix
   if(method=="aitchison") {
     clr_d <- compositions::clr(d+1)
     dist <- vegdist(d,method="euclidean")
   } else {
      dist <- vegdist(d,method=method)
   }
   
   out <- list("dist"=dist,
               "meta"=met1,
               "data"=d)
   return(out)
   
}
 
estimate_alpha <- function(dat,tax="Sequence",meta=c("author")) {
  
  dat <- lapply(dat,function(x) aggregate_taxa(x,tax))
  res <- lapply(dat,function(x) {data.frame("richness"=hillR::hill_taxa(comm=otu_table(x),MARGIN=2,q=0),
                                            "shannon"=round(hillR::hill_taxa(comm=otu_table(x),MARGIN=2,q=1)),
                                            "inv.simpson"=round(hillR::hill_taxa(comm=otu_table(x),MARGIN=2,q=2)))})
  met <- lapply(dat,function(x) sample_data(x) %>% data.frame() %>% dplyr::select(author) %>% rownames_to_column("ID"))
  met1 <- reduce(met,full_join)
  ls <- unlist(lapply(dat,function(x)  rowSums(t(otu_table(x))))) %>% data.frame("N"=.) %>% rownames_to_column("ID")
  
  res <- foreach::foreach(i=1:length(dat)) %do%  {res[[i]] %>% dplyr::mutate(ID=rownames(sample_data(dat[[i]])))}
  
  res <- foreach::foreach(i=1:length(dat)) %do% {dat[[i]] %>% phyloseq::sample_data() %>% plyr::mutate(ID=rownames(.)) %>% data.frame() %>% dplyr::select(ID,status) %>% left_join(res[[i]]) %>% plyr::mutate(study=i)}
  
  out <- res %>% reduce(full_join)
  out <- out %>% left_join(ls,by="ID")
  out <- out %>% left_join(met1,by="ID")
  out <- out %>% mutate(author = fct_reorder(author,study))
  return(out)
  
}


extrapolate_alpha <- function(dat,tax="Sequence",meta=c("author")) {
  
  dat <- lapply(dat,function(x) aggregate_taxa(x,tax))
  res <- lapply(dat,function(x) {data.frame("richness"=hillR::hill_taxa(comm=otu_table(x),MARGIN=2,q=0),
                                            "shannon"=round(hillR::hill_taxa(comm=otu_table(x),MARGIN=2,q=1),0),
                                            "inv.simpson"=round(hillR::hill_taxa(comm=otu_table(x),MARGIN=2,q=2)),0)})
  met <- lapply(dat,function(x) sample_data(x) %>% data.frame() %>% dplyr::select(author) %>% rownames_to_column("ID"))
  met1 <- reduce(met,full_join)
  ls <- unlist(lapply(dat,function(x)  rowSums(t(otu_table(x))))) %>% data.frame("N"=.) %>% rownames_to_column("ID")
  
  res <- foreach::foreach(i=1:length(dat)) %do%  {res[[i]] %>% dplyr::mutate(ID=rownames(sample_data(dat[[i]])))}
  
  res <- foreach::foreach(i=1:length(dat)) %do% {dat[[i]] %>% phyloseq::sample_data() %>% plyr::mutate(ID=rownames(.)) %>% data.frame() %>% dplyr::select(ID,status) %>% left_join(res[[i]]) %>% plyr::mutate(study=i)}
  
  out <- res %>% reduce(full_join)
  out <- out %>% left_join(ls,by="ID")
  out <- out %>% left_join(met1,by="ID")
  out <- out %>% mutate(author = fct_reorder(author,study))
  return(out)
  
}



#################################
#### Data cleaning functions ####
#################################
 
 
## Currently, the function only filters based on relative abundance, not taxonomic information, as taxonomic information is present in almost all ASVs
 
 
clean <- function(data,tax,min_tax="Family",length=150, filter=c("abs","rel","quant"),abs=100,rel=0.005) {
  
  id <- data %>% dplyr::select(ID)
  
  ## Filter kingdoms
  tot_asv <- nrow(tax) # Track full number of ASVs
  
  ft_king <- tax %>% filter(Kingdom=="Bacteria") %>% dplyr::select(ID) %>% pull # Create Filter
  ft_king_nr <- tot_asv - length(ft_king) # Calculate nr. of filtered ASVs
  
  tax <- tax %>% filter(Kingdom=="Bacteria") # Apply filter to taxonomy 
  new_asv <- nrow(tax) # new numbeer of ASVs
  
  d <- data %>% dplyr::select(any_of(ft_king)) # Apply filter to data
  
  ## Filter Chloroplasts
  ft_chloroplast <- tax %>% filter(Class!="Chloroplast") %>% dplyr::select(ID) %>% pull # Create Filter
  ft_chloroplast_nr <- new_asv - length(ft_chloroplast) # Calculate nr. of filtered ASVs
  
  tax <- tax %>% filter(Class!="Chloroplast") # Apply filter to taxonomy 
  new_asv <- nrow(tax) # new numbeer of ASVs
  d <- d %>% dplyr::select(any_of(ft_chloroplast)) # Apply filter to data
  
  ## Filter based on low abundance AND/OR low taxonomic information
  
  abund <- apply(d,2,sum) # Calculate Abundance
  #samp <- t(apply(d,1,function(x) x/sum(x)))
  #
  #filter_samp <- function(x,r=rel,prop=0.3) {
  #  n <- length(x)
  #  s <- sum(x>r)
  #  filt <- ifelse(s/n < 0.1,TRUE,FALSE)
  #  return(filt)
  #}
  
  if(filter=="quant") {ft_abund <- tax %>% filter(abund<summary(abund)[[3]]) %>% dplyr::select(ID) %>% pull # Create filter based on quantile
                       thresh <- summary(abund)[[3]]}
  if(filter == "abs") {ft_abund <- tax %>% filter(abund<abs) %>% dplyr::select(ID) %>% pull # Create filter based on absolute reads
                       thresh <- abs}
  if(filter == "rel") {ind <- (abund/sum(d))*100<rel                                        # Create filter based on relative abundance
                       ind <- names(ind[ind==TRUE])
                       ft_abund <- tax %>% filter(ID %in% ind) %>% dplyr::select(ID) %>% pull
                       thresh <- paste0(rel,"%")}
  
  # Apply combined filter -> only if abundance threshold AND missing annotation are are present, ASV will be filtered !: MAYBE BAD IDEA?
  # %>% filter(is.na(!!as.symbol(min_tax)))
  
  ft_mintax <- tax %>% filter(ID%in%ft_abund|is.na(!!as.symbol(min_tax))) %>% dplyr::select(ID) %>% pull
  ft_mintax_nr <- length(ft_mintax)
  
  ft_abund_nr <- length(ft_abund)
  
  tax <- tax %>% filter(ID %notin% ft_mintax)
  new_asv <- nrow(tax)
  d <- d %>% dplyr::select(-any_of(ft_mintax))
  
  
  ## Filter based on sequence length
  
  seq <- tax %>% dplyr::select(Sequence) %>% convert(chr(Sequence)) %>% pull
  seq_length <- sapply(str_split(seq,''),length)
  ft_seqlength <- tax %>% filter(seq_length > length) %>% dplyr::select(ID) %>% pull
  ft_seqlength_nr <- nrow(tax)-length(ft_seqlength)
  tax <- tax %>% filter(ID %in% ft_seqlength)
  new_asv <- nrow(tax)
  
  d <- d %>% dplyr::select(any_of(ft_seqlength))
  
  tot_filt <- ft_king_nr + ft_chloroplast_nr + ft_mintax_nr + ft_seqlength_nr
  
  d <- data.frame(id,d)
  
  cat('Total number of ASVs:',tot_asv,"\n",
      ft_king_nr,'ASVs filtered due to wrong Kingdom\n',
      ft_chloroplast_nr,'ASVs of Class Chloroplast filtered\n',
      ft_mintax_nr,'ASVs filtered due to low abundance (<',thresh,')\n',
      #"Number of ASVs below threshold, but not filtered:",ft_abund_nr,"\n",
      ft_seqlength_nr,'ASVs filtered due to low sequence length < 150 \n',
      'Number of filtered ASVs:',tot_filt,'\n',
      'New number of ASVs:',new_asv,"\n")
  return(d)
}


###################################
#### Data exploration function ####
###################################


## calculate number of ASVs assigned to one Phylum
features_per_phylum <- function(phyl) {table(tax_table(phyl)[, "Phylum"], exclude = NULL)}


## Calculate average abundance of all ASVs in a Phylum and total abundance
prevalence_per_phylum <- function(phyl) {
  ps <- subset_taxa(phyl,!is.na(Phylum))
  # Compute prevalence of each feature, store as data.frame
  prevdf = apply(X = otu_table(ps),
                 MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps),
                      tax_table(ps))
  prevdf <- prevdf %>% arrange(Phylum)
  prev <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
  colnames(prev) <- c("Phylum","Avg. Prevalence", "Total Prevalence")
  return(prev)
}


## Plot the fraction of samples in which a taxon appears against its abundance by Phylum
prevalence_by_abundance <- function(phyl,thrs=0.05) {
  
  n <- phyloseq::nsamples(phyl)
  threshold <- thrs*n
  ps <- subset_taxa(phyl,!is.na(Phylum))
  prevdf = apply(X = otu_table(ps),
                 MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps),
                      tax_table(ps))
  
  keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= threshold)]
  
  ps0 <- prune_taxa(keepTaxa,ps)
  prevdf = apply(X = otu_table(ps0),
                 MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps0),
                      tax_table(ps0))
  
  
  prevdf <- prevdf %>%
    arrange(Phylum) %>%
    mutate(prev.frac = Prevalence/n)
  
  ggplot(aes(x=TotalAbundance,y=prev.frac,color=Phylum),data=prevdf) +
    geom_jitter(alpha=1/2) +
    facet_wrap(~Phylum,scales = "free") +
    theme_fira() +
    #scale_color_fira() +
    scale_x_log10() +
    xlab("Total abundance[log]") +
    ylab("Prevalence (fraction of samples)") +
    theme(legend.position="none")
}


#### Inspect sparsity

sparsityHeatmap <- function(phyl,level) {
  
  id <- rownames(data.frame(t(otu_table(aggregate_taxa(phyl,level=level)))))  
  d <- as.data.frame(t(otu_table(aggregate_taxa(phyl,level=level))))
  d$N <- rowSums(d)
  
  #rownames(d) <- id
  
  
  # Prepare zero-inflation indicator matrix
  d.mat <- d
  d.mat[d.mat>0] <- 1
  d.mat <- d.mat %>% dplyr::select(-N)
  d.mat <- d.mat[order(rowSums(d.mat),decreasing = T),order(colSums(d.mat),decreasing = T)]
  
  # Prepare sequencing depth vector
  N <- d[order(rowSums(d.mat),decreasing = T),"N"]
  
  rowDat <- data.frame(N=N)
  
  rowAnn <- HeatmapAnnotation(Depth = anno_barplot(as.matrix(rowDat$N), border = F,axis_param = list(direction = "reverse"),width = unit(4, "cm")),
                              which = "row")
  
  ComplexHeatmap::Heatmap(as.matrix(d.mat),
                          col = circlize::colorRamp2(breaks = c(0,1),colors = c("#D60C00FF","#00468BFF"),space="sRGB"),
                          cluster_columns = F,
                          cluster_rows = F, 
                          show_column_names = F,
                          show_row_names = T,
                          left_annotation = rowAnn,
                          heatmap_legend_param = list(
                            title = "Observed", at = c(0,1), 
                            labels = c("No","Yes")),
                          rect_gp = gpar(col = "black", lwd = .002),
                          row_names_gp = gpar(fontsize = 6))
  
}


##########################
#### Helper functions ####
##########################


addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

g.mean <- function(x,seq.depth=1,weighted=F, na.rm=TRUE){
  w <- scale(N)
  if(weighted == F) {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  } else {
    exp(sum((w*log(x[x > 0])), na.rm=na.rm) / sum(w))
  }
}
