#############################################
#### Data preparation and quality checks ####
#############################################


#### load packages ####
source('packages.R')

#### load functions ####
source("functions.R")


#### Load data ####
studies <- list("A_Pietrucci_2019",
                "B_Qian_2018",
                "C_Keshavarzian_2015",
                "E_Scheperjans_2019",
                "F_Buschart_2018",
                "H_Weis_2019",
                "I_Barichella_2019",
                "L_Wallen_2020_1",
                "L_Wallen_2020_2")


path <- "Data/count_tables/"

# Load ASV data
A_data <- readr::read_delim(paste0(path,studies[[1]],".txt"),delim="\t")
B_data <- readr::read_delim(paste0(path,studies[[2]],".txt"),delim="\t")
C_data <- readr::read_delim(paste0(path,studies[[3]],".txt"),delim="\t")
E_data <- readr::read_delim(paste0(path,studies[[4]],".txt"),delim="\t")
F_data <- readr::read_delim(paste0(path,studies[[5]],".txt"),delim="\t")
H_data <- readr::read_delim(paste0(path,studies[[6]],".txt"),delim="\t")
I_data <- readr::read_delim(paste0(path,studies[[7]],".txt"),delim="\t")
#J_data <- readr::read_delim(paste0(path,studies[[8]],".txt"),delim="\t")
L1_data <- readr::read_delim(paste0(path,studies[[8]],".txt"),delim="\t")
L2_data <- readr::read_delim(paste0(path,studies[[9]],".txt"),delim="\t")


### Transform into tidy format
A_dat <- tidy_up(A_data,trim=F,indicator = "sp")
B_dat <- tidy_up(B_data,trim=T,indicator = "sp")
C_dat <- tidy_up(C_data,indicator = "sp")
E_dat <- tidy_up(E_data,trim=T,indicator = "sp")
F_dat <- tidy_up(F_data,indicator = "sp")
H_dat <- tidy_up(H_data,indicator = "sp")
I_dat <- tidy_up(I_data,indicator = "sp")  #%>% mutate(ID = str_replace_all(ID,"-","."))
#J_dat <- tidy_up(J_data)
L1_dat <- tidy_up(L1_data,indicator = "sp")
L2_dat <- tidy_up(L2_data,indicator = "sp")


### Prepare taxonomy based on V4-region trimmed sequences sequences
author.names <- c("Pietrucci et al., 2019",
                  "Qian et al., 2018",
                  "Keshavarzian et al., 2015",
                  "Aho et al., 2019",
                  "Heintz-Buschart et al., 2018",
                  "Weis et al., 2019",
                  "Barichella et al., 2019",
                  "Wallen et al., 2020a",
                  "Wallen et al., 2020b")

#project_path <- c("G:\\Projects\\ParkinsonPooledAnalysis\\")
project_path <- c("C:\\Users\\kleinebardenhorst\\Documents\\Projects\\ParkinsonPooledAnalysis\\")



## Prepare taxonomy based on same region data

path <- paste0(project_path,"Data\\")
seq_same <- paste("sequences_slicedV4\\")
tax_same <- paste("sequences_slicedV4_classified\\") 

taxa <- lapply(studies,function(x) read_tsv(paste0(path,tax_same,"Classified_Unaligned_Sliced_Sequences_",x,".align.fasta.txt"),quote="",col_names = F))
seq <- lapply(studies, function(x) ((seqinr::read.fasta(paste0(path,seq_same,"Unaligned_Sliced_Sequences_",x,".align.fasta"),as.string = T,set.attributes = F))))

seq <- foreach(i=seq_along(seq)) %do% data.frame("ID"=paste0("sp",sub("^[^_]*_", "", names(unlist(seq[[i]])))),"Sequence"=unlist(seq[[i]]))

taxon <- lapply(taxa,prepare_taxonomy_sameregion)  # Here, I trim the study indication - ASV1 in study A may not be ASV1 in study B. Identification between studies happen via exact sequence

A_tax <- taxon[[1]] %>% left_join(seq[[1]]) %>% apply(2,function(x) gsub('\"', "", x, fixed = TRUE)) %>% as_tibble() %>% arrange(ID)
B_tax <- taxon[[2]] %>% left_join(seq[[2]]) %>% apply(2,function(x) gsub('\"', "", x, fixed = TRUE)) %>% as_tibble() %>% arrange(ID)
C_tax <- taxon[[3]] %>% left_join(seq[[3]]) %>% apply(2,function(x) gsub('\"', "", x, fixed = TRUE)) %>% as_tibble() %>% arrange(ID)
E_tax <- taxon[[4]] %>% left_join(seq[[4]]) %>% apply(2,function(x) gsub('\"', "", x, fixed = TRUE)) %>% as_tibble() %>% arrange(ID)
F_tax <- taxon[[5]] %>% left_join(seq[[5]]) %>% apply(2,function(x) gsub('\"', "", x, fixed = TRUE)) %>% as_tibble() %>% arrange(ID)
H_tax <- taxon[[6]] %>% left_join(seq[[6]]) %>% apply(2,function(x) gsub('\"', "", x, fixed = TRUE)) %>% as_tibble() %>% arrange(ID)
I_tax <- taxon[[7]] %>% left_join(seq[[7]]) %>% apply(2,function(x) gsub('\"', "", x, fixed = TRUE)) %>% as_tibble() %>% arrange(ID)
L1_tax <- taxon[[8]] %>% left_join(seq[[8]]) %>% apply(2,function(x) gsub('\"', "", x, fixed = TRUE)) %>% as_tibble() %>% arrange(ID)
L2_tax <- taxon[[9]] %>% left_join(seq[[9]]) %>% apply(2,function(x) gsub('\"', "", x, fixed = TRUE)) %>% as_tibble() %>% arrange(ID)

## Prepare taxonomy based on full data
A_tax_FS <- prepare_taxonomy(A_data,indicator = "sp")
B_tax_FS <- prepare_taxonomy(B_data,indicator = "sp")
C_tax_FS <- prepare_taxonomy(C_data,indicator = "sp")
E_tax_FS <- prepare_taxonomy(E_data,indicator = "sp")
F_tax_FS <- prepare_taxonomy(F_data,indicator = "sp")
H_tax_FS <- prepare_taxonomy(H_data,indicator = "sp")
I_tax_FS <- prepare_taxonomy(I_data,indicator = "sp")
L1_tax_FS <- prepare_taxonomy(L1_data,indicator = "sp")
L2_tax_FS <- prepare_taxonomy(L2_data,indicator = "sp")


## Same region
A_clean <- clean(A_dat,A_tax,filter="rel",rel=0.005)
B_clean <- clean(B_dat,B_tax,filter="rel",rel=0.005)
C_clean <- clean(C_dat,C_tax,filter="rel",rel=0.005)
E_clean <- clean(E_dat,E_tax,filter="rel",rel=0.005)
F_clean <- clean(F_dat,F_tax,filter="rel",rel=0.005)
H_clean <- clean(H_dat,H_tax,filter="rel",rel=0.005)
I_clean <- clean(I_dat,I_tax,filter="rel",rel=0.005)
L1_clean <- clean(L1_dat,L1_tax,filter="rel",rel=0.005)
L2_clean <- clean(L2_dat,L2_tax,filter="rel",rel=0.005)

## Full sequences
A_clean_FS <- clean(A_dat,A_tax_FS,filter="rel",rel=0.005) 
B_clean_FS <- clean(B_dat,B_tax_FS,filter="rel",rel=0.005) 
C_clean_FS <- clean(C_dat,C_tax_FS,filter="rel",rel=0.005) 
E_clean_FS <- clean(E_dat,E_tax_FS,filter="rel",rel=0.005) 
F_clean_FS <- clean(F_dat,F_tax_FS,filter="rel",rel=0.005) 
H_clean_FS <- clean(H_dat,H_tax_FS,filter="rel",rel=0.005) 
I_clean_FS <- clean(I_dat,I_tax_FS,filter="rel",rel=0.005) 
L1_clean_FS <- clean(L1_dat,L1_tax_FS,filter="rel",rel=0.005)
L2_clean_FS <- clean(L2_dat,L2_tax_FS,filter="rel",rel=0.005)



## Read trees

trees <- lapply(studies,function(x) read_tree(paste0(path,"\\trees\\",x,".nwk")))

A_tree <- trees[[1]]
B_tree <- trees[[2]]
C_tree <- trees[[3]]
E_tree <- trees[[4]]
F_tree <- trees[[5]]
H_tree <- trees[[6]]
I_tree <- trees[[7]]
L1_tree <- trees[[8]]
L2_tree <- trees[[9]]

###########################
#### Prepare meta-data ####
###########################

# Load ENA metadata
A_Meta_ENA <- loadENAMeta(path="Data/meta_data/A_meta_ENA.txt") %>% dplyr::select(run_accession,library_name)
B_Meta_ENA <- loadENAMeta(path="Data/meta_data/B_meta_ENA.txt")
#C_Meta <- loadENAMeta(path="Data/C_meta.txt")
#D_Meta <- loadENAMeta(path="Data/D_meta.txt")
E_Meta_ENA <- loadENAMeta(path="Data/meta_data/E_meta_ENA.txt") %>% 
   dplyr::select(run_accession,sample_title) %>%
   dplyr::slice(1:266) %>% # Slice blank controls
   mutate(Timepoint = str_split(.$sample_title," ",simplify = TRUE)[,2],
          Subject = str_split(.$sample_title," ",simplify = TRUE)[,1]) %>%
   dplyr::select(-sample_title)

F_Meta_ENA <- loadENAMeta(path="Data/meta_data/F_meta_ENA.txt") %>% dplyr::select(run_accession, library_name,experiment_alias)
G_Meta_ENA <- loadENAMeta(path="Data/meta_data/G_meta_ENA.txt")
H_Meta_ENA <- loadENAMeta(path="Data/meta_data/H_meta_ENA.txt") %>%
   dplyr::select(run_accession,sample_alias) %>%
   dplyr::slice(1:256) %>% # Slice ion torrent
   mutate(ID = run_accession,
          Sample = str_split(.$sample_alias,"_",simplify = T)[,1]) %>%
   dplyr::select(ID,Sample)


#I_Meta <- loadENAMeta(path="Data/I_meta.txt")
J_Meta_ENA <- loadENAMeta(path="Data/meta_data/J_meta_ENA.txt")

# Load detailed metadata

C_meta1 <- xlsx::read.xlsx(paste0(project_path,"Data/meta_data/C_meta.xlsx"),sheetIndex = 1, header = T,startRow =2) %>%
   mutate(gender = Gender,
          age = Age,.keep = "unused")
names <- str_remove(colnames(C_meta1),"PD.")
colnames(C_meta1) <- names

C_meta2 <- xlsx::read.xlsx(paste0(project_path,"Data/meta_data/C_meta.xlsx"),sheetIndex = 2, header = T,startRow =2) %>%
   mutate(gender = Gender,
          age = Age,.keep = "unused")
names <- str_remove(colnames(C_meta2),"HC.")
colnames(C_meta2) <- names

C_Meta <- C_meta1 %>% full_join(C_meta2) %>% filter(Feces %in% toupper(C_dat$ID) | Tissue %in% toupper(C_dat$ID))

E_meta <- readr::read_csv(paste0(project_path,"Data/meta_data/E_meta.csv")) %>% 
   left_join(E_Meta_ENA) %>% dplyr::select(run_accession,everything()) %>%
   mutate(age = age_at_stool_collection,
          gender = ifelse(gender=="M","Male","Female"),.keep = "unused")

H_meta <- xlsx::read.xlsx(paste0(project_path,"Data/meta_data/H_meta.xlsx"),
                          sheetIndex = 1, header = T,startRow =1)  %>% # Other drugs variable not used and deleted for easier import
   mutate(gender = ifelse(Sex=="m","Male","Female"),
          age = Age, .keep="unused")


I_meta <- xlsx::read.xlsx(paste0(project_path,"Data/meta_data/I_meta.xls"), sheetIndex = 1,header=T,startRow = 2) %>% dplyr::slice(1:306) %>%
   mutate(gender = ifelse(gender==1,"Male","Female"),.keep="unused")

J_meta <- readr::read_csv(paste0(project_path,"Data/meta_data/J_meta.csv")) %>%
   dplyr::select(cc,ID,Age,Sex,Batch,CC,Onset,'UPDRS-I','UPDRS-III',Stool_freq,Constipation,Vegetarian,Nicotine,Alcohol,Caffeine,Ferm_Milk,Age_DBS,L_Dopa,Benserazid,Carbidopa,Dopa_Agon,COMT_Inhibitor,MAO_B_Inh,CV_disease) %>%
   mutate(age = Age,.keep="unused")

## Study L
library(XML)

L_xml <- xmlParse(file = "Data/meta_Data/L_meta_ENA.xml") ## Load XML file
L_meta_ENA <- read.delim("Data/meta_data/L_meta_ENA.txt") %>%
   dplyr::select(study_accession,sample_accession,experiment_accession,run_accession,sample_alias)

L_xml <- xmlToList(L_xml) ## Convert XML to nested list

names(L_xml) <- paste0(names(L_xml),"_",1:840) ## Create new biosample indicator



## Rename list elements of attributes by attribute name
for (i in seq_along(L_xml)) {
   
   for (j in seq_along(L_xml[[i]][[6]])) {
      names(L_xml[[i]][[6]])[[j]] <- L_xml[[i]][[6]][[j]][[2]][[1]]
      L_xml[[i]][[6]][[j]] <- L_xml[[i]][[6]][[j]][[1]]
   }
}

## Bring attributes in nice format
out <- list()
ids <- list()
for(i in seq_along(L_xml)) {
   out[[i]] <- unlist(L_xml[[i]]$Attributes)
   
}

## Save attributes as data.frame
L_meta <- as.data.frame(purrr::reduce(out,rbind)) %>%
   dplyr::mutate(sample_alias = isolate,.keep="unused",.before="host") %>%
   dplyr::left_join(L_meta_ENA,by="sample_alias") %>%
   mutate(status = ifelse(health_state=="Control","HC","PD"),
          gender = ifelse(host_sex=="male","Male","Female"),
          age = host_age,
          .keep="unused")

table(L_meta$dataset)

L1_meta <- L_meta %>% filter(dataset == "dataset_2")
L2_meta <- L_meta %>% filter(dataset == "dataset_1")

##################################
#### Prepare phyloseq objects ####
##################################

# Combine Taxonomies, Metadata and ASV-Data into phyloseq objects

# Study A


## OTU
OTU_A <- as.matrix(A_clean[,-1])
rownames(OTU_A) <- A_dat$ID
OTU_A <- otu_table(t(OTU_A),taxa_are_rows = T)

## Tax
tax_A <- as.matrix(A_tax[,-1])
rownames(tax_A) <- A_tax$ID

tax_A <- tax_table(tax_A)

## Meta
A_sample <- A_Meta_ENA %>%
   mutate(pd = substr(.$library_name,1,3),
          author = "Pietrucci et al., 2019") %>%
   mutate(status = ifelse(pd=="CON","HC","PD"),
          region = "V3-V4",
          location = "Italy") %>%
   column_to_rownames(var="run_accession") %>% sample_data()

## Tree

## Combine
A_phyl <- phyloseq(OTU_A,tax_A,A_sample,A_tree)


## Study A - Original
## OTU
OTU_A_FS <- as.matrix(A_clean_FS[,-1])
rownames(OTU_A_FS) <- A_dat$ID
OTU_A_FS <- otu_table(t(OTU_A_FS),taxa_are_rows = T)

## Tax
tax_A_FS <- as.matrix(A_tax_FS[,-1])
rownames(tax_A_FS) <- A_tax_FS$ID

tax_A_FS <- tax_table(tax_A_FS)

## Combine

A_phyl_FS <- phyloseq(OTU_A_FS,tax_A_FS,A_sample)


# Study B

## OTU
OTU_B <- as.matrix(B_clean[,-1])
rownames(OTU_B) <- B_dat$ID
OTU_B <- otu_table(t(OTU_B),taxa_are_rows = T)

## Tax
tax_B <- as.matrix(B_tax[,-1])
rownames(tax_B) <- B_tax$ID
tax_B <- tax_table(tax_B)

## Meta
B_sample <- B_Meta_ENA %>%
   dplyr::select(run_accession,library_name) %>%
   mutate(status = case_when(startsWith(library_name,"GC") ~ "HC",
                             startsWith(library_name,"GP") ~ "PD"),
          author = "Qian et al., 2018",
          region = "V3-V4",
          location = "China") %>%
   filter(not.na(status)) %>% # Filter out blood metagenome
   column_to_rownames(var="run_accession") %>%
   sample_data()

## Combine
B_phyl <- phyloseq(OTU_B,tax_B,B_sample,B_tree)


## Study B - Full sequence
## OTU
OTU_B_FS <- as.matrix(B_clean_FS[,-1])
rownames(OTU_B_FS) <- B_dat$ID
OTU_B_FS <- otu_table(t(OTU_B_FS),taxa_are_rows = T)

## Tax
tax_B_FS <- as.matrix(B_tax_FS[,-1])
rownames(tax_B_FS) <- B_tax_FS$ID
tax_B_FS <- tax_table(tax_B_FS)


## Combine
B_phyl_FS <- phyloseq(OTU_B_FS,tax_B_FS,B_sample)


# Study C
## OTU
OTU_C <- as.matrix(C_clean[,-1])
rownames(OTU_C) <- toupper(C_dat$ID)
OTU_C <- otu_table(t(OTU_C),taxa_are_rows = T)

## Tax
tax_C <- as.matrix(C_tax[,-1])
#rownames(tax_C) <- toupper(C_tax$ID)

tax_C <- tax_table(tax_C)

## Meta
C_sample <- C_Meta %>%
   dplyr::select(Case.IDs,Tissue,Feces,everything()) %>%
   tidyr::gather(Sample,ID,2:3) %>%
   dplyr::select(Case.IDs,Sample,ID,everything()) %>%
   filter(ID != 'NA') %>%
   filter(Sample=="Feces") %>%
   mutate(status = ifelse(startsWith(Cases,"Parkinson"),"PD","HC"),
          author = "Keshavarzian et al., 2015",
          region = "V4",
          location = "USA") %>%
   dplyr::select(status,everything()) %>%
   column_to_rownames(var='ID') %>%
   sample_data

## Combine
C_phyl <- phyloseq(OTU_C,tax_C,C_sample,C_tree)


## Study C - Full sequence
## OTU
OTU_C_FS <- as.matrix(C_clean_FS[,-1])
rownames(OTU_C_FS) <- C_dat$ID
OTU_C_FS <- otu_table(t(OTU_C_FS),taxa_are_rows = T)

## Tax
tax_C_FS <- as.matrix(C_tax_FS[,-1])
rownames(tax_C_FS) <- C_tax_FS$ID
tax_C_FS <- tax_table(tax_C_FS)


## Combine
C_phyl_FS <- phyloseq(OTU_C_FS,tax_C_FS,C_sample)



# Study E
## OTU
OTU_E <- as.matrix(E_clean[,-1])
rownames(OTU_E) <- E_dat$ID
OTU_E <- otu_table(t(OTU_E),taxa_are_rows = T)

## Tax
tax_E <- as.matrix(E_tax[,-1])
rownames(tax_E) <- E_tax$ID

tax_E <- tax_table(tax_E)

## Meta
E_sample <- E_meta %>%
   filter(Timepoint == "followup") %>%
   mutate(status = ifelse(Parkinson=="control","HC","PD"),
          author = "Aho et al., 2019",
          region = "V3-V4",
          location = "Finland") %>%
   dplyr::select(status,everything()) %>%
   column_to_rownames(var='run_accession') %>%
   sample_data()

## Combine
E_phyl <- phyloseq(OTU_E,tax_E,E_sample,E_tree)


## Study E - Full sequence
## OTU
OTU_E_FS <- as.matrix(E_clean_FS[,-1])
rownames(OTU_E_FS) <- E_dat$ID
OTU_E_FS <- otu_table(t(OTU_E_FS),taxa_are_rows = T)

## Tax
tax_E_FS <- as.matrix(E_tax_FS[,-1])
rownames(tax_E_FS) <- E_tax_FS$ID
tax_E_FS <- tax_table(tax_E_FS)


## Combine
E_phyl_FS <- phyloseq(OTU_E_FS,tax_E_FS,E_sample)



# Study F
## OTU
OTU_F <- as.matrix(F_clean[,-1])
rownames(OTU_F) <- F_dat$ID
OTU_F <- otu_table(t(OTU_F),taxa_are_rows = T)

## Tax
tax_F <- as.matrix(F_tax[,-1])
rownames(tax_F) <- F_tax$ID

tax_F <- tax_table(tax_F)

## Meta
F_sample <- F_Meta_ENA %>%
   filter(endsWith(library_name,".f")) %>%
   filter(startsWith(library_name,"PD")|startsWith(library_name,"HC")) %>%
   mutate(status = ifelse(startsWith(library_name,"PD"),"PD","HC"),
          author = "Heintz-Buschart et al., 2018",
          region = "V4",
          location = "Germany") %>%
   column_to_rownames(var="run_accession") %>% sample_data


## Combine
F_phyl <- phyloseq(OTU_F,tax_F,F_sample,F_tree)


## Study F - Full sequence
## OTU
OTU_F_FS <- as.matrix(F_clean_FS[,-1])
rownames(OTU_F_FS) <- F_dat$ID
OTU_F_FS <- otu_table(t(OTU_F_FS),taxa_are_rows = T)

## Tax
tax_F_FS <- as.matrix(F_tax_FS[,-1])
rownames(tax_F_FS) <- F_tax_FS$ID
tax_F_FS <- tax_table(tax_F_FS)


## Combine
F_phyl_FS <- phyloseq(OTU_F_FS,tax_F_FS,F_sample)


# Study H
## OTU
OTU_H <- H_clean %>%
   left_join(H_Meta_ENA) %>%
   group_by(Sample) %>%
   summarize_if(is.numeric,funs(sum)) %>%
   arrange(desc(Sample))

id <- substr(OTU_H$Sample,4,8)
id[40:length(id)] <- paste0("P",id[40:length(id)])

OTU_H$Sample[1:39] <-  paste0('IfM ',substr(id[1:39],1,1),"0",substr(id[1:39],2,4))
OTU_H$Sample[40:nrow(OTU_H)] <- paste0('IfM ',substr(id[40:length(id)],1,2),"0",substr(id[40:length(id)],3,5))

OTU_H <- OTU_H %>% column_to_rownames(var="Sample") %>% as.matrix %>% t %>% otu_table(taxa_are_rows = T)
## Tax
tax_H <- as.matrix(H_tax[,-1])
rownames(tax_H) <- H_tax$ID
tax_H <- tax_table(tax_H)

## Meta
H_sample <- H_meta %>%
   mutate(status = ifelse(Group =="PD","PD","HC"),
          author = "Weis et al., 2019",
          region = "V4-V5",
          location = "Germany") %>%
   dplyr::select(status,everything(),-Group) %>%
   column_to_rownames(var="SampleID") %>%
   sample_data

## Combine
H_phyl <- phyloseq(OTU_H,tax_H,H_sample,H_tree)


## Study H - Full sequence
## OTU
OTU_H_FS <- H_clean_FS %>%
   left_join(H_Meta_ENA) %>%
   group_by(Sample) %>%
   summarize_if(is.numeric,funs(sum)) %>%
   arrange(desc(Sample))

id <- substr(OTU_H_FS$Sample,4,8)
id[40:length(id)] <- paste0("P",id[40:length(id)])

OTU_H_FS$Sample[1:39] <-  paste0('IfM ',substr(id[1:39],1,1),"0",substr(id[1:39],2,4))
OTU_H_FS$Sample[40:nrow(OTU_H_FS)] <- paste0('IfM ',substr(id[40:length(id)],1,2),"0",substr(id[40:length(id)],3,5))

OTU_H_FS <- OTU_H_FS %>% column_to_rownames(var="Sample") %>% as.matrix %>% t %>% otu_table(taxa_are_rows = T)

## Tax
tax_H_FS <- as.matrix(H_tax_FS[,-1])
rownames(tax_H_FS) <- H_tax_FS$ID
tax_H_FS <- tax_table(tax_H_FS)


## Combine
H_phyl_FS <- phyloseq(OTU_H_FS,tax_H_FS,H_sample)


# Study I

## OTU
OTU_I <- as.matrix(I_clean[,-1])
rownames(OTU_I) <- str_replace_all(I_dat$ID,"-",".")
OTU_I <- otu_table(t(OTU_I),taxa_are_rows = T)

## Tax
tax_I <- as.matrix(I_tax[,-1])
rownames(tax_I) <- I_tax$ID

tax_I <- tax_table(tax_I)

## Meta
I_sample <- I_meta %>%
   mutate(status = ifelse(status==0,"HC","PD"),
          author = "Barichella et al., 2019",
          region = "V3-V4",
          location = "Italy")  %>%
   column_to_rownames(var="ID") %>%
   sample_data


## Combine
I_phyl <- phyloseq(OTU_I,tax_I,I_sample,I_tree)


## Study I - Full sequence
## OTU
OTU_I_FS <- as.matrix(I_clean_FS[,-1])
rownames(OTU_I_FS) <- str_replace_all(I_dat$ID,"-",".")
OTU_I_FS <- otu_table(t(OTU_I_FS),taxa_are_rows = T)

## Tax
tax_I_FS <- as.matrix(I_tax_FS[,-1])
rownames(tax_I_FS) <- I_tax_FS$ID
tax_I_FS <- tax_table(tax_I_FS)


## Combine
I_phyl_FS <- phyloseq(OTU_I_FS,tax_I_FS,I_sample)


## Study L1

## Trimmed seq
## OTU
OTU_L1 <- as.matrix(L1_clean[,-1])
rownames(OTU_L1) <- str_replace_all(L1_dat$ID,"-",".")
OTU_L1 <- otu_table(t(OTU_L1),taxa_are_rows = T)

## Tax
tax_L1 <- as.matrix(L1_tax[,-1])
rownames(tax_L1) <- L1_tax$ID
tax_L1 <- tax_table(tax_L1)

# Sample
L1_sample <- L1_meta %>%
   column_to_rownames("run_accession") %>%
   mutate(author = "Wallen et al., 2020a",
          region = "V4",
          location = "USA") %>%
   sample_data

# Combine
L1_phyl <- phyloseq(OTU_L1,tax_L1,L1_sample,L1_tree)

## Full seq
## OTU
OTU_L1_FS <- as.matrix(L1_clean_FS[,-1])
rownames(OTU_L1_FS) <- str_replace_all(L1_dat$ID,"-",".")
OTU_L1_FS <- otu_table(t(OTU_L1_FS),taxa_are_rows = T)

## Tax
tax_L1_FS <- as.matrix(L1_tax_FS[,-1])
rownames(tax_L1_FS) <- L1_tax_FS$ID
tax_L1_FS <- tax_table(tax_L1_FS)

L1_phyl_FS <- phyloseq(OTU_L1_FS,tax_L1_FS,L1_sample)


## Study L2

## Trimmed seq
## OTU
OTU_L2 <- as.matrix(L2_clean[,-1])
rownames(OTU_L2) <- str_replace_all(L2_dat$ID,"-",".")
OTU_L2 <- otu_table(t(OTU_L2),taxa_are_rows = T)

## Tax
tax_L2 <- as.matrix(L2_tax[,-1])
rownames(tax_L2) <- L2_tax$ID
tax_L2 <- tax_table(tax_L2)

# Sample
L2_sample <- L2_meta %>%
   column_to_rownames("run_accession") %>%
   mutate(author = "Wallen et al., 2020b",
          region = "V4",
          location = "USA") %>%
   sample_data


# Combine
L2_phyl <- phyloseq(OTU_L2,tax_L2,L2_sample,L2_tree)

## Full seq
## OTU
OTU_L2_FS <- as.matrix(L2_clean_FS[,-1])
rownames(OTU_L2_FS) <- str_replace_all(L2_dat$ID,"-",".")
OTU_L2_FS <- otu_table(t(OTU_L2_FS),taxa_are_rows = T)

## Tax
tax_L2_FS <- as.matrix(L2_tax_FS[,-1])
rownames(tax_L2_FS) <- L2_tax_FS$ID
tax_L2_FS <- tax_table(tax_L2_FS)

L2_phyl_FS <- phyloseq(OTU_L2_FS,tax_L2_FS,L2_sample)





# Same region (V4) sequences
A_phyl
B_phyl
C_phyl
E_phyl
F_phyl
H_phyl
I_phyl
L1_phyl
L2_phyl

# Full sequences
A_phyl_FS
B_phyl_FS
C_phyl_FS
E_phyl_FS
F_phyl_FS
H_phyl_FS
I_phyl_FS
L1_phyl_FS
L2_phyl_FS


##############################
#### Filtering and checks ####
##############################

# Filter samples below 10.000 reads

A_phyl <- prune_samples(sample_sums(A_phyl)>=10000, A_phyl)
B_phyl <- prune_samples(sample_sums(B_phyl)>=10000, B_phyl) ## 20 Samples filtered
C_phyl <- prune_samples(sample_sums(C_phyl)>=10000, C_phyl) ## 2 samples filtered
E_phyl <- prune_samples(sample_sums(E_phyl)>=10000, E_phyl) ## 2 Samples filtered
F_phyl <- prune_samples(sample_sums(F_phyl)>=10000, F_phyl) ## 1 sample filtered
H_phyl <- prune_samples(sample_sums(H_phyl)>=10000, H_phyl)
I_phyl <- prune_samples(sample_sums(I_phyl)>=10000, I_phyl) ## 47 samples filtered
L1_phyl <- prune_samples(sample_sums(L1_phyl)>=10000, L1_phyl) ## 0 samples filtered
L2_phyl <- prune_samples(sample_sums(L2_phyl)>=10000, L2_phyl) ## 186 samples filtered



A_phyl_FS <- prune_samples(sample_sums(A_phyl_FS)>=10000, A_phyl_FS)
B_phyl_FS <- prune_samples(sample_sums(B_phyl_FS)>=10000, B_phyl_FS) # 20 Samples filtered
C_phyl_FS <- prune_samples(sample_sums(C_phyl_FS)>=10000, C_phyl_FS) # 2 samples filtered
E_phyl_FS <- prune_samples(sample_sums(E_phyl_FS)>=10000, E_phyl_FS) # 2 samples filtered
F_phyl_FS <- prune_samples(sample_sums(F_phyl_FS)>=10000, F_phyl_FS) # 1 sample filtered
H_phyl_FS <- prune_samples(sample_sums(H_phyl_FS)>=10000, H_phyl_FS)
I_phyl_FS <- prune_samples(sample_sums(I_phyl_FS)>=10000, I_phyl_FS) # 62 samples
L1_phyl_FS <- prune_samples(sample_sums(L1_phyl_FS)>=10000, L1_phyl_FS) ## 0 samples filtered
L2_phyl_FS <- prune_samples(sample_sums(L2_phyl_FS)>=10000, L2_phyl_FS) ## 196 samples filtered


#### Save data ####
dat <- list(A_phyl,
            B_phyl,
            C_phyl,
            E_phyl,
            F_phyl,
            H_phyl,
            I_phyl,
            L1_phyl,
            L2_phyl)

#saveRDS(dat,"phyl_list_unfiltered.RDS") # ran all the code without filtering steps to produce this unfiltered dataset
saveRDS(dat,"Data/processed_data/phyl_list.RDS")

dat <- readRDS("Data/processed_data/phyl_list.RDS")

dat_FS <- list(A_phyl_FS,
               B_phyl_FS,
               C_phyl_FS,
               E_phyl_FS,
               F_phyl_FS,
               H_phyl_FS,
               I_phyl_FS,
               L1_phyl_FS,
               L2_phyl_FS)

saveRDS(dat_FS,"Data/processed_data/phyl_list_FS.RDS")


