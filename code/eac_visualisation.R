
# Quick visualisation of results
# for viral diversity and transmission clusters in EAC


set.seed(777)

library(ape)
library(seqinr)
library(stringi)
library(lubridate)
library(tidyverse)
library(ggtree)
library(igraph)


# Load sequence data -----------


# Burundi

A_BI <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/A_Burundi_protease.Fasta")
C_BI <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/C_Burundi_protease.Fasta")
D_BI <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/D_Burundi_protease.Fasta")

A_BI_seqs <- names(A_BI) # 25
C_BI_seqs <- names(C_BI) # 289
D_BI_seqs <- names(D_BI) # 5

A_BI_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.BI.", "", A_BI_seqs), "\\d{4}"), format = "%Y"))# stringi # lubricate
C_BI_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.BI.", "", C_BI_seqs), "\\d{4}"), format = "%Y"))
D_BI_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.BI.", "", D_BI_seqs), "\\d{4}"), format = "%Y"))



# DRC

A_DRC <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/A_DRC_protease.Fasta")
C_DRC <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/C_DRC_protease.Fasta")
D_DRC <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/D_DRC_protease.Fasta")

A_DRC_seqs <- names(A_DRC) #  43
C_DRC_seqs <- names(C_DRC) #  22
D_DRC_seqs <- names(D_DRC) # 28

A_DRC_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.DC.", "", A_DRC_seqs), "\\d{4}"), format = "%Y"))# stringi # lubricate
C_DRC_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.CD.", "", C_DRC_seqs), "\\d{4}"), format = "%Y"))
D_DRC_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.CD.", "", D_DRC_seqs), "\\d{4}"), format = "%Y"))


# Kenya

A_KE <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/A_Kenya_protease.Fasta")
C_KE <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/C_Kenya_protease.Fasta")
D_KE <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/D_Kenya_protease.Fasta")

A_KE_seqs <- names(A_KE) # 1295
C_KE_seqs <- names(C_KE) # 104
D_KE_seqs <- names(D_KE) # 145

A_KE_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.DC.", "", A_KE_seqs), "\\d{4}"), format = "%Y"))# stringi # lubricate
C_KE_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.CD.", "", C_KE_seqs), "\\d{4}"), format = "%Y"))
D_KE_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.CD.", "", D_KE_seqs), "\\d{4}"), format = "%Y"))


# Uganda

A_UGA <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/A_Uganda_protease.Fasta")
C_UGA <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/C_Uganda_protease.Fasta")
D_UGA <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/D_Uganda_protease.Fasta")


A_UGA_seqs <- names(A_UGA) # 1180
C_UGA_seqs <- names(C_UGA) # 66
D_UGA_seqs <- names(D_UGA) # 1099

A_UGA_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.DC.", "", A_UGA_seqs), "\\d{4}"), format = "%Y"))# stringi # lubricate
C_UGA_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.CD.", "", C_UGA_seqs), "\\d{4}"), format = "%Y"))
D_UGA_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.CD.", "", D_UGA_seqs), "\\d{4}"), format = "%Y"))


# Tanzania

A_TZ <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/A_Tanzania_protease.Fasta")
C_TZ <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/C_Tanzania_protease.Fasta")
D_TZ <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/D_Tanzania_protease.Fasta")


A_TZ_seqs <- names(A_TZ) # 363
C_TZ_seqs <- names(C_TZ) # 388
D_TZ_seqs <- names(D_TZ) # 163

A_TZ_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.DC.", "", A_TZ_seqs), "\\d{4}"), format = "%Y"))# stringi # lubricate
C_TZ_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.CD.", "", C_TZ_seqs), "\\d{4}"), format = "%Y"))
D_TZ_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.CD.", "", D_TZ_seqs), "\\d{4}"), format = "%Y"))



# Rwanda

A_RW <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/A_Rwanda_protease.Fasta")
C_RW <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/C_Rwanda_protease.Fasta")
D_RW <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/D_Rwanda_protease.Fasta")


A_RW_seqs <- names(A_RW) # 189
C_RW_seqs <- names(C_RW) # 24
D_RW_seqs <- names(D_RW) # 5

A_RW_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.DC.", "", A_RW_seqs), "\\d{4}"), format = "%Y"))# stringi # lubricate
C_RW_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.CD.", "", C_RW_seqs), "\\d{4}"), format = "%Y"))
D_RW_seqs_dates <- year(as.Date(stri_extract_first_regex(gsub("A1.CD.", "", D_RW_seqs), "\\d{4}"), format = "%Y"))


# Tables for sequences numbers and sampling dates ------------


# Accession numbers 

access_number_eac_seqs <- c(A_BI_seqs, C_BI_seqs, D_BI_seqs, 
                            A_DRC_seqs, C_DRC_seqs, D_DRC_seqs, 
                            A_KE_seqs, C_KE_seqs, D_KE_seqs, 
                            A_UGA_seqs, C_UGA_seqs, D_UGA_seqs, 
                            A_TZ_seqs, C_TZ_seqs, D_TZ_seqs, 
                            A_RW_seqs, C_RW_seqs, D_RW_seqs)

write(access_number_eac_seqs, file = "data/access_number_seqs_eac.txt")



# Table of sequences per subtype and country

A <- c(length(A_BI_seqs), length(A_DRC_seqs), length(A_KE_seqs), length(A_UGA_seqs), length(A_TZ_seqs), length(A_RW_seqs))
C <- c(length(C_BI_seqs), length(C_DRC_seqs), length(C_KE_seqs), length(C_UGA_seqs), length(C_TZ_seqs), length(C_RW_seqs))
D <- c(length(D_BI_seqs), length(D_DRC_seqs), length(D_KE_seqs), length(D_UGA_seqs), length(D_TZ_seqs), length(D_RW_seqs))

countries <- c("Burundi", "Congo", "Kenya", "Uganda", "Tanzania", "Rwanda")

count_sequences_ACD <- data.frame(countries, A, C, D)

names(count_sequences_ACD) <- c("Country", "Subtype A", "Subtype C", "Subtype D")


write.csv(count_sequences_ACD, file = "results/table_count_sequences_ACD.csv")



# Min and Max sampling dates per subtype per country -

min_max <- function(vec){
  min.i <- min(vec)
  max.i <- max(vec)
  i <- paste0(min.i,"_", max.i)
  i <- paste(i)
  return(i)
}


A_years <- c(paste(min_max(A_BI_seqs_dates)), paste(min_max(A_DRC_seqs_dates)), 
             paste(min_max(A_KE_seqs_dates)), paste(min_max(A_UGA_seqs_dates)), 
             paste(min_max(A_TZ_seqs_dates)), paste(min_max(A_RW_seqs_dates)))


C_years <- c(paste(min_max(C_BI_seqs_dates)), paste(min_max(C_DRC_seqs_dates)), 
             paste(min_max(C_KE_seqs_dates)), paste(min_max(C_UGA_seqs_dates)), 
             paste(min_max(C_TZ_seqs_dates)), paste(min_max(C_RW_seqs_dates)))


D_years <- c(paste(min_max(D_BI_seqs_dates)), paste(min_max(D_DRC_seqs_dates)), 
             paste(min_max(D_KE_seqs_dates)), paste(min_max(D_UGA_seqs_dates)), 
             paste(min_max(D_TZ_seqs_dates)), paste(min_max(D_RW_seqs_dates)))



countries <- c("Burundi", "Congo", "Kenya", "Uganda", "Tanzania", "Rwanda")

sampling_years_ACD <- data.frame(countries, A_years, C_years, D_years)

names(sampling_years_ACD) <- c("Country", "Subtype A", "Subtype C", "Subtype D")


write.csv(sampling_years_ACD, file = "results/table_sampling_years_ACD.csv")




# Trend of internal nodes per subtype and country --------------

# Read tables of internal nodes per year per subtype and country

# Raw numbers for internal nodes 

# Per country we compute yearly proportion of each subtype:

# prop_sub_A = n_sub_A/(n_sub_A + n_sub_C + n_sub_D)
# prop_sub_C = n_sub_A/(n_sub_A + n_sub_C + n_sub_D)
# prop_sub_D = n_sub_A/(n_sub_A + n_sub_C + n_sub_D)

A_BI_nodes.long.df_raw <- readRDS("results/A_BI_nodes.long.df_raw.RDS")
C_BI_nodes.long.df_raw <- readRDS("results/C_BI_nodes.long.df_raw.RDS")
D_BI_nodes.long.df_raw <- readRDS("results/D_BI_nodes.long.df_raw.RDS")

BI_df_raw <- rbind(A_BI_nodes.long.df_raw, C_BI_nodes.long.df_raw, D_BI_nodes.long.df_raw)
BI_df_raw$prop.intern.nodes <- BI_df_raw$intern.nodes/sum(BI_df_raw$intern.nodes)
BI_type_raw <- c(rep("A", nrow(A_BI_nodes.long.df_raw)), rep("C", nrow(C_BI_nodes.long.df_raw)), rep("D", nrow(D_BI_nodes.long.df_raw)))

A_DRC_nodes.long.df_raw <- readRDS("results/A_DRC_nodes.long.df_raw.RDS")
C_DRC_nodes.long.df_raw <- readRDS("results/C_DRC_nodes.long.df_raw.RDS")
D_DRC_nodes.long.df_raw <- readRDS("results/D_DRC_nodes.long.df_raw.RDS")

DRC_df_raw <- rbind(A_DRC_nodes.long.df_raw, C_DRC_nodes.long.df_raw, D_DRC_nodes.long.df_raw)
DRC_df_raw$prop.intern.nodes <- DRC_df_raw$intern.nodes/sum(DRC_df_raw$intern.nodes)
DRC_type_raw <- c(rep("A", nrow(A_DRC_nodes.long.df_raw)), rep("C", nrow(C_DRC_nodes.long.df_raw)), rep("D", nrow(D_DRC_nodes.long.df_raw)))


A_KE_nodes.long.df_raw <- readRDS("results/A_KE_nodes.long.df_raw.RDS")
C_KE_nodes.long.df_raw <- readRDS("results/C_KE_nodes.long.df_raw.RDS")
D_KE_nodes.long.df_raw <- readRDS("results/D_KE_nodes.long.df_raw.RDS")

KE_df_raw <- rbind(A_KE_nodes.long.df_raw, C_KE_nodes.long.df_raw, D_KE_nodes.long.df_raw)
KE_df_raw$prop.intern.nodes <- KE_df_raw$intern.nodes/sum(KE_df_raw$intern.nodes)
KE_type_raw <- c(rep("A", nrow(A_KE_nodes.long.df_raw)), rep("C", nrow(C_KE_nodes.long.df_raw)), rep("D", nrow(D_KE_nodes.long.df_raw)))


A_UGA_nodes.long.df_raw <- readRDS("results/A_UGA_nodes.long.df_raw.RDS")
C_UGA_nodes.long.df_raw <- readRDS("results/C_UGA_nodes.long.df_raw.RDS")
D_UGA_nodes.long.df_raw <- readRDS("results/D_UGA_nodes.long.df_raw.RDS")

UGA_df_raw <- rbind(A_UGA_nodes.long.df_raw, C_UGA_nodes.long.df_raw, D_UGA_nodes.long.df_raw)
UGA_df_raw$prop.intern.nodes <- UGA_df_raw$intern.nodes/sum(UGA_df_raw$intern.nodes)
UGA_type_raw <- c(rep("A", nrow(A_UGA_nodes.long.df_raw)), rep("C", nrow(C_UGA_nodes.long.df_raw)), rep("D", nrow(D_UGA_nodes.long.df_raw)))


A_TZ_nodes.long.df_raw <- readRDS("results/A_TZ_nodes.long.df_raw.RDS")
C_TZ_nodes.long.df_raw <- readRDS("results/C_TZ_nodes.long.df_raw.RDS")
D_TZ_nodes.long.df_raw <- readRDS("results/D_TZ_nodes.long.df_raw.RDS")

TZ_df_raw <- rbind(A_TZ_nodes.long.df_raw, C_TZ_nodes.long.df_raw, D_TZ_nodes.long.df_raw)
TZ_df_raw$prop.intern.nodes <- TZ_df_raw$intern.nodes/sum(TZ_df_raw$intern.nodes)
TZ_type_raw <- c(rep("A", nrow(A_TZ_nodes.long.df_raw)), rep("C", nrow(C_TZ_nodes.long.df_raw)), rep("D", nrow(D_TZ_nodes.long.df_raw)))


A_RW_nodes.long.df_raw <- readRDS("results/A_RW_nodes.long.df_raw.RDS")
C_RW_nodes.long.df_raw <- readRDS("results/C_RW_nodes.long.df_raw.RDS")
D_RW_nodes.long.df_raw <- readRDS("results/D_RW_nodes.long.df_raw.RDS")

RW_df_raw <- rbind(A_RW_nodes.long.df_raw, C_RW_nodes.long.df_raw, D_RW_nodes.long.df_raw)
RW_df_raw$prop.intern.nodes <- RW_df_raw$intern.nodes/sum(RW_df_raw$intern.nodes)
RW_type_raw <- c(rep("A", nrow(A_RW_nodes.long.df_raw)), rep("C", nrow(C_RW_nodes.long.df_raw)), rep("D", nrow(D_RW_nodes.long.df_raw)))


name_country_raw <- c(rep("BI", nrow(BI_df_raw)), rep("DRC", nrow(DRC_df_raw)),
                      rep("KE", nrow(KE_df_raw)), rep("UGA", nrow(UGA_df_raw)),
                      rep("TZ", nrow(TZ_df_raw)), rep("RW", nrow(RW_df_raw)))

HIV_Subtype_raw <- c(BI_type_raw, DRC_type_raw, KE_type_raw, UGA_type_raw, TZ_type_raw, RW_type_raw)

EAC_df_raw <- rbind(BI_df_raw, DRC_df_raw, KE_df_raw, UGA_df_raw, TZ_df_raw, RW_df_raw)

EAC_df_raw$country <- name_country_raw

EAC_df_raw$HIV_Subtype <- HIV_Subtype_raw



plot.EAC_raw_prop_df_col <- ggplot(EAC_df_raw, aes(x=calendaryear, y=prop.intern.nodes, colour=HIV_Subtype, group = HIV_Subtype)) + 
  # geom_errorbar(aes(ymin=lower.Q1, ymax=upper.Q3), width=.1) +
  geom_line(size=0.5) +
  geom_point(aes(shape=HIV_Subtype), size=1) + 
  facet_wrap(vars(country), nrow = 6)+
  xlab("Year") + ylab("Proportion of internal node")+
  theme(legend.position="bottom")



ggsave(filename = "plot.EAC_raw_prop_df_col.JPEG",
       plot = plot.EAC_raw_prop_df_col,
       path = "/home/david/Dropbox/eac_phylo_analysis/results",
       width = 17, height = 20, units = "cm")





# Build a transmission network from transmission clusters per subtype in EAC --------


# Read the clusters' files for subtype A

sub.dir.rename <- "/home/david/Dropbox/eac_phylo_analysis/results/Transmission Clusters/"

dd <- list.files(path = paste0(sub.dir.rename), 
                 pattern = paste0("EAC_HIV_1_A_EAC_HIV_1_A_clusterPicks_cluster"),
                 all.files = FALSE,
                 full.names = FALSE, recursive = FALSE)

# Transmission clusters.

d <- dd

clust.size_a <- vector() # size of each cluster # table.simpact.trans.net.adv

name_clust.size <- vector()

tips_names <- list()

for (i in 1:length(d)) {
  
  transm.df.cl.dat <- NULL
  
  clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
  
  clust.size_a <- c(clust.size_a, nrow(clus.read))
  
  name_clust.size <- c(name_clust.size, d[i])
  
  tips_names[[i]] <- as.character(clus.read[,1])
  
}





# Analysis for Subtype C -------



# Read the clusters' files for subtype C


dd <- list.files(path = paste0(sub.dir.rename), 
                 pattern = paste0("EAC_HIV_1_C_EAC_HIV_1_C_clusterPicks_cluster"),
                 all.files = FALSE,
                 full.names = FALSE, recursive = FALSE)

# Transmission clusters.

d <- dd

clust.size_c <- vector() # size of each cluster # table.simpact.trans.net.adv

name_clust.size <- vector()

tips_names <- list()

for (i in 1:length(d)) {
  
  transm.df.cl.dat <- NULL
  
  clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
  
  clust.size_c <- c(clust.size_c, nrow(clus.read))
  
  name_clust.size <- c(name_clust.size, d[i])
  
  tips_names[[i]] <- as.character(clus.read[,1])
  
}






# Analysis for Subtype D -------



# Read the clusters' files for subtype C


dd <- list.files(path = paste0(sub.dir.rename), 
                 pattern = paste0("EAC_HIV_1_D_EAC_HIV_1_D_clusterPicks_cluster"),
                 all.files = FALSE,
                 full.names = FALSE, recursive = FALSE)

# Transmission clusters.

d <- dd

clust.size_d <- vector() # size of each cluster # table.simpact.trans.net.adv

name_clust.size <- vector()

tips_names <- list()

for (i in 1:length(d)) {
  
  transm.df.cl.dat <- NULL
  
  clus.read <- read.table(file = paste0(paste0(sub.dir.rename,"/"),d[i]), header = FALSE) # Ids of each cluster
  
  clust.size_d <- c(clust.size_d, nrow(clus.read))
  
  name_clust.size <- c(name_clust.size, d[i])
  
  tips_names[[i]] <- as.character(clus.read[,1])
  
}





# Boxplot: summary of transmission clusters sizes --------------

par(mfrow=c(3,1))
boxplot(clust.size_a, main = "Histogram of cluster zise for subtype A")
boxplot(clust.size_c, main = "Histogram of cluster zise for subtype C")
boxplot(clust.size_d, main = "Histogram of cluster zise for subtype D")





# Transmission network accross countries per subtype -------------



# Subtype size per country: available seqs

# subtype A

A_BI_nbr <- length(names(A_BI)) # 25
A_DRC_nbr <- length(names(A_DRC)) #  43
A_KE_nbr <- length(names(A_KE)) # 1295
A_UGA_nbr <- length(names(A_UGA)) # 1180
A_TZ_nbr <- length(names(A_TZ)) # 363
A_RW_nbr <- length(names(A_RW)) # 189



# subtype C

C_BI_nbr <- length(names(C_BI)) # 289
C_DRC_nbr <- length(names(C_DRC)) #  22
C_KE_nbr <- length(names(C_KE)) # 104
C_UGA_nbr <- length(names(C_UGA)) # 66
C_TZ_nbr <- length(names(C_TZ)) # 388
C_RW_nbr <- length(names(C_RW)) # 24


# subtype D

D_BI_nbr <- length(names(D_BI)) # 5
D_DRC_nbr <- length(names(D_DRC)) # 28
D_KE_nbr <- length(names(D_KE)) # 145
D_UGA_nbr <- length(names(D_UGA)) # 1099
D_TZ_nbr <- length(names(D_TZ)) # 163
D_RW_nbr <- length(names(D_RW)) # 5



# Connectedness: number of transmission clusters wherein we could find seqs from two-by-two countries 

seq_A <- matrix(c(0,0,0,0,0,4, # column 
                  0,0,0,0,0,0, 
                  0,0,0,0,0,0, 
                  0,0,0,0,0,0, 
                  0,0,0,0,0,0,
                  0,0,0,0,0,0), nrow = 6)

colnames(seq_A) = rownames(seq_A) = c("BI", "KE", "TZ", "RW", "UG", "CD")



# Connectedness: number of transmission clusters wherein we could find seqs from two-by-two countries 

seq_C <- matrix(c(0,0,0,0,0,0, # column
                  0,0,0,0,0,0, 
                  0,0,0,5,8,3, 
                  0,0,0,0,1,1, 
                  0,0,0,0,0,0,
                  0,0,0,0,0,0), nrow = 6)

colnames(seq_C) = rownames(seq_C) = c("BI", "KE", "TZ", "RW", "UG", "CD")


# Connectedness: number of transmission clusters wherein we could find seqs from two-by-two countries 

seq_D <- matrix(c(0,0,1,0,1,0, # column
                  0,0,0,0,0,0, 
                  0,0,0,0,19,0, 
                  0,0,0,0,3,0, 
                  0,0,0,0,0,0,
                  0,0,0,0,0,0), nrow = 6)

colnames(seq_D) = rownames(seq_D) = c("BI", "KE", "TZ", "RW", "UG", "CD")


network_A <- graph_from_adjacency_matrix(seq_A, mode = "undirected", weighted = TRUE, diag = FALSE)
network_C <- graph_from_adjacency_matrix(seq_C, mode = "undirected", weighted = TRUE, diag = FALSE)
network_D <- graph_from_adjacency_matrix(seq_D, mode = "undirected", weighted = TRUE, diag = FALSE)




# plot for subtype A

V(network_A)
E(network_A)
vertex_attr(network_A)
edge_attr(network_A)

colorsnet <-  c("#FF0000", "#00FF00", "#FFFF00", "#00FFFF", "#FF00FF", "#0000FF")

vertex_attr(network_A)$color <- colorsnet
# vertex_attr(network_A)$size <- c(A_BI_nbr, A_KE_nbr, A_TZ_nbr, A_RW_nbr, A_UGA_nbr, A_DRC_nbr)


# plot for subtype C

V(network_C)
E(network_C)
vertex_attr(network_C)
edge_attr(network_C)


vertex_attr(network_C)$color <- colorsnet
# vertex_attr(network_C)$size <- c(C_BI_nbr, C_KE_nbr, C_TZ_nbr, C_RW_nbr, C_UGA_nbr, C_DRC_nbr)


# plot for subtype D

V(network_D)
E(network_D)
vertex_attr(network_D)
edge_attr(network_D)


vertex_attr(network_D)$color <- colorsnet
# vertex_attr(network_D)$size <- c(D_BI_nbr, D_KE_nbr, D_TZ_nbr, D_RW_nbr, D_UGA_nbr, D_DRC_nbr)



# Plot network

par(mfrow=c(1,3)) # c(3,1)
plot(network_A, 
     vertex.color=vertex_attr(network_A)$color,
     vertex.label=NA,
     # vertex.size=0.1*vertex_attr(network_A)$size,
     edge.width=edge_attr(network_A)$weight,
     main = "Subtype A"
)
plot(network_C, 
     vertex.color=vertex_attr(network_C)$color,
     vertex.label=NA,
     # vertex.size=0.1*vertex_attr(network_A)$size,
     edge.width=edge_attr(network_C)$weight,
     main = "Subtype C"
)
legend(x=-1.5, y=-1.1, c("Burundi","Kenya", "Tanzania", "Rwanda", "Uganda", "Dem Rep Congo"), pch=21,
       pt.bg=vertex_attr(network_A)$color, pt.cex=2, cex=.8, bty="n", ncol=3)
plot(network_D, 
     vertex.color=vertex_attr(network_D)$color,
     vertex.label=NA,
     # vertex.size=0.1*vertex_attr(network_A)$size,
     edge.width=edge_attr(network_D)$weight,
     main = "Subtype D"
)

