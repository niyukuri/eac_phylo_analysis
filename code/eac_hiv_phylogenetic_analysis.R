

# Full analysis of viral diversity and transmission clusters in EAC

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


# Count sequences per subtype per country -----


# A_BI_seqs <- names(A_BI) # 25
# C_BI_seqs <- names(C_BI) # 289
# D_BI_seqs <- names(D_BI) # 5
# 
# A_DRC_seqs <- names(A_DRC) #  43
# C_DRC_seqs <- names(C_DRC) #  22
# D_DRC_seqs <- names(D_DRC) # 28
# 
# A_KE_seqs <- names(A_KE) # 295
# C_KE_seqs <- names(C_KE) # 104
# D_KE_seqs <- names(D_KE) # 145
# 
# A_UGA_seqs <- names(A_UGA) # 1180
# C_UGA_seqs <- names(C_UGA) # 66
# D_UGA_seqs <- names(D_UGA) # 1099
# 
# A_TZ_seqs <- names(A_TZ) # 363
# C_TZ_seqs <- names(C_TZ) # 388
# D_TZ_seqs <- names(D_TZ) # 163
# 
# 
# A_RW_seqs <- names(A_RW) # 189
# C_RW_seqs <- names(C_RW) # 24
# D_RW_seqs <- names(D_RW) # 5

# Accession numbers ---------------

access_number_eac_seqs <- c(A_BI_seqs, C_BI_seqs, D_BI_seqs, 
                            A_DRC_seqs, C_DRC_seqs, D_DRC_seqs, 
                            A_KE_seqs, C_KE_seqs, D_KE_seqs, 
                            A_UGA_seqs, C_UGA_seqs, D_UGA_seqs, 
                            A_TZ_seqs, C_TZ_seqs, D_TZ_seqs, 
                            A_RW_seqs, C_RW_seqs, D_RW_seqs)

write(access_number_eac_seqs, file = "data/access_number_seqs_eac.txt")



# Table of sequences per subtype and country ----------------

A <- c(length(A_BI_seqs), length(A_DRC_seqs), length(A_KE_seqs), length(A_UGA_seqs), length(A_TZ_seqs), length(A_RW_seqs))
C <- c(length(C_BI_seqs), length(C_DRC_seqs), length(C_KE_seqs), length(C_UGA_seqs), length(C_TZ_seqs), length(C_RW_seqs))
D <- c(length(D_BI_seqs), length(D_DRC_seqs), length(D_KE_seqs), length(D_UGA_seqs), length(D_TZ_seqs), length(D_RW_seqs))

countries <- c("Burundi", "Congo", "Kenya", "Uganda", "Tanzania", "Rwanda")

count_sequences_ACD <- data.frame(countries, A, C, D)

names(count_sequences_ACD) <- c("Country", "Subtype A", "Subtype C", "Subtype D")


write.csv(count_sequences_ACD, file = "results/table_count_sequences_ACD.csv")



# Min and Max sampling dates per subtype per country --------------

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



# Compute age of internal nodes per subtype and country --------


# Substituion rates: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6007745/

# Landmark of HIV genome: https://www.hiv.lanl.gov/content/sequence/HIV/MAP/landmark.html 


# Subtype A for Burundi --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/A_Burundi_protease.Fasta > A_Burundi_protease.Fasta.nwk")


A_BI_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/A_Burundi_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

A_BI_tree.tips <- as.character(A_BI_tree.const$tip.label)

A_BI_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("A.BI.", "", A_BI_tree.tips), "\\d{4}"), format = "%Y"))

names(A_BI_tree.tips_dates) <- A_BI_tree.tips


A_BI_dater.tree <- treedater::dater(A_BI_tree.const, A_BI_tree.tips_dates, s=297, omega0 = 0.0013) 
# s is the length of sequence
# omega0 susbtitution rate: use POL value because for
# sequence sampled (1990–2011) +POL The genomic region encoding the viral enzymes protease, reverse transcriptase, and integrase)

class(A_BI_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
A_BI_dater.tree_sim.start.year <- round(A_BI_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
A_BI_dater.tree_first.transmission <- min(A_BI_dater.tree$sts)
A_BI_mrsd <- max(A_BI_dater.tree$sts)

A_BI_dates <- format(date_decimal(c(A_BI_mrsd, A_BI_dater.tree_sim.start.year)), "%Y-%m-%d")
A_BI_dater.tree$root.edge <-  - A_BI_dater.tree_sim.start.year


phylotree.plot_A_BI <- ggtree(A_BI_dater.tree,
                              mrsd = A_BI_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(A_BI_dater.tree_sim.start.year-1, A_BI_mrsd+1),
                     breaks = seq(from = A_BI_dater.tree_sim.start.year-1,
                                  to = A_BI_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_A_BI) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_A_BI.pdf",
       phylotree.plot_A_BI,
       width = 25,
       height = 15,
       units = "cm")


write.tree(A_BI_dater.tree, file = "phylotree.plot_A_BI.tree")

saveRDS(A_BI_dater.tree, file = "A_BI_dater.tree.RDS")

# A_BI_dater.tree <- readRDS("A_BI_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(A_BI_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


A_BI_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)

A_BI_nodes.long.df_raw <- A_BI_nodes.long.df
saveRDS(A_BI_nodes.long.df_raw, file = "A_BI_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

A_BI_nodes.long.df$intern.nodes <- A_BI_nodes.long.df$intern.nodes/sum(A_BI_nodes.long.df$intern.nodes)


saveRDS(A_BI_nodes.long.df, file = "A_BI_nodes.long.df.RDS")


# A_BI_nodes.long.df <- readRDS("A_BI_nodes.long.df.RDS")


A_BI_nodes.plot <- ggplot(data = A_BI_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(), 
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1932, 2000),
                     breaks = seq(from = 1934,
                                  to = 2000,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(A_BI_nodes.plot)

ggsave(filename = "A_BI_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype C for Burundi --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/C_Burundi_protease.Fasta > C_Burundi_protease.Fasta.nwk")


C_BI_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/C_Burundi_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

C_BI_tree.tips <- as.character(C_BI_tree.const$tip.label)

C_BI_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("C.BI.", "", C_BI_tree.tips), "\\d{4}"), format = "%Y"))

names(C_BI_tree.tips_dates) <- C_BI_tree.tips


C_BI_dater.tree <- treedater::dater(C_BI_tree.const, C_BI_tree.tips_dates, s=309, omega0 = 0.0013) # s is the length of sequence


class(C_BI_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
C_BI_dater.tree_sim.start.year <- round(C_BI_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
C_BI_dater.tree_first.transmission <- min(C_BI_dater.tree$sts)
C_BI_mrsd <- max(C_BI_dater.tree$sts)

C_BI_dates <- format(date_decimal(c(C_BI_mrsd, C_BI_dater.tree_sim.start.year)), "%Y-%m-%d")
C_BI_dater.tree$root.edge <-  - C_BI_dater.tree_sim.start.year


phylotree.plot_C_BI <- ggtree(C_BI_dater.tree,
                              mrsd = C_BI_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(C_BI_dater.tree_sim.start.year-1, C_BI_mrsd+1),
                     breaks = seq(from = C_BI_dater.tree_sim.start.year-1,
                                  to = C_BI_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_C_BI) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_C_BI.pdf",
       phylotree.plot_C_BI,
       width = 25,
       height = 15,
       units = "cm")


write.tree(C_BI_dater.tree, file = "phylotree.plot_C_BI.tree")

saveRDS(C_BI_dater.tree, file = "C_BI_dater.tree.RDS")

# C_BI_dater.tree <- readRDS("C_BI_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(C_BI_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


C_BI_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)

C_BI_nodes.long.df_raw <- C_BI_nodes.long.df
saveRDS(C_BI_nodes.long.df_raw, file = "C_BI_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)

# normalize internal nodes 

C_BI_nodes.long.df$intern.nodes <- C_BI_nodes.long.df$intern.nodes/sum(C_BI_nodes.long.df$intern.nodes)


saveRDS(C_BI_nodes.long.df, file = "C_BI_nodes.long.df.RDS")


# C_BI_nodes.long.df <- readRDS("C_BI_nodes.long.df.RDS")


C_BI_nodes.plot <- ggplot(data = C_BI_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1980, 2008),
                     breaks = seq(from = 1981,
                                  to = 2008,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(C_BI_nodes.plot)

ggsave(filename = "C_BI_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype D for Burundi --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/D_Burundi_protease.Fasta > D_Burundi_protease.Fasta.nwk")


D_BI_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/D_Burundi_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

D_BI_tree.tips <- as.character(D_BI_tree.const$tip.label)

D_BI_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("D.BI.", "", D_BI_tree.tips), "\\d{4}"), format = "%Y"))

names(D_BI_tree.tips_dates) <- D_BI_tree.tips


D_BI_dater.tree <- treedater::dater(D_BI_tree.const, D_BI_tree.tips_dates, s=309, omega0 = 0.0013) # s is the length of sequence


class(D_BI_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
D_BI_dater.tree_sim.start.year <- round(D_BI_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
D_BI_dater.tree_first.transmission <- min(D_BI_dater.tree$sts)
D_BI_mrsd <- max(D_BI_dater.tree$sts)

D_BI_dates <- format(date_decimal(c(D_BI_mrsd, D_BI_dater.tree_sim.start.year)), "%Y-%m-%d")
D_BI_dater.tree$root.edge <-  - D_BI_dater.tree_sim.start.year


phylotree.plot_D_BI <- ggtree(D_BI_dater.tree,
                              mrsd = D_BI_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(D_BI_dater.tree_sim.start.year-1, D_BI_mrsd+1),
                     breaks = seq(from = D_BI_dater.tree_sim.start.year-1,
                                  to = D_BI_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_D_BI) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_D_BI.pdf",
       phylotree.plot_D_BI,
       width = 25,
       height = 15,
       units = "cm")


write.tree(D_BI_dater.tree, file = "phylotree.plot_D_BI.tree")

saveRDS(D_BI_dater.tree, file = "D_BI_dater.tree.RDS")

# D_BI_dater.tree <- readRDS("D_BI_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(D_BI_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


D_BI_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)

D_BI_nodes.long.df_raw <- D_BI_nodes.long.df
saveRDS(D_BI_nodes.long.df_raw, file = "D_BI_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)

# normalize internal nodes 

D_BI_nodes.long.df$intern.nodes <- D_BI_nodes.long.df$intern.nodes/sum(D_BI_nodes.long.df$intern.nodes)


saveRDS(D_BI_nodes.long.df, file = "D_BI_nodes.long.df.RDS")


# D_BI_nodes.long.df <- readRDS("D_BI_nodes.long.df.RDS")


D_BI_nodes.plot <- ggplot(data = D_BI_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1813, 1993),
                     breaks = seq(from = 1814,
                                  to = 1993,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(D_BI_nodes.plot)

ggsave(filename = "D_BI_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype A for Congo --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/A_DRC_protease.Fasta > A_DRC_protease.Fasta.nwk")


A_DRC_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/A_DRC_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

A_DRC_tree.tips <- as.character(A_DRC_tree.const$tip.label)

A_DRC_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("A.DRC.", "", A_DRC_tree.tips), "\\d{4}"), format = "%Y"))

names(A_DRC_tree.tips_dates) <- A_DRC_tree.tips


A_DRC_dater.tree <- treedater::dater(A_DRC_tree.const, A_DRC_tree.tips_dates, s=297, omega0 = 0.0013) 
# s is the length of sequence
# omega0 susbtitution rate: use POL value because for
# sequence sampled (1990–2011) +POL The genomic region encoding the viral enzymes protease, reverse transcriptase, and integrase)

class(A_DRC_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
A_DRC_dater.tree_sim.start.year <- round(A_DRC_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
A_DRC_dater.tree_first.transmission <- min(A_DRC_dater.tree$sts)
A_DRC_mrsd <- max(A_DRC_dater.tree$sts)

A_DRC_dates <- format(date_decimal(c(A_DRC_mrsd, A_DRC_dater.tree_sim.start.year)), "%Y-%m-%d")
A_DRC_dater.tree$root.edge <-  - A_DRC_dater.tree_sim.start.year


phylotree.plot_A_DRC <- ggtree(A_DRC_dater.tree,
                               mrsd = A_DRC_dates[1],
                               size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(A_DRC_dater.tree_sim.start.year-1, A_DRC_mrsd+1),
                     breaks = seq(from = A_DRC_dater.tree_sim.start.year-1,
                                  to = A_DRC_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_A_DRC) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_A_DRC.pdf",
       phylotree.plot_A_DRC,
       width = 25,
       height = 15,
       units = "cm")


write.tree(A_DRC_dater.tree, file = "phylotree.plot_A_DRC.tree")

saveRDS(A_DRC_dater.tree, file = "A_DRC_dater.tree.RDS")

# A_DRC_dater.tree <- readRDS("A_DRC_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(A_DRC_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


A_DRC_nodes.long.df <- data.frame(calendaryear = i.vec,
                                  intern.nodes = int.node.vec)


A_DRC_nodes.long.df_raw <- A_DRC_nodes.long.df
saveRDS(A_DRC_nodes.long.df_raw, file = "A_DRC_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

A_DRC_nodes.long.df$intern.nodes <- A_DRC_nodes.long.df$intern.nodes/sum(A_DRC_nodes.long.df$intern.nodes)


saveRDS(A_DRC_nodes.long.df, file = "A_DRC_nodes.long.df.RDS")


# A_DRC_nodes.long.df <- readRDS("A_DRC_nodes.long.df.RDS")


A_DRC_nodes.plot <- ggplot(data = A_DRC_nodes.long.df,
                           aes(x = calendaryear,
                               y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(), 
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1957, 2000),
                     breaks = seq(from = 1958,
                                  to = 2000,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(A_DRC_nodes.plot)

ggsave(filename = "A_DRC_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype C for Congo --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/C_DRC_protease.Fasta > C_DRC_protease.Fasta.nwk")


C_DRC_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/C_DRC_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

C_DRC_tree.tips <- as.character(C_DRC_tree.const$tip.label)

C_DRC_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("C.DRC.", "", C_DRC_tree.tips), "\\d{4}"), format = "%Y"))

names(C_DRC_tree.tips_dates) <- C_DRC_tree.tips


C_DRC_dater.tree <- treedater::dater(C_DRC_tree.const, C_DRC_tree.tips_dates, s=297, omega0 = 0.0013) # s is the length of sequence


class(C_DRC_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
C_DRC_dater.tree_sim.start.year <- round(C_DRC_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
C_DRC_dater.tree_first.transmission <- min(C_DRC_dater.tree$sts)
C_DRC_mrsd <- max(C_DRC_dater.tree$sts)

C_DRC_dates <- format(date_decimal(c(C_DRC_mrsd, C_DRC_dater.tree_sim.start.year)), "%Y-%m-%d")
C_DRC_dater.tree$root.edge <-  - C_DRC_dater.tree_sim.start.year


phylotree.plot_C_DRC <- ggtree(C_DRC_dater.tree,
                               mrsd = C_DRC_dates[1],
                               size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(C_DRC_dater.tree_sim.start.year-1, C_DRC_mrsd+1),
                     breaks = seq(from = C_DRC_dater.tree_sim.start.year-1,
                                  to = C_DRC_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_C_DRC) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_C_DRC.pdf",
       phylotree.plot_C_DRC,
       width = 25,
       height = 15,
       units = "cm")


write.tree(C_DRC_dater.tree, file = "phylotree.plot_C_DRC.tree")

saveRDS(C_DRC_dater.tree, file = "C_DRC_dater.tree.RDS")

# C_DRC_dater.tree <- readRDS("C_DRC_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(C_DRC_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


C_DRC_nodes.long.df <- data.frame(calendaryear = i.vec,
                                  intern.nodes = int.node.vec)


C_DRC_nodes.long.df_raw <- C_DRC_nodes.long.df
saveRDS(C_DRC_nodes.long.df_raw, file = "C_DRC_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

C_DRC_nodes.long.df$intern.nodes <- C_DRC_nodes.long.df$intern.nodes/sum(C_DRC_nodes.long.df$intern.nodes)


saveRDS(C_DRC_nodes.long.df, file = "C_DRC_nodes.long.df.RDS")


# C_DRC_nodes.long.df <- readRDS("C_DRC_nodes.long.df.RDS")


C_DRC_nodes.plot <- ggplot(data = C_DRC_nodes.long.df,
                           aes(x = calendaryear,
                               y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1956, 1992),
                     breaks = seq(from = 1956,
                                  to = 1992,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(C_DRC_nodes.plot)

ggsave(filename = "C_DRC_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype D for Congo --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/D_DRC_protease.Fasta > D_DRC_protease.Fasta.nwk")


D_DRC_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/D_DRC_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

D_DRC_tree.tips <- as.character(D_DRC_tree.const$tip.label)

D_DRC_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("D.DRC.", "", D_DRC_tree.tips), "\\d{4}"), format = "%Y"))

names(D_DRC_tree.tips_dates) <- D_DRC_tree.tips


D_DRC_dater.tree <- treedater::dater(D_DRC_tree.const, D_DRC_tree.tips_dates, s=309, omega0 = 0.0013) # s is the length of sequence


class(D_DRC_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
D_DRC_dater.tree_sim.start.year <- round(D_DRC_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
D_DRC_dater.tree_first.transmission <- min(D_DRC_dater.tree$sts)
D_DRC_mrsd <- max(D_DRC_dater.tree$sts)

D_DRC_dates <- format(date_decimal(c(D_DRC_mrsd, D_DRC_dater.tree_sim.start.year)), "%Y-%m-%d")
D_DRC_dater.tree$root.edge <-  - D_DRC_dater.tree_sim.start.year


phylotree.plot_D_DRC <- ggtree(D_DRC_dater.tree,
                               mrsd = D_DRC_dates[1],
                               size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(D_DRC_dater.tree_sim.start.year-1, D_DRC_mrsd+1),
                     breaks = seq(from = D_DRC_dater.tree_sim.start.year-1,
                                  to = D_DRC_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_D_DRC) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_D_DRC.pdf",
       phylotree.plot_D_DRC,
       width = 25,
       height = 15,
       units = "cm")


write.tree(D_DRC_dater.tree, file = "phylotree.plot_D_DRC.tree")

saveRDS(D_DRC_dater.tree, file = "D_DRC_dater.tree.RDS")

# D_DRC_dater.tree <- readRDS("D_DRC_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(D_DRC_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


D_DRC_nodes.long.df <- data.frame(calendaryear = i.vec,
                                  intern.nodes = int.node.vec)




D_DRC_nodes.long.df_raw <- D_DRC_nodes.long.df
saveRDS(D_DRC_nodes.long.df_raw, file = "D_DRC_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

D_DRC_nodes.long.df$intern.nodes <- D_DRC_nodes.long.df$intern.nodes/sum(D_DRC_nodes.long.df$intern.nodes)


saveRDS(D_DRC_nodes.long.df, file = "D_DRC_nodes.long.df.RDS")


# D_DRC_nodes.long.df <- readRDS("D_DRC_nodes.long.df.RDS")


D_DRC_nodes.plot <- ggplot(data = D_DRC_nodes.long.df,
                           aes(x = calendaryear,
                               y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1960, 1990),
                     breaks = seq(from = 1960,
                                  to = 1990,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(D_DRC_nodes.plot)

ggsave(filename = "D_DRC_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype A for Kenya --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/A_Kenya_protease.Fasta > A_Kenya_protease.Fasta.nwk")


A_KE_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/A_Kenya_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

A_KE_tree.tips <- as.character(A_KE_tree.const$tip.label)

A_KE_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("A.KE.", "", A_KE_tree.tips), "\\d{4}"), format = "%Y"))

names(A_KE_tree.tips_dates) <- A_KE_tree.tips


A_KE_dater.tree <- treedater::dater(A_KE_tree.const, A_KE_tree.tips_dates, s=315, omega0 = 0.0013) 
# s is the length of sequence
# omega0 susbtitution rate: use POL value because for
# sequence sampled (1990–2011) +POL The genomic region encoding the viral enzymes protease, reverse transcriptase, and integrase)

class(A_KE_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
A_KE_dater.tree_sim.start.year <- round(A_KE_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
A_KE_dater.tree_first.transmission <- min(A_KE_dater.tree$sts)
A_KE_mrsd <- max(A_KE_dater.tree$sts)

A_KE_dates <- format(date_decimal(c(A_KE_mrsd, A_KE_dater.tree_sim.start.year)), "%Y-%m-%d")
A_KE_dater.tree$root.edge <-  - A_KE_dater.tree_sim.start.year


phylotree.plot_A_KE <- ggtree(A_KE_dater.tree,
                              mrsd = A_KE_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(A_KE_dater.tree_sim.start.year-1, A_KE_mrsd+1),
                     breaks = seq(from = A_KE_dater.tree_sim.start.year-1,
                                  to = A_KE_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_A_KE) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_A_KE.pdf",
       phylotree.plot_A_KE,
       width = 25,
       height = 15,
       units = "cm")


write.tree(A_KE_dater.tree, file = "phylotree.plot_A_KE.tree")

saveRDS(A_KE_dater.tree, file = "A_KE_dater.tree.RDS")

# A_KE_dater.tree <- readRDS("A_KE_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(A_KE_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


A_KE_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)


A_KE_nodes.long.df_raw <- A_KE_nodes.long.df
saveRDS(A_KE_nodes.long.df_raw, file = "A_KE_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

A_KE_nodes.long.df$intern.nodes <- A_KE_nodes.long.df$intern.nodes/sum(A_KE_nodes.long.df$intern.nodes)


saveRDS(A_KE_nodes.long.df, file = "A_KE_nodes.long.df.RDS")


# A_KE_nodes.long.df <- readRDS("A_KE_nodes.long.df.RDS")


A_KE_nodes.plot <- ggplot(data = A_KE_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(), 
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1973, 2015),
                     breaks = seq(from = 1973,
                                  to = 2015,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(A_KE_nodes.plot)

ggsave(filename = "A_KE_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype C for Kenya --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/C_Kenya_protease.Fasta > C_Kenya_protease.Fasta.nwk")


C_KE_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/C_Kenya_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

C_KE_tree.tips <- as.character(C_KE_tree.const$tip.label)

C_KE_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("C.KE.", "", C_KE_tree.tips), "\\d{4}"), format = "%Y"))

names(C_KE_tree.tips_dates) <- C_KE_tree.tips


C_KE_dater.tree <- treedater::dater(C_KE_tree.const, C_KE_tree.tips_dates, s=300, omega0 = 0.0013) # s is the length of sequence


class(C_KE_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
C_KE_dater.tree_sim.start.year <- round(C_KE_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
C_KE_dater.tree_first.transmission <- min(C_KE_dater.tree$sts)
C_KE_mrsd <- max(C_KE_dater.tree$sts)

C_KE_dates <- format(date_decimal(c(C_KE_mrsd, C_KE_dater.tree_sim.start.year)), "%Y-%m-%d")
C_KE_dater.tree$root.edge <-  - C_KE_dater.tree_sim.start.year


phylotree.plot_C_KE <- ggtree(C_KE_dater.tree,
                              mrsd = C_KE_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(C_KE_dater.tree_sim.start.year-1, C_KE_mrsd+1),
                     breaks = seq(from = C_KE_dater.tree_sim.start.year-1,
                                  to = C_KE_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_C_KE) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_C_KE.pdf",
       phylotree.plot_C_KE,
       width = 25,
       height = 15,
       units = "cm")


write.tree(C_KE_dater.tree, file = "phylotree.plot_C_KE.tree")

saveRDS(C_KE_dater.tree, file = "C_KE_dater.tree.RDS")

# C_KE_dater.tree <- readRDS("C_KE_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(C_KE_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


C_KE_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)



C_KE_nodes.long.df_raw <- C_KE_nodes.long.df
saveRDS(C_KE_nodes.long.df_raw, file = "C_KE_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

C_KE_nodes.long.df$intern.nodes <- C_KE_nodes.long.df$intern.nodes/sum(C_KE_nodes.long.df$intern.nodes)


saveRDS(C_KE_nodes.long.df, file = "C_KE_nodes.long.df.RDS")


# C_KE_nodes.long.df <- readRDS("C_KE_nodes.long.df.RDS")


C_KE_nodes.plot <- ggplot(data = C_KE_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1971, 2015),
                     breaks = seq(from = 1971,
                                  to = 2015,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(C_KE_nodes.plot)

ggsave(filename = "C_KE_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype D for Kenya --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/D_Kenya_protease.Fasta > D_Kenya_protease.Fasta.nwk")


D_KE_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/D_Kenya_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

D_KE_tree.tips <- as.character(D_KE_tree.const$tip.label)

D_KE_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("D.KE.", "", D_KE_tree.tips), "\\d{4}"), format = "%Y"))

names(D_KE_tree.tips_dates) <- D_KE_tree.tips


D_KE_dater.tree <- treedater::dater(D_KE_tree.const, D_KE_tree.tips_dates, s=312, omega0 = 0.0013) # s is the length of sequence


class(D_KE_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
D_KE_dater.tree_sim.start.year <- round(D_KE_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
D_KE_dater.tree_first.transmission <- min(D_KE_dater.tree$sts)
D_KE_mrsd <- max(D_KE_dater.tree$sts)

D_KE_dates <- format(date_decimal(c(D_KE_mrsd, D_KE_dater.tree_sim.start.year)), "%Y-%m-%d")
D_KE_dater.tree$root.edge <-  - D_KE_dater.tree_sim.start.year


phylotree.plot_D_KE <- ggtree(D_KE_dater.tree,
                              mrsd = D_KE_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(D_KE_dater.tree_sim.start.year-1, D_KE_mrsd+1),
                     breaks = seq(from = D_KE_dater.tree_sim.start.year-1,
                                  to = D_KE_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_D_KE) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_D_KE.pdf",
       phylotree.plot_D_KE,
       width = 25,
       height = 15,
       units = "cm")


write.tree(D_KE_dater.tree, file = "phylotree.plot_D_KE.tree")

saveRDS(D_KE_dater.tree, file = "D_KE_dater.tree.RDS")

# D_KE_dater.tree <- readRDS("D_KE_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(D_KE_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


D_KE_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)


D_KE_nodes.long.df_raw <- D_KE_nodes.long.df
saveRDS(D_KE_nodes.long.df_raw, file = "D_KE_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

D_KE_nodes.long.df$intern.nodes <- D_KE_nodes.long.df$intern.nodes/sum(D_KE_nodes.long.df$intern.nodes)


saveRDS(D_KE_nodes.long.df, file = "D_KE_nodes.long.df.RDS")


# D_KE_nodes.long.df <- readRDS("D_KE_nodes.long.df.RDS")


D_KE_nodes.plot <- ggplot(data = D_KE_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1970, 2012),
                     breaks = seq(from = 1970,
                                  to = 2012,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(D_KE_nodes.plot)

ggsave(filename = "D_KE_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype A for Uganda --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/A_Uganda_protease.Fasta > A_Uganda_protease.Fasta.nwk")


A_UGA_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/A_Uganda_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

A_UGA_tree.tips <- as.character(A_UGA_tree.const$tip.label)

A_UGA_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("A.UGA.", "", A_UGA_tree.tips), "\\d{4}"), format = "%Y"))

names(A_UGA_tree.tips_dates) <- A_UGA_tree.tips


A_UGA_dater.tree <- treedater::dater(A_UGA_tree.const, A_UGA_tree.tips_dates, s=321, omega0 = 0.0015, minblen = 0.005) 
# s is the length of sequence
# omega0 susbtitution rate: use POL value because for
# sequence sampled (1990–2011) +POL The genomic region encoding the viral enzymes protease, reverse transcriptase, and integrase)

class(A_UGA_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
A_UGA_dater.tree_sim.start.year <- round(A_UGA_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
A_UGA_dater.tree_first.transmission <- min(A_UGA_dater.tree$sts)
A_UGA_mrsd <- max(A_UGA_dater.tree$sts)

A_UGA_dates <- format(date_decimal(c(A_UGA_mrsd, A_UGA_dater.tree_sim.start.year)), "%Y-%m-%d")
A_UGA_dater.tree$root.edge <-  - A_UGA_dater.tree_sim.start.year


phylotree.plot_A_UGA <- ggtree(A_UGA_dater.tree,
                               mrsd = A_UGA_dates[1],
                               size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(A_UGA_dater.tree_sim.start.year-1, A_UGA_mrsd+1),
                     breaks = seq(from = A_UGA_dater.tree_sim.start.year-1,
                                  to = A_UGA_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_A_UGA) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_A_UGA.pdf",
       phylotree.plot_A_UGA,
       width = 25,
       height = 15,
       units = "cm")


write.tree(A_UGA_dater.tree, file = "phylotree.plot_A_UGA.tree")

saveRDS(A_UGA_dater.tree, file = "A_UGA_dater.tree.RDS")

# A_UGA_dater.tree <- readRDS("A_UGA_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(A_UGA_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


A_UGA_nodes.long.df <- data.frame(calendaryear = i.vec,
                                  intern.nodes = int.node.vec)


A_UGA_nodes.long.df_raw <- A_UGA_nodes.long.df
saveRDS(A_UGA_nodes.long.df_raw, file = "A_UGA_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

A_UGA_nodes.long.df$intern.nodes <- A_UGA_nodes.long.df$intern.nodes/sum(A_UGA_nodes.long.df$intern.nodes)


saveRDS(A_UGA_nodes.long.df, file = "A_UGA_nodes.long.df.RDS")


# A_UGA_nodes.long.df <- readRDS("A_UGA_nodes.long.df.RDS")


A_UGA_nodes.plot <- ggplot(data = A_UGA_nodes.long.df,
                           aes(x = calendaryear,
                               y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(), 
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1958, 2012),
                     breaks = seq(from = 1958,
                                  to = 2012,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(A_UGA_nodes.plot)

ggsave(filename = "A_UGA_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype C for Uganda --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/C_Uganda_protease.Fasta > C_Uganda_protease.Fasta.nwk")


C_UGA_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/C_Uganda_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

C_UGA_tree.tips <- as.character(C_UGA_tree.const$tip.label)

C_UGA_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("C.UGA.", "", C_UGA_tree.tips), "\\d{4}"), format = "%Y"))

names(C_UGA_tree.tips_dates) <- C_UGA_tree.tips


C_UGA_dater.tree <- treedater::dater(C_UGA_tree.const, C_UGA_tree.tips_dates, s=297, omega0 = 0.0013) # s is the length of sequence

# s=321, omega0 = 0.0015, minblen = 0.005

class(C_UGA_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
C_UGA_dater.tree_sim.start.year <- round(C_UGA_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
C_UGA_dater.tree_first.transmission <- min(C_UGA_dater.tree$sts)
C_UGA_mrsd <- max(C_UGA_dater.tree$sts)

C_UGA_dates <- format(date_decimal(c(C_UGA_mrsd, C_UGA_dater.tree_sim.start.year)), "%Y-%m-%d")
C_UGA_dater.tree$root.edge <-  - C_UGA_dater.tree_sim.start.year


phylotree.plot_C_UGA <- ggtree(C_UGA_dater.tree,
                               mrsd = C_UGA_dates[1],
                               size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(C_UGA_dater.tree_sim.start.year-1, C_UGA_mrsd+1),
                     breaks = seq(from = C_UGA_dater.tree_sim.start.year-1,
                                  to = C_UGA_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_C_UGA) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_C_UGA.pdf",
       phylotree.plot_C_UGA,
       width = 25,
       height = 15,
       units = "cm")


write.tree(C_UGA_dater.tree, file = "phylotree.plot_C_UGA.tree")

saveRDS(C_UGA_dater.tree, file = "C_UGA_dater.tree.RDS")

# C_UGA_dater.tree <- readRDS("C_UGA_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(C_UGA_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


C_UGA_nodes.long.df <- data.frame(calendaryear = i.vec,
                                  intern.nodes = int.node.vec)


C_UGA_nodes.long.df_raw <- C_UGA_nodes.long.df
saveRDS(C_UGA_nodes.long.df_raw, file = "C_UGA_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

C_UGA_nodes.long.df$intern.nodes <- C_UGA_nodes.long.df$intern.nodes/sum(C_UGA_nodes.long.df$intern.nodes)


saveRDS(C_UGA_nodes.long.df, file = "C_UGA_nodes.long.df.RDS")


# C_UGA_nodes.long.df <- readRDS("C_UGA_nodes.long.df.RDS")


C_UGA_nodes.plot <- ggplot(data = C_UGA_nodes.long.df,
                           aes(x = calendaryear,
                               y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1936, 2012),
                     breaks = seq(from = 1936,
                                  to = 2012,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(C_UGA_nodes.plot)

ggsave(filename = "C_UGA_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype D for Uganda --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/D_Uganda_protease.Fasta > D_Uganda_protease.Fasta.nwk")


D_UGA_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/D_Uganda_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

D_UGA_tree.tips <- as.character(D_UGA_tree.const$tip.label)

D_UGA_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("D.UGA.", "", D_UGA_tree.tips), "\\d{4}"), format = "%Y"))

names(D_UGA_tree.tips_dates) <- D_UGA_tree.tips


D_UGA_dater.tree <- treedater::dater(D_UGA_tree.const, D_UGA_tree.tips_dates, s=300, omega0 = 0.0013) # s is the length of sequence

# s=321, omega0 = 0.0015, minblen = 0.005

class(D_UGA_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
D_UGA_dater.tree_sim.start.year <- round(D_UGA_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
D_UGA_dater.tree_first.transmission <- min(D_UGA_dater.tree$sts)
D_UGA_mrsd <- max(D_UGA_dater.tree$sts)

D_UGA_dates <- format(date_decimal(c(D_UGA_mrsd, D_UGA_dater.tree_sim.start.year)), "%Y-%m-%d")
D_UGA_dater.tree$root.edge <-  - D_UGA_dater.tree_sim.start.year


phylotree.plot_D_UGA <- ggtree(D_UGA_dater.tree,
                               mrsd = D_UGA_dates[1],
                               size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(D_UGA_dater.tree_sim.start.year-1, D_UGA_mrsd+1),
                     breaks = seq(from = D_UGA_dater.tree_sim.start.year-1,
                                  to = D_UGA_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_D_UGA) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_D_UGA.pdf",
       phylotree.plot_D_UGA,
       width = 25,
       height = 15,
       units = "cm")


write.tree(D_UGA_dater.tree, file = "phylotree.plot_D_UGA.tree")

saveRDS(D_UGA_dater.tree, file = "D_UGA_dater.tree.RDS")

# D_UGA_dater.tree <- readRDS("D_UGA_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(D_UGA_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


D_UGA_nodes.long.df <- data.frame(calendaryear = i.vec,
                                  intern.nodes = int.node.vec)



D_UGA_nodes.long.df_raw <- D_UGA_nodes.long.df
saveRDS(D_UGA_nodes.long.df_raw, file = "D_UGA_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

D_UGA_nodes.long.df$intern.nodes <- D_UGA_nodes.long.df$intern.nodes/sum(D_UGA_nodes.long.df$intern.nodes)


saveRDS(D_UGA_nodes.long.df, file = "D_UGA_nodes.long.df.RDS")


# D_UGA_nodes.long.df <- readRDS("D_UGA_nodes.long.df.RDS")


D_UGA_nodes.plot <- ggplot(data = D_UGA_nodes.long.df,
                           aes(x = calendaryear,
                               y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1948, 2012),
                     breaks = seq(from = 1948,
                                  to = 2012,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(D_UGA_nodes.plot)

ggsave(filename = "D_UGA_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")






# Subtype A for Tanzania --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/A_Tanzania_protease.Fasta > A_Tanzania_protease.Fasta.nwk")


A_TZ_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/A_Tanzania_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

A_TZ_tree.tips <- as.character(A_TZ_tree.const$tip.label)

A_TZ_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("A.TZ.", "", A_TZ_tree.tips), "\\d{4}"), format = "%Y"))

names(A_TZ_tree.tips_dates) <- A_TZ_tree.tips


A_TZ_dater.tree <- treedater::dater(A_TZ_tree.const, A_TZ_tree.tips_dates, s=306, omega0 = 0.0015, minblen = 0.005) 
# s is the length of sequence
# omega0 susbtitution rate: use POL value because for
# sequence sampled (1990–2011) +POL The genomic region encoding the viral enzymes protease, reverse transcriptase, and integrase)

class(A_TZ_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
A_TZ_dater.tree_sim.start.year <- round(A_TZ_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
A_TZ_dater.tree_first.transmission <- min(A_TZ_dater.tree$sts)
A_TZ_mrsd <- max(A_TZ_dater.tree$sts)

A_TZ_dates <- format(date_decimal(c(A_TZ_mrsd, A_TZ_dater.tree_sim.start.year)), "%Y-%m-%d")
A_TZ_dater.tree$root.edge <-  - A_TZ_dater.tree_sim.start.year


phylotree.plot_A_TZ <- ggtree(A_TZ_dater.tree,
                              mrsd = A_TZ_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(A_TZ_dater.tree_sim.start.year-1, A_TZ_mrsd+1),
                     breaks = seq(from = A_TZ_dater.tree_sim.start.year-1,
                                  to = A_TZ_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_A_TZ) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_A_TZ.pdf",
       phylotree.plot_A_TZ,
       width = 25,
       height = 15,
       units = "cm")


write.tree(A_TZ_dater.tree, file = "phylotree.plot_A_TZ.tree")

saveRDS(A_TZ_dater.tree, file = "A_TZ_dater.tree.RDS")

# A_TZ_dater.tree <- readRDS("A_TZ_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(A_TZ_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


A_TZ_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)


A_TZ_nodes.long.df_raw <- A_TZ_nodes.long.df
saveRDS(A_TZ_nodes.long.df_raw, file = "A_TZ_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)



# normalize internal nodes 

A_TZ_nodes.long.df$intern.nodes <- A_TZ_nodes.long.df$intern.nodes/sum(A_TZ_nodes.long.df$intern.nodes)


saveRDS(A_TZ_nodes.long.df, file = "A_TZ_nodes.long.df.RDS")


# A_TZ_nodes.long.df <- readRDS("A_TZ_nodes.long.df.RDS")


A_TZ_nodes.plot <- ggplot(data = A_TZ_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(), 
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1978, 2012),
                     breaks = seq(from = 1978,
                                  to = 2012,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(A_TZ_nodes.plot)

ggsave(filename = "A_TZ_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype C for Tanzania --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/C_Tanzania_protease.Fasta > C_Tanzania_protease.Fasta.nwk")


C_TZ_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/C_Tanzania_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

C_TZ_tree.tips <- as.character(C_TZ_tree.const$tip.label)

C_TZ_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("C.TZ.", "", C_TZ_tree.tips), "\\d{4}"), format = "%Y"))

names(C_TZ_tree.tips_dates) <- C_TZ_tree.tips


C_TZ_dater.tree <- treedater::dater(C_TZ_tree.const, C_TZ_tree.tips_dates, s=297, omega0 = 0.0013) # s is the length of sequence

# s=321, omega0 = 0.0015, minblen = 0.005

class(C_TZ_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
C_TZ_dater.tree_sim.start.year <- round(C_TZ_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
C_TZ_dater.tree_first.transmission <- min(C_TZ_dater.tree$sts)
C_TZ_mrsd <- max(C_TZ_dater.tree$sts)

C_TZ_dates <- format(date_decimal(c(C_TZ_mrsd, C_TZ_dater.tree_sim.start.year)), "%Y-%m-%d")
C_TZ_dater.tree$root.edge <-  - C_TZ_dater.tree_sim.start.year


phylotree.plot_C_TZ <- ggtree(C_TZ_dater.tree,
                              mrsd = C_TZ_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(C_TZ_dater.tree_sim.start.year-1, C_TZ_mrsd+1),
                     breaks = seq(from = C_TZ_dater.tree_sim.start.year-1,
                                  to = C_TZ_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_C_TZ) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_C_TZ.pdf",
       phylotree.plot_C_TZ,
       width = 25,
       height = 15,
       units = "cm")


write.tree(C_TZ_dater.tree, file = "phylotree.plot_C_TZ.tree")

saveRDS(C_TZ_dater.tree, file = "C_TZ_dater.tree.RDS")

# C_TZ_dater.tree <- readRDS("C_TZ_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(C_TZ_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


C_TZ_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)


C_TZ_nodes.long.df_raw <- C_TZ_nodes.long.df
saveRDS(C_TZ_nodes.long.df_raw, file = "C_TZ_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

C_TZ_nodes.long.df$intern.nodes <- C_TZ_nodes.long.df$intern.nodes/sum(C_TZ_nodes.long.df$intern.nodes)


saveRDS(C_TZ_nodes.long.df, file = "C_TZ_nodes.long.df.RDS")


# C_TZ_nodes.long.df <- readRDS("C_TZ_nodes.long.df.RDS")


C_TZ_nodes.plot <- ggplot(data = C_TZ_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1966, 2012),
                     breaks = seq(from = 1966,
                                  to = 2012,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(C_TZ_nodes.plot)

ggsave(filename = "C_TZ_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype D for Tanzania --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/D_Tanzania_protease.Fasta > D_Tanzania_protease.Fasta.nwk")


D_TZ_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/D_Tanzania_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

D_TZ_tree.tips <- as.character(D_TZ_tree.const$tip.label)

D_TZ_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("D.TZ.", "", D_TZ_tree.tips), "\\d{4}"), format = "%Y"))

names(D_TZ_tree.tips_dates) <- D_TZ_tree.tips


D_TZ_dater.tree <- treedater::dater(D_TZ_tree.const, D_TZ_tree.tips_dates, s=297, omega0 = 0.0013) # s is the length of sequence

# s=321, omega0 = 0.0015, minblen = 0.005

class(D_TZ_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
D_TZ_dater.tree_sim.start.year <- round(D_TZ_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
D_TZ_dater.tree_first.transmission <- min(D_TZ_dater.tree$sts)
D_TZ_mrsd <- max(D_TZ_dater.tree$sts)

D_TZ_dates <- format(date_decimal(c(D_TZ_mrsd, D_TZ_dater.tree_sim.start.year)), "%Y-%m-%d")
D_TZ_dater.tree$root.edge <-  - D_TZ_dater.tree_sim.start.year


phylotree.plot_D_TZ <- ggtree(D_TZ_dater.tree,
                              mrsd = D_TZ_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(D_TZ_dater.tree_sim.start.year-1, D_TZ_mrsd+1),
                     breaks = seq(from = D_TZ_dater.tree_sim.start.year-1,
                                  to = D_TZ_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_D_TZ) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_D_TZ.pdf",
       phylotree.plot_D_TZ,
       width = 25,
       height = 15,
       units = "cm")


write.tree(D_TZ_dater.tree, file = "phylotree.plot_D_TZ.tree")

saveRDS(D_TZ_dater.tree, file = "D_TZ_dater.tree.RDS")

# D_TZ_dater.tree <- readRDS("D_TZ_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(D_TZ_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


D_TZ_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)


D_TZ_nodes.long.df_raw <- D_TZ_nodes.long.df
saveRDS(D_TZ_nodes.long.df_raw, file = "D_TZ_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

D_TZ_nodes.long.df$intern.nodes <- D_TZ_nodes.long.df$intern.nodes/sum(D_TZ_nodes.long.df$intern.nodes)


saveRDS(D_TZ_nodes.long.df, file = "D_TZ_nodes.long.df.RDS")


# D_TZ_nodes.long.df <- readRDS("D_TZ_nodes.long.df.RDS")


D_TZ_nodes.plot <- ggplot(data = D_TZ_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1965, 2009),
                     breaks = seq(from = 1965,
                                  to = 2009,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(D_TZ_nodes.plot)

ggsave(filename = "D_TZ_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")




# Subtype A for Rwanda --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/A_Rwanda_protease.Fasta > A_Rwanda_protease.Fasta.nwk")


A_RW_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/A_Rwanda_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

A_RW_tree.tips <- as.character(A_RW_tree.const$tip.label)

A_RW_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("A.RW.", "", A_RW_tree.tips), "\\d{4}"), format = "%Y"))

names(A_RW_tree.tips_dates) <- A_RW_tree.tips


A_RW_dater.tree <- treedater::dater(A_RW_tree.const, A_RW_tree.tips_dates, s=297, omega0 = 0.0015, minblen = 0.02) 
# s is the length of sequence
# omega0 susbtitution rate: use POL value because for
# sequence sampled (1990–2011) +POL The genomic region encoding the viral enzymes protease, reverse transcriptase, and integrase)

class(A_RW_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
A_RW_dater.tree_sim.start.year <- round(A_RW_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
A_RW_dater.tree_first.transmission <- min(A_RW_dater.tree$sts)
A_RW_mrsd <- max(A_RW_dater.tree$sts)

A_RW_dates <- format(date_decimal(c(A_RW_mrsd, A_RW_dater.tree_sim.start.year)), "%Y-%m-%d")
A_RW_dater.tree$root.edge <-  - A_RW_dater.tree_sim.start.year


phylotree.plot_A_RW <- ggtree(A_RW_dater.tree,
                              mrsd = A_RW_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(A_RW_dater.tree_sim.start.year-1, A_RW_mrsd+1),
                     breaks = seq(from = A_RW_dater.tree_sim.start.year-1,
                                  to = A_RW_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_A_RW) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_A_RW.pdf",
       phylotree.plot_A_RW,
       width = 25,
       height = 15,
       units = "cm")


write.tree(A_RW_dater.tree, file = "phylotree.plot_A_RW.tree")

saveRDS(A_RW_dater.tree, file = "A_RW_dater.tree.RDS")

# A_RW_dater.tree <- readRDS("A_RW_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(A_RW_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


A_RW_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)


A_RW_nodes.long.df_raw <- A_RW_nodes.long.df
saveRDS(A_RW_nodes.long.df_raw, file = "A_RW_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

A_RW_nodes.long.df$intern.nodes <- A_RW_nodes.long.df$intern.nodes/sum(A_RW_nodes.long.df$intern.nodes)


saveRDS(A_RW_nodes.long.df, file = "A_RW_nodes.long.df.RDS")


# A_RW_nodes.long.df <- readRDS("A_RW_nodes.long.df.RDS")


A_RW_nodes.plot <- ggplot(data = A_RW_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x=element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(), 
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1920, 2010),
                     breaks = seq(from = 1920,
                                  to = 2010,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(A_RW_nodes.plot)

ggsave(filename = "A_RW_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype C for Rwanda --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/C_Rwanda_protease.Fasta > C_Rwanda_protease.Fasta.nwk")


C_RW_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/C_Rwanda_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

C_RW_tree.tips <- as.character(C_RW_tree.const$tip.label)

C_RW_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("C.RW.", "", C_RW_tree.tips), "\\d{4}"), format = "%Y"))

names(C_RW_tree.tips_dates) <- C_RW_tree.tips


C_RW_dater.tree <- treedater::dater(C_RW_tree.const, C_RW_tree.tips_dates, s=348, omega0 = 0.0013) # s is the length of sequence

# s=321, omega0 = 0.0015, minblen = 0.005

class(C_RW_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
C_RW_dater.tree_sim.start.year <- round(C_RW_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
C_RW_dater.tree_first.transmission <- min(C_RW_dater.tree$sts)
C_RW_mrsd <- max(C_RW_dater.tree$sts)

C_RW_dates <- format(date_decimal(c(C_RW_mrsd, C_RW_dater.tree_sim.start.year)), "%Y-%m-%d")
C_RW_dater.tree$root.edge <-  - C_RW_dater.tree_sim.start.year


phylotree.plot_C_RW <- ggtree(C_RW_dater.tree,
                              mrsd = C_RW_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(C_RW_dater.tree_sim.start.year-1, C_RW_mrsd+1),
                     breaks = seq(from = C_RW_dater.tree_sim.start.year-1,
                                  to = C_RW_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_C_RW) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_C_RW.pdf",
       phylotree.plot_C_RW,
       width = 25,
       height = 15,
       units = "cm")


write.tree(C_RW_dater.tree, file = "phylotree.plot_C_RW.tree")

saveRDS(C_RW_dater.tree, file = "C_RW_dater.tree.RDS")

# C_RW_dater.tree <- readRDS("C_RW_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(C_RW_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


C_RW_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)


C_RW_nodes.long.df_raw <- C_RW_nodes.long.df
saveRDS(C_RW_nodes.long.df_raw, file = "C_RW_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)




# normalize internal nodes 

C_RW_nodes.long.df$intern.nodes <- C_RW_nodes.long.df$intern.nodes/sum(C_RW_nodes.long.df$intern.nodes)


saveRDS(C_RW_nodes.long.df, file = "C_RW_nodes.long.df.RDS")


# C_RW_nodes.long.df <- readRDS("C_RW_nodes.long.df.RDS")


C_RW_nodes.plot <- ggplot(data = C_RW_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1968, 2012),
                     breaks = seq(from = 1968,
                                  to = 2012,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(C_RW_nodes.plot)

ggsave(filename = "C_RW_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")





# Subtype D for Rwanda --------


system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/D_Rwanda_protease.Fasta > D_Rwanda_protease.Fasta.nwk")


D_RW_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/D_Rwanda_protease.Fasta.nwk")


# Match dates and phylogenetic tree leaves

D_RW_tree.tips <- as.character(D_RW_tree.const$tip.label)

D_RW_tree.tips_dates <- year(as.Date(stri_extract_first_regex(gsub("D.RW.", "", D_RW_tree.tips), "\\d{4}"), format = "%Y"))

names(D_RW_tree.tips_dates) <- D_RW_tree.tips


D_RW_dater.tree <- treedater::dater(D_RW_tree.const, D_RW_tree.tips_dates, s=297, omega0 = 0.0013) # s is the length of sequence

# s=321, omega0 = 0.0015, minblen = 0.005

class(D_RW_dater.tree) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
D_RW_dater.tree_sim.start.year <- round(D_RW_dater.tree$timeOfMRCA) # akantu dater.tree$timeToMRCA
D_RW_dater.tree_first.transmission <- min(D_RW_dater.tree$sts)
D_RW_mrsd <- max(D_RW_dater.tree$sts)

D_RW_dates <- format(date_decimal(c(D_RW_mrsd, D_RW_dater.tree_sim.start.year)), "%Y-%m-%d")
D_RW_dater.tree$root.edge <-  - D_RW_dater.tree_sim.start.year


phylotree.plot_D_RW <- ggtree(D_RW_dater.tree,
                              mrsd = D_RW_dates[1],
                              size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(D_RW_dater.tree_sim.start.year-1, D_RW_mrsd+1),
                     breaks = seq(from = D_RW_dater.tree_sim.start.year-1,
                                  to = D_RW_mrsd+1,
                                  by = 2)) +
  xlab("Time") +
  ylab("Proportion of internal nodes")

print(phylotree.plot_D_RW) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("phylotree.plot_D_RW.pdf",
       phylotree.plot_D_RW,
       width = 25,
       height = 15,
       units = "cm")


write.tree(D_RW_dater.tree, file = "phylotree.plot_D_RW.tree")

saveRDS(D_RW_dater.tree, file = "D_RW_dater.tree.RDS")

# D_RW_dater.tree <- readRDS("D_RW_dater.tree.RDS")


# Internal nodels

# Node age with picante package

N <- picante::node.age(D_RW_dater.tree)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti

max.val <- round(max(int.node.age))+1
min.val <- round(min(int.node.age))-1
step.int <- 1

d <- (max.val-min.val)/step.int

dt.node.age.dt <- int.node.age


i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  
  inf <- min.val-1+i
  sup <- min.val+i
  
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
  
}


D_RW_nodes.long.df <- data.frame(calendaryear = i.vec,
                                 intern.nodes = int.node.vec)


D_RW_nodes.long.df_raw <- D_RW_nodes.long.df
saveRDS(D_RW_nodes.long.df_raw, file = "D_RW_nodes.long.df_raw.RDS") # raw number of internal nodes (un-scale)


# normalize internal nodes 

D_RW_nodes.long.df$intern.nodes <- D_RW_nodes.long.df$intern.nodes/sum(D_RW_nodes.long.df$intern.nodes)


saveRDS(D_RW_nodes.long.df, file = "D_RW_nodes.long.df.RDS")


# D_RW_nodes.long.df <- readRDS("D_RW_nodes.long.df.RDS")


D_RW_nodes.plot <- ggplot(data = D_RW_nodes.long.df,
                          aes(x = calendaryear,
                              y = intern.nodes)) + 
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1982, 2010),
                     breaks = seq(from = 1982,
                                  to = 2010,
                                  by = 2)) +
  xlab("Time")+
  ylab("Proportion of internal nodes")

print(D_RW_nodes.plot)

ggsave(filename = "D_RW_nodes.plot.pdf",
       plot = nodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis",
       width = 20, height = 10, units = "cm")



# Trend of internal nodes per subtype and country --------------


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


# saveRDS(EAC_df_raw, file="EAC_df_raw.RDS")


EAC_df_raw <- readRDS("EAC_df_raw.RDS")




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





# Transmission clusters per subtype in EAC --------


# Dealing with seqs lengths differences:

# Each subtype per country has its own length 
# Combine seqs per subtyes for all countries (done with MEGA)

# length(A_TZ$A1.TZ.1997.97TZ02.AF361872)
# length(A_BI$A1.BI.2002.02BU_U1808.AM260224)
# length(A_RW$A1.RW.1992.92RW008.AB253421)
# length(A_KE$A1.KE.1986.ML013_10.AY322184)
# length(A_TZ$A1.TZ.1997.97TZ02.AF361872)
# length(A_DRC$A1.CD.2002.02CD_LBS024.AM040996)
# 
# 
# length(C_TZ$C.TZ.1996.6567.AY036408)
# length(C_BI$C.BI.2002.02BU_B1027.AM260260)
# length(C_RW$`C.RW.2005.8222-PD.KC513726`)
# length(C_KE$C.KE.1991.KNH1268.AY945738)
# length(C_TZ$C.TZ.1996.6567.AY036408)
# length(C_DRC$C.CD.2002.02CD_KSS055.AM041050)
# 
# 
# length(D_TZ$D.TZ.1996.2975.AY036303)
# length(D_BI$D.BI.2002.02BU_U1916.AM260230)
# length(D_RW$`D.RW.2005.105-7101.JQ364514`)
# length(D_KE$D.KE.1996.96KE006.AY492783)
# length(D_TZ$D.TZ.1996.2975.AY036303)
# length(D_DRC$D.CD.1983.ELI.K03454)


# We have same length per subtype: EAC_HIV_1_A.fas, EAC_HIV_1_C.fas, and EAC_HIV_1_D.fas

A_EAC <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/EAC_HIV_1_A.fas") # n = 2095, len = 348
C_EAC <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/EAC_HIV_1_C.fas") # n = 893, len = 309
D_EAC <- read.fasta("/home/david/Dropbox/eac_phylo_analysis/data/EAC_HIV_1_D.fas") # n = 1445, len = 312

# length(A_EAC$A1.TZ.1997.97TZ02.AF361872)
# length(A_EAC$A1.BI.2002.02BU_U1808.AM260224)
# length(A_EAC$A1.RW.1992.92RW008.AB253421)
# length(A_EAC$A1.KE.1986.ML013_10.AY322184)
# length(A_EAC$A1.TZ.1997.97TZ02.AF361872)
# length(A_EAC$A1.CD.2002.02CD_LBS024.AM040996)
# 
# 
# length(C_EAC$C.TZ.1996.6567.AY036408)
# length(C_EAC$C.BI.2002.02BU_B1027.AM260260)
# length(C_EAC$`C.RW.2005.8222-PD.KC513726`)
# length(C_EAC$C.KE.1991.KNH1268.AY945738)
# length(C_EAC$C.TZ.1996.6567.AY036408)
# length(C_EAC$C.CD.2002.02CD_KSS055.AM041050)
# 
# 
# length(D_EAC$D.TZ.1996.2975.AY036303)
# length(D_EAC$D.BI.2002.02BU_U1916.AM260230)
# length(D_EAC$`D.RW.2005.105-7101.JQ364514`)
# length(D_EAC$D.KE.1996.96KE006.AY492783)
# length(D_EAC$D.TZ.1996.2975.AY036303)
# length(D_EAC$D.CD.1983.ELI.K03454)




# Analysis for Subtype A -------


# Phylogenetic tree with FastTree

system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/EAC_HIV_1_A.fas > EAC_HIV_1_A.fas.nwk")


A_EAC_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/EAC_HIV_1_A.fas.nwk")



# Compute transmission clusters


work.dir <- getwd()

# run ClusterPicker which is in the working directory
# we take sequence file in data folder and newik tree file in the working directory where we saved it before

system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), 
                                 paste0(work.dir,"/data/", "EAC_HIV_1_A.fas"), paste0(work.dir,"/","EAC_HIV_1_A.fas.nwk"),  paste0("0.8 0.7 0.045 2 gap"))))

# Read clusters' files

dd <- list.files(path = paste0(work.dir), 
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



# Phylogenetic tree with FastTree

system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/EAC_HIV_1_C.fas > EAC_HIV_1_C.fas.nwk")


C_EAC_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/EAC_HIV_1_C.fas.nwk")



# Compute transmission clusters


work.dir <- getwd()


# run ClusterPicker which is in the working directory
# we take sequence file in data folder and newik tree file in the working directory where we saved it before

system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), 
                                 paste0(work.dir,"/data/", "EAC_HIV_1_C.fas"), paste0(work.dir,"/","EAC_HIV_1_C.fas.nwk"),  paste0("0.8 0.7 0.045 2 gap"))))

# Read clusters' files

dd <- list.files(path = paste0(work.dir), 
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


# Phylogenetic tree with FastTree

system("/home/david/Dropbox/eac_phylo_analysis/FastTree -gtr -nt < /home/david/Dropbox/eac_phylo_analysis/data/EAC_HIV_1_D.fas > EAC_HIV_1_D.fas.nwk")


D_EAC_tree.const <- read.tree("/home/david/Dropbox/eac_phylo_analysis/EAC_HIV_1_D.fas.nwk")



# Compute transmission clusters


work.dir <- getwd()
sub.dir.rename <- "/home/david/Dropbox/eac_phylo_analysis/results/Transmission Clusters"


# run ClusterPicker which is in the working directory
# we take sequence file in data folder and newik tree file in the working directory where we saved it before

system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), 
                                 paste0(work.dir,"/data/", "EAC_HIV_1_D.fas"), paste0(work.dir,"/","EAC_HIV_1_D.fas.nwk"),  paste0("0.8 0.7 0.045 2 gap"))))

# Read clusters' files

dd <- list.files(path = paste0(work.dir), 
                 pattern = paste0("EAC_HIV_1_D_EAC_HIV_1_A_clusterPicks_cluster"),
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




# Transmission clusters sizes --------------


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

