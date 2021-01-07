
# Before running the code below, make sure you have followed the instructions for installing all the required software


# Set the working directory

setwd("/home/david/Dropbox/eac_phylo_analysis/simpact_simulation")

# setwd("/path/to/your/working_directory/") # Change this to your own working directory

# Make sure you put Seq-gen, FastTree, ClusterPicker_1.2.3.jar tools and root sequence in the working directory
#



# Make sure compiled tools (Seq-Gen and FastTree) are in same working directory

## Load required packages
library(RSimpactCyan)
library(RSimpactHelper)
library(phangorn)
library(treedater)
library(picante)
library(igraph)
library(geomnet)
library(lubridate)
library(ggtree)
library(ggplot2)
library(tidyr)
library(dplyr)

# Sourcing the IDs.Seq.Random.skew function and renaming IDS in transmission table 
# to match tips names in the tree



inputvector <- c(777, -0.52, -0.05, 2, 10, 5, 0.25, -0.3, -0.1,
                 -1, -90, 0.5, 0.05, -0.14, 5, 7, 12, -3.0) # -1.7 




# Step 1: Run Simpact ------------


## Run Simpact for specific parameter combination



age.distr <- agedistr.creator(shape = 5, scale = 65)

cfg.list <- input.params.creator(population.eyecap.fraction = 0.2,
                                 population.simtime = 40, 
                                 population.nummen = 5000, 
                                 population.numwomen = 5000,
                                 hivseed.time = 10, 
                                 hivseed.type = "amount",
                                 hivseed.amount = 10, 
                                 hivseed.age.min = 20,
                                 hivseed.age.max = 50,
                                 formation.hazard.agegapry.meanage = -0.025,
                                 debut.debutage = 15
)

# Assumption of nature of sexual network


cfg.list["population.msm"] = "no"


# Sexual behaviour

seedid <- inputvector[1]

cfg.list["dissolution.alpha_0"] <- inputvector[2] # [1] # -0.52 c("unif", -1, 0)
cfg.list["dissolution.alpha_4"] <- inputvector [3] # [2] # -0.05 c("unif", -0.5, 0)
cfg.list["formation.hazard.agegapry.baseline"] <- inputvector[4] # [3] # 2 c("unif", 1, 3)
cfg.list["person.agegap.man.dist.normal.mu"] <- inputvector[5] # [4] # 0 c("unif", -0.5, 0.5)
cfg.list["person.agegap.woman.dist.normal.mu"] <- inputvector[5] # [4] # 0
cfg.list["person.agegap.man.dist.normal.sigma"] <- inputvector[6] # [5] # 3 c("unif", 2, 4)
cfg.list["person.agegap.woman.dist.normal.sigma"] <- inputvector[6] # [5] # 3 
cfg.list["formation.hazard.agegapry.gap_agescale_man"] <- inputvector[7] # [6] # 0.25 c("unif", 0, 1)
cfg.list["formation.hazard.agegapry.gap_agescale_woman"] <- inputvector[7] # [6] # 0.25
cfg.list["formation.hazard.agegapry.numrel_man"] <- inputvector[8] # [7] # -0.3 c("unif", -1, 0)
cfg.list["formation.hazard.agegapry.numrel_woman"] <- inputvector[8] # [7] # -0.3
cfg.list["formation.hazard.agegapry.numrel_diff"] <- inputvector[9] # [8] # -0.1 c("unif", -0.9, 0)


# HIV transmission

cfg.list["hivtransmission.param.a"] <- inputvector[10] # [10] # -1 c("unif", -2, 0)
cfg.list["hivtransmission.param.b"] <- inputvector[11] # [11] # -90 c("unif", -100, -80)
cfg.list["hivtransmission.param.c"] <- inputvector[12] # [12] # 0.5 c("unif", 0, 1)
cfg.list["hivtransmission.param.f1"] <- inputvector[13] # [13] # 0.04879016 c("unif", 0, 0.5)
cfg.list["hivtransmission.param.f2"] <- inputvector[14] # [14] # -0.1386294 c("unif", -0.5, 0)

# Disease progression > may be remove in parameter to estimates

cfg.list["person.vsp.toacute.x"] <- inputvector[15] # [15] # 5 c("unif", 3, 7)
cfg.list["person.vsp.toaids.x"] <- inputvector[16] # [16] # 7 c("unif", 5, 9)
cfg.list["person.vsp.tofinalaids.x"] <- inputvector[17] # [17] # 12 c("unif", 10, 14)



# Demographic

cfg.list["conception.alpha_base"] <- inputvector[18] # [18] # -2.7 c("unif", -3.5, -1.7)


# Assumptions to avoid negative branch lengths

# + sampling == start ART
# when someone start ART, he/she is sampled and becomes non-infectious

cfg.list["monitoring.fraction.log_viralload"] <- 0 # very important for transmission tree and sequence simulation

# Note: If treatment is started, the person’s set-point viral load value will be lowered according 
# to the setting in monitoring.fraction.log_viralload, if we set it to 0, it means the person is no more infectious

#
# ## Add-ons
#
### BEGIN Add-on
cfg.list["mortality.aids.survtime.C"] <- 65
cfg.list["mortality.aids.survtime.k"] <- -0.2
cfg.list["dropout.interval.dist.type"] <- "uniform"
cfg.list["dropout.interval.dist.uniform.min"] <- 1000
cfg.list["dropout.interval.dist.uniform.max"] <- 2000

cfg.list["person.survtime.logoffset.dist.type"] <- "normal"
cfg.list["person.survtime.logoffset.dist.normal.mu"] <- 0
cfg.list["person.survtime.logoffset.dist.normal.sigma"] <- 0.1

cfg.list["person.agegap.man.dist.type"] <- "normal" #fixed
#cfg.list["person.agegap.man.dist.fixed.value"] <- -6
cfg.list["person.agegap.woman.dist.type"] <- "normal" #"fixed"
#cfg.list["person.agegap.woman.dist.fixed.value"] <- -6


cfg.list["person.eagerness.man.dist.gamma.a"] <- 0.23 # 0.23
cfg.list["person.eagerness.woman.dist.gamma.a"] <- 0.23 # 0.23
cfg.list["person.eagerness.man.dist.gamma.b"] <- 45 # 45
cfg.list["person.eagerness.woman.dist.gamma.b"] <- 45 # 45



cfg.list["monitoring.cd4.threshold"] <- 0 # 0 means nobody qualifies for ART
cfg.list["person.art.accept.threshold.dist.fixed.value"] <- 0.75
cfg.list["diagnosis.baseline"] <- -99999 # -2


#### END Add-ons


# # ART intervention



# Let's introduce ART, and evaluate whether the HIV prevalence drops less  rapidly

# Introducing ART
art.intro <- list()
art.intro["time"] <- 23    # ~2000
art.intro["diagnosis.baseline"] <- -2
art.intro["monitoring.cd4.threshold"] <- 100

art.intro1 <- list()
art.intro1["time"] <- 25     # ~2002
art.intro1["diagnosis.baseline"] <- -1.8
art.intro1["monitoring.cd4.threshold"] <- 150

art.intro2 <- list()
art.intro2["time"] <- 28     # ~2005
art.intro2["diagnosis.baseline"] <- -1.5
art.intro2["monitoring.cd4.threshold"] <- 200

art.intro3 <- list()
art.intro3["time"] <- 33     # ~2010
art.intro3["diagnosis.baseline"] <- -1
art.intro3["monitoring.cd4.threshold"] <- 350

art.intro4 <- list()
art.intro4["time"] <- 36     # ~2013
art.intro4["monitoring.cd4.threshold"] <- 500

art.intro5 <- list()
art.intro5["time"] <- 39     # ~2016
art.intro5["monitoring.cd4.threshold"] <- 700 # 


interventionlist <- list(art.intro, art.intro1, art.intro2, art.intro3, art.intro4, art.intro5)

intervention <- interventionlist


# Events
cfg.list["population.maxevents"] <- as.numeric(cfg.list["population.simtime"][1]) * as.numeric(cfg.list["population.nummen"][1]) * 3


#
# # # # Run Simpact
results <- simpact.run(configParams = cfg.list,
                       destDir = "temp",
                       agedist = age.distr,
                       seed = seedid,
                       intervention = intervention)

datalist <- readthedata(results)



# save(datalist, file="/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/datalist.RData")

datalist <- get(load("/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/datalist.RData"))


# Epidemic characteristics: prevalence


# Women

hiv.prev.lt25.women <- prevalence.calculator(datalist = datalist,
                                             agegroup = c(15, 25),
                                             timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()

hiv.prev.25.30.women <- prevalence.calculator(datalist = datalist,
                                              agegroup = c(25, 30),
                                              timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()

hiv.prev.30.35.women <- prevalence.calculator(datalist = datalist,
                                              agegroup = c(30, 35),
                                              timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()

hiv.prev.35.40.women <- prevalence.calculator(datalist = datalist,
                                              agegroup = c(35, 40),
                                              timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()

hiv.prev.40.45.women <- prevalence.calculator(datalist = datalist,
                                              agegroup = c(40, 45),
                                              timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()

hiv.prev.45.50.women <- prevalence.calculator(datalist = datalist,
                                              agegroup = c(45, 50),
                                              timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(2) %>%
  as.numeric()



# Men

hiv.prev.lt25.men <- prevalence.calculator(datalist = datalist,
                                           agegroup = c(15, 25),
                                           timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()

hiv.prev.25.30.men <- prevalence.calculator(datalist = datalist,
                                            agegroup = c(25, 30),
                                            timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()

hiv.prev.30.35.men <- prevalence.calculator(datalist = datalist,
                                            agegroup = c(30, 35),
                                            timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()

hiv.prev.35.40.men <- prevalence.calculator(datalist = datalist,
                                            agegroup = c(35, 40),
                                            timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()

hiv.prev.40.45.men <- prevalence.calculator(datalist = datalist,
                                            agegroup = c(40, 45),
                                            timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()

hiv.prev.45.50.men <- prevalence.calculator(datalist = datalist,
                                            agegroup = c(45, 50),
                                            timepoint = 40) %>%
  dplyr::select(pointprevalence) %>%
  dplyr::slice(1) %>%
  as.numeric()





agegroup <- c("15-24", "25-29", "30-34", "35-39", "40-44", "45-49")

men.prev <- c(hiv.prev.lt25.men, hiv.prev.25.30.men, hiv.prev.30.35.men, hiv.prev.35.40.men, hiv.prev.40.45.men, hiv.prev.45.50.men)
women.prev <- c(hiv.prev.lt25.women, hiv.prev.25.30.women, hiv.prev.30.35.women, hiv.prev.35.40.women, hiv.prev.40.45.women, hiv.prev.45.50.women)


men.prev.dat <- data.frame(agegroup, men.prev)
men.prev.dat$gender <- "men"
names(men.prev.dat) <- c("age", "prev", "gender")

women.prev.dat <- data.frame(agegroup, women.prev)
women.prev.dat$gender <- "women"
names(women.prev.dat) <- c("age", "prev", "gender")

prev_data <- rbind(men.prev.dat, women.prev.dat)

plot.prev.men.women <- ggplot(prev_data, aes(x=age, y=prev, colour=gender, group = gender)) + 
  # geom_errorbar(aes(ymin=L, ymax=U), width=.1) +
  geom_line(size=.3) +
  geom_point() + 
  xlab("Age groups") + ylab("HIV prevalence")




# Step 2: Construct transmission networks ----------


# Produce a list of transmission networks in epi object format

simpact.trans.net <- transmission.network.builder(datalist = datalist, endpoint = 40)


# save(simpact.trans.net, file="/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/simpact.trans.net.RData")

simpact.trans.net <- get(load("/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/simpact.trans.net.RData"))


# (i) Transmission network


# Entire transmission network 

trans.net <- simpact.trans.net

trans.net <- data.table::rbindlist(simpact.trans.net)

trans.net <- dplyr::filter(trans.net, trans.net$parent != "-1") # remove universal seed individuals


graph.build <- as.data.frame(trans.net)

graph.build[,4] <- as.character(graph.build$parent) # donors
graph.build[,3] <- as.character(graph.build$id) # recipients
gag = as.matrix(graph.build)
gag = gag[-1,] # remove universall seed -1
ga.graph = graph.edgelist(gag[,4:3])

V(ga.graph)$color <- "red"

transmi.network.full <- ga.graph





# Step 3: Sequence simulation ------------

# --- DONE ONCE
# with polimerase fragment and with protease fragment

dirseqgen <- "/home/david/Dropbox/eac_phylo_analysis/simpact_simulation"
dirfasttree <- "/home/david/Dropbox/eac_phylo_analysis/simpact_simulation" 
sub.dir.rename <- "/home/david/Dropbox/eac_phylo_analysis/simpact_simulation"




sequence.simulation.seqgen.par(dir.seq = dirseqgen,
                               sub.dir.rename = sub.dir.rename,
                               simpact.trans.net = simpact.trans.net, # simpact.trans.net,
                               seq.gen.tool = "seq-gen",
                               seeds.num = 777,
                               endpoint = 40,
                               limitTransmEvents = 7, # no less than 7
                               hiv.seq.file = "hiv.seq.C.pol.j.fasta", # hiv.seq.C.protease.j.fasta  # hiv.seq.C.pol.j.fasta
                               clust = FALSE) # hiv.seq.file lodged in work.dir

# Output sequence file: C.Epidemic_seed.seq.bis.sim.nwk.fasta 

# Transform the sequence format to be handled by ClusterPicker
sequ.dna <- read.dna(file = paste0(sub.dir.rename,"/C.Epidemic_seed.seq.bis.sim.nwk.fasta"), format = "interleaved")

write.dna(sequ.dna, file = paste0(sub.dir.rename,"/C.Epidemic.pol.fas") , format = "fasta") # C.Epidemic.pol.fas

# write.dna(sequ.dna, file = paste0(sub.dir.rename,"/C.Epidemic.protease.fas") , format = "fasta") # C.Epidemic.protease.fas



# Step 4: Construct time stamped phylogenetic trees  ----------


# Tree with polimerase -----------

dirfasttree <- "/home/david/Dropbox/eac_phylo_analysis/simpact_simulation"
sub.dir.rename <- "/home/david/Dropbox/eac_phylo_analysis/simpact_simulation"

dir.tree <- dirfasttree
sub.dir.rename <- sub.dir.rename
fasttree.tool <- "FastTree"
calendar.dates <- "samplingtimes.all.csv"
simseqfile <- "C.Epidemic.pol.fas" # C.Epidemic.protease.fas  # C.Epidemic.pol.fas
count.start <- 1977
endsim <- 40


dates.Transform.NamedVector  <- function(dates=dates){
  
  dates.val <- endsim - dates$V2 + count.start # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  
  return(dates.val)
}



out.fast.tree.file <- paste0(sub.dir.rename,"/", simseqfile ,".nwk")


system(paste0("./", paste0(paste(fasttree.tool, "-gtr")," -nt < ", paste0(sub.dir.rename,"/",simseqfile), paste0("> ", out.fast.tree.file))))


samp.dates <- read.csv(paste0(sub.dir.rename, "/", calendar.dates))

tree.const <- read.tree(out.fast.tree.file)

time.samp <- dates.Transform.NamedVector(dates=samp.dates) # name the dates


# Match dates and phylogenetic tree leaves 

tree.tips <- as.character(tree.const$tip.label)

ord.tree.dates <- vector()

for(i in 1:length(tree.tips)){
  for(j in 1:length(time.samp)){
    if(tree.tips[i] == samp.dates$V1[j]){
      ord.tree.dates <- c(ord.tree.dates, time.samp[j])
    }
  }
}


dater.tree.pol <- treedater::dater(tree.const, 
                                   ord.tree.dates, 
                                   s = 3000) # C.Epidemic.pol.fas.nwk



rootToTipRegressionPlot(dater.tree.pol) # regression (figure 8 * 6)


editedRootToTipRegressionPlot(dater.tree.pol) # to add summary in the plot


# save(dater.tree.pol, file="/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/dater.tree.pol.RData")

dater.tree.pol <- get(load("/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/dater.tree.pol.RData"))


tree.calib <- dater.tree.pol



# Node age with picante package

N <- node.age(tree.calib)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti


latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date


# Entire transmission network 

trans.net <- simpact.trans.net

trans.net <- data.table::rbindlist(simpact.trans.net)

trans.net <- dplyr::filter(trans.net, trans.net$parent != "-1") # remove universal seed individuals

# We considered the epidemic to start at 10 and individuals infected up to 40 simulation time: 30 years of epidemic
# With a seed sequence sampled in 1989, we assume it existed two years before (1987)
# It means that the simulation started in 1977, and the infection in 1987 for 30 years in 2017

trans.net$dtimes <- abs(trans.net$dtimes-40)+1977
trans.net$itimes <- abs(trans.net$itimes-40)+1977

min.val = 1977
max.val = round(max(trans.net$itimes))



# (i) Number of internal nodes & transmission events

step.int=1 # in one year

d <- (max.val-min.val)/step.int

dat.f.trans <- as.data.frame(trans.net)
dt.node.age.dt <- int.node.age

num.tra <- vector() # initialize transmission events
i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  inf <- 1976+i
  sup <- 1977+i
  dat.f.trans.i <- dat.f.trans[which(dat.f.trans$itimes <= sup & dat.f.trans$itimes  > inf),]
  num.i <- nrow(dat.f.trans.i)
  num.tra <- c(num.tra, num.i)
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
}


int.node.pol.vec <- int.node.vec



# Tree with protease ------------



dirfasttree <- "/home/david/Dropbox/eac_phylo_analysis/simpact_simulation"
sub.dir.rename <- "/home/david/Dropbox/eac_phylo_analysis/simpact_simulation"

dir.tree <- dirfasttree
sub.dir.rename <- sub.dir.rename
fasttree.tool <- "FastTree"
calendar.dates <- "samplingtimes.all.csv"
simseqfile <- "C.Epidemic.protease.fas" # C.Epidemic.protease.fas  # C.Epidemic.pol.fas
count.start <- 1977
endsim <- 40


dates.Transform.NamedVector  <- function(dates=dates){
  
  dates.val <- endsim - dates$V2 + count.start # dates
  names(dates.val) <- as.character(dates$V1) # names are the names of the tips
  
  return(dates.val)
}



out.fast.tree.file <- paste0(sub.dir.rename,"/", simseqfile ,".nwk")


system(paste0("./", paste0(paste(fasttree.tool, "-gtr")," -nt < ", paste0(sub.dir.rename,"/",simseqfile), paste0("> ", out.fast.tree.file))))


samp.dates <- read.csv(paste0(sub.dir.rename, "/", calendar.dates))

tree.const <- read.tree(out.fast.tree.file)

time.samp <- dates.Transform.NamedVector(dates=samp.dates) # name the dates


# Match dates and phylogenetic tree leaves 

tree.tips <- as.character(tree.const$tip.label)

ord.tree.dates <- vector()

for(i in 1:length(tree.tips)){
  for(j in 1:length(time.samp)){
    if(tree.tips[i] == samp.dates$V1[j]){
      ord.tree.dates <- c(ord.tree.dates, time.samp[j])
    }
  }
}


dater.tree.prot <- treedater::dater(tree.const, 
                                    ord.tree.dates, 
                                    s = 297) # C.Epidemic.protease.fas.nwk


fixedRootToTipRegressionPlot(dater.tree.prot)


editedRootToTipRegressionPlot(dater.tree.prot) # to add summary in the plot


# save(dater.tree.prot, file="/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/dater.tree.prot.RData")

dater.tree.prot <- get(load("/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/dater.tree.prot.RData"))


tree.calib <- dater.tree.prot



# Node age with picante package

N <- node.age(tree.calib)

# Time to MRCA: internal nodes ages

int.node.age <- N$Ti


latest.samp <- N$timeToMRCA+N$timeOfMRCA # latest sampling date


# Entire transmission network 

trans.net <- simpact.trans.net

trans.net <- data.table::rbindlist(simpact.trans.net)

trans.net <- dplyr::filter(trans.net, trans.net$parent != "-1") # remove universal seed individuals

# We considered the epidemic to start at 10 and individuals infected up to 40 simulation time: 30 years of epidemic
# With a seed sequence sampled in 1989, we assume it existed two years before (1987)
# It means that the simulation started in 1977, and the infection in 1987 for 30 years in 2017

trans.net$dtimes <- abs(trans.net$dtimes-40)+1977
trans.net$itimes <- abs(trans.net$itimes-40)+1977

min.val = 1977
max.val = round(max(trans.net$itimes))



# (i) Number of internal nodes & transmission events

step.int=1 # in one year

d <- (max.val-min.val)/step.int

dat.f.trans <- as.data.frame(trans.net)
dt.node.age.dt <- int.node.age

num.tra <- vector() # initialize transmission events
i.vec <- vector() # initialise time intervals
int.node.vec <- vector() # initialize internal nodes

for (i in 1:d) {
  inf <- 1976+i
  sup <- 1977+i
  dat.f.trans.i <- dat.f.trans[which(dat.f.trans$itimes <= sup & dat.f.trans$itimes  > inf),]
  num.i <- nrow(dat.f.trans.i)
  num.tra <- c(num.tra, num.i)
  i.vec <- c(i.vec, sup)
  int.node.age.i <- int.node.age[int.node.age <= sup & dt.node.age.dt > inf]
  int.node.vec <- c(int.node.vec,length(int.node.age.i))
}


int.node.prot.vec <- int.node.vec



# All outputs ------------

SimpactSimulation <- list()
SimpactSimulation$dater.tree.pol <- dater.tree.pol
SimpactSimulation$dater.tree.prot <- dater.tree.prot
SimpactSimulation$interv.year <- i.vec
SimpactSimulation$int.node.pol.vec <- int.node.pol.vec
SimpactSimulation$int.node.prot.vec <- int.node.prot.vec
SimpactSimulation$num.transm <- num.tra
SimpactSimulation$transmi.network.full <- transmi.network.full


# save(SimpactSimulation, file = "/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/SimpactSimulation.RData")

SimpactSimulation <- get(load("/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/SimpactSimulation.RData"))


# Visualise transmission events versus internal nodes -------------

calendaryear.isrelevant <- SimpactSimulation$interv.year >= 1987
calendaryear <- SimpactSimulation$interv.year[calendaryear.isrelevant]
intern.nodes.pol <- SimpactSimulation$int.node.pol.vec[calendaryear.isrelevant]
intern.nodes.prot <- SimpactSimulation$int.node.prot.vec[calendaryear.isrelevant]
trans.events <- SimpactSimulation$num.transm[calendaryear.isrelevant]

trans.and.nodes.df <- data.frame(calendaryear = calendaryear,
                                 intern.nodes.pol = intern.nodes.pol,
                                 intern.nodes.prot = intern.nodes.prot,
                                 trans.events = trans.events)

trans.and.nodes.df$intern.nodes.pol <- trans.and.nodes.df$intern.nodes.pol/sum(trans.and.nodes.df$intern.nodes.pol)
trans.and.nodes.df$intern.nodes.prot <- trans.and.nodes.df$intern.nodes.prot/sum(trans.and.nodes.df$intern.nodes.prot)
trans.and.nodes.df$trans.events <- trans.and.nodes.df$trans.events/sum(trans.and.nodes.df$trans.events)


trans.and.nodes.long.df <- gather(trans.and.nodes.df,
                                  key = "Events",
                                  value = "Proportion",
                                  intern.nodes.pol, intern.nodes.prot, trans.events,
                                  factor_key = TRUE)


transandnodes.plot <- ggplot(data = trans.and.nodes.long.df, # trans.and.nodes.long.df,
                             aes(x = calendaryear,
                                 y = Proportion,
                                 colour = factor(Events))) +
  geom_point() +
  scale_color_brewer(palette="Set1",
                     name = "",
                     labels = c("Internal nodes pol",
                                "Internal nodes prot",
                                "Transmission events")) +
  geom_line() +
  theme(axis.line.x = element_line(),
        legend.position=c(0.75, 0.9),
        legend.key = element_blank(),
        legend.background = element_blank()) +
  scale_x_continuous(limits = c(1985, 2020),
                     breaks = seq(from = 1985,
                                  to = 2020,
                                  by = 5)) +
  xlab("Time")
print(transandnodes.plot)


ggsave(filename = "internalnodels_transmission_events.pdf",
       plot = transandnodes.plot,
       path = "/home/david/Dropbox/eac_phylo_analysis/simpact_simulation",
       width = 18, height = 10, units = "cm")



# Compute error between transmission events and internal nodes --------

intern.nodes.pol <- SimpactSimulation$int.node.pol.vec[calendaryear.isrelevant]
intern.nodes.prot <- SimpactSimulation$int.node.prot.vec[calendaryear.isrelevant]
trans.events <- SimpactSimulation$num.transm[calendaryear.isrelevant]

error_transm_nodes_pol <- Metrics::rmse(trans.events, intern.nodes.pol)

error_transm_nodes_prot <- Metrics::rmse(trans.events, intern.nodes.prot)






# Visualize the phylogenetic tree -------------


# Tree with polimerase ---------


dater.tree.pol <- SimpactSimulation$dater.tree.pol
  
tree.calib <- dater.tree.pol

class(tree.calib) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
tree.calib_sim.start.year <- round(tree.calib$timeOfMRCA) # akantu dater.tree$timeToMRCA
tree.calib_first.transmission <- min(tree.calib$sts)
mrsd <- max(tree.calib$sts)

dates <- format(date_decimal(c(mrsd, tree.calib_sim.start.year)), "%Y-%m-%d")
tree.calib$root.edge <-  - tree.calib_sim.start.year


phylotree.pol.plot <- ggtree(tree.calib,
                         mrsd = dates[1],
                         size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(tree.calib_sim.start.year-1, mrsd+1),
                     breaks = seq(from = tree.calib_sim.start.year-1,
                                  to = mrsd+1,
                                  by = 2)) +
  xlab("Time") 


print(phylotree.pol.plot) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/phylotree.pol.plot.pdf",
       phylotree.pol.plot,
       width = 25,
       height = 15,
       units = "cm")






# Tree with protease ---------


dater.tree.prot <- SimpactSimulation$dater.tree.prot

tree.calib <- dater.tree.prot


class(tree.calib) <- "phylo" # Removing "treedater" as one of the classes that this object belongs to.
tree.calib_sim.start.year <- round(tree.calib$timeOfMRCA) # akantu dater.tree$timeToMRCA
tree.calib_first.transmission <- min(tree.calib$sts)
mrsd <- max(tree.calib$sts)

dates <- format(date_decimal(c(mrsd, tree.calib_sim.start.year)), "%Y-%m-%d")
tree.calib$root.edge <-  - tree.calib_sim.start.year


phylotree.prot.plot <- ggtree(tree.calib,
                         mrsd = dates[1],
                         size = 0.2) + 
  theme_tree2() +
  theme_grey() +
  theme(axis.text.x = element_text(angle =- 90, vjust = 0.5),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(limits = c(tree.calib_sim.start.year-1, mrsd+1),
                     breaks = seq(from = tree.calib_sim.start.year-1,
                                  to = mrsd+1,
                                  by = 2)) +
  xlab("Time") 


print(phylotree.prot.plot) # more @ https://4va.github.io/biodatasci/r-ggtree.html


ggsave("/home/david/Dropbox/eac_phylo_analysis/simpact_simulation/phylotree.prot.plot.pdf",
       phylotree.prot.plot,
       width = 25,
       height = 15,
       units = "cm")



# Compare trees  -----


tree_polimerase <- dater.tree.pol$intree
tree_protease <- dater.tree.prot$intree

ape::all.equal.phylo(tree_polimerase, tree_protease)

identical(tree_polimerase, tree_protease)


ape::comparePhylo(tree_polimerase, tree_protease, plot = TRUE, force.rooted = FALSE,
           use.edge.length = FALSE)

# => Comparing dater.tree.pol$intree with dater.tree.prot$intree.
# Both trees have the same number of tips: 810.
# Both trees have the same tip labels.
# Both trees have the same number of nodes: 809.
# Both trees are rooted.
# Both trees are not ultrametric.
# 475 clades in dater.tree.pol$intree not in dater.tree.prot$intree.
# 475 clades in dater.tree.prot$intree not in dater.tree.pol$intree.
# Node labels of clades in common between both trees: see ..$NODES
# (node number in parentheses).


dist.topo(raw_tree_protease, raw_tree_polimerase, method = "PH85")


dist.topo(tree_polimerase, tree_protease, method = "score")


cherry(tree_polimerase)

cherry(tree_protease)



# Compare transmission clusters ---------


# Transmission clusters with polimerase ---------

work.dir <- getwd()

# run ClusterPicker which is in the working directory
# we take sequence file in data folder and newik tree file in the working directory where we saved it before

system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), 
                                 paste0(work.dir, "/C.Epidemic.pol.fas"), 
                                 paste0(work.dir, "/C.Epidemic.pol.fas.nwk"),  
                                 paste0("0.9 0.9 4.5 5 gap"))))



pol <- list.files(path = paste0(work.dir, "/"),
                 pattern = paste0("C.Epidemic.pol_C.Epidemic.pol_cluster"),
                 all.files = FALSE,
                 full.names = FALSE, recursive = FALSE) 

d <- pol

clust.size_pol <- vector() # size of each cluster # table.simpact.trans.net.adv

name_clust.size_pol <- vector()

tips_names_pol <- list()

for (i in 1:length(d)) {
  
  transm.df.cl.dat <- NULL
  
  clus.read <- read.table(file = paste0(paste0(work.dir, "/"),d[i]), header = FALSE) # Ids of each cluster
  
  clust.size_pol <- c(clust.size_pol, nrow(clus.read))
  
  name_clust.size_pol <- c(name_clust.size_pol, d[i])
  
  tips_names_pol[[i]] <- as.character(clus.read[,1])
  
}



# Transmission clusters with protease ---------


work.dir <- getwd()

# run ClusterPicker which is in the working directory
# we take sequence file in data folder and newik tree file in the working directory where we saved it before

system(paste("java -jar ", paste(paste0(work.dir,"/ClusterPicker_1.2.3.jar"), 
                                 paste0(work.dir, "/C.Epidemic.protease.fas"), 
                                 paste0(work.dir, "/C.Epidemic.protease.fas.nwk"),  
                                 paste0("0.9 0.9 4.5 5 gap"))))


prot <- list.files(path = paste0(work.dir, "/"),
                  pattern = paste0("C.Epidemic.protease_C.Epidemic.protease_cluster"),
                  all.files = FALSE,
                  full.names = FALSE, recursive = FALSE) 

d <- prot

clust.size_prot <- vector() # size of each cluster # table.simpact.trans.net.adv

name_clust.size_prot <- vector()

tips_names_prot <- list()

for (i in 1:length(d)) {
  
  transm.df.cl.dat <- NULL
  
  clus.read <- read.table(file = paste0(paste0(work.dir, "/"),d[i]), header = FALSE) # Ids of each cluster
  
  clust.size_prot <- c(clust.size_prot, nrow(clus.read))
  
  name_clust.size_prot <- c(name_clust.size_prot, d[i])
  
  tips_names_prot[[i]] <- as.character(clus.read[,1])
  
}





# Compare IDs in the bigest transmission clusters in polimerase and protease -------


clust.size_pol

tips_names_pol

ids_pol <- which(clust.size_pol >= 300)

seq_pol3 <- sort(tips_names_pol[[3]]) # ids_pol = 3



# %

clust.size_prot

tips_names_prot

ids_prot <- which(clust.size_prot >= 300)

seq_prot3 <- sort(tips_names_prot[[3]]) # ids_prot = 3


intersect_3_3 <- intersect(seq_pol3, seq_prot3) # length(seq_pol3)=324, length(seq_prot3)=320
length(intersect_3_3) # = 317





