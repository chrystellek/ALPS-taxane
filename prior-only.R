# Using ALPS for taxane metabolic pathway
# Chrystelle Kiang
# Created June 10, 2024

# This is prior only run for prior probabilities 
# The prior only run results also used to create cpsi grid 
library(here)
library(tidyverse)

# setwd("..")
# import prior forest
prior_forest <- read_csv("./data/priorforest.csv", show_col_types = FALSE)
prior.forest <- subset(prior_forest, select = c("snp1","snp2"))

# import ALPS files
source(file = paste0("./ALPS/alps2.R"))

# TODO update with new dataset
# import scrambled dataset
simgeno <- read_csv(file = "./data/simpath_cohort_scramCK.csv", show_col_types = FALSE)
# only keeping snp names we are using
simgeno <- simgeno %>%
  select(-abc10, -c(abc7:abc9), -c(cyp10:cyp14), -c(cyp6:cyp9), -sul5, -ugt2, -ugt3, -sul2)

# setwd("output")
# set up run
ds <- simgeno %>% arrange(time)
time.var <- ds$time
status.var <- simgeno$event
dos <- ds %>%
  mutate(labid = NULL,
         time = NULL,
         event = NULL,
         X = NULL,
         xn = NULL)
dos <- t(as.matrix(dos))

# ALPS needs pointer to wd 
wd <- here()
set.seed(404)
# start at random spot
spot <- prior.forest[sample(nrow(prior.forest),1),]
spot.tree <- paste0("(", spot$snp1, ",", spot$snp2, ");")
curtree <- read.tree(text = spot.tree)

# TODO confirm if initializing psi as NA is an OK solve
# was getting error otherwise that object 'psi' not found
# alternatively was considering all 0's?: psi <- numeric(100)

psi <- rep(NA, 100)
# starting with low number to check feasibility 
system.time(fitalps(iter = 1000,
                    initipsi = 19,
                    curtree,
                    normpotts = TRUE,
                    prioronly = TRUE,
                    lik = "coxph-ties",
                    prefix="/output/prior"))

# code below is from cpsi.R by TPA (github)
norm <- read.table("./output/prior-tree-parameters.txt", header = TRUE, sep = "\t")
# summary(norm$distance) 

maxpsi <- 3
npsi <- 100  
# overwriting psi 
psi <- (maxpsi*1:npsi)/npsi
cpsi <- array(dim = npsi)
for (i in 1:npsi){
  cpsi[i] <- sum(exp(-psi[i]*norm$distance))/nrow(norm)
}

# export
save(cpsi, maxpsi, npsi, psi, file = "./data/cpsi_taxane.RData")
