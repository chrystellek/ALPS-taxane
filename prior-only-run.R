# Using ALPS for taxane metabolic pathway
# Chrystelle Kiang
# last updated October 16, 2024

# This is prior only run for prior probabilities 
# These results are used to create cpsi grid 
library(here)
library(tidyverse)

# import prior forest
# TODO change file location
prior.forest <- read.csv("./data/priorforeststructure.csv") %>% select("snp1","snp2")

# import ALPS files
# TODO change file location
source(file = paste0("./alps2.R"))

# import dataset
load("./taxane_imp_agg.RData")

#### retrofitting scram data for testing 
genodata <- read.csv("./data/simpath_cohort_scram_tax.csv") 
genodata <- genodata %>% 
  rename(id = cpr) %>%
  select(id, 
         abc2, abc3, abc4, abc11, abc12,
         cyp13, cyp14, cyp15, cyp16, cyp17, 
         gst1, slc1, slc2, slc3,
         event, time)
# TODO delete for actual run 
##########


# set up run - sort dataset by time
ds <- genodata %>% 
  arrange(time)
time.var <- ds$time
status.var <- ds$event
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

# note that I initialized psi as NA 
# was getting error that object 'psi' not found
psi <- rep(NA, 100)

# starting with low number to check feasibility 
# TODO does this number of iter have to = number of ALPS iter?  
system.time(fitalps(iter = 1000,
                    initipsi = 19,
                    curtree,
                    normpotts = TRUE,
                    prioronly = TRUE,
                    lik = "coxph-ties",
                    # TODO ensure this folder exists inside current wd
                    prefix="/output-test/prior"))

# code below adapted from cpsi.R by TPA
norm <- read.table("./output-test/prior-tree-parameters.txt", header = TRUE, sep = "\t")
summary(norm$distance)  # should not be NA

maxpsi <- 3
npsi <- 100  
# overwriting psi 
psi <- (maxpsi*1:npsi)/npsi
cpsi <- array(dim = npsi)
for (i in 1:npsi){
  cpsi[i] <- sum(exp(-psi[i]*norm$distance))/nrow(norm)
}

# export
# TODO file location 
save(cpsi, maxpsi, npsi, psi, file = "./output-test/cpsi_taxane.RData")
