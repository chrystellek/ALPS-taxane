# Using ALPS for taxane metabolic pathway
# Chrystelle Kiang
# last updated 6 November 2024

# This is prior only run for prior probabilities 
# These results are used to create cpsi grid 
library(here)
library(tidyverse)

# Files: there are some that need to be read in, ideally in same folder as this Rproject and/or script. tagged with TODO if there are changes that need to be made

# import prior forest
# TODO can change file location
prior.forest <- read.csv("./priorforeststructure.csv") %>% 
  select("snp1","snp2")

# import ALPS files
# TODO file location - I imagine you might already have this saved
source(file = paste0("./alps2.R"))

# import dataset created from missing imputation 
# TODO file location 
#imp_agg <- read.csv("./data/simpath_cohort_scram_tax.csv")
imp_agg <- readRDS("./data/genodata.rds")
# set up run - sort dataset by time
ds <- imp_agg %>% 
  # only keeping cpr, genes, event, and time 
  select(cpr, abc2, abc3, abc4, abc11, abc12, 
         cyp13, cyp14, cyp15, cyp16, cyp17, 
         gst1, slc1, slc2, slc3, event, time) %>%
  arrange(time)
time.var <- ds$time
status.var <- ds$event
dos <- ds %>%
  mutate(cpr = NULL,
         time = NULL,
         event = NULL,
         X = NULL,
         xn = NULL)
dos <- t(as.matrix(dos))

# ALPS needs pointer to wd 
wd <- here()
# TODO if you have trouble or prefer to save output in another folder change below
# wd <- "H:/home"
set.seed(4)
# start at random spot
spot <- prior.forest[sample(nrow(prior.forest),1),]
spot.tree <- paste0("(", spot$snp1, ",", spot$snp2, ");")
curtree <- read.tree(text = spot.tree)

# note that I initialized psi as NA 
# was getting error that object 'psi' not found
psi <- rep(1, 1000)
start_time <- Sys.time()
# This should create 4 .txt files that start with "prior-"
system.time(fitalps(iter = 1000,
                    initipsi = 19,
                    curtree,
                    normpotts = TRUE,
                    prioronly = TRUE,
                    lik = "coxph-ties",
            # if you want them in a new folder, can change to "/foldername/prior"
            # might need to create folder before running this 
                    prefix="/prior"))
# NOTE: moved to output folder since running this 
end_time <- Sys.time()
time_elapsed <- end_time - start_time


# code below adapted from cpsi.R by TPA
# TODO file location here should match above in prefix
norm <- read.table("./prior-tree-parameters.txt", header = TRUE, sep = "\t")
summary(norm$distance)  
# should not be NA

maxpsi <- 3
npsi <- 100  
# overwriting psi 
psi <- (maxpsi*1:npsi)/npsi
cpsi <- array(dim = npsi)
for (i in 1:npsi){
  cpsi[i] <- sum(exp(-psi[i]*norm$distance))/nrow(norm)
}
new_cpsi <- cpsi

# export
# TODO file location 
# save(cpsi, maxpsi, npsi, psi, file = "./cpsi_taxane.RData")
######################
# save to compare later
cpsi_100 <- cpsi
maxpsi_100 <- maxpsi
npsi_100 <- npsi
psi_100 <- psi
save(cpsi_100, maxpsi_100, npsi_100, psi_100, file = "./cpsi_taxane_old.RData")
load("./data/cpsi_taxane.RData")

cpsi_comp <- data.frame(matrix(nrow = 100, ncol = 2))
cpsi_comp$new_cpsi <- new_cpsi
cpsi_comp$old_cpsi <- cpsi
identical(cpsi_comp$new_cpsi, cpsi_comp$old_cpsi)
