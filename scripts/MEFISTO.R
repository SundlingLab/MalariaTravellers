## Script used to train MEFISTO model
## by Maximilian Julius Lautenbach

## libraries
library(MOFA2)
library(tidyverse)

## load data
mefisto.input <- readRDS("data/MEFISTO_input_data.Rds")
sampleTable <- read_csv("data/TravellerCohort_SampleTable.csv")

## create mofa object
MOFAobject <- create_mofa(data = mefisto.input %>% select(sample,feature,value,group,view))# %>% mutate_if(is.character, stri_trans_general,id = "latin-ascii") )

samples_metadata(MOFAobject) <- samples_metadata(MOFAobject) %>% 
  left_join(sampleTable %>% 
              select(sampleID,weeks.po, days.po, days.po.log) %>%
              rename(sample = sampleID), by="sample") %>% 
  distinct() %>% 
  filter(!duplicated(sample))
#%>% 
#  filter(sample !="2016009|Y1")

## set "weeks post" onset to as covariate
MOFAobject <- set_covariates(MOFAobject, covariates = "weeks.po" )

## get data options
data_opts <- get_default_data_options(MOFAobject)

## set model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 3

## set train options
train_opts <- get_default_training_options(MOFAobject)
train_opts$drop_factor_threshold <- 0.05
train_opts$maxiter <- 10000
train_opts$convergence_mode <-"slow" # "slow" # use "fast" for faster training

## set mefisto options
mefisto_opts <- get_default_mefisto_options(MOFAobject)
mefisto_opts$warping <- TRUE
mefisto_opts$warping_ref <- "primary_infected"
mefisto_opts$new_values <- seq(0, 70, 1)

## prepare mefisto object
MEFISTO.untrained <- prepare_mofa(MOFAobject, 
                                  model_options = model_opts,
                                  mefisto_options = mefisto_opts,
                                  training_options = train_opts,
                                  data_options = data_opts)

## train and save mefisto object
outfile = file.path("data/MEFISTO_model_updated.hdf5")
MEFISTO.trained <- run_mofa(MEFISTO.untrained, outfile, use_basilisk = T)


