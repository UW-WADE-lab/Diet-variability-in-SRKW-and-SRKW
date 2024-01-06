###Amplification bias correction
###SRKW and ARKW fecal samples
###Amy Van Cise, adapted from Shelton et al. 2022 (Towards Quantitative Metabarcoding)
###Includes edits to allow species in environmental samples to = 0,min, or max alpha (amplification efficiency)

library(dplyr)
library(compositions)
library(tidyverse)
library(rstan)
library(ggsci)

### Read and format input files ----------------------------------------------------------

#1. mock_input_data: longform dataframe with input proportions of each species in each mixture. 1 column for each species
#1 row for each mixture replicate. first three columns, before species, specify the mixutre "community" ID,
#the replicate number, and the number of PCR cycles.

#2. mock_seq_data: longform dataframe same as above with read count per species for each mixture.

#3. sample_seq_data: dataframe same as above with read count per species for each unknown sample.

#x <- gsub("\\\\", "/", readClipboard())
setwd("G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding")

## Import and format mock mixture input data
mock_input_data <- read.csv("Mockmix_control_proportion.csv") %>% 
  select(-X) %>% 
  pivot_longer(-c(community,tech_rep,Cycles), names_to = "Species", values_to = "true_prop")

mock_seq_data <- read.csv("Mockmix_control_reads.csv") %>% 
  select(-X) %>% 
  pivot_longer(-c(community,tech_rep,Cycles), names_to = "Species", values_to = "nReads") %>% 
  filter(Species %in% mock_input_data$Species) #removes species not in original input

#combine mock_input and mock_seq
# get rid of errors in the species that are zero
mock_dat <- full_join(mock_input_data,mock_seq_data) %>% 
              mutate(nReads = ifelse(true_prop==0,0,nReads))


## Import and format unknown mixture sample data
sample_seq_data <- read.csv("KWfecal_techrep_reads_allspecies_noreps.csv") %>% 
  select(-MiSeq.run, -FieldID, -Pop, -year) %>% 
  separate(Sample, sep = "_", into = c("community",NA)) %>% 
  mutate(community = str_remove(community, paste0(c("rep", "b", "B", "-dup"), collapse = "|"))) %>% 
  group_by(community) %>% 
  mutate(tech_rep = row_number(), .after = community) %>% 
  ungroup() %>% 
  mutate(Cycles = 32, .after = tech_rep) %>% 
  pivot_longer(-c(community,tech_rep,Cycles), names_to = "Species", values_to = "nReads") %>% 
  #filter(Species %in% mock_input_data$Species) %>% #removes species not in original input
  group_by(Species) %>% 
  mutate(tot.reads.species = sum(nReads)) %>% 
  filter(tot.reads.species > 0) %>% 
  select(-tot.reads.species) %>% 
  ungroup() %>% 
  group_by(community,tech_rep,Cycles) %>% 
  mutate(total.reads = sum(nReads)) %>% 
  ungroup() %>% 
  group_by(community) %>% 
  filter(total.reads == max(total.reads)) %>%  #for replicated samples, keep replicate with highest number of reads
  select(-total.reads) 
  

# Make a simple data structure that will only be used for making posterior predictions 
sample_data_small <- sample_seq_data #%>% filter(tech_rep==1)

#assign most common species to be the reference species
mostCommon <- sample_seq_data %>% 
  group_by(Species) %>% 
  tally(nReads > 0) %>%
  arrange(desc(n)) %>% 
  head(1) %>% 
  pull(Species)

sample_seq_data$Species[sample_seq_data$Species == mostCommon] <- paste0("zRefSpecies_", mostCommon) #make the most common species the reference species

sample_seq_data <- sample_seq_data %>% 
  arrange(Species, community) %>% 
  group_by(Species)

mock_dat$Species[mock_dat$Species == mostCommon] <- paste0("zRefSpecies_", mostCommon)

# limit sample species to those in that are non-zero in mock community data set 
mock_sp <- mock_dat %>% filter(true_prop>0) %>% distinct(Species) %>% pull(Species)

# also make a data frame with species that are not in the mock community data set, for later use
sample_seq_data_nomock <- sample_seq_data %>% 
  filter(!(Species %in% mock_sp)) %>% 
  arrange(Species, community)

sample_seq_data <- sample_seq_data %>%
  filter(Species %in% mock_sp)%>% #limit to species occurring in mock community data set
  arrange(Species, community)  

# Make a single species list:
        #sp.list   <- data.frame(Species = sort(unique(mock_dat$Species)) ) %>% mutate(sp_idx =1:length(Species))
#sp.list <- sample_seq_data %>% select(Species, sp_idx) %>% distinct(.keep_all = TRUE)
sp.list <- sort(unique(c(mock_dat$Species, sample_seq_data$Species, sample_seq_data_nomock$Species)) ) %>% 
  as.data.frame() %>% 
  rename("Species" = 1) %>% 
  mutate(sp_idx = 1:length(Species))

sp.list.mock <- sp.list %>% 
  filter(Species %in% mock_dat$Species)

sp.list.nomock <- sp.list %>% 
  filter(!(Species %in% mock_dat$Species))
  
#N_species <- nrow(sp.list)
N_species <- length(mock_sp)

# Make a list of samples/communities in the mock data set and the environmental data set
comm.mock.list <- mock_dat %>% group_by(community, tech_rep,Cycles) %>% summarise(n=length(tech_rep)) %>%
  ungroup() %>% mutate(id=1:length(n))
comm.env.list   <- sample_seq_data %>% group_by(community, tech_rep,Cycles) %>% summarise(n=length(tech_rep)) %>%
  ungroup() %>% mutate(id=1:length(n))

# make a list of species that are in mock community but not environment, 
# expand grid to make it so the the environmental samples get padded with all the
# missing species for all of the communities and all tech replicates.

sp.comm.mc  <- expand_grid(Species = sp.list.mock$Species, id = comm.mock.list$id) %>% 
  left_join(.,sp.list %>% select(Species,sp_idx)) %>%
  left_join(.,comm.mock.list %>% select(community,tech_rep,Cycles,id) ) %>% select(-id)

sp.comm.env <- expand_grid(Species = sp.list.mock$Species, id = comm.env.list$id) %>% 
  left_join(.,sp.list %>% select(Species,sp_idx)) %>%
  left_join(.,comm.env.list %>% select(community,tech_rep,Cycles,id) ) %>% select(-id)

# make a similar list of species that are NOT in the mock community but are in the environment
sp.comm.env.nomock <- expand_grid(Species = sp.list.nomock$Species, id = comm.env.list$id) %>% 
  left_join(.,sp.list %>% select(Species,sp_idx)) %>%
  left_join(.,comm.env.list %>% select(community,tech_rep,Cycles,id) ) %>% select(-id)

# convert to matrices
# merge in species and indices first to make pivoting more efficient.
mc  <- left_join(sp.comm.mc,mock_dat) %>%
  mutate(nReads = ifelse(is.na(nReads),0,nReads)) 

env <- left_join(sp.comm.env,sample_seq_data) %>%
  mutate(nReads = ifelse(is.na(nReads),0,nReads)) 

env.nomock <- left_join(sp.comm.env.nomock, sample_seq_data_nomock) %>%
  mutate(nReads = ifelse(is.na(nReads),0,nReads)) 

sample_data <- env %>% 
  ungroup() %>% 
  dplyr::select(community, sp_idx, nReads,tech_rep, Cycles) %>% 
  arrange(sp_idx) %>% 
  pivot_wider(names_from = "sp_idx", values_from = "nReads", values_fill = 0) 

sample_data_small <- sample_data #%>% filter(tech_rep==1)

sample_data_all <- bind_rows(env, env.nomock) %>% 
  ungroup() %>% 
  dplyr::select(community, sp_idx, nReads,tech_rep, Cycles) %>% 
  arrange(sp_idx) %>% 
  pivot_wider(names_from = "sp_idx", values_from = "nReads", values_fill = 0)

sample_data_all_prop <- sample_data_all %>% 
  pivot_longer(-c(community,tech_rep,Cycles), names_to = "sp_idx", values_to = "nReads") %>% 
  group_by(community) %>% 
  mutate(prop.reads = nReads/sum(nReads))

sample_data_all_small <- sample_data_all #%>% filter(tech_rep==1)

mock_data_fin <- mc %>% 
  ungroup() %>% 
  dplyr::select(community, sp_idx, nReads,tech_rep, Cycles) %>% 
  arrange(sp_idx) %>% 
  pivot_wider(names_from = "sp_idx", values_from = "nReads", values_fill = 0)

# add proportion data to mock community matrix

# proportions
p_mock <- mc %>%
  ungroup() %>% 
  dplyr::select(community, sp_idx, true_prop,tech_rep, Cycles) %>% 
  arrange(sp_idx) %>% 
  pivot_wider(names_from = "sp_idx", values_from = "true_prop", values_fill = 0)

# calculate additive log ratios 
alr_mock_true_prop <- p_mock[,4:(ncol(p_mock)-1)]*0
for(i in 1:nrow(p_mock)){
  alr_mock_true_prop[i,] <- alr(p_mock[i,4:(ncol(p_mock))] + 1e-12)
}
names.temp <- names(alr_mock_true_prop)
alr_mock_true_prop[,N_species] <- 0 #adding explicit reference species column

final.sp.name <- mc %>% 
                ungroup() %>% 
                filter(sp_idx == max(sp_idx)) %>% 
                slice(1) %>% 
                mutate(sp_idx = as.character(sp_idx)) %>% 
                pull(sp_idx)

names(alr_mock_true_prop) <- c(names.temp, final.sp.name)

rm(names.temp, final.sp.name)
### Create Design Matrices ---------------------------------------------------------------------------

## First, create design matrices for mock communities
# species compositions (betas)

N_pcr_mock <- mock_seq_data$Cycles

if(length(unique(mock_seq_data$community))==1){
  formula_b <- Cycles ~ 1  # what is on the left side of the equation doesn't matter.
} else {
  formula_b <- Cycles ~ community # what is on the left side of the equation doesn't matter.
}
model_frame <- model.frame(formula_b, mock_seq_data)
model_matrix_b_mock <- model.matrix(formula_b, model_frame)

# set up framework to estimate alpha (amplification efficiency)

formula_a <- community ~ Cycles -1

model_frame <- model.frame(formula_a, mock_data_fin)
model_vector_a_mock <- model.matrix(formula_a, model_frame) %>% as.numeric()

N_obs_mock       <- nrow(mock_data_fin)
N_b_mock_col     <- ncol(model_matrix_b_mock)  

## Second, create design matrices for unknown communities
# species compositions (betas)

N_pcr_samp <- sample_seq_data$Cycles

if(length(unique(sample_seq_data$community))==1){
  formula_b <- Cycles ~ 1  
} else {
  formula_b <- Cycles ~ community
}
model_frame <- model.frame(formula_b, sample_data_all)
model_matrix_b_samp <- model.matrix(formula_b, model_frame)

model_frame <- model.frame(formula_b, sample_data_all_small)
model_matrix_b_samp_small <- model.matrix(formula_b, model_frame)

# set up framework to estimate alpha (amplification efficiency)

formula_a <- community ~ Cycles -1

model_frame <- model.frame(formula_a, sample_data_all)
model_vector_a_samp <- model.matrix(formula_a, model_frame) %>% as.numeric()

model_frame <- model.frame(formula_a, sample_data_all_small)
model_vector_a_samp_small <- model.matrix(formula_a, model_frame) %>% as.numeric()

#### Fill in species that are in the unknown samples but not the mock communities ---------------------------------------

# Provide the species that are in the unknown samples but not in the mock communities
sp.missing <- 
  sp.comm.env.nomock %>% 
  distinct(Species) %>%
  as.data.frame() %>% 
  rename("species" = 1) %>% 
  mutate(species = case_when(species == "Oncorhynchus.tshawytscha" ~ "zRefSpecies_Oncorhynchus.tshawytscha", TRUE ~ species)) %>% 
  filter(!(species %in% mock_sp)) %>% 
  pull(species)
 
# invoke one of four possible scenarios.  DO NOT SET MORE THAN ONE TO TRUE (ALL THREE MAY BE FALSE)
MAX.ALPHA = c(TRUE, FALSE, FALSE) # Listed species get the maximum alpha estimated
MIN.ALPHA = c(FALSE, TRUE, FALSE) # Listed species get the minimum alpha estimated
ZERO.ALPHA = c(FALSE, FALSE, TRUE) # Listed species get an alpha equivalent to the reference species.

#list object for stan model output

output.list <- list()

for (a in 1:3) {
# Create a species index with numbers for all species
SP_IDX <- rbind(env, env.nomock) %>% distinct(Species, sp_idx) %>% arrange(sp_idx)

if(MAX.ALPHA[a] == TRUE){
  SP_IDX <- SP_IDX %>% mutate(max.alpha  = ifelse(Species %in% sp.missing,1,0),
                              min.alpha  = 0,
                              zero.alpha = ifelse(!Species %in% sp.missing,1,0))
}else if(MIN.ALPHA[a] == TRUE){
  SP_IDX <- SP_IDX %>% mutate(max.alpha  = 0,
                              min.alpha  = ifelse(Species %in% sp.missing,1,0),
                              zero.alpha = ifelse(!Species %in% sp.missing,1,0))
}else if(ZERO.ALPHA[a] == TRUE){
  SP_IDX <- SP_IDX %>% mutate(max.alpha  = 0,
                              min.alpha  = 0,
                              zero.alpha = ifelse(!Species %in% sp.missing,1,0))
}else{
  SP_IDX <- SP_IDX %>% mutate(max.alpha  = 0,
                              min.alpha  = 0,
                              zero.alpha = 1)
}

# total number of species in both datasets
N_species_all = nrow(SP_IDX)

# set up counters 
N_obs_samp_small <- nrow(model_matrix_b_samp_small)
N_obs_samp       <- nrow(sample_data)
N_b_samp_col     <- ncol(model_matrix_b_samp)  

# create list of data objects to send to stan model
stan_data <- list(
  N_species = N_species,   # Number of species in mock data
  N_species_all = N_species_all, # Number of species in environmental dataset
  N_obs_samp = N_obs_samp, # Number of observed samples 
  N_obs_mock = N_obs_mock, # Number of observed mock samples
  N_obs_samp_small = N_obs_samp_small, # Number of observed samples 
  
  # Observed data of community matrices
  sample_data = sample_data_all[,4:ncol(sample_data_all)],
  mock_data   = mock_data_fin[,4:ncol(mock_data_fin)],
  
  # True proportions for mock community
  alr_mock_true_prop = alr_mock_true_prop,
  
  # vectors of PCR numbers
  N_pcr_samp = as.array(N_pcr_samp),
  N_pcr_mock = N_pcr_mock,
  
  # Design matrices: field samples
  N_b_samp_col = N_b_samp_col,
  model_matrix_b_samp = model_matrix_b_samp,
  model_matrix_b_samp_small = model_matrix_b_samp_small,
  model_vector_a_samp = as.array(model_vector_a_samp),
  model_vector_a_samp_small = as.array(model_vector_a_samp_small),
  
  # Design matrices: mock community samples
  N_b_mock_col = N_b_mock_col,
  model_matrix_b_mock = model_matrix_b_mock,
  model_vector_a_mock = model_vector_a_mock,
  
  # Vectors for assigning min, max or zero for the alpha vector
  max_alpha  = SP_IDX$max.alpha,
  min_alpha  = SP_IDX$min.alpha,
  zero_alpha = SP_IDX$zero.alpha,
  
  # Priors
  alpha_prior = c(0,0.1),  # normal prior
  beta_prior = c(0,5),    # normal prior
  tau_prior = c(1,1)   # gamma prior
)


stan_pars <- c(
  "alpha",
  "beta",
  "mu_samp",
  "mu_mock",
  "eta_mock", # turn off for model without overdispersion
  "tau", # turn off for model without overdispersion
  "int_samp_small"
  )

stan_init_f2 <- function(n.chain,N_species){
  A <- list()
  for(i in 1:n.chain){
    A[[i]] <- list(
      # tau = runif(N_species-1,0.1,0.5),
      alpha_raw = runif(N_species-1,-0.5,0.5)
    )
  }
  return(A)
}


### Run stan Model ---------------------------------------------------------------------------------------

N_CHAIN = 1
Warm = 50
Iter = 50
Treedepth = 12
Adapt_delta = 0.8

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanMod = stan(file = "quant_metabar_no_overdispersion_fixed_alpha.stan", data = stan_data,
               verbose = FALSE, chains = N_CHAIN, thin = 1,
               warmup = Warm, iter = Warm + Iter,
               control = list(adapt_init_buffer = 175,
                              max_treedepth=Treedepth,
                              stepsize=0.01,
                              adapt_delta=Adapt_delta,
                              metric="diag_e"),
               pars = stan_pars,
               refresh = 10,
               boost_lib = NULL,
               init = stan_init_f2(n.chain=N_CHAIN,N_species=N_species),
               sample_file = paste0("./tmpB.csv")
)



### Collect and save model output data ----------------------------------------------------------------
#extract model parameters
pars <- rstan::extract(stanMod, permuted = TRUE)
samp_params <- get_sampler_params(stanMod)

#create list of output objects from the model
stanMod_summary <- list()
stanMod_summary[["alpha"]] <- summary(stanMod,pars="alpha")$summary
stanMod_summary[["beta"]] <- summary(stanMod,pars="beta")$summary
stanMod_summary[["mu_samp"]] <- summary(stanMod,pars="mu_samp")$summary
stanMod_summary[["mu_mock"]] <- summary(stanMod,pars="mu_mock")$summary
stanMod_summary[["tau"]] <- summary(stanMod,pars="tau")$summary
stanMod_summary[["int_samp_small"]] <- summary(stanMod,pars="int_samp_small")$summary

Output <- list(
  
  # input data
  p_true = p_mock,
  p_samp_all = sample_data_all,
  p_mock_all = mock_data_fin,
  env = env, #this should be the same as sample_input_data
  mc = mc, #this should be the same as mock_input_data
  
  # stan input objects
  stan_data = stan_data,
  Warm=Warm,
  Iter=Iter,
  
  # Fitted Objects
  stanMod = stanMod, # Full stan Model fitted object
  pars = pars, # MCMC output
  samp_params=samp_params, # Sampler information
  stanMod_summary = stanMod_summary # posterior summaries.
)

output.list[[a]] <- Output

}

#save(Output,file=paste0("C:/Users/Amy/Google Drive/00 NWFSC/01 Oo fecal DNA/04 Data analysis/prey metabarcoding/quantmeta_QAQC/data/Oorca_no_overdisp",".Rdata"))

### Preliminary Output Plots ------------------------------------------------------------------------------------

# plot model traces for posterior probability and alpha
#rstan::traceplot(stanMod,pars=c("lp__","alpha"),inc_warmup=FALSE) #which alpha is for mock and which is for samples? 
#rstan::traceplot(stanMod,pars=c("lp__","tau"),inc_warmup=FALSE) #which alpha is for mock and which is for samples? 

# extract alphas
alpha.list <- list()

for (a in 1:3) {
alpha.list[[a]] <- data.frame(
  mean_log_alpha = output.list[[a]]$stanMod_summary[["alpha"]][, 1],
  lower_bound_alpha = output.list[[a]]$stanMod_summary[["alpha"]][, 4],
  upper_bound_alpha = output.list[[a]]$stanMod_summary[["alpha"]][, 8],
  sp_idx = 1:length(stanMod_summary[["alpha"]][, 1])
) %>% 
  left_join(bind_rows(sp.list, sp.list.nomock), by = c("sp_idx" = "sp_idx")) %>% 
  mutate(Species = str_replace_all(Species, "\\.", " ")) %>% 
  mutate(Species = str_remove(Species, "zRefSpecies_"))

}

# re-plot alphas with species names
alpha.plot <- list()

for (a in 1:3) {
alpha.plot[[a]] <- ggplot(alpha.list[[a]], aes(x = mean_log_alpha, y = Species)) +
  geom_point() +
  geom_segment(aes(y = Species, yend = Species, x = lower_bound_alpha, xend = upper_bound_alpha)) +
  xlim(-0.5,0.2)

if(a > 1) {
  alpha.plot[[a]] <- alpha.plot[[a]] +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}
}

library(patchwork)

alpha.plot[[1]] + alpha.plot[[2]] + alpha.plot[[3]]

# extract posterior estimations of species proportions for each unknown sample
posterior.list <- list()

for (a in 1:3){
posterior.list[[a]] <- output.list[[a]]$stanMod_summary[["int_samp_small"]][, c(1,4:8)] %>%
  as.data.frame() %>%
  mutate(community = rep(unique(env$community), each = N_species_all)) %>%
  mutate(sp_idx = rep(1:N_species_all, times = length(unique(env$community)))) %>% 
  left_join(SP_IDX[,1:2], by = ("sp_idx" = "sp_idx")) %>% 
  mutate(Species = str_replace_all(Species, "\\.", " ")) %>% 
  mutate(Species = str_remove(Species, "zRefSpecies_"))
}

posterior.plots <- list()

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colorCount = length(unique(posterior.list[[1]]$Species))

for (a in 1:3){
  
posterior.plots[[a]] <- ggplot(posterior.list[[a]], aes(x = community, y = mean, fill = Species)) +
  geom_col()  + 
  scale_fill_manual(values = getPalette(colorCount))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ylab(ifelse(a==1, "max alpha", ifelse(a==2, "min alpha", "zero alpha")))

if (a < 3) {
  posterior.plots[[a]] <- posterior.plots[[a]] +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  
  
} else {
  posterior.plots[[a]] <- posterior.plots[[a]] +
    theme(legend.position = "bottom")
}
}

sample_data_all_prop <- sample_data_all_prop %>% 
  mutate(sp_idx = as.numeric(sp_idx)) %>% 
  left_join(SP_IDX[,1:2], by = ("sp_idx" = "sp_idx")) %>% 
  mutate(Species = str_replace_all(Species, "\\.", " ")) %>% 
  mutate(Species = str_remove(Species, "zRefSpecies_"))

input.plot <- ggplot(sample_data_all_prop, aes(x = community, y = prop.reads, fill = Species)) +
  geom_col()  + 
  scale_fill_manual(values = getPalette(colorCount))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  ylab("input")+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

input.plot / posterior.plots[[1]] / posterior.plots[[2]] / posterior.plots[[3]]
