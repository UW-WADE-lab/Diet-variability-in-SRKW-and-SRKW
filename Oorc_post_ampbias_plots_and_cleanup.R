library(tidyverse)
library(ggsci)
library(patchwork)
library(RColorBrewer)
library(ggpubr)
library(hrbrthemes)
library(phyloseq)
library(microbiomeutilities)

### load data ------------------------------------------------------------------
setwd("G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding")

load("QAQC_allsamps_overdisp_tax95_newref_halibut_08312023.Rdata")

list2env(Output, .GlobalEnv)
N_species <- nrow(sp.list)

### Preliminary Output Plots ---------------------------------------------------

# plot model traces for posterior probability and alpha
#rstan::traceplot(stanMod,pars=c("lp__","alpha"),inc_warmup=FALSE)  
#rstan::traceplot(stanMod,pars=c("lp__","tau"),inc_warmup=FALSE) 

# extract alphas
alphas <- data.frame(
  mean_log_alpha = stanMod_summary[["alpha"]][, 1],
  lower_bound_alpha = stanMod_summary[["alpha"]][, 4],
  upper_bound_alpha = stanMod_summary[["alpha"]][, 8],
  sp_idx = 1:length(stanMod_summary[["alpha"]][, 1])
) %>% 
  left_join(sp.list, by = c("sp_idx" = "sp_idx")) %>% 
  mutate(Species = str_replace_all(Species, "\\.", " ")) %>% 
  mutate(Species = str_remove(Species, "zRefSpecies_")) %>% 
  mutate(common_name = c("Pacific herring", "Pacific halibut", "pink", "chum", "coho", "steelhead", "sockeye", "lingcod", "Atlantic salmon", "Chinook"))

# save data for control plot figure
save(alphas, file = "ABC_alpha_vals.Rdata")

### Posterior estimations of species proportions for each unknown sample -------
posteriorProportions <- stanMod_summary[["int_samp_small"]][, c(1,4:8)] %>%
  as.data.frame() %>%
  mutate(community = rep(unique(env$community), each = N_species)) %>%
  mutate(sp_idx = rep(1:N_species, times = length(unique(env$community))))

# add metadata
sample_metadata <- read.csv("KW_fecal_raw_readcount_allsamps_allspecies_tax95_newref_halibut_08312023.csv") %>% 
  select(Sample, Population, Pod) %>% 
  separate(Sample, sep = "_", into = c("community",NA))# %>% 
  #mutate(community = str_remove(community, paste0(c("rep", "b", "B", "-dup"), collapse = "|")))

SP_IDX <- SP_IDX %>% 
  mutate(Species = str_remove(Species,"zRefSpecies_"))

### Make a dataframe of common names for fish species --------------------------

sci.to.common <- data.frame(sci_name = c("Oncorhynchus.mykiss","Oncorhynchus.tshawytscha", "Ophiodon.elongatus", "Orcinus.orca","Oncorhynchus.nerka",      
                                         "Oncorhynchus.keta", "Parophrys.vetulus", "Raja.binoculata", "Hippoglossus.stenolepis",  "Oncorhynchus.kisutch",  
                                         "Atheresthes.stomias", "Clupea.pallasii", "Oncorhynchus.gorbuscha", "Squalus.acanthias",  "Anoplopoma.fimbria",      
                                         "Zaprora.silenus", "Eumicrotremus.orbis", "Thaleichthys.pacificus", "Microstomus.pacificus", "Liparis.sp..BOLD.AAB4898",
                                         "Citharichthys.sordidus"),
                            common_name = c("steelhead", "Chinook", "lingcod", "killer whale", "sockeye",
                                            "chum", "English sole", "big skate", "Pacific halibut", "coho",
                                            "arrowtooth flounder", "Pacific herring", "pink salmon", "spiny dogfish", "sablefish",
                                            "prowfish", "Pacific spiny lumpsucker", "eulachon", "Pacific sole", "snailfish (genus)",
                                            "Pacific sanddab"))

### Convert posterior proportions to reads for each sample ----------------------
  
sample_posterior_reads <- sample_seq_data_input %>% 
  mutate(Species = str_remove(Species,"zRefSpecies_")) %>% 
  left_join(SP_IDX, by = c("Species" = "Species")) %>% 
  left_join(posteriorProportions, by = c("sp_idx" = "sp_idx", "community" = "community")) %>% 
  group_by(community) %>% 
  mutate(total_community_reads = sum(nReads[!(is.na(sp_idx))])) %>% 
  mutate(post.reads = mean*total_community_reads) %>% 
  select(community, Species, nReads, mean, post.reads) %>% 
  ungroup() %>% 
  mutate(nReads = as.numeric(nReads)) %>% 
  mutate(post.prop.reads = case_when(is.na(post.reads) ~ nReads, TRUE ~ post.reads)) %>% 
  left_join(sample_metadata, by = c("community" = "community")) %>% 
  distinct(.keep_all = TRUE) %>% 
  left_join(sci.to.common, by = c("Species" = "sci_name"))

sample_posterior_prop <- sample_posterior_reads %>% 
  group_by(community) %>% 
  mutate(total_community_reads = sum(post.prop.reads)) %>% 
  mutate(post.prop.readprop = post.prop.reads/total_community_reads) %>% 
  filter(post.prop.readprop > 0.01) 

sample_posterior_reads_ampbias <- sample_posterior_reads %>% 
  filter(!is.na(mean))

ARKWsample_posterior_reads_ampbias <- sample_posterior_reads_ampbias %>% filter(Population=="ARKW")
SRKWsample_posterior_reads_ampbias <- sample_posterior_reads_ampbias %>% filter(Population=="SRKW")

ARKWsample_posterior_prop <- sample_posterior_prop %>% filter(Population=="ARKW")
SRKWsample_posterior_prop <- sample_posterior_prop %>% filter(Population=="SRKW")

### Plot proportions before and after correction -------------------------------

# set color palette and legend labels for ARKW and SRKW
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
ARKWcolorCount = length(unique(c(ARKWsample_posterior_prop$common_name, alphas$common_name)))
SRKWcolorCount = length(unique(c(SRKWsample_posterior_prop$common_name, alphas$common_name)))

ARKWspecies <- unique(c(ARKWsample_posterior_prop$common_name, alphas$common_name))
SRKWspecies <- unique(c(SRKWsample_posterior_prop$common_name, alphas$common_name))

ARKWcolor_map <- set_names(getPalette(ARKWcolorCount), ARKWspecies)
SRKWcolor_map <- set_names(getPalette(SRKWcolorCount), SRKWspecies)

# plot alphas with species names
ARKWalpha_plot <- ggplot(alphas, aes(x = mean_log_alpha, y = common_name, color = common_name, fill = common_name)) +
  geom_point(size=1) +
  geom_segment(aes(y = common_name, yend = common_name, x = lower_bound_alpha, xend = upper_bound_alpha))+
  xlab("Mean log(alpha)") +
  ylab("Species") +
  theme_minimal() +
  scale_fill_manual(values = ARKWcolor_map) +
  scale_color_manual(values = ARKWcolor_map) +
  theme(axis.text.y = element_text(angle = 45, hjust = 0.5)) + 
  theme(legend.position = "none")

SRKWalpha_plot <- ggplot(alphas, aes(x = mean_log_alpha, y = common_name, color = common_name, fill = common_name)) +
  geom_point(size=1) +
  geom_segment(aes(y = common_name, yend = common_name, x = lower_bound_alpha, xend = upper_bound_alpha))+
  xlab("Mean log(alpha)") +
  ylab("Species") +
  theme_minimal() +
  scale_fill_manual(values = SRKWcolor_map) +
  scale_color_manual(values = SRKWcolor_map) +
  theme(axis.text.y = element_text(angle = 45, hjust = 0.5)) + 
  theme(legend.position = "none")

# calculate proportions before correction
sample_seq_prop <- sample_seq_data_input %>% 
  mutate(Species = str_remove(Species,"zRefSpecies_")) %>% 
  filter(Species %in% SP_IDX$Species) %>% 
  group_by(community) %>% 
  mutate(total_community_reads = sum(nReads)) %>% 
  mutate(prop.reads = nReads/total_community_reads) %>% 
  left_join(sample_metadata, by = c("community" = "community")) %>% 
  distinct(.keep_all = TRUE) %>% 
  left_join(sci.to.common, by = c("Species" = "sci_name"))

ARKWsample_seq_prop <- sample_seq_prop %>% filter(Population=="ARKW") 
SRKWsample_seq_prop <- sample_seq_prop %>% filter(Population=="SRKW")

# ARKW plots
ARKWafter <- ggplot(ARKWsample_posterior_reads_ampbias, aes(x=community, y=mean, fill = common_name)) +
  geom_col() +
  scale_fill_manual(values = ARKWcolor_map)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels = NULL) +
  theme_classic() +
  ylab("Proportion of reads") +
  labs(fill = "Species common name") +
  theme(legend.position = "none")

ARKWbefore <- ggplot(ARKWsample_seq_prop, aes(x=community, y=prop.reads, fill = common_name)) +
  geom_col() +
  scale_fill_manual(values = ARKWcolor_map)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels = NULL) +
  theme_classic() +
  ylab("") +
  labs(fill = "Species common name") +
  theme(legend.position = "none")
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ARKWafter_all <- ggplot(ARKWsample_posterior_prop, aes(x=community, y=post.prop.readprop, fill = common_name)) +
  geom_col() +
  scale_fill_manual(values = ARKWcolor_map) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels = NULL) +
  theme_classic() +
  ylab("") +
  labs(fill = "Species common name")

pdf(file = "ARKW_ampbias_correction.pdf", width = 10, height = 6)
(ARKWalpha_plot + (ARKWbefore / ARKWafter / ARKWafter_all)) + 
  plot_layout(guides = "collect") +
  plot_layout(widths = c(1,3))

dev.off()

# SRKW plots
SRKWafter <- ggplot(SRKWsample_posterior_reads_ampbias, aes(x=community, y=mean, fill = common_name)) +
  geom_col() +
  scale_fill_manual(values = SRKWcolor_map)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels = NULL) +
  theme_classic() +
  ylab("Proportion of reads") +
  theme(legend.position = "none")

SRKWbefore <- ggplot(SRKWsample_seq_prop, aes(x=community, y=prop.reads, fill = common_name)) +
  geom_col() +
  scale_fill_manual(values = SRKWcolor_map)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels = NULL) +
  theme_classic() +
  ylab("Proportion of reads") +
  theme(legend.position = "none")

SRKWafter_all <- ggplot(SRKWsample_posterior_prop, aes(x=community, y=post.prop.readprop, fill = common_name)) +
  geom_col() +
  scale_fill_manual(values = SRKWcolor_map) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels = NULL) +
  theme_classic() +
  ylab("Proportion of reads") +
  labs(fill = "Species common name")

pdf(file = "SRKW_ampbias_correction.pdf", width = 10, height = 6)
(SRKWalpha_plot + (SRKWbefore / SRKWafter / SRKWafter_all)) + 
  plot_layout(guides = "collect") +
  plot_layout(widths = c(1,3))

dev.off()

### Convert back to phyloseq object --------------------------------------------

# load dada2 output
load("Oorc_2016-2021_dada2_QAQC_tax95_newref_halibut_output.Rdata")
load("ASVreference_allsamps_tax95_newref_halibut_08312023.Rdata")

#new otu table
sample_otucount <- read.csv("KW_fecal_raw_readcount_allsamps_allspecies_tax95_newref_halibut_08312023.csv") %>% 
  left_join(read.csv("KW_fecal_raw_readcount_allsamps_asvs_tax95_newref_halibut_08312023.csv"), 
            by = c("Sample","LabID","FieldID","bio_rep","Population","year","month","Pod")) %>% 
  separate(Sample, sep = "_", into = c("community",NA), remove = FALSE) %>% 
  select(-LabID, -FieldID, -bio_rep,-Population, -year,-month,-Pod) %>% 
  pivot_longer(!Sample:community, names_to = "Species", values_to = "nReads") %>% 
  left_join(sample_posterior_reads, by = c("community" = "community", "Species"="Species")) %>% 
  rename("nReads"=nReads.x) %>% 
  select(Sample, Species, nReads) %>% 
  mutate(Species = str_replace(Species, "\\.", " ")) %>% 
  mutate(Species = case_when(grepl("Liparis sp",Species) ~ "Liparis sp", grepl("Osmeridae sp",Species) ~ "Osmeridae sp", TRUE ~ Species)) %>%
  pivot_wider(id_cols = Sample, names_from = Species, values_from = nReads) %>% 
  column_to_rownames("Sample") %>% 
  as.matrix()

#new tax table
tax_convert <- as.data.frame(taxa) %>% 
  rownames_to_column(var="DNAseq") %>% 
  left_join(reference.table, by = "DNAseq") %>% 
  mutate(Species = case_when(grepl("Liparis sp",Species) ~ "Liparis sp", Species == "Osmeridae sp. BOLD:AAC0320" ~ "Osmeridae sp", TRUE ~ Species)) %>% 
  mutate(Species.asv = case_when(is.na(Species) ~ ASVname, TRUE ~ Species)) %>% 
  distinct(Species.asv, .keep_all = TRUE) %>% 
  remove_rownames() %>% 
  column_to_rownames("Species.asv")

new_tax_table <- tax_convert %>% 
  select(-DNAseq,-ASVname) %>% 
  as.matrix()

new_reference_DNA <- tax_convert %>% 
  pull(DNAseq)

names(new_reference_DNA) <- rownames(tax_convert)

# sample metadata
samdf <- read.csv("new_prey_meta_7.26.22_CE_KP.csv") %>% 
  mutate(year = as.factor(year)) %>% 
  mutate(month = as.factor(month)) %>% 
  unite("Month.Year", month:year, sep = ".", remove = FALSE) %>% 
  unite("Pop.year", c(Population,year), sep = ".", remove = FALSE) %>% 
  unite("Pop.month", c(Population,month), sep = ".", remove = FALSE) %>% 
  unite("Pop.month.year", c(Population,month,year), sep = ".", remove = FALSE) %>% 
  unite("Region.month", c(Region,month), sep = ".", remove = FALSE) %>% 
  filter(Population != "NRKW") %>% 
  column_to_rownames("Sample")

### create master phyloseq object
ps.ampbias <- phyloseq(otu_table(sample_otucount, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(new_tax_table))

save(ps.ampbias, sample_posterior_prop, new_reference_DNA, file = "psconvert_QAQC_allsamps_overdisp_tax95_newref_halibut_08312023.Rdata")

     