### High-level summary of data
## Amy Van Cise
## March 2022

library(tidyverse)
library(ggsci)
library(phyloseq)
library(phyloseqCompanion)
library(microViz)
library(vegan)

### Set up environment ---------------------------------------------------------------------------
wd <- "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding"
setwd(wd)

load("KWfecal_2016-2019_allsamps_tax95_newref_halibut_ampbias_overdisp_QAQC_08312023.Rdata")

#filter out AOKW population
ps.sp <- ps.sp %>% ps_filter(Population != "AOKW")
nsamples(ps.sp)
ps.spcount <- ps.spcount %>% ps_filter(Population != "AOKW")
nsamples(ps.sp)

### Minor taxa -----------------------------------------------------------------

ps_sanddab <- ps.sp %>% 
  prune_taxa("Citharichthys sordidus",.) %>% 
  prune_samples(sample_sums(.) > 0.01,.) 

dabs <- c("Pacific sanddab",nsamples(ps_sanddab), round(mean(sample_sums(ps_sanddab)*100), digits = 1), 
          ps_sanddab@sam_data$Population, as.character(ps_sanddab@sam_data$month))

ps_sockeye <- ps.sp %>% 
  prune_taxa("Oncorhynchus nerka",.) %>% 
  prune_samples(sample_sums(.) > 0.01,.) 

sockeye <- c("sockeye salmon", nsamples(ps_sockeye), round(mean(sample_sums(ps_sockeye)*100),digits = 1), 
             unique(ps_sockeye@sam_data$Population), paste0(as.character(ps_sockeye@sam_data$month),collapse = ","))

ps_skate <- ps.sp %>% 
  prune_taxa("Raja binoculata",.) %>% 
  prune_samples(sample_sums(.) > 0.01,.) 

skate <- c("big skate", nsamples(ps_skate), round(mean(sample_sums(ps_skate)*100),digits = 1), 
             unique(ps_skate@sam_data$Population), paste0(unique(as.character(ps_skate@sam_data$month)),collapse = ","))

ps_prowfish <- ps.sp %>% 
  prune_taxa("Zaprora silenus",.) %>% 
  prune_samples(sample_sums(.) > 0.01,.) 

prowfish <- c("prowfish", nsamples(ps_prowfish), round(mean(sample_sums(ps_prowfish)*100),digits = 1), 
              unique(ps_prowfish@sam_data$Population), paste0(unique(as.character(ps_prowfish@sam_data$month)),collapse = ","))

minor.prey <- rbind(dabs, sockeye, skate, prowfish) %>% 
  as.data.frame() %>% 
  rename("Species" = 1,
         "No. samples" = 2,
         "Mean proportion (%)" = 3,
         "Population" = 4,
         "Month(s)" = 5)

### Major taxa -----------------------------------------------------------------
#repeat filter to remove low count taxa in diet: 
#must be at least 1% of diet in 4 or more samples in both populations
#we're repeating this because the last subset allowed 1% of diet in 1 more more samples
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps.sp, f1, A=4)
ps.sp <- prune_taxa(lowcount.filt, ps.sp)
nsamples(ps.sp)

# also remove low count taxa from count dataset
ps.spcount <- prune_taxa(lowcount.filt, ps.spcount) 
nsamples(ps.spcount)

# PERMANOVA of major taxa
allsamps.otu <- as.data.frame(ps.sp@otu_table) 

allsamps.env <- sample.data.frame(ps.sp)

pop_difference <- adonis2(allsamps.otu ~ Population + month + year, data = allsamps.env)
pop_difference

# convert ps.sp to longform dataframe
allsamps.long <- as.data.frame(ps.sp@otu_table) %>% 
  rownames_to_column(var = "SampleID") %>% 
  left_join(sample.data.frame(ps.sp) %>% rownames_to_column(var = "SampleID"), 
            by = c("SampleID" = "SampleID")) %>% 
  pivot_longer(cols = 2:9, names_to = "species", values_to = "proportion") %>% 
  mutate(proportion = round(proportion*100, digits = 2)) %>% 
  filter(LabID != "")

# remove sample/species combos if proportion is less than 1%
allsamps.nozero <- allsamps.long %>% filter(proportion >= 1)
  
### plot count of samples by species in each population
prey.present.count <- ggplot(data = allsamps.nozero, aes(x = species, fill = Population)) +
  geom_histogram(stat = "count", position = position_dodge2(width = 0.9,preserve = "single"), alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(limits = c("Oncorhynchus tshawytscha", "Oncorhynchus keta", 
                              "Oncorhynchus kisutch", "Hippoglossus stenolepis",
                              "Atheresthes stomias", "Anoplopoma fimbria", 
                              "Oncorhynchus mykiss", "Ophiodon elongatus"),
                   labels = c("Chinook", "chum", "coho", "halibut", "arrowtooth", "sablefish", "steelhead", "lingcod")) +
  theme_bw() +
  scale_fill_manual(values = c("cornflowerblue","purple"))+
  theme(panel.border = element_blank()) +
  ylab("Number of samples with species present")+
  xlab("")

### make table of proportion of samples with each species for each population

prop.species.pop <- allsamps.nozero %>% 
  group_by(Population,species) %>% 
  mutate(Count.samples = n()) %>% 
  ungroup() %>% 
  group_by(Population) %>% 
  mutate(Pop.samples = length(unique(LabID))) %>% 
  ungroup() %>% 
  mutate(Prop.sample.species = Count.samples/Pop.samples*100) %>% 
  group_by(Population,species) %>% 
  slice_head() %>% 
  ungroup()

### plot proportion of samples with each species present in each population

prey.present.prop <- ggplot(data = prop.species.pop, aes(x = species, y = Prop.sample.species, fill = Population)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9,preserve = "single"), alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylim(0,90)+
  scale_x_discrete(limits = c("Oncorhynchus tshawytscha", "Oncorhynchus keta", 
                              "Oncorhynchus kisutch",
                              "Atheresthes stomias", "Anoplopoma fimbria", 
                              "Oncorhynchus mykiss", "Ophiodon elongatus"),
                   labels = c("Chinook", "chum", "coho", "arrowtooth", "sablefish", "steelhead", "lingcod")) +
  theme_bw() +
  scale_fill_manual(values = c("cornflowerblue","purple"))+
  theme(panel.border = element_blank()) +
  ylab("Species present in sample (%)")+
  xlab("")+
  theme(legend.position = "none")

### make table of proportion of samples with each species > 50% per population 

prey.50.pop <- allsamps.nozero %>% 
  filter(proportion > 50) %>% 
  group_by(Population,species) %>% 
  mutate(Count.samples = n()) %>% 
  ungroup() %>% 
  group_by(Population) %>% 
  mutate(Pop.samples = length(unique(LabID))) %>% 
  ungroup() %>% 
  mutate(Prop.sample.species = Count.samples/Pop.samples*100) %>% 
  group_by(Population,species) %>% 
  slice_head() %>% 
  ungroup()

### plot prop of samples by 50% species in each population

prey.50.plot <- ggplot(data = prey.50.pop, aes(x = species, y = Prop.sample.species)) +
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9,preserve = "single"), aes(fill = Population), alpha = 0.6) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylim(0,90)+
  scale_x_discrete(limits = c("Oncorhynchus tshawytscha", "Oncorhynchus keta", 
                              "Oncorhynchus kisutch",
                              "Atheresthes stomias", "Anoplopoma fimbria", 
                              "Oncorhynchus mykiss", "Ophiodon elongatus"),
                   labels = c("Chinook", "chum", "coho", "arrowtooth", "sablefish", "steelhead", "lingcod")) +
  theme_bw() +
  scale_fill_manual(values = c("cornflowerblue","purple"))+
  theme(panel.border = element_blank()) +
  ylab("Species > 50% sample (%)")+
  xlab("")+
  theme(legend.position = "top")

presence.legend <- ggpubr::get_legend(prey.50.plot)

prey.50.plot <- prey.50.plot +
  theme(legend.position = "none")

pdf(file = "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/05 Manuscript/figure_pdf/species_presence.pdf",
    width = 10, height = 8)
plot_grid(presence.legend,
          prey.present.prop,prey.50.plot, 
          ncol=1, rel_heights = c(1,4,4,4))
dev.off()

prey_presence_plot <- plot_grid(presence.legend,
                                prey.present.prop,prey.50.plot, 
                                ncol=1, rel_heights = c(1,4,4,4))

save(ps.sp, ps.spcount, pop_difference, prey_presence_plot, minor.prey, file = "Oorc_95tax_newref_halibut_phyloseq_majortaxfilt.Rdata")
save.image(file = "Oorc_95tax_newref_halibut_phyloseq_majortaxfilt_prey_hist_minor_species.Rdata")

ggsave("Oorc_prey_presence_plot.png", prey_presence_plot, height = 4, width = 8, units = "in")