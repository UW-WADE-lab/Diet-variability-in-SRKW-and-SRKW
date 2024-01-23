## analysis of Oorc_fecal technical replicates
## 1/26/2022
## Amy Van Cise

### set up working environment --------------------------------------------------

library(tidyverse)
library(phyloseq)
library(hrbrthemes)
library(ggsci)

wd <- "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding"

setwd(wd)

# load dada2 output
load("Oorc_2016-2021_dada2_QAQC_tax95_newref_halibut_output.Rdata")

# species palatte
species.palette = c(pal_uchicago(alpha = 0.6)(9),pal_jama(alpha = 0.6)(5))
# names(species.palette) <- c("Atheresthes stomias", "Hippoglossus stenolepis", "Oncorhynchus keta", 
#                             "Oncorhynchus kisutch", "Oncorhynchus tshawytscha",
#                             "Oncorhynchus mykiss", "Anoplopoma fimbria","Ophiodon elongatus", "Raja binoculata","Parophrys vetulus",
#                             "Clupea pallasii", "Salmo salar",
#                             "Oncorhynchus nerka", "Oncorhynchus gorbuscha")
names(species.palette) <- c("arrowtooth flounder", "Pacific halibut", "chum",
                            "coho", "Chinook",
                            "steelhead", "sablefish", "lingcod", "big skate", "English sole",
                            "Pacific herring", "Atlantic salmon", 
                            "sockeye", "pink")

# species common names
sci.to.common <- data.frame(sci_name = c("Oncorhynchus mykiss","Oncorhynchus tshawytscha", "Ophiodon elongatus", "Orcinus orca","Oncorhynchus nerka",      
                                         "Oncorhynchus keta", "Parophrys vetulus", "Raja binoculata", "Hippoglossus stenolepis",  "Oncorhynchus kisutch",  
                                         "Atheresthes stomias", "Clupea pallasii", "Oncorhynchus gorbuscha", "Squalus acanthias",  "Anoplopoma fimbria",      
                                         "Zaprora silenus", "Eumicrotremus orbis", "Thaleichthys pacificus", "Microstomus pacificus", "Liparis.sp..BOLD.AAB4898",
                                         "Citharichthys sordidus", "Salmo salar"),
                            common_name = c("steelhead", "Chinook", "lingcod", "killer whale", "sockeye",
                                            "chum", "English sole", "big skate", "Pacific halibut", "coho",
                                            "arrowtooth flounder", "Pacific herring", "pink", "spiny dogfish", "sablefish",
                                            "prowfish", "Pacific spiny lumpsucker", "eulachon", "Pacific sole", "snailfish (genus)",
                                            "Pacific sanddab", "Atlantic salmon"))


# get sample metadata
samdf <- read.csv("G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding/preybarcode_metadata.csv") %>% 
  mutate(year = as.factor(year)) %>% 
  unite("Month.Year", month:year, sep = ".", remove = FALSE) %>% 
  unite("Pop.year", c(Population,year), sep = ".", remove = FALSE) %>% 
  column_to_rownames("Sample")

### create master phyloseq object --------------------------------------------------
ps.raw <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf),
                   tax_table(taxa))

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.raw))
names(dna) <- taxa_names(ps.raw)
ps.all <- merge_phyloseq(ps.raw, dna)
taxa_names(ps.raw) <- paste0("ASV", seq(ntaxa(ps.raw)))

### glom taxa to species at master level, change read counts to proportion, remove Oorc sequences
ps.prop.sp <- tax_glom(ps.raw, "Species") %>% 
  subset_taxa(Genus != "Orcinus") %>% 
  transform_sample_counts(function(x) x / sum(x)) 

# filter to remove low count taxa 
# must be at least 1% of diet in 1 or more samples
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps.prop.sp, f1, A=1)

ps.prop.sp <- prune_taxa(lowcount.filt, ps.prop.sp)

### technical replicate sample subset
techrep.index <- grepl("rep|dup", sample_names(ps.prop.sp))

ps.reps <- ps.prop.sp %>% 
  prune_samples(techrep.index, .)

rep.samps <- strsplit(sample_names(ps.reps), "_") %>% 
  sapply("[", 1) %>% 
  str_remove("rep") %>% 
  str_remove("-dup") %>% 
  c("51530-008", "51530-009", "52357-05")

ps.techreps <- ps.prop.sp %>% 
  subset_samples(LabID %in% rep.samps) 

# plot proportions
techrep.prop.plot <- 
  plot_bar(ps.techreps, fill = "Species") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  scale_fill_manual(values = species.palette, labels = c("sablefish", "arrowtooth flounder",
                                                         "Pacific herring", "Pacific halibut",
                                                         "pink", "chum",
                                                         "coho", "steelhead", "sockeye",
                                                         "Chinook", "lingcod", "Atlantic salmon")) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.2)) +
  ylab("Prop. obs. reads") +
  xlab("Expected read proportion") +
  facet_wrap(LabID~., scales = "free_x", ncol = 15) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  theme(strip.text = element_blank(),
        panel.margin.y = unit(0, "lines"))+
  #theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(0.2,1,0.2), expand = c(0,0)) 

#plot raw reads
common.species <- ps.techreps@tax_table@.Data[,7]

ps.raw.techrep <- tax_glom(ps.raw, "Species") %>% 
  subset_taxa(Genus != "Orcinus") %>% 
  subset_taxa(Species %in% common.species) %>% 
  subset_samples(LabID %in% rep.samps)

techrep.count.plot <- plot_bar(ps.raw.techrep, fill = "Species") +
  theme_ipsum() +
  scale_fill_futurama() +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.2)) +
  ylab("Number of reads")

save(techrep.count.plot,techrep.prop.plot, file = "techrep_plots.Rdata")

ggsave("Oorc_techrep_prop_plot.png", techrep.prop.plot, height = 4, width = 8, units = "in")