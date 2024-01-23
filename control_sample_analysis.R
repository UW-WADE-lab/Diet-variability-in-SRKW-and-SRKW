## analysis of Oorc_fecal control samples
## 1/26/2022
## Amy Van Cise

### set up working environment -------------------------------------------------

library(tidyverse)
library(phyloseq)
library(hrbrthemes)
library(ggsci)
library(patchwork)

wd <- "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding"

setwd(wd)

# load dada2 output
load("Oorc_2016-2021_dada2_QAQC_tax95_newref_halibut_output.Rdata")

# load alphas
load("ABC_alpha_vals.Rdata")

# species color palette

# names(species.palette) <- c("Atheresthes stomias", "Hippoglossus stenolepis", "Oncorhynchus keta", 
#                             "Oncorhynchus kisutch", "Oncorhynchus tshawytscha",
#                             "Oncorhynchus mykiss", "Anoplopoma fimbria","Ophiodon elongatus", "Raja binoculata","Parophrys vetulus",
#                             "Clupea pallasii", "Salmo salar",
#                             "Oncorhynchus nerka", "Oncorhynchus gorbuscha")

species.palette = c(pal_uchicago(alpha = 0.6)(9)[c(1,3:9)],pal_jama(alpha = 0.6)(6))

names(species.palette) <- c("arrowtooth", "chum",
                            "coho", "Chinook",
                            "steelhead", "sablefish", "lingcod",  "Pacific halibut","big skate", "English sole",
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

### create master phyloseq object ---------------------------------------------
ps.raw <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   tax_table(taxa))

# shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.raw))
names(dna) <- taxa_names(ps.raw)
ps.all <- merge_phyloseq(ps.raw, dna)
taxa_names(ps.raw) <- paste0("ASV", seq(ntaxa(ps.raw)))

# glom taxa to species at master level, change read counts to proportion, remove Oorc sequences
ps.prop.sp <- tax_glom(ps.raw, "Species") %>% 
  subset_taxa(Genus != "Orcinus") %>% 
  transform_sample_counts(function(x) x / sum(x))

# control sample subset
control.index <- grepl("Control|control", sample_names(ps.prop.sp))

ps.control <- ps.prop.sp %>% 
  prune_samples(control.index, .)

# extract proportions from ps object to data.frame

control.observed <- ps.control@otu_table %>% 
  as.data.frame()
names(control.observed) <- ps.control@tax_table[,7]

control.observed.long <- control.observed %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(cols = -sample, names_to = "species", values_to = "prop.obs")


## add metadata for control samples (both as phyloseq and data.frame) and for raw phyloseq

# read and format control data
control.meta <- read.csv("Oorc_control_metadata.csv", stringsAsFactors = FALSE) %>% 
  filter(!grepl("AmpBias", LabID)) %>% 
  mutate(year = as.factor(year)) %>% 
  column_to_rownames("sample")

# add to phyloseq
sample_data(ps.control) <- control.meta

# reformat metadata
control.meta.long <- control.meta %>% 
  rownames_to_column("sample") %>% 
  pivot_longer(cols = -c(sample, LabID, control.type, year, plate), names_to = "species", values_to = "prop.expected") %>% 
  mutate(species = str_replace(species, "\\.", " ")) %>% 
  replace_na(list(prop.expected = 0))

# add to data frame
control.exp.obs <- control.observed.long %>% 
  left_join(control.meta.long, by = c("sample" = "sample", "species" = "species")) %>% 
  #filter(prop.obs > 0) %>% 
  mutate(prop.obs = prop.obs*100)

# read and format metadata for all samples
samdf <- read.csv("G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding/Oorc_fecal_metadata.csv") %>% 
  mutate(Year = as.factor(Year)) %>% 
  unite("Month.Year", Month:Year, sep = ".", remove = FALSE) %>% 
  column_to_rownames("Sample")

# add to raw data phyloseq
sample_data(ps.raw) <- samdf

### plot observed vs. expected proportion for controls 2 and 3

control23 <- control.exp.obs %>% 
  filter(year %in% c("2018","2019","2021")) %>% 
  filter(control.type %in% c(2,3)) %>% 
  replace_na(list(prop.expected = 0)) %>% 
  left_join(sci.to.common, by = c("species" = "sci_name"))
  
control.type.labs = c("Control 1", "Control 2")
names(control.type.labs) = c("2","3")

control.prop.plot <- ggplot(data = control23, aes(x = prop.expected, y = prop.obs, color = common_name, shape = year)) +
  geom_point(size = 3) +
  scale_color_manual(values = species.palette) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(~control.type, labeller = labeller(control.type = control.type.labs), ncol = 1) +
  theme_bw() +
  theme(panel.border = element_blank()) +
  ylab("Proportion observed reads")+
  xlab("Porportion input DNA")+
  theme(legend.position="bottom", legend.title = element_blank(),
        plot.margin=unit(x=c(0,0,0,1),units="mm"),
        legend.margin=margin(t = 0, r = 0, b = 0, l = 0, unit='pt')) 

legend <- ggpubr::get_legend(control.prop.plot)

control.prop.plot <- control.prop.plot +
  theme(legend.position = "none")

# alpha plot

alpha_plot <- ggplot(alphas, aes(x = mean_log_alpha, y = common_name, color = common_name, fill = common_name)) +
  geom_point(size=1) +
  geom_segment(aes(y = common_name, yend = common_name, x = lower_bound_alpha, xend = upper_bound_alpha))+
  xlab("Mean log(alpha)") +
  ylab("Species") +
  theme_minimal() +
  scale_fill_manual(values = species.palette) +
  scale_color_manual(values = species.palette) +
  theme(axis.text.y = element_text(angle = 45, hjust = 0.5)) + 
  theme(legend.position = "none",
        plot.margin=unit(x=c(0,2,0,0),units="mm"))

controlABCplot <- (alpha_plot + control.prop.plot) / legend +
  plot_layout(heights = c(5,2))

save(controlABCplot, control23, file = "control_sample_alpha_plot.Rdata")
ggsave("Oorc_control_sample_alpha_plot.png", controlABCplot, height = 4, width = 8, units = "in")

### compare species proportion by plate to difference from expected proportion

# full phyloseq merged by plate
ps.raw.merged <- tax_glom(ps.raw, "Species") %>% 
  merge_samples("MiSeq.run") %>% 
  transform_sample_counts(function(x) x / sum(x))

# export data frame
allsamps.prop <- ps.raw.merged@otu_table %>% 
  as.data.frame() %>% 
  dplyr::slice(2:n())
names(allsamps.prop) <- ps.raw.merged@tax_table[,7] 

allsamps.prop.long <- allsamps.prop %>% 
  rownames_to_column("plate") %>% 
  filter(plate != "") %>% 
  pivot_longer(cols = -plate, names_to = "species", values_to = "prop.obs.plate") %>% 
  mutate(prop.obs.plate = prop.obs.plate * 100)

# calcuate obs-exp for control samples
control.diff <- control.exp.obs %>% 
  replace_na(list(prop.expected = 0)) %>% 
  mutate(diff.exp.obs = prop.obs - prop.expected) %>% 
  left_join(allsamps.prop.long, by = c("plate" = "plate", "species" = "species")) %>% 
  filter(!is.na(prop.obs.plate)) 

# plot expected-observed vs. proportion on plate

control.diff23 <- control.diff %>% 
  #filter(year %in% c("2019","2021")) %>% 
  filter(control.type %in% c(2,3)) #%>% 
  #filter(abs(diff.exp.obs) > 1)

control_names <- list("2" = "Control 1", "3" = "Control 2")
control_labeller <- function(variable,value){
  return(control_names[value])
}
control.diff.plot <- ggplot(control.diff23, aes(x = prop.obs.plate, y = diff.exp.obs)) +
  geom_point(size = 3, aes(color = species, shape = plate)) + 
  geom_smooth(method = "lm", color = "grey50", fill = NA) +
  scale_color_futurama() +
  facet_wrap(~control.type,labeller=control_labeller) +
  theme_ipsum() +
  geom_hline(aes(yintercept=0)) +
  ylim(-10,16) +
  xlab("Proportion on plate") +
  ylab("Difference expected - observed")
  
save(control.prop.plot, control.diff.plot, file = "control_plots.Rdata")

