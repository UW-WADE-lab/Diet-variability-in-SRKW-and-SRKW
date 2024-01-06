### QAQC Master phyloseq object after amplification bias correction
## Amy Van Cise
## August 2023

library(tidyverse)
library(ggsci)
library(phyloseq)
library(microViz)
library(phyloseqCompanion)
library(patchwork)

#### set up environment ------------------------------------------------------- 

wd <- "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding"

setwd(wd)

load("psconvert_QAQC_allsamps_overdisp_tax95_newref_halibut_08312023.Rdata")

mycolors <- pal_futurama()(12)[c(1:4,8,9)]

# create prey species factor object
species.palette = c(pal_uchicago(alpha = 0.5)(9),pal_jama(alpha = 0.8)(6))
names(species.palette) <- c("Atheresthes stomias", "Hippoglossus stenolepis", "Oncorhynchus keta", 
                            "Oncorhynchus kisutch", "Oncorhynchus tshawytscha",
                            "Oncorhynchus mykiss", "Anoplopoma fimbria","Ophiodon elongatus", "Raja binoculata","Parophrys vetulus",
                            "Clupea pallasii", "Salmo salar",
                            "Oncorhynchus nerka", "Oncorhynchus gorbuscha", NA)

allsamps <- as.data.frame(as.matrix(ps.ampbias@sam_data)) %>% 
  rownames_to_column("SampleID") %>% 
  pull(SampleID)

#### dedup LabID -------------------------------------------------------------
# look at duplicate labID
LabID_dup <- as.data.frame(as.matrix(ps.ampbias@sam_data)) %>% 
  rownames_to_column("SampleID") %>% 
  group_by(LabID) %>% 
  mutate(n.obs = n()) %>% 
  ungroup() %>% 
  filter(n.obs > 1) %>% 
  pull(SampleID)

ps.labdups <- prune_samples(ps.ampbias, samples = LabID_dup) 

labdup_sample_data <- ps.labdups@sam_data %>% 
  rownames_to_column(var = "sample_name") %>% 
  as.data.frame() %>% 
  group_by(LabID) %>% 
  mutate(Lab_label = cur_group_id()) %>% 
  column_to_rownames("sample_name") 

sample_data(ps.labdups) <- labdup_sample_data

ps.labdups.prop <- ps.labdups %>% #tax_glom(ps.inddups, "Species", NArm = FALSE) %>% 
  subset_taxa(Genus != "Orcinus") %>% 
  transform_sample_counts(function(x) x / sum(x))

# plot dup individual/date samples
labdup.prop.plot <- 
  plot_bar(ps.labdups.prop, fill = "Species") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  scale_fill_manual(values = species.palette, labels = c("sablefish", "arrowtooth", "herring", "halibut", "pink", "chum", "coho", "steelhead", "sockeye","Chinook", 
                                                         "lingcod", "big skate", "Atlantic salmon", "NA")) +
  ylab("Prop. reads by species") +
  xlab("Laboratory replicates") +
  facet_wrap(Lab_label~., scales = "free_x", ncol = 10) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  theme(strip.text = element_blank(),
        panel.margin.y = unit(0, "lines"))+
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(0.2,1,0.2), expand = c(0,0))

# remove duplicate labID
LabID_keep <- as.data.frame(as.matrix(ps.ampbias@sam_data)) %>% 
  rownames_to_column("SampleID") %>% 
  mutate(readcount = phyloseq::sample_sums(ps.ampbias)) %>% 
  group_by(LabID) %>% 
  slice_max(readcount) %>% 
  ungroup() %>% 
  pull(SampleID)

ps.raw <- prune_samples(ps.ampbias, samples = LabID_keep)
nsamples(ps.raw)

# remove samples with too few reads
#ps.raw <- prune_samples(sample_sums(ps.raw) >= 25000,ps.raw)

### dedup by Individual/Date -------------------------------------------------
# look at dup individual/date samples for consistency
Ind.date_dups <- as.data.frame(as.matrix(ps.raw@sam_data)) %>% 
  rownames_to_column("SampleID") %>% 
  mutate(readcount = phyloseq::sample_sums(ps.raw)) %>% 
  group_by(Individual.ID, year, month, day) %>% 
  mutate(numsamps = n()) %>% 
  filter(numsamps > 1) %>% 
  ungroup() %>% 
  pull(SampleID)

ps.inddups <- prune_samples(ps.raw, samples = Ind.date_dups)
nsamples(ps.inddups)

ps.indups.prop <- ps.inddups %>% #tax_glom(ps.inddups, "Species", NArm = FALSE) %>% 
  subset_taxa(Genus != "Orcinus") %>% 
  transform_sample_counts(function(x) x / sum(x)) 

sample_data(ps.indups.prop)$Inddate <- paste0(sample_data(ps.indups.prop)$Individual.ID,
                                              sample_data(ps.indups.prop)$year,
                                              sample_data(ps.indups.prop)$month,
                                              sample_data(ps.indups.prop)$day)

# plot dup individual/date samples
indrep.prop.plot <- 
  plot_bar(ps.indups.prop, fill = "Species") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  scale_fill_manual(values = species.palette, labels = c("sablefish", "arrowtooth", "herring", "halibut", "pink", "chum", "coho", "steelhead", "sockeye","Chinook", 
                                                         "lingcod", "big skate", "Atlantic salmon", "NA")) +
  ylab("Prop. reads by species") +
  xlab("Sample duplicates") +
  facet_wrap(Inddate~., scales = "free_x", ncol = 15) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  theme(strip.text = element_blank(),
        panel.margin.y = unit(0, "lines"))+
  theme(legend.position = "bottom") +
  scale_y_continuous(breaks = seq(0.2,1,0.2), expand = c(0,0)) 

 pdf(file = "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/05 Manuscript/figure_pdf/ind_rep_comparison.pdf",
     width = 10, height = 6)
 indrep.prop.plot
 dev.off()
 
 ## Combine replicate plots for Oorc_diet paper --------------------------------
 
 Oorc_lab_field_dups <- labdup.prop.plot/indrep.prop.plot
 
 ggsave(file = "Oorc_RKW_lab_field_dups.png", Oorc_lab_field_dups, height = 6, width = 7, units = "in", bg = "white")
 

### dedup by Ind and date -----------------------------------------------------
Ind.date_keep <- as.data.frame(as.matrix(ps.raw@sam_data)) %>% 
  rownames_to_column("SampleID") %>% 
  mutate(readcount = phyloseq::sample_sums(ps.raw)) %>% 
  group_by(Individual.ID, year, month, day) %>% 
  slice_max(readcount) %>% 
  ungroup() %>% 
  pull(SampleID)

ps.raw <- prune_samples(ps.raw, samples = Ind.date_keep)
nsamples(ps.raw)

# track sample loss across both dedup steps
sample.track <- as.data.frame(allsamps) %>% 
  left_join(as.data.frame(LabID_keep) %>% mutate(keep = "Y"), 
            by = c("allsamps" = "LabID_keep")) %>% 
  mutate(keep = case_when(is.na(keep) ~ "N", TRUE ~ keep)) %>% 
  rename("LabID_dedup" = keep) %>% 
  left_join(as.data.frame(Ind.date_keep) %>% mutate(keep = "Y"),
            by = c("allsamps" = "Ind.date_keep")) %>% 
  mutate(keep = case_when(is.na(keep) ~ "N", TRUE ~ keep)) %>% 
  rename("Ind.date_dedup" = keep) %>% 
  rename("SampleID" = allsamps)

write.csv(sample.track, file = "Sample_loss_during_QAQC_newref_halibut_08312023.csv", row.names=FALSE)

# Plot all taxa by family before filtering taxa
order.palette <- c(pal_uchicago(alpha = 0.5)(9),pal_jama(alpha = 0.8)(5),"#FFFFFFCC")
order.plot <- 
  plot_bar(ps.raw %>% transform_sample_counts(function(x) x / sum(x)), fill = "Order") +
  theme_bw() +
  theme(panel.border = element_blank()) +
  #scale_fill_manual(values = order.palette) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.2)) +
  ylab("Prop. reads by order") +
  xlab("Samples") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  theme(strip.text = element_blank(),
        panel.margin.y = unit(0, "lines"))+
  theme(legend.position = "bottom") +
  scale_y_continuous(breaks = seq(0.2,1,0.2), expand = c(0,0)) 

### Filter taxa ----------------------------------------------------------------

# First, figure out how many NA taxa there are and what proportion of sequences they are
num.speciesID <- ps.raw %>% subset_taxa(is.na(Species)) %>% ntaxa()

num.NAspecies <- ps.raw %>% subset_taxa(!is.na(Species)) %>% ntaxa()

prop.NAspecies <- ps.raw %>% transform_sample_counts(function(x) x / sum(x)) %>% 
  subset_taxa(is.na(Species)) %>% 
  sample_sums() %>% as.data.frame() %>% 
  rename("prop" = 1) %>% 
  summarize(mean(prop), min(prop), max(prop), sd(prop))

# glom taxa to species at master level, remove Oorc sequences, convert reads to proportion
# this removes NA species from the dataset
ps.sp <- tax_glom(ps.raw, "Species") %>% 
  subset_taxa(Genus != "Orcinus") %>% 
  transform_sample_counts(function(x) x / sum(x)) 

# filter to remove low count taxa 
# must be at least 1% of diet in 1 or more samples
f1 <- filterfun_sample(function(x) x >= 0.01)
lowcount.filt <- genefilter_sample(ps.sp, f1, A=1)

ps.sp <- prune_taxa(lowcount.filt, ps.sp)

#glom taxa again but keep count data
# remove low count taxa according to above filter
ps.spcount <- tax_glom(ps.raw, "Species") %>% 
  subset_taxa(Genus != "Orcinus")

ps.spcount <- prune_taxa(lowcount.filt, ps.spcount) %>% 
  prune_samples(sample_names(ps.sp),.) %>% 
  subset_samples(Population != "AOKW")

#plot sampling effort by population, month, year
sampling.plot <- ggplot(data = ps.spcount@sam_data, aes(x = month, fill  = year)) +
  geom_bar(aes(y=(..count..)), 
           position = "stack") +
  facet_wrap(~Population) +
  scale_fill_aaas(alpha = 0.6) +
  ylab("Number of samples") +
  theme_minimal() +
  labs(fill="Year") +
  theme(legend.position = "bottom")

pdf(file = "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/05 Manuscript/figure_pdf/Oorc_sampling_effort.pdf",
    width = 10, height = 6)
sampling.plot
dev.off()

# Add sampling map to annual sampling effort figure 

Oorc_sampling_effort <- RKW_sample_map / sampling.plot +
  plot_layout(widths = c(5,4))

ggsave(file = "Oorc_RKW_sampling_effort.png", Oorc_sampling_effort, height = 8, width = 4.5, units = "in", bg = "white")


#repeat for ARKW
#sampling effort by Region, month, year
ARKWsampling.plot <- ggplot(data = subset(sample.data.frame(ps.sp), Population == "ARKW") %>% 
                              mutate(Region = case_when(Region == "West PWS" ~ "PWS",
                                                        TRUE ~ Region)), aes(x = month, fill  = year)) +
  geom_bar(aes(y=(..count..)), 
           position = "stack") +
  facet_wrap(~Region) +
  scale_fill_aaas(alpha = 0.6) +
  ylab("Number of samples") +
  theme_bw()

pdf(file = "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/05 Manuscript/figure_pdf/ARKWsampling_effort.pdf",
    width = 10, height = 6)
ARKWsampling.plot
dev.off()

save(ps.ampbias, ps.raw, ps.sp, ps.spcount, 
     file = "KWfecal_2016-2019_allsamps_tax95_newref_halibut_ampbias_overdisp_QAQC_08312023.Rdata")

save(indrep.prop.plot, sampling.plot,
     file = "Oorc_indrep_sampling_plots_tax95_newref_halibut_ampbias_overdisp_QAQC_08312023.Rdata")

ggsave("Oorc_sampling_plot.png", sampling.plot, width = 8, height = 6, units = "in")

ggsave("Oorc_ind_rep_plot.png", indrep.prop.plot, width = 8, height = 6, units = "in")

ggsave("Oorc_lab_rep_plot.png", labdup.prop.plot, width = 8, height = 6, units = "in")
