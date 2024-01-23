### Average diet proportions over time, by month and year
## Amy Van Cise
## March 2022

library(tidyverse)
library(ggsci)
library(phyloseq)
library(phyloseqCompanion)
library(vegan)
library(cowplot)

### Set up environment ---------------------------------------------------------------------------
wd <- "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding"
setwd(wd)

load("Oorc_95tax_newref_halibut_phyloseq_majortaxfilt.Rdata")

mycolors <- pal_futurama()(12)[c(1:4,8,9)]

### Convert phlyoseq object to dataframe ------------------------------------------
species_prop <- ps.sp@otu_table %>% as.data.frame()
names(species_prop) <- ps.sp@tax_table[,7]
samdf <- sample.data.frame(ps.sp)

# 10/26 new sample data with 6 new pod IDs for SRKW
samdf <- read.csv("new_prey_meta_10.26.23_pedPod.csv") %>% 
  select(-X) %>% 
  mutate(year = as.factor(year)) %>% 
  mutate(month = as.factor(month)) %>% 
  unite("Month.Year", month:year, sep = ".", remove = FALSE) %>% 
  unite("Pop.year", c(Population,year), sep = ".", remove = FALSE) %>% 
  unite("Pop.month", c(Population,month), sep = ".", remove = FALSE) %>% 
  unite("Pop.month.year", c(Population,month,year), sep = ".", remove = FALSE) %>% 
  unite("Region.month", c(Region,month), sep = ".", remove = FALSE) %>% 
  filter(Population != "NRKW") %>% 
  mutate(newpedPod = case_when(is.na(pedPod)~Pod,
                               TRUE~pedPod), .after = pedPod) %>% 
  select(-pedPod) %>% 
  rename("pedPod" = newpedPod)

species_prop_meta <- species_prop %>% 
  rownames_to_column("Sample") %>% 
  left_join(samdf, by = c("Sample" = "Sample")) %>% 
  select(-Region.month, -Pop.year, -Pop.month, -Pop.month.year, -dom.sp, -prey50) %>% 
  relocate(10:length(.), .after = Sample)

species_prop_meta_long <- species_prop_meta %>% 
  pivot_longer(18:length(.), names_to = "Species", values_to = "Proportion") %>% 
  mutate(species_f = factor(Species, levels = c("Oncorhynchus tshawytscha", "Oncorhynchus keta",
                                                  "Oncorhynchus kisutch", "Oncorhynchus mykiss", 
                                                  "Hippoglossus stenolepis","Atheresthes stomias",
                                                  "Anoplopoma fimbria","Ophiodon elongatus", 
                                                "Raja binoculata","Parophrys vetulus")))
#facet label names
month.labs <- c("January", "February","March","April","May", "June", "July", "August","Sept","Oct","Nov","Dec")
names(month.labs) <- as.factor(c(1,2,3,4,5,6,7,8,9,10,11,12))
# species.labs <- c("arrowtooth", "halibut", "chum", "coho", "Chinook","steelhead","sablefish","lingcod","big skate", "English sole")
# names(species.labs) <- c("Atheresthes stomias", "Hippoglossus stenolepis", "Oncorhynchus keta", 
#                          "Oncorhynchus kisutch", "Oncorhynchus tshawytscha",
#                          "Oncorhynchus mykiss", "Anoplopoma fimbria","Ophiodon elongatus", "Raja binoculata","Parophrys vetulus") 


species.labs <- c("arrowtooth", "chum", "coho", "Chinook","steelhead","sablefish","lingcod", "halibut")
names(species.labs) <- c("Atheresthes stomias", "Oncorhynchus keta", 
                         "Oncorhynchus kisutch", "Oncorhynchus tshawytscha",
                         "Oncorhynchus mykiss", "Anoplopoma fimbria","Ophiodon elongatus", "Hippoglossus stenolepis") 


# plot yearly proportional abundance
annual.prop.plot <- ggplot(species_prop_meta_long, aes(x = as.numeric(as.character(year)), y = Proportion, color = Population, fill = Population)) +
  geom_point() +
  geom_smooth(method = "glm", alpha = 0.6)+
  facet_grid(Species~month, labeller = labeller(Species = species.labs, month = month.labs), scales = 'free_y')+
  xlab("Year") +
  ylab("Proportion of Diet") +
  theme_bw()+
  scale_fill_manual(values =  c("cornflowerblue","purple"))+
  scale_color_manual(values =  c("cornflowerblue","purple"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5)) 

# plot monthly proportional abundance
monthly.prop.plot <- ggplot(species_prop_meta_long, aes(x = as.numeric(as.character(month)), y = Proportion, color = Population, fill = Population)) +
  geom_point() +
  geom_smooth(method = "glm", alpha = 0.6)+
  facet_grid(Species~year, labeller = labeller(Species = species.labs), scales = 'free_y')+
  xlab("Month") +
  ylab("Proportion of Diet") +
  theme_bw()+
  scale_color_manual(values=mycolors) +
  scale_fill_manual(values=mycolors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))+
  scale_x_continuous(name = "Month", breaks = c(3,5,6,7,8,9,10,11,12))

### plot monthly proportional abundance aggregated across all years -----------------------
salmon_species = c("Oncorhynchus keta", 
                   "Oncorhynchus kisutch", "Oncorhynchus tshawytscha",
                   "Oncorhynchus mykiss")

nonsalmon_species = c("Atheresthes stomias", "Hippoglossus stenolepis", "Anoplopoma fimbria","Ophiodon elongatus", "Raja binoculata","Parophrys vetulus")

species_prop_salmon <- species_prop_meta_long %>% 
  filter(species_f %in% salmon_species)

species_prop_nonsalmon <- species_prop_meta_long %>% 
  filter(species_f %in% nonsalmon_species)
  
monthly.prop.agg.salmon <- ggplot(species_prop_salmon, aes(x = as.numeric(as.character(month)), y = Proportion, color = Population, fill = Population)) +
  geom_jitter(width = 0.1, size = 2) +
  geom_smooth(method = "loess", span = 0.85, alpha = 0.2)+
  facet_grid(species_f~., labeller = labeller(species_f = species.labs), scales = 'free_y')+
  xlab("Month") +
  ylab("Proportion of Diet") +
  theme_bw()+
  scale_color_manual(values=mycolors[3:4]) +
  scale_fill_manual(values=mycolors[3:4]) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))+
  scale_x_continuous(name = "Month", breaks = c(1,2,3,4,5,6,7,8,9,10,11,12))+
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(breaks = seq(0.25,1,0.5)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

seasonal.trends.legend <- ggpubr::get_legend(monthly.prop.agg.salmon)

monthly.prop.agg.salmon <- monthly.prop.agg.salmon +
  theme(legend.position = "none")

monthly.prop.agg.nonsalmon <- ggplot(species_prop_nonsalmon, aes(x = as.numeric(as.character(month)), y = Proportion, color = Population, fill = Population)) +
  geom_jitter(width = 0.1, size = 2) +
  geom_smooth(method = "loess", span = 0.65, alpha = 0.2)+
  facet_grid(species_f~., labeller = labeller(species_f = species.labs), scales = 'free_y')+
  xlab("Month") +
  ylab("") +
  theme_bw()+
  scale_color_manual(values=mycolors[3:4]) +
  scale_fill_manual(values=mycolors[3:4]) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))+
  scale_x_continuous(name = "Month", breaks = c(1,2,3,4,5,6,7,8,9,10,11,12))+
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(breaks = seq(0.25,1,0.5)) +
  theme(legend.position = "none")

seasonal_plots <- plot_grid(monthly.prop.agg.salmon, monthly.prop.agg.nonsalmon, ncol = 2)

seasonal_plots_legend <- plot_grid(seasonal.trends.legend, seasonal_plots, ncol = 1, rel_heights = c(0.5,5))


pdf(file = "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/05 Manuscript/figure_pdf/Oorc_temporal_prop_population.pdf",
    width = 10, height = 4)
seasonal_plots_legend
dev.off()

### PERMANOVA of diet composition across pods -----------------------------------
podsamps.env <- samdf %>% 
  filter(pedPod %in% c("J", "K", "L")) %>% 
  filter(Sample %in% sample_names(ps.sp)) %>% 
  #filter(month %in% c(10,11,12,1,2,3)) %>% 
  column_to_rownames("Sample")

podsamps.otu <- as.data.frame(ps.sp@otu_table) %>% 
  filter(rownames(.) %in% rownames(podsamps.env))

pod_difference <- adonis2(podsamps.otu ~ pedPod + month + year, data = podsamps.env)
pod_difference

### Plot monthly proportional abundance aggregated by year across all months by pod ------------

species_prop_meta_long_SRKW <- species_prop_meta_long %>% 
  filter(Population == "SRKW") %>% 
  filter(pedPod != "unknown") %>% 
  filter(pedPod %in% c("J", "K", "L"))

species_prop_salmon_SRKW <- species_prop_meta_long_SRKW %>% 
  filter(species_f %in% salmon_species)

species_prop_nonsalmon_SRKW <- species_prop_meta_long_SRKW %>% 
  filter(species_f %in% nonsalmon_species)

monthly.prop.agg.salmon.SRKW <- ggplot(species_prop_salmon_SRKW, aes(x = as.numeric(as.character(month)), y = Proportion, color = pedPod, fill = pedPod)) +
  geom_jitter(width = 0.1, size = 2) +
  geom_smooth(method = "loess", span = 0.85, alpha = 0.2)+
  facet_grid(species_f~., labeller = labeller(species_f = species.labs), scales = 'free_y')+
  xlab("Month") +
  ylab("Proportion of Diet") +
  theme_bw()+
  scale_color_manual(values=mycolors[c(3,5,8)]) +
  scale_fill_manual(values=mycolors[c(3,5,8)]) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))+
  scale_x_continuous(name = "Month", breaks = c(1,2,3,4,5,6,7,8,9,10,11,12))+
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(breaks = seq(0.25,1,0.5)) +
  theme(legend.position="top", legend.title = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

seasonal.trends.legend.SRKW <- ggpubr::get_legend(monthly.prop.agg.salmon.SRKW)

monthly.prop.agg.salmon.SRKW <- monthly.prop.agg.salmon.SRKW +
  theme(legend.position = "none")

monthly.prop.agg.nonsalmon.SRKW <- ggplot(species_prop_nonsalmon_SRKW, aes(x = as.numeric(as.character(month)), y = Proportion, color = pedPod, fill = pedPod)) +
  geom_jitter(width = 0.1, size = 2) +
  geom_smooth(method = "loess", span = 0.85, alpha = 0.2)+
  facet_grid(species_f~., labeller = labeller(species_f = species.labs), scales = 'free_y')+
  xlab("Month") +
  ylab("Proportion of Diet") +
  theme_bw()+
  scale_color_manual(values=mycolors[c(3,5,8)]) +
  scale_fill_manual(values=mycolors[c(3,5,8)]) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))+
  scale_x_continuous(name = "Month", breaks = c(1,2,3,4,5,6,7,8,9,10,11,12))+
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(breaks = seq(0.25,1,0.5)) +
  theme(legend.position = "none")

seasonal_plots_SRKWpod <- plot_grid(monthly.prop.agg.salmon.SRKW, monthly.prop.agg.nonsalmon.SRKW, ncol = 2)

seasonal_plots_legend_SRKWpod <- plot_grid(seasonal.trends.legend.SRKW, seasonal_plots_SRKWpod, ncol = 1, rel_heights = c(0.5,5))


pdf(file = "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/05 Manuscript/figure_pdf/Oorc_temporal_prop_SRKW_pod.pdf",
    width = 10, height = 4)
seasonal_plots_legend_SRKWpod
dev.off()

### Plot proportional abundance across years by species for the months with most years sampled -------------------------

#filtered dataset
species_prop_meta_long_annual <- species_prop_meta_long %>% 
  filter(Population == "SRKW" & as.character(month) == "9" |
           Population == "SRKW" & as.character(month) == "11" |
           Population == "ARKW" & as.character(month) == "6") %>% 
  filter(!is.na(species_f)) 

# PERMANOVA test of differences among years
SRKW_annual_sept <- species_prop_meta_long_annual %>% 
  select(-species_f) %>% 
  pivot_wider(names_from = "Species", values_from = "Proportion") %>% 
  filter(as.character(month) == "9") %>% 
  select(-(2:17)) %>% 
  column_to_rownames("Sample")

SRKW_annual_sept_env <- species_prop_meta_long_annual %>% 
  select(-species_f) %>% 
  pivot_wider(names_from = "Species", values_from = "Proportion") %>% 
  filter(as.character(month) == "9") %>% 
  select((1:17)) %>% 
  column_to_rownames("Sample")

sept_difference <- adonis2(SRKW_annual_sept ~ year, data = SRKW_annual_sept_env)

SRKW_annual_nov <- species_prop_meta_long_annual %>% 
  select(-species_f) %>% 
  pivot_wider(names_from = "Species", values_from = "Proportion") %>% 
  filter(as.character(month) == "11") %>% 
  select(-(2:17)) %>% 
  column_to_rownames("Sample")

SRKW_annual_nov_env <- species_prop_meta_long_annual %>% 
  select(-species_f) %>% 
  pivot_wider(names_from = "Species", values_from = "Proportion") %>% 
  filter(as.character(month) == "11") %>% 
  select((1:17)) %>% 
  column_to_rownames("Sample")

nov_difference <- adonis2(SRKW_annual_nov ~ year, data = SRKW_annual_nov_env)

ARKW_annual_jun <- species_prop_meta_long_annual %>% 
  select(-species_f) %>% 
  pivot_wider(names_from = "Species", values_from = "Proportion") %>% 
  filter(as.character(month) == "6") %>% 
  select(-(2:17)) %>% 
  column_to_rownames("Sample")

ARKW_annual_jun_env <- species_prop_meta_long_annual %>% 
  select(-species_f) %>% 
  pivot_wider(names_from = "Species", values_from = "Proportion") %>% 
  filter(as.character(month) == "6") %>% 
  select((1:17)) %>% 
  column_to_rownames("Sample")

jun_difference <- adonis2(ARKW_annual_jun ~ year, data = ARKW_annual_jun_env)

# Bar plot proportion among years
species_prop_meta_long_annual_summary <- species_prop_meta_long_annual %>% 
  group_by(Population, year, month, species_f) %>% 
  summarize(mean_prop = mean(Proportion)) %>% 
  ungroup()

facet_labs <- c("ARKW, June", "SRKW, September", "SRKW, November")
names(facet_labs) <- as.factor(c(6,9,11))

annual.prop.agg.plot.monthly <- ggplot(species_prop_meta_long_annual_summary, aes(x = as.character(year), y = mean_prop, color = species_f, fill = species_f)) +
  geom_col(position="stack") +
  xlab("Year") +
  ylab("Proportion of Diet") +
  theme_bw()+
  facet_grid(~month, scales = "free_x", labeller = labeller(month = facet_labs)) + 
  scale_color_manual(limits = names(species.labs), labels = species.labs, values=c(pal_uchicago(alpha = 0.6)(9)[c(1,3:9)],pal_jama(alpha = 0.6)(6))) +
  scale_fill_manual(limits = names(species.labs), labels = species.labs, values=c(pal_uchicago(alpha = 0.6)(9)[c(1,3:9)],pal_jama(alpha = 0.6)(6))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.5))+
  scale_x_discrete(name = "Year", breaks = c(2014, 2016,2017,2018,2019,2020,2021))+
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(breaks = seq(0.25,1,0.5)) +
  theme(legend.position = "bottom") +
  labs(color = "", fill = "")

pdf(file = "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/05 Manuscript/figure_pdf/Oorc_annual_sept_prop_SRKW.pdf",
    width = 10, height = 4)
annual.prop.agg.plot.monthly
dev.off()

save(pod_difference, seasonal_plots_legend, seasonal_plots_legend_SRKWpod, annual.prop.agg.plot.monthly, 
     jun_difference, SRKW_annual_sept_env, SRKW_annual_nov_env, ARKW_annual_jun_env,
     nov_difference, sept_difference, species_prop_meta_long_SRKW, file = "Oorc_temporal_prop_plots.Rdata")

ggsave("Oorc_seasonal_plots_legend.png",seasonal_plots_legend, height = 4, width = 8, units = "in")
ggsave("Oorc_seasonal_plots_legend_SRKW.png", seasonal_plots_legend_SRKWpod, height = 4, width = 8, units = "in")
ggsave("Oorc_annual_plots_Sept_SRKW.png", annual.prop.agg.plot.monthly, height = 4, width = 8, units = "in")
