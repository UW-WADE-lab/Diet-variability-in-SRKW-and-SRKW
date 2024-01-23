library(tidyverse)
library(ggmap)
library(ggOceanMaps)
library(ggspatial)
library(ggsci)
library(patchwork)

setwd("G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding/")

options(ggOceanMaps.userpath = "C:/Users/Amy Van Cise/Downloads/gebco_2023_sub_ice_topo/GEBCO_2023_sub_ice_topo.nc")

# species color palette
species.palette = c(pal_uchicago(alpha = 0.6)(9)[c(1,3:9)],pal_jama(alpha = 0.6)(6))

names(species.palette) <- c("arrowtooth", "chum",
                            "coho", "chinook",
                            "steelhead", "sablefish", "lingcod",  "halibut","big skate", "English sole",
                            "Pacific herring", 
                            "sockeye", "pink", "None")

species.shapes = c(3,15,16,17,8,25,23,13)
names(species.shapes) <- c("arrowtooth", "chum",
                            "coho", "chinook",
                            "sablefish", "halibut","big skate", "None")

months.palette = c(pal_uchicago(alpha = 0.6)(9)[c(1,3:9)],pal_jama(alpha = 0.6)(4))
names(months.palette) <- (1:12)

# get data --------------------------------------------------------------
ARKW_sample_locations <- read.csv("new_prey_meta_7.26.22.csv") %>% 
  select(-latitudeN, -longitudeW) %>% 
  left_join(read.csv("new_prey_meta_7.26.22_CE_KP.csv") %>% select(Sample, latitudeN, longitudeW), 
            by = c("Sample" = "Sample")) %>% 
  filter(Population == "ARKW") %>% 
  filter(!(is.na(longitudeW))) %>% 
  filter(!(is.na(latitudeN))) %>% 
  mutate(longitudeW =longitudeW*(-1))

SRKW_sample_locations <- read.csv("new_prey_meta_7.26.22.csv") %>% 
  select(-latitudeN, -longitudeW) %>% 
  left_join(read.csv("new_prey_meta_7.26.22_CE_KP.csv") %>% select(Sample, latitudeN, longitudeW), 
            by = c("Sample" = "Sample")) %>% 
  filter(Population == "SRKW") %>% 
  filter(!(is.na(latitudeN))) %>% 
  mutate(longitudeW = longitudeW*(-1))

### ARKW basemap ----------------------------------------------------------------
ARKW_sample_map <- basemap(limits = c(min(ARKW_sample_locations$longitudeW)-0.5,
                                      max(ARKW_sample_locations$longitudeW)+0.5,
                                      min(ARKW_sample_locations$latitudeN)-0.5,
                                      max(ARKW_sample_locations$latitudeN)+0.5),
                           bathy.style = "rub", rotate = TRUE, grid.col = NA) +
  ggspatial::geom_spatial_point(data = ARKW_sample_locations, aes(x = longitudeW, 
                                                                  y = latitudeN, 
                                                                  color = as.factor(month), 
                                                                  shape = prey50),
                                size = 3) +
  geom_jitter()+
  scale_color_manual(values = months.palette, name = "Month") + 
  scale_shape_manual(values=species.shapes) +
  labs(shape = "Dominant prey species (>50%)") +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank()) +
  ylab("Latitude") +
  xlab("Longitude") +
  #scale_fill_viridis_c("Water depth (m)") +
  ggspatial::annotation_scale(location = "br") + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true") +
  guides(color = guide_legend(override.aes = list(fill = NA, order = NA), order = 1),
         shape = guide_legend(override.aes = list(fill = NA, order = NA), order = 2)) +
  theme(legend.key = element_rect(fill = "white"))

ggsave(file = "Oorc_ARKW_sample_map.png", ARKW_sample_map, height = 6, width = 8, units = "in", bg = "white")


### SRKW sample map -----------------------------------------------------------------------
SRKW_sample_map <- basemap(limits = c(min(SRKW_sample_locations$longitudeW)-0.5,
                                      max(SRKW_sample_locations$longitudeW)+0.5,
                                      min(SRKW_sample_locations$latitudeN)-0.5,
                                      max(SRKW_sample_locations$latitudeN)+0.5),
                           bathy.style = "rub", rotate = TRUE, grid.col = NA) +
  ggspatial::geom_spatial_point(data = SRKW_sample_locations, aes(x = longitudeW, y = latitudeN, 
                                                                  color = as.factor(month), shape = prey50),
                                size = 3) +
  geom_jitter() +
  scale_color_manual(values = months.palette, name = "Month") +
  scale_shape_manual(values=species.shapes) +
  labs(shape = "Dominant prey species (>50%)") +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank()) +
  ylab("Latitude") +
  xlab("Longitude") +
  #scale_fill_viridis_c("Water depth (m)") +
  ggspatial::annotation_scale(location = "bl") + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true") +
  guides(color = guide_legend(override.aes = list(fill = NA, order = NA), order = 1),
         shape = guide_legend(override.aes = list(fill = NA, order = NA), order = 2)) +
  theme(legend.key = element_rect(fill = "white"))

ggsave(file = "Oorc_SRKW_sample_map.png", SRKW_sample_map, height = 6, width = 8, units = "in", bg = "white")

### Combine plots -------------------------------------------------------------

RKW_dom_species_map <- ARKW_sample_map + SRKW_sample_map

ggsave(file = "Oorc_RKW_dom_species.png", RKW_dom_species_map, height = 7, width = 16, units = "in", bg = "white")
