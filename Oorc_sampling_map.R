library(tidyverse)
library(ggmap)
library(ggOceanMaps)
library(ggspatial)
library(ggsci)
library(sf)
library(ggnewscale)
library(patchwork)


setwd("G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding/")

options(ggOceanMaps.userpath = "C:/Users/Amy/Downloads/gebco_2023_sub_ice_topo/GEBCO_2023_sub_ice_topo.nc")

# species color palette
species.palette = c(pal_uchicago(alpha = 0.6)(9)[c(1,3:9)],pal_jama(alpha = 0.6)(6))

names(species.palette) <- c("arrowtooth", "chum",
                            "coho", "chinook",
                            "steelhead", "sablefish", "lingcod",  "halibut","big skate", "English sole",
                            "Pacific herring", 
                            "sockeye", "pink", "None")

# get data --------------------------------------------------------------
ARKW_sample_locations <- read.csv("new_prey_meta_7.26.22.csv") %>% 
  select(-latitudeN, -longitudeW) %>% 
  left_join(read.csv("new_prey_meta_7.26.22_CE_KP.csv") %>% select(Sample, latitudeN, longitudeW), 
            by = c("Sample" = "Sample")) %>% 
  filter(Population == "ARKW") %>% 
  filter(!(is.na(longitudeW))) %>% 
  filter(!(is.na(latitudeN))) %>% 
  mutate(longitudeW =longitudeW*(-1))

ARKW_limits = data.frame(long = c(min(ARKW_sample_locations$longitudeW)-0.5,
                                  min(ARKW_sample_locations$longitudeW)-0.5,
                                  max(ARKW_sample_locations$longitudeW)+0.5,
                                  max(ARKW_sample_locations$longitudeW)+0.5),
                         lat = c(min(ARKW_sample_locations$latitudeN)-0.5,
                                 max(ARKW_sample_locations$latitudeN)+0.5,
                                 max(ARKW_sample_locations$latitudeN)+0.5,
                                 min(ARKW_sample_locations$latitudeN)-0.5
                                 ))

SRKW_sample_locations <- read.csv("new_prey_meta_7.26.22.csv") %>% 
  select(-latitudeN, -longitudeW) %>% 
  left_join(read.csv("new_prey_meta_7.26.22_CE_KP.csv") %>% select(Sample, latitudeN, longitudeW), 
            by = c("Sample" = "Sample")) %>% 
  filter(Population == "SRKW") %>% 
  filter(!(is.na(latitudeN))) %>% 
  mutate(longitudeW = longitudeW*(-1))

SRKW_limits = data.frame(long = c(min(SRKW_sample_locations$longitudeW)-0.5,
                                  min(SRKW_sample_locations$longitudeW)-0.5,
                                  max(SRKW_sample_locations$longitudeW)+0.5,
                                  max(SRKW_sample_locations$longitudeW)+0.5),
                         lat = c(min(SRKW_sample_locations$latitudeN)-0.5,
                                 max(SRKW_sample_locations$latitudeN)+0.5,
                                 max(SRKW_sample_locations$latitudeN)+0.5,
                                 min(SRKW_sample_locations$latitudeN)-0.5
                                 ))

RKW_dist <- st_read("G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding/RKW spatial polygons/resident_all-alb27.shp") %>% 
  mutate(COMMONCODE = case_when(COMMONCODE %in% c("PWKW", "SEKW") ~ "ARKW", TRUE ~ COMMONCODE)) %>% 
  filter(COMMONCODE %in% c("ARKW", "SRKW"))

### RKW basemap -----------------------------------------------------------------

RKW_dist_map <- basemap(RKW_dist, rotate = TRUE, bathy.style = "rcb", grid.col = NA) +
  new_scale_fill() +
  geom_sf(data=RKW_dist, color = "white", aes(fill = COMMONCODE)) +
  scale_fill_aaas() +
  ggspatial::geom_spatial_polygon(
    data = ARKW_limits, 
    aes(x = long, y = lat), fill = NA, color = "black", size = 1) +
  ggspatial::geom_spatial_polygon(
    data = SRKW_limits, 
    aes(x = long, y = lat), fill = NA, color = "black", size = 1) +
  theme(legend.position = "none")

  
### ARKW sample map ----------------------------------------------------------------
ARKW_sample_map <- basemap(limits = c(min(ARKW_sample_locations$longitudeW)-1,
                    max(ARKW_sample_locations$longitudeW)+1,
                    min(ARKW_sample_locations$latitudeN)-1,
                    max(ARKW_sample_locations$latitudeN)+1),
        bathy.style = "rcb", rotate = TRUE, grid.col = NA) +
  ggspatial::geom_spatial_point(data = ARKW_sample_locations, aes(x = longitudeW, y = latitudeN), size = 2) +
  geom_jitter()+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.margin = unit(c(0, 0, 0, 0), "null"),
        panel.border = element_rect(linewidth = 1.5))
  


### SRKW sample map -----------------------------------------------------------------------
SRKW_sample_map <- basemap(limits = c(min(SRKW_sample_locations$longitudeW)-0.5,
                                      max(SRKW_sample_locations$longitudeW)+0.5,
                                      min(SRKW_sample_locations$latitudeN)-0.5,
                                      max(SRKW_sample_locations$latitudeN)+0.5),
                           bathy.style = "rcb", rotate = TRUE, grid.col = NA) +
  ggspatial::geom_spatial_point(data = SRKW_sample_locations, aes(x = longitudeW, y = latitudeN),
                                size = 2) +
  geom_jitter() +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),  
        axis.ticks.y=element_blank(),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        panel.margin = unit(c(0, 0, 0, 0), "null"),
        panel.border = element_rect(linewidth = 1.5))




### final sampling map --------------------------------------------------------
RKW_sample_map <- RKW_dist_map + inset_element(SRKW_sample_map, 0.5, 0, 0.8, 0.5) + inset_element(ARKW_sample_map, 0.28, 0.5, 0.58, 0.8)

ggsave(file = "Oorc_RKW_sample_map.png", RKW_sample_map, height = 6, width = 8, units = "in", bg = "white")

save(RKW_sample_map, file = "RKW_sample_map.Rdata")

### combine with sampling effort -----------------------------------------------

load("Oorc_indrep_sampling_plots_tax95_newref_halibut_ampbias_overdisp_QAQC_08312023.Rdata")

Oorc_sampling_effort <- RKW_sample_map / sampling.plot +
  plot_layout(widths = c(5,4))

ggsave(file = "Oorc_RKW_sampling_effort.png", Oorc_sampling_effort, height = 8, width = 4.5, units = "in", bg = "white")
