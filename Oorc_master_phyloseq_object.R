### Create master phyloseq object
## Amy Van Cise
## March 2022

library(tidyverse)
library(ggsci)
library(phyloseq)
library(microViz)

wd <- "G:/My Drive/02 NWFSC postdoc research/01 Oo fecal DNA/04 Data analysis/prey metabarcoding"

setwd(wd)

mycolors <- pal_futurama()(12)[c(1:4,8,9)]

### load dada2 output
load("Oorc_2016-2021_dada2_QAQC_tax95_newref_halibut_output.Rdata")

### get sample metadata
samdf <- read.csv("new_prey_meta_7.26.22.csv") %>% 
  mutate(year = as.factor(year)) %>% 
  mutate(month = as.factor(month)) %>% 
  unite("Month.Year", month:year, sep = ".", remove = FALSE) %>% 
  unite("Pop.year", c(Population,year), sep = ".", remove = FALSE) %>% 
  unite("Pop.month", c(Population,month), sep = ".", remove = FALSE) %>% 
  unite("Pop.month.year", c(Population,month,year), sep = ".", remove = FALSE) %>% 
  unite("Region.month", c(Region,month), sep = ".", remove = FALSE) %>% 
  filter(Population != "NRKW") %>% 
  column_to_rownames("Sample")

### create prey species factor object
species.palette = c(pal_uchicago(alpha = 0.6)(9),pal_jama(alpha = 0.6)(5))
names(species.palette) <- c("Atheresthes stomias", "Hippoglossus stenolepis", "Oncorhynchus keta", 
                  "Oncorhynchus kisutch", "Oncorhynchus tshawytscha",
                  "Oncorhynchus mykiss", "Anoplopoma fimbria","Ophiodon elongatus", "Raja binoculata","Parophrys vetulus",
                  "Clupea pallasii", "Salmo salar",
                  "Oncorhynchus nerka", "Oncorhynchus gorbuscha")

### create master phyloseq object
ps.raw <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))

### shorten ASV seq names, store sequences as reference
dna <- Biostrings::DNAStringSet(taxa_names(ps.raw))
names(dna) <- taxa_names(ps.raw)
ps.raw <- merge_phyloseq(ps.raw, dna)
taxa_names(ps.raw) <- paste0("ASV", seq(ntaxa(ps.raw)))

nsamples(ps.raw)

### export csv for ampbias correction
readcount.table <- as.data.frame(otu_table(ps.raw))
taxon.table <- as.data.frame(tax_table(ps.raw)) %>% rownames_to_column(var = "ASV")
metadata.table <- samdf %>% rownames_to_column(var = "Sample")
reference.table <- as.data.frame(refseq(ps.raw)) %>% rownames_to_column("ASVname") %>% 
  rename("DNAseq" = x)

readcount.table.long <- readcount.table %>% 
  rownames_to_column(var = "Sample") %>% 
  pivot_longer(-Sample, names_to = "ASV", values_to = "read.count") %>% 
  left_join(metadata.table, by = c("Sample" = "Sample")) %>% 
  select(Sample, LabID, Sample.ID, Population, year, month, Pod, read.count, ASV) %>% 
  left_join(taxon.table, by = c("ASV" = "ASV")) %>% 
  select(Sample, LabID, Sample.ID, Population, year, month, Pod,read.count,Species, ASV) %>% 
  mutate(Species = case_when(is.na(Species) ~ ASV, TRUE ~ Species)) %>% 
  select(-ASV) %>% 
  group_by(Sample, Species) %>% 
  mutate(sp.read.count = sum(read.count)) %>% 
  slice_head() %>% 
  select(-read.count)

readcount.table.species <- readcount.table.long %>% 
  filter(!grepl("ASV",Species)) %>% 
  pivot_wider(names_from = "Species", values_from = sp.read.count) %>% 
  rename("FieldID" = Sample.ID) %>% 
  mutate("bio_rep" = 1, .after = FieldID)

#write.csv(readcount.table.species, file = "testtax50.csv")

readcount.table.asv <- readcount.table.long %>% 
  filter(grepl("ASV",Species)) %>% 
  pivot_wider(names_from = "Species", values_from = sp.read.count) %>% 
  rename("FieldID" = Sample.ID) %>% 
  mutate("bio_rep" = 1, .after = FieldID)

# Export for amplification bias correction
write.csv(readcount.table.species, "KW_fecal_raw_readcount_allsamps_allspecies_tax95_newref_halibut_08312023.csv", row.names = FALSE)
write.csv(readcount.table.asv, "KW_fecal_raw_readcount_allsamps_asvs_tax95_newref_halibut_08312023.csv", row.names = FALSE)
save(reference.table, file = "ASVreference_allsamps_tax95_newref_halibut_08312023.Rdata")
