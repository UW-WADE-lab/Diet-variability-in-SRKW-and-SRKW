## dada2 QAQC of Oorc fecal fastq sequences
## 1/3/2022
## Amy Van Cise

### set up working environment

#devtools::install_github("benjjneb/dada2", ref="v1.16")
#BiocManager::install("dada2")
library(dada2)
library(tidyverse)
library(ggplot2)
library(seqinr)

fecal.seqs.file <- "/scratch/avancise/Oorc_fecals/fastq_files"
taxref <- "ref_test_add_halibut.fasta"


### read fastq files in working directory
fnFs <- sort(list.files(fecal.seqs.file, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(fecal.seqs.file, pattern="_R2_001.fastq", full.names = TRUE))
sample.names1 <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names2 <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)
sample.names <- paste(sample.names1, sample.names2, sep = "_")

### vizualize read quality profiles
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])

### Name filtered files in filtered/subdirectory
filtFs <- file.path(fecal.seqs.file, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(fecal.seqs.file, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


### Filter and Trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
### Dereplicate
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

### Learn Error Rates
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out

#plotErrors(dadaFs.lrn[[1]], nominalQ=TRUE)

### Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

### Merge Paired Reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

### Construct sequence table
seqtab <- makeSequenceTable(mergers)

### Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

freq.nochim <- sum(seqtab.nochim)/sum(seqtab)

### Track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

### Create Taxonomy file
# library(taxonomizr)
# 
## prepareDatabase('accessionTaxa.sql') #this downloads all of the taxonomy on NCBI. delete afterward
# 
# 
# blastResults<-read.delim("C:/Users/Amy.M.VanCise/Desktop/PreySpeciesandMarMamandSharksdbase.fasta",header=FALSE,stringsAsFactors=FALSE)
# 
# blast_results_tax <- data.frame(Accession = blastResults$V1[c(TRUE, FALSE)], 
#                                  Sequence = blastResults$V1[c(FALSE, TRUE)]) %>% 
#   separate(Accession, into = c(NA,NA,NA,"Accession","Species"), sep = "\\|") %>% 
#   mutate(taxID = accessionToTaxa(Accession,"accessionTaxa.sql")) %>% 
#   bind_cols(as.data.frame(getTaxonomy(blast_results$taxID,"accessionTaxa.sql")))
# 
# ref_fasta <- blast_results_tax %>% 
#   mutate(Taxonomy = paste(superkingdom, phylum, class, order, family, genus, species, sep = ";"), 
#          .before = Sequence, .keep = "unused") %>% 
#   select(-Accession,-Species,-taxID)
# 
# write.fasta(sequences = as.list(ref_fasta$Sequence), as.string = FALSE, names = ref_fasta$Taxonomy, file.out = "ref_test.fasta")
  
### Assign Taxonomy
taxa <- assignTaxonomy(seqtab.nochim, taxref, tryRC = TRUE, minBoot = 95, outputBootstraps = FALSE)

### Save data
save(seqtab.nochim, freq.nochim, track, taxa, file = "Oorc_2016-2021_dada2_QAQC_tax95_newref_halibut_output.Rdata")
