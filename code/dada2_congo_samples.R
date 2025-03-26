library(pacman)
p_load(tidyverse, patchwork, here, Biostrings, ShortRead, dada2, doParallel, DECIPHER, phangorn, update = F)

setwd("~/Desktop/drc_cvmb_paper")

files_path <- "/raw_reads"
head(list.files(files_path))

fnFs <- sort(list.files(files_path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(files_path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- fnFs %>% basename() %>% str_split("_") %>% map_chr(., 1) 

plotQualityProfile(fnFs[70:81])
plotQualityProfile(fnRs[1:10])
# make a string to subset
n <- 1:length(fnFs)
# Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ subdirectory
dir.create("filtered")
filtFs <- file.path(files_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(files_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# run the function on good sequences
filtered_out <- filterAndTrim(
  fnFs,
  # “forward_reads” - input
  filtFs,
  # “forward_reads” - output filtered
  fnRs,
  # “reverse_reads” - input
  filtRs,
  # “forward_reads” - output filtered
  maxEE =c(5,5),
  # quality filtering threshold
  rm.phix = TRUE,
  # removes any reads that match the PhiX bacteriophage genome
  #   minLen=30, # minimum length reads we want to keep after trimming
  #truncLen = c(300,300), #skipped as read quality was good! also this maintained the most amount of reads after filtering
  # minimum size to trim the forward and reverse reads (keep the QS above 30 overall) - Shorter sequences are discarded
  trimLeft = c(20,20), #remo those trailing primers at the beggining of the sequences
  compress = TRUE,
  # gzipped fastq files
  truncQ = 2,
  #  trims all bases after the first quality score of 2
  multithread = TRUE,
  # TRUE if working in a HPC, otherwise FALSE (e.g., Windows machine)
  verbose = T)

head(filtered_out)
df_reads3 <- filtered_out  %>% as.data.frame() %>% rownames_to_column("sampleID") %>% mutate(perc = (reads.out/reads.in)*100)

# create new folder
path.dada2 <- file.path(files_path, "dada2")
if(!dir.exists(path.dada2)) dir.create(path.dada2)
rm(path.dada2)

# save
saveRDS(filtered_out, paste0(files_path, "/dada2/filtered_trim.RDS"))
plotQualityProfile(filtFs[1:10])
plotQualityProfile(filtRs[1:10])

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# save RDS
saveRDS(errF, paste0(files_path, "/dada2/errF.RDS"))
saveRDS(errF, paste0(files_path, "/dada2/errR.RDS"))

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#derepliation
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

# save RDS
saveRDS(derep_forward, paste0(files_path, "/dada2/derep_forward.RDS"))
saveRDS(derep_reverse, paste0(files_path, "/dada2/derep_reverse.RDS"))

#ASV inference
dada_forward <- dada(derep_forward, err=errF, 
                     #  pool="pseudo", 
                     multithread=TRUE)
# save RDS
saveRDS(dada_forward, paste0(files_path, "/dada2/dada_forward.RDS"))
#dada_forward <- read_rds(paste0(files_path, "/dada2/dada_forward.RDS"))

dada_reverse <- dada(derep_reverse, err=errR, 
                     #pool="pseudo", 
                     multithread=TRUE)
# save RDS
saveRDS(dada_reverse, paste0(files_path, "/dada2/dada_reverse.RDS"))
#dada_reverse <- read_rds(paste0(files_path, "/dada2/dada_reverse.RDS"))

dada_forward[[6]]
dada_reverse[[6]]

#merging
merged_amplicons.default <- mergePairs(dada_forward, derep_forward, 
                                       dada_reverse, derep_reverse, 
                                       #  trimOverhang = T,
                                       #  minOverlap=30, 
                                       verbose = T)


merged_amplicons <- merged_amplicons.default
# merged_amplicons <- merged_amplicons.overhang.T
# merged_amplicons <- merged_amplicons.min30

head(merged_amplicons[[4]])
# this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
class(merged_amplicons) # list
length(merged_amplicons) # 81 elements in this list, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list 


# remove big objects
rm("derep_forward", 'derep_reverse')

#make the sequence table
seqtab.default <- makeSequenceTable(merged_amplicons)


# inspect the number of ASVs in each seq table
dim(seqtab.default) # samples 4534 ASVs

# get the chosen one
seqtab <- seqtab.default
class(seqtab) # matrix

# save RDS
saveRDS(merged_amplicons, paste0(files_path, "/dada2/merged_amplicons.RDS"))

#merge with sample 003
seqtab2 <- mergeSequenceTables(seqtab, seqtab_003)
#chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim) 
sum(seqtab.nochim)/sum(seqtab)  # 97%

saveRDS(seqtab.nochim, paste0(files_path, "/dada2/seqtab.nochim.RDS"))

getN <- function(x) sum(getUniques(x))


# making a little table


getN <- function(x) sum(getUniques(x))

summary_tab <- data.frame(row.names = sample.names, 
                          dada2_input = filtered_out[,1],
                          filtered = filtered_out[,2], 
                          input_percent = round((filtered_out[,2]/filtered_out[,1])*100,2),
                          dada_f = sapply(dada_forward, getN),
                          dada_r = sapply(dada_reverse, getN), 
                          merged = sapply(merged_amplicons, getN),
                          input_perc = round((sapply(merged_amplicons, getN)/filtered_out[,1])*100,2),
                          filtered_perc = round((sapply(merged_amplicons, getN)/filtered_out[,2])*100, 2),
                          nonchim = rowSums(seqtab.nochim),
                          final_perc_reads_retained = round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 2),
                          merged_perc = round(rowSums(seqtab.nochim)/sapply(merged_amplicons, getN)*100, 2))

summary_tab

# save
summary_tab %>% 
  rownames_to_column("sample_id") %>% 
  write_csv(paste0(getwd(), "/outputs/summary_tab_dada2.csv"))

#Assign taxonomy
tax_file <- "/Users/ndlovu.3/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Possible thesis projects/HIV project/Vaginal study w Ann&Ricardo/16S_data_vaginal_DNA/Attempt_2/scripts/greengenes2_trainset.fa.gz"
taxHS <- assignTaxonomy(seqtab.nochim, tax_file, taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), multithread = T, verbose = TRUE)
colnames(taxHS) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species")

# Save sequence and tax tables as TSV files
write.table(taxHS, file = "16S_taxa_SV_gg.tsv", sep = "\t", quote = FALSE)

# Save taxHS as RDS file
saveRDS(taxHS, "taxHS_gg2.RDS")
rownames(seqtab.nochim) <- sprintf("V%03d", as.integer(gsub("\\D", "", rownames(seqtab.nochim))))

meta <- read_excel("~/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Possible thesis projects/HIV project/Vaginal study w Ann&Ricardo/16S_data_vaginal_DNA/Attempt_2/final_otu_tables/metadata/metadata_congosamples.xlsx")
rownames(meta) <- meta$sampleID
# Create phyloseq object
ps_congo <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), tax_table(taxHS), sample_data(meta))

# Save phyloseq object as RDS fi
saveRDS(ps_congo, "ps_rds/congo2/ps_congo_raw.rds")


#Filtering phyloseq
ps_filt <- ps_congo_raw  %>% subset_taxa(Family != "f__mitochondria" | is.na(Family) & Class != "c__Chloroplast" | is.na(Class))

#Filter to remove tax with less than 10 reads
ps_filt2  <- prune_taxa(taxa_sums(ps_filt) > 10, ps_filt)
# Remove samples with less than 100 reads
ps_filt3  <- prune_samples(sample_sums(ps_filt2) >= 100, ps_filt2) #1770 samples
#remove taxa in less tna 0.5% abund across all samples
ps_filt4 <- prune_taxa(taxa_sums(ps_filt3) >= 0.005 * sample_sums(ps_filt3), ps_filt3)
#remove taxa not present more than 2 times in at least 10% of samples
ps_filt5<- filter_taxa(ps_filt4, function(x) sum(x > 2) > (0.1*length(x)), TRUE)
#. ps_spp3 <- prune_taxa(taxa_sums(ps_spp2) >= 0.01 * sample_sums(ps_spp2), ps_spp2)
# Remove doubletons and singletons
ps_filt6 <- filter_taxa(ps_filt5, function(x) { sum(x > 0) > 2 }, prune = TRUE)
# Remove taxa not assigned at phylum level
ps_filt7 <- subset_taxa(ps_filt6 , Phylum != "NA")
read_per_samplefilt <- as.data.frame(sample_sums(ps_filt7))
write.table(read_per_samplefilt, file = "readsperssample_congofilt.csv", sep = ",", quote = F)

#add refseq and phy_tree
seqs <- getSequences(otu_table(ps_filt7))
dups <- duplicated(seqs) #no duplicates
dups_indices <- which(dups)

names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA, verbose = F)

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

phy_tree(ps_filt7) <- fitGTR$tree
#add refseq
dna <- Biostrings::DNAStringSet(taxa_names(ps_filt7))
names(dna) <- taxa_names(ps_filt7)
ps_filt8 <- merge_phyloseq(ps_filt7, dna)
taxa_names(ps_filt8) <- paste0("ASV", seq(ntaxa(ps_filt8)))

#rename properly
otu_filt <- as.data.frame(otu_table(ps_filt8))
rownames(otu_filt) <- sprintf("V%03d", as.integer(gsub("\\D", "", rownames(otu_filt))))
otu_filt <- otu_table(otu_filt, taxa_are_rows = F)
otu_table(ps_filt8) <- otu_filt
#save RDS
saveRDS(ps_filt8, "ps_congo_filt.rds")

#load metadata file and clean up
metadata_congosamples <- read_excel("metadata/metadata_congosamples.xlsx")

metadata_congosamples[metadata_congosamples=="NA" ] <- NA
metadata_congosamples <- as.data.frame(metadata_congosamples)
rownames(metadata_congosamples)<- metadata_congosamples$sampleID 

cols_to_convert <- grep("^enrol", names(metadata_congosamples), value = TRUE)
cols_to_convert2 <- c("Gravidity", "SES_quartile")

metadata_congosamples <- metadata_congosamples %>%
  mutate_at(vars(cols_to_convert), as.numeric)

metadata_congosamples <- metadata_congosamples %>%
  mutate_at(vars(cols_to_convert2), as.character)
#add meta to ps object
sample_data(ps_congo_filt) <- metadata_congosamples        

#source functions
source("functions_drc_cvmb.R")
#agglomerate and save
ps_congo_filt_genus <- tax_glom(ps_congo_filt, "Genus", NArm = T)
ps_congo_filt_spp <- tax_glom(ps_congo_filt, "Species", NArm = T)
saveRDS(ps_congo_filt_genus, "ps_rds/ps_congo_filt_genus.rds")
saveRDS(ps_congo_filt_spp, "ps_rds/ps_congo_filt_spp.rds")
#transform all
transform_list <- c("ps_congo_filt_genus", "ps_congo_filt_spp", "ps_congo_filt")

for (ps_name in transform_list) {
  
  ps <- get(ps_name)
  
  ps_rclr <- robust_clr_transform(ps)
  rclr_ps_name <- paste0(ps_name, "_rclr")
  saveRDS(ps_rclr, file = paste0("ps_rds/transformed/" ,rclr_ps_name, ".rds"))
}



