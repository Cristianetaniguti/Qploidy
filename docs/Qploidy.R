## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# # install.packages("devtools")
# devtools::install_github("cristianetaniguti/Qploidy")

## -----------------------------------------------------------------------------
library(Qploidy)

# Simulate a VCF file with 500 markers and 60 samples
simulate_vcf(
  seed = 1234, file_path = "vcf_example_simulated.vcf",
  n_tetraploid = 35, n_diploid = 5, n_triploid = 10,
  n_markers = 500
)

## -----------------------------------------------------------------------------
data <- qploidy_read_vcf(vcf_file = "vcf_example_simulated.vcf")
head(data)

## ----echo=FALSE---------------------------------------------------------------
temp1 <- tempfile(fileext = ".txt")
simulate_axiom_summary(seed = 1234, file_path = temp1)
temp1 <- read.table(temp1, header = TRUE)
head(temp1)

## ----eval=FALSE---------------------------------------------------------------
# data <- read_axiom("AxiomGT1_summary.txt")
# head(data)

## ----eval=FALSE---------------------------------------------------------------
# data <- read_illumina_array(
#   "data/illumina/Set1.txt", # This is a example code, data not available
#   "data/illumina/Set2.txt",
#   "data/illumina/Set3.txt"
# )
# head(input_data)

## -----------------------------------------------------------------------------
# All samples
all_samples <- unique(data$SampleName)
all_samples
tetraploid_samples <- all_samples[grep("Tetraploid", all_samples)]
tetraploid_samples

## ----eval=FALSE---------------------------------------------------------------
# data_reference <- data

## -----------------------------------------------------------------------------
genos <- qploidy_read_vcf("vcf_example_simulated.vcf", geno = TRUE)
dim(genos)

genos.pos <- qploidy_read_vcf("vcf_example_simulated.vcf", geno.pos = TRUE)
head(genos.pos)

## ----eval=FALSE---------------------------------------------------------------
# genos.pos$MarkerName <- paste0(genos.pos$Chromosome, "_", genos.pos$Position)
# head(genos.pos)

## -----------------------------------------------------------------------------
genos <- genos[which(genos$SampleName %in% tetraploid_samples), ]

## -----------------------------------------------------------------------------
library(tidyr)

# Prepare inputs for updog
## Approach (1)
# tobe_genotyped <- data[which(data$SampleName %in% tetraploid_samples),]

## Approach (2)
tobe_genotyped <- data

ref <- pivot_wider(tobe_genotyped[, 1:3], names_from = SampleName, values_from = X)
ref <- as.matrix(ref)
rownames_ref <- ref[, 1]
ref <- ref[, -1]
ref <- apply(ref, 2, as.numeric)
rownames(ref) <- rownames_ref

size <- pivot_wider(tobe_genotyped[, c(1, 2, 5)], names_from = SampleName, values_from = R)
size <- as.matrix(size)
rownames_size <- size[, 1]
size <- size[, -1]
size <- apply(size, 2, as.numeric)
rownames(size) <- rownames_size

size[1:10, 1:10]
ref[1:10, 1:10]

## -----------------------------------------------------------------------------
library(updog)
multidog_obj <- multidog(
  refmat = ref,
  sizemat = size,
  model = "norm", # It depends of your population structure
  ploidy = 4, # Most common ploidy in your population
  nc = 6
) # Change the parameters accordingly!!

genos <- data.frame(
  MarkerName = multidog_obj$inddf$snp,
  SampleName = multidog_obj$inddf$ind,
  geno = multidog_obj$inddf$geno,
  prob = multidog_obj$inddf$maxpostprob
)

head(genos)

## ----eval=FALSE---------------------------------------------------------------
# library(fitPoly)
# saveMarkerModels(
#   ploidy = 4, # Most common ploidy among Texas Garden Roses collection
#   data = data_reference, # output of the previous function
#   filePrefix = "fitpoly_output/texas_roses", # Define the path and prefix for the result files
#   ncores = 2
# ) # Change it accordingly to your system!!!!
# 
# library(vroom) # Package to speed up reading files
# fitpoly_scores <- vroom("fitpoly_output/texas_roses_scores.dat")
# 
# genos <- data.frame(
#   MarkerName = fitpoly_scores$MarkerName,
#   SampleName = fitpoly_scores$SampleName,
#   geno = fitpoly_scores$maxgeno,
#   prob = fitpoly_scores$maxP
# )
# head(genos)

## ----eval=FALSE---------------------------------------------------------------
# # File containing markers genomic positions
# genos.pos_file <- read.table("geno.pos_roses_texas.txt", header = T)
# head(genos.pos)

## ----eval=FALSE---------------------------------------------------------------
# # Edit for Qploidy input format
# genos.pos <- data.frame(
#   MarkerName = genos.pos_file$probes,
#   Chromosome = genos.pos_file$chr,
#   Position = genos.pos_file$pos
# )
# head(genos.pos)

## ----echo=FALSE---------------------------------------------------------------
temp_file <- tempfile(fileext = ".tsv.gz") # saving in a temporary file
simu_data_standardized <- standardize(
  data = data,
  genos = genos,
  geno.pos = genos.pos,
  ploidy.standardization = 4,
  threshold.n.clusters = 5,
  n.cores = 1,
  out_filename = temp_file,
  verbose = TRUE
)

simu_data_standardized

## ----eval=FALSE---------------------------------------------------------------
# simu_data_standardized <- standardize(
#   data = data,
#   genos = genos,
#   geno.pos = genos.pos,
#   ploidy.standardization = 4,
#   threshold.n.clusters = 5,
#   n.cores = 1,
#   out_filename = "my_standardized_data.tsv.gz",
#   verbose = TRUE
# )
# 
# simu_data_standardized

## ----echo=FALSE---------------------------------------------------------------
simu_data_standardized <- read_qploidy_standardization(temp_file)

## ----eval=FALSE---------------------------------------------------------------
# simu_data_standardized <- read_qploidy_standardization("my_standardized_data.tsv.gz")

## -----------------------------------------------------------------------------
# Select a sample to be evaluated
samples <- unique(simu_data_standardized$data$SampleName)
head(samples)

## -----------------------------------------------------------------------------
sample <- "Tetraploid1"

# Proportion of heterozygous loci, BAF (Qploidy standardized ratio), and zscore by genomic positions
p <- plot_qploidy_standardization(
  x = simu_data_standardized,
  sample = sample,
  type = c("het", "BAF", "zscore"),
  dot.size = 0.05,
  chr = 1:2
)
p

## -----------------------------------------------------------------------------
# Heterozygous frequency, BAF (Qploidy standardized ratio), and zscore by genomic positions
# centromere positions added
p <- plot_qploidy_standardization(
  x = simu_data_standardized,
  sample = sample,
  type = c("het", "BAF", "zscore"),
  dot.size = 0.05,
  chr = 1:2,
  add_centromeres = TRUE,
  centromeres = c("chr1" = 15000000, "chr2" = 19000000)
)

p

## -----------------------------------------------------------------------------
# Raw allele intensity or read count ratio and BAF (Qploidy standardized ratio) histograms
# combining all markers in the sample (sample level resolution)
p <- plot_qploidy_standardization(
  x = simu_data_standardized,
  sample = sample,
  type = c("Ratio_hist_overall", "BAF_hist_overall"),
  chr = 1:2,
  ploidy = 4,
  add_expected_peaks = TRUE
)
p

## -----------------------------------------------------------------------------
# BAF histograms (chromosome level resolution) and Z score
p <- plot_qploidy_standardization(
  x = simu_data_standardized,
  sample = sample,
  type = c("BAF_hist", "zscore"),
  chr = 1:2,
  add_expected_peaks = TRUE,
  ploidy = c(4, 4)
)
p

## -----------------------------------------------------------------------------
# BAF histograms combining all markers in the sample (chromosome-arm level resolution) and Z score
p <- plot_qploidy_standardization(
  x = simu_data_standardized,
  sample = sample,
  type = c("BAF_hist", "zscore"),
  chr = 1:2,
  ploidy = c(4, 4, 4, 4), # Provide ploidy for each arm
  add_expected_peaks = TRUE,
  add_centromeres = TRUE,
  centromeres = c("chr1" = 15000000, "chr2" = 19000000)
)

p

## -----------------------------------------------------------------------------
# If sample level resolution
estimated_ploidies_sample <- area_estimate_ploidy(
  qploidy_standardization = simu_data_standardized, # standardization result object
  samples = "all", # Samples or "all" to estimate all samples
  level = "sample", # Resolution level
  ploidies = c(2, 3, 4, 5)
) # Ploidies to be tested
estimated_ploidies_sample

## -----------------------------------------------------------------------------
head(estimated_ploidies_sample$ploidy)

## -----------------------------------------------------------------------------
# If chromosome resolution
estimated_ploidies_chromosome <- area_estimate_ploidy(
  qploidy_standardization = simu_data_standardized,
  samples = "all",
  level = "chromosome",
  ploidies = c(2, 3, 4, 5)
)
estimated_ploidies_chromosome

## -----------------------------------------------------------------------------
head(estimated_ploidies_chromosome$ploidy)

## -----------------------------------------------------------------------------
# If chromosome-arm resolution
estimated_ploidies_chromosome_arm <- area_estimate_ploidy(
  qploidy_standardization = simu_data_standardized,
  samples = "all",
  level = "chromosome-arm",
  ploidies = c(2, 3, 4, 5),
  centromeres = c("chr1" = 15000000, "chr2" = 19000000)
)

estimated_ploidies_chromosome_arm

## -----------------------------------------------------------------------------
head(estimated_ploidies_chromosome_arm$ploidy)

## -----------------------------------------------------------------------------
estimated_ploidies_format <- merge_arms_format(estimated_ploidies_chromosome_arm)
estimated_ploidies_format

## -----------------------------------------------------------------------------
head(estimated_ploidies_format$ploidy)

## ----eval=FALSE---------------------------------------------------------------
# write.csv(estimated_ploidies_format$ploidy, file = "ploidies_before_visual_evaluation.csv")

## ----eval=FALSE---------------------------------------------------------------
# dir.create("figures")
# 
# for (i in seq_along(samples)) {
#   print(paste0("Part:", i, "/", length(samples)))
#   print(paste("Generating figures for", samples[i], "..."))
#   # This function will generate figures for sample, chromosome and chromosome-arm resolution level evaluation
#   all_resolutions_plots(
#     data_standardized = simu_data_standardized,
#     sample = samples[i],
#     ploidy = estimated_ploidies_chromosome$ploidy[i, 1:2],
#     chr = 1:2,
#     centromeres = c("chr1" = 15000000, "chr2" = 19000000),
#     file_name = paste0("figures/", samples[i])
#   )
#   print(paste("Done!"))
# }

## -----------------------------------------------------------------------------
ploidy_corrected <- estimated_ploidies_chromosome$ploidy

## ----eval=FALSE---------------------------------------------------------------
# # i) Select the highest resolution plots from the first standardization round.
# tetraploids <- apply(ploidy_corrected, 1, function(x) all(x == 4))
# tetraploids
# 
# data_reference2 <- rownames(ploidy_corrected)[which(tetraploids)]
# length(data_reference2)
# 
# # ii) Filter these samples from the `data` object.
# genos_round2 <- genos[which(genos$SampleName %in% data_reference2), ]
# 
# # iii) Run standardization again using this subset as the reference.
# simu_data_standardized_round2 <- standardize(
#   data = data,
#   genos = genos_round2,
#   geno.pos = genos.pos,
#   ploidy.standardization = 4,
#   threshold.n.clusters = 5,
#   n.cores = 1,
#   out_filename = "my_round2_standardization.tsv.gz",
#   verbose = TRUE
# )
# simu_data_standardized_round2

