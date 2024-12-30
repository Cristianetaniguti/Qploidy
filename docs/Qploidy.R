## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# #install.packages("devtools")
# devtools::install_github("cristianetaniguti/Qploidy")
# 
# library(Qploidy)

## ----eval=FALSE---------------------------------------------------------------
# data <- read_axiom("Qploidy_AxiomGT1_summary.txt")
# head(data)

## ----eval=FALSE---------------------------------------------------------------
# data <- read_illumina_array("data/illumina/Set1.txt", # This is a example code, data not available
#                             "data/illumina/Set2.txt",
#                             "data/illumina/Set3.txt")
# head(input_data)

## ----eval=FALSE---------------------------------------------------------------
# data <- qploidy_read_vcf("my_file.vcf") # This is a example code, data not available
# head(data)
# 

## ----eval=FALSE---------------------------------------------------------------
# tetraploids_subset <- read.csv("tetraploids.csv")
# data_reference <- data[which(data$SampleName %in% tetraploids_subset$x),]

## ----eval=FALSE---------------------------------------------------------------
# data_reference <- data

## ----eval=FALSE---------------------------------------------------------------
# library(fitPoly)
# saveMarkerModels(ploidy=4,                                # Most common ploidy among Texas Garden Roses collection
#                  data=data_reference,                     # output of the previous function
#                  filePrefix="fitpoly_output/texas_roses", # Define the path and prefix for the result files
#                  ncores=8)                                # Change it accordingly to your system!!!!
# 
# library(vroom) # Package to speed up reading files
# fitpoly_scores <- vroom("fitpoly_output/texas_roses_scores.dat")
# 
# genos <- data.frame(MarkerName = fitpoly_scores$MarkerName,
#                     SampleName = fitpoly_scores$SampleName,
#                     geno = fitpoly_scores$maxgeno,
#                     prob = fitpoly_scores$maxP)
# head(genos)

## ----eval=FALSE---------------------------------------------------------------
# genos.pos_file <- read.table("geno.pos_roses_texas.txt", header = T) # File containing markers genomic positions
# head(genos.pos)

## ----eval=FALSE---------------------------------------------------------------
# # Edit for Qploidy input format
# genos.pos <- data.frame(MarkerName = genos.pos_file$probes,
#                         Chromosome = genos.pos_file$chr,
#                         Position = genos.pos_file$pos)
# head(genos.pos)

## ----eval=FALSE---------------------------------------------------------------
# genos <- qploidy_read_vcf("my_file.vcf", geno = TRUE) # This is a example code, data not available
# dim(genos)
# 
# genos.pos <- qploidy_read_vcf("my_file.vcf", geno.pos = TRUE) # This is a example code, data not available
# head(genos.pos)

## ----eval=FALSE---------------------------------------------------------------
# genos.pos$MarkerName=paste0(genos.pos$Chromosome,"_",genos.pos$Position)
# head(genos.pos)

## ----eval=FALSE---------------------------------------------------------------
# genos <- genos[which(genos$SampleName %in% tetraploids_subset),]

## ----eval=FALSE---------------------------------------------------------------
# # Prepare inputs for updog
# ref <- pivot_wider(data[,1:3], names_from = SampleName, values_from = X) # This is a example code, data not available
# ref <- as.matrix(ref)
# rownames_ref <- ref[,1]
# ref<- ref[,-1]
# ref <- apply(ref, 2, as.numeric)
# rownames(ref) <- rownames_ref
# 
# size <- pivot_wider(data[,c(1,2,5)], names_from = SampleName, values_from = R)
# size <- as.matrix(size)
# rownames_size <- size[,1]
# size<- size[,-1]
# size <- apply(size, 2, as.numeric)
# rownames(size) <- rownames_size
# 
# library(updog)
# multidog_obj <- multidog(refmat = ref,
#                          sizemat = size,
#                          model = "norm", # It depends of your population structure
#                          ploidy = 6,     # Most common ploidy in your population
#                          nc = 6)         # Change the parameters accordingly!!
# 
# genos <- data.frame(MarkerName = multidog_obj$inddf$snp,
#                     SampleName = multidog_obj$inddf$ind,
#                     geno = multidog_obj$inddf$geno,
#                     prob = multidog_obj$inddf$maxpostprob)
# 
# head(genos)

## ----eval=FALSE---------------------------------------------------------------
# pos <- strsplit(multidog_obj$snpdf$snp, "_") # This is a example code, data not available
# 
# genos.pos <- data.frame(MarkerName = multidog_obj$snpdf$snp,
#                         Chromosome = sapply(pos, "[[", 1),
#                         Position = as.numeric(sapply(pos, "[[", 2)))
# head(geno.pos)

## ----eval=FALSE---------------------------------------------------------------
# roses_data_standardized <- standardize(data = data,
#                                        genos = genos,
#                                        geno.pos = genos.pos,
#                                        ploidy.standardization = 4,
#                                        threshold.n.clusters = 5,
#                                        n.cores =8,
#                                        out_filename = "standardization_results/roses_texas_5clust.tsv.gz",
#                                        type = "intensities",
#                                        verbose = TRUE)
# 
# roses_data_standardized

## ----eval=FALSE---------------------------------------------------------------
# roses_data_standardized <- read_qploidy_standardization("standardization_results/roses_texas_5clust.tsv.gz")

## ----eval=FALSE---------------------------------------------------------------
# # Select a sample to be evaluated
# samples <- unique(roses_data_standardized$data$SampleName)
# head(samples)

## ----eval=FALSE---------------------------------------------------------------
# sample <- "262_97_4"
# 
# # Proportion of heterozygous loci, BAF (Qploidy standardized ratio), and zscore by genomic positions
# p <- plot_qploidy_standardization(x = roses_data_standardized,
#                                   sample = sample,
#                                   type = c("het", "BAF", "zscore"),
#                                   dot.size = 0.05,
#                                   chr = 2:8)
# ggsave(p, filename = "fig1.png")

## ----eval=FALSE---------------------------------------------------------------
# # Heterozygous frequency, BAF (Qploidy standardized ratio), and zscore by genomic positions  - centromere positions added
# p <- plot_qploidy_standardization(x = roses_data_standardized,
#                                   sample = sample,
#                                   type = c("het", "BAF", "zscore"),
#                                   dot.size = 0.05,
#                                   chr = 2:8,
#                                   add_centromeres = TRUE,
#                                   centromeres = c("1" = 22000000, "2" = 36000000, "3" = 4000000,
#                                                   "4" = 20000000, "5" = 52000000, "6" = 32000000, "7" = 20000000))
# 
# ggsave(p, filename = "fig2.png")
# 

## ----eval=FALSE---------------------------------------------------------------
# # Raw allele intensity or read count ratio and BAF (Qploidy standardized ratio) histograms combining all markers in the sample (sample level resolution)
# p <- plot_qploidy_standardization(x = roses_data_standardized,
#                                   sample = sample,
#                                   type = c("Ratio_hist_overall", "BAF_hist_overall"),
#                                   chr = 2:8,
#                                   ploidy = c(4,4,5,4,4,4,4),
#                                   add_expected_peaks = TRUE)
# 
# ggsave(p, filename = "fig3.png")

## ----eval=FALSE---------------------------------------------------------------
# # BAF histograms (chromosome level resolution) and Z score
# p <- plot_qploidy_standardization(x = roses_data_standardized,
#                                   sample = sample,
#                                   type = c("BAF_hist", "zscore"),
#                                   chr = 2:8,
#                                   add_expected_peaks = TRUE,
#                                   ploidy = c(4,4,5,4,4,4,4))
# ggsave(p, filename = "fig4.png")
# 

## ----eval=FALSE---------------------------------------------------------------
# # BAF histograms combining all markers in the sample (chromosome-arm level resolution) and Z score
# p <- plot_qploidy_standardization(x = roses_data_standardized,
#                                   sample = sample,
#                                   type = c("BAF_hist", "zscore"),
#                                   chr = 2:8,
#                                   ploidy = rep(c(4,4,5,4,4,4,4), each = 2), # Provide ploidy for each arm
#                                   add_expected_peaks = TRUE,
#                                   add_centromeres = TRUE,
#                                   centromeres = c("1" = 22000000, "2" = 36000000, "3" = 4000000,
#                                                   "4" = 20000000, "5" = 52000000, "6" = 32000000, "7" = 20000000))
# 
# ggsave(p, filename = "fig5.png")

## ----eval=FALSE---------------------------------------------------------------
# # If sample level resolution
# estimated_ploidies_sample <- area_estimate_ploidy(qploidy_standardization = roses_data_standardized, # standardization result object
#                                                   samples = "all",                                  # Samples or "all" to estimate all samples
#                                                   level = "sample",                                 # Resolution level
#                                                   ploidies = c(2,5))                                # Ploidy range to investigate
# estimated_ploidies_sample

## ----eval=FALSE---------------------------------------------------------------
# head(estimated_ploidies_sample$ploidy)

## ----eval=FALSE---------------------------------------------------------------
# # If chromosome resolution
# estimated_ploidies_chromosome <- area_estimate_ploidy(qploidy_standardization = roses_data_standardized,
#                                                       samples = "all",
#                                                       level = "chromosome",
#                                                       ploidies = c(2,5))
# estimated_ploidies_chromosome

## ----eval=FALSE---------------------------------------------------------------
# head(estimated_ploidies_chromosome$ploidy)

## ----eval=FALSE---------------------------------------------------------------
# # If chromosome-arm resolution
# estimated_ploidies_chromosome_arm <- area_estimate_ploidy(qploidy_standardization = roses_data_standardized,
#                                                           samples = "all",
#                                                           level = "chromosome-arm",
#                                                           ploidies = c(2,5),
#                                                           centromeres = c("1" = 22000000, "2" = 36000000, "3" = 4000000,
#                                                                           "4" = 20000000, "5" = 52000000, "6" = 32000000,
#                                                                           "7" = 20000000))
# estimated_ploidies_chromosome_arm

## ----eval=FALSE---------------------------------------------------------------
# head(estimated_ploidies_chromosome_arm$ploidy)

## ----eval=FALSE---------------------------------------------------------------
# estimated_ploidies_format <- merge_arms_format(estimated_ploidies_chromosome_arm)
# estimated_ploidies_format

## ----eval=FALSE---------------------------------------------------------------
# head(estimated_ploidies_format$ploidy)

## ----eval=FALSE---------------------------------------------------------------
# write.csv(estimated_ploidies_format$ploidy, file = "ploidies_before_visual_evaluation.csv")

## ----eval=FALSE---------------------------------------------------------------
# dir.create("figures")
# 
# for(i in 1:length(samples)){
#   print(paste0("Part:",i,"/",length(samples)))
#   print(paste("Generating figures for", samples[i],"..."))
#   # This function will generate figures for sample, chromosome and chromosome-arm resolution level evaluation
#   all_resolutions_plots(data_standardized = roses_data_standardized,
#                         sample = samples[i],
#                         ploidy = estimated_ploidies_chromosome$ploidy[i,2:8],
#                         chr = 2:8,
#                         centromeres = c("1" = 22000000, "2" = 36000000, "3" = 4000000,
#                                         "4" = 20000000, "5" = 52000000, "6" = 32000000, "7" = 20000000),
#                         file_name = paste0("figures/",samples[i]))
#   print(paste("Done!"))
# }

## ----eval=FALSE---------------------------------------------------------------
# ploidy_corrected <- read.csv("ploidies_after_visual_inspection.csv")
# 

## ----eval=FALSE---------------------------------------------------------------
# 
# #i) Select the highest resolution plots from the first standardization round.
# data_reference2 <- ploidy_corrected$ID_array[which(ploidy_corrected$ploidy == 4 & ploidy_corrected$Resolution == "chromosome-arm")]
# length(data_reference2) # Total of 42 individuals, more is better
# 
# data_reference2 <- ploidy_corrected$ID_array[which(ploidy_corrected$ploidy == 4 & ploidy_corrected$Resolution == "chromosome")]
# length(data_reference2) # With chromosome resolution we found 219, a better sample size
# 
# #ii) Filter these samples from the `data` object.
# genos_round2 <- genos[which(genos$SampleName %in% data_reference2),]
# 
# #iii) Run standardization again using this subset as the reference.
# roses_data_standardized_round2 <- standardize(data = data,
#                                               genos = genos_round2,
#                                               geno.pos = genos.pos,
#                                               ploidy.standardization = 4,
#                                               threshold.n.clusters = 5,
#                                               n.cores = 8,
#                                               out_filename = "standardization_results/roses_texas_5clust_2round.tsv.gz",
#                                               type = "intensities",
#                                               verbose = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# #iv) Reevaluate the plots and estimate the ploidy once more.
# 
# sample <- "262_97_4"
# # Proportion of heterozygous loci, BAF (Qploidy standardized ratio), and zscore by genomic positions
# p <- plot_qploidy_standardization(x = roses_data_standardized_round2,
#                                   sample = sample,
#                                   type = c("het", "BAF", "zscore"),
#                                   dot.size = 0.05,
#                                   chr = 2:8,
#                                   ploidy = c(4,4,5,4,4,4,4))
# 
# ggsave(p, filename = "fig7.png")
# 

## ----eval=FALSE---------------------------------------------------------------
# # BAF histograms combining all markers in the sample (chromosome-arm level resolution) and Z score
# p <- plot_qploidy_standardization(x = roses_data_standardized_round2,
#                                   sample = sample,
#                                   type = c("BAF_hist", "zscore"),
#                                   chr = 2:8,
#                                   ploidy = c(4,4,5,4,4,4,4),
#                                   add_expected_peaks = TRUE)
# 
# ggsave(p, filename = "fig6.png")

## ----eval=FALSE---------------------------------------------------------------
# # If chromosome-arm resolution
# estimated_ploidies_chromosome_arm <- area_estimate_ploidy(qploidy_standardization = roses_data_standardized_round2,
#                                                           samples = "all",
#                                                           level = "chromosome-arm",
#                                                           ploidies = c(2,5),
#                                                           centromeres = c("1" = 22000000, "2" = 36000000, "3" = 4000000,
#                                                                           "4" = 20000000, "5" = 52000000, "6" = 32000000,
#                                                                           "7" = 20000000))
# 
# estimated_ploidies_format <- merge_arms_format(estimated_ploidies_chromosome_arm)
# 
# 
# write.csv(estimated_ploidies_format$ploidy, file = "ploidies_before_visual_evaluation_round2.csv")
# 
# for(i in 1:length(samples)){
#   print(paste0("Part:",i,"/",length(samples)))
#   print(paste("Generating figures for", samples[i],"..."))
#   # This function will generate figures for sample, chromosome and chromosome-arm resolution level evaluation
#   all_resolutions_plots(data_standardized = roses_data_standardized_round2,
#                         sample = samples[i],
#                         ploidy = estimated_ploidies_chromosome$ploidy[i,2:8],
#                         chr = 2:8,
#                         centromeres = c("1" = 22000000, "2" = 36000000, "3" = 4000000,
#                                         "4" = 20000000, "5" = 52000000, "6" = 32000000, "7" = 20000000),
#                         file_name = paste0("figures/",samples[i],"_round2"))
#   print(paste("Done!"))
# }

