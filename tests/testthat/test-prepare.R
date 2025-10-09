library(testthat)

# Test for qploidy_read_vcf
# Assuming a sample VCF file exists for testing
# Replace 'sample.vcf' with the actual path to a test VCF file

test_that("qploidy_read_vcf handles VCF input correctly", {
  temp <- tempfile(fileext = ".vcf")
  simulate_vcf(temp, seed = 1213)

  checks <- vcf_sanity_check(temp)
  expect_true(all(checks$checks[-c(12,13,14,15)]))

  result <- qploidy_read_vcf(temp, geno = FALSE, geno.pos = FALSE)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("MarkerName", "SampleName", "X", "Y", "R", "ratio") %in% colnames(result)))

  result_geno <- qploidy_read_vcf(temp, geno = TRUE, geno.pos = FALSE)
  expect_s3_class(result_geno, "data.frame")
  expect_true(all(c("MarkerName", "SampleName", "geno", "prob") %in% colnames(result_geno)))

  result_geno_pos <- qploidy_read_vcf(temp, geno = FALSE, geno.pos = TRUE)
  expect_s3_class(result_geno_pos, "data.frame")
  expect_true(all(c("MarkerName", "Chromosome", "Position") %in% colnames(result_geno_pos)))
})

# Test for read_illumina_array
# Assuming sample Illumina array files exist for testing
# Replace 'sample1.txt' and 'sample2.txt' with actual paths to test files

test_that("read_illumina_array processes Illumina array files correctly", {
  temp1 <- tempfile(fileext = ".txt")
  temp2 <- tempfile(fileext = ".txt")
  simulate_illumina_file(file = temp1, seed = 1213, num_snps = 100, num_samples = 50, mk_id = "MK1-")
  simulate_illumina_file(temp2, seed = 2312, num_snps = 100, num_samples = 30, mk_id = "MK2-")

  result <- read_illumina_array(temp1, temp2)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("MarkerName", "SampleName", "X", "Y", "R", "ratio") %in% colnames(result)))
})

# Test for read_axiom
# Assuming a sample Axiom summary file exists for testing
# Replace 'summary.txt' with the actual path to a test summary file

test_that("read_axiom processes Axiom summary files correctly", {
  temp <- tempfile(fileext = ".txt")
  simulate_axiom_summary(temp, seed = 1213, n_probes = 100, n_samples = 50)

  # Create the data frame
  ind_names_example <- data.frame(
    Plate_Name = paste0("Sample", 1:50),
    Sample_Name = paste0("NewName", 1:50)
  )

  temp_name <- tempfile(fileext = ".txt")
  # Save the data frame to a file for testing
  write.table(
    ind_names_example,
    file = temp_name,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  result <- read_axiom(summary_file = temp, atan = FALSE, ind_names = temp_name)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("MarkerName", "SampleName", "X", "Y", "R", "ratio") %in% colnames(result)))
})
