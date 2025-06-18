library(testthat)
library(stringr)
library(tools)

test_that("extract.markers runs without errors and correctly detects file type", {
   expect_equal(detect_file_type("test.vcf"), "VCF")
   expect_equal(detect_file_type("test.bcf"), "BCF")
   expect_equal(detect_file_type("test.bed", "test.bed", "test.bim", "test.fam"), "PLINK")
   
   # Ensure function stops on unsupported file type
   expect_error(detect_file_type("test.txt"))
})

test_that("create_output_directory creates the expected directory", {
   test_dir <- "test_output"
   create_output_directory(test_dir)
   expect_true(dir.exists(test_dir))
   
   # Cleanup
   unlink(test_dir, recursive = TRUE)
})

test_that("construct_plink_command generates correct rsID extraction commands", {
   cmd_vcf <- construct_plink_command("VCF", "test.vcf", snps.list = "test_snps.txt", output.dir = "output")
   expect_match(cmd_vcf, "--vcf test.vcf")
   
   cmd_bcf <- construct_plink_command("BCF", "test.bcf", snps.list = "test_snps.txt", output.dir = "output")
   expect_match(cmd_bcf, "--bcf test.bcf")
   
   cmd_plink <- construct_plink_command("PLINK", "test.bed", snps.list = "test_snps.txt", bed.file = "test.bed", bim.file = "test.bim", fam.file = "test.fam", output.dir = "output")
   expect_match(cmd_plink, "--bfile --bed test.bed --bim test.bim --fam test.fam")
})

test_that("construct_plink_command generates correct position-based extraction commands", {
   pos_list <- data.frame(Chr = c(1, 2), Start = c(1000, 2000), End = c(5000, 6000), Output = c("marker1", "marker2"))
   
   commands <- construct_plink_command("VCF", "test.vcf", pos_list = pos_list, output.dir = "output")
   expect_type(commands, "list")
   expect_true(all(grepl("--vcf test.vcf", commands)))
})

test_that("execute_plink_commands runs without error", {
   expect_silent(execute_plink_commands("echo Test Command"))
})

test_that("merge_vcf_files correctly identifies VCF files for merging", {
   test_dir <- "test_merge_output"
   dir.create(test_dir)
   
   # Create dummy VCF files
   file.create(file.path(test_dir, "test1.vcf"))
   file.create(file.path(test_dir, "test2.vcf"))
   
   # Test merging function
   expect_silent(merge_vcf_files(test_dir, "merged.vcf"))
   
   # Cleanup
   unlink(test_dir, recursive = TRUE)
})