#' Extracts a set of markers from PLINK, VCF, or BCF files
#' @description
#' The 'extract_markers' function extracts a set of markers, either using their rsID or chromosome position, from PLINK, VCF, or BCF files. 
#' In using rsIDs for extraction, a text file containing one marker per line is accepted. When using GRCh38/GRCh37 positions, the text file should contain 4 columns:
#' the first column should contain the chromosome number only (e.g. "1" and not "chr1"); the second column should contain the starting bp position; 
#' the third column the final bp position; and the fourth column should specify the desired output name with no file extension (e.g. "snpsfromchr1").
#' @param input.file (character) The file path to the genetic dataset. Accepts vcf or bcf files.
#' @param snps.list (character) The file path to the text file containing the list of markers (rsID) to be extracted. Each line will be read as one marker. Should be a text file.
#' @param pos.list (character) The file path to the text file containing the list of markers to be extracted. If using rsIDs, each line should correspond to one marker. If using POS (positions), there should be 4 columns: [1] chromosome numner (integer), [2] starting base position, [3] final base position, and [4] output file name. 
#' @param bed.file (character) The file path for the bed file. Should be provided if plink.files = TRUE
#' @param bim.file (character) The file path for the bim file. Should be provided if plink.files = TRUE
#' @param fam.file (character) The file path for the fam file. Should be provided if plink.files = TRUE
#' @param plink_args (character) The set of plink-specific arguments and their values.
#' @import stringr
#' @import tools
#' @examples
#' extractmarkers(input.file = "test.vcf", snps.list = "markers.txt", args = "--maf 0.1 --geno 0.2")
#' @examples extractmarkers(bed.file = "data.bed", bim.file = "data.bim", fam.file = "data.fam", pos.list = "markers.txt")
#' @export
extract_markers <- function(input.file, 
                            snps.list = NULL, 
                            pos.list = NULL, 
                            bed.file = NULL, 
                            bim.file = NULL, 
                            fam.file = NULL, 
                            plink_args = NULL,
                            output.dir = output.dir, 
                            merged.file = "final_merged.vcf",
                            plink_path = plink_path) {
  
  file_type <- detect_file_type(input.file, bed.file, bim.file, fam.file)
  
  commands <- extraction(file_type, 
                         input.file, 
                         snps.list, 
                         pos.list, 
                         bed.file, 
                         bim.file, 
                         fam.file, 
                         plink_args,
                         output.dir, 
                         merged.file,
                         plink_path)
  
  return(list(
    extracted_files = list.files(output.dir, full.names = FALSE)
  ))
}

# TO DO #1: ADDITIONAL ARGUMENTS FOR PLINK - OPTIONAL FILTERING PROCEDURES
detect_file_type <- function(input.file, bed.file = NULL, bim.file = NULL, fam.file = NULL) {
  file_extension <- tools::file_ext(input.file)
  
  if (!is.null(bed.file) && !is.null(bim.file) && !is.null(fam.file)) {
    return("PLINK")
  } else if (grepl("\\.vcf\\.gz$", input.file)) {
    return("VCF_GZ")
  } else if (file_extension == "vcf") {
    return("VCF")
  } else if (file_extension == "bcf") {
    return("BCF")
  } else {
    stop("Unsupported file type. Please provide a VCF, VCF.GZ, BCF, or PLINK file.")
  }
}

merge_vcf_files <- function(output.dir, merged.file) {
  vcf_files <- list.files(output.dir, pattern = "*.vcf$", full.names = TRUE)
  
  if (length(vcf_files) == 0) {
    stop("No extracted VCF files found for merging.")
  }
  
  # July 30 note: replaced str_c with system2 to source out bcftools
  merge_command <- system2("bcftools",
                           args = c("concat", "-o", merged.file),
                           stdout = paste(vcf_files, collapse = " "))
  
  system(merge_command)
  print(paste("Merged VCF file created:", merged.file))
}

extraction <- function(file_type, 
                       input.file, 
                       snps.list = NULL, 
                       pos.list = NULL, 
                       bed.file = NULL, 
                       bim.file = NULL, 
                       fam.file = NULL, 
                       plink_args = NULL,
                       output.dir, 
                       merged.file,
                       plink_path) {
  
  if (!dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  if (!is.null(snps.list)) {
    args <- if (!is.null(plink_args)) paste(plink_args, collapse = " ") else ""
    
    # using rsID
    # added args 4 Aug 2025
    command <- switch(file_type,
                      "VCF" = stringr::str_c(plink_path, 
                                             " --vcf ", 
                                             input.file, 
                                             " --const-fid 0 --cow --extract ", 
                                             snps.list, 
                                             " ", args,
                                             " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                             file.path(output.dir, "rsid_extracted")),
                      "VCF_GZ" = stringr::str_c(plink_path, 
                                                " --vcf ", input.file, 
                                                " --const-fid 0 --cow --extract ", snps.list, 
                                                " ", args,
                                                " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                                file.path(output.dir, "rsid_extracted")),
                      "BCF" = stringr::str_c(plink_path, 
                                             " --bcf ", 
                                             input.file, 
                                             " --const-fid 0 --cow --extract ", 
                                             snps.list, 
                                             " ", args,
                                             " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                             file.path(output.dir, "rsid_extracted")),
                      "PLINK" = stringr::str_c(plink_path, 
                                               " --bfile ", 
                                               sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file), 
                                               " --const-fid 0 --cow --extract ", 
                                               snps.list, 
                                               " ", args,
                                               " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                               file.path(output.dir, "rsid_extracted"))
    )
    system(command)  # Execute PLINK
  } else if (!is.null(pos.list)) {
    args <- if (!is.null(plink_args)) paste(plink_args, collapse = " ") else ""
    
    # Extract markers using positions (creates multiple output files)
    command_list <- lapply(seq_len(nrow(pos.list)), function(i) {
      chr_num <- pos.list[i, 1]
      start_bp <- pos.list[i, 2]
      end_bp <- pos.list[i, 3]
      
      # Subtract 1 from start position for filename
      filename_base <- paste0("chr", chr_num, "_", start_bp - 1)
      output_file <- file.path(output.dir, paste0(filename_base, ".vcf"))
      
      cmd <- stringr::str_c(
        plink_path, " --", tolower(file_type), " ", input.file,
        " --cow --chr ", chr_num,
        " --from-bp ", start_bp,
        " --to-bp ", end_bp,
        " ", args,
        " --recode vcf --keep-allele-order --out ", tools::file_path_sans_ext(output_file)
      )
      
      system(cmd)  # Run immediately
      return(cmd)
    })
    
    merge_vcf_files(output.dir, merged.file)  # Merge extracted files
    
    return(command_list)
  } else {
    stop("Either `snps.list` or `pos.list` must be provided.")
  }
  
  return(command)
}
