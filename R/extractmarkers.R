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
#' @import stringr
#' @import tools
#' @examples
#' extractmarkers(input.file = "test.vcf", snps.list = "markers.txt")
#' @examples extractmarkers(bed.file = "data.bed", bim.file = "data.bim", fam.file = "data.fam", pos.list = "markers.txt")
#' @export
extract_markers <- function(input.file, 
                            snps.list = NULL, 
                            pos.list = NULL, 
                            bed.file = NULL, 
                            bim.file = NULL, 
                            fam.file = NULL, 
                            output.dir = "marker_outputs", 
                            merged.file = "final_merged.vcf") {
  # Step 1: Detect file type
  file_type <- detect_file_type(input.file, bed.file, bim.file, fam.file)
  
  # Step 2: Generate extraction commands & execute PLINK
  commands <- construct_and_execute_plink_command(file_type, 
                                                  input.file, 
                                                  snps.list, 
                                                  pos.list, 
                                                  bed.file, 
                                                  bim.file, 
                                                  fam.file, 
                                                  output.dir, 
                                                  merged.file)
  
  # Step 3: Print completion message
  print("Marker extraction and merging completed successfully!")
}


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

get_plink_path <- function() {
  # Default to 'plink' assuming it's in the Linux container's PATH
  plink_path <- Sys.which("plink")
  
  if (plink_path == "") {
    stop("PLINK is not available in the PATH. Ensure PLINK is installed and accessible.")
  }
  
  return(plink_path)
}

create_output_directory <- function(directory) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
    print(paste("Created directory:", directory))
  } else {
    print(paste("Directory already exists:", directory))
  }
}

merge_vcf_files <- function(output.dir, merged.file) {
  vcf_files <- list.files(output.dir, pattern = "*.vcf$", full.names = TRUE)
  
  if (length(vcf_files) == 0) {
    stop("No extracted VCF files found for merging.")
  }
  
  merge_command <- stringr::str_c("bcftools concat -o ", merged.file, " ", paste(vcf_files, collapse = " "))
  
  system(merge_command)
  print(paste("Merged VCF file created:", merged.file))
}

construct_and_execute_plink_command <- function(file_type, 
                                                input.file, 
                                                snps.list = NULL, 
                                                pos.list = NULL, 
                                                bed.file = NULL, 
                                                bim.file = NULL, 
                                                fam.file = NULL, 
                                                output.dir, 
                                                merged.file) {
  create_output_directory(output.dir)  # Ensure directory exists
  plink_path <- get_plink_path()  # Get correct PLINK path
  
  if (!is.null(snps.list)) {
    # Extract markers using rsID
    command <- switch(file_type,
                      "VCF" = stringr::str_c(plink_path, 
                                             " --vcf ", 
                                             input.file, 
                                             " --const-fid 0 --cow --extract ", 
                                             snps.list, 
                                             " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                             file.path(output.dir, "rsid_extracted")),
                      "VCF_GZ" = stringr::str_c(plink_path, 
                                            " --vcf ", input.file, 
                                            " --const-fid 0 --cow --extract ", snps.list, 
                                            " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                            file.path(output.dir, "rsid_extracted")),
                      "BCF" = stringr::str_c(plink_path, 
                                             " --bcf ", 
                                             input.file, 
                                             " --const-fid 0 --cow --extract ", 
                                             snps.list, 
                                             " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                             file.path(output.dir, "rsid_extracted")),
                      "PLINK" = stringr::str_c(plink_path, 
                                               " --bfile ", 
                                               sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file), 
                                               " --const-fid 0 --cow --extract ", 
                                               snps.list, 
                                               " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                               file.path(output.dir, "rsid_extracted"))
    )
    system(command)  # Execute PLINK
  } else if (!is.null(pos.list)) {
    # Extract markers using positions (creates multiple output files)
    command_list <- lapply(seq_len(nrow(pos.list)), function(i) {
      output_file <- file.path(output.dir, paste0("marker_", i, ".vcf"))
      
      cmd <- stringr::str_c(plink_path, " --", tolower(file_type), " ", input.file,
                            " --cow --chr ", pos.list[i, 1],
                            " --from-bp ", pos.list[i, 2],
                            " --to-bp ", pos.list[i, 3],
                            " --recode vcf --keep-allele-order --out ", output_file
      )
      system(cmd)  # Execute PLINK immediately
      return(cmd)  # Store for reference
    })
    
    merge_vcf_files(output.dir, merged.file)  # Merge extracted files
    
    return(command_list)
  } else {
    stop("Either `snps.list` or `pos.list` must be provided.")
  }
  
  return(command)
}



