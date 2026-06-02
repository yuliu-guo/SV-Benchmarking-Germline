#!/usr/bin/env Rscript

#
#library(VariantAnnotation)
#library(stringr)
#library(GenomicRanges)
#
## Args: input folder
#args <- commandArgs(trailingOnly = TRUE)
#if (length(args) != 1) {
#  stop("Usage: Rscript annotate_vcfs.R /path/to/vcfs")
#}
#input_dir <- args[1]
#
## Simple SV type classifier based on available INFO fields
#simpleEventType <- function(vcf) {
#  info_df <- info(vcf)
#
#  # Try to extract INSLEN and SVLEN; default to NA if missing
#  inslen <- suppressWarnings(as.numeric(info_df$INSLEN))
#  svlen <- suppressWarnings(as.numeric(info_df$SVLEN))
#  svtype <- as.character(info_df$SVTYPE)
#
#  n <- length(svtype)
#  result <- rep("BND", n)
#
#  for (i in seq_len(n)) {
#    if (!is.na(svtype[i])) {
#      if (svtype[i] == "DEL") {
#        result[i] <- "DEL"
#      } else if (svtype[i] == "DUP") {
#        result[i] <- "DUP"
#      } else if (svtype[i] == "INV") {
#        result[i] <- "INV"
#      } else if (!is.na(inslen[i]) && !is.na(svlen[i]) && abs(inslen[i]) > abs(svlen[i]) * 0.7) {
#        result[i] <- "INS"
#      } else {
#        result[i] <- svtype[i]
#      }
#    }
#  }
#
#  return(result)
#}
#
#
## Process all .vcf.gz files in the folder
#vcf_files <- list.files(input_dir, pattern = "\\.vcf$", full.names = TRUE)
#
#for (vcf_file in vcf_files) {
#  cat("Processing:", vcf_file, "\n")
#
#  vcf <- readVcf(vcf_file, "hg38")
#
#  # Add SIMPLE_TYPE info field to header
#  hdr <- header(vcf)
#  new_info <- DataFrame(
#    Number = "1",
#    Type = "String",
#    Description = "Simple event type annotation based on SVTYPE/SVLEN/INSLEN",
#    row.names = "SIMPLE_TYPE"
#  )
#  info(hdr) <- rbind(info(hdr), new_info)
#  header(vcf) <- hdr
#
#  # Generate annotations
#  simple_type <- simpleEventType(vcf)
#  info(vcf)$SIMPLE_TYPE <- simple_type
#
#  # Output filename
#  out_base <- sub("\\.vcf$", "", basename(vcf_file))
#  output_file <- file.path(input_dir, paste0(out_base, "_annotated.vcf"))
#  writeVcf(vcf, output_file)
#
#  cat("Written annotated VCF to:", output_file, "\n")
#}
#

#!/usr/bin/env Rscript

suppressMessages({
  library(VariantAnnotation)
})

# Define SV type classifier
simpleEventType <- function(gr) {
  pgr <- partner(gr)
  ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomosal
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS",
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL", "DUP")
      )
    )
  )
}

# Read input directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript batch_svtype_annotation.R /path/to/vcf_folder")
}
input_dir <- args[1]

vcf_files <- list.files(input_dir, pattern = "\\.vcf$", full.names = TRUE)
if (length(vcf_files) == 0) {
  stop("No VCF files found.")
}

cat("Found", length(vcf_files), "VCF file(s)\n")

for (vcf_file in vcf_files) {
  cat("Processing:", basename(vcf_file), "\n")
  vcf <- readVcf(vcf_file, "hg38")

  # Add missing headers for new INFO fields
  hdr <- header(vcf)

  # Add SIMPLE_TYPE
  if (!"SIMPLE_TYPE" %in% rownames(info(hdr))) {
    new_simple <- DataFrame(
      Number = "1",
      Type = "String",
      Description = "Simple event type based on breakend pairs",
      row.names = "SIMPLE_TYPE"
    )
    info(hdr) <- rbind(info(hdr), new_simple)
  }

  # Add SVTYPE
  if (!"SVTYPE" %in% rownames(info(hdr))) {
    new_svtype <- DataFrame(
      Number = "1",
      Type = "String",
      Description = "Inferred structural variant type",
      row.names = "SVTYPE"
    )
    info(hdr) <- rbind(info(hdr), new_svtype)
  }

  # Add SVLEN
  if (!"SVLEN" %in% rownames(info(hdr))) {
    new_svlen <- DataFrame(
      Number = "1",
      Type = "Integer",
      Description = "Inferred structural variant length",
      row.names = "SVLEN"
    )
    info(hdr) <- rbind(info(hdr), new_svlen)
  }

  header(vcf) <- hdr

  # Annotate breakpoints
  gr <- breakpointRanges(vcf)
  gr <- gr[gr$sourceId %in% names(vcf)]  # ensure valid mapping
  svtype <- simpleEventType(gr)

  # Assign annotations to matching VCF rows
  info(vcf)$SIMPLE_TYPE <- NA_character_
  info(vcf)$SVTYPE <- NA_character_
  info(vcf)$SVLEN <- NA_integer_

  idx <- gr$sourceId
  info(vcf)[idx, "SIMPLE_TYPE"] <- svtype
  info(vcf)[idx, "SVTYPE"] <- svtype
  info(vcf)[idx, "SVLEN"] <- gr$svLen

  # Output filename
  out_base <- sub("\\.vcf(\\.gz)?$", "", basename(vcf_file))
  output_file <- file.path(input_dir, paste0(out_base, "_annotated.vcf"))
  writeVcf(vcf, output_file)

  cat("Annotated VCF saved to:", output_file, "\n\n")
}

cat("All files processed successfully")
