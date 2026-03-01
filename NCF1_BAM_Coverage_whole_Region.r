# ===============================
# SNP Coverage vs. Position Plot from BAMs - Right and Working Script for line graph showing coverage across the whole NCF1 region between different Thresholds
# ===============================

# Load required packages
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(ggplot2)
library(dplyr)

# 1. Load BED file with SNP positions ------------------------
bed_file <- "/home/kiliankapitza/Desktop/NCF1_KK/NCF1_start_to_end.bed"

bed <- read.table(bed_file, sep = "\t", header = FALSE,
                  col.names = c("Chromosome", "Start", "End", "Name"),
                  stringsAsFactors = FALSE)


# 2. List BAM files ------------------------------------------
bam_dir <- "/home/kiliankapitza/Desktop/NCF1_KK/Ergebnisse_Stand_031225/NCF1_20250520_AS_NP_TP5/19-mer/Alignment"
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

# Check BAM index files
missing_index <- bam_files[!file.exists(paste0(bam_files, ".bai"))]
if (length(missing_index) > 0) {
  stop("❌ Missing BAM index (.bai) for:\n", paste(missing_index, collapse = "\n"))
} else {
  message("✅ All BAM index files found.")
}

# 3. Adjust BED chromosome names to match BAM ----------------
bam_chr <- names(scanBamHeader(bam_files[1])[[1]]$targets)

if (!all(bed$Chromosome %in% bam_chr)) {
  if (all(paste0("chr", bed$Chromosome) %in% bam_chr)) {
    bed$Chromosome <- paste0("chr", bed$Chromosome)
    message("✅ Added 'chr' prefix to BED chromosomes.")
  } else if (all(sub("^chr","", bed$Chromosome) %in% bam_chr)) {
    bed$Chromosome <- sub("^chr","", bed$Chromosome)
    message("✅ Removed 'chr' prefix from BED chromosomes.")
  } else {
    stop("❌ Chromosome names in BED and BAM do not match.")
  }
} else {
  message("✅ Chromosomes already match BAM.")
}

# Create GRanges object for SNP positions
regions <- GRanges(seqnames = bed$Chromosome,
                   ranges = IRanges(start = bed$Start, end = bed$End),
                   names = bed$Name)


# 4. Function to get coverage at every position ----------------
get_coverage_in_region <- function(bam_file, regions) {
  message("Processing: ", bam_file)
  
  # Read alignments only from the defined regions
  aln <- readGAlignments(bam_file, param = ScanBamParam(which = regions))
  
  # If no alignments found at all
  if (length(aln) == 0) {
    warning("⚠️ No reads found in defined regions for file: ", bam_file)
    return(NULL)
  }
  
  cov <- coverage(aln)
  result_list <- list()
  
  for (i in seq_along(regions)) {
    chr <- as.character(seqnames(regions[i]))
    start_pos <- start(regions[i])
    end_pos <- end(regions[i])
    
    if (!chr %in% names(cov)) {
      warning("⚠️ Chromosome ", chr, " not found in BAM coverage: ", bam_file)
      next
    }
    
    cov_chr <- cov[[chr]]
    chr_len <- length(cov_chr)
    
    # Limit to actual coverage range
    start_pos <- max(1, start_pos)
    end_pos <- min(end_pos, chr_len)
    
    region_length <- end_pos - start_pos + 1
    if (region_length <= 0) {
      warning("⚠️ Invalid or empty region for ", chr, " in ", bam_file)
      next
    }
    
    # Extract coverage safely
    cov_values <- as.numeric(cov_chr[start_pos:end_pos])
    
    # If coverage vector is empty, fill with zeros
    if (length(cov_values) == 0) {
      cov_values <- rep(0, region_length)
      message("ℹ️ No reads covering ", chr, ":", start_pos, "-", end_pos, " — filled with zeros.")
    }
    
    result_list[[i]] <- data.frame(
      Chromosome = chr,
      Position = start_pos:end_pos,
      Coverage = cov_values,
      Region = ifelse(is.null(names(regions)[i]), paste0(chr, ":", start_pos, "-", end_pos), names(regions)[i]),
      Filename = tools::file_path_sans_ext(basename(bam_file))
    )
  }
  
  if (length(result_list) == 0) return(NULL)
  do.call(rbind, result_list)
}




# 5. Compute coverage for all BAMs ----------------------------
coverage_list <- lapply(bam_files, get_coverage_in_region, regions = regions)
df_final <- do.call(rbind, coverage_list)
df_final <- df_final[order(df_final$Chromosome, df_final$Position), ]

# 6. Extract and simplify Threshold from filenames ------------
df_final$Threshold <- as.numeric(sub(".*_TH([0-9]+)_.*", "\\1", df_final$Filename))
df_final$Threshold <- factor(df_final$Threshold, levels = sort(unique(df_final$Threshold)))

# 7. Plot Coverage vs. Position -------------------------------
# The x-axis: SNP position
# The y-axis: coverage
# Lines colored by threshold
library(scales)

# Define custom stretch transform
stretch_transform <- trans_new(
  name = "stretch",
  transform = function(x) ifelse(x <= 100, x * 2.5, 250 + (x - 100) * 0.625),
  inverse  = function(x) ifelse(x <= 250, x / 2.5, 100 + (x - 250) / 0.625)
)

# Determine start and end coordinates from BED
x_start <- min(bed$Start)
x_end   <- max(bed$End)

# Generate pretty breaks within the range (e.g., ~15 ticks)
x_breaks <- scales::pretty_breaks(n = 15)(c(x_start, x_end))

p <- ggplot(df_final, aes(x = Position, y = Coverage, color = Threshold, group = Threshold)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 0.2) +  
  geom_hline(yintercept = 40, linetype = "dashed", color = "black", linewidth = 0.4) +  # threshold line
  labs(
    title = "Coverage Across NCF1 Region",
    x = "Genomic Coordinate",
    y = "Coverage",
    color = "Threshold"
  ) +
  # More coordinate labels along X-axis
   scale_x_continuous(
    limits = c(x_start, x_end),
    breaks = x_breaks
  ) +
  # Custom Y scaling and breaks
  scale_y_continuous(
    trans = stretch_transform,
    breaks = c(0, 20, 40, 60, 80, 100, 200, 300, 400),
    limits = c(0, 400)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Display plot
print(p)

# 8. Save plot to file ----------------------------------------
output_plot <- file.path(bam_dir, "NCF1_KM19_NP_TP5_SNP_Coverage_whole_region.png")
ggsave(output_plot, plot = p, width = 10, height = 6, dpi = 300)
message("✅ Plot saved to: ", output_plot)

