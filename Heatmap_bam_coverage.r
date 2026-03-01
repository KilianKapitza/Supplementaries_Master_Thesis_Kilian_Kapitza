# ===============================
# SNP Coverage Heatmap from BAMs
# ===============================
# Author: Kilian Kapitza
# Creates Coverage Heatmap of NCF1 specific SNPs
# conda environment Coverage_env has to be activated

# Load required packages
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(ggplot2)

# 1. Load BED file with SNP positions ------------------------
bed_file <- "/NCF1AB_SNP_UTF8.bed"		# NCF1AB_SNP_UTF8_bed as input to give coordinates of NCF1 specific SNPs (has to be UTF8, NOT UTF16)

bed <- read.table(bed_file, sep = "\t", header = FALSE,
                  col.names = c("Chromosome", "Position", "REF", "ALT"),
                  stringsAsFactors = FALSE)

# 2. List BAM files ------------------------------------------
bam_dir <- "/path/to/your/bam/folder"  # folder with aligned bam files as input
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

# Automatically create BAM index files if missing
for (bam in bam_files) {
  bai <- paste0(bam, ".bai")
  
  if (!file.exists(bai)) {
    message("⚠️  Index missing for: ", basename(bam))
    message("➡️  Creating index...")

    tryCatch(
      {
        Rsamtools::indexBam(bam)
        message("✅ Created index: ", basename(bai))
      },
      error = function(e) {
        stop("❌ Failed to index BAM file: ", bam, "\nError: ", e$message)
      }
    )
  }
}

message("✅ All BAM files correctly indexed.")

# 3. Automatically adjust BED chromosome names to match BAM ----
bam_chr <- names(scanBamHeader(bam_files[1])[[1]]$targets)

if (!all(bed$Chromosome %in% bam_chr)) {
  if (all(paste0("chr", bed$Chromosome) %in% bam_chr)) {
    message("✅ Adding 'chr' prefix to BED chromosomes to match BAM")
    bed$Chromosome <- paste0("chr", bed$Chromosome)
  } else if (all(sub("^chr","", bed$Chromosome) %in% bam_chr)) {
    message("✅ Removing 'chr' prefix from BED chromosomes to match BAM")
    bed$Chromosome <- sub("^chr","", bed$Chromosome)
  } else {
    stop("❌ Chromosome names in BED and BAM do not match. Please check your files.")
  }
} else {
  message("✅ Chromosomes in BED already match BAM.")
}

# Create GRanges object for SNP positions
snps <- GRanges(seqnames = bed$Chromosome,
                ranges = IRanges(start = bed$Position, end = bed$Position))

# 4. Function to get coverage at SNP positions ----------------
get_coverage_at_snps <- function(bam_file, snps) {
  message("Processing: ", bam_file)
  
  # Read alignments
  aln <- readGAlignments(bam_file, use.names = TRUE)
  
  # Compute coverage
  cov <- coverage(aln)
  
  # Extract coverage at each SNP safely
  counts <- sapply(seq_along(snps), function(i) {
    chr <- as.character(seqnames(snps[i]))
    pos <- start(snps[i])
    if (!chr %in% names(cov)) return(0)  # chromosome not present
    as.numeric(cov[[chr]][pos])
  })
  
  # Return as data frame
  data.frame(
    Chromosome = as.character(seqnames(snps)),
    Position = start(snps),
    Coverage = counts,
    Filename = tools::file_path_sans_ext(basename(bam_file))
  )
}

# 5. Apply function to all BAM files --------------------------
coverage_list <- lapply(bam_files, get_coverage_at_snps, snps = snps)
df_final <- do.call(rbind, coverage_list)

# Sort SNPs by chromosome and position
df_final <- df_final[order(df_final$Chromosome, df_final$Position), ]

# 6. Simplify sample names and order by Threshold -----------
# Extract numeric Threshold from filenames
df_final$Threshold <- as.numeric(sub(".*_TH([0-9]+)_.*", "\\1", df_final$Filename))

# Use only the numeric Threshold as label (e.g., "50", "60", ...)
df_final$Filename <- as.character(df_final$Threshold)

# Order y-axis by numeric Threshold (ascending)
df_final$Filename <- factor(df_final$Filename,
                            levels = as.character(sort(unique(df_final$Threshold))))

# 7. Plot as heatmap ------------------------------------------

library(scales)

# Preserve ordering by Threshold
df_final$Threshold <- as.numeric(sub(".*_TH([0-9]+)_.*", "\\1", df_final$Filename))
df_final$Filename <- factor(paste0("Threshold ", df_final$Threshold),
                            levels = paste0("Threshold ", sort(unique(df_final$Threshold))))

# Define segment boundaries
max_cov <- max(df_final$Coverage, na.rm = TRUE)
if (is.na(max_cov) || max_cov < 500) max_cov <- 500  # set at least to 500 for consistent scale

# Segmented color palette:
# 0–40 (black→gray), 40–100 (yellow→red), 100–max (light blue→dark blue)
cols <- c("#EE6677", "#CCBB44", "#66CCEE", "#4477AA", "#55A868", "#1B7837")
values <- c(0, 40, 40.001, 100, 100.001, max_cov)

# Normalize segment positions (between 0 and 1)
norm_values <- (values - min(values)) / (max(values) - min(values))

p <- ggplot(df_final, aes(x = factor(Position, levels = unique(Position)),
                          y = Filename, fill = Coverage)) +
  geom_tile(colour = "grey90", linewidth = 0.1) +
  scale_fill_gradientn(
    colours = cols,
    values = norm_values,
    limits = c(0, max_cov),
    breaks = c(0, 40, 100, 200, 300, 400, 500),
    labels = c(0, 40, 100, 200, 300, 400, 500),
    name = "Coverage"
  ) +
  scale_y_discrete(labels = function(x) sub("^Threshold\\s+", "", x)) +
  labs(title = "SNP Coverage Across Thresholds",
       x = "SNP Position",
       y = "Threshold") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1, size = 6),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8)
  )




# Display plot
print(p)

# 8. Save plot to file ----------------------------------------
output_plot <- file.path(bam_dir, "SNP_Coverage_Heatmap_Thresholds_AB_KM19_NP_TP5.png")
ggsave(output_plot, plot = p, width = 10, height = 6, dpi = 300)
message("✅ Heatmap saved to: ", output_plot)

