# Adjustment of k-mer count per read threshold of fastq files
# Author: Kilian Kapitza
# Über gene_detection_fastq.sh gefilterte Fastq-Dateien erhalten einen match score in ihrem Header. 
# Dieses Skript ist dafür da den ursprünglich genutzten Threshold nachträglich zu erhöhen. Dabei werden alle Reads, deren match score geringer als der neu angelegte threshold ist, aussortiert.
# Ausführung über Ausführung_FilterThreshold_Adjustment.sh
# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
})

# Define CLI options
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input FASTQ file"),
  make_option(c("-o", "--output"), type = "character", help = "Output FASTQ file"),
  make_option(c("-t", "--threshold"), type = "integer", help = "Match threshold"),
  make_option(c("-p", "--plot"), type = "character", default = "histogram.png", help = "Output plot file [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check for missing required arguments
if (is.null(opt$input) || is.null(opt$output) || is.null(opt$threshold)) {
  print_help(opt_parser)
  stop("❌ Please provide --input, --output, and --threshold", call. = FALSE)
}

# Initialize
input_fastq <- opt$input
output_fastq <- opt$output
threshold <- opt$threshold
plot_file <- opt$plot

cat("🔍 Reading:", input_fastq, "\n")
cat("🔎 Threshold:", threshold, "\n")
cat("📤 Writing filtered FASTQ to:", output_fastq, "\n")
cat("📊 Saving histogram to:", plot_file, "\n")

# Open input/output files
con_in <- file(input_fastq, open = "r")
con_out <- file(output_fastq, open = "w")

match_scores <- c()
count <- 0
kept <- 0

repeat {
  lines <- readLines(con_in, n = 4)
  if (length(lines) < 4) break

  header <- lines[1]
  seq <- lines[2]
  plus <- lines[3]
  qual <- lines[4]

  # Extract match score
  match_score <- as.numeric(sub(".*matches=([0-9]+).*", "\\1", header))
  match_scores <- c(match_scores, match_score)
  count <- count + 1

  if (!is.na(match_score) && match_score >= threshold) {
    writeLines(c(header, seq, plus, qual), con_out)
    kept <- kept + 1
  }
}

close(con_in)
close(con_out)

cat("✅ Total reads processed:", count, "\n")
cat("✅ Reads above threshold:", kept, "\n")

# Plot histogram
plot_data <- data.frame(matches = match_scores)

p <- ggplot(plot_data, aes(x = matches)) +
  geom_histogram(binwidth = 10, fill = "#69b3a2", color = "black") +
  geom_vline(xintercept = threshold, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = threshold, y = max(table(cut(match_scores, breaks = 50))) * 0.9,
           label = paste("Threshold =", threshold), color = "red", hjust = -0.1, size = 4) +
  theme_minimal() +
  labs(title = "Distribution of Matches per Read",
       x = "Number of Matches",
       y = "Read Count")

ggsave(plot_file, plot = p, width = 8, height = 5)

cat("📈 Histogram saved to:", plot_file, "\n")

