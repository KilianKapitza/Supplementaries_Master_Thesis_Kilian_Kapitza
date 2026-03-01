## ---------------------------
##
## Script name: Heatmap pseudogene frequency
##
## Purpose of script:Compares variant data from multiple VCF-derived files to a 
## set of NCF1-specific SNP positions and visualizes allele frequencies in 
## a heatmap.
##
## Author: Annika Vogt, Kilian Kapitza
##
## Date Created: 2025-08-14
##
## ---------------------------
##
## Notes:
##  - The functioning and running script
##  - Merges variant data with a predefined SNP reference list.
##  - Filters out mismatches and keeps the most reliable allele per position/sample.
##  - Missing or invalid frequencies are set to zero.
##  - use NCF1AB_SNP or NCF1AC_SNP bed file
##  
## ---------------------------

library(ggplot2)
library(reshape2)
library(dplyr)
library(tools)
library(forcats)

#1. import the reference BED file containing RhesusBox-specific SNPs (position + nucleotide modification)
bed_file <- file.choose()
bed <- read.table(bed_file, sep = "\t", header = FALSE,
                  col.names = c("Chromosome", "Position", "REF", "ALT"),
                  fill = TRUE, stringsAsFactors = FALSE,
                  fileEncoding = "UTF-16LE")

#2. select all variant files
variant_files <- choose.files(caption = "Select variant files including RhesusBox controls.")
merged_list <- list()

#3. read and merge each variant file with the reference BED list
for (file in variant_files) {
  variant_data <- read.csv(file, sep = ";", header = TRUE, check.names = FALSE,
                           stringsAsFactors = FALSE)
  filename <- file_path_sans_ext(basename(file))
  
  # sicherstellen, dass die Join-Keys gleiche Klasse haben
  bed$Position <- as.character(bed$Position)
  variant_data$Region <- as.character(variant_data$Region)
  
  merged_data <- merge(bed, variant_data[, c("Region", "Reference", "Allele", "Frequency")],
                       by.x = "Position", by.y = "Region", all.x = TRUE)
  
  
  # check if the called allele matches the expected ALT; if not, set frequency to NA
  mismatch <- merged_data$ALT != merged_data$Allele
  merged_data$Frequency[mismatch] <- NA
  
  # convert frequency to numeric (handle commas as decimal points)
  merged_data$Frequency <- as.numeric(gsub(",", ".", merged_data$Frequency))
  
  # store sample file name
  merged_data$Filename <- filename
  merged_list[[file]] <- merged_data
}

#4. combine all merged entries into a single dataframe
df <- bind_rows(merged_list)

#5. rename columns for clarity
colnames(df) <- c("Position", "Chromosome", "REF", "ALT", "Reference", "Allele", "Frequency", "Filename")

#6. filter for valid frequency data and resolve duplicates
df_clean <- df %>%
  group_by(Position, Filename) %>%
  # step 1: Keep only rows with Frequency > 0
  filter(Frequency > 0) %>%
  # step 2: If multiple rows exist, prefer those where Reference == REF
  mutate(match_ref = Reference == REF) %>%
  arrange(desc(match_ref)) %>%     # Prioritize matching reference alleles
  slice(1) %>%                     # Keep only the top-ranked row per group
  ungroup()

#7. identify remaining entries not covered by df_clean (e.g., missing or mismatched variants)
df_rest <- anti_join(df, df_clean, 
                     by = c("Position", "Filename", "Reference", "Allele", "Frequency"))

#8. set frequency to 0 for unmatched or uncalled variants
df_rest$Frequency <- 0

  # identify valid Position+Sample combinations from the cleaned dataset
valid_keys <- df_clean %>% select(Position, Filename) %>% distinct()

  # remove any redundant entries already present in valid_keys
df_rest_filtered <- anti_join(df_rest, valid_keys, by = c("Position", "Filename"))

#9. combine cleaned and remaining filtered data
df_final <- bind_rows(df_clean, df_rest_filtered)
df_final$Filename <- factor(df_final$Filename, levels = sapply(variant_files, function(f) file_path_sans_ext(basename(f))))


#10 extract numeric TH value
df_final$TH_value <- as.numeric(gsub(".*TH([0-9]+)_.*", "\\1", df_final$Filename))

# reorder factor DESCENDING so heatmap is top-to-bottom
df_final$TH_value <- forcats::fct_reorder(
  factor(df_final$TH_value),
  df_final$TH_value,
  .desc = TRUE
)


#11. generate heatmap visualization of allele frequencies
ggplot(df_final, aes(x = factor(Position), y = TH_value, fill = Frequency)) +
  geom_tile(colour = "grey", linewidth = 0.1) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "SNP-Positions", y = "k-mer count per read") +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed(ratio = 6) +
  theme_grey(base_size = 8) +
  scale_x_discrete(breaks = levels(factor(df_final$Position))[seq(1, length(levels(factor(df_final$Position))), by = 5)]) +
  theme(
    legend.key.size = unit(0.35, 'cm'),
    legend.title = element_text(size=5),
    legend.text = element_text(size=4),
    axis.text = element_text(face = "bold"),
    axis.ticks = element_blank(), # element_line(linewidth = 0.2),
    axis.text.x = element_blank(), # element_text(angle = 55, vjust = 1, hjust = 1, size = 3),
    axis.title.x = element_text(margin = margin(t = 6), size=7),
    axis.text.y = element_text(size = 8),
    plot.background = element_blank(),
    panel.border = element_blank()

  )
