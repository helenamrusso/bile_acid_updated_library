library(readr)
library(tidyr)
library(ape)
library(taxize)
library(purrr)
library(dplyr)


# Make Tree with Species

# Declare Variables
# .tsv file that contains all organisms whose genera will be in newick file
input_file <- "C:/Users/megan/OneDrive/Desktop/DorresteinLab/TaxaTree/all_sampleinformation.tsv"
# save_file and final_file contain information extracted from NCBI for each batch
save_file <- "ncbi_classifications_species_partial.rds"   # intermediate save
final_file <- "ncbi_classifications_species_final.rds"
# batch size for processing samples
batch_size <- 5            
# number of retries for failed chunks
max_retries <- 3            
# api key to get information from NCBI (domain, kingdom, clade, etc.)
Sys.setenv(ENTREZ_KEY = "c985981b0a043989665e4dfc0f42edd2ac08")

# Load data to build newick file
df <- read_tsv(input_file)

# Split NCBITaxonomy into TaxaID and OrganismName
df_separated <- df %>%
  separate(NCBITaxonomy, into = c("TaxaID", "OrganismName"), sep = "\\|", fill = "right", extra = "merge")

# Filter valid entries 
df_filtered <- df_separated[df_separated$TaxaID != "missing value" & !is.na(df_separated$OrganismName), ]

# Extract genus (first word only)
#df_filtered$Genus <- sub(" .*", "", df_filtered$OrganismName)

# Remove repeats (one TaxaID per genus)
df_unique <- df_filtered[!duplicated(df_filtered$TaxaID), ]
taxa_ids <- as.numeric(df_unique$TaxaID)
names(taxa_ids) <- df_unique$OrganismName

# Split into chunks and batch size
# split function split(dataframe, splitting into batches)
# the ceiling function divides the length of taxa_ids by batch_size and puts the remainder into the last batch (e.g. last batch might just have 3 values)
taxa_chunks <- split(taxa_ids, ceiling(seq_along(taxa_ids) / batch_size))
# total_chunks represents the total number of batches
total_chunks <- length(taxa_chunks)

# Load previous progress if it exists (e.g. some batches were already processed)
if (file.exists(save_file)) {
  message("Resuming from saved file...")
  # this line reads the RDS file (the type of file that save_file is) back into R
  # so that you can pick up where you left off
  all_results <- readRDS(save_file)
} else {
  # this line of code lists the number of batches so that you know it starts at beginning
  all_results <- vector("list", total_chunks)
}

# Function to process one batch safely
# I am guessing taxa_chunks is chunk and i is each value or batch in it
process_chunk <- function(chunk, i) {
  message("Processing chunk ", i, " of ", total_chunks, " (", length(chunk), " taxa)")
  for (attempt in 1:max_retries) {
    # runs NCBI query to get information for each batch
    result <- tryCatch(
      classification(chunk, db = "ncbi"),
      error = function(e) {
        message("  Error in chunk ", i, " attempt ", attempt, ": ", e$message)
        return(NULL)
      }
    )
    if (!is.null(result)) {
      message("  Chunk ", i, " completed successfully")
      return(result)
    }
    Sys.sleep(2)  #  delay time before retry
  }
  message("  Skipping chunk ", i, " after ", max_retries, " attempts")
  return(NULL)
}

# Main loop for batches
for (i in seq_len(total_chunks)) {
  if (!is.null(all_results[[i]])) {
    message("Skipping chunk ", i, " (already done)")
    next
  }
  all_results[[i]] <- process_chunk(taxa_chunks[[i]], i)
  saveRDS(all_results, save_file)
  Sys.sleep(1)
}

# Retry failed batches at the end
skipped <- which(sapply(all_results, is.null))
if (length(skipped) > 0) {
  message("Retrying ", length(skipped), " skipped chunks...")
  for (i in skipped) {
    all_results[[i]] <- process_chunk(taxa_chunks[[i]], i)
    saveRDS(all_results, save_file)
    Sys.sleep(1)
  }
}

# Finalize results
all_results_combined <- do.call(c, all_results)
saveRDS(all_results_combined, final_file)
message("All results saved as ", final_file)

# Clean results & build tree
NCBIClass_clean <- all_results_combined[!sapply(all_results_combined, function(x) is.null(x) || nrow(x) == 0)]
#rank_counts <- sapply(NCBIClass_clean, nrow)
#max_rank <- as.numeric(names(which.max(table(rank_counts))))
#NCBIClass_filtered <- NCBIClass_clean[rank_counts == max_rank]

# Build genus-level tree
#tree_data <- class2tree(NCBIClass_filtered)
# Remove duplicates from NCBIClass_clean by keeping the first occurrence for each name
NCBIClass_nodup <- NCBIClass_clean[!duplicated(names(NCBIClass_clean))]

### Now only selecting for organisms that are in the animal kingdom/Metazoa
# Extract data that belongs to domain, kingdom, phylum, class, order, family, genus, and species
standard_ranks <- c("domain", "kingdom", "phylum", "class", 
                    "order", "family", "genus", "species")

# Filter each organism to standard_ranks
NCBIClass_filtered <- lapply(NCBIClass_nodup, function(df) {
  df[df$rank %in% standard_ranks, , drop = FALSE]
})

# Keep only lineages with all standard_ranks values
NCBIClass_complete <- NCBIClass_filtered[sapply(NCBIClass_filtered, function(x) {
  all(standard_ranks %in% x$rank)
})]

# 🔹 NEW: Filter to only keep organisms in kingdom Metazoa
NCBIClass_metazoa <- NCBIClass_complete[sapply(NCBIClass_complete, function(x) {
  any(x$rank == "kingdom" & x$name == "Metazoa")
})]

cat("Original taxa:", length(NCBIClass_nodup), "\n")
cat("With all standard ranks:", length(NCBIClass_complete), "\n")
cat("With kingdom = Metazoa:", length(NCBIClass_metazoa), "\n")

# Build the tree
tree_data <- class2tree(NCBIClass_metazoa)
phylo_tree <- tree_data$phylo

# Replace spaces with underscores before saving
phylo_tree$tip.label <- gsub(" ", "_", phylo_tree$tip.label)
phylo_tree$tip.label <- gsub("\\.x$|\\.y$", "", phylo_tree$tip.label)

dupes <- phylo_tree$tip.label[duplicated(phylo_tree$tip.label)]
if (length(dupes) > 0) {
  phylo_tree <- drop.tip(phylo_tree, which(duplicated(phylo_tree$tip.label)))
}

bad_tips <- which(is.na(phylo_tree$tip.label))
if (length(bad_tips) > 0) {
  phylo_tree <- drop.tip(phylo_tree, bad_tips)
}

# Save tree
write.tree(phylo_tree, file = "tree_of_life_metazoa_species.nw")
message("Species-level tree saved as tree_of_life_metazoa_species.nw")

# Plot the tree
png("tree_of_life_metazoa_species.png", width = 8000, height = 8000, res = 300)  
plot(
  phylo_tree,
  cex = 1.0,
  type = "fan",
  label.offset = 0.5,
  no.margin = TRUE
)
title("Species-level Tree of Life (Metazoa)", cex.main = 2)
dev.off()
