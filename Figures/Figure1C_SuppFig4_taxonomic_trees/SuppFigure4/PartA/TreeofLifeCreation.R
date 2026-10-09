library(readr)
library(tidyr)
library(ape)
library(taxize)
library(purrr)
library(dplyr)


# Make Tree with Species & Map Species by Kingdom

# Declare Variables
# .tsv file that contains all organisms whose genera will be in newick file
input_file <- "C:/Users/megan/OneDrive/Desktop/DorresteinLab/TaxaTree/all_sampleinformation.tsv"
# save_file and final_file contain information extracted from NCBI for each batch
save_file <- "ncbi_classifications_species_partial.rds"   # intermediate save
final_file <- "ncbi_classifications_species_final.rds"
# batch size for processing samples - make it small to reduce errors
batch_size <- 5            
# number of retries for failed chunks
max_retries <- 3            
# api key to get information from NCBI (domain, kingdom, clade, etc.)
# get your own key to get information from NCBI
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
# taxa_chunks is chunk and i is each value or batch in it
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

NCBIClass_clean <- all_results_combined[!sapply(all_results_combined, function(x) {
  is.null(x) || nrow(x) == 0
})]

# Remove duplicates
NCBIClass_nodup <- NCBIClass_clean[!duplicated(names(NCBIClass_clean))]

# Standard taxonomic ranks needed for all species
standard_ranks <- c("domain", "kingdom", "phylum", "class", 
                    "order", "family", "genus", "species")

# Keep only these ranks for all speices (some species contain others that will be removed for so that the tree can be created)
NCBIClass_filtered <- lapply(NCBIClass_nodup, function(df) {
  df[df$rank %in% standard_ranks, , drop = FALSE]
})

# Keep lineages with all standard ranks
NCBIClass_complete <- NCBIClass_filtered[sapply(NCBIClass_filtered, function(x) {
  all(standard_ranks %in% x$rank)
})]

cat("Original taxa:", length(NCBIClass_nodup), "\n")
cat("With all standard ranks:", length(NCBIClass_complete), "\n")

# The code above retrieves the taxonomic information from NCBI
# the code below will map each organism to a kingdom in a colored kingdom ring around the tree and save the tree as a .nwk and .png file



# Build species → kingdom mapping

# Build species_to_kingdom first
species_to_kingdom <- sapply(NCBIClass_complete, function(df) {
  org <- df[df$rank == "species", "name"]   # species name
  k   <- df[df$rank == "kingdom", "name"]   # kingdom name
  org <- gsub(" ", "_", org)                # match tree format
  org <- gsub("\\.x$|\\.y$", "", org)
  if (length(org) > 0 && length(k) > 0) {
    setNames(k[1], paste0(df$tax_id[1], ".", org[1]))
  } else {
    NULL
  }
})


# Convert species_to_kingdom into a data.frame
species_to_kingdom_df <- data.frame(
  TaxaID       = sub("\\..*", "", names(species_to_kingdom)),           
  OrganismName = sub("^[0-9]+\\.", "", names(species_to_kingdom)),      
  Kingdom      = unname(species_to_kingdom),
  stringsAsFactors = FALSE
)

# Collapse subgroups into major categories
species_to_kingdom_df$Kingdom <- dplyr::recode(species_to_kingdom_df$Kingdom,
                                               "Bacillati"        = "Bacteria",
                                               "Pseudomonadati"   = "Bacteria",
                                               "Thermotogati"     = "Bacteria",
                                               "Fusobacteriati"   = "Bacteria",
                                               "Methanobacteriati"= "Archaea",
                                               "Thermoproteati"   = "Archaea",
                                               .default = species_to_kingdom_df$Kingdom
)

# Clean OrganismName to match tip labels
species_to_kingdom_df$OrganismName <- gsub(" ", "_", species_to_kingdom_df$OrganismName)
species_to_kingdom_df$OrganismName <- gsub("\\.x$|\\.y$", "", species_to_kingdom_df$OrganismName)
species_to_kingdom_df$OrganismName <- gsub("^\\.", "", species_to_kingdom_df$OrganismName)

head(species_to_kingdom_df)
# Build the tree

tree_data <- class2tree(NCBIClass_complete) # takes a LONG time
phylo_tree <- tree_data$phylo

grep("symbiosum|scindens|lactaris|innocuum|papulosa", phylo_tree$tip.label, value = TRUE)


phylo_tree$tip.label <- gsub(" ", "_", phylo_tree$tip.label)
phylo_tree$tip.label <- gsub("\\.x$|\\.y$", "", phylo_tree$tip.label)
phylo_tree$tip.label <- gsub("\\.(x|y)\\.\\d+$", "", phylo_tree$tip.label)

dupes <- phylo_tree$tip.label[duplicated(phylo_tree$tip.label)]
if (length(dupes) > 0) {
  message("Dropping ", length(dupes), " duplicate tips")
  phylo_tree <- drop.tip(phylo_tree, which(duplicated(phylo_tree$tip.label)))
}

bad_tips <- which(is.na(phylo_tree$tip.label))
if (length(bad_tips) > 0) {
  phylo_tree <- drop.tip(phylo_tree, bad_tips)
}

# Match kingdoms to tree tips
tip_kingdoms <- setNames(
  species_to_kingdom_df$Kingdom[match(phylo_tree$tip.label, species_to_kingdom_df$OrganismName)],
  phylo_tree$tip.label
)

# Replace any missing matches
tip_kingdoms[is.na(tip_kingdoms)] <- "Other"

# See distribution
print(table(tip_kingdoms))

# Assign colors
library(RColorBrewer)
all_kingdoms <- unique(tip_kingdoms)
#kingdom_colors <- setNames(brewer.pal(length(all_kingdoms), "Set3"), all_kingdoms)


kingdom_colors <- c(
  "Archaea"   = "#80B1D3",
  "Bacteria"  = "#BEBADA",
  "Fungi"    = "#D9D9D9",
  "Metazoa"   = "#FB8072",
  "Viridiplantae" = "#CCEBC5"
)
# Apply coloring
#edge_colors <- assign_edge_colors(phylo_tree, tip_kingdoms)



# Save and plot
# Write newick file to be called later
write.tree(phylo_tree, file = "tree_of_life_species_kingdom_combined_color_noother_check.nw")
message("Species-level tree saved as tree_of_life_species_kingdom_combined_color_noother_check.nw")

png("tree_of_life_species_kingdom_combined_background_noother.png", width = 6000, height = 6000, res = 300)

# First get coordinates for tips (no drawing)
plot(phylo_tree, type = "fan", cex = 0.3, show.tip.label = FALSE, plot = FALSE)
coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tree_radius <- max(sqrt(coords$xx^2 + coords$yy^2))

# Angles for tips
n_tips <- length(phylo_tree$tip.label)
tip_x <- coords$xx[1:n_tips]
tip_y <- coords$yy[1:n_tips]

tip_angles <- (atan2(tip_y, tip_x) + 2*pi) %% (2*pi)

tip_data <- data.frame(
  species = phylo_tree$tip.label,
  angle   = tip_angles,
  kingdom = tip_kingdoms[phylo_tree$tip.label],
  stringsAsFactors = FALSE
)
tip_data <- tip_data[order(tip_data$angle), ]

# Function to draw outer-ring wedges
draw_kingdom_ring <- function(tip_data, tree_radius, kingdom_colors, ring_width = 0.05) {
  kingdoms <- unique(tip_data$kingdom)
  for (k in kingdoms) {
    subset <- tip_data[tip_data$kingdom == k, ]
    if (nrow(subset) == 0) next
    start_angle <- min(subset$angle)
    end_angle   <- max(subset$angle)
    theta <- seq(start_angle, end_angle, length.out = 500)
    
    # make sure the kingdom ring is outside the tree radius
    r_inner <- 1.05 * tree_radius
    r_outer <- (1.05 + ring_width) * tree_radius
    
    polygon(c(r_inner*cos(theta), rev(r_outer*cos(theta))),
            c(r_inner*sin(theta), rev(r_outer*sin(theta))),
            col = adjustcolor(kingdom_colors[k], alpha.f = 0.6),
            border = NA)
  }
}

par(oma = c(1, 1, 6, 1))
# Draw base tree with tighter limits (so kigndom ring doesn’t push tree off page)
plot(phylo_tree, type = "fan", cex = 0.3,
     show.tip.label = FALSE, no.margin = TRUE,
     x.lim = c(-1.2, 1.2) * tree_radius,
     y.lim = c(-1.2, 1.2) * tree_radius)

# draw kingdom ring
draw_kingdom_ring(tip_data, tree_radius, kingdom_colors, ring_width = 0.1)

# Title and legend
mtext("Tree of Life", 
      side = 3, outer = TRUE, line = 2, cex = 5)
legend("bottomleft", legend = names(kingdom_colors),
       fill = kingdom_colors, border = "black", cex = 2, bty = "n", title = "Kingdoms")

dev.off()