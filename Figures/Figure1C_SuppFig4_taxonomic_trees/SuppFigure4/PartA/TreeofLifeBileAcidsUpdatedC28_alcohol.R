library(ape)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(RColorBrewer)
# FIX: removed library(fields) -- unused in this script (no image.plot() call here)

setwd('C:/Users/megan/OneDrive/Desktop/bile-acid-library-update')
# Load files
tree_file <- "C:/Users/megan/OneDrive/Desktop/bile-acid-library-update/TreeofLife/tree_of_life_species_kingdom_combined_color_noother.nw"
rds_file  <- "C:/Users/megan/OneDrive/Desktop/bile-acid-library-update/TreeofLife/ncbi_classifications_species_final.rds"
ba_file   <- "C:/Users/megan/OneDrive/Desktop/bile-acid-library-update/TreeofLife/BA_scan_counts_per_taxonomy_no_artifacts_alkamine.tsv"


phylo_tree <- read.tree(tree_file)
cat("Loaded tree with", length(phylo_tree$tip.label), "tips\n")


NCBIClass_complete <- readRDS(rds_file)
head(phylo_tree$tip.label)
# Safely build species to kingdom mapping
species_to_kingdom_list <- lapply(NCBIClass_complete, function(df) {
  org <- df[df$rank == "species", "name"]
  k   <- df[df$rank == "kingdom", "name"]
  if (length(org) > 0 && length(k) > 0) {
    org <- gsub(" ", "_", org)
    org <- gsub("\\.x$|\\.y$", "", org)
    tibble(
      TaxaID       = df$tax_id[1],
      OrganismName = org[1],
      Kingdom      = k[1]
    )
  } else {
    NULL
  }
})
species_to_kingdom_df <- bind_rows(species_to_kingdom_list)

# Standardize kingdom names
species_to_kingdom_df$Kingdom <- dplyr::recode(species_to_kingdom_df$Kingdom,
                                               "Bacillati"        = "Bacteria",
                                               "Pseudomonadati"   = "Bacteria",
                                               "Thermotogati"     = "Bacteria",
                                               "Fusobacteriati"   = "Bacteria",
                                               "Methanobacteriati"= "Archaea",
                                               "Thermoproteati"   = "Archaea",
                                               .default = species_to_kingdom_df$Kingdom
)

species_to_kingdom_df$OrganismName <- gsub(" ", "_", species_to_kingdom_df$OrganismName)
species_to_kingdom_df$OrganismName <- gsub("\\.x$|\\.y$", "", species_to_kingdom_df$OrganismName)

cat("Kingdom table built with", nrow(species_to_kingdom_df), "rows\n")

# Filter to match tree tips (robust fix)
species_to_kingdom_df <- species_to_kingdom_df %>%
  filter(OrganismName %in% phylo_tree$tip.label) %>%
  distinct(OrganismName, .keep_all = TRUE)

cat("Filtered kingdom table now has", nrow(species_to_kingdom_df), "rows (should match ~",
    length(phylo_tree$tip.label), ")\n")

# Match kingdoms to tree tips
tip_kingdoms <- setNames(
  species_to_kingdom_df$Kingdom[match(phylo_tree$tip.label,
                                      species_to_kingdom_df$OrganismName)],
  phylo_tree$tip.label
)
tip_kingdoms[is.na(tip_kingdoms)] <- "Other"

# Remove tips with no kingdom assignment
other_tips <- phylo_tree$tip.label[tip_kingdoms == "Other"]
cat("Removing", length(other_tips), "tips with 'Other' kingdom:\n")
print(other_tips)
phylo_tree <- drop.tip(phylo_tree, other_tips)

# Rebuild tip_kingdoms to match pruned tree
tip_kingdoms <- tip_kingdoms[names(tip_kingdoms) %in% phylo_tree$tip.label]

# Sanity checks
cat("Tree tips:", length(phylo_tree$tip.label), "\n")
cat("Mapped kingdoms:", length(tip_kingdoms), "\n")
stopifnot(length(tip_kingdoms) == length(phylo_tree$tip.label))

cat("Kingdom distribution:\n")
print(table(tip_kingdoms))

# Define colors
kingdom_colors <- c(
  "Archaea"        = "#80B1D3",
  "Bacteria"       = "#BEBADA",
  "Fungi"          = "#D9D9D9",
  "Metazoa"        = "#FB8072",
  "Viridiplantae"  = "#CCEBC5"
)

# Set which core_std group to analyze (change this each run)
ba_group <- "C28_alcohol"  # options: "C23", "C24", "C24_alkamine", "C26", "C27", "C28_alcohol"

# FIX: ba_levels and ba_colors are now defined TOGETHER per ba_group, so they
# can't drift out of sync the way the old "## Remember to change these"
# comment risked. This also replaces the old grepl("^ba_group_", core_std)
# prefix match, which would silently double-count: "C24_" as a prefix also
# matches "C24_alkamine_..." entries, so setting ba_group <- "C24" previously
# would have pulled in C24_alkamine rows too. Add the remaining groups'
# core_std levels/colors here as you build their trees.
#
# NOTE: C28 alcohols only have 1-3 side-chain hydroxyl groups (not the
# nonhydroxy->pentahydroxy 6-level scheme used for C23/C24/C26/C27), so this
# group has 3 levels instead of 6. Naming confirmed from a direct hit in your
# own delta-counts TSV ("C28_alcohol_3OH"); "C28_alcohol_1OH" and
# "C28_alcohol_2OH" are inferred from that pattern -- please verify they
# match your actual TSV before running, since a name mismatch here will
# silently filter to zero rows rather than error.
ba_group_config <- list(
  "C27" = list(
    levels = c("C27_nonhydroxy","C27_monohydroxy","C27_dihydroxy",
               "C27_trihydroxy","C27_tetrahydroxy","C27_pentahydroxy"),
    colors = c(
      "C27_nonhydroxy"   = "#1F77B4",
      "C27_monohydroxy"  = "#2CA02C",
      "C27_dihydroxy"    = "#FF7F0E",
      "C27_trihydroxy"   = "#FFDD57",
      "C27_tetrahydroxy" = "#9467BD",
      "C27_pentahydroxy" = "#7F7F7F"
    )
  ),
  "C28_alcohol" = list(
    levels = c("C28_alcohol_1OH", "C28_alcohol_2OH", "C28_alcohol_3OH"),
    colors = c(
      "C28_alcohol_1OH"  = "#E377C2",
      "C28_alcohol_2OH"    = "#17BECF",
      "C28_alcohol_3OH"   = "#DDB892"
    )
  )
)

if (!ba_group %in% names(ba_group_config)) {
  stop("ba_levels/ba_colors not yet defined for ba_group = '", ba_group,
       "'. Add its core_std levels and colors to ba_group_config above.")
}
ba_levels <- ba_group_config[[ba_group]]$levels
ba_colors <- ba_group_config[[ba_group]]$colors


# Load and prepare bile acid data
ba_class <- read_tsv(ba_file)

ba_class_filtered <- ba_class %>%
  filter(!(NCBITaxonomy %in% c("missing value", "N/A")) & !is.na(NCBITaxonomy)) %>%
  filter(core_std %in% ba_levels) %>%   # FIX: exact match against ba_levels, was grepl(paste0("^", ba_group, "_"), core_std)
  separate(NCBITaxonomy, into = c("TaxaID", "OrganismName"),
           sep = "\\|", fill = "right", extra = "merge")

# QC: warn early if the level names above don't actually match anything in
# this TSV, since the filter above fails silently (returns zero rows) rather
# than erroring on a name mismatch.
if (nrow(ba_class_filtered) == 0) {
  warning("No rows matched core_std %in% ba_levels for ba_group = '", ba_group,
          "'. Double-check the exact core_std spelling in ", ba_file, ".")
}

ba_class_filtered$OrganismName <- gsub(" ", "_", ba_class_filtered$OrganismName)
ba_class_filtered$OrganismName <- gsub("\\.x$|\\.y$", "", ba_class_filtered$OrganismName)

ba_class_summary <- ba_class_filtered %>%
  group_by(OrganismName, core_std) %>%
  summarise(total = sum(scan_counts), .groups = "drop")

# Factor ordering
ba_class_summary$core_std <- factor(
  ba_class_summary$core_std,
  levels = names(ba_colors),
  ordered = TRUE
)


species_totals <- split(ba_class_summary, ba_class_summary$OrganismName)

# Draw kingdom ring
draw_kingdom_ring <- function(tip_data, tree_radius, kingdom_colors, ring_width = 0.1) {
  r_inner <- 1.05 * tree_radius
  r_outer <- (1.05 + ring_width) * tree_radius
  
  # tip_data is already sorted by angle
  n <- nrow(tip_data)
  
  # Find contiguous runs of the same kingdom
  run_kingdom <- tip_data$kingdom[1]
  run_start   <- 1
  
  for (i in 2:(n + 1)) {
    current <- if (i <= n) tip_data$kingdom[i] else ""
    if (current != run_kingdom || i == n + 1) {
      # Draw arc for this run
      k <- run_kingdom
      if (k %in% names(kingdom_colors)) {
        start_angle <- tip_data$angle[run_start]
        end_angle   <- tip_data$angle[min(i - 1, n)]
        theta <- seq(start_angle, end_angle, length.out = 500)
        polygon(c(r_inner*cos(theta), rev(r_outer*cos(theta))),
                c(r_inner*sin(theta), rev(r_outer*sin(theta))),
                col = adjustcolor(kingdom_colors[k], alpha.f = 0.6),
                border = NA)
      }
      run_kingdom <- current
      run_start   <- i
    }
  }
}


# --- Compute layout BEFORE opening the PDF device ---
# This avoids creating blank pages from plot.new()/plot(..., plot=FALSE)

plot(phylo_tree, type = "fan", cex = 0.4, show.tip.label = FALSE, plot = FALSE)

coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tree_radius <- max(sqrt(coords$xx^2 + coords$yy^2))

n_tips <- length(phylo_tree$tip.label)
tip_x <- coords$xx[1:n_tips]
tip_y <- coords$yy[1:n_tips]
tip_angles <- (atan2(tip_y, tip_x) + 2*pi) %% (2*pi)

names(tip_kingdoms) <- NULL

tip_data <- data.frame(
  species = phylo_tree$tip.label,
  angle   = tip_angles,
  kingdom = tip_kingdoms,
  stringsAsFactors = FALSE
)
tip_data <- tip_data[order(tip_data$angle), ]

# --- NOW open the PDF and draw only the real plot ---
pdf("C28_alcohol_TreeofLife.pdf", width = 30, height = 30)

# Draw base tree
plot(phylo_tree, type = "fan", cex = 0.4,
     show.tip.label = FALSE, no.margin = TRUE,
     x.lim = c(-1.8, 1.8) * tree_radius,
     y.lim = c(-1.8, 1.8) * tree_radius)

# Add kingdom ring
draw_kingdom_ring(tip_data, tree_radius, kingdom_colors, ring_width = 0.1)

# Overlay bile acid bars
# FIX (v2): the previous version summed log10(count_j + 1) across classes per
# species -- but sum(log10(x_i+1)) is NOT proportional to total signal: a
# species with the same total spread across more classes gets an inflated
# bar purely from summing more log terms (e.g. one class at 50 -> log10(51)
# ~= 1.71; five classes at 10 each, same total of 50 -> 5*log10(11) ~= 5.21).
# My prior fix only recalibrated the scale factor, which didn't address this
# -- the relative distortion between species was still there regardless of
# scaling. Correct approach (matching the proportional-split pattern already
# used for the BA bars in the Metazoa tree scripts): compute ONE
# log10(total+1) per species for the bar's total length, then split that
# length proportionally by each class's share of the raw (unlogged) total.
species_log_totals <- ba_class_summary %>%
  group_by(OrganismName) %>%
  summarise(raw_total = sum(total), log_total = log10(sum(total) + 1), .groups = "drop")

max_log_total <- max(species_log_totals$log_total, na.rm = TRUE)
scale_factor  <- 0.5 * tree_radius / max_log_total
ring_outer <- (1.05 + 0.1) * tree_radius

log_total_lookup <- setNames(species_log_totals$log_total, species_log_totals$OrganismName)
raw_total_lookup  <- setNames(species_log_totals$raw_total,  species_log_totals$OrganismName)

for (i in seq_along(phylo_tree$tip.label)) {
  sp <- phylo_tree$tip.label[i]
  if (sp %in% names(species_totals)) {
    dat <- species_totals[[sp]]
    dat <- dat[order(dat$core_std), ]
    x <- coords$xx[i]
    y <- coords$yy[i]
    angle <- atan2(y, x)
    
    bar_total_len <- log_total_lookup[[sp]] * scale_factor
    dat$prop <- dat$total / raw_total_lookup[[sp]]   # proportional share of raw (unlogged) total
    
    r_start <- ring_outer
    for (j in seq_len(nrow(dat))) {
      seg_len <- dat$prop[j] * bar_total_len
      r_end <- r_start + seg_len
      segments(r_start*cos(angle), r_start*sin(angle),
               r_end*cos(angle),   r_end*sin(angle),
               col = ba_colors[dat$core_std[j]],
               lwd = 2)
      r_start <- r_end
    }
  }
}

# Add legends and title
mtext("Tree of Life with C28 Bile Alcohols",
      side = 3, outer = TRUE, line = 2, cex = 4)

legend("bottomleft", legend = names(kingdom_colors),
       fill = kingdom_colors, border = "black", cex = 2.5, bty = "n", title = "Kingdoms")

legend("bottomright", legend = names(ba_colors),
       fill = ba_colors, border = "black", cex = 2.5, bty = "n", title = "Bile Acid Classes")



dev.off()

cat("Figure saved \n")