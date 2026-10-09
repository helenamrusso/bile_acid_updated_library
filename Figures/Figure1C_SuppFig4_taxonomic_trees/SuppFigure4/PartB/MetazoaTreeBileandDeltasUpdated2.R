### Metazoa Tree + Phyla Ring + Normalized Bile Acids + Log2 Delta Masses + top 10 organisms printed
### USE THIS ONE
library(readr)
library(tidyr)
library(ape)
library(purrr)
library(dplyr)
library(fields)
library(extrafont)
library(RColorBrewer)   # FIX: was missing -- brewer.pal() below would error without this

setwd('C:/Users/megan/OneDrive/Desktop/bile-acid-library-update')

# Select BA Group
#### You have to manually change the BA group for each tree
ba_group <- "C27"  # options: "C23", "C24", "C24_alkamine", "C26", "C27", "C28_alcohol"

# FIX: ba_levels is now derived from ba_group instead of being hardcoded, so
# switching ba_group above actually changes what gets plotted. Extend this
# list if/when you reuse the script for the other five libraries -- each
# core's core_std naming convention needs to match what's in your TSVs.
ba_levels <- switch(
  ba_group,
  "C27" = c("C27_nonhydroxy", "C27_monohydroxy", "C27_dihydroxy",
            "C27_trihydroxy", "C27_tetrahydroxy", "C27_pentahydroxy"),
  stop("ba_levels not yet defined for ba_group = '", ba_group, "'. Add its core_std levels above.")
)

# Load Metazoa tree (created in MakeMetazoaTree.R)
tree_life <- read.tree("C:/Users/megan/OneDrive/Desktop/DorresteinLab/TaxaTree/tree_of_life_metazoa_species.nw")

# FIX: multi2di() randomly resolves polytomies (multifurcating nodes) into
# bifurcations each time it's called -- this changes internal node/edge
# order, which changes tip angular position in the fan layout below. The
# seed MUST be set immediately before this call to make tip placement
# reproducible across runs. (Previously set.seed() was called much later,
# after multi2di() had already consumed randomness, so it had no effect.)
set.seed(123)
tree_life <- multi2di(tree_life)

rds_file  <- "C:/Users/megan/OneDrive/Desktop/bile-acid-library-update/MetazoaTree/ncbi_classifications_species_final.rds"
NCBIClass_metazoa <- readRDS(rds_file)

# Get coordinates of tree
plot(tree_life, type = "fan", cex = 0.3, show.tip.label = FALSE, plot = FALSE)
coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tree_radius <- max(sqrt(coords$xx^2 + coords$yy^2))

# Angles for tips
n_tips <- length(tree_life$tip.label)
tip_x <- coords$xx[1:n_tips]
tip_y <- coords$yy[1:n_tips]
tip_angles <- (atan2(tip_y, tip_x) + 2*pi) %% (2*pi)


# Assign Phyla info for Metazoa

species_to_phylum_df <- do.call(rbind, lapply(NCBIClass_metazoa, function(df) {
  sp <- df[df$rank == "species", "name"]
  ph <- df[df$rank == "phylum", "name"]
  if (length(sp) > 0 && length(ph) > 0) {
    data.frame(
      OrganismName = gsub(" ", "_", sp[1]),
      Phylum = ph[1],
      stringsAsFactors = FALSE
    )
  } else NULL
}))

species_to_phylum_df <- species_to_phylum_df %>%
  distinct(OrganismName, .keep_all = TRUE)

tip_phyla <- setNames(
  species_to_phylum_df$Phylum[match(tree_life$tip.label, species_to_phylum_df$OrganismName)],
  tree_life$tip.label
)
tip_phyla[is.na(tip_phyla)] <- "Other"

# Manual override: Pycnonotus melanicterus is absent from NCBIClass_metazoa
# entirely (confirmed -- not a name-formatting mismatch), so it falls into
# "Other" by default. Correcting it here since we know from taxonomy it's
# Chordata (Aves).
tip_phyla["Pycnonotus_melanicterus"] <- "Chordata"

# QC: which tips have no phylum assignment and will fall into "Other"
other_tips <- names(tip_phyla)[tip_phyla == "Other"]
cat(length(other_tips), "tips assigned to 'Other' phylum:\n")
print(other_tips)

# Phylum colors (your palette)
phylum_colors <- c(
  "Annelida"      = "#007A99",
  "Arthropoda"    = "#00A9CC",
  "Chordata"      = "#65EFFF",
  "Cnidaria"      = "#CCFDFF",
  "Echinodermata" = "#F2DACD",   # FIX: was misspelled "Enchinodermata" -- this string must also
  # match however NCBI/your classification data spells the phylum,
  # or matched tips will fall through to "Other" below
  "Mollusca"      = "#D8AF97",
  "Nematoda"      = "#662F00",
  "Porifera"      = "#331900",
  "Other"         = "#B0B0B0"    # FIX: explicit fallback color so any tip without full phylum
  # info draws a visible (grey) ring segment instead of silently
  # drawing nothing (adjustcolor(NA, ...) was invisible before)
)


# Prepare Bile Acids
ba_class <- read_tsv("C:/Users/megan/OneDrive/Desktop/bile-acid-library-update/MetazoaTree/BA_scan_counts_per_taxonomy_no_artifacts_alkamine.tsv")

ba_colors <- c(
  "C27_nonhydroxy"   = "#1F77B4",
  "C27_monohydroxy"  = "#2CA02C",
  "C27_dihydroxy"    = "#FF7F0E",
  "C27_trihydroxy"   = "#FFDD57",
  "C27_tetrahydroxy" = "#9467BD",
  "C27_pentahydroxy" = "#7F7F7F"
)

ba_class_filtered <- ba_class %>%
  filter(!(NCBITaxonomy %in% c("missing value", "N/A")) & !is.na(NCBITaxonomy)) %>%
  filter(core_std %in% ba_levels) %>%
  separate(NCBITaxonomy, into = c("TaxaID", "OrganismName"),
           sep = "\\|", fill = "right", extra = "merge")

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

# NOTE: the >3-per-class eligibility threshold block that used to live here has
# been removed. It computed `eligible_species` but was never applied to
# `ba_class_summary` (the filter line was commented out) -- if you want to
# reintroduce a threshold, add it back deliberately and apply it here.

species_totals <- split(ba_class_summary, ba_class_summary$OrganismName)


# Load Delta Masses File
deltas <- read_tsv("C:/Users/megan/OneDrive/Desktop/bile-acid-library-update/MetazoaTree/BA_delta_counts_per_taxonomy_no_artifacts_alkamine.tsv")

deltas_filtered <- deltas %>%
  filter(!(NCBITaxonomy %in% c("missing value", "N/A")) & !is.na(NCBITaxonomy)) %>%
  filter(core_std %in% ba_levels) %>%
  separate(NCBITaxonomy, into = c("TaxaID", "OrganismName"),
           sep = "\\|", fill = "right", extra = "merge")

deltas_filtered$OrganismName <- gsub(" ", "_", deltas_filtered$OrganismName)
deltas_filtered$OrganismName <- gsub("\\.x$|\\.y$", "", deltas_filtered$OrganismName)

deltas_summary <- deltas_filtered %>%
  group_by(OrganismName, core_std) %>%
  summarise(num_deltas = sum(delta_counts), .groups = "drop")

# FIX -- THE MAIN BUG: previously this was
#   species_totals_deltas <- setNames(deltas_summary$num_deltas, deltas_summary$OrganismName)
# which produced a vector with REPEATED names (one entry per core_std class per
# organism). Downstream, species_totals_deltas[sp] silently returns only the
# FIRST matching entry for that name, not the organism's true total across all
# six C24 hydroxylation classes. Confirmed against your TSV: 2,020 of 2,641
# organisms (76.5%) have counts split across more than one class, and for your
# top species this meant plotting a fraction of the true value (e.g.
# M. musculus showed 13% of its true total, H. sapiens 47%). Fixed by summing
# across core_std per organism BEFORE building the named vector:
species_totals_deltas <- deltas_summary %>%
  group_by(OrganismName) %>%
  summarise(total = sum(num_deltas), .groups = "drop") %>%
  { setNames(.$total, .$OrganismName) }

# log2 transformation of delta masses data
log_species_totals_deltas <- log2(species_totals_deltas + 1)
max_log_val <- max(log_species_totals_deltas, na.rm = TRUE)

# Palette for delta masses
delta_palette <- colorRampPalette(brewer.pal(9, "Reds"))(100)


# Plot Metazoa tree

pdf("C27_check.pdf",  width = 30, height = 35,
    family = "Helvetica")

par(oma = c(1,1,6,6), family = "Helvetica")

# Fix limits
# FIX: widened from ±1.8 to ±2.3*tree_radius. With the aggregation bug fixed,
# t_max now reflects true per-species totals, so the tallest delta-mass bars
# (e.g. Mus_musculus, Homo_sapiens) reach close to their full L_max length.
# Ring_outer + bile-acid bars + separator + L_max sums to roughly
# 2.16*tree_radius, which exceeded the old 1.8 limit and clipped the tallest
# bars against the plot boundary.
plot(tree_life, type = "fan", cex = 0.3, show.tip.label = FALSE,
     no.margin = TRUE,
     x.lim = c(-2.3,2.3)*tree_radius,
     y.lim = c(-2.3,2.3)*tree_radius)
coords <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Angles for tips
tip_angles <- (atan2(coords$yy[1:n_tips], coords$xx[1:n_tips]) + 2*pi) %% (2*pi)

tip_data <- data.frame(
  species = tree_life$tip.label,
  angle   = tip_angles,
  phylum  = tip_phyla[tree_life$tip.label],
  stringsAsFactors = FALSE
)
tip_data <- tip_data[order(tip_data$angle), ]

# Draw Phyla ring
draw_phylum_ring <- function(tip_data, tree_radius, phylum_colors, ring_width = 0.1) {
  phyla <- unique(tip_data$phylum)
  for (p in phyla) {
    subset <- tip_data[tip_data$phylum == p, ]
    if (nrow(subset) == 0) next
    
    start_angle <- min(subset$angle, na.rm = TRUE)
    end_angle   <- max(subset$angle, na.rm = TRUE)
    theta <- seq(start_angle, end_angle, length.out = 500)
    
    r_inner <- 1.05 * tree_radius
    r_outer <- (1.05 + ring_width) * tree_radius
    
    polygon(c(r_inner*cos(theta), rev(r_outer*cos(theta))),
            c(r_inner*sin(theta), rev(r_outer*sin(theta))),
            col = adjustcolor(phylum_colors[p], alpha.f = 0.6),
            border = NA)
  }
}

ring_width <- 0.1
draw_phylum_ring(tip_data, tree_radius, phylum_colors, ring_width)
ring_outer <- (1.05+ring_width)*tree_radius

# Normalized bile acid bars
bar_thickness <- 2
bar_total_len <- 0.15*tree_radius

for (i in seq_along(tree_life$tip.label)) {
  sp <- tree_life$tip.label[i]
  if (sp %in% names(species_totals)) {
    dat <- species_totals[[sp]]
    dat <- dat[order(factor(dat$core_std, levels=ba_levels)), ]
    x <- coords$xx[i]; y <- coords$yy[i]; angle <- atan2(y,x)
    r_start <- ring_outer
    dat$prop <- dat$total/sum(dat$total)
    for (j in seq_len(nrow(dat))) {
      seg_len <- dat$prop[j]*bar_total_len
      r_end <- r_start+seg_len
      segments(r_start*cos(angle), r_start*sin(angle),
               r_end*cos(angle), r_end*sin(angle),
               col = ba_colors[as.character(dat$core_std[j])],
               lwd=bar_thickness)
      r_start <- r_end
    }
  }
}

# Separator ring between bile acids and delta masses
separator_inner <- ring_outer+bar_total_len
separator_outer <- separator_inner+0.01*tree_radius
theta <- seq(0,2*pi,length.out=2000)
polygon(c(separator_inner*cos(theta), rev(separator_outer*cos(theta))),
        c(separator_inner*sin(theta), rev(separator_outer*sin(theta))),
        col="black", border=NA)

# Delta mass bars (log2 scale, colored gradient)
# t = log2(count + 1) for every species; bar length and color both driven by t
t_all <- log_species_totals_deltas  # log2(species_totals_deltas + 1), computed earlier
t_max <- max_log_val                # max(t_all, na.rm = TRUE), computed earlier

L_max <- 0.85 * tree_radius
L_min <- 0.02 * tree_radius
delta_scale <- (L_max - L_min) / t_max

for (i in seq_along(tree_life$tip.label)) {
  sp <- tree_life$tip.label[i]
  if (sp %in% names(species_totals_deltas)) {
    val <- species_totals_deltas[sp]
    t   <- log2(val + 1)
    
    x <- coords$xx[i]; y <- coords$yy[i]; angle <- atan2(y, x)
    r_start <- separator_outer
    bar_len <- L_min + t * delta_scale
    r_end   <- r_start + bar_len
    
    # color by same transform:
    col_idx <- round((t / t_max) * (length(delta_palette) - 1)) + 1
    segments(r_start*cos(angle), r_start*sin(angle),
             r_end*cos(angle),   r_end*sin(angle),
             col = delta_palette[col_idx], lwd = 3)
  }
}

# Legend for delta masses (same log2 scale as the bars)
image.plot(
  legend.only=TRUE,
  zlim=c(0,t_max),
  col=delta_palette,
  legend.width=1.5,
  legend.shrink=0.6,
  legend.mar=12,
  legend.args=list(text="log2(# Unique Delta Masses + 1)", side=4, line=6, cex=2),
  axis.args=list(cex.axis=1.8)
)

# Highlight top 10 species by delta masses that are actually in the tree
N <- 10  # how many species to label

# Restrict to overlap
overlap_species <- intersect(names(species_totals_deltas), tree_life$tip.label)

# Rank only within overlap
overlap_deltas <- species_totals_deltas[overlap_species]
top_overlap_species <- names(sort(overlap_deltas, decreasing = TRUE))[1:N]

cat("Top", N, "species by delta masses (in tree):\n")
print(top_overlap_species)

# Add labels at bar tips
for (sp in top_overlap_species) {
  i <- which(tree_life$tip.label == sp)
  if (length(i) == 0) next
  
  x <- coords$xx[i]; y <- coords$yy[i]; angle <- atan2(y, x)
  val <- species_totals_deltas[sp]
  t   <- log2(val + 1)
  bar_len <- L_min + t * delta_scale
  r_end   <- separator_outer + bar_len
  
  text(r_end * cos(angle) * 1.05,   # slightly outside the bar tip
       r_end * sin(angle) * 1.05,
       labels = sp,
       cex = 1.2, font = 3)
}

# Titles & legends
mtext("Metazoa Tree with Phylum Ring, C27 Bile Acids & Unique Delta Masses",
      side=3, outer=TRUE, line=2, cex=3.5)

legend("bottomleft", legend=names(phylum_colors),
       fill=phylum_colors, border="black", cex=2.5, bty="n", title="Phylum")

legend("bottomright", legend=names(ba_colors),
       fill=ba_colors, border="black", cex=2.5, bty="n", title="Bile Acids")

dev.off()