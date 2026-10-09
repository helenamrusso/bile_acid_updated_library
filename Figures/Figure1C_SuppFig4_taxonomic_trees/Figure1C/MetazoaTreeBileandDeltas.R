### Metazoa Tree + Phyla Ring + Normalized Bile Acids + Log 2 Delta Masses + top 10 organisms printed 
### USE THIS ONE
library(readr)
library(tidyr)
library(ape)
library(taxize)
library(purrr)
library(dplyr)
library(fields)
library(extrafont)

# Load Metazoa tree (created in MakeMetazoaTree.R)

tree_life <- read.tree("C:/Users/megan/OneDrive/Desktop/DorresteinLab/TaxaTree/tree_of_life_metazoa_species.nw")
set.seed(123)
tree_life <- multi2di(tree_life)


# Get coordinates of tree
# set.seed(123)
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

# Phylum colors (your palette)
phylum_colors <- c(
  "Annelida"      = "#007A99",
  "Arthropoda"    = "#00A9CC",
  "Chordata"      = "#65EFFF",
  "Cnidaria"      = "#CCFDFF",
  "Enchinodermata"= "#F2DACD",
  "Mollusca"      = "#D8AF97",
  "Nematoda"      = "#662F00",
  "Porifera"      = "#331900"
)


# Prepare Bile Acids
ba_class <- read_tsv("C:/Users/megan/OneDrive/Desktop/DorresteinLab/TaxaTree/NCBI_group_BA_class.tsv")

ba_class_filtered <- ba_class %>%
  filter(!(NCBITaxonomy %in% c("missing value", "N/A")) & !is.na(NCBITaxonomy)) %>%
  separate(NCBITaxonomy, into = c("TaxaID", "OrganismName"),
           sep = "\\|", fill = "right", extra = "merge")

ba_class_filtered$OrganismName <- gsub(" ", "_", ba_class_filtered$OrganismName)

ba_levels <- c("C24_nonhydroxy","C24_monohydroxy","C24_dihydroxy",
               "C24_trihydroxy","C24_tetrahydroxy","C24_pentahydroxy")

ba_colors <- c(
  "C24_nonhydroxy"   = "blue",
  "C24_monohydroxy"  = "forestgreen",
  "C24_dihydroxy"    = "orange",
  "C24_trihydroxy"   = "yellow",
  "C24_tetrahydroxy" = "purple",
  "C24_pentahydroxy" = "darkgray"
)

ba_class_summary <- ba_class_filtered %>%
  group_by(OrganismName, BA_class_clean) %>%
  summarise(total = sum(count), .groups = "drop")

ba_class_summary$BA_class_clean <- factor(ba_class_summary$BA_class_clean,
                                          levels = ba_levels, ordered = TRUE)


# Apply threshold: keep organisms where *every* bile acid class > 3
# The line that applies this to the tree is commented out, so it should not affect the tree
# Summarize by organism × class
ba_class_check <- ba_class_summary %>%
  group_by(OrganismName, BA_class_clean) %>%
  summarise(total = sum(total), .groups = "drop")

# Find organisms that have *all* classes present and > 3
eligible_species <- ba_class_check %>%
  group_by(OrganismName) %>%
  filter(all(ba_levels %in% BA_class_clean)) %>%       # must have all BA classes
  summarise(all_above_3 = all(total > 3)) %>%          # every class > 3
  filter(all_above_3) %>%
  pull(OrganismName)

cat("Keeping", length(eligible_species), "species with all BA classes > 3\n")

# Filter the original summary table
#ba_class_summary <- ba_class_summary %>%
  #filter(OrganismName %in% eligible_species)



species_totals <- split(ba_class_summary, ba_class_summary$OrganismName)


# Load Delta Masses File
deltas <- read_tsv("C:/Users/megan/OneDrive/Desktop/DorresteinLab/TaxaTree/NCBI_group_deltas.tsv")

deltas_filtered <- deltas %>%
  filter(!(NCBITaxonomy %in% c("missing value", "N/A")) & !is.na(NCBITaxonomy)) %>%
  separate(NCBITaxonomy, into = c("TaxaID", "OrganismName"),
           sep = "\\|", fill = "right", extra = "merge")

deltas_filtered$OrganismName <- gsub(" ", "_", deltas_filtered$OrganismName)

deltas_summary <- deltas_filtered %>%
  group_by(OrganismName) %>%
  summarise(num_deltas = n_distinct(delta_mass_round), .groups = "drop")

species_totals_deltas <- setNames(deltas_summary$num_deltas, deltas_summary$OrganismName)

# log2 transformation of delta masses data
log_species_totals_deltas <- log2(species_totals_deltas + 1)
max_log_val <- max(log_species_totals_deltas, na.rm = TRUE)

# Palette for delta masses
x <- seq(0, 1, length.out = 100)^1.5   # exponent >1 compresses low end

delta_palette <-  colorRampPalette(brewer.pal(9, "Reds"))(100) #viridisLite::turbo(100)
#colorRampPalette(c("deeppink", "red", "darkred"))(100)


# Plot Metazoa tree

pdf("metazoa_tree_phylaring_BA_deltas_log2_vector_reds_text_threshold.pdf",  width = 30, height = 35,
    family = "Helvetica")

par(oma = c(1,1,6,6), family = "Helvetica")

# Fix limits
plot(tree_life, type = "fan", cex = 0.3, show.tip.label = FALSE,
     no.margin = TRUE,
     x.lim = c(-1.8,1.8)*tree_radius,
     y.lim = c(-1.8,1.8)*tree_radius)
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
    dat <- dat[order(factor(dat$BA_class_clean, levels=ba_levels)), ]
    x <- coords$xx[i]; y <- coords$yy[i]; angle <- atan2(y,x)
    r_start <- ring_outer
    dat$prop <- dat$total/sum(dat$total)
    for (j in seq_len(nrow(dat))) {
      seg_len <- dat$prop[j]*bar_total_len
      r_end <- r_start+seg_len
      segments(r_start*cos(angle), r_start*sin(angle),
               r_end*cos(angle), r_end*sin(angle),
               col = ba_colors[as.character(dat$BA_class_clean[j])],
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
# target overall max height and a small floor so tiny values are still visible
k <- as.numeric(stats::quantile(species_totals_deltas, 0.8, na.rm = TRUE)) # tweak 0.7–0.9
t_all <- asinh(species_totals_deltas / max(k, 1))
t_max <- max(t_all, na.rm = TRUE)

L_max <- 0.85 * tree_radius
L_min <- 0.02 * tree_radius
delta_scale <- (L_max - L_min) / t_max

for (i in seq_along(tree_life$tip.label)) {
  sp <- tree_life$tip.label[i]
  if (sp %in% names(species_totals_deltas)) {
    val <- species_totals_deltas[sp]
    t   <- asinh(val / max(k, 1))
    
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

# Legend for delta masses
image.plot(
  legend.only=TRUE,
  zlim=c(0,max_log_val),
  col=delta_palette,
  legend.width=1.5,
  legend.shrink=0.6,
  legend.mar=12,
  legend.args=list(text="# Unique Delta Masses Scale", side=4, line=6, cex=2),
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
  t   <- asinh(val / max(k, 1))
  bar_len <- L_min + t * delta_scale
  r_end   <- separator_outer + bar_len
  
  text(r_end * cos(angle) * 1.05,   # slightly outside the bar tip
       r_end * sin(angle) * 1.05,
       labels = sp,
       cex = 1.2, font = 3)
}

# Titles & legends
mtext("Metazoa Tree with Phylum Ring, C24 Bile Acids & Unique Delta Masses",
      side=3, outer=TRUE, line=2, cex=3.5)

legend("bottomleft", legend=names(phylum_colors),
       fill=phylum_colors, border="black", cex=2.5, bty="n", title="Phylum")

legend("bottomright", legend=names(ba_colors),
       fill=ba_colors, border="black", cex=2.5, bty="n", title="Bile Acids")

dev.off()

