#!/usr/bin/env Rscript

# =============================================================================
# SUPPLEMENTARY FIGURES 6–7
# Dated circular phylogenies with inversion-trait annotation rings
# Required R packages:
#   ape, dplyr, readr, scales
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20261022)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit) == 0) {
    return(normalizePath(getwd(), mustWork = TRUE))
  }
  normalizePath(dirname(sub("^--file=", "", hit[[1]])), mustWork = TRUE)
}

ROOT <- get_script_dir()
setwd(ROOT)

FILES <- list(
  per_genome = "per_genome.csv",
  per_inversion = "per_inversion.csv",
  animal = "animal_mcmctree_augmented_verified.nwk",
  plant =  "plant_mcmctree_sanitized.nwk")

OUTDIR <- file.path(ROOT, "outputs", "supplementary")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# IMPORTANT:
# These MCMCTree trees are treated as:
#   1 branch-length unit = 100 million years.
# Keep this conversion explicit rather than embedding it inside plotting code.
TIME_MULTIPLIER_MA <- 100
TIME_SCALE_MA <- 100

# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------

  library(ape)
  library(dplyr)
  library(readr)
  library(scales)

# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

require_columns <- function(df, cols, label) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(
      label,
      " is missing required column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
}

prune_and_order <- function(tree, dat) {

  if (!"tip" %in% names(dat)) {
    stop("Analysis data must contain column 'tip'.", call. = FALSE)
  }

  if (anyNA(dat$tip) || any(dat$tip == "")) {
    stop("Analysis data contain empty/NA tree-tip labels.", call. = FALSE)
  }

  if (anyDuplicated(dat$tip)) {
    stop("Duplicate tip labels in analysis data.", call. = FALSE)
  }

  keep <- intersect(tree$tip.label, dat$tip)

  if (length(keep) < 4) {
    stop(
      "Fewer than four matched tips remain after tree/data matching.",
      call. = FALSE
    )
  }

  tr <- ape::keep.tip(tree, keep)

  d <- dat[
    match(tr$tip.label, dat$tip),
    ,
    drop = FALSE
  ]

  if (anyNA(d$tip)) {
    stop("Tree/data matching failure after pruning.", call. = FALSE)
  }

  if (!identical(tr$tip.label, d$tip)) {
    stop("Tree/data ordering failure after pruning.", call. = FALSE)
  }

  rownames(d) <- d$tip

  list(tree = tr, data = d)
}

# -----------------------------------------------------------------------------
# Read and audit dated trees
# -----------------------------------------------------------------------------

read_tree_strict <- function(path, label) {

  tr <- ape::read.tree(file = path)

  if (!inherits(tr, "phylo")) {
    stop("Cannot read tree: ", label, call. = FALSE)
  }

  tr$tip.label <- trimws(tr$tip.label)

  if (length(tr$tip.label) < 4) {
    stop("Tree has <4 tips: ", label, call. = FALSE)
  }

  if (anyNA(tr$tip.label) || any(tr$tip.label == "")) {
    stop("Empty/NA tip label: ", label, call. = FALSE)
  }

  if (anyDuplicated(tr$tip.label)) {
    stop("Duplicated tip label: ", label, call. = FALSE)
  }

  if (is.null(tr$edge.length)) {
    stop("Missing branch lengths: ", label, call. = FALSE)
  }

  if (any(!is.finite(tr$edge.length))) {
    stop("Non-finite branch length: ", label, call. = FALSE)
  }

  if (any(tr$edge.length < -1e-10)) {
    stop("Negative branch length: ", label, call. = FALSE)
  }

  # MCMCTree trees should be very close to ultrametric.
  dep <- ape::node.depth.edgelength(tr)[seq_len(ape::Ntip(tr))]
  rel_spread <- (max(dep) - min(dep)) / max(dep)

  if (rel_spread > 1e-4) {
    warning(
      label,
      " is not closely ultrametric (relative root-to-tip spread = ",
      signif(rel_spread, 4),
      ")."
    )
  }

  tr
}

tree_audit <- function(tr, label, time_multiplier_Ma = TIME_MULTIPLIER_MA) {

  dep <- ape::node.depth.edgelength(tr)[seq_len(ape::Ntip(tr))]

  tibble::tibble(
    tree = label,
    n_tips = ape::Ntip(tr),
    n_nodes = ape::Nnode(tr),
    duplicate_tips = sum(duplicated(tr$tip.label)),
    negative_edges = sum(tr$edge.length < -1e-10),
    zero_edges = sum(abs(tr$edge.length) <= 1e-12),
    root_to_tip_min_raw = min(dep),
    root_to_tip_max_raw = max(dep),
    root_to_tip_spread_raw = max(dep) - min(dep),
    relative_spread = (max(dep) - min(dep)) / max(dep),
    time_multiplier_Ma = time_multiplier_Ma,
    inferred_root_age_Ma = max(dep) * time_multiplier_Ma
  )
}

# -----------------------------------------------------------------------------
# Tree ↔ genome mapping
# -----------------------------------------------------------------------------

map_tree_to_genomes <- function(tree, genome_df, kingdom_name) {

  dat <- genome_df %>%
    filter(.data$kingdom == kingdom_name)

  tips <- tree$tip.label

  # Priority 1: exact sample_id match.
  exact_idx <- match(tips, dat$sample_id)
  exact_ok <- !is.na(exact_idx)

  exact <- tibble::tibble(
    tip = tips[exact_ok],
    row_id = dat$.row_id[exact_idx[exact_ok]],
    mapping_method = "sample_id_exact"
  )

  # Priority 2: exact species match only when that species occurs once
  # in the empirical table and the row has not already been used.
  used_rows <- exact$row_id
  remaining_tips <- setdiff(tips, exact$tip)

  sp_counts <- dat %>%
    count(.data$species_tip, name = "n_rows")

  unique_sp <- dat %>%
    left_join(sp_counts, by = "species_tip") %>%
    filter(
      .data$n_rows == 1,
      !(.data$.row_id %in% used_rows)
    )

  idx <- match(remaining_tips, unique_sp$species_tip)
  ok <- !is.na(idx)

  fallback <- tibble::tibble(
    tip = remaining_tips[ok],
    row_id = unique_sp$.row_id[idx[ok]],
    mapping_method = "unique_species_fallback"
  )

  mp <- bind_rows(exact, fallback)

  if (anyDuplicated(mp$tip)) {
    stop("A tree tip was mapped more than once.", call. = FALSE)
  }

  if (anyDuplicated(mp$row_id)) {
    stop("A genome row was mapped to >1 tree tip.", call. = FALSE)
  }

  mapped <- mp %>%
    left_join(dat, by = c("row_id" = ".row_id"))

  unmapped_tree_tips <- setdiff(tips, mapped$tip)

  unmapped_genomes <- dat %>%
    filter(!(.data$.row_id %in% mapped$row_id))

  list(
    mapped = mapped,
    unmapped_tree_tips = unmapped_tree_tips,
    unmapped_genomes = unmapped_genomes
  )
}

# -----------------------------------------------------------------------------
# Ring-colour helpers
# -----------------------------------------------------------------------------

continuous_cols <- function(x, palette = c("white", "black")) {

  x <- suppressWarnings(as.numeric(x))

  ok <- is.finite(x)

  out <- rep("#D9D9D9", length(x))

  if (sum(ok) <= 1) {
    return(out)
  }

  z <- scales::rescale(x[ok], to = c(0, 1))

  pal <- grDevices::colorRampPalette(palette)(101)

  out[ok] <- pal[
    pmax(
      1,
      pmin(101, 1 + floor(z * 100))
    )
  ]

  out
}

clade_palette <- function(x) {

  x <- as.character(x)

  lev <- sort(
    unique(
      x[!is.na(x) & trimws(x) != ""]
    )
  )

  n <- length(lev)

  if (n == 0) {
    return(character(0))
  }

  cols <- if (n <= 8) {
    grDevices::hcl.colors(n, "Dark 3")
  } else {
    grDevices::hcl.colors(n, "Dynamic")
  }

  setNames(cols, lev)
}

# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

pg <- readr::read_csv(
  FILES$per_genome,
  show_col_types = FALSE,
  progress = FALSE
)

require_columns(
  pg,
  c(
    "sample_id",
    "Species",
    "kingdom",
    "major_clade"
  ),
  "per_genome.csv"
)

if (anyDuplicated(pg$sample_id)) {
  stop(
    "per_genome.csv contains duplicated sample_id values; ",
    "tree mapping would be ambiguous.",
    call. = FALSE
  )
}

# Build the variables that the tree-mapping and ring-plotting functions require.
pg <- pg %>%
  mutate(
    .row_id = row_number(),
    species_tip = gsub(
      "[[:space:]]+",
      "_",
      trimws(.data$Species)
    )
  )

# Inversion count/Gb: use existing transformed column if present;
# otherwise derive it from inv_count_per_Gb.
if (!"log10_inv_count_per_Gb" %in% names(pg)) {

  if (!"inv_count_per_Gb" %in% names(pg)) {
    stop(
      "per_genome.csv must contain either 'log10_inv_count_per_Gb' ",
      "or 'inv_count_per_Gb'.",
      call. = FALSE
    )
  }

  pg <- pg %>%
    mutate(
      log10_inv_count_per_Gb = if_else(
        is.finite(as.numeric(.data$inv_count_per_Gb)) &
          as.numeric(.data$inv_count_per_Gb) > 0,
        log10(as.numeric(.data$inv_count_per_Gb)),
        NA_real_
      )
    )
}

# Median inversion length: the full phylogenetic script uses
# median_log10_inv_length. If it is absent, derive it from a recognized
# raw median-length column.
if (!"median_log10_inv_length" %in% names(pg)) {

  raw_length_candidates <- c(
    "median_inv_length_bp",
    "median_inversion_length_bp",
    "median_inv_length",
    "inv_length_median_bp"
  )

  hit <- raw_length_candidates[
    raw_length_candidates %in% names(pg)
  ]

  if (length(hit) == 0) {
    stop(
      "per_genome.csv is missing 'median_log10_inv_length' and no recognized ",
      "raw median-inversion-length column was found.",
      call. = FALSE
    )
  }

  raw_col <- hit[[1]]

  x <- suppressWarnings(as.numeric(pg[[raw_col]]))

  pg$median_log10_inv_length <- ifelse(
    is.finite(x) & x > 0,
    log10(x),
    NA_real_
  )
}

require_columns(
  pg,
  c(
    ".row_id",
    "species_tip",
    "log10_inv_count_per_Gb",
    "median_log10_inv_length"
  ),
  "prepared per_genome data"
)

# -----------------------------------------------------------------------------
# Read trees
# -----------------------------------------------------------------------------

animal_tree <- read_tree_strict(
  FILES$animal,
  "Metazoa"
)

plant_tree <- read_tree_strict(
  FILES$plant,
  "Viridiplantae"
)

tree_qc <- bind_rows(
  tree_audit(animal_tree, "Metazoa"),
  tree_audit(plant_tree, "Viridiplantae")
)

readr::write_csv(
  tree_qc,
  file.path(
    OUTDIR,
    "SuppFig6-7_tree_QC.csv"
  )
)

# -----------------------------------------------------------------------------
# Map trees to empirical genomes
# -----------------------------------------------------------------------------

animal_map <- map_tree_to_genomes(
  animal_tree,
  pg,
  "Metazoa"
)

plant_map <- map_tree_to_genomes(
  plant_tree,
  pg,
  "Viridiplantae"
)

animal_dat <- animal_map$mapped
plant_dat <- plant_map$mapped

coverage_tbl <- tibble::tibble(
  kingdom = c("Metazoa", "Viridiplantae"),
  empirical_genomes = c(
    sum(pg$kingdom == "Metazoa", na.rm = TRUE),
    sum(pg$kingdom == "Viridiplantae", na.rm = TRUE)
  ),
  tree_tips = c(
    ape::Ntip(animal_tree),
    ape::Ntip(plant_tree)
  ),
  mapped_genomes = c(
    nrow(animal_dat),
    nrow(plant_dat)
  ),
  unmapped_tree_tips = c(
    length(animal_map$unmapped_tree_tips),
    length(plant_map$unmapped_tree_tips)
  )
) %>%
  mutate(
    mapping_coverage_empirical =
      mapped_genomes / empirical_genomes
  )

readr::write_csv(
  coverage_tbl,
  file.path(
    OUTDIR,
    "SuppFig6-7_tree_mapping_coverage.csv"
  )
)

mapping_detail <- bind_rows(
  animal_dat %>%
    transmute(
      kingdom = "Metazoa",
      tip,
      sample_id,
      Species,
      major_clade,
      mapping_method
    ),

  plant_dat %>%
    transmute(
      kingdom = "Viridiplantae",
      tip,
      sample_id,
      Species,
      major_clade,
      mapping_method
    )
)

readr::write_csv(
  mapping_detail,
  file.path(
    OUTDIR,
    "SuppFig6-7_tree_mapping_detail.csv"
  )
)

# -----------------------------------------------------------------------------
# Publication-grade circular dated tree
# -----------------------------------------------------------------------------

plot_tree_with_rings <- function(
    tree,
    dat,
    kingdom_name,
    out_pdf,
    time_multiplier_Ma = TIME_MULTIPLIER_MA,
    scale_Ma = TIME_SCALE_MA,
    show_title = TRUE
) {

  # ---------------------------------------------------------------------------
  # 1. Match/prune tree and data
  # ---------------------------------------------------------------------------

  po <- prune_and_order(tree, dat)

  tr <- po$tree
  d <- po$data

  ntip <- ape::Ntip(tr)

  # ---------------------------------------------------------------------------
  # 2. Verify dated-tree geometry
  # ---------------------------------------------------------------------------

  tip_depth <- ape::node.depth.edgelength(tr)[
    seq_len(ntip)
  ]

  tree_height_raw <- max(
    tip_depth,
    na.rm = TRUE
  )

  root_age_Ma <- tree_height_raw *
    time_multiplier_Ma

  rel_spread <- (
    max(tip_depth, na.rm = TRUE) -
      min(tip_depth, na.rm = TRUE)
  ) / tree_height_raw

  if (rel_spread > 1e-4) {
    warning(
      kingdom_name,
      " pruned tree is not closely ultrametric: relative spread = ",
      signif(rel_spread, 4)
    )
  }

  if (
    !is.finite(scale_Ma) ||
      scale_Ma <= 0 ||
      scale_Ma >= root_age_Ma
  ) {
    stop(
      "scale_Ma must be >0 and smaller than the root age (",
      round(root_age_Ma, 1),
      " Ma).",
      call. = FALSE
    )
  }

  

  # ---------------------------------------------------------------------------
  # 3. Prepare ring colours
  # ---------------------------------------------------------------------------

  # Ring 1: major clade
  clade_cols <- clade_palette(
    d$major_clade
  )

  clade_value <- as.character(
    d$major_clade
  )

  c1 <- rep(
    "#D9D9D9",
    nrow(d)
  )

  clade_ok <- (
    !is.na(clade_value) &
      trimws(clade_value) != "" &
      clade_value %in% names(clade_cols)
  )

  c1[clade_ok] <- unname(
    clade_cols[
      clade_value[clade_ok]
    ]
  )

  # Ring 2: inversion count/Gb
  c2 <- continuous_cols(
    d$log10_inv_count_per_Gb,
    palette = c(
      "white",
      "black"
    )
  )

  # Ring 3: median inversion length
  c3 <- continuous_cols(
    d$median_log10_inv_length,
    palette = c(
      "white",
      "black"
    )
  )

  # ---------------------------------------------------------------------------
  # 4. Open vector PDF
  # ---------------------------------------------------------------------------

  grDevices::pdf(
    out_pdf,
    width = 10.5,
    height = 10.5,
    useDingbats = FALSE,
    family = "sans"
  )

  oldpar <- par(
    no.readonly = TRUE
  )

  on.exit(
    {
      par(oldpar)
      grDevices::dev.off()
    },
    add = TRUE
  )

  par(
    mar = c(
      0.35,
      0.35,
      1.45,
      0.35
    ),
    xpd = NA,
    lend = "butt"
  )

  # Leave enough room outside the tree for three rings,
  # title, scale bar and legends.
  plot_lim <- 1.36 *
    tree_height_raw

  # ---------------------------------------------------------------------------
  # 5. Plot dated fan phylogeny
  # ---------------------------------------------------------------------------

  ape::plot.phylo(
    tr,
    type = "fan",
    show.tip.label = FALSE,
    edge.color = "grey58",
    edge.width = ifelse(ntip > 500, 0.22, 0.30),
    no.margin = FALSE,
    label.offset = 0,
    x.lim = c(-plot_lim, plot_lim),
    y.lim = c(-plot_lim, plot_lim)
  )

  # Retrieve the exact coordinates used by ape.
  plot_env <- get(
    ".PlotPhyloEnv",
    envir = asNamespace("ape")
  )

  lp <- get(
    "last_plot.phylo",
    envir = plot_env
  )

  xx <- lp$xx[
    seq_len(ntip)
  ]

  yy <- lp$yy[
    seq_len(ntip)
  ]

  rr <- sqrt(
    xx^2 + yy^2
  )

  good <- (
    is.finite(rr) &
      rr > 0
  )

  if (!all(good)) {
    warning(
      kingdom_name,
      ": some plotted tip coordinates were non-finite."
    )
  }

  ux <- rep(
    NA_real_,
    ntip
  )

  uy <- rep(
    NA_real_,
    ntip
  )

  ux[good] <- xx[good] /
    rr[good]

  uy[good] <- yy[good] /
    rr[good]

  # For an ultrametric fan tree all tips should lie on nearly the same radius.
  r0 <- median(
    rr[good],
    na.rm = TRUE
  )

  # ---------------------------------------------------------------------------
  # 6. Add annotation rings
  # ---------------------------------------------------------------------------

  ring_r <- r0 * c(
    1.035,
    1.078,
    1.121
  )

  ring_cols <- list(
    c1,
    c2,
    c3
  )

  ring_cex <- if (
    ntip > 700
  ) {
    0.29
  } else if (
    ntip > 300
  ) {
    0.36
  } else {
    0.48
  }

  # Thin guide circles delimit the three annotation tracks.
  for (r in ring_r) {

    graphics::symbols(
      0,
      0,
      circles = r,
      inches = FALSE,
      add = TRUE,
      fg = "grey90",
      bg = NA,
      lwd = 0.22
    )
  }

  # Draw annotation cells.
  for (j in seq_along(ring_r)) {

    graphics::points(
      ux * ring_r[j],
      uy * ring_r[j],
      pch = 15,
      cex = ring_cex,
      col = ring_cols[[j]]
    )
  }

  # ---------------------------------------------------------------------------
  # 7. True branch-length-calibrated time scale
  # ---------------------------------------------------------------------------

  # Convert biological time into the coordinate system actually used by
  # plot.phylo. This remains correct even if the plotting region is expanded.
  plot_units_per_tree_unit <-
    r0 / tree_height_raw

  scale_tree_units <-
    scale_Ma / time_multiplier_Ma

  scale_plot_length <-
    scale_tree_units *
    plot_units_per_tree_unit

  # Keep the scale outside the outer ring but away from the corner legends.
  scale_y <- -1.235 * r0

  scale_x0 <- -scale_plot_length / 2
  scale_x1 <-  scale_plot_length / 2

  graphics::segments(
    x0 = scale_x0,
    y0 = scale_y,
    x1 = scale_x1,
    y1 = scale_y,
    lwd = 1.15,
    col = "black"
  )

  tick_height <- 0.016 *
    r0

  graphics::segments(
    x0 = scale_x0,
    y0 = scale_y - tick_height,
    x1 = scale_x0,
    y1 = scale_y + tick_height,
    lwd = 0.9,
    col = "black"
  )

  graphics::segments(
    x0 = scale_x1,
    y0 = scale_y - tick_height,
    x1 = scale_x1,
    y1 = scale_y + tick_height,
    lwd = 0.9,
    col = "black"
  )

  graphics::text(
    x = 0,
    y = scale_y - 0.038 * r0,
    labels = paste0(
      scale_Ma,
      " Ma"
    ),
    cex = 0.67,
    col = "black"
  )

  
  # ---------------------------------------------------------------------------
  #  Annotation-track key
  # ---------------------------------------------------------------------------

  graphics::legend(
    "bottomleft",
    inset = c(
      0.012,
      0.020
    ),
    legend = c(
      "Inner ring: major clade",
      "Middle ring: inversion count/Gb",
      "Outer ring: median inversion length"
    ),
    pch = 15,
    col = c(
      "grey35",
      "grey45",
      "grey70"
    ),
    pt.cex = 0.70,
    bty = "n",
    cex = 0.56,
    y.intersp = 1.02
  )

  # ---------------------------------------------------------------------------
  # Major-clade colour key
  # ---------------------------------------------------------------------------

  n_clades <- length(
    clade_cols
  )

  if (n_clades > 0) {

    graphics::legend(
      "bottomright",
      inset = c(
        0.004,
        0.012
      ),
      legend = names(
        clade_cols
      ),
      pch = 15,
      col = unname(
        clade_cols
      ),
      pt.cex = 0.62,
      bty = "n",
      cex = ifelse(
        n_clades > 14,
        0.39,
        ifelse(
          n_clades > 9,
          0.46,
          0.53
        )
      ),
      ncol = ifelse(
        n_clades > 9,
        2,
        1
      ),
      y.intersp = 0.88,
      x.intersp = 0.48
    )
  }

  invisible(
    list(
      kingdom = kingdom_name,
      n_tips = ntip,
      tree_height_raw = tree_height_raw,
      time_multiplier_Ma = time_multiplier_Ma,
      root_age_Ma = root_age_Ma,
      relative_ultrametric_spread = rel_spread,
      scale_Ma = scale_Ma,
      out_pdf = out_pdf
    )
  )
}

# -----------------------------------------------------------------------------
# Generate Supplementary Figures 6 and 7
# -----------------------------------------------------------------------------

plot_tree_with_rings(
  animal_tree,
  animal_dat,
  "Metazoa",
  file.path(
    OUTDIR,
    "SuppFig6_Tree_Metazoa_FINAL.pdf"
  ),
  time_multiplier_Ma = TIME_MULTIPLIER_MA,
  scale_Ma = TIME_SCALE_MA
)

plot_tree_with_rings(
  plant_tree,
  plant_dat,
  "Viridiplantae",
  file.path(
    OUTDIR,
    "SuppFig7_Tree_Viridiplantae_FINAL.pdf"
  ),
  time_multiplier_Ma = TIME_MULTIPLIER_MA,
  scale_Ma = TIME_SCALE_MA
)


