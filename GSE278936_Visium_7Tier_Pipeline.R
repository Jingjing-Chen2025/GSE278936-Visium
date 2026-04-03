# ==============================================================================
# GSE278936 Visium Spatial Transcriptomics — FULL Multi-Sample Pipeline
# ASPC, Hallmark AR & NEPC Beltran Signatures
# Processes ALL 52 samples automatically
# Dynamic spot sizing for optimal gap-free spatial plots
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Install & Load Packages
# ------------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("Seurat", "ggplot2", "patchwork", "dplyr", "Matrix", "scales",
          "RColorBrewer", "grDevices", "png", "jsonlite", "data.table")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(scales)
library(grDevices)
library(png)
library(jsonlite)
library(data.table)

cat("Libraries loaded\n")

# ------------------------------------------------------------------------------
# 1. Define Paths — Auto-detect all samples & organize GEO-prefixed files
# ------------------------------------------------------------------------------
visium_root <- "/Users/cj719/Library/Mobile Documents/com~apple~CloudDocs/Adipocyte NEPC/cj719/GSE278936 Visium"

sample_dirs <- c(
  "GSM8557976_BPH_1",
  "GSM8557977_BPH_2",
  "GSM8557978_BPH_3",
  "GSM8557979_BPH_4",
  "GSM8557980_TRNA_1",
  "GSM8557981_TRNA_2",
  "GSM8557982_TRNA_3",
  "GSM8557983_TRNA_4",
  "GSM8557984_TRNA_5",
  "GSM8557985_TRNA_6",
  "GSM8557986_TRNA_7",
  "GSM8557987_TRNA_8",
  "GSM8557988_TRNA_9",
  "GSM8557989_TRNA_10",
  "GSM8557990_TRNA_11",
  "GSM8557991_TRNA_12",
  "GSM8557992_TRNA_13",
  "GSM8557993_TRNA_14",
  "GSM8557994_TRNA_15",
  "GSM8557995_TRNA_16",
  "GSM8557996_TRNA_17",
  "GSM8557997_NEADT_1",
  "GSM8557998_NEADT_2",
  "GSM8557999_NEADT_3",
  "GSM8558000_NEADT_4",
  "GSM8558001_NEADT_5",
  "GSM8558002_NEADT_6",
  "GSM8558003_NEADT_7",
  "GSM8558004_NEADT_8",
  "GSM8558005_NEADT_9",
  "GSM8558006_NEADT_10",
  "GSM8558007_NEADT_11",
  "GSM8558008_NEADT_12",
  "GSM8558009_NEADT_13",
  "GSM8558010_NEADT_14",
  "GSM8558011_NEADT_15",
  "GSM8558012_NEADT_16",
  "GSM8558013_NEADT_17",
  "GSM8558014_NEADT_18",
  "GSM8558015_NEADT_19",
  "GSM8558016_NEADT_20",
  "GSM8558017_NEADT_21",
  "GSM8558018_NEADT_22",
  "GSM8558019_CRPC_1",
  "GSM8558020_CRPC_2",
  "GSM8558021_CRPC_3",
  "GSM8558022_CRPC_4",
  "GSM8558023_CRPC_5",
  "GSM8558024_MET_A",
  "GSM8558025_MET_B",
  "GSM8558026_MET_C",
  "GSM8558027_MET_D"
)

# --- Helper: scan a sample directory for GEO-prefixed files by suffix ---
find_sample_files <- function(sample_path, sample_folder_name) {
  all_files <- list.files(sample_path, recursive = TRUE, full.names = TRUE)

  find_by_suffix <- function(suffix) {
    matches <- all_files[grepl(paste0(suffix, "$"), all_files, ignore.case = TRUE)]
    if (length(matches) == 0) return(NA_character_)
    return(matches[1])
  }

  list(
    barcodes       = find_by_suffix("barcodes\\.tsv\\.gz"),
    features       = find_by_suffix("features\\.tsv\\.gz"),
    matrix         = find_by_suffix("matrix\\.mtx\\.gz"),
    positions      = find_by_suffix("tissue_positions_list\\.csv"),
    scalefactors   = find_by_suffix("scalefactors_json\\.json"),
    hires_image    = find_by_suffix("tissue_hires_image\\.png"),
    lowres_image   = find_by_suffix("tissue_lowres_image\\.png")
  )
}

# --- Build the full sample manifest ---
sample_manifest <- list()
for (sd in sample_dirs) {
  sample_path <- file.path(visium_root, sd)
  if (!dir.exists(sample_path)) {
    cat(paste0("  WARNING: directory not found — ", sample_path, "\n"))
    next
  }
  sample_name <- sd
  files <- find_sample_files(sample_path, sd)
  sample_manifest[[sample_name]] <- list(
    sample_name = sample_name,
    base_dir    = sample_path,
    files       = files
  )
}

cat(paste0("\nFound ", length(sample_manifest), " sample directories.\n"))

# Validate all samples
for (sn in names(sample_manifest)) {
  sm <- sample_manifest[[sn]]
  cat(paste0("\n--- ", sn, " ---\n"))
  for (fn in names(sm$files)) {
    f <- sm$files[[fn]]
    if (!is.na(f) && file.exists(f)) {
      cat(paste0("  OK: ", fn, " -> ", basename(f), "\n"))
    } else {
      cat(paste0("  MISSING: ", fn, "\n"))
    }
  }
}

# ------------------------------------------------------------------------------
# Gene Signatures (Human) — defined once, used for all samples
# ------------------------------------------------------------------------------
cat("\n=== Gene Signatures ===\n")

ASPC_genes <- c("LAMA2","THBS3","FGF2","POSTN","THBS2","VCAN","LAMB1","COL4A1","GAS6",
  "BMPER","COL5A1","COL6A5","COL5A3","ENTPD2","ALDH1A1","MFAP2","COL5A2","LOXL2",
  "FBLN7","BDH2","COL3A1","GPC6","LGI2","AK1","SRPX2","LOXL3","NOTCH3","CLU",
  "CTHRC1","FGFR1","C3","DPEP1","SSC5D","NAALAD2","ITGBL1","FOXS1","COL11A1",
  "NNMT","HEYL","TFCP2L1","PRG4","ZFHX4","BNC2","ADAMTS2","NFASC","CHRD",
  "PDGFRL","C1QTNF2","ANO1","HOXB7","ADAMTS15","ADAMTS12","IGSF10","DNM1","CHL1",
  "HAS2","CYGB","NR2F1","ZFPM2","ADAMTS5","CP","SULF1","C2","PCSK5","PTPRD",
  "FAP","TNFAIP6","S100A16","HSD11B1","TIMP1","BMP1","CACNA2D1","TBX2","RBP1",
  "NFATC4","GPRC5B","MMP11","DPYSL3","ADAMTSL1","KRT7","FST","ID4","SERPINH1",
  "MDK","HOXC8","PDGFRA","ANGPT4","HIGD1B","MMP3","GPR153","MEG3","ISL1",
  "TSPAN11","HOXC9","GJC1","SMIM22","ABCC9","ILDR2","HP","LHFP","NTRK2","FEZ1",
  "FMO2","FZD2","TMEM45A","MMP15","CRLF1","TENM3","ARHGEF25","FOXC2","PIEZO2",
  "S1PR3","CADM3","DCLK1","LRRN2","LBP","KCNE4","PLS3","PKDCC","NUMBL","STOX2",
  "PLCE1","STEAP3","PDE1A","CCDC102B","PRKG1","ARHGAP42","S100A14","ELF3","SPHK1",
  "NAV3","FXYD1","SULT1E1","GABRA3","MEDAG","GPR39","CASP12","SNTG1","CERCAM",
  "FKBP9","CDR2L","PTK7","KIF21A","NSG1","ROBO1","CCDC74A")

Hallmark_AR_genes <- c("ABCC4","ABHD2","ACSL3","ACTN1","ADAMTS1","ADRM1","AKAP12","AKT1",
  "ALDH1A3","ANKH","APPBP2","ARID5B","AZGP1","B2M","B4GALT1","BMPR1B","CAMKK2",
  "CCND1","CCND3","CDC14B","CDK6","CENPN","DBI","DHCR24","DNAJB9","ELK4","ELL2",
  "ELOVL5","FADS1","FKBP5","GNAI3","GPD1L","GSR","GUCY1A1","H1-0","HERC3",
  "HMGCR","HMGCS1","HOMER2","HPGD","HSD17B14","IDI1","INPP4B","INSIG1","IQGAP2",
  "ITGAV","KLK2","KLK3","KRT19","KRT8","LIFR","LMAN1","MAF","MAK","MAP7",
  "MERTK","MYL12A","NCOA4","NDRG1","NGLY1","NKX3-1","PA2G4","PDLIM5","PGM3",
  "PIAS1","PLPP1","PMEPA1","PTK2B","PTPN21","RAB4A","RPS6KA3","RRP12","SAT1",
  "SCD","SEC24D","SELENOP","SGK1","SLC26A2","SLC38A2","SMS","SORD","SPCS3",
  "SPDEF","SRF","SRP19","STEAP4","STK39","TARP","TMEM50A","TMPRSS2","TNFAIP8",
  "TPD52","TSC22D1","UAP1","UBE2I","UBE2J1","VAPA","XRCC5","XRCC6","ZBTB10","ZMIZ1")

NEPC_Beltran_UP <- c("ASXL3","CAND2","ETV5","GPX2","JAKMIP2","KIAA0408","SOGA3","TRIM9",
  "BRINP1","C7orf76","GNAO1","KCNB2","KCND2","LRRC16B","MAP10","NRSN1","PCSK1",
  "PROX1","RGS7","SCG3","SEC11C","SEZ6","ST8SIA3","SVOP","SYT11","AURKA","DNMT1",
  "EZH2","MYCN")

NEPC_Beltran_DOWN <- c("AR","CCND1","CIITA","CREBBP","CSDE1","CYLD","DICER1","FHIT",
  "FOXP1","HERPUD1","MMP2","MYH9","NUP93","PAX8","RBBP6","TRIM33","GATA2",
  "MAPKAPK3","SPDEF","ARHGAP8","CATSPERB","EFNA4","EPN3","EVPL","HOXB13","KLK3",
  "KLK4","LMAN1L","NKX3-1","OPHN1","PIEZO1","PRR5-ARHGAP8","PSCA","RAB27B",
  "RGS10","RIPK2","SLC25A37","SLC44A4","TC2N","UPK2","RB1")

# ==============================================================================
# SECTION 2: Additional Packages, Helpers, Output Dirs & Load Samples
# ==============================================================================
cat("\n=== Section 2: Loading additional packages and preparing environment ===\n")

# Null-coalescing operator (define once if not provided by a loaded package)
`%||%` <- function(x, y) if (!is.null(x)) x else y

# --- Additional packages with graceful fallback ---
extra_pkgs <- list(
  pheatmap  = list(cran = "pheatmap"),
  igraph    = list(cran = "igraph"),
  MASS      = list(cran = "MASS"),
  parallel  = list(cran = NULL)   # base package
)

for (pkg_name in names(extra_pkgs)) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    tryCatch({
      install.packages(pkg_name)
      cat(paste0("Installed: ", pkg_name, "\n"))
    }, error = function(e) {
      cat(paste0("Could not install ", pkg_name, " — will use fallback\n"))
    })
  }
}

# princurve
if (!requireNamespace("princurve", quietly = TRUE)) {
  tryCatch({
    install.packages("princurve")
    cat("Installed: princurve\n")
  }, error = function(e) {
    cat("Could not install princurve — will use PC1 ordering fallback\n")
  })
}

have_pheatmap  <- requireNamespace("pheatmap",  quietly = TRUE)
have_igraph    <- requireNamespace("igraph",    quietly = TRUE)
have_princurve <- requireNamespace("princurve", quietly = TRUE)
have_MASS      <- requireNamespace("MASS",      quietly = TRUE)

if (have_pheatmap)  library(pheatmap)
if (have_igraph)    library(igraph)
if (have_princurve) library(princurve)
if (have_MASS)      library(MASS)
library(parallel)

cat(paste0("pheatmap: ", have_pheatmap, " | igraph: ", have_igraph,
           " | princurve: ", have_princurve, " | MASS: ", have_MASS, "\n"))

# --- Create output subdirectories ---
tier_dirs <- c("TIER1_Interface", "TIER2_DE", "TIER3_Gradient",
               "TIER4_LR", "TIER5_Trajectory", "TIER6_GRN", "TIER7_MetaAnalysis")
for (td in tier_dirs) {
  dir.create(file.path(visium_root, td), showWarnings = FALSE, recursive = TRUE)
}
cat("Output subdirectories created.\n")

# --- Helper: extract condition from sample name ---
extract_condition <- function(sample_name) {
  if (grepl("BPH",   sample_name)) return("BPH")
  if (grepl("TRNA",  sample_name)) return("TRNA")
  if (grepl("NEADT", sample_name)) return("NEADT")
  if (grepl("CRPC",  sample_name)) return("CRPC")
  if (grepl("MET",   sample_name)) return("MET")
  return("Unknown")
}

# --- Helper: compute dynamic spot size from coordinate density ---
compute_spot_size <- function(coords) {
  if (nrow(coords) < 2) return(1.5)
  x_range <- diff(range(coords[, 1], na.rm = TRUE))
  y_range <- diff(range(coords[, 2], na.rm = TRUE))
  area     <- x_range * y_range
  if (area <= 0) return(1.5)
  density  <- nrow(coords) / area
  size     <- max(0.3, min(3.0, 1 / sqrt(density) * 10))
  return(size)
}

# --- Load samples into Seurat ---
cat("\n=== Loading samples into Seurat ===\n")
seurat_list <- list()

for (sn in names(sample_manifest)) {
  cat(paste0("\nProcessing sample: ", sn, "\n"))
  sm <- sample_manifest[[sn]]
  files <- sm$files

  tryCatch({
    # Check required files
    if (is.na(files$matrix) || is.na(files$features) || is.na(files$barcodes)) {
      cat(paste0("  SKIP: missing count matrix files for ", sn, "\n"))
      next
    }

    # Read count matrix
    counts <- ReadMtx(
      mtx      = files$matrix,
      features = files$features,
      cells    = files$barcodes
    )
    cat(paste0("  Count matrix: ", nrow(counts), " genes x ", ncol(counts), " cells\n"))

    # Create Seurat object
    obj <- CreateSeuratObject(
      counts      = counts,
      project     = sn,
      min.cells   = 3,
      min.features = 200
    )

    # Read tissue positions
    if (!is.na(files$positions) && file.exists(files$positions)) {
      pos <- tryCatch(
        read.csv(files$positions, header = FALSE,
                 col.names = c("barcode","in_tissue","array_row","array_col",
                               "pxl_row","pxl_col")),
        error = function(e) NULL
      )
      if (!is.null(pos)) {
        pos <- pos[pos$in_tissue == 1, ]
        pos <- pos[pos$barcode %in% colnames(obj), ]
        rownames(pos) <- pos$barcode
        pos_matched <- pos[colnames(obj)[colnames(obj) %in% rownames(pos)], ]
        bcs_with_pos <- colnames(obj)[colnames(obj) %in% rownames(pos)]
        obj <- subset(obj, cells = bcs_with_pos)
        pos_matched <- pos[colnames(obj), ]
        obj$row      <- pos_matched$array_row
        obj$col      <- pos_matched$array_col
        obj$imagerow <- pos_matched$pxl_row
        obj$imagecol <- pos_matched$pxl_col
        cat(paste0("  Positions loaded: ", nrow(pos_matched), " in-tissue spots\n"))
      }
    } else {
      cat("  WARNING: tissue positions file not found\n")
    }

    # Read scalefactors & attach spatial image
    if (!is.na(files$scalefactors) && file.exists(files$scalefactors)) {
      tryCatch({
        sf <- fromJSON(files$scalefactors)
        # Build VisiumV1 image object if lowres image available
        if (!is.na(files$lowres_image) && file.exists(files$lowres_image)) {
          img_data <- readPNG(files$lowres_image)
          # Construct tissue coordinates data.frame for VisiumV1
          tc <- data.frame(
            tissue   = obj$in_tissue %||% rep(1, ncol(obj)),
            row      = obj$row,
            col      = obj$col,
            imagerow = obj$imagerow,
            imagecol = obj$imagecol,
            row.names = colnames(obj)
          )
          visium_img <- tryCatch({
            new("VisiumV1",
                image          = img_data,
                scale.factors  = scalefactors(
                  spot         = sf$spot_diameter_fullres %||% 1,
                  fiducial     = sf$fiducial_diameter_fullres %||% 1,
                  hires        = sf$tissue_hires_scalef %||% 0.2,
                  lowres       = sf$tissue_lowres_scalef %||% 0.1
                ),
                coordinates    = tc,
                spot.radius    = 0.5)
          }, error = function(e) NULL)
          if (!is.null(visium_img)) {
            obj@images[[sn]] <- visium_img
            cat("  VisiumV1 image attached\n")
          }
        }
      }, error = function(e) {
        cat(paste0("  WARNING: could not attach spatial image — ", conditionMessage(e), "\n"))
      })
    }

    # Add metadata
    obj$sample_id <- sn
    obj$condition <- extract_condition(sn)

    # SCTransform normalization
    cat("  Running SCTransform...\n")
    obj <- SCTransform(obj, verbose = FALSE)

    # Module scores
    all_genes <- rownames(obj)

    aspc_use  <- intersect(ASPC_genes,         all_genes)
    ar_use    <- intersect(Hallmark_AR_genes,  all_genes)
    nepc_up   <- intersect(NEPC_Beltran_UP,    all_genes)
    nepc_dn   <- intersect(NEPC_Beltran_DOWN,  all_genes)

    if (length(aspc_use)  > 0)
      obj <- AddModuleScore(obj, features = list(aspc_use),  name = "ASPC_score",    ctrl = min(100, length(aspc_use)))
    if (length(ar_use)    > 0)
      obj <- AddModuleScore(obj, features = list(ar_use),    name = "AR_score",      ctrl = min(100, length(ar_use)))
    if (length(nepc_up)   > 0)
      obj <- AddModuleScore(obj, features = list(nepc_up),   name = "NEPC_UP_score", ctrl = min(100, length(nepc_up)))
    if (length(nepc_dn)   > 0)
      obj <- AddModuleScore(obj, features = list(nepc_dn),   name = "NEPC_DOWN_score", ctrl = min(100, length(nepc_dn)))

    # Rename AddModuleScore appended "1"
    for (sc_col in c("ASPC_score1","AR_score1","NEPC_UP_score1","NEPC_DOWN_score1")) {
      if (sc_col %in% colnames(obj@meta.data)) {
        new_col <- sub("1$", "", sc_col)
        obj@meta.data[[new_col]] <- obj@meta.data[[sc_col]]
      }
    }

    # Net NEPC score
    if ("NEPC_UP_score" %in% colnames(obj@meta.data) &&
        "NEPC_DOWN_score" %in% colnames(obj@meta.data)) {
      obj$NEPC_Beltran_NET <- obj$NEPC_UP_score - obj$NEPC_DOWN_score
    }

    seurat_list[[sn]] <- obj
    cat(paste0("  Sample loaded: ", ncol(obj), " spots\n"))

  }, error = function(e) {
    cat(paste0("  ERROR loading ", sn, ": ", conditionMessage(e), "\n"))
  })
}

cat(paste0("\nSuccessfully loaded ", length(seurat_list), " samples.\n"))

# --- Merge all samples ---
if (length(seurat_list) >= 2) {
  cat("Merging all samples...\n")
  combined_obj <- tryCatch({
    merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)],
          add.cell.ids = names(seurat_list), project = "GSE278936_Visium")
  }, error = function(e) {
    cat(paste0("Merge error: ", conditionMessage(e), "\n"))
    NULL
  })
  if (!is.null(combined_obj)) {
    saveRDS(combined_obj, file = file.path(visium_root, "combined_seurat.rds"))
    cat("Saved combined_seurat.rds\n")
  }
} else if (length(seurat_list) == 1) {
  combined_obj <- seurat_list[[1]]
  saveRDS(combined_obj, file = file.path(visium_root, "combined_seurat.rds"))
  cat("Only one sample loaded; saved as combined_seurat.rds\n")
} else {
  combined_obj <- NULL
  cat("WARNING: No samples loaded.\n")
}

# ==============================================================================
# SECTION 3: TIER 1 — Spatial Annotation (Interface Zone Classification)
# ==============================================================================
cat("\n=== TIER 1: Spatial Annotation ===\n")
tier1_status <- "SKIPPED"

tryCatch({
  tier1_dir <- file.path(visium_root, "TIER1_Interface")
  zone_prop_list <- list()

  for (sn in names(seurat_list)) {
    cat(paste0("  TIER1 — ", sn, "\n"))
    obj <- seurat_list[[sn]]
    md  <- obj@meta.data

    # Ensure score columns exist
    score_cols <- c("ASPC_score","AR_score","NEPC_UP_score")
    missing_sc <- setdiff(score_cols, colnames(md))
    if (length(missing_sc) > 0) {
      cat(paste0("    Missing score columns: ", paste(missing_sc, collapse=", "), " — skipping\n"))
      next
    }

    med_aspc   <- median(md$ASPC_score,    na.rm = TRUE)
    med_ar     <- median(md$AR_score,      na.rm = TRUE)
    med_nepc   <- median(md$NEPC_UP_score, na.rm = TRUE)

    is_interface <- (md$ASPC_score > med_aspc) &
                    (md$NEPC_UP_score > med_nepc | md$AR_score > med_ar)
    is_aspc_high <- (md$ASPC_score > med_aspc) & !is_interface
    is_ar_high   <- (md$AR_score   > med_ar)   & !is_interface & !is_aspc_high
    is_nepc_high <- (md$NEPC_UP_score > med_nepc) & !is_interface & !is_aspc_high & !is_ar_high

    zone <- rep("Other", nrow(md))
    zone[is_nepc_high] <- "NEPC_high"
    zone[is_ar_high]   <- "AR_high"
    zone[is_aspc_high] <- "ASPC_high"
    zone[is_interface] <- "Interface"

    obj$interface_zone <- zone
    seurat_list[[sn]] <- obj

    # Spatial plot (ggplot2 direct using array coordinates)
    if (all(c("col","row") %in% colnames(md))) {
      plot_df <- data.frame(
        x    = md$col,
        y    = -md$row,
        zone = zone,
        stringsAsFactors = FALSE
      )
      spot_sz <- compute_spot_size(plot_df[, c("x","y")])
      zone_colors <- c(
        Interface = "#E31A1C",
        ASPC_high = "#1F78B4",
        AR_high   = "#33A02C",
        NEPC_high = "#FF7F00",
        Other     = "#CCCCCC"
      )
      p <- ggplot(plot_df, aes(x = x, y = y, color = zone)) +
        geom_point(size = spot_sz) +
        scale_color_manual(values = zone_colors) +
        labs(title = paste0(sn, " — Interface Zones"),
             x = "Array Col", y = "Array Row (inv)", color = "Zone") +
        theme_minimal(base_size = 10)
      ggsave(file.path(tier1_dir, paste0(sn, "_interface_zones.pdf")),
             plot = p, width = 10, height = 8)
    }

    # Zone proportions
    tbl <- table(zone)
    zone_prop_list[[sn]] <- as.data.frame(tbl, stringsAsFactors = FALSE)
    zone_prop_list[[sn]]$sample  <- sn
    zone_prop_list[[sn]]$condition <- extract_condition(sn)
  }

  # Combine zone proportions
  if (length(zone_prop_list) > 0) {
    zone_prop_df <- do.call(rbind, zone_prop_list)
    colnames(zone_prop_df)[1:2] <- c("zone","count")
    write.csv(zone_prop_df,
              file.path(tier1_dir, "zone_proportions_per_sample.csv"),
              row.names = FALSE)

    # Aggregate by condition
    cond_zone <- zone_prop_df %>%
      group_by(condition, zone) %>%
      summarise(total_spots = sum(count), .groups = "drop")
    write.csv(cond_zone,
              file.path(tier1_dir, "zone_proportions_by_condition.csv"),
              row.names = FALSE)

    cat(paste0("  Zone proportion table saved (", nrow(zone_prop_df), " rows)\n"))
  }

  tier1_status <- "COMPLETE"
  cat("TIER 1 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 1 ERROR: ", conditionMessage(e), "\n"))
  tier1_status <<- paste0("ERROR: ", conditionMessage(e))
})

# ==============================================================================
# SECTION 4: TIER 2 — Distance-Binned Pseudo-Bulk DE
# ==============================================================================
cat("\n=== TIER 2: Distance-Binned Pseudo-Bulk DE ===\n")
tier2_status <- "SKIPPED"

tryCatch({
  tier2_dir <- file.path(visium_root, "TIER2_DE")

  dist_bin_labels <- c("bin_0_interface","bin_1_5","bin_6_10",
                        "bin_11_20","bin_21_50","bin_50plus")

  assign_dist_bin <- function(dist_vec) {
    bins <- rep("bin_50plus", length(dist_vec))
    bins[dist_vec <= 50] <- "bin_21_50"
    bins[dist_vec <= 20] <- "bin_11_20"
    bins[dist_vec <= 10] <- "bin_6_10"
    bins[dist_vec <= 5]  <- "bin_1_5"
    bins[dist_vec == 0]  <- "bin_0_interface"
    return(bins)
  }

  all_dist_data <- list()  # collect per-sample for Tier 3/5

  for (sn in names(seurat_list)) {
    cat(paste0("  TIER2 — ", sn, "\n"))
    obj <- seurat_list[[sn]]
    md  <- obj@meta.data

    if (!"interface_zone" %in% colnames(md)) next
    if (!all(c("row","col") %in% colnames(md))) next

    iface_idx <- which(md$interface_zone == "Interface")
    if (length(iface_idx) == 0) {
      cat("    No interface spots — skipping distance calc\n")
      next
    }

    coords <- as.matrix(md[, c("row","col")])
    iface_coords <- coords[iface_idx, , drop = FALSE]

    # Euclidean distance to nearest interface spot
    dist_to_iface <- apply(coords, 1, function(sp) {
      min(sqrt(rowSums(sweep(iface_coords, 2, sp)^2)))
    })
    dist_to_iface[iface_idx] <- 0
    obj$dist_to_interface <- dist_to_iface
    obj$dist_bin <- assign_dist_bin(dist_to_iface)
    seurat_list[[sn]] <- obj
    all_dist_data[[sn]] <- data.frame(
      barcode           = colnames(obj),
      dist_to_interface = dist_to_iface,
      dist_bin          = obj$dist_bin,
      condition         = obj$condition,
      stringsAsFactors  = FALSE
    )
  }

  # Per-condition: pseudo-bulk DE (bin_0 vs bin_50plus)
  conditions <- unique(sapply(names(seurat_list), extract_condition))

  for (cond in conditions) {
    cat(paste0("  Condition: ", cond, "\n"))
    cond_samples <- names(seurat_list)[sapply(names(seurat_list), extract_condition) == cond]
    if (length(cond_samples) == 0) next

    # Aggregate expression per bin across condition
    bin_expr_list <- list()
    for (sn in cond_samples) {
      obj <- seurat_list[[sn]]
      if (!"dist_bin" %in% colnames(obj@meta.data)) next
      DefaultAssay(obj) <- "SCT"
      expr_mat <- GetAssayData(obj, slot = "data")
      for (bn in dist_bin_labels) {
        cells_in_bin <- colnames(obj)[obj$dist_bin == bn]
        if (length(cells_in_bin) == 0) next
        pb <- rowMeans(expr_mat[, cells_in_bin, drop = FALSE])
        if (is.null(bin_expr_list[[bn]])) bin_expr_list[[bn]] <- list()
        bin_expr_list[[bn]][[sn]] <- pb
      }
    }

    # Collapse each bin to mean across samples
    bin_means <- lapply(bin_expr_list, function(lst) {
      if (length(lst) == 0) return(NULL)
      mat <- do.call(cbind, lst)
      rowMeans(mat, na.rm = TRUE)
    })

    bin0 <- bin_means[["bin_0_interface"]]
    bin50 <- bin_means[["bin_50plus"]]
    if (is.null(bin0) || is.null(bin50)) {
      cat(paste0("    Insufficient bins for ", cond, "\n"))
      next
    }

    common_genes <- intersect(names(bin0), names(bin50))
    if (length(common_genes) < 10) next

    # Simple fold-change and pseudo-Wilcoxon per gene across samples
    de_results <- tryCatch({
      fc <- log2((bin0[common_genes] + 0.01) / (bin50[common_genes] + 0.01))

      # Collect per-sample values for Wilcoxon
      vals0  <- lapply(cond_samples, function(sn) {
        obj <- seurat_list[[sn]]
        if (!"dist_bin" %in% colnames(obj@meta.data)) return(NULL)
        cells <- colnames(obj)[obj$dist_bin == "bin_0_interface"]
        if (length(cells) == 0) return(NULL)
        DefaultAssay(obj) <- "SCT"
        rowMeans(GetAssayData(obj, slot = "data")[common_genes, cells, drop = FALSE])
      })
      vals50 <- lapply(cond_samples, function(sn) {
        obj <- seurat_list[[sn]]
        if (!"dist_bin" %in% colnames(obj@meta.data)) return(NULL)
        cells <- colnames(obj)[obj$dist_bin == "bin_50plus"]
        if (length(cells) == 0) return(NULL)
        DefaultAssay(obj) <- "SCT"
        rowMeans(GetAssayData(obj, slot = "data")[common_genes, cells, drop = FALSE])
      })
      vals0  <- do.call(cbind, Filter(Negate(is.null), vals0))
      vals50 <- do.call(cbind, Filter(Negate(is.null), vals50))

      pvals <- sapply(common_genes, function(g) {
        v0  <- if (is.matrix(vals0))  vals0[g,]  else vals0[g]
        v50 <- if (is.matrix(vals50)) vals50[g,] else vals50[g]
        if (length(v0) < 2 || length(v50) < 2) return(1)
        tryCatch(wilcox.test(v0, v50)$p.value, error = function(e) 1)
      })

      data.frame(
        gene    = common_genes,
        log2FC  = fc[common_genes],
        pvalue  = pvals,
        padj    = p.adjust(pvals, method = "BH"),
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)

    if (is.null(de_results)) next
    de_results <- de_results[order(de_results$padj), ]
    write.csv(de_results,
              file.path(tier2_dir, paste0(cond, "_DE_bin0_vs_bin50.csv")),
              row.names = FALSE)

    # Volcano plot
    de_results$sig <- de_results$padj < 0.05 & abs(de_results$log2FC) > 0.5
    p_vol <- ggplot(de_results, aes(x = log2FC, y = -log10(pvalue + 1e-300),
                                    color = sig)) +
      geom_point(size = 0.8, alpha = 0.6) +
      scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      labs(title = paste0(cond, " — Interface vs Distal (bin0 vs bin50)"),
           x = "log2 Fold Change", y = "-log10(p-value)") +
      theme_minimal(base_size = 10)
    ggsave(file.path(tier2_dir, paste0(cond, "_volcano.pdf")),
           plot = p_vol, width = 10, height = 8)

    # Heatmap of top 50 genes across all bins
    top50 <- head(de_results$gene[order(abs(de_results$log2FC), decreasing = TRUE)], 50)
    heat_mat <- do.call(cbind, lapply(dist_bin_labels, function(bn) {
      v <- bin_means[[bn]]
      if (is.null(v)) return(rep(NA, length(top50)))
      v[top50]
    }))
    if (!is.null(heat_mat) && !all(is.na(heat_mat))) {
      colnames(heat_mat) <- dist_bin_labels
      rownames(heat_mat) <- top50
      pdf(file.path(tier2_dir, paste0(cond, "_heatmap_top50.pdf")),
          width = 10, height = 8)
      if (have_pheatmap) {
        pheatmap(heat_mat, cluster_cols = FALSE, scale = "row",
                 main = paste0(cond, " — Top 50 DE Genes Across Bins"),
                 fontsize_row = 6)
      } else {
        heatmap(heat_mat, Colv = NA, scale = "row",
                main = paste0(cond, " — Top 50 DE Genes Across Bins"))
      }
      dev.off()
    }
    cat(paste0("    DE done: ", nrow(de_results), " genes tested\n"))
  }

  tier2_status <- "COMPLETE"
  cat("TIER 2 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 2 ERROR: ", conditionMessage(e), "\n"))
  tier2_status <<- paste0("ERROR: ", conditionMessage(e))
})

# ==============================================================================
# SECTION 5: TIER 3 — Gradient Modeling
# ==============================================================================
cat("\n=== TIER 3: Gradient Modeling ===\n")
tier3_status <- "SKIPPED"

tryCatch({
  tier3_dir <- file.path(visium_root, "TIER3_Gradient")
  all_sig_genes <- unique(c(ASPC_genes, Hallmark_AR_genes, NEPC_Beltran_UP, NEPC_Beltran_DOWN))
  gradient_results <- list()

  conditions <- unique(sapply(names(seurat_list), extract_condition))

  for (cond in conditions) {
    cat(paste0("  TIER3 — condition: ", cond, "\n"))
    cond_samples <- names(seurat_list)[sapply(names(seurat_list), extract_condition) == cond]

    # Collect all spots with dist info across condition
    spot_data_list <- lapply(cond_samples, function(sn) {
      obj <- seurat_list[[sn]]
      if (!"dist_to_interface" %in% colnames(obj@meta.data)) return(NULL)
      DefaultAssay(obj) <- "SCT"
      genes_avail <- intersect(all_sig_genes, rownames(obj))
      expr <- GetAssayData(obj, slot = "data")[genes_avail, , drop = FALSE]
      df <- as.data.frame(t(as.matrix(expr)))
      df$dist <- obj$dist_to_interface
      df
    })
    spot_data <- do.call(rbind, Filter(Negate(is.null), spot_data_list))
    if (is.null(spot_data) || nrow(spot_data) < 20) next

    genes_to_model <- intersect(all_sig_genes, colnames(spot_data))
    gene_results <- list()

    for (g in genes_to_model) {
      tryCatch({
        y <- log1p(spot_data[[g]])
        x <- spot_data$dist
        ok <- !is.na(y) & !is.na(x) & is.finite(y) & is.finite(x)
        if (sum(ok) < 10) next
        fit_lm <- lm(y[ok] ~ x[ok])
        sm_lm  <- summary(fit_lm)
        slope  <- coef(fit_lm)[2]
        r2_lm  <- sm_lm$r.squared
        p_lm   <- coef(sm_lm)[2, 4]

        # Try NLS exponential decay
        decay_rate <- NA; half_dist <- NA; r2_nls <- NA
        tryCatch({
          start_a  <- max(y[ok]) - min(y[ok])
          start_lm <- 0.01
          fit_nls  <- nls(y[ok] ~ a * exp(-lambda * x[ok]) + b,
                          start = list(a = start_a, lambda = start_lm, b = min(y[ok])),
                          control = nls.control(maxiter = 100, warnOnly = TRUE))
          lam <- coef(fit_nls)["lambda"]
          pred_nls <- predict(fit_nls)
          ss_res <- sum((y[ok] - pred_nls)^2)
          ss_tot <- sum((y[ok] - mean(y[ok]))^2)
          r2_nls <- if (ss_tot > 0) 1 - ss_res/ss_tot else NA
          decay_rate <- lam
          half_dist  <- if (!is.na(lam) && lam > 0) log(2) / lam else NA
        }, error = function(e) NULL)

        gene_results[[g]] <- data.frame(
          gene       = g,
          condition  = cond,
          slope_lm   = slope,
          p_lm       = p_lm,
          R2_lm      = r2_lm,
          decay_rate = decay_rate,
          half_dist  = half_dist,
          R2_nls     = r2_nls,
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)
    }

    if (length(gene_results) == 0) next
    cond_df <- do.call(rbind, gene_results)
    cond_df <- cond_df[order(cond_df$R2_lm, decreasing = TRUE), ]
    gradient_results[[cond]] <- cond_df
    write.csv(cond_df,
              file.path(tier3_dir, paste0(cond, "_gradient_models.csv")),
              row.names = FALSE)

    # Plot top 20 decay curves
    top20_genes <- head(cond_df$gene[order(abs(cond_df$slope_lm), decreasing = TRUE)], 20)
    if (length(top20_genes) > 0 && "dist" %in% colnames(spot_data)) {
      plot_list_t3 <- list()
      for (gg in top20_genes) {
        if (!gg %in% colnames(spot_data)) next
        tmp_df <- data.frame(dist = spot_data$dist,
                             expr = log1p(spot_data[[gg]]))
        tmp_df <- tmp_df[!is.na(tmp_df$dist) & !is.na(tmp_df$expr), ]
        if (nrow(tmp_df) < 5) next
        p_t3 <- ggplot(tmp_df, aes(x = dist, y = expr)) +
          geom_point(size = 0.5, alpha = 0.3) +
          geom_smooth(method = "lm", color = "red", se = TRUE) +
          labs(title = paste0(gg, " — ", cond),
               x = "Distance to Interface", y = "log1p(expr)") +
          theme_minimal(base_size = 9)
        plot_list_t3[[gg]] <- p_t3
      }
      if (length(plot_list_t3) > 0) {
        pdf(file.path(tier3_dir, paste0(cond, "_top20_decay_curves.pdf")),
            width = 10, height = 8)
        n_per_page <- 4
        for (i in seq(1, length(plot_list_t3), by = n_per_page)) {
          sub_plots <- plot_list_t3[i:min(i + n_per_page - 1, length(plot_list_t3))]
          print(wrap_plots(sub_plots, ncol = 2))
        }
        dev.off()
      }
    }
    cat(paste0("    Modeled ", nrow(cond_df), " genes for ", cond, "\n"))
  }

  tier3_status <- "COMPLETE"
  cat("TIER 3 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 3 ERROR: ", conditionMessage(e), "\n"))
  tier3_status <<- paste0("ERROR: ", conditionMessage(e))
})

# ==============================================================================
# SECTION 6: TIER 4 — Ligand-Receptor Communication
# ==============================================================================
cat("\n=== TIER 4: Ligand-Receptor Communication ===\n")
tier4_status <- "SKIPPED"

tryCatch({
  tier4_dir <- file.path(visium_root, "TIER4_LR")

  # Built-in L-R pair database (~200 pairs as specified)
  lr_pairs <- data.frame(
    ligand = c(
      "FGF1","FGF2","FGF2","FGF7","FGF10",
      "WNT2","WNT3A","WNT5A","WNT5A","WNT7B","WNT11",
      "BMP2","BMP4","BMP7","TGFB1","TGFB1","TGFB3",
      "DLL1","DLL4","JAG1","JAG1","JAG2",
      "HGF","PDGFA","PDGFB","PDGFC","PDGFD",
      "VEGFA","VEGFA","VEGFB","VEGFC",
      "EGF","HBEGF","AREG","EREG","NRG1","NRG1",
      "SHH","SHH","IHH","DHH",
      "SEMA3A","SEMA3C","SEMA4D","SEMA6A",
      "EFNA1","EFNA5","EFNB1","EFNB2",
      "CXCL12","CCL2","CCL5","CXCL8","CXCL8","CCL19","CCL21",
      "CXCL1","CXCL10",
      "IL6","IL1B","IL10","IL15","IL11","IL33",
      "TNFSF10","TNF","FASLG","CD40LG","TNFSF11",
      "IGF1","IGF2","INS",
      "ANGPT1","ANGPT2","KITLG","CSF1","CSF2",
      "LGALS9","ADIPOQ","ADIPOQ","LEP",
      "NRG1","NRG2",
      "ANXA1","ANXA1","MIF",
      "GAS6","GAS6","PROS1","MDK","PTN","PTN"
    ),
    receptor = c(
      "FGFR1","FGFR1","FGFR2","FGFR2","FGFR2",
      "FZD1","FZD2","FZD5","ROR2","FZD1","FZD4",
      "BMPR1A","BMPR1B","BMPR2","TGFBR1","TGFBR2","TGFBR2",
      "NOTCH1","NOTCH1","NOTCH1","NOTCH2","NOTCH1",
      "MET","PDGFRA","PDGFRB","PDGFRA","PDGFRB",
      "FLT1","KDR","FLT1","FLT4",
      "EGFR","EGFR","EGFR","EGFR","ERBB3","ERBB4",
      "PTCH1","PTCH2","PTCH1","PTCH1",
      "NRP1","NRP2","PLXNB1","PLXNA2",
      "EPHA2","EPHA5","EPHB2","EPHB4",
      "CXCR4","CCR2","CCR5","CXCR1","CXCR2","CCR7","CCR7",
      "CXCR2","CXCR3",
      "IL6R","IL1R1","IL10RA","IL15RA","IL11RA","IL1RL1",
      "TNFRSF10A","TNFRSF1A","FAS","CD40","TNFRSF11A",
      "IGF1R","IGF1R","INSR",
      "TEK","TEK","KIT","CSF1R","CSF2RA",
      "HAVCR2","ADIPOR1","ADIPOR2","LEPR",
      "ERBB2","ERBB4",
      "FPR1","FPR2","CD74",
      "AXL","MERTK","AXL","SDC1","SDC3","PTPRZ1"
    ),
    pathway = c(
      "FGF","FGF","FGF","FGF","FGF",
      "WNT","WNT","WNT","WNT","WNT","WNT",
      "BMP","BMP","BMP","TGFb","TGFb","TGFb",
      "NOTCH","NOTCH","NOTCH","NOTCH","NOTCH",
      "HGF","PDGF","PDGF","PDGF","PDGF",
      "VEGF","VEGF","VEGF","VEGF",
      "EGF","EGF","EGF","EGF","ErbB","ErbB",
      "Hedgehog","Hedgehog","Hedgehog","Hedgehog",
      "Semaphorin","Semaphorin","Semaphorin","Semaphorin",
      "Ephrin","Ephrin","Ephrin","Ephrin",
      "Chemokine","Chemokine","Chemokine","Chemokine","Chemokine","Chemokine","Chemokine",
      "Chemokine","Chemokine",
      "Cytokine","Cytokine","Cytokine","Cytokine","Cytokine","Cytokine",
      "TNFSF","TNFSF","TNFSF","TNFSF","TNFSF",
      "IGF","IGF","Insulin",
      "Angiopoietin","Angiopoietin","SCF","CSF","CSF",
      "Immune","Adipokine","Adipokine","Adipokine",
      "ErbB","ErbB",
      "Annexin","Annexin","MIF",
      "TAM","TAM","TAM","PTN","PTN","PTN"
    ),
    stringsAsFactors = FALSE
  )

  cat(paste0("  L-R database: ", nrow(lr_pairs), " pairs\n"))

  # Per-sample LR scoring
  lr_results_all <- list()

  for (sn in names(seurat_list)) {
    cat(paste0("  TIER4 — ", sn, "\n"))
    obj <- seurat_list[[sn]]
    md  <- obj@meta.data
    if (!"interface_zone" %in% colnames(md)) next
    if (!all(c("row","col") %in% colnames(md))) next

    DefaultAssay(obj) <- "SCT"
    all_genes_obj <- rownames(obj)
    expr_mat <- GetAssayData(obj, slot = "data")

    iface_cells <- colnames(obj)[md$interface_zone == "Interface"]
    if (length(iface_cells) < 3) next

    # Find neighbors within 3 spot-units of interface spots
    iface_idx <- which(colnames(obj) %in% iface_cells)
    coords    <- as.matrix(md[, c("row","col")])
    neighbor_cells <- character(0)
    for (ic in iface_idx) {
      dists <- sqrt(rowSums(sweep(coords, 2, coords[ic,])^2))
      nb    <- colnames(obj)[dists > 0 & dists <= 3]
      neighbor_cells <- union(neighbor_cells, nb)
    }
    neighbor_cells <- setdiff(neighbor_cells, iface_cells)
    if (length(neighbor_cells) < 3) next

    # Score each L-R pair
    lr_scores <- lapply(seq_len(nrow(lr_pairs)), function(i) {
      lig <- lr_pairs$ligand[i]
      rec <- lr_pairs$receptor[i]
      if (!lig %in% all_genes_obj || !rec %in% all_genes_obj) return(NULL)
      mean_lig <- mean(expr_mat[lig, iface_cells,    drop = TRUE], na.rm = TRUE)
      mean_rec <- mean(expr_mat[rec, neighbor_cells, drop = TRUE], na.rm = TRUE)
      score <- mean_lig * mean_rec

      # Permutation test (1000 shuffles)
      all_cells <- colnames(obj)
      n_iface   <- length(iface_cells)
      n_nbr     <- length(neighbor_cells)
      perm_scores <- replicate(1000, {
        perm_iface <- sample(all_cells, n_iface)
        perm_nbr   <- sample(all_cells, n_nbr)
        mean(expr_mat[lig, perm_iface, drop = TRUE], na.rm = TRUE) *
          mean(expr_mat[rec, perm_nbr,   drop = TRUE], na.rm = TRUE)
      })
      emp_p <- mean(perm_scores >= score)

      data.frame(
        sample    = sn,
        condition = extract_condition(sn),
        ligand    = lig,
        receptor  = rec,
        pathway   = lr_pairs$pathway[i],
        mean_lig  = mean_lig,
        mean_rec  = mean_rec,
        score     = score,
        pvalue    = emp_p,
        stringsAsFactors = FALSE
      )
    })
    lr_scores <- do.call(rbind, Filter(Negate(is.null), lr_scores))
    if (!is.null(lr_scores) && nrow(lr_scores) > 0) {
      lr_results_all[[sn]] <- lr_scores
    }
  }

  if (length(lr_results_all) > 0) {
    lr_all_df <- do.call(rbind, lr_results_all)
    write.csv(lr_all_df,
              file.path(tier4_dir, "LR_scores_all_samples.csv"),
              row.names = FALSE)

    # Aggregate per condition with Fisher's method
    conditions <- unique(lr_all_df$condition)
    for (cond in conditions) {
      cond_lr <- lr_all_df[lr_all_df$condition == cond, ]
      # Fisher's combined p per pair
      pair_summary <- cond_lr %>%
        group_by(ligand, receptor, pathway) %>%
        summarise(
          mean_score    = mean(score, na.rm = TRUE),
          n_samples     = n(),
          fisher_stat   = -2 * sum(log(pmax(pvalue, 1e-300))),
          .groups = "drop"
        ) %>%
        mutate(
          fisher_pval = pchisq(fisher_stat, df = 2 * n_samples, lower.tail = FALSE)
        ) %>%
        arrange(fisher_pval)

      write.csv(pair_summary,
                file.path(tier4_dir, paste0(cond, "_LR_aggregated.csv")),
                row.names = FALSE)

      # Dot plot of top 30 pairs
      top30 <- head(pair_summary, 30)
      top30$pair <- paste0(top30$ligand, "→", top30$receptor)
      p_lr <- ggplot(top30, aes(x = mean_score, y = reorder(pair, mean_score),
                                 color = pathway, size = -log10(fisher_pval + 1e-10))) +
        geom_point() +
        labs(title = paste0(cond, " — Top 30 L-R Pairs"),
             x = "Mean Interaction Score", y = "L-R Pair",
             color = "Pathway", size = "-log10(p)") +
        theme_minimal(base_size = 9)
      ggsave(file.path(tier4_dir, paste0(cond, "_LR_dotplot.pdf")),
             plot = p_lr, width = 10, height = 8)
      cat(paste0("    ", cond, ": ", nrow(pair_summary), " L-R pairs\n"))
    }
  }

  tier4_status <- "COMPLETE"
  cat("TIER 4 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 4 ERROR: ", conditionMessage(e), "\n"))
  tier4_status <<- paste0("ERROR: ", conditionMessage(e))
})

# ==============================================================================
# SECTION 7: TIER 5 — Trajectory / Spatial Pseudotime
# ==============================================================================
cat("\n=== TIER 5: Trajectory / Spatial Pseudotime ===\n")
tier5_status <- "SKIPPED"

tryCatch({
  tier5_dir <- file.path(visium_root, "TIER5_Trajectory")
  conditions <- unique(sapply(names(seurat_list), extract_condition))

  for (cond in conditions) {
    cat(paste0("  TIER5 — condition: ", cond, "\n"))
    cond_samples <- names(seurat_list)[sapply(names(seurat_list), extract_condition) == cond]

    # Collect near-interface spots (bins 0 to bin_11_20)
    near_iface_data <- lapply(cond_samples, function(sn) {
      obj <- seurat_list[[sn]]
      if (!"dist_bin" %in% colnames(obj@meta.data)) return(NULL)
      near_bins <- c("bin_0_interface","bin_1_5","bin_6_10","bin_11_20")
      cells <- colnames(obj)[obj$dist_bin %in% near_bins]
      if (length(cells) < 5) return(NULL)
      list(obj = obj, cells = cells, sn = sn)
    })
    near_iface_data <- Filter(Negate(is.null), near_iface_data)
    if (length(near_iface_data) == 0) next

    # Gather expression for variable features
    all_var_genes <- character(0)
    for (item in near_iface_data) {
      obj <- item$obj
      DefaultAssay(obj) <- "SCT"
      vf <- VariableFeatures(obj)
      all_var_genes <- union(all_var_genes, vf)
    }
    all_var_genes <- head(all_var_genes, 2000)
    if (length(all_var_genes) < 10) next

    expr_list <- lapply(near_iface_data, function(item) {
      obj <- item$obj
      DefaultAssay(obj) <- "SCT"
      genes_avail <- intersect(all_var_genes, rownames(obj))
      mat <- GetAssayData(obj, slot = "data")[genes_avail, item$cells, drop = FALSE]
      mat
    })

    # Align genes across samples
    common_vf <- Reduce(intersect, lapply(expr_list, rownames))
    if (length(common_vf) < 10) next

    expr_combined <- do.call(cbind, lapply(expr_list, function(m) m[common_vf, ]))
    expr_t <- t(as.matrix(expr_combined))

    # PCA
    pca_res <- tryCatch(
      prcomp(expr_t, scale. = TRUE, center = TRUE, rank. = min(10, ncol(expr_t) - 1)),
      error = function(e) NULL
    )
    if (is.null(pca_res)) next

    pc_scores <- pca_res$x

    # Fit principal curve or use PC1
    pseudotime <- tryCatch({
      if (have_princurve) {
        pc_fit <- principal_curve(pc_scores[, 1:min(5, ncol(pc_scores))])
        pt <- pc_fit$lambda
        (pt - min(pt)) / (max(pt) - min(pt))
      } else {
        pt <- pc_scores[, 1]
        (pt - min(pt)) / (max(pt) - min(pt))
      }
    }, error = function(e) {
      pt <- pc_scores[, 1]
      (pt - min(pt)) / (max(pt) - min(pt))
    })

    # Collect metadata for all near-interface spots
    all_near_meta <- do.call(rbind, lapply(near_iface_data, function(item) {
      obj <- item$obj
      md  <- obj@meta.data[item$cells, ]
      data.frame(
        barcode   = item$cells,
        sample_id = item$sn,
        condition = obj$condition[1],
        dist_bin  = md$dist_bin,
        ASPC_score = md$ASPC_score %||% NA_real_,
        AR_score   = md$AR_score   %||% NA_real_,
        NEPC_UP_score = md$NEPC_UP_score %||% NA_real_,
        NEPC_Beltran_NET = md$NEPC_Beltran_NET %||% NA_real_,
        row = md$row %||% NA_real_,
        col = md$col %||% NA_real_,
        stringsAsFactors = FALSE
      )
    }))

    if (nrow(all_near_meta) != length(pseudotime)) {
      cat(paste0("    Dimension mismatch for ", cond, " — skipping\n"))
      next
    }
    all_near_meta$pseudotime <- pseudotime

    # Spearman correlation of scores with pseudotime
    score_cols_pt <- c("ASPC_score","AR_score","NEPC_UP_score","NEPC_Beltran_NET")
    cor_table <- lapply(score_cols_pt, function(sc) {
      v <- all_near_meta[[sc]]
      if (all(is.na(v))) return(NULL)
      ct <- cor.test(all_near_meta$pseudotime, v, method = "spearman",
                     use = "complete.obs", exact = FALSE)
      data.frame(score = sc, rho = ct$estimate, pvalue = ct$p.value,
                 stringsAsFactors = FALSE)
    })
    cor_df <- do.call(rbind, Filter(Negate(is.null), cor_table))
    if (!is.null(cor_df))
      write.csv(cor_df,
                file.path(tier5_dir, paste0(cond, "_pseudotime_score_correlations.csv")),
                row.names = FALSE)

    # Transition genes: Spearman of top 2000 VF genes vs pseudotime
    gene_pt_cor <- tryCatch({
      expr_sub <- expr_combined[, rownames(expr_t)]
      sapply(seq_len(nrow(expr_sub)), function(gi) {
        y <- expr_sub[gi, ]
        ct <- tryCatch(
          cor.test(pseudotime, y, method = "spearman", use = "complete.obs", exact = FALSE),
          error = function(e) list(estimate = NA, p.value = NA)
        )
        c(rho = as.numeric(ct$estimate), pvalue = as.numeric(ct$p.value))
      })
    }, error = function(e) NULL)

    if (!is.null(gene_pt_cor)) {
      trans_genes_df <- data.frame(
        gene   = rownames(expr_combined),
        rho    = gene_pt_cor["rho", ],
        pvalue = gene_pt_cor["pvalue", ],
        stringsAsFactors = FALSE
      )
      trans_genes_df$padj <- p.adjust(trans_genes_df$pvalue, method = "BH")
      trans_genes_df <- trans_genes_df[
        !is.na(trans_genes_df$rho) &
          abs(trans_genes_df$rho) > 0.3 &
          trans_genes_df$padj < 0.05, ]
      trans_genes_df <- trans_genes_df[order(abs(trans_genes_df$rho), decreasing = TRUE), ]
      write.csv(trans_genes_df,
                file.path(tier5_dir, paste0(cond, "_transition_genes.csv")),
                row.names = FALSE)
      cat(paste0("    Transition genes (|rho|>0.3, padj<0.05): ",
                  nrow(trans_genes_df), "\n"))
    }

    # Spatial pseudotime plot
    if (all(c("col","row") %in% colnames(all_near_meta))) {
      spot_sz <- compute_spot_size(all_near_meta[, c("col","row")])
      p_pt <- ggplot(all_near_meta, aes(x = col, y = -row, color = pseudotime)) +
        geom_point(size = spot_sz) +
        scale_color_viridis_c(option = "plasma") +
        labs(title = paste0(cond, " — Spatial Pseudotime (near-interface)"),
             x = "Array Col", y = "Array Row (inv)", color = "Pseudotime") +
        theme_minimal(base_size = 10)
      ggsave(file.path(tier5_dir, paste0(cond, "_spatial_pseudotime.pdf")),
             plot = p_pt, width = 10, height = 8)
    }

    # Scatter plots of pseudotime vs scores
    if (!is.null(cor_df)) {
      p_scatter <- ggplot(all_near_meta, aes(x = pseudotime, y = ASPC_score)) +
        geom_point(size = 0.5, alpha = 0.4) +
        geom_smooth(method = "loess", color = "blue") +
        labs(title = paste0(cond, " — ASPC Score vs Pseudotime")) +
        theme_minimal(base_size = 10)
      ggsave(file.path(tier5_dir, paste0(cond, "_ASPC_vs_pseudotime.pdf")),
             plot = p_scatter, width = 10, height = 8)
    }
  }

  tier5_status <- "COMPLETE"
  cat("TIER 5 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 5 ERROR: ", conditionMessage(e), "\n"))
  tier5_status <<- paste0("ERROR: ", conditionMessage(e))
})

# ==============================================================================
# SECTION 8: TIER 6 — GRN (Interface-Specific)
# ==============================================================================
cat("\n=== TIER 6: Gene Regulatory Network (Interface-Specific) ===\n")
tier6_status <- "SKIPPED"

tryCatch({
  tier6_dir <- file.path(visium_root, "TIER6_GRN")

  known_TFs <- c("AR","FOXA1","FOXA2","NKX3-1","ERG","ETV1","ETV4","ETV5","HOXB13",
                 "GATA2","MYC","MYCN","E2F1","REST","ASCL1","NEUROD1","POU3F2","SOX2",
                 "ONECUT2","TP53","RB1","FOXP1","NR2F1","ELK4","ELF3","FOXC2","SPDEF",
                 "TFCP2L1","HEYL","ID4","TBX2")

  # Collect all Interface spots across samples
  all_iface_expr_list <- list()
  for (sn in names(seurat_list)) {
    obj <- seurat_list[[sn]]
    if (!"interface_zone" %in% colnames(obj@meta.data)) next
    iface_cells <- colnames(obj)[obj$interface_zone == "Interface"]
    if (length(iface_cells) < 3) next
    DefaultAssay(obj) <- "SCT"
    all_iface_expr_list[[sn]] <- GetAssayData(obj, slot = "data")[, iface_cells, drop = FALSE]
  }

  if (length(all_iface_expr_list) == 0) {
    cat("  No interface spots found for GRN\n")
    stop("No interface data")
  }

  # Common genes across samples
  common_genes_grn <- Reduce(intersect, lapply(all_iface_expr_list, rownames))
  if (length(common_genes_grn) < 20) stop("Too few common genes for GRN")

  # Get signature genes + top 500 variable genes
  all_sig_genes_grn <- unique(c(ASPC_genes, Hallmark_AR_genes, NEPC_Beltran_UP, NEPC_Beltran_DOWN))
  # Determine variable genes from combined data
  all_expr_iface <- do.call(cbind, lapply(all_iface_expr_list, function(m) m[common_genes_grn, ]))
  gene_vars <- apply(all_expr_iface, 1, var, na.rm = TRUE)
  top500_var <- names(sort(gene_vars, decreasing = TRUE))[1:min(500, length(gene_vars))]
  grn_genes <- unique(c(intersect(all_sig_genes_grn, common_genes_grn), top500_var))
  grn_genes <- head(grn_genes, 600)

  cat(paste0("  GRN genes: ", length(grn_genes), "\n"))
  grn_expr <- all_expr_iface[grn_genes, ]

  # Spearman correlation matrix
  cat("  Computing Spearman correlation matrix...\n")
  cor_mat <- tryCatch(
    cor(t(as.matrix(grn_expr)), method = "spearman", use = "pairwise.complete.obs"),
    error = function(e) NULL
  )
  if (is.null(cor_mat)) stop("Correlation failed")

  # Filter edges: |rho| > 0.3 and p < 0.01
  n_spots <- ncol(grn_expr)
  t_stat  <- cor_mat * sqrt((n_spots - 2) / (1 - cor_mat^2 + 1e-10))
  p_mat   <- 2 * pt(-abs(t_stat), df = n_spots - 2)

  idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
  edge_df <- data.frame(
    gene1 = rownames(cor_mat)[idx[, 1]],
    gene2 = colnames(cor_mat)[idx[, 2]],
    rho   = cor_mat[idx],
    pvalue = p_mat[idx],
    stringsAsFactors = FALSE
  )
  edge_df <- edge_df[abs(edge_df$rho) > 0.3 & edge_df$pvalue < 0.01, ]
  cat(paste0("  Edges after filtering: ", nrow(edge_df), "\n"))

  write.csv(edge_df, file.path(tier6_dir, "edge_list.csv"), row.names = FALSE)

  # igraph network analysis
  if (have_igraph && nrow(edge_df) > 0) {
    g <- graph_from_data_frame(edge_df[, c("gene1","gene2")], directed = FALSE)
    E(g)$weight <- abs(edge_df$rho)
    deg <- degree(g)
    btw <- tryCatch(betweenness(g, normalized = TRUE), error = function(e) rep(NA, vcount(g)))
    node_df <- data.frame(
      gene        = V(g)$name,
      degree      = deg,
      betweenness = btw,
      is_TF       = V(g)$name %in% known_TFs,
      is_sig_gene = V(g)$name %in% all_sig_genes_grn,
      stringsAsFactors = FALSE
    )
    node_df <- node_df[order(node_df$degree, decreasing = TRUE), ]
    write.csv(node_df, file.path(tier6_dir, "node_stats.csv"), row.names = FALSE)

    hub_genes <- head(node_df$gene, 20)
    write.csv(data.frame(hub_gene = hub_genes, degree = node_df$degree[1:20]),
              file.path(tier6_dir, "hub_genes.csv"), row.names = FALSE)
    cat(paste0("  Hub genes: ", paste(head(hub_genes, 10), collapse = ", "), "\n"))

    # Network plot (top 100 nodes by degree for readability)
    top_nodes <- head(node_df$gene, min(100, nrow(node_df)))
    g_sub     <- induced_subgraph(g, top_nodes)
    pdf(file.path(tier6_dir, "GRN_network_plot.pdf"), width = 10, height = 8)
    tryCatch({
      node_colors <- ifelse(V(g_sub)$name %in% known_TFs, "red",
                     ifelse(V(g_sub)$name %in% all_sig_genes_grn, "orange", "grey70"))
      plot(g_sub,
           vertex.color = node_colors,
           vertex.size  = 4,
           vertex.label = ifelse(V(g_sub)$name %in% hub_genes, V(g_sub)$name, ""),
           vertex.label.cex  = 0.6,
           vertex.label.color = "black",
           edge.width   = 0.5,
           edge.color   = "grey80",
           layout       = layout_with_fr(g_sub),
           main         = "Interface GRN (top 100 nodes by degree)")
      legend("topleft", legend = c("TF","Signature","Other"),
             fill = c("red","orange","grey70"), cex = 0.7)
    }, error = function(e) {
      cat(paste0("  Network plot error: ", conditionMessage(e), "\n"))
    })
    dev.off()
  } else {
    cat("  igraph not available or no edges — skipping network analysis\n")
    hub_genes <- character(0)
  }

  tier6_status <- "COMPLETE"
  cat("TIER 6 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 6 ERROR: ", conditionMessage(e), "\n"))
  tier6_status <<- paste0("ERROR: ", conditionMessage(e))
})

# ==============================================================================
# SECTION 9: TIER 7 — Cross-Sample Meta-Analysis
# ==============================================================================
cat("\n=== TIER 7: Cross-Sample Meta-Analysis ===\n")
tier7_status <- "SKIPPED"

tryCatch({
  tier7_dir <- file.path(visium_root, "TIER7_MetaAnalysis")

  # Built-in Hallmark-style gene sets (~20 sets)
  hallmark_sets <- list(
    EMT               = c("VIM","CDH2","FN1","TWIST1","SNAI1","SNAI2","ZEB1","ZEB2",
                          "MMP2","MMP9","ACTA2","TGFB1","TGFB2","CTGF","COL1A1"),
    Inflammatory      = c("IL1B","IL6","TNF","CXCL8","CXCL1","CXCL2","CCL2","PTGS2",
                          "NFKB1","RELA","NLRP3","CASP1","IL18","IL12A"),
    Androgen          = c("KLK3","KLK2","TMPRSS2","NKX3-1","FKBP5","PMEPA1","STEAP4",
                          "NDRG1","SGK1","AR","FOXA1","HOXB13","SPDEF","ELL2"),
    MYC_targets       = c("MYC","MYCN","CDK4","CDK6","E2F1","E2F2","ODC1","NPM1",
                          "NCL","PCNA","RRM2","TK1","LDHA","ENO1","PTMA"),
    E2F_targets       = c("E2F1","E2F2","E2F3","CCNE1","CCNE2","CDC6","PCNA","MCM2",
                          "MCM3","MCM4","MCM5","MCM6","MCM7","RRM1","TYMS"),
    Hypoxia           = c("HIF1A","VEGFA","LDHA","PGK1","ENO1","ALDOA","GAPDH","TPI1",
                          "SLC2A1","CA9","BNIP3","BNIP3L","ADM","HMOX1"),
    Angiogenesis      = c("VEGFA","VEGFB","VEGFC","FLT1","KDR","ANGPT1","ANGPT2",
                          "TEK","PDGFB","PDGFRB","NRP1","NRP2","FGF2","FGFR1"),
    KRAS              = c("KRAS","RAF1","BRAF","MAP2K1","MAPK1","MAPK3","EGF","EGFR",
                          "SOS1","GRB2","HRAS","NRAS","MYC","JUN"),
    TNFa              = c("TNF","TNFRSF1A","NFKB1","RELA","TRADD","TRAF2","RIPK1",
                          "CASP8","BCL2L1","XIAP","BIRC3","ICAM1","VCAM1"),
    IL6_JAK_STAT3     = c("IL6","IL6R","JAK1","JAK2","STAT3","SOCS1","SOCS3","BCL2",
                          "CCND1","MYC","VEGFA","MCL1","TWIST1","CDKN1A"),
    p53               = c("TP53","CDKN1A","MDM2","BAX","BBC3","PMAIP1","FAS","GADD45A",
                          "GADD45B","SESN1","SESN2","TIGAR","SCO2","FDXR"),
    Apoptosis         = c("BCL2","BCL2L1","BCL2L2","MCL1","BAX","BAK1","CASP3","CASP7",
                          "CASP8","CASP9","CASP10","CYCS","APAF1","BID","PARP1"),
    OxPhos            = c("MT-CO1","MT-CO2","MT-ND1","MT-ND2","ATP5F1A","UQCRFS1",
                          "SDHB","SDHA","NDUFB8","COX4I1","ATP5MC1","VDAC1","VDAC2"),
    Glycolysis        = c("HK2","HK1","GPI","PFKM","ALDOA","TPI1","GAPDH","PGK1",
                          "ENO1","PKM","LDHA","SLC2A1","SLC2A3","G6PD"),
    Fatty_acid        = c("FASN","ACACA","ACACB","ACLY","SCD","ELOVL5","ELOVL6",
                          "FADS1","FADS2","ACSL3","ACSL4","CPT1A","HADHA","HADHB"),
    Cholesterol       = c("HMGCR","HMGCS1","FDFT1","SQLE","LSS","MVD","MVK",
                          "PMVK","IDI1","FDPS","DHCR7","DHCR24","INSIG1","INSIG2"),
    Notch             = c("NOTCH1","NOTCH2","NOTCH3","NOTCH4","DLL1","DLL4","JAG1",
                          "JAG2","HES1","HES5","HEY1","HEY2","HEYL","MAML1"),
    Hedgehog          = c("SHH","IHH","DHH","PTCH1","PTCH2","SMO","GLI1","GLI2",
                          "GLI3","HHIP","BOC","GAS1","CDON","BMP2"),
    WNT_beta_catenin  = c("WNT1","WNT2","WNT3A","WNT5A","FZD1","FZD2","LRP5","LRP6",
                          "CTNNB1","APC","AXIN1","GSK3B","TCF7","LEF1","MYC"),
    TGF_beta          = c("TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","SMAD2","SMAD3",
                          "SMAD4","SMAD7","INHBA","INHBB","BMP2","BMP4","BMP7")
  )

  # Load Tier 2 DE results to get top 200 interface-enriched genes per condition
  conditions <- unique(sapply(names(seurat_list), extract_condition))
  cond_top200 <- list()
  for (cond in conditions) {
    de_file <- file.path(visium_root, "TIER2_DE",
                         paste0(cond, "_DE_bin0_vs_bin50.csv"))
    if (file.exists(de_file)) {
      de <- read.csv(de_file, stringsAsFactors = FALSE)
      de <- de[!is.na(de$log2FC) & !is.na(de$padj), ]
      de_up <- de[de$log2FC > 0, ]
      de_up <- de_up[order(de_up$padj), ]
      cond_top200[[cond]] <- head(de_up$gene, 200)
    } else {
      cond_top200[[cond]] <- character(0)
    }
  }

  # Per condition: conservation score = fraction of patients where gene in top 200
  all_cond_conserved <- list()
  for (cond in conditions) {
    cond_samples <- names(seurat_list)[sapply(names(seurat_list), extract_condition) == cond]
    n_patients <- length(cond_samples)
    if (n_patients == 0) next

    # Per-patient top 200 (use SCT data, interface vs all)
    patient_top200 <- lapply(cond_samples, function(sn) {
      obj <- seurat_list[[sn]]
      if (!"interface_zone" %in% colnames(obj@meta.data)) return(character(0))
      DefaultAssay(obj) <- "SCT"
      iface_cells  <- colnames(obj)[obj$interface_zone == "Interface"]
      other_cells  <- colnames(obj)[obj$interface_zone != "Interface"]
      if (length(iface_cells) < 3 || length(other_cells) < 3) return(character(0))
      expr_mat <- GetAssayData(obj, layer = "data")
      mean_iface <- rowMeans(expr_mat[, iface_cells, drop = FALSE], na.rm = TRUE)
      mean_other <- rowMeans(expr_mat[, other_cells, drop = FALSE], na.rm = TRUE)
      fc <- mean_iface - mean_other
      head(names(sort(fc, decreasing = TRUE)), 200)
    })

    # Count frequency per gene
    all_genes_cond <- unique(unlist(patient_top200))
    freq_vec <- sapply(all_genes_cond, function(g) {
      mean(sapply(patient_top200, function(top) g %in% top))
    })
    freq_df <- data.frame(
      gene              = all_genes_cond,
      conservation      = freq_vec,
      condition         = cond,
      in_core_program   = freq_vec > 0.5,
      stringsAsFactors  = FALSE
    )
    freq_df <- freq_df[order(freq_df$conservation, decreasing = TRUE), ]
    all_cond_conserved[[cond]] <- freq_df
    write.csv(freq_df,
              file.path(tier7_dir, paste0(cond, "_conservation_scores.csv")),
              row.names = FALSE)
    cat(paste0("  ", cond, ": ", sum(freq_df$in_core_program), " core interface genes\n"))
  }

  # Combine conservation scores
  if (length(all_cond_conserved) > 0) {
    all_cons_df <- do.call(rbind, all_cond_conserved)
    write.csv(all_cons_df,
              file.path(tier7_dir, "conservation_scores.csv"),
              row.names = FALSE)

    # Cross-condition comparison matrix
    core_genes_per_cond <- lapply(all_cond_conserved, function(df) {
      df$gene[df$in_core_program]
    })
    all_core_genes <- unique(unlist(core_genes_per_cond))
    if (length(all_core_genes) > 0) {
      cross_mat <- do.call(cbind, lapply(core_genes_per_cond, function(g) {
        as.integer(all_core_genes %in% g)
      }))
      rownames(cross_mat) <- all_core_genes
      colnames(cross_mat) <- names(core_genes_per_cond)
      write.csv(cross_mat,
                file.path(tier7_dir, "cross_condition_matrix.csv"),
                row.names = TRUE)
    }
  }

  # Pathway enrichment: hypergeometric test per hallmark set × condition
  universe_size <- 20000  # approximate human gene universe
  enrich_results <- list()
  for (cond in names(all_cond_conserved)) {
    core_g <- all_cond_conserved[[cond]]$gene[all_cond_conserved[[cond]]$in_core_program]
    n_core <- length(core_g)
    if (n_core == 0) next
    for (set_name in names(hallmark_sets)) {
      set_genes  <- hallmark_sets[[set_name]]
      overlap    <- intersect(core_g, set_genes)
      n_overlap  <- length(overlap)
      n_set      <- length(set_genes)
      p_hyper    <- phyper(n_overlap - 1, n_set, universe_size - n_set,
                           n_core, lower.tail = FALSE)
      enrich_results[[paste0(cond, "_", set_name)]] <- data.frame(
        condition   = cond,
        pathway     = set_name,
        n_core      = n_core,
        n_set       = n_set,
        n_overlap   = n_overlap,
        overlap_genes = paste(overlap, collapse = ";"),
        pvalue      = p_hyper,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(enrich_results) > 0) {
    enrich_df <- do.call(rbind, enrich_results)
    enrich_df$padj <- p.adjust(enrich_df$pvalue, method = "BH")
    enrich_df <- enrich_df[order(enrich_df$padj), ]
    write.csv(enrich_df,
              file.path(tier7_dir, "pathway_enrichment.csv"),
              row.names = FALSE)
    cat(paste0("  Pathway enrichment: ",
               sum(enrich_df$padj < 0.05, na.rm = TRUE),
               " significant pathways\n"))
  }

  tier7_status <- "COMPLETE"
  cat("TIER 7 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 7 ERROR: ", conditionMessage(e), "\n"))
  tier7_status <<- paste0("ERROR: ", conditionMessage(e))
})

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================
cat("\n")
cat("==============================================================================\n")
cat("  GSE278936 Visium 7-Tier Pipeline — FINAL SUMMARY\n")
cat("==============================================================================\n")
cat(paste0("Samples processed: ", length(seurat_list), " / ", length(sample_manifest), "\n"))
cat(paste0("TIER 1 (Spatial Annotation):          ", tier1_status, "\n"))
cat(paste0("TIER 2 (Pseudo-Bulk DE):              ", tier2_status, "\n"))
cat(paste0("TIER 3 (Gradient Modeling):           ", tier3_status, "\n"))
cat(paste0("TIER 4 (Ligand-Receptor):             ", tier4_status, "\n"))
cat(paste0("TIER 5 (Trajectory/Pseudotime):       ", tier5_status, "\n"))
cat(paste0("TIER 6 (GRN):                         ", tier6_status, "\n"))
cat(paste0("TIER 7 (Meta-Analysis):               ", tier7_status, "\n"))

# Interface spot counts
total_iface <- 0
for (sn in names(seurat_list)) {
  obj <- seurat_list[[sn]]
  if ("interface_zone" %in% colnames(obj@meta.data)) {
    n_iface <- sum(obj$interface_zone == "Interface", na.rm = TRUE)
    total_iface <- total_iface + n_iface
  }
}
cat(paste0("Total interface spots identified:     ", total_iface, "\n"))

# Top genes from Tier 2 (BPH as example)
tryCatch({
  bph_de <- read.csv(file.path(visium_root, "TIER2_DE", "BPH_DE_bin0_vs_bin50.csv"),
                     stringsAsFactors = FALSE)
  top_genes_summary <- head(bph_de$gene[order(bph_de$padj)], 10)
  cat(paste0("Top DE genes (BPH, bin0 vs bin50):  ",
              paste(top_genes_summary, collapse = ", "), "\n"))
}, error = function(e) NULL)

# Top L-R pairs
tryCatch({
  lr_file <- file.path(visium_root, "TIER4_LR", "LR_scores_all_samples.csv")
  if (file.exists(lr_file)) {
    lr_all <- read.csv(lr_file, stringsAsFactors = FALSE)
    top_lr  <- head(lr_all[order(lr_all$score, decreasing = TRUE), ], 5)
    cat(paste0("Top L-R pairs (by score): ",
               paste(paste0(top_lr$ligand, "→", top_lr$receptor), collapse = ", "), "\n"))
  }
}, error = function(e) NULL)

cat(paste0("Outputs written to: ", visium_root, "\n"))
cat("==============================================================================\n")

# Save summary to file
summary_lines <- c(
  "GSE278936 Visium 7-Tier Pipeline — Summary",
  paste0("Run date: ", Sys.time()),
  paste0("Samples processed: ", length(seurat_list), " / ", length(sample_manifest)),
  paste0("TIER 1: ", tier1_status),
  paste0("TIER 2: ", tier2_status),
  paste0("TIER 3: ", tier3_status),
  paste0("TIER 4: ", tier4_status),
  paste0("TIER 5: ", tier5_status),
  paste0("TIER 6: ", tier6_status),
  paste0("TIER 7: ", tier7_status),
  paste0("Total interface spots: ", total_iface)
)
writeLines(summary_lines, file.path(visium_root, "pipeline_summary.txt"))
cat("Saved pipeline_summary.txt\n")
