# ==============================================================================
# GSE278936 Visium Spatial Transcriptomics — FULL Multi-Sample Pipeline
# ASPC, Hallmark AR & NEPC Beltran Signatures
# Processes ALL 52 samples automatically
# Dynamic spot sizing for optimal gap-free spatial plots
# MEMORY-OPTIMIZED: process-save-clear per sample, reload as needed
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

# Increase memory limit if possible (macOS R)
tryCatch({
  if (.Platform$OS.type == "unix") {
    # Try to set mem.maxVSize to available system memory
    mem_info <- tryCatch(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)) / 1e9,
                         error = function(e) NA)
    if (!is.na(mem_info) && mem_info > 16) {
      cat(paste0("System memory: ", round(mem_info, 1), " GB\n"))
    }
  }
}, error = function(e) NULL)

# ------------------------------------------------------------------------------
# 1. Define Paths — Auto-detect all samples & organize GEO-prefixed files
# ------------------------------------------------------------------------------
visium_root <- "/Users/cj719/Library/Mobile Documents/com~apple~CloudDocs/Adipocyte NEPC/cj719/GSE278936 Visium"

# Directory for individual sample RDS files (memory management)
rds_dir <- file.path(visium_root, "sample_rds")
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

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

# Null-coalescing operator
`%||%` <- function(x, y) if (!is.null(x)) x else y

# --- Helper: select best available assay (SCT preferred, RNA fallback) ---
set_best_assay <- function(obj) {
  if ("SCT" %in% Assays(obj)) {
    DefaultAssay(obj) <- "SCT"
  } else {
    DefaultAssay(obj) <- "RNA"
  }
  return(obj)
}

# --- Helper: safely get assay data from object (strips images first) ---
safe_GetAssayData <- function(obj, layer = "data") {
  img_backup <- obj@images
  obj@images <- list()
  mat <- tryCatch(
    GetAssayData(obj, layer = layer),
    error = function(e) NULL
  )
  # No need to restore — obj is not modified in caller scope
  mat
}

# --- Helper: load a sample from RDS, return NULL if missing ---
load_sample <- function(sn) {
  rds_path <- file.path(rds_dir, paste0(sn, ".rds"))
  if (file.exists(rds_path)) {
    return(readRDS(rds_path))
  }
  return(NULL)
}

# --- Additional packages with graceful fallback ---
extra_pkgs <- list(
  pheatmap  = list(cran = "pheatmap"),
  igraph    = list(cran = "igraph"),
  MASS      = list(cran = "MASS"),
  parallel  = list(cran = NULL)
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

# ==============================================================================
# SECTION 2b: Load Samples — Process-Save-Clear Strategy
# ==============================================================================
# Strategy: process one sample at a time, save to RDS, keep only the
# sample name list in memory. Each tier reloads samples on demand.
# This prevents memory accumulation across 52 samples.
# ==============================================================================
cat("\n=== Loading samples into Seurat (process-save-clear) ===\n")
successfully_loaded <- character(0)

for (sn in names(sample_manifest)) {
  cat(paste0("\nProcessing sample: ", sn, "\n"))

  # Skip if already processed in a prior run
  rds_path <- file.path(rds_dir, paste0(sn, ".rds"))
  if (file.exists(rds_path)) {
    cat("  Already processed (RDS exists) — skipping\n")
    successfully_loaded <- c(successfully_loaded, sn)
    next
  }

  sm <- sample_manifest[[sn]]
  files <- sm$files

  # Force garbage collection before each sample
  invisible(gc())

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

    # Free counts immediately
    rm(counts); invisible(gc())

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
        bcs_with_pos <- colnames(obj)[colnames(obj) %in% rownames(pos)]
        obj <- subset(obj, cells = bcs_with_pos)
        pos_matched <- pos[colnames(obj), ]
        obj$in_tissue <- pos_matched$in_tissue
        obj$row      <- pos_matched$array_row
        obj$col      <- pos_matched$array_col
        obj$imagerow <- pos_matched$pxl_row
        obj$imagecol <- pos_matched$pxl_col
        cat(paste0("  Positions loaded: ", nrow(pos_matched), " in-tissue spots\n"))
        rm(pos, pos_matched, bcs_with_pos); invisible(gc())
      }
    } else {
      cat("  WARNING: tissue positions file not found\n")
    }

    # Add metadata BEFORE normalization
    obj$sample_id <- sn
    obj$condition <- extract_condition(sn)

    # JoinLayers on RNA assay (no image attached so safe)
    DefaultAssay(obj) <- "RNA"
    tryCatch({
      obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
      cat("  RNA layers joined\n")
    }, error = function(e) {
      cat(paste0("  Note: JoinLayers skipped (", conditionMessage(e), ")\n"))
    })

    # SCTransform normalization with log-norm fallback, with ncells
    # subsampling to reduce peak memory usage
    cat("  Running SCTransform...\n")
    n_cells <- ncol(obj)
    obj <- tryCatch({
      # Subsample for SCT model fitting if large to save memory
      sct_ncells <- min(3000, n_cells)
      SCTransform(obj, verbose = FALSE, ncells = sct_ncells)
    }, error = function(e) {
      cat(paste0("  SCTransform failed (", conditionMessage(e),
                 ") — falling back to log-normalization\n"))
      DefaultAssay(obj) <- "RNA"
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
      obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
      obj
    })

    # Ensure RNA assay is normalized for module scoring
    DefaultAssay(obj) <- "RNA"
    tryCatch({
      obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
    }, error = function(e) {
      cat(paste0("  Note: RNA NormalizeData skipped (", conditionMessage(e), ")\n"))
    })

    # Module scores on RNA assay with individual error handling
    all_genes <- rownames(obj[["RNA"]])

    aspc_use  <- intersect(ASPC_genes,         all_genes)
    ar_use    <- intersect(Hallmark_AR_genes,  all_genes)
    nepc_up   <- intersect(NEPC_Beltran_UP,    all_genes)
    nepc_dn   <- intersect(NEPC_Beltran_DOWN,  all_genes)

    if (length(aspc_use) > 0)
      obj <- tryCatch(
        AddModuleScore(obj, features = list(aspc_use), name = "ASPC_score",
                       ctrl = min(5, length(aspc_use)), assay = "RNA"),
        error = function(e) { cat(paste0("  ASPC score failed: ", conditionMessage(e), "\n")); obj })
    if (length(ar_use) > 0)
      obj <- tryCatch(
        AddModuleScore(obj, features = list(ar_use), name = "AR_score",
                       ctrl = min(5, length(ar_use)), assay = "RNA"),
        error = function(e) { cat(paste0("  AR score failed: ", conditionMessage(e), "\n")); obj })
    if (length(nepc_up) > 0)
      obj <- tryCatch(
        AddModuleScore(obj, features = list(nepc_up), name = "NEPC_UP_score",
                       ctrl = min(5, length(nepc_up)), assay = "RNA"),
        error = function(e) { cat(paste0("  NEPC_UP score failed: ", conditionMessage(e), "\n")); obj })
    if (length(nepc_dn) > 0)
      obj <- tryCatch(
        AddModuleScore(obj, features = list(nepc_dn), name = "NEPC_DOWN_score",
                       ctrl = min(5, length(nepc_dn)), assay = "RNA"),
        error = function(e) { cat(paste0("  NEPC_DOWN score failed: ", conditionMessage(e), "\n")); obj })

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

    # Build VisiumV1 image AFTER normalization, attach, then save
    if (!is.na(files$scalefactors) && file.exists(files$scalefactors)) {
      tryCatch({
        sf <- fromJSON(files$scalefactors)
        if (!is.na(files$lowres_image) && file.exists(files$lowres_image)) {
          img_data <- readPNG(files$lowres_image)
          in_tissue_vals <- if ("in_tissue" %in% colnames(obj@meta.data)) {
            obj$in_tissue
          } else {
            rep(1, ncol(obj))
          }
          tc <- data.frame(
            tissue   = in_tissue_vals,
            row      = if ("row" %in% colnames(obj@meta.data)) obj$row else rep(0, ncol(obj)),
            col      = if ("col" %in% colnames(obj@meta.data)) obj$col else rep(0, ncol(obj)),
            imagerow = if ("imagerow" %in% colnames(obj@meta.data)) obj$imagerow else rep(0, ncol(obj)),
            imagecol = if ("imagecol" %in% colnames(obj@meta.data)) obj$imagecol else rep(0, ncol(obj)),
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
          rm(img_data, tc); invisible(gc())
        }
      }, error = function(e) {
        cat(paste0("  WARNING: could not read spatial files — ", conditionMessage(e), "\n"))
      })
    }

    # Strip SCT assay before saving to reduce RDS size (can be
    # recomputed if needed; RNA + metadata is sufficient for all tiers)
    if ("SCT" %in% Assays(obj)) {
      DefaultAssay(obj) <- "RNA"
      obj[["SCT"]] <- NULL
    }

    # Save to disk and clear from memory
    saveRDS(obj, rds_path)
    cat(paste0("  Saved to disk: ", ncol(obj), " spots\n"))
    successfully_loaded <- c(successfully_loaded, sn)

    rm(obj); invisible(gc())

  }, error = function(e) {
    cat(paste0("  ERROR loading ", sn, ": ", conditionMessage(e), "\n"))
    # Clean up on error
    if (exists("obj", inherits = FALSE)) rm(obj)
    if (exists("counts", inherits = FALSE)) rm(counts)
    invisible(gc())
  })
}

cat(paste0("\nSuccessfully loaded ", length(successfully_loaded), " samples.\n"))

# ==============================================================================
# SECTION 3: TIER 1 — Spatial Annotation (Interface Zone Classification)
# Saves each sample as individual RDS + defines helpers for downstream tiers
# ==============================================================================
cat("\n=== TIER 1: Spatial Annotation ===\n")
tier1_status <- "SKIPPED"

tryCatch({
  tier1_dir <- file.path(visium_root, "TIER1_Interface")
  dir.create(tier1_dir, showWarnings = FALSE, recursive = TRUE)

  # ================================================================
  # Create RDS directory for per-sample objects (memory-efficient)
  # ================================================================
  rds_dir <- file.path(visium_root, "sample_rds")
  dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

  # ================================================================
  # Define helper functions used by Tier 1-7
  # ================================================================

  # Load a single sample from RDS
  load_sample <- function(sample_name) {
    rds_path <- file.path(rds_dir, paste0(sample_name, ".rds"))
    if (!file.exists(rds_path)) return(NULL)
    tryCatch(readRDS(rds_path), error = function(e) {
      cat(paste0("    Error loading ", sample_name, ": ",
                 conditionMessage(e), "\n"))
      NULL
    })
  }

  # Safe GetAssayData — handles Seurat v5 (layers) and v4 (slots)
  safe_GetAssayData <- function(obj, layer = "data") {
    tryCatch({
      # Try SCT first, fall back to RNA
      for (assay_name in c("SCT", "RNA")) {
        if (assay_name %in% Assays(obj)) {
          DefaultAssay(obj) <- assay_name
          mat <- tryCatch(
            GetAssayData(obj, layer = layer),
            error = function(e) {
              tryCatch(
                GetAssayData(obj, slot = layer),
                error = function(e2) NULL
              )
            }
          )
          if (!is.null(mat) && ncol(mat) > 0 && nrow(mat) > 0)
            return(mat)
        }
      }
      NULL
    }, error = function(e) NULL)
  }

  # ================================================================
  # Process each sample: score, classify zones, save RDS
  # ================================================================
  zone_prop_list <- list()
  successfully_loaded <- character(0)

  for (sn in names(seurat_list)) {
    cat(paste0("  TIER1 — ", sn, "\n"))

    tryCatch({
      obj <- seurat_list[[sn]]
      md  <- obj@meta.data

      # --- Verify score columns exist ---
      score_cols <- c("ASPC_score", "AR_score", "NEPC_UP_score")
      missing_sc <- setdiff(score_cols, colnames(md))
      if (length(missing_sc) > 0) {
        cat(paste0("    Missing: ",
                   paste(missing_sc, collapse = ", "), " — skipping\n"))
        next
      }

      # --- Compute medians for thresholding ---
      med_aspc <- median(md$ASPC_score,    na.rm = TRUE)
      med_ar   <- median(md$AR_score,      na.rm = TRUE)
      med_nepc <- median(md$NEPC_UP_score, na.rm = TRUE)

      # --- Classify interface zones ---
      # Interface: high ASPC AND (high NEPC OR high AR)
      is_interface <- (md$ASPC_score > med_aspc) &
        (md$NEPC_UP_score > med_nepc | md$AR_score > med_ar)
      is_aspc_high <- (md$ASPC_score > med_aspc) & !is_interface
      is_ar_high   <- (md$AR_score > med_ar) &
        !is_interface & !is_aspc_high
      is_nepc_high <- (md$NEPC_UP_score > med_nepc) &
        !is_interface & !is_aspc_high & !is_ar_high

      zone <- rep("Other", nrow(md))
      zone[is_nepc_high] <- "NEPC_high"
      zone[is_ar_high]   <- "AR_high"
      zone[is_aspc_high] <- "ASPC_high"
      zone[is_interface] <- "Interface"

      obj$interface_zone <- zone

      # --- Zone proportions ---
      tbl <- table(zone)
      zp <- as.data.frame(tbl, stringsAsFactors = FALSE)
      colnames(zp) <- c("zone", "count")
      zp$sample    <- sn
      zp$condition <- extract_condition(sn)
      zone_prop_list[[sn]] <- zp

      cat(paste0("    Zones: ",
                 paste(paste0(names(tbl), "=", tbl), collapse = ", "),
                 "\n"))

      # --- Spatial zone plot ---
      if (all(c("col", "row") %in% colnames(md))) {
        tryCatch({
          plot_df <- data.frame(
            x    = md$col,
            y    = -md$row,
            zone = zone,
            stringsAsFactors = FALSE
          )
          spot_sz <- compute_spot_size(plot_df[, c("x", "y")])
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
                 x = "Array Col", y = "Array Row (inv)",
                 color = "Zone") +
            theme_minimal(base_size = 10) +
            coord_fixed()
          ggsave(file.path(tier1_dir,
                           paste0(sn, "_interface_zones.pdf")),
                 plot = p, width = 10, height = 8)
        }, error = function(e) {
          cat(paste0("    Plot error: ",
                     conditionMessage(e), "\n"))
        })
      }

      # --- Score distribution plots per sample ---
      tryCatch({
        score_df <- data.frame(
          ASPC    = md$ASPC_score,
          AR      = md$AR_score,
          NEPC_UP = md$NEPC_UP_score,
          zone    = zone,
          stringsAsFactors = FALSE
        )
        p_scores <- ggplot(
          tidyr::pivot_longer(score_df, cols = c("ASPC","AR","NEPC_UP"),
                              names_to = "score_type",
                              values_to = "value"),
          aes(x = zone, y = value, fill = zone)) +
          geom_boxplot(outlier.size = 0.5) +
          facet_wrap(~score_type, scales = "free_y") +
          scale_fill_manual(values = c(
            Interface = "#E31A1C", ASPC_high = "#1F78B4",
            AR_high = "#33A02C", NEPC_high = "#FF7F00",
            Other = "#CCCCCC")) +
          labs(title = paste0(sn, " — Score Distributions by Zone"),
               x = "Zone", y = "Score") +
          theme_minimal(base_size = 9) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        ggsave(file.path(tier1_dir,
                         paste0(sn, "_score_distributions.pdf")),
               plot = p_scores, width = 12, height = 5)
      }, error = function(e) NULL)

      # --- Save individual RDS ---
      saveRDS(obj, file.path(rds_dir, paste0(sn, ".rds")))
      successfully_loaded <- c(successfully_loaded, sn)

      # Update in seurat_list too
      seurat_list[[sn]] <- obj

    }, error = function(e) {
      cat(paste0("    ERROR processing ", sn, ": ",
                 conditionMessage(e), "\n"))
    })
  }

  cat(paste0("\n  Samples processed: ",
             length(successfully_loaded), "/",
             length(names(seurat_list)), "\n"))

  # ================================================================
  # Aggregate zone proportions across all samples
  # ================================================================
  if (length(zone_prop_list) > 0) {
    zone_prop_df <- do.call(rbind, zone_prop_list)
    write.csv(zone_prop_df,
              file.path(tier1_dir,
                        "zone_proportions_per_sample.csv"),
              row.names = FALSE)

    # By condition
    cond_zone <- zone_prop_df %>%
      group_by(condition, zone) %>%
      summarise(total_spots = sum(count), .groups = "drop") %>%
      group_by(condition) %>%
      mutate(fraction = total_spots / sum(total_spots)) %>%
      ungroup()
    write.csv(cond_zone,
              file.path(tier1_dir,
                        "zone_proportions_by_condition.csv"),
              row.names = FALSE)

    cat(paste0("  Zone proportions saved (",
               nrow(zone_prop_df), " rows)\n"))

    # ================================================================
    # Summary PDF with all conditions
    # ================================================================
    tryCatch({
      pdf(file.path(tier1_dir, "TIER1_summary.pdf"),
          width = 12, height = 8)

      # --- Page 1: Stacked bar of zone proportions by condition ---
      tryCatch({
        p_bar <- ggplot(cond_zone,
                        aes(x = condition, y = fraction,
                            fill = zone)) +
          geom_bar(stat = "identity", position = "stack") +
          scale_fill_manual(values = c(
            Interface = "#E31A1C", ASPC_high = "#1F78B4",
            AR_high = "#33A02C", NEPC_high = "#FF7F00",
            Other = "#CCCCCC")) +
          labs(title = "Zone Proportions by Condition",
               x = "Condition", y = "Fraction of Spots",
               fill = "Zone") +
          theme_minimal(base_size = 12)
        print(p_bar)
      }, error = function(e) NULL)

      # --- Page 2: Total spot counts by zone per condition ---
      tryCatch({
        p_counts <- ggplot(cond_zone,
                           aes(x = condition, y = total_spots,
                               fill = zone)) +
          geom_bar(stat = "identity", position = "dodge") +
          scale_fill_manual(values = c(
            Interface = "#E31A1C", ASPC_high = "#1F78B4",
            AR_high = "#33A02C", NEPC_high = "#FF7F00",
            Other = "#CCCCCC")) +
          labs(title = "Total Spots per Zone by Condition",
               x = "Condition", y = "Spot Count",
               fill = "Zone") +
          theme_minimal(base_size = 12)
        print(p_counts)
      }, error = function(e) NULL)

      # --- Page 3: Per-sample zone proportions heatmap ---
      tryCatch({
        # Pivot to matrix
        samp_zone_wide <- zone_prop_df %>%
          group_by(sample) %>%
          mutate(frac = count / sum(count)) %>%
          ungroup() %>%
          tidyr::pivot_wider(id_cols = sample,
                             names_from = zone,
                             values_from = frac,
                             values_fill = 0)
        mat <- as.matrix(samp_zone_wide[, -1])
        rownames(mat) <- samp_zone_wide$sample

        if (have_pheatmap && nrow(mat) >= 2 && ncol(mat) >= 2) {
          # Condition annotation
          cond_annot <- data.frame(
            Condition = sapply(rownames(mat), extract_condition),
            row.names = rownames(mat))
          cond_colors <- list(
            Condition = c(BPH = "#66C2A5", TRNA = "#FC8D62",
                          NEADT = "#8DA0CB", CRPC = "#E78AC3",
                          MET = "#A6D854"))
          pheatmap(mat, cluster_cols = FALSE,
                   annotation_row = cond_annot,
                   annotation_colors = cond_colors,
                   main = "Zone Fraction per Sample",
                   fontsize_row = 5, fontsize_col = 8,
                   color = colorRampPalette(
                     c("white", "steelblue", "darkblue"))(100))
        }
      }, error = function(e) NULL)

      # --- Page 4: Interface fraction by condition (boxplot) ---
      tryCatch({
        iface_frac <- zone_prop_df %>%
          group_by(sample) %>%
          mutate(frac = count / sum(count)) %>%
          ungroup() %>%
          filter(zone == "Interface")
        if (nrow(iface_frac) > 0) {
          p_iface <- ggplot(iface_frac,
                            aes(x = condition, y = frac,
                                fill = condition)) +
            geom_boxplot(outlier.shape = 16) +
            geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
            labs(title = "Interface Fraction by Condition",
                 x = "Condition",
                 y = "Fraction of Spots = Interface") +
            theme_minimal(base_size = 12) +
            theme(legend.position = "none")
          print(p_iface)
        }
      }, error = function(e) NULL)

      # --- Pages 5+: Score distributions by zone (one per condition) ---
      for (cond_name in unique(sapply(successfully_loaded,
                                      extract_condition))) {
        tryCatch({
          cond_samps <- successfully_loaded[
            sapply(successfully_loaded, extract_condition) == cond_name]
          score_all <- list()
          for (sn in cond_samps) {
            obj <- load_sample(sn)
            if (is.null(obj)) next
            md <- obj@meta.data
            if (!all(c("ASPC_score", "AR_score",
                       "NEPC_UP_score", "interface_zone") %in%
                     colnames(md))) {
              rm(obj); next
            }
            score_all[[sn]] <- data.frame(
              ASPC    = md$ASPC_score,
              AR      = md$AR_score,
              NEPC_UP = md$NEPC_UP_score,
              zone    = md$interface_zone,
              sample  = sn,
              stringsAsFactors = FALSE)
            rm(obj); invisible(gc())
          }
          if (length(score_all) > 0) {
            sd_all <- do.call(rbind, score_all)
            sd_long <- tidyr::pivot_longer(
              sd_all, cols = c("ASPC", "AR", "NEPC_UP"),
              names_to = "score_type",
              values_to = "value")
            p_sc <- ggplot(sd_long,
                           aes(x = zone, y = value, fill = zone)) +
              geom_violin(scale = "width", alpha = 0.7) +
              geom_boxplot(width = 0.15, outlier.size = 0.3) +
              facet_wrap(~score_type, scales = "free_y") +
              scale_fill_manual(values = c(
                Interface = "#E31A1C", ASPC_high = "#1F78B4",
                AR_high = "#33A02C", NEPC_high = "#FF7F00",
                Other = "#CCCCCC")) +
              labs(title = paste0(cond_name,
                                  " — Score Distributions by Zone"),
                   x = "Zone", y = "Score") +
              theme_minimal(base_size = 10) +
              theme(axis.text.x = element_text(angle = 45,
                                               hjust = 1))
            print(p_sc)
            rm(sd_all, sd_long, score_all)
          }
        }, error = function(e) NULL)
      }

      dev.off()
      cat("  TIER1_summary.pdf saved\n")
    }, error = function(e) {
      cat(paste0("  Summary PDF error: ",
                 conditionMessage(e), "\n"))
      tryCatch(dev.off(), error = function(e2) NULL)
    })
  }

  # ================================================================
  # Print summary
  # ================================================================
  total_iface <- 0
  total_spots <- 0
  for (sn in successfully_loaded) {
    obj <- load_sample(sn)
    if (!is.null(obj) &&
        "interface_zone" %in% colnames(obj@meta.data)) {
      total_iface <- total_iface +
        sum(obj$interface_zone == "Interface", na.rm = TRUE)
      total_spots <- total_spots + ncol(obj)
    }
    rm(obj); invisible(gc())
  }
  cat(paste0("  Total spots: ", total_spots,
             " | Interface spots: ", total_iface,
             " (", round(100 * total_iface / max(total_spots, 1), 1),
             "%)\n"))

  tier1_status <- "COMPLETE"
  cat("TIER 1 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 1 ERROR: ", conditionMessage(e), "\n"))
  tier1_status <<- paste0("ERROR: ", conditionMessage(e))
})

# ==============================================================================
# TIER 2 PART 1 — Setup, Gene Sets, Distance Computation, Bin Assignment
# Source this file first, then source TIER2_Part2.R in the same session.
# ==============================================================================
cat("\n=== TIER 2: Distance-Binned Pseudo-Bulk DE + GSEA ===\n")
tier2_status <- "SKIPPED"

tryCatch({
  tier2_dir <- file.path(visium_root, "TIER2_DE")
  dir.create(tier2_dir, showWarnings = FALSE, recursive = TRUE)

  # ================================================================
  # Helpers
  # ================================================================
  if (!exists("rds_dir"))
    rds_dir <- file.path(visium_root, "sample_rds")
  dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

  if (!exists("load_sample", mode = "function")) {
    load_sample <- function(sample_name) {
      p <- file.path(rds_dir, paste0(sample_name, ".rds"))
      if (!file.exists(p)) return(NULL)
      tryCatch(readRDS(p), error = function(e) NULL)
    }
  }

  if (!exists("safe_GetAssayData", mode = "function")) {
    safe_GetAssayData <- function(obj, layer = "data") {
      tryCatch({
        for (an in c("SCT", "RNA")) {
          if (an %in% Assays(obj)) {
            DefaultAssay(obj) <- an
            mat <- tryCatch(
              GetAssayData(obj, layer = layer),
              error = function(e) {
                tryCatch(GetAssayData(obj, slot = layer),
                         error = function(e2) NULL)
              })
            if (!is.null(mat) && ncol(mat) > 0 && nrow(mat) > 0)
              return(mat)
          }
        }
        NULL
      }, error = function(e) NULL)
    }
  }

  # ================================================================
  # Determine sample list
  # ================================================================
  if (!exists("successfully_loaded") ||
      length(successfully_loaded) == 0) {
    rds_files <- list.files(rds_dir, pattern = "\\.rds$",
                            full.names = FALSE)
    if (length(rds_files) > 0) {
      successfully_loaded <- sub("\\.rds$", "", rds_files)
      cat(paste0("  Found ", length(successfully_loaded),
                 " RDS files\n"))
    } else if (exists("seurat_list") && length(seurat_list) > 0) {
      cat("  Saving seurat_list to RDS...\n")
      successfully_loaded <- character(0)
      for (sn in names(seurat_list)) {
        tryCatch({
          saveRDS(seurat_list[[sn]],
                  file.path(rds_dir, paste0(sn, ".rds")))
          successfully_loaded <- c(successfully_loaded, sn)
        }, error = function(e) NULL)
      }
    } else {
      stop("No samples available")
    }
  }

  cat("  Validating samples...\n")
  valid_samples <- character(0)
  for (sn in successfully_loaded) {
    obj <- load_sample(sn)
    if (is.null(obj)) next
    md <- obj@meta.data
    ok <- "interface_zone" %in% colnames(md) &&
      all(c("row", "col") %in% colnames(md)) &&
      sum(md$interface_zone == "Interface", na.rm = TRUE) > 0
    if (ok) valid_samples <- c(valid_samples, sn)
    else cat(paste0("    ", sn, ": skipped\n"))
    rm(obj); invisible(gc())
  }
  cat(paste0("  Valid: ", length(valid_samples), "\n"))
  if (length(valid_samples) == 0)
    stop("No valid samples with interface_zone")

  # ================================================================
  # Install/load enrichment packages
  # ================================================================
  for (pkg in c("msigdbr", "fgsea")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      tryCatch({
        if (pkg == "fgsea")
          BiocManager::install("fgsea", ask = FALSE, update = FALSE)
        else install.packages(pkg)
      }, error = function(e) NULL)
  }
  have_msigdbr <- requireNamespace("msigdbr", quietly = TRUE)
  have_fgsea   <- requireNamespace("fgsea",   quietly = TRUE)
  if (have_msigdbr) library(msigdbr)
  if (have_fgsea)   library(fgsea)
  cat(paste0("  msigdbr: ", have_msigdbr,
             " | fgsea: ", have_fgsea, "\n"))

  # ================================================================
  # MSigDB toggles
  # ================================================================
  MSIGDB_USE_H  <- TRUE
  MSIGDB_USE_C2 <- TRUE
  MSIGDB_USE_C5 <- TRUE
  MSIGDB_USE_C6 <- TRUE
  MSIGDB_USE_C7 <- TRUE
  MIN_GS <- 15
  MAX_GS <- 500

  # ================================================================
  # MSigDB loader helpers
  # ================================================================
  if (have_msigdbr) {
    .msig_args <- names(formals(msigdbr::msigdbr))
    .use_new   <- "collection" %in% .msig_args
  }

  load_msigdb_collection <- function(cat_val, subcat_val = NULL,
                                      label = cat_val,
                                      min_size = MIN_GS,
                                      max_size = MAX_GS) {
    if (!have_msigdbr) return(list())
    tryCatch({
      if (.use_new) {
        if (!is.null(subcat_val))
          df <- msigdbr(species = "Homo sapiens",
                        collection = cat_val,
                        subcollection = subcat_val)
        else
          df <- msigdbr(species = "Homo sapiens",
                        collection = cat_val)
      } else {
        if (!is.null(subcat_val))
          df <- msigdbr(species = "Homo sapiens",
                        category = cat_val,
                        subcategory = subcat_val)
        else
          df <- msigdbr(species = "Homo sapiens",
                        category = cat_val)
      }
      gs <- split(df$gene_symbol, df$gs_name)
      sz <- sapply(gs, length)
      gs <- gs[sz >= min_size & sz <= max_size]
      cat(paste0("      ", label, ": ",
                 length(gs), " sets\n"))
      gs
    }, error = function(e) {
      cat(paste0("      ", label, " error: ",
                 conditionMessage(e), "\n"))
      list()
    })
  }

  load_by_pattern <- function(cat_val, pattern_list,
                               min_size = MIN_GS,
                               max_size = MAX_GS) {
    if (!have_msigdbr) return(list())
    result <- list()
    for (pl in pattern_list) {
      tryCatch({
        sub_gs <- load_msigdb_collection(
          cat_val, subcat_val = pl$pattern,
          label = pl$label,
          min_size = min_size, max_size = max_size)
        result <- c(result, sub_gs)
      }, error = function(e) {
        tryCatch({
          if (.use_new)
            df <- msigdbr(species = "Homo sapiens",
                          collection = cat_val)
          else
            df <- msigdbr(species = "Homo sapiens",
                          category = cat_val)
          gs_all <- split(df$gene_symbol, df$gs_name)
          matched <- gs_all[grepl(pl$pattern, names(gs_all),
                                  ignore.case = TRUE)]
          sz <- sapply(matched, length)
          matched <- matched[sz >= min_size & sz <= max_size]
          cat(paste0("      ", pl$label, " (grep): ",
                     length(matched), "\n"))
          result <<- c(result, matched)
        }, error = function(e2) NULL)
      })
    }
    result
  }

  # ================================================================
  # Build gene set collection
  # ================================================================
  cat("  Building gene sets...\n")
  gene_set_list <- list()

  # ==============================================================
  # A) CUSTOM SIGNATURES — 7 original + 7 cell death + 4 mito
  # ==============================================================
  cat("    Custom signatures:\n")

  ASPC_SIGNATURE <- if (exists("ASPC_genes")) ASPC_genes else c()
  NEPC_BELTRAN_UP <- if (exists("NEPC_Beltran_UP")) NEPC_Beltran_UP else c()
  NEPC_BELTRAN_DOWN <- if (exists("NEPC_Beltran_DOWN")) NEPC_Beltran_DOWN else c()
  HALLMARK_AR_SIGNATURE <- if (exists("Hallmark_AR_genes")) Hallmark_AR_genes else c()

  if (!exists("TNFSF12_SIGNATURE"))
    TNFSF12_SIGNATURE <- c(
      "TNFSF12","TNFRSF12A","TNFRSF25","TRAF2","TRAF5",
      "NFKB1","RELA","MAP3K14","BIRC2","BIRC3",
      "RIPK1","CFLAR","CASP8","MMP9","CXCL8",
      "CCL2","ICAM1","VCAM1","VEGFA","FGF2",
      "ANGPT2","IL6","STAT3","AKT1","MAPK1")

  if (!exists("MDK_SIGNATURE"))
    MDK_SIGNATURE <- c(
      "MDK","PTN","SDC1","SDC3","SDC4",
      "GPC1","GPC2","GPC3","GPC4","GPC6",
      "PTPRZ1","ALK","LRP1","ITGA4","ITGA6",
      "ITGB1","NCL","NOTCH2","NTRK2","NTRK3",
      "FGFR1","VEGFA","STAT3","AKT1","MAPK1")

  if (!exists("HALLMARK_ADIPOGENESIS_CUSTOM"))
    HALLMARK_ADIPOGENESIS_CUSTOM <- c(
      "PPARG","CEBPA","CEBPB","CEBPD","FABP4",
      "ADIPOQ","LEP","LPL","PLIN1","PLIN2",
      "FASN","SCD","DGAT1","DGAT2","ACACA",
      "ACLY","LIPE","PNPLA2","CD36","SLC27A1",
      "GPD1","PCK1","CIDEC","CIDEA","RETN",
      "CFD","ADIPOR1","ADIPOR2","IRS1","IRS2",
      "AQP7","SLC2A4","GYK","ACSL1","ACSL5",
      "ELOVL6","MGLL","AKR1B1","SCD5","THRSP")

  # --- Cell death signatures ---
  NECROPTOSIS_SIG <- c(
    "RIPK1","RIPK3","MLKL","TNFRSF1A","TNF",
    "TRADD","FADD","CASP8","CFLAR","CYLD",
    "BIRC2","BIRC3","XIAP","TRAF2","TRAF5",
    "TAB2","TAB3","MAP3K7","CHUK","IKBKB",
    "NFKBIA","ZBP1","DAI","ADAR","IFNAR1",
    "IFNAR2","JAK1","TYK2","STAT1","STAT2",
    "IRF9","USP21","SPATA2","HOIP","SHARPIN",
    "OTULIN","A20","TNFAIP3","LUBAC","RNF31")

  FERROPTOSIS_SIG <- c(
    "GPX4","SLC7A11","SLC3A2","ACSL4","LPCAT3",
    "TFRC","FTH1","FTL","HMOX1","NRF2","NFE2L2",
    "SLC40A1","NCOA4","IREB2","CISD1","VDAC2",
    "VDAC3","RSL3","ALOX15","ALOX12","ALOX5",
    "POR","CBS","FSP1","AIFM2","DHODH",
    "GCH1","BH4","GCLC","GCLM","GSS",
    "GSR","SLC11A2","DMT1","STEAP3","PCBP1",
    "PCBP2","PROM2","ATG5","ATG7","BECN1",
    "HSPA5","HSPB1","FANCD2","TP53","SAT1")

  PYROPTOSIS_SIG <- c(
    "NLRP3","NLRC4","NLRP1","AIM2","PYCARD",
    "CASP1","CASP4","CASP5","CASP11","GSDMD",
    "GSDME","GSDMA","GSDMB","GSDMC","IL1B",
    "IL18","IL33","HMGB1","P2RX7","PANX1",
    "NEK7","TXNIP","BRCC3","TLR4","MYD88",
    "IRAK1","IRAK4","TRAF6","NFKB1","RELA",
    "NLRP6","NLRP7","MEFV","CARD8","NAIP")

  CUPROPTOSIS_SIG <- c(
    "FDX1","LIAS","LIPT1","DLD","DLAT",
    "PDHA1","PDHB","MTF1","GLS","CDKN2A",
    "SLC31A1","ATP7A","ATP7B","SOD1","CCS",
    "MT1A","MT1B","MT1E","MT1F","MT1G",
    "MT1H","MT1M","MT1X","MT2A","COX11",
    "COX17","SCO1","SCO2","COA6","ATOX1",
    "COMMD1","XIAP","PARK7","CP","STEAP4",
    "LOX","LOXL2","SLC25A3","DBT","GCSH")

  APOPTOSIS_INTRINSIC_SIG <- c(
    "BAX","BAK1","BCL2","BCL2L1","MCL1",
    "BID","BIM","BCL2L11","BAD","PUMA",
    "BBC3","NOXA","PMAIP1","CYCS","APAF1",
    "CASP9","CASP3","CASP7","SMAC","DIABLO",
    "XIAP","BIRC5","AIF","AIFM1","ENDOG",
    "TP53","MDM2","CDKN1A","FAS","FASLG",
    "CASP8","CASP10","TRADD","RIPK1","CFLAR",
    "PARP1","PARP2","DFFA","DFFB","ICAD")

  APOPTOSIS_EXTRINSIC_SIG <- c(
    "FAS","FASLG","TNFRSF10A","TNFRSF10B","TRAIL",
    "TNFSF10","TNF","TNFRSF1A","TNFRSF1B","TRADD",
    "FADD","CASP8","CASP10","CFLAR","BID",
    "RIPK1","TRAF2","CASP3","CASP7","XIAP",
    "BIRC2","BIRC3","DIABLO","SMAC","CYCS",
    "APAF1","CASP9","BCL2","BCL2L1","MCL1")

  AUTOPHAGY_SIG <- c(
    "BECN1","ATG5","ATG7","ATG12","ATG16L1",
    "ATG3","ATG4B","ATG9A","ATG14","ATG101",
    "ULK1","ULK2","FIP200","RB1CC1","PIK3C3",
    "VPS34","PIK3R4","AMBRA1","UVRAG","LC3B",
    "MAP1LC3B","MAP1LC3A","GABARAP","GABARAPL1",
    "GABARAPL2","SQSTM1","NBR1","OPTN","NDP52",
    "CALCOCO2","BNIP3","BNIP3L","FUNDC1","PINK1",
    "PRKN","PARKIN","LAMP1","LAMP2","CTSD","TFEB")

  # --- Mitochondrial function signatures ---
  MITO_OXPHOS_SIG <- c(
    "MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L",
    "MT-ND5","MT-ND6","NDUFA1","NDUFA2","NDUFA4",
    "NDUFA5","NDUFA9","NDUFA10","NDUFA13","NDUFB3",
    "NDUFB5","NDUFB8","NDUFB10","NDUFC2","NDUFS1",
    "NDUFS2","NDUFS3","NDUFS4","NDUFS6","NDUFS7",
    "NDUFS8","NDUFV1","NDUFV2","SDHA","SDHB",
    "SDHC","SDHD","MT-CYB","UQCRC1","UQCRC2",
    "UQCRFS1","UQCRB","UQCRQ","MT-CO1","MT-CO2",
    "MT-CO3","COX4I1","COX5A","COX5B","COX6A1",
    "COX6B1","COX6C","COX7A2","COX7C","COX8A",
    "MT-ATP6","MT-ATP8","ATP5F1A","ATP5F1B","ATP5F1C",
    "ATP5F1D","ATP5F1E","ATP5MC1","ATP5MC2","ATP5MC3",
    "ATP5ME","ATP5MF","ATP5PB","ATP5PD","ATP5PO")

  MITO_DYNAMICS_SIG <- c(
    "MFN1","MFN2","OPA1","DRP1","DNM1L",
    "FIS1","MFF","MIEF1","MIEF2","MID49",
    "MID51","PINK1","PRKN","PARKIN","BNIP3",
    "BNIP3L","FUNDC1","PHB","PHB2","RHOT1",
    "RHOT2","TRAK1","TRAK2","KIF5B","MIRO1",
    "MIRO2","TOMM20","TOMM40","TOMM70","TIMM23",
    "TIMM44","TIMM50","TIMM13","SAMM50","MTX1",
    "MTX2","VDAC1","VDAC2","VDAC3","SLC25A4")

  MITO_BIOGENESIS_SIG <- c(
    "PPARGC1A","PGC1A","PPARGC1B","NRF1","GABPA",
    "TFAM","TFB1M","TFB2M","POLG","POLG2",
    "TWNK","SSBP1","MTERF1","MTERF2","MTERF3",
    "MTERF4","MRPL11","MRPL12","MRPL15","MRPL16",
    "MRPS12","MRPS15","MRPS22","MRPS25","MRPS27",
    "COX10","COX15","SURF1","SCO1","SCO2",
    "COA3","COA5","COA6","COA7","BCS1L",
    "UQCC1","UQCC2","NDUFAF1","NDUFAF2","ACAD9")

  MITO_APOPTOSIS_SIG <- c(
    "BAX","BAK1","BCL2","BCL2L1","MCL1",
    "BID","BIM","BCL2L11","BAD","PUMA",
    "BBC3","NOXA","PMAIP1","CYCS","APAF1",
    "CASP9","CASP3","DIABLO","SMAC","HTRA2",
    "AIF","AIFM1","ENDOG","VDAC1","VDAC2",
    "ANT1","SLC25A4","ANT2","SLC25A5","TSPO",
    "PPIF","CYPD","MOMP","BOK","BMF",
    "HRK","BIK","BNIP3","BNIP3L","AIFM2")

  # --- Additional biology ---
  SENESCENCE_SIG <- c(
    "CDKN1A","CDKN2A","CDKN2B","TP53","RB1",
    "SERPINE1","IGFBP3","IGFBP5","IGFBP7","MMP1",
    "MMP3","MMP10","IL1A","IL1B","IL6",
    "IL8","CXCL8","CXCL1","CXCL2","CXCL3",
    "CCL2","CCL3","CCL5","CCL20","VEGFA",
    "FGF2","HGF","AREG","EREG","MIF",
    "HMGA1","HMGA2","LMNB1","H2AFX","GLB1",
    "SA_BGAL","GADD45A","GADD45B","MDM2","BCL2L1")

  LIPID_METABOLISM_SIG <- c(
    "FASN","ACACA","ACACB","ACLY","SCD",
    "ELOVL5","ELOVL6","FADS1","FADS2","ACSL1",
    "ACSL3","ACSL4","ACSL5","CPT1A","CPT1B",
    "CPT2","HADHA","HADHB","ACADM","ACADL",
    "ACADVL","ECHS1","EHHADH","HSD17B4","ABCA1",
    "ABCG1","NPC1","NPC2","LIPA","LIPG",
    "DGAT1","DGAT2","AGPAT1","GPAT4","PLIN1",
    "PLIN2","PLIN3","CD36","FABP4","FABP5",
    "LPL","HMGCR","HMGCS1","SQLE","LDLR")

  NEURO_DIFF_SIG <- c(
    "SYP","SYN1","SNAP25","VAMP2","STX1A",
    "CHGA","CHGB","ENO2","NCAM1","TUBB3",
    "MAP2","MAPT","NEFH","NEFL","NEFM",
    "RBFOX3","NEUROD1","NEUROD2","NEUROG1","NEUROG2",
    "ASCL1","INSM1","POU3F2","SOX2","SOX11",
    "DLL3","HES6","SRRM4","REST","RCOR2",
    "EZH2","BRN2","AURKA","MYCN","DLG4",
    "GRIA2","GRIN1","SCG3","PCSK1","DNER")

  # === Assemble custom_sigs ===
  custom_sigs <- list(
    # Original 7
    "CUSTOM_ASPC_SIGNATURE"             = ASPC_SIGNATURE,
    "CUSTOM_NEPC_BELTRAN_UP"            = NEPC_BELTRAN_UP,
    "CUSTOM_NEPC_BELTRAN_DOWN"          = NEPC_BELTRAN_DOWN,
    "CUSTOM_TNFSF12_SIGNATURE"          = TNFSF12_SIGNATURE,
    "CUSTOM_MDK_SIGNATURE"              = MDK_SIGNATURE,
    "CUSTOM_HALLMARK_ADIPOGENESIS"      = HALLMARK_ADIPOGENESIS_CUSTOM,
    "CUSTOM_HALLMARK_ANDROGEN_RESPONSE" = HALLMARK_AR_SIGNATURE,
    # Cell death (7)
    "CUSTOM_NECROPTOSIS"                = NECROPTOSIS_SIG,
    "CUSTOM_FERROPTOSIS"                = FERROPTOSIS_SIG,
    "CUSTOM_PYROPTOSIS"                 = PYROPTOSIS_SIG,
    "CUSTOM_CUPROPTOSIS"                = CUPROPTOSIS_SIG,
    "CUSTOM_APOPTOSIS_INTRINSIC"        = APOPTOSIS_INTRINSIC_SIG,
    "CUSTOM_APOPTOSIS_EXTRINSIC"        = APOPTOSIS_EXTRINSIC_SIG,
    "CUSTOM_AUTOPHAGY"                  = AUTOPHAGY_SIG,
    # Mitochondrial (4)
    "CUSTOM_MITO_OXPHOS"                = MITO_OXPHOS_SIG,
    "CUSTOM_MITO_DYNAMICS"              = MITO_DYNAMICS_SIG,
    "CUSTOM_MITO_BIOGENESIS"            = MITO_BIOGENESIS_SIG,
    "CUSTOM_MITO_APOPTOSIS"             = MITO_APOPTOSIS_SIG,
    # Additional biology (3)
    "CUSTOM_SENESCENCE"                 = SENESCENCE_SIG,
    "CUSTOM_LIPID_METABOLISM"           = LIPID_METABOLISM_SIG,
    "CUSTOM_NEURO_DIFFERENTIATION"      = NEURO_DIFF_SIG
  )

  for (nm in names(custom_sigs)) {
    gene_set_list[[nm]] <- custom_sigs[[nm]]
    cat(paste0("      ", nm, ": ",
               length(custom_sigs[[nm]]), " genes\n"))
  }

  # Group custom sigs for organized plotting
  custom_groups <- list(
    "Original"     = c("CUSTOM_ASPC_SIGNATURE",
                       "CUSTOM_NEPC_BELTRAN_UP",
                       "CUSTOM_NEPC_BELTRAN_DOWN",
                       "CUSTOM_TNFSF12_SIGNATURE",
                       "CUSTOM_MDK_SIGNATURE",
                       "CUSTOM_HALLMARK_ADIPOGENESIS",
                       "CUSTOM_HALLMARK_ANDROGEN_RESPONSE"),
    "Cell_Death"   = c("CUSTOM_NECROPTOSIS",
                       "CUSTOM_FERROPTOSIS",
                       "CUSTOM_PYROPTOSIS",
                       "CUSTOM_CUPROPTOSIS",
                       "CUSTOM_APOPTOSIS_INTRINSIC",
                       "CUSTOM_APOPTOSIS_EXTRINSIC",
                       "CUSTOM_AUTOPHAGY"),
    "Mitochondria" = c("CUSTOM_MITO_OXPHOS",
                       "CUSTOM_MITO_DYNAMICS",
                       "CUSTOM_MITO_BIOGENESIS",
                       "CUSTOM_MITO_APOPTOSIS"),
    "Biology"      = c("CUSTOM_SENESCENCE",
                       "CUSTOM_LIPID_METABOLISM",
                       "CUSTOM_NEURO_DIFFERENTIATION")
  )

  # ==============================================================
  # B) MSigDB collections
  # ==============================================================
  cat("    MSigDB:\n")
  h_list <- list()
  if (MSIGDB_USE_H) {
    h_list <- load_msigdb_collection("H", label = "Hallmark",
                                      min_size = 5, max_size = 500)
    gene_set_list <- c(gene_set_list, h_list)
  }

  c2_cp_list <- list()
  if (MSIGDB_USE_C2) {
    c2_cp_list <- load_by_pattern("C2", list(
      list(pattern = "KEGG",         label = "KEGG"),
      list(pattern = "REACTOME",     label = "Reactome"),
      list(pattern = "WIKIPATHWAYS", label = "WikiPathways"),
      list(pattern = "PID",          label = "PID"),
      list(pattern = "BIOCARTA",     label = "BioCarta")))
    gene_set_list <- c(gene_set_list, c2_cp_list)
  }

  c5_go_list <- list()
  if (MSIGDB_USE_C5) {
    c5_go_list <- load_by_pattern("C5", list(
      list(pattern = "GO:BP", label = "GO BP"),
      list(pattern = "GO:MF", label = "GO MF"),
      list(pattern = "GO:CC", label = "GO CC")))
    gene_set_list <- c(gene_set_list, c5_go_list)
  }

  c6_onc_list <- list()
  if (MSIGDB_USE_C6) {
    c6_onc_list <- load_msigdb_collection("C6", label = "Oncogenic")
    gene_set_list <- c(gene_set_list, c6_onc_list)
  }

  c7_imm_list <- list()
  if (MSIGDB_USE_C7) {
    c7_imm_list <- load_msigdb_collection("C7", label = "Immunologic")
    gene_set_list <- c(gene_set_list, c7_imm_list)
  }

  cat(paste0("    Total gene sets: ",
             length(gene_set_list), "\n"))

  gs_source <- setNames(rep("OTHER", length(gene_set_list)),
                        names(gene_set_list))
  gs_source[grepl("^CUSTOM_", names(gs_source))] <- "CUSTOM"
  gs_source[names(gs_source) %in% names(h_list)]      <- "HALLMARK"
  gs_source[names(gs_source) %in% names(c2_cp_list)]   <- "C2_CP"
  gs_source[names(gs_source) %in% names(c5_go_list)]   <- "C5_GO"
  gs_source[names(gs_source) %in% names(c6_onc_list)]  <- "C6_ONCO"
  gs_source[names(gs_source) %in% names(c7_imm_list)]  <- "C7_IMMUNO"
  rm(h_list, c2_cp_list, c5_go_list, c6_onc_list, c7_imm_list)
  invisible(gc())

  # ================================================================
  # Phase 1: Compute distances
  # ================================================================
  cat("\n  Phase 1: Distances...\n")
  samples_with_dist <- character(0)
  all_distances <- list()

  for (sn in valid_samples) {
    obj <- load_sample(sn)
    if (is.null(obj)) next
    md <- obj@meta.data
    if (!"dist_to_interface" %in% colnames(md)) {
      iface_idx <- which(md$interface_zone == "Interface")
      if (length(iface_idx) == 0) {
        rm(obj); invisible(gc()); next
      }
      coords <- as.matrix(md[, c("row", "col")])
      ic <- coords[iface_idx, , drop = FALSE]
      ns <- nrow(coords); dti <- numeric(ns)
      for (ci in seq(1, ns, by = 1000)) {
        ce <- min(ci + 999, ns)
        ch <- coords[ci:ce, , drop = FALSE]
        if (nrow(ic) <= 5000) {
          rd <- outer(ch[,1], ic[,1], "-")
          cd <- outer(ch[,2], ic[,2], "-")
          dti[ci:ce] <- apply(sqrt(rd^2+cd^2), 1, min)
          rm(rd, cd)
        } else {
          for (j in seq_len(nrow(ch))) {
            dd <- sweep(ic, 2, ch[j,])
            dti[ci+j-1] <- min(sqrt(rowSums(dd^2)))
          }
        }
      }
      dti[iface_idx] <- 0
      obj$dist_to_interface <- dti
      saveRDS(obj, file.path(rds_dir, paste0(sn, ".rds")))
      cat(paste0("    ", sn, ": computed (",
                 sum(dti==0), " iface)\n"))
      rm(coords, ic)
    } else {
      cat(paste0("    ", sn, ": present\n"))
    }
    cn <- extract_condition(sn)
    if (is.null(all_distances[[cn]]))
      all_distances[[cn]] <- numeric(0)
    all_distances[[cn]] <- c(all_distances[[cn]],
                              obj$dist_to_interface)
    samples_with_dist <- c(samples_with_dist, sn)
    rm(obj); invisible(gc())
  }
  cat(paste0("  With dist: ",length(samples_with_dist),"\n"))
  if (length(samples_with_dist) == 0)
    stop("No samples with distance data")

  # Adaptive bins
  cat("\n  Adaptive bins:\n")
  cond_bin_info <- list()
  for (cn in names(all_distances)) {
    d <- all_distances[[cn]]; d_ni <- d[d > 0]
    if (length(d_ni)==0) { cond_bin_info[[cn]] <- NULL; next }
    d_max <- max(d_ni); d_med <- median(d_ni)
    d_q75 <- quantile(d_ni, 0.75)
    cb <- sort(unique(c(0, 1.5, 3, d_q75, d_max)))
    fb <- cb[1]
    for (i in 2:length(cb))
      if (cb[i]-fb[length(fb)] >= 0.3) fb <- c(fb, cb[i])
    if (length(fb)<3) fb <- sort(unique(c(0, d_med, d_max)))
    if (length(fb)<3) fb <- c(0, d_max)
    bl <- character(length(fb)-1)
    for (i in seq_along(bl))
      bl[i] <- paste0("dist_",round(fb[i],1),"_to_",round(fb[i+1],1))
    cond_bin_info[[cn]] <- list(breaks=fb, labels=bl)
    cat(paste0("    ",cn,": ",paste(bl,collapse=" | "),"\n"))
  }
  rm(all_distances); invisible(gc())

  # Assign bins
  cat("  Assigning bins...\n")
  for (sn in samples_with_dist) {
    obj <- load_sample(sn)
    if (is.null(obj)) next
    cn <- extract_condition(sn)
    bi <- cond_bin_info[[cn]]
    d <- obj$dist_to_interface
    res <- rep(NA_character_, length(d))
    res[d==0] <- "bin_0_interface"
    if (!is.null(bi) && length(bi$breaks)>=2) {
      ni <- which(d>0)
      if (length(ni)>0) {
        dn <- pmin(d[ni], max(bi$breaks))
        dn <- pmax(dn, min(bi$breaks[bi$breaks>0]))
        ct <- cut(dn, breaks=bi$breaks, include.lowest=TRUE,
                  labels=bi$labels, right=TRUE)
        res[ni] <- as.character(ct)
        na_i <- ni[is.na(res[ni])]
        if (length(na_i)>0) res[na_i] <- bi$labels[length(bi$labels)]
      }
    } else res[d>0] <- "bin_noninterface"
    obj$dist_bin_adaptive <- res
    saveRDS(obj, file.path(rds_dir, paste0(sn, ".rds")))
    rm(obj); invisible(gc())
  }

  cat("\n  Part 1 complete. Source TIER2_Part2.R now.\n")

    # ================================================================
  # TIER 2 PART 2 — DE + GSEA + Plots + PDFs
  # Continues the tryCatch from Part 1.
  # ================================================================
  cat("\n  Phase 2: Per-condition DE + GSEA...\n")
  conditions <- unique(sapply(samples_with_dist, extract_condition))

  for (cond in conditions) {
    cat(paste0("\n  ===== ", cond, " =====\n"))
    cond_samples <- samples_with_dist[
      sapply(samples_with_dist, extract_condition) == cond]
    if (length(cond_samples)==0) next

    bi <- cond_bin_info[[cond]]
    if (is.null(bi)) { cat("    No bins\n"); next }
    abl <- c("bin_0_interface", bi$labels)
    dbin <- bi$labels[length(bi$labels)]
    cat(paste0("    DE: bin_0_interface vs ", dbin, "\n"))

    # --- Pseudo-bulk ---
    bel <- list(); sgl <- list(); sok <- character(0)
    for (sn in cond_samples) {
      obj <- load_sample(sn)
      if (is.null(obj) ||
          !"dist_bin_adaptive" %in% colnames(obj@meta.data)) {
        if (!is.null(obj)) rm(obj); invisible(gc()); next
      }
      em <- safe_GetAssayData(obj, layer="data")
      if (is.null(em)) { rm(obj); invisible(gc()); next }
      sgl[[sn]] <- rownames(em); hd <- FALSE
      for (bn in abl) {
        cc <- colnames(obj)[!is.na(obj$dist_bin_adaptive) &
                              obj$dist_bin_adaptive==bn]
        cc <- cc[cc %in% colnames(em)]
        if (length(cc)==0) next
        if (is.null(bel[[bn]])) bel[[bn]] <- list()
        bel[[bn]][[sn]] <- rowMeans(em[,cc,drop=FALSE])
        hd <- TRUE
      }
      if (hd) sok <- c(sok, sn)
      rm(obj, em); invisible(gc())
    }
    if (length(sok)==0) next

    ug <- if (length(sgl)>=2) Reduce(intersect, sgl) else sgl[[1]]
    if (length(ug)<100) { rm(bel,sgl); invisible(gc()); next }

    col_b <- function(l, g) {
      if (length(l)==0) return(NULL)
      m <- do.call(cbind, lapply(l, function(v) v[g]))
      if (is.matrix(m)) rowMeans(m,na.rm=TRUE) else m
    }
    bm <- lapply(bel, function(l) col_b(l, ug))

    b0 <- bm[["bin_0_interface"]]; bf <- bm[[dbin]]
    if (is.null(bf))
      for (a in rev(bi$labels))
        if (!is.null(bm[[a]])) { dbin <- a; bf <- bm[[a]]; break }
    if (is.null(b0)||is.null(bf)) {
      rm(bel,bm,sgl); invisible(gc()); next
    }
    cg <- intersect(names(b0),names(bf))
    cg <- cg[!is.na(cg)]
    if (length(cg)<10) { rm(bel,bm,sgl); invisible(gc()); next }

    # --- DE ---
    de <- tryCatch({
      fc <- log2((b0[cg]+0.01)/(bf[cg]+0.01))
      v0L <- list(); vFL <- list()
      for (sn in sok) {
        obj <- load_sample(sn)
        if (is.null(obj) ||
            !"dist_bin_adaptive" %in% colnames(obj@meta.data)) {
          if (!is.null(obj)) rm(obj); invisible(gc()); next
        }
        em <- safe_GetAssayData(obj,layer="data")
        if (is.null(em)) { rm(obj); invisible(gc()); next }
        gh <- intersect(cg,rownames(em))
        if (length(gh)<10) { rm(obj,em); invisible(gc()); next }
        c0 <- colnames(obj)[!is.na(obj$dist_bin_adaptive) &
                              obj$dist_bin_adaptive=="bin_0_interface"]
        c0 <- c0[c0 %in% colnames(em)]
        cF <- colnames(obj)[!is.na(obj$dist_bin_adaptive) &
                              obj$dist_bin_adaptive==dbin]
        cF <- cF[cF %in% colnames(em)]
        if (length(c0)>0) {
          v <- rowMeans(em[gh,c0,drop=FALSE])
          fv <- setNames(rep(NA_real_,length(cg)),cg)
          fv[gh] <- v[gh]; v0L[[sn]] <- fv
        }
        if (length(cF)>0) {
          v <- rowMeans(em[gh,cF,drop=FALSE])
          fv <- setNames(rep(NA_real_,length(cg)),cg)
          fv[gh] <- v[gh]; vFL[[sn]] <- fv
        }
        rm(obj,em); invisible(gc())
      }
      v0 <- if(length(v0L)>0) do.call(cbind,v0L) else NULL
      vF <- if(length(vFL)>0) do.call(cbind,vFL) else NULL
      pv <- sapply(cg, function(g) {
        a2 <- if(!is.null(v0)){if(is.matrix(v0))as.numeric(v0[g,])
          else as.numeric(v0[g])} else numeric(0)
        b2 <- if(!is.null(vF)){if(is.matrix(vF))as.numeric(vF[g,])
          else as.numeric(vF[g])} else numeric(0)
        a2 <- a2[!is.na(a2)]; b2 <- b2[!is.na(b2)]
        if(length(a2)<2||length(b2)<2) return(1)
        tryCatch(wilcox.test(a2,b2,exact=FALSE)$p.value,
                 error=function(e)1)
      })
      rm(v0L,vFL,v0,vF); invisible(gc())
      data.frame(gene=cg, log2FC=fc[cg], pvalue=pv,
                 padj=p.adjust(pv,"BH"),
                 mean_iface=b0[cg], mean_distal=bf[cg],
                 stringsAsFactors=FALSE)
    }, error=function(e) {
      cat(paste0("    DE error: ",conditionMessage(e),"\n")); NULL
    })

    if (is.null(de)||nrow(de)==0) {
      rm(bel,bm,sgl); invisible(gc()); next
    }
    de$cscore <- -log10(de$padj+1e-300)*abs(de$log2FC)
    de$sig <- !is.na(de$padj) & de$padj<0.05 & abs(de$log2FC)>0.5
    de <- de[order(de$cscore,decreasing=TRUE),]
    tryCatch(write.csv(de,
      file.path(tier2_dir,paste0(cond,"_DE_interface_vs_distal.csv")),
      row.names=FALSE), error=function(e) NULL)
    nsig <- sum(de$sig,na.rm=TRUE)
    cat(paste0("    DE: ",nrow(de)," genes, ",nsig," sig\n"))
    sgr <- de[de$sig,]
    t15 <- head(sgr$gene,15)
    t50 <- head(sgr$gene,50)
    if(length(t50)<50)
      t50 <- c(t50,head(setdiff(de$gene,t50),50-length(t50)))
    t10 <- head(if(nrow(sgr)>=10) sgr$gene else de$gene, 10)

    # ================================================================
    # fgsea — all collections
    # ================================================================
    cat("    fgsea...\n")
    fres <- NULL
    if (have_fgsea && length(gene_set_list)>0) {
      rnk <- setNames(de$log2FC, de$gene)
      rnk <- sort(rnk[!is.na(rnk)&is.finite(rnk)], decreasing=TRUE)
      gsf <- lapply(gene_set_list, function(gs)
        intersect(gs,names(rnk)))
      gsf <- gsf[sapply(gsf,length)>=5]
      cat(paste0("      Sets >=5: ",length(gsf),"\n"))
      if (length(gsf)>0) {
        fres <- tryCatch(
          fgsea(pathways=gsf, stats=rnk,
                minSize=5, maxSize=500, nPermSimple=10000),
          error=function(e) {
            cat(paste0("      fgsea error: ",
                       conditionMessage(e),"\n")); NULL })
        if (!is.null(fres) && nrow(fres)>0) {
          fres <- fres[order(fres$pval),]
          fres$source <- gs_source[fres$pathway]
          fres$source[is.na(fres$source)] <- "OTHER"
          fsv <- as.data.frame(fres)
          fsv$leadingEdge <- vapply(fsv$leadingEdge,
            function(x) paste(head(x,50),collapse=";"),
            character(1))
          tryCatch(write.csv(
            fsv[!is.na(fsv$padj)&fsv$padj<0.25,],
            file.path(tier2_dir,
                      paste0(cond,"_fgsea_all_padj025.csv")),
            row.names=FALSE), error=function(e) NULL)
          for (src in unique(fsv$source)) {
            sr <- fsv[fsv$source==src &
                        !is.na(fsv$padj) & fsv$padj<0.25,]
            if(nrow(sr)>0) tryCatch(write.csv(sr,
              file.path(tier2_dir,
                        paste0(cond,"_fgsea_",src,".csv")),
              row.names=FALSE), error=function(e) NULL)
          }
          cat("      Per source:\n")
          for (src in c("CUSTOM","HALLMARK","C2_CP","C5_GO",
                        "C6_ONCO","C7_IMMUNO")) {
            nt <- sum(fres$source==src,na.rm=TRUE)
            ns2 <- sum(fres$source==src & !is.na(fres$padj) &
                         fres$padj<0.05, na.rm=TRUE)
            if(nt>0) cat(paste0("        ",src,": ",nt,
                                " tested, ",ns2," sig\n"))
          }
          rm(fsv)
        }
      }
      rm(rnk,gsf); invisible(gc())
    }

    # ================================================================
    # Custom enrichment (Fisher + NES)
    # ================================================================
    cat("    Custom enrichment...\n")
    ce_rows <- lapply(names(custom_sigs), function(sn) {
      sg <- custom_sigs[[sn]]
      it <- intersect(sg,de$gene)
      di <- sgr$gene[sgr$gene %in% sg]
      up <- di[sgr$log2FC[match(di,sgr$gene)]>0]
      dn <- di[sgr$log2FC[match(di,sgr$gene)]<0]
      a <- length(di); b <- nrow(sgr)-a
      cc2 <- length(it)-a; dd <- nrow(de)-nrow(sgr)-cc2
      fp <- tryCatch(fisher.test(
        matrix(c(max(a,0),max(b,0),max(cc2,0),max(dd,0)),nrow=2),
        alternative="greater")$p.value, error=function(e)NA)
      mfc <- mean(de$log2FC[de$gene %in% it], na.rm=TRUE)
      nes <- NA_real_; fpa <- NA_real_
      if (!is.null(fres)) {
        idx <- which(fres$pathway==sn)
        if (length(idx)>0) {
          nes <- fres$NES[idx[1]]; fpa <- fres$padj[idx[1]]
        }
      }
      data.frame(condition=cond, signature=sn,
                 n_in=length(it), n_sig=a,
                 n_up=length(up), n_down=length(dn),
                 fisher_p=fp, mean_FC=mfc,
                 NES=nes, fgsea_padj=fpa,
                 stringsAsFactors=FALSE)
    })
    cedf <- do.call(rbind, ce_rows)
    tryCatch(write.csv(cedf,
      file.path(tier2_dir,paste0(cond,"_custom_enrichment.csv")),
      row.names=FALSE), error=function(e) NULL)
    for (i in seq_len(nrow(cedf))) {
      r <- cedf[i,]
      cat(paste0("      ",sub("^CUSTOM_","",r$signature),": ",
                 r$n_sig,"/",r$n_in,
                 ifelse(!is.na(r$NES),
                        paste0(" NES=",round(r$NES,2)),""),"\n"))
    }

    # Keyword pathway filter
    cgp_pro <- NULL
    if (!is.null(fres)) {
      tryCatch({
        kw <- c("PROSTATE","ANDROGEN","NEUROENDOCRINE","NEPC",
                "BELTRAN","CRPC","CASTRATION","AR_","RB1","TP53",
                "MYCN","EZH2","SMALL_CELL","NEURAL","NEURO",
                "EPITHELIAL_MESENCHYMAL","STEM_CELL","PLASTICITY",
                "ADIPOGEN","MDK","TWEAK","TNFSF12",
                "NECROP","FERROP","PYROP","CUPROP",
                "APOPTOS","AUTOPHAGY","MITOCHOND","OXPHOS",
                "SENESCEN","LIPID","FATTY_ACID","CHOLESTEROL")
        af <- as.data.frame(fres)
        hits <- af[grepl(paste(kw,collapse="|"),af$pathway,
                         ignore.case=TRUE) &
                     !grepl("^CUSTOM_",af$pathway),]
        if (nrow(hits)>0) {
          hits$leadingEdge <- vapply(hits$leadingEdge,
            function(x) if(is.character(x))
              paste(head(x,30),collapse=";") else "",
            character(1))
          hits <- hits[order(hits$pval),]
          cgp_pro <- hits
          tryCatch(write.csv(hits,
            file.path(tier2_dir,
                      paste0(cond,"_keyword_pathways.csv")),
            row.names=FALSE), error=function(e) NULL)
          cat(paste0("      Keyword pathways: ",nrow(hits),"\n"))
        }
      }, error=function(e) NULL)
    }

    # ================================================================
    # PLOTS
    # ================================================================

    # Volcano
    de$label <- ifelse(de$gene %in% t15, de$gene, "")
    p_vol <- ggplot(de, aes(x=log2FC,y=-log10(pvalue+1e-300),
                            color=sig)) +
      geom_point(size=0.8,alpha=0.6) +
      scale_color_manual(values=c("FALSE"="grey60","TRUE"="red")) +
      geom_vline(xintercept=c(-0.5,0.5),linetype="dashed") +
      geom_hline(yintercept=-log10(0.05),linetype="dashed") +
      labs(title=paste0(cond," — Interface vs Distal"),
           subtitle=paste0(nsig," sig"),
           x="log2FC",y="-log10(p)") +
      theme_minimal(base_size=10)
    if (requireNamespace("ggrepel",quietly=TRUE) &&
        any(nchar(de$label)>0))
      p_vol <- p_vol + ggrepel::geom_text_repel(
        aes(label=label),size=2.5,max.overlaps=20,color="black")

    # Heatmap
    tryCatch({
      hc <- list()
      for (bn in abl) {
        v <- bm[[bn]]
        if (is.null(v))
          hc[[bn]] <- setNames(rep(NA_real_,length(t50)),t50)
        else {
          vv <- setNames(rep(NA_real_,length(t50)),t50)
          m <- t50[t50 %in% names(v)]; vv[m] <- v[m]
          hc[[bn]] <- vv
        }
      }
      hm <- do.call(cbind,hc); colnames(hm) <- names(hc)
      rownames(hm) <- t50
      gr2 <- apply(hm,1,function(r)!all(is.na(r)))
      gc2 <- apply(hm,2,function(c2)!all(is.na(c2)))
      hm <- hm[gr2,gc2,drop=FALSE]
      if (nrow(hm)>=2 && ncol(hm)>=2) {
        hm[is.na(hm)] <- 0
        pdf(file.path(tier2_dir,paste0(cond,"_heatmap_top50.pdf")),
            width=10,height=max(8,nrow(hm)*0.2))
        if (have_pheatmap)
          pheatmap(hm,cluster_cols=FALSE,scale="row",
                   main=paste0(cond," — Top 50 DE"),
                   fontsize_row=6,fontsize_col=8)
        else heatmap(hm,Colv=NA,scale="row")
        dev.off()
      }
    }, error=function(e) tryCatch(dev.off(),error=function(e2)NULL))

    # Gradient
    p_grad <- NULL
    tryCatch({
      gl <- list()
      for (bn in abl) {
        v <- bm[[bn]]; if(is.null(v)) next
        for (gg in t10) if(gg %in% names(v))
          gl[[length(gl)+1]] <- data.frame(
            gene=gg,bin=bn,expr=v[gg],stringsAsFactors=FALSE)
      }
      if (length(gl)>0) {
        gd <- do.call(rbind,gl)
        gd$bin <- factor(gd$bin,levels=abl)
        p_grad <- ggplot(gd,aes(x=bin,y=expr,
                                color=gene,group=gene))+
          geom_line(linewidth=0.8)+geom_point(size=2)+
          labs(title=paste0(cond," — Top 10 Gradient"),
               x="Bin",y="Expr")+
          theme_minimal(base_size=10)+
          theme(axis.text.x=element_text(angle=45,hjust=1))
      }
    }, error=function(e) NULL)

    # Combined score bar
    p_cs <- NULL
    tryCatch({
      cd <- head(de[de$cscore>0,],40)
      if(nrow(cd)>0) {
        cd$gene <- factor(cd$gene,levels=rev(cd$gene))
        cd$dir <- ifelse(cd$log2FC>0,"Up","Down")
        p_cs <- ggplot(cd,aes(x=gene,y=cscore,fill=dir))+
          geom_bar(stat="identity")+coord_flip()+
          scale_fill_manual(values=c(Up="firebrick",
                                     Down="steelblue"))+
          labs(title=paste0(cond," — Top 40 Score"),
               x="",y="Score")+
          theme_minimal(base_size=10)+
          theme(axis.text.y=element_text(size=7))
      }
    }, error=function(e) NULL)

    # --- NES bars per custom group ---
    grp_nes_plots <- list()
    tryCatch({
      if (!is.null(fres)) {
        for (grp in names(custom_groups)) {
          gn <- custom_groups[[grp]]
          cf <- as.data.frame(fres[fres$pathway %in% gn,])
          if (nrow(cf)==0) next
          cf$pw <- sub("^CUSTOM_","",cf$pathway)
          cf$sl <- ifelse(!is.na(cf$padj)&cf$padj<0.05,"*","")
          cf$d <- ifelse(cf$NES>0,"Interface","Distal")
          grp_nes_plots[[grp]] <- ggplot(
            cf, aes(x=reorder(pw,NES),y=NES,fill=d))+
            geom_bar(stat="identity")+coord_flip()+
            scale_fill_manual(values=c(Interface="firebrick",
                                       Distal="steelblue"))+
            geom_hline(yintercept=0,linetype="dashed")+
            geom_text(aes(label=sl),hjust=-0.3,size=5)+
            labs(title=paste0(cond," — ",grp," NES"),
                 subtitle="*=padj<0.05",
                 x="",y="NES",fill="")+
            theme_minimal(base_size=11)
        }
      }
    }, error=function(e) NULL)

    # Enrichment count bar (all 21 sigs)
    p_ebar <- NULL
    tryCatch({
      el <- data.frame(
        sig=rep(sub("^CUSTOM_","",cedf$signature),2),
        dir=rep(c("Up","Down"),each=nrow(cedf)),
        n=c(cedf$n_up,cedf$n_down),stringsAsFactors=FALSE)
      p_ebar <- ggplot(el,aes(x=sig,y=n,fill=dir))+
        geom_bar(stat="identity",position="dodge")+
        scale_fill_manual(values=c(Up="firebrick",
                                   Down="steelblue"))+
        labs(title=paste0(cond," — DE per Signature"),
             x="",y="Count")+
        theme_minimal(base_size=8)+
        theme(axis.text.x=element_text(angle=70,hjust=1))
    }, error=function(e) NULL)

    # --- Per-source MSigDB pathway bars ---
    src_plots <- list()
    for (src in c("HALLMARK","C2_CP","C5_GO",
                  "C6_ONCO","C7_IMMUNO")) {
      tryCatch({
        if (!is.null(fres)) {
          sf <- as.data.frame(fres[fres$source==src,])
          ms <- sf[!is.na(sf$padj)&sf$padj<0.05,]
          if (nrow(ms)>0) {
            tp <- head(ms[ms$NES>0,][
              order(-ms$NES[ms$NES>0]),],15)
            tn <- head(ms[ms$NES<0,][
              order(ms$NES[ms$NES<0]),],15)
            tb <- rbind(tp,tn)
            if (nrow(tb)>0) {
              tb$pw <- substring(gsub(
                paste0("^HALLMARK_|^GOBP_|^GOCC_|^GOMF_|",
                       "^REACTOME_|^KEGG_|^WP_|^PID_|^BIOCARTA_"),
                "",tb$pathway),1,50)
              tb$d <- ifelse(tb$NES>0,"Interface","Distal")
              src_plots[[src]] <- ggplot(
                tb,aes(x=reorder(pw,NES),y=NES,fill=d))+
                geom_bar(stat="identity")+coord_flip()+
                scale_fill_manual(
                  values=c(Interface="firebrick",
                           Distal="steelblue"))+
                geom_hline(yintercept=0)+
                labs(title=paste0(cond," — ",src," (padj<0.05)"),
                     subtitle=paste0(nrow(ms)," sig"),
                     x="",y="NES")+
                theme_minimal(base_size=9)+
                theme(axis.text.y=element_text(size=6))
            }
          }
        }
      }, error=function(e) NULL)
    }

    # Keyword pathway bar
    p_cgp <- NULL
    tryCatch({
      if (!is.null(cgp_pro) && nrow(cgp_pro)>0) {
        ct <- head(cgp_pro,30)
        ct$pw <- substring(ct$pathway,1,50)
        ct$d <- ifelse(ct$NES>0,"Interface","Distal")
        p_cgp <- ggplot(ct,aes(x=reorder(pw,NES),y=NES,fill=d))+
          geom_bar(stat="identity")+coord_flip()+
          scale_fill_manual(values=c(Interface="firebrick",
                                     Distal="steelblue"))+
          labs(title=paste0(cond," — Keyword Pathways"),
               x="",y="NES")+
          theme_minimal(base_size=9)+
          theme(axis.text.y=element_text(size=6))
      }
    }, error=function(e) NULL)

    # fgsea enrichment curves (all 21 custom sigs)
    fcurves <- list()
    if (have_fgsea && !is.null(fres)) {
      rk <- setNames(de$log2FC, de$gene)
      rk <- sort(rk[!is.na(rk)&is.finite(rk)], decreasing=TRUE)
      for (sn2 in names(custom_sigs)) {
        tryCatch({
          gg <- intersect(custom_sigs[[sn2]], names(rk))
          if (length(gg)>=5)
            fcurves[[sn2]] <- plotEnrichment(gg,rk) +
              labs(title=paste0(cond," — ",
                                sub("^CUSTOM_","",sn2)))
        }, error=function(e) NULL)
      }
      rm(rk)
    }

    # ================================================================
    # Summary PDF
    # ================================================================
    tryCatch({
      pdf(file.path(tier2_dir,
                    paste0(cond,"_TIER2_summary.pdf")),
          width=11,height=8.5)

      # Distance histogram
      tryCatch({
        dv <- list()
        for (sn in cond_samples) {
          obj <- load_sample(sn)
          if (is.null(obj) ||
              !"dist_to_interface" %in% colnames(obj@meta.data)) {
            if (!is.null(obj)) rm(obj); invisible(gc()); next
          }
          dv[[sn]] <- data.frame(d=obj$dist_to_interface,
                                 stringsAsFactors=FALSE)
          rm(obj); invisible(gc())
        }
        if (length(dv)>0) {
          dd <- do.call(rbind,dv)
          print(ggplot(dd,aes(x=d))+
                  geom_histogram(bins=50,fill="steelblue",
                                 color="white")+
                  geom_vline(xintercept=bi$breaks,
                             linetype="dashed",color="red")+
                  labs(title=paste0(cond," — Distances"),
                       x="Distance",y="Count")+
                  theme_minimal(base_size=10))
          rm(dd,dv)
        }
      }, error=function(e) NULL)

      # Core plots
      print(p_vol)
      if (!is.null(p_cs)) print(p_cs)
      if (!is.null(p_grad)) print(p_grad)

      # MA plot
      tryCatch({
        de$me <- (de$mean_iface+de$mean_distal)/2
        print(ggplot(de,aes(x=log10(me+0.01),y=log2FC,color=sig))+
                geom_point(size=0.6,alpha=0.5)+
                scale_color_manual(values=c("FALSE"="grey60",
                                            "TRUE"="red"))+
                geom_hline(yintercept=c(-0.5,0.5),linetype="dashed")+
                labs(title=paste0(cond," — MA"),
                     x="log10(Mean)",y="log2FC")+
                theme_minimal(base_size=10))
      }, error=function(e) NULL)

      # NES by group
      for (grp in names(grp_nes_plots))
        print(grp_nes_plots[[grp]])

      if (!is.null(p_ebar)) print(p_ebar)

      # MSigDB source bars
      for (src in names(src_plots))
        print(src_plots[[src]])

      if (!is.null(p_cgp)) print(p_cgp)

      # fgsea curves (grouped)
      for (grp in names(custom_groups)) {
        for (sn2 in custom_groups[[grp]]) {
          if (!is.null(fcurves[[sn2]]))
            print(fcurves[[sn2]])
        }
      }

      # Signature volcanos (all 21)
      for (sn3 in names(custom_sigs)) {
        tryCatch({
          de$in_s <- de$gene %in% custom_sigs[[sn3]]
          if (sum(de$in_s)<3) next
          de$sl <- ifelse(de$in_s & de$sig, de$gene, "")
          pp <- ggplot(de,aes(x=log2FC,y=-log10(pvalue+1e-300)))+
            geom_point(data=de[!de$in_s,],
                       size=0.5,alpha=0.3,color="grey80")+
            geom_point(data=de[de$in_s,],
                       aes(color=sig),size=2,alpha=0.8)+
            scale_color_manual(values=c("FALSE"="orange",
                                        "TRUE"="red"))+
            geom_vline(xintercept=c(-0.5,0.5),linetype="dashed")+
            labs(title=paste0(cond," — ",
                              sub("^CUSTOM_","",sn3)),
                 x="log2FC",y="-log10(p)")+
            theme_minimal(base_size=10)
          if (requireNamespace("ggrepel",quietly=TRUE) &&
              any(nchar(de$sl)>0))
            pp <- pp + ggrepel::geom_text_repel(
              data=de[de$sl!="",],
              aes(label=sl),size=2.5,color="black",
              max.overlaps=20)
          print(pp)
        }, error=function(e) NULL)
      }

      dev.off()
      cat(paste0("    Summary PDF saved\n"))
    }, error=function(e) {
      cat(paste0("    Summary PDF error: ",
                 conditionMessage(e),"\n"))
      tryCatch(dev.off(),error=function(e2) NULL)
    })

    # ================================================================
    # Enrichment PDF
    # ================================================================
    tryCatch({
      pdf(file.path(tier2_dir,
                    paste0(cond,"_TIER2_enrichment.pdf")),
          width=12,height=9)

      # NES by group
      for (grp in names(grp_nes_plots))
        print(grp_nes_plots[[grp]])

      if (!is.null(p_ebar)) print(p_ebar)

      # fgsea curves grouped
      for (grp in names(custom_groups)) {
        for (sn2 in custom_groups[[grp]]) {
          if (!is.null(fcurves[[sn2]]))
            print(fcurves[[sn2]])
        }
      }

      # MSigDB source bars
      for (src in names(src_plots))
        print(src_plots[[src]])

      if (!is.null(p_cgp)) print(p_cgp)

      # Signature volcanos
      for (sn3 in names(custom_sigs)) {
        tryCatch({
          de$in_s <- de$gene %in% custom_sigs[[sn3]]
          if (sum(de$in_s)<3) next
          de$sl <- ifelse(de$in_s & de$sig, de$gene, "")
          pp <- ggplot(de,aes(x=log2FC,y=-log10(pvalue+1e-300)))+
            geom_point(data=de[!de$in_s,],
                       size=0.5,alpha=0.3,color="grey80")+
            geom_point(data=de[de$in_s,],
                       aes(color=sig),size=2,alpha=0.8)+
            scale_color_manual(values=c("FALSE"="orange",
                                        "TRUE"="red"))+
            geom_vline(xintercept=c(-0.5,0.5),linetype="dashed")+
            labs(title=paste0(cond," — ",
                              sub("^CUSTOM_","",sn3)),
                 x="log2FC",y="-log10(p)")+
            theme_minimal(base_size=10)
          if (requireNamespace("ggrepel",quietly=TRUE) &&
              any(nchar(de$sl)>0))
            pp <- pp + ggrepel::geom_text_repel(
              data=de[de$sl!="",],
              aes(label=sl),size=2.5,color="black",
              max.overlaps=20)
          print(pp)
        }, error=function(e) NULL)
      }

      dev.off()
      cat(paste0("    Enrichment PDF saved\n"))
    }, error=function(e) {
      cat(paste0("    Enrichment PDF error: ",
                 conditionMessage(e),"\n"))
      tryCatch(dev.off(),error=function(e2) NULL)
    })

    # ================================================================
    # Clean up condition-level objects
    # ================================================================
    rm(bel, bm, sgl, de, sgr, cedf)
    if (exists("fres",inherits=FALSE)) rm(fres)
    if (exists("fcurves",inherits=FALSE)) rm(fcurves)
    if (exists("cgp_pro",inherits=FALSE)) rm(cgp_pro)
    if (exists("p_vol",inherits=FALSE)) rm(p_vol)
    if (exists("p_grad",inherits=FALSE)) rm(p_grad)
    if (exists("p_cs",inherits=FALSE)) rm(p_cs)
    if (exists("p_ebar",inherits=FALSE)) rm(p_ebar)
    if (exists("p_cgp",inherits=FALSE)) rm(p_cgp)
    if (exists("grp_nes_plots",inherits=FALSE)) rm(grp_nes_plots)
    if (exists("src_plots",inherits=FALSE)) rm(src_plots)
    invisible(gc())

  }  # end for (cond in conditions)

  tier2_status <- "COMPLETE"
  cat("\nTIER 2 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 2 ERROR: ", conditionMessage(e), "\n"))
  tier2_status <<- paste0("ERROR: ", conditionMessage(e))
  tryCatch(dev.off(), error = function(e2) NULL)
})
# ==============================================================================
# SECTION 5: TIER 3 — Gradient Modeling (Distance-Expression Relationships)
# ==============================================================================
cat("\n=== TIER 3: Gradient Modeling ===\n")
tier3_status <- "SKIPPED"

tryCatch({
  tier3_dir <- file.path(visium_root, "TIER3_Gradient")
  all_sig_genes <- unique(c(ASPC_genes, Hallmark_AR_genes,
                            NEPC_Beltran_UP, NEPC_Beltran_DOWN))
  gradient_results <- list()

  conditions <- unique(sapply(samples_with_dist, extract_condition))

  for (cond in conditions) {
    cat(paste0("\n  TIER3 — condition: ", cond, "\n"))
    cond_samples <- samples_with_dist[
      sapply(samples_with_dist, extract_condition) == cond]
    if (length(cond_samples) == 0) next

    # ==============================================================
    # Phase 1: Collect expression + distance for signature genes
    #   Load one sample at a time; store as data.frames
    # ==============================================================
    spot_data_list <- list()
    sample_gene_counts <- list()

    for (sn in cond_samples) {
      obj <- load_sample(sn)
      if (is.null(obj)) next
      if (!"dist_to_interface" %in% colnames(obj@meta.data)) {
        rm(obj); invisible(gc()); next
      }

      expr_mat <- safe_GetAssayData(obj, layer = "data")
      if (is.null(expr_mat)) { rm(obj); invisible(gc()); next }

      genes_avail <- intersect(all_sig_genes, rownames(expr_mat))
      if (length(genes_avail) == 0) { rm(obj, expr_mat); invisible(gc()); next }

      sample_gene_counts[[sn]] <- genes_avail

      # Extract expression for signature genes only (memory efficient)
      expr_sub <- expr_mat[genes_avail, , drop = FALSE]
      df <- as.data.frame(t(as.matrix(expr_sub)))
      df$dist <- obj$dist_to_interface
      df$sample_id <- sn

      # Add score columns if available
      md <- obj@meta.data
      for (sc in c("ASPC_score","AR_score","NEPC_UP_score","NEPC_Beltran_NET")) {
        df[[sc]] <- if (sc %in% colnames(md)) md[[sc]] else NA_real_
      }

      spot_data_list[[sn]] <- df

      rm(obj, expr_mat, expr_sub, df); invisible(gc())
    }

    if (length(spot_data_list) == 0) {
      cat("    No spot data — skipping\n"); next
    }

    # Find common signature genes across all samples in this condition
    if (length(sample_gene_counts) >= 2) {
      common_sig_genes <- Reduce(intersect, sample_gene_counts)
    } else {
      common_sig_genes <- sample_gene_counts[[1]]
    }
    cat(paste0("    Common signature genes: ", length(common_sig_genes), "\n"))

    if (length(common_sig_genes) < 5) {
      cat("    Too few common genes — skipping\n")
      rm(spot_data_list, sample_gene_counts); invisible(gc())
      next
    }

    # Keep only common genes + metadata columns in each data.frame
    meta_cols <- c("dist","sample_id","ASPC_score","AR_score",
                   "NEPC_UP_score","NEPC_Beltran_NET")
    keep_cols <- c(common_sig_genes, meta_cols)

    spot_data_list <- lapply(spot_data_list, function(df) {
      cols_here <- intersect(keep_cols, colnames(df))
      df[, cols_here, drop = FALSE]
    })

    spot_data <- tryCatch(
      do.call(rbind, spot_data_list),
      error = function(e) {
        cat(paste0("    rbind error: ", conditionMessage(e), "\n"))
        NULL
      }
    )
    rm(spot_data_list, sample_gene_counts); invisible(gc())

    if (is.null(spot_data) || nrow(spot_data) < 20) {
      cat("    Too few spots — skipping\n"); next
    }

    cat(paste0("    Combined: ", nrow(spot_data), " spots, ",
               length(common_sig_genes), " genes\n"))

    # Report distance range
    d_range <- range(spot_data$dist, na.rm = TRUE)
    cat(paste0("    Distance range: [", round(d_range[1], 2), ", ",
               round(d_range[2], 2), "]\n"))

    # ==============================================================
    # Phase 2: Fit linear + exponential models for each gene
    # ==============================================================
    genes_to_model <- common_sig_genes[common_sig_genes %in% colnames(spot_data)]
    gene_results <- list()

    for (g in genes_to_model) {
      tryCatch({
        y <- log1p(spot_data[[g]])
        x <- spot_data$dist
        ok <- !is.na(y) & !is.na(x) & is.finite(y) & is.finite(x)
        if (sum(ok) < 10) next

        # --- Linear model ---
        fit_lm <- lm(y[ok] ~ x[ok])
        sm_lm  <- summary(fit_lm)
        slope  <- coef(fit_lm)[2]
        r2_lm  <- sm_lm$r.squared
        p_lm   <- tryCatch(coef(sm_lm)[2, 4], error = function(e) NA)

        # --- Exponential decay NLS ---
        decay_rate <- NA; half_dist <- NA; r2_nls <- NA
        tryCatch({
          y_range <- max(y[ok]) - min(y[ok])
          if (y_range > 0.01) {
            fit_nls <- nls(y[ok] ~ a * exp(-lambda * x[ok]) + b,
                           start = list(a = y_range, lambda = 0.1,
                                        b = min(y[ok])),
                           control = nls.control(maxiter = 200,
                                                 warnOnly = TRUE))
            lam <- coef(fit_nls)["lambda"]
            pred_nls <- predict(fit_nls)
            ss_res <- sum((y[ok] - pred_nls)^2)
            ss_tot <- sum((y[ok] - mean(y[ok]))^2)
            r2_nls <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA
            decay_rate <- lam
            half_dist  <- if (!is.na(lam) && lam > 0) log(2) / lam else NA
          }
        }, error = function(e) NULL)

        # --- Classify gene gradient ---
        direction <- if (!is.na(slope) && slope > 0) "increasing" else "decreasing"

        # Signature membership
        in_sig <- c()
        if (g %in% ASPC_genes)         in_sig <- c(in_sig, "ASPC")
        if (g %in% Hallmark_AR_genes)  in_sig <- c(in_sig, "AR")
        if (g %in% NEPC_Beltran_UP)    in_sig <- c(in_sig, "NEPC_UP")
        if (g %in% NEPC_Beltran_DOWN)  in_sig <- c(in_sig, "NEPC_DOWN")

        gene_results[[g]] <- data.frame(
          gene       = g,
          condition  = cond,
          signature  = paste(in_sig, collapse = ";"),
          direction  = direction,
          slope_lm   = slope,
          p_lm       = p_lm,
          R2_lm      = r2_lm,
          decay_rate = decay_rate,
          half_dist  = half_dist,
          R2_nls     = r2_nls,
          n_spots    = sum(ok),
          stringsAsFactors = FALSE
        )
      }, error = function(e) NULL)
    }

    if (length(gene_results) == 0) {
      cat("    No genes modeled\n")
      rm(spot_data); invisible(gc())
      next
    }

    cond_df <- do.call(rbind, gene_results)
    cond_df <- cond_df[order(cond_df$R2_lm, decreasing = TRUE), ]
    gradient_results[[cond]] <- cond_df

    write.csv(cond_df,
              file.path(tier3_dir, paste0(cond, "_gradient_models.csv")),
              row.names = FALSE)

    n_sig_grad <- sum(cond_df$p_lm < 0.05 & cond_df$R2_lm > 0.01, na.rm = TRUE)
    cat(paste0("    Modeled ", nrow(cond_df), " genes, ",
               n_sig_grad, " with significant gradient (p<0.05, R2>0.01)\n"))

    # ==============================================================
    # Phase 3: Multi-page PDF output
    # ==============================================================
    cat(paste0("    Generating PDF for ", cond, "...\n"))

    pdf(file.path(tier3_dir, paste0(cond, "_TIER3_gradient_results.pdf")),
        width = 11, height = 8.5)

    # ---- Page 1: Summary bar plot — R² distribution ----
    cond_df$sig_gradient <- !is.na(cond_df$p_lm) &
      cond_df$p_lm < 0.05 & cond_df$R2_lm > 0.01
    p_r2 <- ggplot(cond_df, aes(x = R2_lm, fill = sig_gradient)) +
      geom_histogram(bins = 40, alpha = 0.8) +
      scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "steelblue"),
                        labels = c("NS", "Sig (p<0.05, R²>0.01)")) +
      labs(title = paste0(cond, " — Linear Model R² Distribution"),
           subtitle = paste0(nrow(cond_df), " genes modeled, ",
                             n_sig_grad, " significant"),
           x = "R² (linear)", y = "Count", fill = "") +
      theme_minimal(base_size = 11)
    print(p_r2)

    # ---- Page 2: Top 30 genes ranked by |slope| ----
    top30 <- head(cond_df[order(abs(cond_df$slope_lm), decreasing = TRUE), ], 30)
    top30$gene <- factor(top30$gene, levels = rev(top30$gene))
    p_bar <- ggplot(top30, aes(x = gene, y = slope_lm, fill = signature)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = paste0(cond, " — Top 30 Genes by Gradient Slope"),
           x = "", y = "Slope (linear model)", fill = "Signature") +
      theme_minimal(base_size = 10) +
      theme(axis.text.y = element_text(size = 7))
    print(p_bar)

    # ---- Page 3: Slope vs R² scatter (all genes) ----
    p_sr <- ggplot(cond_df, aes(x = slope_lm, y = R2_lm,
                                color = signature)) +
      geom_point(size = 1.5, alpha = 0.6) +
      geom_hline(yintercept = 0.01, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      labs(title = paste0(cond, " — Slope vs R² for All Signature Genes"),
           x = "Slope", y = "R² (linear)", color = "Signature") +
      theme_minimal(base_size = 10)

    # Label top 10 by R²
    top10_r2 <- head(cond_df[order(cond_df$R2_lm, decreasing = TRUE), ], 10)
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p_sr <- p_sr +
        ggrepel::geom_text_repel(
          data = top10_r2,
          aes(label = gene), size = 2.5, color = "black",
          max.overlaps = 15
        )
    }
    print(p_sr)

    # ---- Page 4: Per-signature boxplot of slopes ----
    if (any(nchar(cond_df$signature) > 0)) {
      # Expand multi-signature genes to separate rows
      sig_expanded <- do.call(rbind, lapply(seq_len(nrow(cond_df)), function(i) {
        sigs <- strsplit(cond_df$signature[i], ";")[[1]]
        if (length(sigs) == 0 || all(sigs == "")) sigs <- "none"
        data.frame(gene = cond_df$gene[i],
                   sig = sigs,
                   slope = cond_df$slope_lm[i],
                   R2 = cond_df$R2_lm[i],
                   stringsAsFactors = FALSE)
      }))
      p_box <- ggplot(sig_expanded, aes(x = sig, y = slope, fill = sig)) +
        geom_boxplot(outlier.size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(title = paste0(cond, " — Gradient Slopes by Signature"),
             x = "Signature", y = "Slope") +
        theme_minimal(base_size = 10) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 30, hjust = 1))
      print(p_box)
    }

    # ---- Page 5: Score ~ distance scatter plots ----
    score_cols_t3 <- c("ASPC_score","AR_score","NEPC_UP_score","NEPC_Beltran_NET")
    score_plots <- list()
    for (sc in score_cols_t3) {
      if (!sc %in% colnames(spot_data) || all(is.na(spot_data[[sc]]))) next
      tmp <- data.frame(dist = spot_data$dist, score = spot_data[[sc]])
      tmp <- tmp[!is.na(tmp$dist) & !is.na(tmp$score), ]
      if (nrow(tmp) < 10) next

      rho <- tryCatch({
        ct <- cor.test(tmp$dist, tmp$score, method = "spearman",
                       use = "complete.obs", exact = FALSE)
        paste0("rho=", round(ct$estimate, 3),
               ", p=", formatC(ct$p.value, format = "e", digits = 2))
      }, error = function(e) "")

      score_plots[[sc]] <- ggplot(tmp, aes(x = dist, y = score)) +
        geom_point(size = 0.3, alpha = 0.15, color = "grey40") +
        geom_smooth(method = "loess", color = "blue",
                    linewidth = 1, se = TRUE) +
        labs(title = paste0(sc, "  ", rho),
             x = "Distance to Interface", y = sc) +
        theme_minimal(base_size = 9)
    }
    if (length(score_plots) > 0) {
      print(wrap_plots(score_plots, ncol = 2) +
              plot_annotation(
                title = paste0(cond, " — Signature Scores vs Distance")))
    }

    # ---- Pages 6+: Top 20 gene decay curves ----
    top20_genes <- head(
      cond_df$gene[order(abs(cond_df$slope_lm), decreasing = TRUE)], 20)
    top20_genes <- top20_genes[top20_genes %in% colnames(spot_data)]

    if (length(top20_genes) > 0) {
      gene_curve_plots <- list()
      for (gg in top20_genes) {
        tmp_df <- data.frame(
          dist = spot_data$dist,
          expr = log1p(spot_data[[gg]])
        )
        tmp_df <- tmp_df[!is.na(tmp_df$dist) & !is.na(tmp_df$expr), ]
        if (nrow(tmp_df) < 5) next

        # Get model info for subtitle
        gi <- cond_df[cond_df$gene == gg, ]
        sub_txt <- paste0("slope=", round(gi$slope_lm, 4),
                          "  R²=", round(gi$R2_lm, 3),
                          "  p=", formatC(gi$p_lm, format = "e", digits = 2))

        p_curve <- ggplot(tmp_df, aes(x = dist, y = expr)) +
          geom_point(size = 0.3, alpha = 0.15, color = "grey50") +
          geom_smooth(method = "lm", color = "red", linewidth = 0.8,
                      se = TRUE) +
          geom_smooth(method = "loess", color = "blue", linewidth = 0.8,
                      se = FALSE, linetype = "dashed") +
          labs(title = gg, subtitle = sub_txt,
               x = "Distance", y = "log1p(expr)") +
          theme_minimal(base_size = 8)

        gene_curve_plots[[gg]] <- p_curve
      }

      if (length(gene_curve_plots) > 0) {
        n_per_page <- 6
        for (i in seq(1, length(gene_curve_plots), by = n_per_page)) {
          sub <- gene_curve_plots[i:min(i + n_per_page - 1,
                                        length(gene_curve_plots))]
          print(wrap_plots(sub, ncol = 3) +
                  plot_annotation(
                    title = paste0(cond,
                                   " — Gene Expression vs Distance to Interface")))
        }
      }
    }

    # ---- Decay rate summary page (genes with successful NLS) ----
    nls_ok <- cond_df[!is.na(cond_df$decay_rate) & !is.na(cond_df$R2_nls), ]
    if (nrow(nls_ok) > 0) {
      nls_ok <- nls_ok[order(nls_ok$R2_nls, decreasing = TRUE), ]
      top_nls <- head(nls_ok, 30)
      top_nls$gene <- factor(top_nls$gene, levels = rev(top_nls$gene))

      p_decay <- ggplot(top_nls, aes(x = gene, y = decay_rate, fill = R2_nls)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_viridis_c(option = "plasma") +
        labs(title = paste0(cond, " — Top 30 Exponential Decay Rates"),
             subtitle = "Higher decay rate = steeper drop from interface",
             x = "", y = "Decay Rate (lambda)", fill = "NLS R²") +
        theme_minimal(base_size = 10) +
        theme(axis.text.y = element_text(size = 7))
      print(p_decay)

      # Half-distance summary
      nls_hd <- nls_ok[!is.na(nls_ok$half_dist) & nls_ok$half_dist > 0 &
                          is.finite(nls_ok$half_dist), ]
      if (nrow(nls_hd) > 0) {
        top_hd <- head(nls_hd[order(nls_hd$half_dist), ], 30)
        top_hd$gene <- factor(top_hd$gene, levels = rev(top_hd$gene))
        p_hd <- ggplot(top_hd, aes(x = gene, y = half_dist, fill = signature)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          labs(title = paste0(cond, " — Spatial Half-Distance (exp decay)"),
               subtitle = "Distance at which gene expression drops to 50%",
               x = "", y = "Half-Distance (array units)", fill = "Signature") +
          theme_minimal(base_size = 10) +
          theme(axis.text.y = element_text(size = 7))
        print(p_hd)
      }
    }

    # ---- Heatmap: top genes binned by distance quantiles ----
    # Bin spots into 5 equal-count groups by distance
    spot_data$dist_quantile <- tryCatch({
      as.character(cut(spot_data$dist,
                       breaks = quantile(spot_data$dist,
                                         probs = seq(0, 1, 0.2),
                                         na.rm = TRUE),
                       include.lowest = TRUE, labels = FALSE))
    }, error = function(e) NA)

    if (!all(is.na(spot_data$dist_quantile))) {
      top40 <- head(cond_df$gene[order(abs(cond_df$slope_lm),
                                       decreasing = TRUE)], 40)
      top40 <- top40[top40 %in% colnames(spot_data)]

      if (length(top40) >= 3) {
        quant_means <- sapply(sort(unique(spot_data$dist_quantile)), function(q) {
          idx <- spot_data$dist_quantile == q
          if (sum(idx, na.rm = TRUE) < 3) return(rep(NA, length(top40)))
          colMeans(spot_data[idx, top40, drop = FALSE], na.rm = TRUE)
        })

        if (is.matrix(quant_means) && nrow(quant_means) >= 3 &&
            ncol(quant_means) >= 2) {
          rownames(quant_means) <- top40
          colnames(quant_means) <- paste0("Q", sort(unique(spot_data$dist_quantile)))

          quant_means[is.na(quant_means)] <- 0
          quant_scaled <- t(scale(t(quant_means)))
          quant_scaled[quant_scaled > 2]  <-  2
          quant_scaled[quant_scaled < -2] <- -2
          quant_scaled[is.na(quant_scaled)] <- 0

          if (have_pheatmap) {
            pheatmap(quant_scaled, cluster_cols = FALSE,
                     main = paste0(cond,
                                   " — Top 40 Gradient Genes (distance quintiles)"),
                     fontsize_row = 7, fontsize_col = 10)
          } else {
            heatmap(quant_scaled, Colv = NA, scale = "none",
                    main = paste0(cond, " — Gradient Genes"))
          }
        }
      }
    }

    dev.off()
    cat(paste0("    PDF saved: ", cond, "_TIER3_gradient_results.pdf\n"))

    # Also save standalone decay curve plots
    if (length(top20_genes) > 0) {
      pdf(file.path(tier3_dir, paste0(cond, "_top20_decay_curves.pdf")),
          width = 11, height = 8)
      if (exists("gene_curve_plots", inherits = FALSE) &&
          length(gene_curve_plots) > 0) {
        n_per_page <- 6
        for (i in seq(1, length(gene_curve_plots), by = n_per_page)) {
          sub <- gene_curve_plots[i:min(i + n_per_page - 1,
                                        length(gene_curve_plots))]
          print(wrap_plots(sub, ncol = 3))
        }
      }
      dev.off()
    }

    rm(spot_data, gene_results, cond_df)
    if (exists("gene_curve_plots", inherits = FALSE)) rm(gene_curve_plots)
    if (exists("sig_expanded", inherits = FALSE)) rm(sig_expanded)
    invisible(gc())
  }

  # ==============================================================
  # Cross-condition gradient comparison (if multiple conditions)
  # ==============================================================
  if (length(gradient_results) >= 2) {
    cat("\n  Cross-condition gradient comparison...\n")

    all_grad <- do.call(rbind, gradient_results)
    write.csv(all_grad,
              file.path(tier3_dir, "all_conditions_gradient_models.csv"),
              row.names = FALSE)

    pdf(file.path(tier3_dir, "cross_condition_gradient_comparison.pdf"),
        width = 12, height = 9)

    # Page 1: slope distribution by condition
    p_cross_slope <- ggplot(all_grad, aes(x = condition, y = slope_lm,
                                          fill = condition)) +
      geom_boxplot(outlier.size = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(title = "Gradient Slopes Across Conditions (All Signature Genes)",
           x = "Condition", y = "Slope (linear model)") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "none")
    print(p_cross_slope)

    # Page 2: slope by condition AND signature
    sig_exp_all <- do.call(rbind, lapply(seq_len(nrow(all_grad)), function(i) {
      sigs <- strsplit(all_grad$signature[i], ";")[[1]]
      if (length(sigs) == 0 || all(sigs == "")) sigs <- "other"
      data.frame(gene = all_grad$gene[i],
                 condition = all_grad$condition[i],
                 sig = sigs,
                 slope = all_grad$slope_lm[i],
                 stringsAsFactors = FALSE)
    }))
    p_cross_sig <- ggplot(sig_exp_all, aes(x = condition, y = slope,
                                           fill = sig)) +
      geom_boxplot(outlier.size = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      facet_wrap(~ sig, scales = "free_y") +
      labs(title = "Gradient Slopes by Signature Across Conditions",
           x = "Condition", y = "Slope") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    print(p_cross_sig)

    # Page 3: Top 20 genes with steepest gradient IN ANY condition
    top20_overall <- head(
      all_grad[order(abs(all_grad$slope_lm), decreasing = TRUE), ], 60)
    top20_unique <- unique(top20_overall$gene)[1:min(20,
                                                     length(unique(top20_overall$gene)))]
    cross_sub <- all_grad[all_grad$gene %in% top20_unique, ]
    cross_sub$gene <- factor(cross_sub$gene, levels = rev(top20_unique))
    p_cross_top <- ggplot(cross_sub, aes(x = gene, y = slope_lm,
                                         fill = condition)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      labs(title = "Top 20 Gradient Genes — Slope Across Conditions",
           x = "", y = "Slope", fill = "Condition") +
      theme_minimal(base_size = 10)
    print(p_cross_top)

    # Page 4: Heatmap of slopes — genes × conditions
    slope_mat_genes <- unique(all_grad$gene)
    slope_mat_conds <- unique(all_grad$condition)
    slope_mat <- matrix(0, nrow = length(slope_mat_genes),
                        ncol = length(slope_mat_conds),
                        dimnames = list(slope_mat_genes, slope_mat_conds))
    for (i in seq_len(nrow(all_grad))) {
      slope_mat[all_grad$gene[i], all_grad$condition[i]] <- all_grad$slope_lm[i]
    }
    # Keep genes with highest variance across conditions
    gene_var_slope <- apply(slope_mat, 1, var, na.rm = TRUE)
    top_var_genes <- names(sort(gene_var_slope, decreasing = TRUE))[
      1:min(50, length(gene_var_slope))]
    slope_sub <- slope_mat[top_var_genes, , drop = FALSE]

    if (have_pheatmap && nrow(slope_sub) >= 3 && ncol(slope_sub) >= 2) {
      pheatmap(slope_sub, scale = "row", cluster_cols = FALSE,
               main = "Top 50 Variable-Gradient Genes Across Conditions",
               fontsize_row = 6)
    }

    dev.off()
    cat("    Cross-condition PDF saved\n")

    rm(all_grad, sig_exp_all); invisible(gc())
  }

  tier3_status <- "COMPLETE"
  cat("\nTIER 3 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 3 ERROR: ", conditionMessage(e), "\n"))
  tier3_status <<- paste0("ERROR: ", conditionMessage(e))
  tryCatch(dev.off(), error = function(e2) NULL)
})
# ==============================================================================
# SECTION 6: TIER 4 — Ligand-Receptor Communication (CellChatDB)
# ==============================================================================
cat("\n=== TIER 4: Ligand-Receptor Communication (CellChatDB) ===\n")
tier4_status <- "SKIPPED"

tryCatch({
  tier4_dir <- file.path(visium_root, "TIER4_LR")

  # --- Install CellChat if needed ---
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    cat("  Installing CellChat...\n")
    tryCatch({
      if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
      devtools::install_github("jinworks/CellChat")
    }, error = function(e) {
      cat(paste0("  CellChat installation failed: ", conditionMessage(e), "\n"))
    })
  }

  have_cellchat <- requireNamespace("CellChat", quietly = TRUE)

  if (!have_cellchat) {
    cat("  CellChat not available — using built-in fallback L-R database\n")
    lr_pairs_db <- data.frame(
      ligand       = c("FGF2","TGFB1","WNT5A","BMP4","PDGFA","VEGFA","HGF","IL6",
                        "CXCL12","TNF","IGF1","ANGPT1","EGF","SHH","DLL1","JAG1"),
      receptor     = c("FGFR1","TGFBR1","FZD5","BMPR1A","PDGFRA","KDR","MET","IL6R",
                        "CXCR4","TNFRSF1A","IGF1R","TEK","EGFR","PTCH1","NOTCH1","NOTCH2"),
      pathway_name = c("FGF","TGFb","WNT","BMP","PDGF","VEGF","HGF","IL6",
                        "CXCL","TNF","IGF","Angiopoietin","EGF","Hedgehog","NOTCH","NOTCH"),
      annotation   = rep("Secreted Signaling", 16),
      stringsAsFactors = FALSE
    )
    cat(paste0("  Fallback L-R database: ", nrow(lr_pairs_db), " pairs\n"))

  } else {
    library(CellChat)
    cat("  CellChat loaded\n")

    CellChatDB <- CellChatDB.human
    interaction_df <- CellChatDB$interaction
    complex_df     <- CellChatDB$complex
    cat(paste0("  CellChatDB.human: ", nrow(interaction_df), " interactions, ",
               nrow(complex_df), " complexes\n"))

    resolve_genes <- function(name_str, complex_lookup) {
      if (is.na(name_str) || name_str == "") return(character(0))
      if (name_str %in% rownames(complex_lookup)) {
        row <- complex_lookup[name_str, ]
        sub_cols <- grep("^subunit", colnames(complex_lookup), value = TRUE)
        subunits <- unlist(row[sub_cols])
        subunits <- subunits[!is.na(subunits) & subunits != ""]
        return(as.character(subunits))
      }
      return(name_str)
    }

    complex_lookup <- complex_df
    if ("complex_name" %in% colnames(complex_lookup)) {
      rownames(complex_lookup) <- complex_lookup$complex_name
    }

    cat("  Resolving interactions to gene-level pairs...\n")
    lr_pairs_list <- list()
    for (i in seq_len(nrow(interaction_df))) {
      lig_name <- interaction_df$ligand[i]
      rec_name <- interaction_df$receptor[i]
      pw_name  <- interaction_df$pathway_name[i]
      annot    <- if ("annotation" %in% colnames(interaction_df)) interaction_df$annotation[i] else "Unknown"

      lig_genes <- resolve_genes(lig_name, complex_lookup)
      rec_genes <- resolve_genes(rec_name, complex_lookup)

      if (length(lig_genes) == 0 || length(rec_genes) == 0) next

      lr_pairs_list[[i]] <- data.frame(
        interaction_name = interaction_df$interaction_name[i],
        ligand_name      = lig_name,
        receptor_name    = rec_name,
        ligand_genes     = paste(lig_genes, collapse = "/"),
        receptor_genes   = paste(rec_genes, collapse = "/"),
        pathway_name     = pw_name,
        annotation       = annot,
        n_lig_subunits   = length(lig_genes),
        n_rec_subunits   = length(rec_genes),
        stringsAsFactors = FALSE
      )
    }
    lr_pairs_db <- do.call(rbind, Filter(Negate(is.null), lr_pairs_list))
    rm(lr_pairs_list); invisible(gc())

    cat(paste0("  Resolved L-R pairs: ", nrow(lr_pairs_db), "\n"))
    write.csv(lr_pairs_db,
              file.path(tier4_dir, "CellChatDB_resolved_LR_pairs.csv"),
              row.names = FALSE)
  }

  compute_complex_expr <- function(gene_str, expr_mat, cells) {
    genes <- strsplit(gene_str, "/")[[1]]
    genes <- genes[genes %in% rownames(expr_mat)]
    if (length(genes) == 0) return(0)
    subunit_means <- sapply(genes, function(g) {
      mean(expr_mat[g, cells], na.rm = TRUE)
    })
    min(subunit_means)
  }

  zones_list <- list(
    Interface_to_ASPC = list(sender = "Interface", receiver = "ASPC_high"),
    Interface_to_AR   = list(sender = "Interface", receiver = "AR_high"),
    Interface_to_NEPC = list(sender = "Interface", receiver = "NEPC_high"),
    ASPC_to_Interface = list(sender = "ASPC_high", receiver = "Interface"),
    AR_to_Interface   = list(sender = "AR_high",   receiver = "Interface"),
    NEPC_to_Interface = list(sender = "NEPC_high", receiver = "Interface")
  )

  lr_results_all <- list()
  result_idx <- 1

  for (sn in successfully_loaded) {
    cat(paste0("  TIER4 — ", sn, "\n"))
    obj <- load_sample(sn)
    if (is.null(obj)) next
    md <- obj@meta.data

    if (!"interface_zone" %in% colnames(md)) { rm(obj); invisible(gc()); next }

    expr_mat <- safe_GetAssayData(obj, layer = "data")
    if (is.null(expr_mat)) { rm(obj); invisible(gc()); next }

    avail_genes <- rownames(expr_mat)

    for (zn in names(zones_list)) {
      zinfo <- zones_list[[zn]]
      sender_cells   <- colnames(obj)[md$interface_zone == zinfo$sender]
      receiver_cells <- colnames(obj)[md$interface_zone == zinfo$receiver]
      if (length(sender_cells) < 3 || length(receiver_cells) < 3) next

      for (i in seq_len(nrow(lr_pairs_db))) {
        if ("ligand_genes" %in% colnames(lr_pairs_db)) {
          lig_str  <- lr_pairs_db$ligand_genes[i]
          rec_str  <- lr_pairs_db$receptor_genes[i]
          pw       <- lr_pairs_db$pathway_name[i]
          annot    <- lr_pairs_db$annotation[i]
          int_name <- lr_pairs_db$interaction_name[i]
          lig_name <- lr_pairs_db$ligand_name[i]
          rec_name <- lr_pairs_db$receptor_name[i]
        } else {
          lig_str  <- lr_pairs_db$ligand[i]
          rec_str  <- lr_pairs_db$receptor[i]
          pw       <- lr_pairs_db$pathway_name[i]
          annot    <- if ("annotation" %in% colnames(lr_pairs_db)) lr_pairs_db$annotation[i] else "Unknown"
          int_name <- paste0(lig_str, "_", rec_str)
          lig_name <- lig_str
          rec_name <- rec_str
        }

        lig_genes_vec <- strsplit(lig_str, "/")[[1]]
        rec_genes_vec <- strsplit(rec_str, "/")[[1]]
        if (!any(lig_genes_vec %in% avail_genes) || !any(rec_genes_vec %in% avail_genes)) next

        lig_expr <- compute_complex_expr(lig_str, expr_mat, sender_cells)
        rec_expr <- compute_complex_expr(rec_str, expr_mat, receiver_cells)

        score <- sqrt(lig_expr * rec_expr)

        if (score > 0) {
          lr_results_all[[result_idx]] <- data.frame(
            sample           = sn,
            condition        = extract_condition(sn),
            zone_pair        = zn,
            interaction_name = int_name,
            ligand_name      = lig_name,
            receptor_name    = rec_name,
            pathway          = pw,
            annotation       = annot,
            lig_expr         = lig_expr,
            rec_expr         = rec_expr,
            score            = score,
            stringsAsFactors = FALSE
          )
          result_idx <- result_idx + 1
        }
      }
    }

    rm(obj, expr_mat); invisible(gc())
  }

  if (length(lr_results_all) > 0) {
    lr_df <- do.call(rbind, lr_results_all)
    lr_df <- lr_df[order(lr_df$score, decreasing = TRUE), ]
    write.csv(lr_df,
              file.path(tier4_dir, "LR_scores_CellChatDB_all_samples.csv"),
              row.names = FALSE)
    rm(lr_results_all); invisible(gc())

    lr_agg <- lr_df %>%
      group_by(condition, zone_pair, interaction_name, ligand_name,
               receptor_name, pathway, annotation) %>%
      summarise(
        mean_score = mean(score, na.rm = TRUE),
        mean_lig   = mean(lig_expr, na.rm = TRUE),
        mean_rec   = mean(rec_expr, na.rm = TRUE),
        n_samples  = n(),
        .groups    = "drop"
      ) %>%
      arrange(desc(mean_score))
    write.csv(lr_agg,
              file.path(tier4_dir, "LR_scores_CellChatDB_aggregated.csv"),
              row.names = FALSE)

    conditions_lr <- unique(lr_agg$condition)
    for (cond in conditions_lr) {
      cond_lr <- lr_agg[lr_agg$condition == cond & grepl("^Interface_to", lr_agg$zone_pair), ]
      if (nrow(cond_lr) == 0) next
      cond_lr$pair_label <- paste0(cond_lr$ligand_name, " -> ", cond_lr$receptor_name)
      cond_lr <- cond_lr[!duplicated(cond_lr$pair_label), ]
      top30 <- head(cond_lr, 30)
      top30$pair_label <- factor(top30$pair_label, levels = rev(top30$pair_label))
      p_lr <- ggplot(top30, aes(x = zone_pair, y = pair_label,
                                 size = mean_score, color = pathway)) +
        geom_point() +
        scale_size_continuous(range = c(2, 8)) +
        labs(title = paste0(cond, " — Top 30 CellChatDB L-R Pairs (Interface -> Zone)"),
             x = "Zone Pair", y = "L-R Pair", size = "Score", color = "Pathway") +
        theme_minimal(base_size = 9) +
        theme(axis.text.y = element_text(size = 7))
      ggsave(file.path(tier4_dir, paste0(cond, "_LR_CellChatDB_dotplot.pdf")),
             plot = p_lr, width = 12, height = 10)
    }

    pw_summary <- lr_agg %>%
      group_by(condition, pathway) %>%
      summarise(total_score = sum(mean_score, na.rm = TRUE), .groups = "drop")

    pw_conds <- unique(pw_summary$condition)
    pw_names <- unique(pw_summary$pathway)
    pw_mat   <- matrix(0, nrow = length(pw_names), ncol = length(pw_conds),
                        dimnames = list(pw_names, pw_conds))
    for (r in seq_len(nrow(pw_summary))) {
      pw_mat[pw_summary$pathway[r], pw_summary$condition[r]] <- pw_summary$total_score[r]
    }

    write.csv(pw_mat,
              file.path(tier4_dir, "CellChatDB_pathway_summary_by_condition.csv"),
              row.names = TRUE)

    if (have_pheatmap && nrow(pw_mat) > 1 && ncol(pw_mat) > 1) {
      pw_totals <- rowSums(pw_mat)
      top_pw    <- names(sort(pw_totals, decreasing = TRUE))[1:min(40, length(pw_totals))]
      pw_mat_sub <- pw_mat[top_pw, , drop = FALSE]
      pdf(file.path(tier4_dir, "CellChatDB_pathway_heatmap.pdf"),
          width = 10, height = 12)
      pheatmap(pw_mat_sub, scale = "row", cluster_cols = FALSE,
               main = "CellChatDB Pathway Scores by Condition",
               fontsize_row = 8)
      dev.off()
    }

    if ("annotation" %in% colnames(lr_agg)) {
      annot_summary <- lr_agg %>%
        group_by(condition, annotation) %>%
        summarise(total_score = sum(mean_score, na.rm = TRUE),
                  n_interactions = n(), .groups = "drop")
      write.csv(annot_summary,
                file.path(tier4_dir, "CellChatDB_annotation_summary.csv"),
                row.names = FALSE)

      p_annot <- ggplot(annot_summary, aes(x = condition, y = total_score,
                                            fill = annotation)) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(title = "CellChatDB — Interaction Types by Condition",
             x = "Condition", y = "Total Score", fill = "Annotation") +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(file.path(tier4_dir, "CellChatDB_annotation_barplot.pdf"),
             plot = p_annot, width = 10, height = 6)
    }

    cat(paste0("  CellChatDB L-R analysis complete: ", nrow(lr_df), " scored interactions\n"))
    rm(lr_df, lr_agg); invisible(gc())
  } else {
    cat("  No L-R interactions scored\n")
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

  # Use samples_with_dist from Tier 2 (these have interface zones + distances)
  # Fall back to successfully_loaded if samples_with_dist doesn't exist
  tier5_samples <- if (exists("samples_with_dist")) samples_with_dist else successfully_loaded
  conditions <- unique(sapply(tier5_samples, extract_condition))

  # --- Helper: get or compute variable features ---
  get_or_compute_vf <- function(obj, n_features = 2000) {
    img_bak <- obj@images
    obj@images <- list()

    for (assay_name in c("SCT", "RNA")) {
      if (!assay_name %in% Assays(obj)) next
      DefaultAssay(obj) <- assay_name
      vf <- tryCatch(VariableFeatures(obj), error = function(e) character(0))
      if (length(vf) > 0) { obj@images <- img_bak; return(head(vf, n_features)) }
    }

    # Compute on the fly
    DefaultAssay(obj) <- "RNA"
    obj <- tryCatch(FindVariableFeatures(obj, nfeatures = n_features, verbose = FALSE),
                    error = function(e) obj)
    vf <- tryCatch(VariableFeatures(obj), error = function(e) character(0))
    if (length(vf) > 0) { obj@images <- img_bak; return(head(vf, n_features)) }

    # Last resort: top-variance genes
    expr_mat <- tryCatch(GetAssayData(obj, layer = "data"), error = function(e) NULL)
    obj@images <- img_bak
    if (!is.null(expr_mat) && ncol(expr_mat) > 1) {
      gene_var <- apply(expr_mat, 1, var, na.rm = TRUE)
      return(head(names(sort(gene_var, decreasing = TRUE)), n_features))
    }
    return(character(0))
  }

  for (cond in conditions) {
    cat(paste0("\n  TIER5 — condition: ", cond, "\n"))
    cond_samples <- tier5_samples[sapply(tier5_samples, extract_condition) == cond]
    if (length(cond_samples) == 0) next

    # ============================================================
    # Phase 1: Identify near-interface spots and variable features
    # ============================================================
    near_iface_info <- list()  # sn -> list(cells)
    all_var_genes <- character(0)

    for (sn in cond_samples) {
      obj <- load_sample(sn)
      if (is.null(obj)) { cat(paste0("    ", sn, ": load failed\n")); next }
      md <- obj@meta.data

      # Determine near-interface spots using adaptive bins or distance
      near_cells <- character(0)

      if ("dist_bin_adaptive" %in% colnames(md)) {
        # Use the first 2-3 adaptive bins (interface + nearest)
        all_bins <- sort(unique(md$dist_bin_adaptive[!is.na(md$dist_bin_adaptive)]))
        # Always include interface; then take next 2 nearest bins
        near_bin_names <- c("bin_0_interface")
        other_bins <- setdiff(all_bins, "bin_0_interface")
        near_bin_names <- c(near_bin_names, head(other_bins, 2))
        near_cells <- colnames(obj)[!is.na(md$dist_bin_adaptive) &
                                      md$dist_bin_adaptive %in% near_bin_names]
        cat(paste0("    ", sn, ": using adaptive bins ",
                   paste(near_bin_names, collapse = ", "), "\n"))

      } else if ("dist_to_interface" %in% colnames(md)) {
        # Fallback: use distance threshold (bottom 60% of distances)
        d <- md$dist_to_interface
        thresh <- quantile(d[d > 0], 0.6, na.rm = TRUE)
        near_cells <- colnames(obj)[!is.na(d) & d <= thresh]
        cat(paste0("    ", sn, ": using distance threshold <= ",
                   round(thresh, 2), "\n"))

      } else if ("interface_zone" %in% colnames(md)) {
        # Fallback: use interface + neighboring zones
        near_cells <- colnames(obj)[md$interface_zone %in%
                                      c("Interface","ASPC_high","AR_high","NEPC_high")]
        cat(paste0("    ", sn, ": using zone-based selection\n"))

      } else {
        cat(paste0("    ", sn, ": no spatial data — skipping\n"))
        rm(obj); invisible(gc()); next
      }

      if (length(near_cells) < 10) {
        cat(paste0("    ", sn, ": only ", length(near_cells),
                   " near-interface spots — skipping\n"))
        rm(obj); invisible(gc()); next
      }

      # Get variable features
      vf <- get_or_compute_vf(obj, n_features = 2000)

      near_iface_info[[sn]] <- list(cells = near_cells)
      all_var_genes <- union(all_var_genes, vf)
      cat(paste0("    ", sn, ": ", length(near_cells), " spots, ",
                 length(vf), " VF\n"))

      rm(obj); invisible(gc())
    }

    cat(paste0("    Samples with near-interface spots: ",
               length(near_iface_info), "\n"))
    cat(paste0("    Union of variable genes: ", length(all_var_genes), "\n"))

    if (length(near_iface_info) == 0) {
      cat(paste0("    No usable samples for ", cond, " — skipping\n"))
      next
    }

    all_var_genes <- head(all_var_genes, 2000)
    if (length(all_var_genes) < 20) {
      cat(paste0("    Too few variable genes for ", cond, " — skipping\n"))
      next
    }

    # ============================================================
    # Phase 2: Collect expression matrices (one sample at a time)
    # ============================================================
    expr_list <- list()
    sample_gene_sets <- list()

    for (sn in names(near_iface_info)) {
      obj <- load_sample(sn)
      if (is.null(obj)) next

      expr_mat <- safe_GetAssayData(obj, layer = "data")
      if (is.null(expr_mat)) { rm(obj); invisible(gc()); next }

      cells <- near_iface_info[[sn]]$cells
      cells <- cells[cells %in% colnames(expr_mat)]
      if (length(cells) < 10) { rm(obj, expr_mat); invisible(gc()); next }

      genes_avail <- intersect(all_var_genes, rownames(expr_mat))
      if (length(genes_avail) < 20) { rm(obj, expr_mat); invisible(gc()); next }

      expr_list[[sn]] <- expr_mat[genes_avail, cells, drop = FALSE]
      sample_gene_sets[[sn]] <- genes_avail

      rm(obj, expr_mat); invisible(gc())
    }

    if (length(expr_list) == 0) {
      cat(paste0("    No expression data for ", cond, " — skipping\n"))
      next
    }

    # Align to common genes across all samples (prevents cbind mismatch)
    if (length(sample_gene_sets) >= 2) {
      common_vf <- Reduce(intersect, sample_gene_sets)
    } else {
      common_vf <- sample_gene_sets[[1]]
    }
    rm(sample_gene_sets); invisible(gc())

    cat(paste0("    Common variable features: ", length(common_vf), "\n"))
    if (length(common_vf) < 20) {
      cat(paste0("    Too few common VF for ", cond, " — skipping\n"))
      rm(expr_list); invisible(gc()); next
    }

    expr_combined <- do.call(cbind, lapply(expr_list, function(m) {
      m[common_vf, , drop = FALSE]
    }))
    total_spots <- ncol(expr_combined)
    rm(expr_list); invisible(gc())

    cat(paste0("    Combined: ", length(common_vf), " genes x ",
               total_spots, " spots\n"))

    expr_t <- t(as.matrix(expr_combined))

    # Remove zero-variance genes (PCA requires nonzero variance)
    gene_vars <- apply(expr_t, 2, var, na.rm = TRUE)
    nonzero <- gene_vars > 0 & !is.na(gene_vars)
    if (sum(nonzero) < 20) {
      cat(paste0("    Too few non-zero-variance genes for ", cond, "\n"))
      rm(expr_combined, expr_t); invisible(gc()); next
    }
    expr_t <- expr_t[, nonzero, drop = FALSE]
    # Keep expr_combined aligned
    expr_combined <- expr_combined[colnames(expr_t), , drop = FALSE]

    cat(paste0("    After variance filter: ", ncol(expr_t), " genes\n"))

    # ============================================================
    # Phase 3: PCA + Pseudotime
    # ============================================================
    n_pcs <- min(10, ncol(expr_t) - 1, nrow(expr_t) - 1)
    if (n_pcs < 2) {
      cat(paste0("    Not enough dims for PCA — skipping ", cond, "\n"))
      rm(expr_combined, expr_t); invisible(gc()); next
    }

    pca_res <- tryCatch(
      prcomp(expr_t, scale. = TRUE, center = TRUE, rank. = n_pcs),
      error = function(e) {
        cat(paste0("    PCA failed: ", conditionMessage(e), "\n"))
        NULL
      }
    )
    if (is.null(pca_res)) { rm(expr_combined, expr_t); invisible(gc()); next }

    pc_scores <- pca_res$x
    pca_var_explained <- summary(pca_res)$importance[2, ]  # proportion of variance
    cat(paste0("    PCA: ", ncol(pc_scores), " PCs, ",
               "PC1 var=", round(pca_var_explained[1] * 100, 1), "%\n"))
    rm(pca_res); invisible(gc())

    # Fit principal curve or use PC1
    pseudotime <- tryCatch({
      if (have_princurve) {
        n_fit <- min(5, ncol(pc_scores))
        pc_fit <- principal_curve(pc_scores[, 1:n_fit])
        pt <- pc_fit$lambda
        pt_r <- max(pt) - min(pt)
        if (pt_r > 0) (pt - min(pt)) / pt_r else rep(0.5, length(pt))
      } else {
        pt <- pc_scores[, 1]
        pt_r <- max(pt) - min(pt)
        if (pt_r > 0) (pt - min(pt)) / pt_r else rep(0.5, length(pt))
      }
    }, error = function(e) {
      cat(paste0("    Principal curve failed — using PC1: ",
                 conditionMessage(e), "\n"))
      pt <- pc_scores[, 1]
      pt_r <- max(pt) - min(pt)
      if (pt_r > 0) (pt - min(pt)) / pt_r else rep(0.5, length(pt))
    })

    cat(paste0("    Pseudotime: ", length(pseudotime), " spots, ",
               "range=[", round(min(pseudotime), 3), ", ",
               round(max(pseudotime), 3), "]\n"))

    # ============================================================
    # Phase 4: Collect metadata (reload one sample at a time)
    # ============================================================
    meta_list <- list()
    for (sn in names(near_iface_info)) {
      obj <- load_sample(sn)
      if (is.null(obj)) next
      cells <- near_iface_info[[sn]]$cells
      cells <- cells[cells %in% colnames(obj)]
      if (length(cells) == 0) { rm(obj); invisible(gc()); next }
      md <- obj@meta.data[cells, ]

      meta_list[[sn]] <- data.frame(
        barcode          = cells,
        sample_id        = sn,
        condition        = extract_condition(sn),
        interface_zone   = if ("interface_zone" %in% colnames(md)) md$interface_zone else NA,
        dist_to_interface = if ("dist_to_interface" %in% colnames(md)) md$dist_to_interface else NA_real_,
        dist_bin_adaptive = if ("dist_bin_adaptive" %in% colnames(md)) md$dist_bin_adaptive else NA,
        ASPC_score       = if ("ASPC_score" %in% colnames(md)) md$ASPC_score else NA_real_,
        AR_score         = if ("AR_score" %in% colnames(md)) md$AR_score else NA_real_,
        NEPC_UP_score    = if ("NEPC_UP_score" %in% colnames(md)) md$NEPC_UP_score else NA_real_,
        NEPC_Beltran_NET = if ("NEPC_Beltran_NET" %in% colnames(md)) md$NEPC_Beltran_NET else NA_real_,
        row = if ("row" %in% colnames(md)) md$row else NA_real_,
        col = if ("col" %in% colnames(md)) md$col else NA_real_,
        stringsAsFactors = FALSE
      )
      rm(obj); invisible(gc())
    }

    all_near_meta <- do.call(rbind, meta_list)
    rm(meta_list); invisible(gc())

    # ============================================================
    # Align pseudotime to metadata by barcode
    # ============================================================
    expr_barcodes <- colnames(expr_combined)
    meta_barcodes <- all_near_meta$barcode

    common_bcs <- intersect(expr_barcodes, meta_barcodes)
    cat(paste0("    Barcode alignment: expr=", length(expr_barcodes),
               ", meta=", nrow(all_near_meta),
               ", common=", length(common_bcs), "\n"))

    if (length(common_bcs) < 20) {
      cat(paste0("    Too few common barcodes — skipping ", cond, "\n"))
      rm(all_near_meta, expr_combined, expr_t, pc_scores, pseudotime)
      invisible(gc()); next
    }

    # Subset and align both
    all_near_meta <- all_near_meta[match(common_bcs, all_near_meta$barcode), ]
    pt_idx <- match(common_bcs, expr_barcodes)
    pseudotime <- pseudotime[pt_idx]
    pc_scores <- pc_scores[pt_idx, , drop = FALSE]
    expr_combined <- expr_combined[, common_bcs, drop = FALSE]

    all_near_meta$pseudotime <- pseudotime

    # Save data table
    write.csv(all_near_meta,
              file.path(tier5_dir, paste0(cond, "_pseudotime_data.csv")),
              row.names = FALSE)

    # ============================================================
    # Phase 5: Spearman correlations of scores with pseudotime
    # ============================================================
    score_cols_pt <- c("ASPC_score","AR_score","NEPC_UP_score","NEPC_Beltran_NET")
    cor_table <- lapply(score_cols_pt, function(sc) {
      v <- all_near_meta[[sc]]
      if (all(is.na(v))) return(NULL)
      ct <- tryCatch(
        cor.test(pseudotime, v, method = "spearman",
                 use = "complete.obs", exact = FALSE),
        error = function(e) NULL
      )
      if (is.null(ct)) return(NULL)
      data.frame(score = sc, rho = ct$estimate, pvalue = ct$p.value,
                 stringsAsFactors = FALSE)
    })
    cor_df <- do.call(rbind, Filter(Negate(is.null), cor_table))
    if (!is.null(cor_df) && nrow(cor_df) > 0) {
      write.csv(cor_df,
                file.path(tier5_dir, paste0(cond, "_pseudotime_score_correlations.csv")),
                row.names = FALSE)
      for (ri in seq_len(nrow(cor_df))) {
        cat(paste0("      ", cor_df$score[ri], ": rho=",
                   round(cor_df$rho[ri], 3), ", p=",
                   formatC(cor_df$pvalue[ri], format = "e", digits = 2), "\n"))
      }
    }

    # ============================================================
    # Phase 6: Transition genes (gene ~ pseudotime correlation)
    # ============================================================
    trans_sig <- NULL
    gene_pt_cor <- tryCatch({
      expr_for_cor <- as.matrix(expr_combined)
      n_genes <- nrow(expr_for_cor)
      cat(paste0("    Computing gene-pseudotime correlations for ",
                 n_genes, " genes...\n"))
      result <- matrix(NA_real_, nrow = 2, ncol = n_genes,
                       dimnames = list(c("rho","pvalue"), rownames(expr_for_cor)))
      for (gi in seq_len(n_genes)) {
        y <- expr_for_cor[gi, ]
        ct <- tryCatch(
          cor.test(pseudotime, y, method = "spearman",
                   use = "complete.obs", exact = FALSE),
          error = function(e) list(estimate = NA, p.value = NA)
        )
        result[1, gi] <- as.numeric(ct$estimate)
        result[2, gi] <- as.numeric(ct$p.value)
      }
      result
    }, error = function(e) {
      cat(paste0("    Gene-PT correlation error: ", conditionMessage(e), "\n"))
      NULL
    })

    if (!is.null(gene_pt_cor)) {
      trans_genes_df <- data.frame(
        gene   = colnames(gene_pt_cor),
        rho    = gene_pt_cor["rho", ],
        pvalue = gene_pt_cor["pvalue", ],
        stringsAsFactors = FALSE
      )
      trans_genes_df$padj <- p.adjust(trans_genes_df$pvalue, method = "BH")

      # Annotate signature membership
      trans_genes_df$in_ASPC    <- trans_genes_df$gene %in% ASPC_genes
      trans_genes_df$in_AR      <- trans_genes_df$gene %in% Hallmark_AR_genes
      trans_genes_df$in_NEPC_UP <- trans_genes_df$gene %in% NEPC_Beltran_UP
      trans_genes_df$in_NEPC_DN <- trans_genes_df$gene %in% NEPC_Beltran_DOWN

      # Save all
      write.csv(trans_genes_df,
                file.path(tier5_dir, paste0(cond, "_transition_genes_all.csv")),
                row.names = FALSE)

      # Filter significant
      trans_sig <- trans_genes_df[
        !is.na(trans_genes_df$rho) &
          abs(trans_genes_df$rho) > 0.2 &
          trans_genes_df$padj < 0.05, ]
      trans_sig <- trans_sig[order(abs(trans_sig$rho), decreasing = TRUE), ]

      write.csv(trans_sig,
                file.path(tier5_dir, paste0(cond, "_transition_genes_sig.csv")),
                row.names = FALSE)
      cat(paste0("    Transition genes (|rho|>0.2, padj<0.05): ",
                 nrow(trans_sig), " / ", nrow(trans_genes_df), "\n"))
    }

    # ============================================================
    # Phase 7: Multi-page PDF output
    # ============================================================
    cat(paste0("    Generating PDF for ", cond, "...\n"))

    pdf(file.path(tier5_dir, paste0(cond, "_TIER5_trajectory_results.pdf")),
        width = 11, height = 8.5)

    # ---- Page 1: Spatial pseudotime map ----
    if (all(c("col","row") %in% colnames(all_near_meta)) &&
        !all(is.na(all_near_meta$col))) {
      spot_sz <- compute_spot_size(all_near_meta[, c("col","row")])
      p_spatial <- ggplot(all_near_meta, aes(x = col, y = -row,
                                             color = pseudotime)) +
        geom_point(size = spot_sz) +
        scale_color_viridis_c(option = "plasma") +
        labs(title = paste0(cond, " — Spatial Pseudotime Map"),
             subtitle = paste0(nrow(all_near_meta), " spots from ",
                               length(unique(all_near_meta$sample_id)),
                               " samples"),
             x = "Array Col", y = "Array Row (inv)",
             color = "Pseudotime") +
        theme_minimal(base_size = 10)
      print(p_spatial)
      cat("      Page: spatial pseudotime map\n")
    }

    # ---- Page 2: PCA colored by pseudotime ----
    if (ncol(pc_scores) >= 2) {
      pca_df <- data.frame(
        PC1 = pc_scores[, 1], PC2 = pc_scores[, 2],
        pseudotime = pseudotime
      )
      p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = pseudotime)) +
        geom_point(size = 0.8, alpha = 0.6) +
        scale_color_viridis_c(option = "plasma") +
        labs(title = paste0(cond, " — PCA (colored by pseudotime)"),
             subtitle = paste0("PC1: ", round(pca_var_explained[1]*100,1),
                               "%, PC2: ", round(pca_var_explained[2]*100,1), "%"),
             x = "PC1", y = "PC2", color = "Pseudotime") +
        theme_minimal(base_size = 10)
      print(p_pca)

      # Also color by interface zone
      if ("interface_zone" %in% colnames(all_near_meta) &&
          !all(is.na(all_near_meta$interface_zone))) {
        pca_df$zone <- all_near_meta$interface_zone
        zone_colors <- c(Interface = "#E31A1C", ASPC_high = "#1F78B4",
                         AR_high = "#33A02C", NEPC_high = "#FF7F00",
                         Other = "#CCCCCC")
        p_pca_zone <- ggplot(pca_df, aes(x = PC1, y = PC2, color = zone)) +
          geom_point(size = 0.8, alpha = 0.5) +
          scale_color_manual(values = zone_colors) +
          labs(title = paste0(cond, " — PCA (colored by zone)"),
               x = "PC1", y = "PC2", color = "Zone") +
          theme_minimal(base_size = 10)
        print(p_pca_zone)
      }
      cat("      Page: PCA plots\n")
    }

    # ---- Page 3: Pseudotime distribution ----
    p_hist <- ggplot(all_near_meta, aes(x = pseudotime)) +
      geom_histogram(bins = 50, fill = "steelblue", color = "white",
                     alpha = 0.8) +
      labs(title = paste0(cond, " — Pseudotime Distribution"),
           x = "Pseudotime", y = "Count") +
      theme_minimal(base_size = 10)
    print(p_hist)

    # Per-sample density
    p_density <- ggplot(all_near_meta, aes(x = pseudotime,
                                           color = sample_id)) +
      geom_density(linewidth = 0.6, alpha = 0.5) +
      labs(title = paste0(cond, " — Pseudotime Density per Sample"),
           x = "Pseudotime", y = "Density") +
      theme_minimal(base_size = 10) +
      theme(legend.text = element_text(size = 5),
            legend.key.size = unit(0.3, "cm"))
    print(p_density)
    cat("      Page: pseudotime distribution\n")

    # ---- Page 4: Pseudotime vs distance to interface ----
    if ("dist_to_interface" %in% colnames(all_near_meta) &&
        !all(is.na(all_near_meta$dist_to_interface))) {
      rho_dist <- tryCatch({
        ct <- cor.test(all_near_meta$pseudotime,
                       all_near_meta$dist_to_interface,
                       method = "spearman", exact = FALSE)
        paste0("rho=", round(ct$estimate, 3),
               ", p=", formatC(ct$p.value, format = "e", digits = 2))
      }, error = function(e) "")

      p_pt_dist <- ggplot(all_near_meta,
                          aes(x = pseudotime, y = dist_to_interface)) +
        geom_point(size = 0.3, alpha = 0.15, color = "grey40") +
        geom_smooth(method = "loess", color = "blue", linewidth = 1) +
        labs(title = paste0(cond, " — Pseudotime vs Distance to Interface"),
             subtitle = rho_dist,
             x = "Pseudotime", y = "Distance to Interface") +
        theme_minimal(base_size = 10)
      print(p_pt_dist)
      cat("      Page: pseudotime vs distance\n")
    }

    # ---- Page 5: Pseudotime by adaptive distance bin ----
    if ("dist_bin_adaptive" %in% colnames(all_near_meta) &&
        !all(is.na(all_near_meta$dist_bin_adaptive))) {
      p_box_bin <- ggplot(all_near_meta,
                          aes(x = dist_bin_adaptive, y = pseudotime,
                              fill = dist_bin_adaptive)) +
        geom_boxplot(outlier.size = 0.3) +
        labs(title = paste0(cond, " — Pseudotime by Distance Bin"),
             x = "Distance Bin", y = "Pseudotime") +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
      print(p_box_bin)
      cat("      Page: pseudotime by bin\n")
    }

    # ---- Pages 6+: Score vs pseudotime scatters ----
    for (sc in score_cols_pt) {
      if (!sc %in% colnames(all_near_meta) ||
          all(is.na(all_near_meta[[sc]]))) next

      rho_str <- ""
      if (!is.null(cor_df) && nrow(cor_df) > 0 && sc %in% cor_df$score) {
        rv <- cor_df$rho[cor_df$score == sc]
        pv <- cor_df$pvalue[cor_df$score == sc]
        rho_str <- paste0("rho=", round(rv, 3),
                          ", p=", formatC(pv, format = "e", digits = 2))
      }

      p_sc <- ggplot(all_near_meta,
                     aes(x = pseudotime, y = .data[[sc]])) +
        geom_point(size = 0.3, alpha = 0.15, color = "grey40") +
        geom_smooth(method = "loess", color = "blue",
                    linewidth = 1.2, se = TRUE) +
        labs(title = paste0(cond, " — ", sc, " vs Pseudotime"),
             subtitle = rho_str,
             x = "Pseudotime", y = sc) +
        theme_minimal(base_size = 10)
      print(p_sc)
    }
    cat("      Page: score vs pseudotime scatters\n")

    # ---- Transition gene heatmap ----
    if (!is.null(trans_sig) && nrow(trans_sig) > 0) {
      top_trans <- head(trans_sig$gene, min(40, nrow(trans_sig)))
      top_trans <- top_trans[top_trans %in% rownames(expr_combined)]

      if (length(top_trans) >= 3) {
        pt_order <- order(pseudotime)
        heat_mat <- as.matrix(expr_combined[top_trans, pt_order, drop = FALSE])
        heat_scaled <- t(scale(t(heat_mat)))
        heat_scaled[heat_scaled > 2]  <-  2
        heat_scaled[heat_scaled < -2] <- -2
        heat_scaled[is.na(heat_scaled)] <- 0

        if (have_pheatmap) {
          # Pseudotime bin annotation
          pt_sorted <- pseudotime[pt_order]
          pt_bins <- cut(pt_sorted, breaks = 5, labels = FALSE)
          annot_col <- data.frame(
            PT_bin = factor(pt_bins),
            row.names = colnames(heat_mat)
          )
          pheatmap(heat_scaled,
                   cluster_cols = FALSE, cluster_rows = TRUE,
                   show_colnames = FALSE,
                   annotation_col = annot_col,
                   main = paste0(cond,
                                 " — Top Transition Genes Along Pseudotime"),
                   fontsize_row = 7)
        } else {
          heatmap(heat_scaled, Colv = NA, scale = "none", labCol = "",
                  main = paste0(cond, " — Transition Genes"))
        }
        cat("      Page: transition gene heatmap\n")
      }
    }

    # ---- Top 12 transition gene individual curves ----
    if (!is.null(trans_sig) && nrow(trans_sig) > 0) {
      top12 <- head(trans_sig$gene, min(12, nrow(trans_sig)))
      top12 <- top12[top12 %in% rownames(expr_combined)]

      if (length(top12) > 0) {
        gene_plots <- list()
        for (gg in top12) {
          gdf <- data.frame(
            pseudotime = pseudotime,
            expr = as.numeric(expr_combined[gg, ])
          )
          rho_g <- trans_sig$rho[trans_sig$gene == gg]
          padj_g <- trans_sig$padj[trans_sig$gene == gg]
          sub_txt <- paste0("rho=", round(rho_g, 3),
                            "  padj=", formatC(padj_g, format = "e", digits = 2))

          gp <- ggplot(gdf, aes(x = pseudotime, y = expr)) +
            geom_point(size = 0.3, alpha = 0.2, color = "grey50") +
            geom_smooth(method = "loess", color = "red",
                        linewidth = 0.8, se = TRUE) +
            labs(title = gg, subtitle = sub_txt,
                 x = "Pseudotime", y = "Expression") +
            theme_minimal(base_size = 8)
          gene_plots[[gg]] <- gp
        }

        n_per_page <- 6
        for (i in seq(1, length(gene_plots), by = n_per_page)) {
          sub <- gene_plots[i:min(i + n_per_page - 1, length(gene_plots))]
          print(wrap_plots(sub, ncol = 3) +
                  plot_annotation(
                    title = paste0(cond,
                                   " — Transition Gene Expression vs Pseudotime")))
        }
        cat("      Page: individual gene curves\n")
      }
    }

    # ---- Transition gene volcano-style plot ----
    if (!is.null(gene_pt_cor) && exists("trans_genes_df")) {
      tg <- trans_genes_df
      tg$sig <- !is.na(tg$padj) & tg$padj < 0.05 & abs(tg$rho) > 0.2
      tg$label <- ""
      top_lbl <- head(tg$gene[tg$sig][order(abs(tg$rho[tg$sig]),
                                            decreasing = TRUE)], 15)
      tg$label[tg$gene %in% top_lbl] <- tg$gene[tg$gene %in% top_lbl]

      p_rho <- ggplot(tg, aes(x = rho, y = -log10(pvalue + 1e-300),
                               color = sig)) +
        geom_point(size = 0.8, alpha = 0.5) +
        scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
        geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed",
                   color = "grey40") +
        labs(title = paste0(cond, " — Gene-Pseudotime Correlation"),
             subtitle = paste0(sum(tg$sig), " significant transition genes"),
             x = "Spearman rho", y = "-log10(p-value)", color = "Sig") +
        theme_minimal(base_size = 10)

      if (requireNamespace("ggrepel", quietly = TRUE) &&
          sum(nchar(tg$label) > 0) > 0) {
        p_rho <- p_rho +
          ggrepel::geom_text_repel(aes(label = label),
                                   size = 2.5, color = "black",
                                   max.overlaps = 15)
      }
      print(p_rho)
      cat("      Page: transition gene volcano\n")
    }

    # ---- Signature gene enrichment in transition genes ----
    if (!is.null(trans_sig) && nrow(trans_sig) > 0) {
      sig_lists <- list(ASPC = ASPC_genes, AR = Hallmark_AR_genes,
                        NEPC_UP = NEPC_Beltran_UP, NEPC_DOWN = NEPC_Beltran_DOWN)
      enrich_data <- lapply(names(sig_lists), function(sl) {
        sl_g <- sig_lists[[sl]]
        n_trans <- sum(trans_sig$gene %in% sl_g)
        n_up    <- sum(trans_sig$gene %in% sl_g & trans_sig$rho > 0)
        n_down  <- sum(trans_sig$gene %in% sl_g & trans_sig$rho < 0)
        data.frame(signature = sl, n_transition = n_trans,
                   n_increasing = n_up, n_decreasing = n_down,
                   stringsAsFactors = FALSE)
      })
      enrich_df <- do.call(rbind, enrich_data)
      write.csv(enrich_df,
                file.path(tier5_dir,
                          paste0(cond, "_transition_signature_enrichment.csv")),
                row.names = FALSE)

      # Bar plot
      enrich_long <- data.frame(
        signature = rep(enrich_df$signature, 2),
        direction = rep(c("Increasing","Decreasing"), each = nrow(enrich_df)),
        count = c(enrich_df$n_increasing, enrich_df$n_decreasing),
        stringsAsFactors = FALSE
      )
      p_enrich <- ggplot(enrich_long,
                         aes(x = signature, y = count, fill = direction)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_fill_manual(values = c(Increasing = "firebrick",
                                     Decreasing = "steelblue")) +
        labs(title = paste0(cond,
                            " — Signature Genes Among Transition Genes"),
             x = "Signature", y = "Count", fill = "Direction") +
        theme_minimal(base_size = 10)
      print(p_enrich)
      cat("      Page: signature enrichment\n")
    }

    dev.off()
    cat(paste0("    PDF saved: ", cond, "_TIER5_trajectory_results.pdf\n"))

    # ---- Also save standalone spatial pseudotime PDF ----
    if (all(c("col","row") %in% colnames(all_near_meta)) &&
        !all(is.na(all_near_meta$col))) {
      spot_sz <- compute_spot_size(all_near_meta[, c("col","row")])
      p_pt_standalone <- ggplot(all_near_meta,
                                aes(x = col, y = -row,
                                    color = pseudotime)) +
        geom_point(size = spot_sz) +
        scale_color_viridis_c(option = "plasma") +
        labs(title = paste0(cond, " — Spatial Pseudotime"),
             x = "Array Col", y = "Array Row (inv)",
             color = "Pseudotime") +
        theme_minimal(base_size = 10)
      ggsave(file.path(tier5_dir, paste0(cond, "_spatial_pseudotime.pdf")),
             plot = p_pt_standalone, width = 10, height = 8)
    }

    # Clean up
    rm(all_near_meta, expr_combined, expr_t, pc_scores, pseudotime,
       near_iface_info)
    if (exists("trans_sig", inherits = FALSE)) rm(trans_sig)
    if (exists("trans_genes_df", inherits = FALSE)) rm(trans_genes_df)
    if (exists("gene_pt_cor", inherits = FALSE)) rm(gene_pt_cor)
    if (exists("cor_df", inherits = FALSE)) rm(cor_df)
    invisible(gc())
  }

  tier5_status <- "COMPLETE"
  cat("\nTIER 5 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 5 ERROR: ", conditionMessage(e), "\n"))
  tier5_status <<- paste0("ERROR: ", conditionMessage(e))
  tryCatch(dev.off(), error = function(e2) NULL)
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

  all_sig_genes_grn <- unique(c(ASPC_genes, Hallmark_AR_genes, NEPC_Beltran_UP, NEPC_Beltran_DOWN))

  cat("  Collecting interface expression...\n")
  iface_expr_list <- list()
  iface_var_genes <- list()

  for (sn in successfully_loaded) {
    obj <- load_sample(sn)
    if (is.null(obj)) next
    if (!"interface_zone" %in% colnames(obj@meta.data)) { rm(obj); invisible(gc()); next }

    iface_cells <- colnames(obj)[obj$interface_zone == "Interface"]
    if (length(iface_cells) < 3) { rm(obj); invisible(gc()); next }

    expr_mat <- safe_GetAssayData(obj, layer = "data")
    if (is.null(expr_mat)) { rm(obj); invisible(gc()); next }

    sig_avail <- intersect(all_sig_genes_grn, rownames(expr_mat))
    iface_sub <- expr_mat[, iface_cells, drop = FALSE]
    gene_var  <- apply(iface_sub, 1, var, na.rm = TRUE)
    top_var   <- names(sort(gene_var, decreasing = TRUE))[1:min(200, length(gene_var))]

    keep_genes <- unique(c(sig_avail, top_var))
    iface_expr_list[[sn]] <- iface_sub[keep_genes, , drop = FALSE]
    iface_var_genes[[sn]] <- top_var

    rm(obj, expr_mat, iface_sub, gene_var); invisible(gc())
  }

  if (length(iface_expr_list) == 0) stop("No interface spots found for GRN")

  common_genes_grn <- Reduce(intersect, lapply(iface_expr_list, rownames))
  if (length(common_genes_grn) < 20) stop("Too few common genes for GRN")

  all_var_union <- unique(unlist(iface_var_genes))
  grn_genes <- unique(c(intersect(all_sig_genes_grn, common_genes_grn),
                         intersect(all_var_union, common_genes_grn)))
  grn_genes <- head(grn_genes, 600)

  cat(paste0("  GRN genes: ", length(grn_genes), "\n"))

  all_expr_iface <- do.call(cbind, lapply(iface_expr_list, function(m) m[grn_genes, ]))

  # *** FIX: compute n_spots_grn BEFORE removing all_expr_iface ***
  n_spots_grn <- ncol(all_expr_iface)
  cat(paste0("  Total interface spots for GRN: ", n_spots_grn, "\n"))

  rm(iface_expr_list, iface_var_genes); invisible(gc())

  # Spearman correlation matrix
  cat("  Computing Spearman correlation matrix...\n")
  cor_mat <- tryCatch(
    cor(t(as.matrix(all_expr_iface)), method = "spearman", use = "pairwise.complete.obs"),
    error = function(e) NULL
  )

  # *** FIX: remove all_expr_iface AFTER correlation and n_spots_grn are computed ***
  rm(all_expr_iface); invisible(gc())

  if (is.null(cor_mat)) stop("Correlation computation failed")

  # Filter edges: |rho| > 0.3 and p < 0.01
  t_stat  <- cor_mat * sqrt((n_spots_grn - 2) / (1 - cor_mat^2 + 1e-10))
  p_mat   <- 2 * pt(-abs(t_stat), df = max(n_spots_grn - 2, 1))

  idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
  edge_df <- data.frame(
    gene1  = rownames(cor_mat)[idx[, 1]],
    gene2  = colnames(cor_mat)[idx[, 2]],
    rho    = cor_mat[idx],
    pvalue = p_mat[idx],
    stringsAsFactors = FALSE
  )
  edge_df <- edge_df[abs(edge_df$rho) > 0.3 & edge_df$pvalue < 0.01, ]
  cat(paste0("  Edges after filtering: ", nrow(edge_df), "\n"))

  write.csv(edge_df, file.path(tier6_dir, "edge_list.csv"), row.names = FALSE)
  rm(cor_mat, t_stat, p_mat, idx); invisible(gc())

  # igraph network analysis
  hub_genes <- character(0)
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
    write.csv(data.frame(hub_gene = hub_genes,
                          degree = node_df$degree[1:min(20, nrow(node_df))]),
              file.path(tier6_dir, "hub_genes.csv"), row.names = FALSE)
    cat(paste0("  Hub genes: ", paste(head(hub_genes, 10), collapse = ", "), "\n"))

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
    rm(g, g_sub, node_df); invisible(gc())
  } else {
    cat("  igraph not available or no edges — skipping network analysis\n")
  }

  rm(edge_df); invisible(gc())
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

  # --- Hallmark gene sets for pathway scoring ---
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
                          "CTNNB1","APC","AXIN1","AXIN2","TCF7","LEF1","MYC"),
    PI3K_AKT_mTOR     = c("PIK3CA","PIK3CB","PIK3R1","AKT1","AKT2","MTOR","RPTOR",
                          "RICTOR","RPS6KB1","EIF4EBP1","PTEN","TSC1","TSC2","RHEB"),
    Interferon_alpha   = c("IFNA1","IFNAR1","IFNAR2","STAT1","STAT2","IRF9","ISG15",
                           "MX1","MX2","OAS1","OAS2","IFIT1","IFIT2","IFI27"),
    Interferon_gamma   = c("IFNG","IFNGR1","IFNGR2","STAT1","IRF1","IRF8","GBP1",
                           "GBP2","CXCL9","CXCL10","CXCL11","IDO1","CIITA","TAP1"),
    TGFb              = c("TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","SMAD2","SMAD3",
                          "SMAD4","SMAD7","SERPINE1","CTGF","COL1A1","FN1","ACTA2"),
    NF_kB             = c("NFKB1","NFKB2","RELA","RELB","REL","NFKBIA","NFKBIB",
                          "IKBKB","IKBKG","CHUK","TNFAIP3","BCL2","BIRC3","ICAM1"),
    Neuroendocrine    = c("SYP","CHGA","CHGB","ENO2","NCAM1","ASCL1","NEUROD1",
                          "SOX2","POU3F2","REST","INSM1","SCG3","BRN2","FOXA2")
  )

  # ============================================================
  # 7A: Per-sample summary table
  # ============================================================
  cat("  7A: Building per-sample summary...\n")
  sample_summary_list <- list()

  for (sn in successfully_loaded) {
    obj <- load_sample(sn)
    if (is.null(obj)) next
    md <- obj@meta.data

    n_spots <- ncol(obj)
    cond    <- extract_condition(sn)

    # Zone counts
    zone_counts <- list(
      Interface = 0, ASPC_high = 0, AR_high = 0, NEPC_high = 0, Other = 0
    )
    if ("interface_zone" %in% colnames(md)) {
      zt <- table(md$interface_zone)
      for (zn in names(zt)) zone_counts[[zn]] <- as.integer(zt[zn])
    }

    # Score medians
    score_medians <- list(
      ASPC_score_med       = if ("ASPC_score" %in% colnames(md)) median(md$ASPC_score, na.rm = TRUE) else NA,
      AR_score_med         = if ("AR_score" %in% colnames(md)) median(md$AR_score, na.rm = TRUE) else NA,
      NEPC_UP_score_med    = if ("NEPC_UP_score" %in% colnames(md)) median(md$NEPC_UP_score, na.rm = TRUE) else NA,
      NEPC_DOWN_score_med  = if ("NEPC_DOWN_score" %in% colnames(md)) median(md$NEPC_DOWN_score, na.rm = TRUE) else NA,
      NEPC_Beltran_NET_med = if ("NEPC_Beltran_NET" %in% colnames(md)) median(md$NEPC_Beltran_NET, na.rm = TRUE) else NA
    )

    row <- c(
      list(sample = sn, condition = cond, n_spots = n_spots),
      zone_counts,
      score_medians
    )
    sample_summary_list[[sn]] <- as.data.frame(row, stringsAsFactors = FALSE)

    rm(obj, md); invisible(gc())
  }

  if (length(sample_summary_list) > 0) {
    sample_summary <- do.call(rbind, sample_summary_list)
    write.csv(sample_summary,
              file.path(tier7_dir, "sample_summary_table.csv"),
              row.names = FALSE)
    cat(paste0("    Sample summary: ", nrow(sample_summary), " samples\n"))
  } else {
    sample_summary <- data.frame()
    cat("    WARNING: no samples for summary\n")
  }

  # ============================================================
  # 7B: Hallmark pathway scoring per condition per zone
  # ============================================================
  cat("  7B: Hallmark pathway scoring per condition x zone...\n")
  pw_zone_results <- list()

  conditions_meta <- unique(sapply(successfully_loaded, extract_condition))
  zone_names <- c("Interface","ASPC_high","AR_high","NEPC_high","Other")

  for (cond in conditions_meta) {
    cond_samples <- successfully_loaded[sapply(successfully_loaded, extract_condition) == cond]

    for (zn in zone_names) {
      # Collect mean expression across all spots of this zone across samples
      zone_means_list <- list()

      for (sn in cond_samples) {
        obj <- load_sample(sn)
        if (is.null(obj)) next
        if (!"interface_zone" %in% colnames(obj@meta.data)) { rm(obj); invisible(gc()); next }

        zone_cells <- colnames(obj)[obj$interface_zone == zn]
        if (length(zone_cells) < 3) { rm(obj); invisible(gc()); next }

        expr_mat <- safe_GetAssayData(obj, layer = "data")
        if (is.null(expr_mat)) { rm(obj); invisible(gc()); next }

        zone_mean <- rowMeans(expr_mat[, zone_cells, drop = FALSE])
        zone_means_list[[sn]] <- zone_mean

        rm(obj, expr_mat); invisible(gc())
      }

      if (length(zone_means_list) == 0) next

      # Average across samples
      zone_mean_combined <- rowMeans(do.call(cbind, zone_means_list), na.rm = TRUE)

      # Score each hallmark set
      for (pw_name in names(hallmark_sets)) {
        pw_genes <- hallmark_sets[[pw_name]]
        pw_avail <- intersect(pw_genes, names(zone_mean_combined))
        if (length(pw_avail) < 3) next

        pw_score <- mean(zone_mean_combined[pw_avail], na.rm = TRUE)
        pw_zone_results[[length(pw_zone_results) + 1]] <- data.frame(
          condition    = cond,
          zone         = zn,
          pathway      = pw_name,
          score        = pw_score,
          n_genes_used = length(pw_avail),
          n_samples    = length(zone_means_list),
          stringsAsFactors = FALSE
        )
      }

      rm(zone_means_list, zone_mean_combined); invisible(gc())
    }
  }

  if (length(pw_zone_results) > 0) {
    pw_zone_df <- do.call(rbind, pw_zone_results)
    write.csv(pw_zone_df,
              file.path(tier7_dir, "hallmark_pathway_scores_by_condition_zone.csv"),
              row.names = FALSE)
    cat(paste0("    Pathway x zone scores: ", nrow(pw_zone_df), " rows\n"))

    # Heatmap: pathway x condition (Interface zone only)
    iface_pw <- pw_zone_df[pw_zone_df$zone == "Interface", ]
    if (nrow(iface_pw) > 0) {
      pw_names_avail <- unique(iface_pw$pathway)
      conds_avail    <- unique(iface_pw$condition)
      iface_mat <- matrix(0, nrow = length(pw_names_avail), ncol = length(conds_avail),
                          dimnames = list(pw_names_avail, conds_avail))
      for (r in seq_len(nrow(iface_pw))) {
        iface_mat[iface_pw$pathway[r], iface_pw$condition[r]] <- iface_pw$score[r]
      }

      if (have_pheatmap && nrow(iface_mat) > 1 && ncol(iface_mat) > 1) {
        pdf(file.path(tier7_dir, "hallmark_interface_heatmap.pdf"),
            width = 10, height = 12)
        pheatmap(iface_mat, scale = "row", cluster_cols = FALSE,
                 main = "Hallmark Pathway Scores — Interface Zone by Condition",
                 fontsize_row = 8, fontsize_col = 10)
        dev.off()
        cat("    Interface hallmark heatmap saved\n")
      }
    }

    # Dot plot: top pathways across conditions for Interface zone
    if (nrow(iface_pw) > 0) {
      p_pw_dot <- ggplot(iface_pw,
                          aes(x = condition, y = reorder(pathway, score),
                              size = score, color = score)) +
        geom_point() +
        scale_color_viridis_c(option = "inferno") +
        scale_size_continuous(range = c(2, 8)) +
        labs(title = "Hallmark Pathways — Interface Zone",
             x = "Condition", y = "Pathway", size = "Score", color = "Score") +
        theme_minimal(base_size = 9) +
        theme(axis.text.y = element_text(size = 7))
      ggsave(file.path(tier7_dir, "hallmark_interface_dotplot.pdf"),
             plot = p_pw_dot, width = 10, height = 12)
    }

    rm(pw_zone_results, pw_zone_df); invisible(gc())
  }

  # ============================================================
  # 7C: Score distributions across conditions (boxplots)
  # ============================================================
  cat("  7C: Score distributions across conditions...\n")

  score_dist_list <- list()
  for (sn in successfully_loaded) {
    obj <- load_sample(sn)
    if (is.null(obj)) next
    md <- obj@meta.data
    cond <- extract_condition(sn)

    score_cols_avail <- intersect(
      c("ASPC_score","AR_score","NEPC_UP_score","NEPC_DOWN_score","NEPC_Beltran_NET"),
      colnames(md)
    )
    if (length(score_cols_avail) == 0) { rm(obj); invisible(gc()); next }

    for (sc in score_cols_avail) {
      score_dist_list[[length(score_dist_list) + 1]] <- data.frame(
        sample    = sn,
        condition = cond,
        score_name = sc,
        value     = md[[sc]],
        stringsAsFactors = FALSE
      )
    }
    rm(obj, md); invisible(gc())
  }

  if (length(score_dist_list) > 0) {
    score_dist_df <- do.call(rbind, score_dist_list)
    rm(score_dist_list); invisible(gc())

    # Boxplots per score faceted by score_name
    p_box <- ggplot(score_dist_df, aes(x = condition, y = value, fill = condition)) +
      geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.3) +
      facet_wrap(~ score_name, scales = "free_y", ncol = 2) +
      labs(title = "Signature Score Distributions by Condition",
           x = "Condition", y = "Score") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    ggsave(file.path(tier7_dir, "score_distributions_by_condition.pdf"),
           plot = p_box, width = 12, height = 10)
    cat("    Score distribution boxplots saved\n")

    # Violin plots
    p_violin <- ggplot(score_dist_df, aes(x = condition, y = value, fill = condition)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
      geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white", alpha = 0.8) +
      facet_wrap(~ score_name, scales = "free_y", ncol = 2) +
      labs(title = "Signature Scores — Violin Plots",
           x = "Condition", y = "Score") +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    ggsave(file.path(tier7_dir, "score_violin_by_condition.pdf"),
           plot = p_violin, width = 12, height = 10)

    rm(score_dist_df); invisible(gc())
  }

  # ============================================================
  # 7D: Kruskal-Wallis tests across conditions per zone
  # ============================================================
  cat("  7D: Statistical tests across conditions...\n")
  kw_results <- list()

  score_names_test <- c("ASPC_score","AR_score","NEPC_UP_score",
                         "NEPC_DOWN_score","NEPC_Beltran_NET")

  for (zn in zone_names) {
    # Collect per-sample median scores for each condition
    sample_zone_scores <- list()
    for (sn in successfully_loaded) {
      obj <- load_sample(sn)
      if (is.null(obj)) next
      md <- obj@meta.data
      cond <- extract_condition(sn)

      if ("interface_zone" %in% colnames(md)) {
        zone_cells <- which(md$interface_zone == zn)
      } else {
        zone_cells <- seq_len(nrow(md))
      }
      if (length(zone_cells) < 3) { rm(obj); invisible(gc()); next }

      for (sc in score_names_test) {
        if (!sc %in% colnames(md)) next
        med_val <- median(md[[sc]][zone_cells], na.rm = TRUE)
        sample_zone_scores[[length(sample_zone_scores) + 1]] <- data.frame(
          sample    = sn,
          condition = cond,
          zone      = zn,
          score     = sc,
          median_value = med_val,
          stringsAsFactors = FALSE
        )
      }
      rm(obj, md); invisible(gc())
    }

    if (length(sample_zone_scores) == 0) next
    szs_df <- do.call(rbind, sample_zone_scores)

    for (sc in score_names_test) {
      sc_data <- szs_df[szs_df$score == sc & !is.na(szs_df$median_value), ]
      if (nrow(sc_data) < 3 || length(unique(sc_data$condition)) < 2) next

      kw_test <- tryCatch(
        kruskal.test(median_value ~ condition, data = sc_data),
        error = function(e) NULL
      )
      if (is.null(kw_test)) next

      kw_results[[length(kw_results) + 1]] <- data.frame(
        zone    = zn,
        score   = sc,
        statistic = kw_test$statistic,
        pvalue  = kw_test$p.value,
        n_groups = length(unique(sc_data$condition)),
        n_samples = nrow(sc_data),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(kw_results) > 0) {
    kw_df <- do.call(rbind, kw_results)
    kw_df$padj <- p.adjust(kw_df$pvalue, method = "BH")
    kw_df <- kw_df[order(kw_df$padj), ]
    write.csv(kw_df,
              file.path(tier7_dir, "kruskal_wallis_tests.csv"),
              row.names = FALSE)
    cat(paste0("    KW tests: ", nrow(kw_df), " tests, ",
               sum(kw_df$padj < 0.05), " significant (padj < 0.05)\n"))
  }

  # ============================================================
  # 7E: Pairwise condition comparisons (Wilcoxon) for Interface zone
  # ============================================================
  cat("  7E: Pairwise Wilcoxon tests for Interface zone...\n")
  pw_results <- list()

  # Collect Interface-zone median scores per sample
  iface_sample_scores <- list()
  for (sn in successfully_loaded) {
    obj <- load_sample(sn)
    if (is.null(obj)) next
    md <- obj@meta.data
    cond <- extract_condition(sn)

    if (!"interface_zone" %in% colnames(md)) { rm(obj); invisible(gc()); next }
    iface_cells <- which(md$interface_zone == "Interface")
    if (length(iface_cells) < 3) { rm(obj); invisible(gc()); next }

    for (sc in score_names_test) {
      if (!sc %in% colnames(md)) next
      med_val <- median(md[[sc]][iface_cells], na.rm = TRUE)
      iface_sample_scores[[length(iface_sample_scores) + 1]] <- data.frame(
        sample = sn, condition = cond, score = sc, median_value = med_val,
        stringsAsFactors = FALSE
      )
    }
    rm(obj, md); invisible(gc())
  }

  if (length(iface_sample_scores) > 0) {
    iss_df <- do.call(rbind, iface_sample_scores)

    cond_pairs <- combn(unique(iss_df$condition), 2, simplify = FALSE)

    for (sc in score_names_test) {
      for (cp in cond_pairs) {
        v1 <- iss_df$median_value[iss_df$score == sc & iss_df$condition == cp[1]]
        v2 <- iss_df$median_value[iss_df$score == sc & iss_df$condition == cp[2]]
        if (length(v1) < 2 || length(v2) < 2) next

        wt <- tryCatch(wilcox.test(v1, v2), error = function(e) NULL)
        if (is.null(wt)) next

        pw_results[[length(pw_results) + 1]] <- data.frame(
          score       = sc,
          condition_1 = cp[1],
          condition_2 = cp[2],
          W_statistic = wt$statistic,
          pvalue      = wt$p.value,
          median_1    = median(v1, na.rm = TRUE),
          median_2    = median(v2, na.rm = TRUE),
          n1          = length(v1),
          n2          = length(v2),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(pw_results) > 0) {
    pw_df <- do.call(rbind, pw_results)
    pw_df$padj <- p.adjust(pw_df$pvalue, method = "BH")
    pw_df <- pw_df[order(pw_df$padj), ]
    write.csv(pw_df,
              file.path(tier7_dir, "pairwise_wilcoxon_interface.csv"),
              row.names = FALSE)
    cat(paste0("    Pairwise tests: ", nrow(pw_df), " comparisons, ",
               sum(pw_df$padj < 0.05), " significant\n"))
  }

  # ============================================================
  # 7F: Interface zone proportion analysis
  # ============================================================
  cat("  7F: Interface zone proportion analysis...\n")

  if (nrow(sample_summary) > 0 && "Interface" %in% colnames(sample_summary)) {
    sample_summary$Interface_prop <- sample_summary$Interface / sample_summary$n_spots

    p_iface_prop <- ggplot(sample_summary,
                            aes(x = condition, y = Interface_prop, fill = condition)) +
      geom_boxplot(outlier.size = 1) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
      labs(title = "Interface Zone Proportion by Condition",
           x = "Condition", y = "Proportion of Spots in Interface Zone") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none")
    ggsave(file.path(tier7_dir, "interface_proportion_by_condition.pdf"),
           plot = p_iface_prop, width = 8, height = 6)
    cat("    Interface proportion plot saved\n")
  }

  # ============================================================
  # 7G: Score correlation across conditions (scatter matrix)
  # ============================================================
  cat("  7G: Score correlation matrix across conditions...\n")

  if (nrow(sample_summary) > 0) {
    score_med_cols <- intersect(
      c("ASPC_score_med","AR_score_med","NEPC_UP_score_med",
        "NEPC_DOWN_score_med","NEPC_Beltran_NET_med"),
      colnames(sample_summary)
    )
    if (length(score_med_cols) >= 2) {
      score_med_mat <- as.matrix(sample_summary[, score_med_cols])
      score_med_mat <- score_med_mat[complete.cases(score_med_mat), , drop = FALSE]

      if (nrow(score_med_mat) > 3) {
        cor_scores <- cor(score_med_mat, method = "spearman", use = "pairwise.complete.obs")
        write.csv(cor_scores,
                  file.path(tier7_dir, "score_correlation_matrix.csv"),
                  row.names = TRUE)

        if (have_pheatmap) {
          pdf(file.path(tier7_dir, "score_correlation_heatmap.pdf"),
              width = 8, height = 7)
          pheatmap(cor_scores,
                   display_numbers = TRUE, number_format = "%.2f",
                   main = "Spearman Correlation of Median Scores Across Samples",
                   fontsize = 10)
          dev.off()
          cat("    Score correlation heatmap saved\n")
        }
      }
    }
  }

  # ============================================================
  # 7H: NEPC vs AR scatter colored by condition
  # ============================================================
  cat("  7H: NEPC vs AR scatter plot...\n")

  if (nrow(sample_summary) > 0 &&
      "AR_score_med" %in% colnames(sample_summary) &&
      "NEPC_Beltran_NET_med" %in% colnames(sample_summary)) {

    ss_plot <- sample_summary[!is.na(sample_summary$AR_score_med) &
                               !is.na(sample_summary$NEPC_Beltran_NET_med), ]

    if (nrow(ss_plot) > 0) {
      p_ar_nepc <- ggplot(ss_plot,
                           aes(x = AR_score_med, y = NEPC_Beltran_NET_med,
                               color = condition, label = sample)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        labs(title = "AR Score vs NEPC Beltran NET Score (per-sample medians)",
             x = "Median AR Score", y = "Median NEPC Beltran NET Score",
             color = "Condition") +
        theme_minimal(base_size = 10)
      ggsave(file.path(tier7_dir, "AR_vs_NEPC_scatter.pdf"),
             plot = p_ar_nepc, width = 10, height = 8)
      cat("    AR vs NEPC scatter saved\n")
    }
  }

  # ============================================================
  # 7I: ASPC vs NEPC scatter colored by condition
  # ============================================================
  cat("  7I: ASPC vs NEPC scatter plot...\n")

  if (nrow(sample_summary) > 0 &&
      "ASPC_score_med" %in% colnames(sample_summary) &&
      "NEPC_Beltran_NET_med" %in% colnames(sample_summary)) {

    ss_plot2 <- sample_summary[!is.na(sample_summary$ASPC_score_med) &
                                !is.na(sample_summary$NEPC_Beltran_NET_med), ]

    if (nrow(ss_plot2) > 0) {
      p_aspc_nepc <- ggplot(ss_plot2,
                             aes(x = ASPC_score_med, y = NEPC_Beltran_NET_med,
                                 color = condition)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        labs(title = "ASPC Score vs NEPC Beltran NET Score (per-sample medians)",
             x = "Median ASPC Score", y = "Median NEPC Beltran NET Score",
             color = "Condition") +
        theme_minimal(base_size = 10)
      ggsave(file.path(tier7_dir, "ASPC_vs_NEPC_scatter.pdf"),
             plot = p_aspc_nepc, width = 10, height = 8)
      cat("    ASPC vs NEPC scatter saved\n")
    }
  }

  # ============================================================
  # 7J: Comprehensive multi-page PDF summary
  # ============================================================
  cat("  7J: Generating comprehensive summary PDF...\n")

  pdf(file.path(tier7_dir, "TIER7_meta_analysis_summary.pdf"),
      width = 12, height = 9)

  # Page 1: Sample counts by condition
  if (nrow(sample_summary) > 0) {
    cond_counts <- sample_summary %>%
      group_by(condition) %>%
      summarise(n_samples = n(), total_spots = sum(n_spots), .groups = "drop")
    p_counts <- ggplot(cond_counts, aes(x = condition, y = n_samples, fill = condition)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = paste0("n=", n_samples, "\n", total_spots, " spots")),
                vjust = -0.2, size = 3) +
      labs(title = "Samples and Spots per Condition",
           x = "Condition", y = "Number of Samples") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "none")
    print(p_counts)
  }

  # Page 2: Zone proportions stacked bar
  zone_prop_file <- file.path(visium_root, "TIER1_Interface", "zone_proportions_per_sample.csv")
  if (file.exists(zone_prop_file)) {
    zp <- read.csv(zone_prop_file, stringsAsFactors = FALSE)
    if (nrow(zp) > 0 && "condition" %in% colnames(zp)) {
      zp_agg <- zp %>%
        group_by(condition, zone) %>%
        summarise(count = sum(count), .groups = "drop") %>%
        group_by(condition) %>%
        mutate(prop = count / sum(count))

      zone_colors <- c(
        Interface = "#E31A1C", ASPC_high = "#1F78B4",
        AR_high   = "#33A02C", NEPC_high = "#FF7F00", Other = "#CCCCCC"
      )
      p_zone_bar <- ggplot(zp_agg, aes(x = condition, y = prop, fill = zone)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = zone_colors) +
        labs(title = "Zone Proportions by Condition",
             x = "Condition", y = "Proportion", fill = "Zone") +
        theme_minimal(base_size = 10)
      print(p_zone_bar)
    }
  }

  # Page 3: Score medians per condition (bar plot)
  if (nrow(sample_summary) > 0) {
    score_med_cols_plot <- intersect(
      c("ASPC_score_med","AR_score_med","NEPC_UP_score_med","NEPC_Beltran_NET_med"),
      colnames(sample_summary)
    )
    if (length(score_med_cols_plot) > 0) {
      ss_long <- tidyr::pivot_longer(
        sample_summary,
        cols = all_of(score_med_cols_plot),
        names_to = "score_type",
        values_to = "median_value"
      )
      ss_long <- ss_long[!is.na(ss_long$median_value), ]

      p_score_cond <- ggplot(ss_long, aes(x = condition, y = median_value, fill = condition)) +
        geom_boxplot(outlier.size = 0.5) +
        facet_wrap(~ score_type, scales = "free_y", ncol = 2) +
        labs(title = "Median Signature Scores by Condition (per-sample)",
             x = "Condition", y = "Median Score") +
        theme_minimal(base_size = 10) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
      print(p_score_cond)
    }
  }

  dev.off()
  cat("    Summary PDF saved\n")

  tier7_status <- "COMPLETE"
  cat("TIER 7 complete.\n")

}, error = function(e) {
  cat(paste0("TIER 7 ERROR: ", conditionMessage(e), "\n"))
  tier7_status <<- paste0("ERROR: ", conditionMessage(e))
  tryCatch(dev.off(), error = function(e2) NULL)
})

# ==============================================================================
# SECTION 10: Final Status Report
# ==============================================================================
cat("\n")
cat("================================================================\n")
cat("  GSE278936 Visium 7-Tier Pipeline — FINAL STATUS REPORT\n")
cat("================================================================\n")
cat(paste0("  Samples processed:  ", length(successfully_loaded), "\n"))
cat(paste0("  TIER 1 (Interface): ", tier1_status, "\n"))
cat(paste0("  TIER 2 (DE):        ", tier2_status, "\n"))
cat(paste0("  TIER 3 (Gradient):  ", tier3_status, "\n"))
cat(paste0("  TIER 4 (L-R):       ", tier4_status, "\n"))
cat(paste0("  TIER 5 (Trajectory):", tier5_status, "\n"))
cat(paste0("  TIER 6 (GRN):       ", tier6_status, "\n"))
cat(paste0("  TIER 7 (Meta):      ", tier7_status, "\n"))
cat("================================================================\n")
cat(paste0("  Output root: ", visium_root, "\n"))
cat("================================================================\n")
cat("Pipeline complete.\n")
