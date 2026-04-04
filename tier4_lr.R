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
        expr_mat <- GetAssayData(obj, layer = "data")

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
            top30$pair <- paste0(top30$ligand, "\u2192", top30$receptor)
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
