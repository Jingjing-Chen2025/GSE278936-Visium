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
            expr_mat <- GetAssayData(obj, layer = "data")
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
                rowMeans(GetAssayData(obj, layer = "data")[common_genes, cells, drop = FALSE])
            })
            vals50 <- lapply(cond_samples, function(sn) {
                obj <- seurat_list[[sn]]
                if (!"dist_bin" %in% colnames(obj@meta.data)) return(NULL)
                cells <- colnames(obj)[obj$dist_bin == "bin_50plus"]
                if (length(cells) == 0) return(NULL)
                DefaultAssay(obj) <- "SCT"
                rowMeans(GetAssayData(obj, layer = "data")[common_genes, cells, drop = FALSE])
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
