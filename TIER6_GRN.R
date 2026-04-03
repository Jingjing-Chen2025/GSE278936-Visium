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
        all_iface_expr_list[[sn]] <- GetAssayData(obj, layer = "data")[, iface_cells, drop = FALSE]
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
