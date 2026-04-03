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
            expr <- GetAssayData(obj, layer = "data")[genes_avail, , drop = FALSE]
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
