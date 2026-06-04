run_full_nma <- function(
    data,
    outcome = "Mortality",
    sample_size_var = "Overall participants ITT",
    effect_modifiers = c("year_group2", "Inpatient or outapatient", "Publication_bias")
) {
  # ================================================================
  # 0. SETUP
  # ================================================================
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gemtc)
  library(igraph)
  library(ggrepel)       
  library(purrr)
  
  # === Base path (confirmed by user) ===
  base_path <- "/Users/kasimallelhenriquez/Dropbox/B_projects/0_UniversityofOxford/GAPI/0_SR_NetworkMA/0_Article_CAP/0_Analyses"
  
  data_name <- deparse(substitute(data))
  
  # Create output folder
  outdir <- file.path(base_path, data_name, sample_size_var, outcome)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  message("Saving all results to: ", outdir)
  
  # ================================================================
  # 1. CLEAN DATA
  # ================================================================
  dat_arm <- data %>%
    mutate(
      responders = as.numeric(.data[[outcome]]),
      sampleSize = as.numeric(.data[[sample_size_var]])
    ) %>%
    transmute(
      study      = ID_study,
      treatment  = ID_abxDuration,
      responders,
      sampleSize
    ) %>%
    filter(!is.na(responders), !is.na(sampleSize))
  
  
  # ================================================================
  # 2. RUN BAYESIAN NMA USING USER'S FUNCTION
  # ================================================================
  message("Running Bayesian NMA for outcome: ", outcome)
  
  summary_res <- run_bayesian_nma(data, outcome = outcome, sample_size_var = sample_size_var)
  network     <- summary_res$network
  
  # Reference treatment = most frequent
  ref_treat <- names(which.max(table(dat_arm$treatment)))
  message("Reference treatment: ", ref_treat)
  
  
  # ================================================================
  # 3. RELATIVE EFFECTS → FINAL TABLE
  # ================================================================
  re <- gemtc::relative.effect(summary_res$results, t1 = ref_treat)
  draws <- as.matrix(re$samples)
  
  summary_df <- data.frame(
    param = colnames(draws),
    mean  = apply(draws, 2, mean),
    sd    = apply(draws, 2, sd),
    q2.5  = apply(draws, 2, quantile, 0.025),
    q97.5 = apply(draws, 2, quantile, 0.975)
  ) %>%
    filter(grepl("^d\\.", param)) %>%
    mutate(
      param_clean = gsub("^d\\.", "", param),
      Reference   = sub("\\..*$", "", param_clean),
      Comparator  = sub("^.*?\\.", "", param_clean)
    ) %>%
    filter(Reference == ref_treat) %>%
    mutate(
      OR      = exp(mean),
      OR_low  = exp(q2.5),
      OR_high = exp(q97.5)
    )
  
  
  # ================================================================
  # 4. ADD MORTALITY / SAMPLE METADATA
  # ================================================================
  arm_info <- summary_res$arm_data %>%
    mutate(
      deaths = responders,
      total  = sampleSize
    )
  
  mort_tbl <- arm_info %>%
    group_by(treatment) %>%
    summarise(
      total_sample = sum(total),
      total_deaths = sum(deaths),
      .groups = "drop"
    )
  
  zero_tbl <- arm_info %>%
    group_by(treatment) %>%
    summarise(
      zero_event_studies = sum(deaths == 0),
      .groups = "drop"
    )
  
  meta_tbl <- mort_tbl %>%
    left_join(zero_tbl, by = "treatment")
  
  final_results <- summary_df %>%
    left_join(meta_tbl, by = c("Comparator" = "treatment")) %>%
    filter(total_deaths > 0) %>%
    arrange(OR)
  
  
  # ================================================================
  # 5. SUCRA + RANK PROBABILITIES
  # ================================================================
  
  rank_prob_res <- try(rank.probability(summary_res$results), silent = TRUE)
  
  if (inherits(rank_prob_res, "try-error") || is.null(rank_prob_res)) {
    
    warning("Rank probabilities could not be computed. Possibly too few treatments or unstable model.")
    
    rank_df  <- NULL
    sucra_df <- NULL
    
  } else {
    
    # ---- Convert to a clean matrix ----
    rank_mat <- as.matrix(rank_prob_res)
    
    # Strip ALL attributes except dimensions + assign clean names
    attributes(rank_mat) <- list(
      dim = dim(rank_mat),
      dimnames = list(
        rownames(rank_mat),
        paste0("Rank_", seq_len(ncol(rank_mat)))
      )
    )
    
    # ---- SUCRA ----
    sucra_values <- gemtc::sucra(rank_mat)
    
    sucra_df <- data.frame(
      Treatment = names(sucra_values),
      SUCRA = as.numeric(sucra_values),
      row.names = NULL
    )
    
    # ---- LONG FORMAT RANK DATA (pivot_longer-proof) ----
    rank_probs2 <- rank_mat %>%
      as.data.frame(check.names = FALSE) %>%   # KEEP "Rank_1", "Rank_2", ...
      tibble::rownames_to_column("Treatment")
    
    # Sanity check: must have rank columns
    rank_cols <- grep("^Rank_", colnames(rank_probs2), value = TRUE)
    
    if (length(rank_cols) == 0) {
      stop("No rank columns detected after cleaning. Something is wrong with rank_mat.")
    }
    
    rank_df <- rank_probs2 %>%
      pivot_longer(
        cols = all_of(rank_cols),
        names_to = "Rank",
        values_to = "Probability"
      ) %>%
      mutate(
        Rank = as.numeric(sub("Rank_", "", Rank)),
        Treatment = factor(Treatment)
      )
  }
  
  sucra_plot <- NULL
  sucra_uncertainty_df <- NULL
  sucra_uncertainty_plot <- NULL
  posterior_best_plot <- NULL
  posterior_rank_heatmap <- NULL
  
  # ================================================================
  # 5B. SUCRA-BASED UNCERTAINTY PLOT (CONSISTENT SCALE)
  # ================================================================
  
  if (!is.null(rank_df) & !is.null(sucra_df)) {
    
    # Number of treatments
    n_trt <- length(unique(rank_df$Treatment))
    
    # ---- Convert rank distribution → SUCRA-like uncertainty ----
    sucra_uncertainty_df <- rank_df %>%
      arrange(Treatment, Rank) %>%
      group_by(Treatment) %>%
      mutate(
        cum_prob = cumsum(Probability),
        # Convert rank to SUCRA scale (0–1)
        sucra_rank = (n_trt - Rank) / (n_trt - 1)
      ) %>%
      summarise(
        SUCRA_mean = sum(sucra_rank * Probability),
        
        # credible interval using cumulative distribution
        SUCRA_low = sucra_rank[min(which(cum_prob >= 0.975))],  # inverted
        SUCRA_high = sucra_rank[min(which(cum_prob >= 0.025))],
        
        .groups = "drop"
      ) %>%
      
      # 🔴 Apply inversion (your convention)
      mutate(
        SUCRA_mean = 1 - SUCRA_mean,
        SUCRA_low  = 1 - SUCRA_low,
        SUCRA_high = 1 - SUCRA_high
      ) %>%
      
      left_join(
        meta_tbl %>%
          transmute(
            Treatment = treatment,
            total_sample,
            total_deaths
          ),
        by = "Treatment"
      ) %>%
      
      mutate(
        Treatment_label = paste0(
          Treatment,
          " (N=", total_sample,
          ", events=", total_deaths, ")"
        )
      ) %>%
      
      arrange(SUCRA_mean) %>%
      mutate(
        Treatment_label = factor(Treatment_label, levels = Treatment_label)
      )
    
    # ------------------------------------------------
    # Plot (NOW ON SUCRA SCALE)
    # ------------------------------------------------
    sucra_uncertainty_plot <- ggplot(
      sucra_uncertainty_df,
      aes(
        x = SUCRA_mean * 100,
        y = Treatment_label
      )
    ) +
      
      geom_vline(
        xintercept = seq(10, 100, 10),
        linetype = "dashed",
        color = "grey85",
        linewidth = 0.6
      ) +
      
      geom_errorbarh(
        aes(
          xmin = SUCRA_low * 100,
          xmax = SUCRA_high * 100
        ),
        height = 0.2,
        linewidth = 0.9,
        color = "grey50"
      ) +
      
      geom_point(
        size = 5.2,
        shape = 21,
        fill = "#3C5488",
        color = "black",
        stroke = 1.2
      ) +
      
      scale_x_continuous(
        limits = c(0, 100),
        breaks = seq(0, 100, 10),
        expand = expansion(mult = c(0, 0.02))
      ) +
      
      labs(
        title = paste0("SUCRA with Uncertainty – ", outcome),
        subtitle = "Mean SUCRA with 95% credible interval derived from rank distribution",
        x = "SUCRA (%)",
        y = "Treatment (N, events)"
      ) +
      
      theme_classic(base_size = 16) +
      theme(
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.line = element_line(color = "grey60", linewidth = 0.7),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black")
      )
    
    tiff(
      file.path(outdir, paste0("SUCRA_Uncertainty_", outcome, ".tiff")),
      width = 14, height = 11, units = "in", res = 500
    )
    print(sucra_uncertainty_plot)
    dev.off()
    
  } else {
    sucra_uncertainty_plot <- NULL
  }
  
  # ================================================================
  # 6. FOREST PLOT (TIFF SAVED)
  # ================================================================
  forest_df <- final_results %>%
    mutate(
      Comparator_label = paste0(
        Comparator,
        " (N=", total_sample,
        ", events=", total_deaths, ")"
      ),
      Comparator_label = factor(Comparator_label, levels = Comparator_label[order(OR)])
    )
  
  # data-driven breaks
  get_log_breaks <- function(x) {
    x <- x[is.finite(x) & x > 0]
    min_x <- min(x)
    max_x <- max(x)
    possible <- 10^(seq(-10, 10, 1))
    br <- possible[possible >= min_x & possible <= max_x]
    sort(unique(c(br, 1)))
  }
  
  format_or <- function(x) {
    sapply(x, function(val) {
      if (is.na(val)) return(NA)
      if (val < 1) format(val, nsmall = 3, trim = TRUE, scientific = FALSE)
      else format(round(val), scientific = FALSE, trim = TRUE)
    })
  }
  
  forest_plot <- ggplot(forest_df, aes(y = Comparator_label, x = OR)) +
    geom_point(size = 3, color = "#1B4F72") +
    geom_errorbarh(aes(xmin = OR_low, xmax = OR_high), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    scale_x_log10(
      breaks = get_log_breaks(c(forest_df$OR_low, forest_df$OR_high)),
      labels = format_or
    ) +
    labs(
      title = paste0("Network Meta-Analysis – Forest Plot (", outcome, ")"),
      subtitle = paste("Reference:", ref_treat),
      x = "Odds Ratio",
      y = "Comparator (N, events)"
    ) +
    theme_minimal(base_size = 17) + theme(
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"))
  
  tiff(file.path(outdir, paste0("Forest_", outcome, ".tiff")),
       width = 10, height = 12, units = "in", res = 500)
  print(forest_plot)
  dev.off()
  


  # ================================================================
  # 7. SUCRA PLOT (USE RAW SUCRA FOR RANKING)
  # ================================================================
  

  if (!is.null(sucra_df)) {
    
    sucra_plot_df <- sucra_df %>%
      mutate(
        SUCRA_raw = SUCRA,
        SUCRA_raw = 1 - SUCRA_raw   # invert probabilities
      ) %>%
      left_join(
        meta_tbl %>%
          transmute(
            Treatment = treatment,
            total_sample,
            total_deaths
          ),
        by = "Treatment"
      ) %>%
      mutate(
        Treatment_label = paste0(
          Treatment,
          " (N=", total_sample,
          ", events=", total_deaths, ")"
        )
      ) %>%
      arrange(SUCRA_raw) %>%
      mutate(
        Treatment_label = factor(Treatment_label, levels = Treatment_label)
      )
    
    sucra_plot <- ggplot(
      sucra_plot_df,
      aes(
        x = SUCRA_raw * 100,
        y = Treatment_label
      )
    ) +
      
      geom_vline(
        xintercept = seq(10, 100, 10),
        linetype = "dashed",
        color = "grey85",
        linewidth = 0.6
      ) +
      
      geom_segment(
        aes(
          x = 0,
          xend = SUCRA_raw * 100,
          yend = Treatment_label
        ),
        linewidth = 0.8,
        color = "grey70"
      ) +
      
      geom_point(
        size = 5.2,
        shape = 21,
        fill = "#3C5488",
        color = "black",
        stroke = 1.2
      ) +
      
      scale_x_continuous(
        limits = c(0, 100),
        breaks = seq(0, 100, 10),
        expand = expansion(mult = c(0, 0.02))
      ) +
      
      labs(
        title = paste0("SUCRA Ranking – ", outcome),
        subtitle = "SUCRA probabilities: higher values indicate better ranking",
        x = "SUCRA (%)",
        y = "Treatment (N, events)"
      ) +
      
      theme_classic(base_size = 16) +
      theme(
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.line = element_line(color = "grey60", linewidth = 0.7),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black")
      )
    
    tiff(
      file.path(outdir, paste0("SUCRA_", outcome, ".tiff")),
      width = 14, height = 11, units = "in", res = 500
    )
    print(sucra_plot)
    dev.off()
  }
  
  # ================================================================
  # 7.5 NON-INFERIORITY ANALYSIS (10% MARGIN) + PLOT
  # ================================================================
  
  message("Computing Bayesian non-inferiority probabilities (10% margin)")
  
  NI_margin_log <- log(1.10)
  
  draws_df <- as.data.frame(draws) |>
    tibble::rownames_to_column("iter") |>
    pivot_longer(
      cols = -iter,
      names_to = "param",
      values_to = "logOR"
    ) |>
    filter(grepl("^d\\.", param)) |>
    mutate(
      param_clean = gsub("^d\\.", "", param),
      Reference   = sub("\\..*$", "", param_clean),
      Comparator  = sub("^.*?\\.", "", param_clean)
    ) |>
    filter(Reference == ref_treat)
  
  ni_df <- draws_df |>
    group_by(Comparator) |>
    summarise(
      n_draws = n(),
      n_NI    = sum(logOR < NI_margin_log),
      p_non_inferior = n_NI / n_draws,
      ni_low  = qbeta(0.025, n_NI + 0.5, n_draws - n_NI + 0.5),
      ni_high = qbeta(0.975, n_NI + 0.5, n_draws - n_NI + 0.5),
      .groups = "drop"
    ) |>
    left_join(
      final_results |> select(Comparator, total_sample, total_deaths),
      by = "Comparator"
    ) |>
    arrange(p_non_inferior) |>
    mutate(
      Comparator_label = paste0(
        Comparator,
        " (N=", total_sample,
        ", events=", total_deaths, ")"
      ),
      Comparator_label = factor(
        Comparator_label,
        levels = Comparator_label
      )
    )
  
  

  
  # ================================================================
  # 7B. POSTERIOR PROBABILITY OF BEING BEST / RANK 1
  # For mortality: lower effect = better, so use Rank == n_trt
  # ================================================================
  
  posterior_best_df <- NULL
  posterior_best_plot <- NULL
  
  if (!is.null(rank_df)) {
    
    n_trt <- length(unique(rank_df$Treatment))
    
    posterior_best_df <- rank_df %>%
      filter(Rank == n_trt) %>%   # best for lower mortality after your inversion convention
      rename(Prob_rank1_best = Probability) %>%
      left_join(
        meta_tbl %>%
          transmute(
            Treatment = treatment,
            total_sample,
            total_deaths
          ),
        by = "Treatment"
      ) %>%
      mutate(
        Prob_rank1_best_percent = Prob_rank1_best * 100,
        Treatment_label = paste0(
          Treatment,
          " (N=", total_sample,
          ", events=", total_deaths, ")"
        )
      ) %>%
      arrange(desc(Prob_rank1_best))
    
    # Sanity check: should be 1
    prob_sum <- sum(posterior_best_df$Prob_rank1_best, na.rm = TRUE)
    
    message("Sum of posterior Prob(rank 1 / best) across treatments = ",
            round(prob_sum, 6))
    
    if (abs(prob_sum - 1) > 0.01) {
      warning("Posterior best probabilities do not sum to 1. Check rank orientation or rank.probability output.")
    }
    
    # Save Excel output
    writexl::write_xlsx(
      posterior_best_df,
      file.path(outdir, paste0("Posterior_Probability_Rank1_Best_", outcome, ".xlsx"))
    )
    
    # Plot
    posterior_best_plot_df <- posterior_best_df %>%
      arrange(Prob_rank1_best) %>%
      mutate(
        Treatment_label = factor(Treatment_label, levels = Treatment_label)
      )
    
    posterior_best_plot <- ggplot(
      posterior_best_plot_df,
      aes(
        x = Prob_rank1_best_percent,
        y = Treatment_label
      )
    ) +
      geom_vline(
        xintercept = seq(10, 100, 10),
        linetype = "dashed",
        color = "grey85",
        linewidth = 0.6
      ) +
      geom_segment(
        aes(
          x = 0,
          xend = Prob_rank1_best_percent,
          yend = Treatment_label
        ),
        linewidth = 0.8,
        color = "grey70"
      ) +
      geom_point(
        size = 5.2,
        shape = 21,
        fill = "#3C5488",
        color = "black",
        stroke = 1.2
      ) +
      scale_x_continuous(
        limits = c(0, max(posterior_best_plot_df$Prob_rank1_best_percent, na.rm = TRUE) * 1.10),
        expand = expansion(mult = c(0, 0.02))
      ) +
      labs(
        title = paste0("Posterior Probability of Ranking First – ", outcome),
        subtitle = paste0("Probabilities sum to ", round(prob_sum, 3), " across treatments"),
        x = "Posterior probability of rank 1 / best (%)",
        y = "Treatment (N, events)"
      ) +
      theme_classic(base_size = 16) +
      theme(
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.line = element_line(color = "grey60", linewidth = 0.7),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black")
      )
    
    tiff(
      file.path(outdir, paste0("Posterior_Probability_Rank1_Best_", outcome, ".tiff")),
      width = 14, height = 11, units = "in", res = 500
    )
    print(posterior_best_plot)
    dev.off()
  }
  

  
#. --- -- -- -- -- -- -- -- -- --  -- -- -- -- -- ##
  # 7C. EXCEL SUMMARY:
  # Prob(rank 1 best), closeness to best, non-inferiority vs priority ref,
  # SUCRA, total sample, total deaths
  #. --- -- -- -- -- -- -- -- -- --  -- -- -- -- -- ##
  
  message("Creating antibiotic-duration summary Excel file")
  
  delta <- 0.10
  NI_margin_log <- log(1 + delta)
  

  # Convert posterior relative-effect draws to wide matrix
  # Reference treatment has logOR = 0 by construction
  # Lower logOR = better for mortality

  effect_wide <- draws_df %>%
    select(iter, Comparator, logOR) %>%
    pivot_wider(
      names_from = Comparator,
      values_from = logOR
    ) %>%
    mutate(!!ref_treat := 0)
  
  treatment_cols <- setdiff(names(effect_wide), "iter")
  
  effect_mat <- effect_wide %>%
    select(all_of(treatment_cols)) %>%
    as.matrix()
  


  # ------------------------------------------------
  # A. Closeness to best based on Prob_rank1_best
  # Best treatment always gets 1
  # Close if Prob_rank1_best >= best_prob * (1 - delta)

  best_prob_rank1 <- max(posterior_best_df$Prob_rank1_best, na.rm = TRUE)
  
  best_prob_rank1_treatment <- posterior_best_df %>%
    filter(Prob_rank1_best == best_prob_rank1) %>%
    slice(1) %>%
    pull(Treatment)
  
  close_to_best_df <- posterior_best_df %>%
    transmute(
      Treatment,
      Prob_close_to_best_delta10 =
        pmin(Prob_rank1_best / best_prob_rank1, 1),
      Best_prob_rank1_reference = best_prob_rank1_treatment,
      Best_prob_rank1_threshold_delta10 = best_prob_rank1 * (1 - delta)
    )
  

  # B. Select reference antibiotic-duration using priority order

  reference_priority <- c(
    "Levofloxacin_S", "Levofloxacin_L",
    "Clarithromycin_S", "Clarithromycin_L",
    "Erythromycin_S", "Erythromycin_L",
    "Moxifloxacin_S", "Moxifloxacin_L"
  )
  
  ni_reference <- reference_priority[
    reference_priority %in% treatment_cols
  ][1]
  
  if (is.na(ni_reference)) {
    warning("None of the preferred NI reference treatments were found.")
    ni_reference <- treatment_cols[which.max(
      posterior_best_df$Prob_rank1_best[
        match(treatment_cols, posterior_best_df$Treatment)
      ]
    )]
  }
  
  message("Non-inferiority reference used: ", ni_reference)
  

  # C. Probability of being non-inferior to selected reference
  # For mortality: treatment is non-inferior if
  # logOR_treatment - logOR_reference <= log(1.10)

  ref_draws <- effect_mat[, ni_reference]
  
  prob_NI_vs_ref <- sapply(treatment_cols, function(trt) {
    if (trt == ni_reference) {
      return(NA_real_)
    }
    
    mean(
      effect_mat[, trt] - ref_draws <= NI_margin_log,
      na.rm = TRUE
    )
  })
  
  ni_vs_ref_df <- data.frame(
    Treatment = names(prob_NI_vs_ref),
    Prob_noninferior_vs_priority_ref_delta10 = as.numeric(prob_NI_vs_ref),
    row.names = NULL
  ) %>%
    mutate(
      Prob_noninferior_vs_priority_ref_delta10 =
        ifelse(
          Treatment == ni_reference,
          "REF",
          as.character(round(Prob_noninferior_vs_priority_ref_delta10, 6))
        )
    )
  

  # D. Build final Excel table

  antibiotic_duration_summary <- posterior_best_df %>%
    transmute(
      Antibiotic_duration = Treatment,
      Prob_rank1_best
    ) %>%
    left_join(
      close_to_best_df,
      by = c("Antibiotic_duration" = "Treatment")
    ) %>%
    left_join(
      ni_vs_ref_df,
      by = c("Antibiotic_duration" = "Treatment")
    ) %>%
    left_join(
      sucra_uncertainty_df %>%
        transmute(
          Antibiotic_duration = as.character(Treatment),
          SUCRA_mean
        ),
      by = "Antibiotic_duration"
    ) %>%
    left_join(
      meta_tbl %>%
        transmute(
          Antibiotic_duration = treatment,
          total_sample,
          total_deaths
        ),
      by = "Antibiotic_duration"
    ) %>%
    arrange(desc(Prob_rank1_best))
  
  writexl::write_xlsx(
    antibiotic_duration_summary,
    file.path(
      outdir,
      paste0("Antibiotic_Duration_Summary_", outcome, "_delta10.xlsx")
    )
  )
  
  
  
  
  # ================================================================
  # NON-INFERIORITY LOLLIPOP PLOT (Refined Lancet Style)
  # ================================================================
  
  ni_plot <- ggplot(
    ni_df,
    aes(
      x = p_non_inferior * 100,
      y = Comparator_label
    )
  ) +
    
    # Light grey dashed vertical grid lines (every 10%)
    geom_vline(
      xintercept = seq(10, 100, 10),
      linetype = "dashed",
      color = "grey85",
      linewidth = 0.6
    ) +
    
    # Lollipop stems
    geom_segment(
      aes(
        x = 0,
        xend = p_non_inferior * 100,
        yend = Comparator_label
      ),
      linewidth = 0.8,
      color = "grey70"
    ) +
    geom_point(
      size = 5.2,
      shape = 21,
      fill = "#3C5488",
      color = "black",
      stroke = 1.2
    ) + scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 10),
      expand = expansion(mult = c(0, 0.02))  # small space on right
    ) +labs(
      title = paste0("Bayesian Non-Inferiority – ", outcome),
      subtitle = paste0(
        "Reference: ",
        ref_treat
      ),
      x = "Posterior probability of non-inferiority using a 10% margin (%)",
      y = "Comparator (N, events)"
    ) +
    
    theme_classic(base_size = 16) +
    theme(
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      
      # Softer axis lines
      axis.line = element_line(color = "grey60", linewidth = 0.7),
      
      # Remove horizontal grid lines
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      
      axis.text = element_text(color = "black")
    )
  
  # Wider export
  tiff(
    file.path(outdir, paste0("NonInferiority_", outcome, ".tiff")),
    width = 14, height = 11, units = "in", res = 500
  )
  print(ni_plot)
  dev.off()
  
  write.csv(
    ni_df,
    file = file.path(outdir, paste0("NonInferiority_", outcome, "_data.csv")),
    row.names = FALSE
  )
  
  
  
  
  # ================================================================
  # 8. NETWORK PLOT (TIFF SAVED)
  # ================================================================
  
  arm_list <- network$data.ab %>%
    distinct(study, treatment) %>%
    group_by(study) %>%
    filter(n_distinct(treatment) >= 2) %>%
    ungroup()
  
  direct_edges <- arm_list %>%
    group_by(study) %>%
    summarise(
      pairs = list(combn(as.character(unique(treatment)), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(pairs) %>%
    transmute(
      source = purrr::map_chr(pairs, 1),
      target = purrr::map_chr(pairs, 2),
      direct = TRUE
    )
  
  treats <- sort(unique(arm_list$treatment))
  
  all_pairs <- expand.grid(
    source = treats,
    target = treats,
    stringsAsFactors = FALSE
  ) %>%
    filter(source < target)
  
  indirect_edges <- all_pairs %>%
    anti_join(direct_edges, by = c("source","target")) %>%
    mutate(direct = FALSE)
  
  edges <- bind_rows(direct_edges, indirect_edges)
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  coords <- layout_with_fr(g)
  rownames(coords) <- V(g)$name
  
  coords_df <- as.data.frame(coords)
  coords_df$Treatment <- V(g)$name
  
  edge_df <- edges %>%
    mutate(
      x    = coords[source, 1],
      y    = coords[source, 2],
      xend = coords[target, 1],
      yend = coords[target, 2]
    )
  
  network_plot <- ggplot() +
    geom_segment(
      data = edge_df %>% filter(!direct),
      aes(x = x, y = y, xend = xend, yend = yend),
      colour = "grey85", linetype = "dashed", linewidth = 0.3
    ) +
    geom_segment(
      data = edge_df %>% filter(direct),
      aes(x = x, y = y, xend = xend, yend = yend),
      colour = "grey85", linetype = "solid", linewidth = 0.6
    ) +
    geom_point(
      data = coords_df,
      aes(x = V1, y = V2),
      size = 4,
      colour = "#1B4F72"
    ) +
    geom_text_repel(
      data = coords_df,
      aes(x = V1, y = V2, label = Treatment),
      size = 3
    ) +
    theme_void() +
    ggtitle(paste0("NMA Network – ", outcome))
  
  tiff(file.path(outdir, paste0("Network_", outcome, ".tiff")),
       width = 10, height = 10, units = "in", res = 500)
  print(network_plot)
  dev.off()
  
  ##NEW ANALYSIS.
  
  
  
  # ================================================================
  # DIRECT COMPARISON COUNTS (CORRECT OBJECT)
  # ================================================================
  
  arm_data_used <- summary_res$arm_data
  
  direct_edges <- arm_data_used %>%
    group_by(study) %>%
    filter(n_distinct(treatment) >= 2) %>%
    summarise(
      pairs = list(combn(as.character(treatment), 2, simplify = FALSE)),
      .groups = "drop"
    ) %>%
    unnest(pairs) %>%
    transmute(
      source = map_chr(pairs, 1),
      target = map_chr(pairs, 2)
    )
  
  direct_counts <- direct_edges %>%
    count(source, target, name = "n_direct")
  
  writexl::write_xlsx(
    direct_counts,
    file.path(outdir, "Direct_Comparison_Counts.xlsx")
  )
  
    
  # ================================================================
  # HETEROGENEITY (ROBUST: VARIABILITY OF BASIC PARAMETERS)
  # ================================================================
  
  samples_mat <- summary_res$results$samples
  
  # basic treatment effects are named like d[treatment]
  d_cols <- grep("^d\\[", colnames(samples_mat), value = TRUE)
  
  if (length(d_cols) > 1) {
    
    # SD of treatment effects per iteration
    hetero_draws <- apply(
      samples_mat[, d_cols, drop = FALSE],
      1,
      sd
    )
    
    heterogeneity_tbl <- data.frame(
      heterogeneity_sd = mean(hetero_draws),
      q2.5 = quantile(hetero_draws, 0.025),
      q97.5 = quantile(hetero_draws, 0.975)
    )
    
  } else {
    
    heterogeneity_tbl <- data.frame(
      heterogeneity_sd = NA_real_,
      q2.5 = NA_real_,
      q97.5 = NA_real_
    )
  }
  
  write.csv(
    heterogeneity_tbl,
    file.path(outdir, "Heterogeneity_treatment_effect_variability.csv"),
    row.names = FALSE
  )
  
  
  # ================================================================
  # 8.5 ASSUMPTION / DIAGNOSTIC CHECKS
  # Homogeneity, similarity/transitivity, consistency, exchangeability
  # ================================================================
  
  suppressPackageStartupMessages({
    library(coda)
    library(writexl)
  })
  
  diagnostics <- list()
  
  # ------------------------------------------------
  # A. MCMC convergence / exchangeability diagnostics
  # ------------------------------------------------
  message("Running MCMC convergence diagnostics")
  
  mcmc_obj <- try(as.mcmc.list(summary_res$results), silent = TRUE)
  
  if (!inherits(mcmc_obj, "try-error")) {
    
    gelman_res <- try(coda::gelman.diag(mcmc_obj, multivariate = FALSE), silent = TRUE)
    geweke_res <- try(coda::geweke.diag(mcmc_obj), silent = TRUE)
    
    if (!inherits(gelman_res, "try-error")) {
      gelman_tbl <- as.data.frame(gelman_res$psrf)
      gelman_tbl$parameter <- rownames(gelman_tbl)
      rownames(gelman_tbl) <- NULL
      
      write.csv(
        gelman_tbl,
        file.path(outdir, "Diagnostics_MCMC_Gelman_Rhat.csv"),
        row.names = FALSE
      )
      
      diagnostics$gelman <- gelman_tbl
    }
    
    pdf(file.path(outdir, "Diagnostics_MCMC_traceplots.pdf"),
        width = 12, height = 8)
    plot(mcmc_obj)
    dev.off()
    
    diagnostics$mcmc <- mcmc_obj
  }
  
  # ------------------------------------------------
  # B. Model fit / residual deviance / DIC
  # ------------------------------------------------
  message("Extracting model fit diagnostics")
  
  model_summary <- try(summary(summary_res$results), silent = TRUE)
  
  if (!inherits(model_summary, "try-error")) {
    capture.output(
      model_summary,
      file = file.path(outdir, "Diagnostics_Model_Summary.txt")
    )
    diagnostics$model_summary <- model_summary
  }
  
  # ------------------------------------------------
  # C. Between-study heterogeneity from NMA
  # ------------------------------------------------
  message("Extracting posterior heterogeneity parameter")
  
  samples_mat <- as.matrix(summary_res$results)
  
  tau_cols <- grep("sd|tau", colnames(samples_mat), value = TRUE, ignore.case = TRUE)
  
  if (length(tau_cols) > 0) {
    
    tau_tbl <- lapply(tau_cols, function(x) {
      vals <- samples_mat[, x]
      data.frame(
        parameter = x,
        mean = mean(vals),
        median = median(vals),
        q2.5 = quantile(vals, 0.025),
        q97.5 = quantile(vals, 0.975)
      )
    }) |>
      bind_rows()
    
    write.csv(
      tau_tbl,
      file.path(outdir, "Diagnostics_NMA_Heterogeneity_Tau.csv"),
      row.names = FALSE
    )
    
    diagnostics$tau <- tau_tbl
  }
  
  # ------------------------------------------------
  # D. Direct-comparison homogeneity
  # ------------------------------------------------
  message("Running direct-comparison heterogeneity checks")
  
  dat_diag <- data %>%
    mutate(
      event = as.numeric(.data[[outcome]]),
      n = as.numeric(.data[[sample_size_var]])
    ) %>%
    filter(!is.na(event), !is.na(n)) %>%
    distinct(ID_study, ID_abxDuration, .keep_all = TRUE)
  
  pw_diag <- meta::pairwise(
    treat = ID_abxDuration,
    event = event,
    n = n,
    studlab = ID_study,
    data = dat_diag,
    sm = "OR",
    add = 1,
    allstudies = TRUE
  ) %>%
    filter(!is.na(TE), !is.na(seTE), seTE > 0) %>%
    mutate(
      comparison = paste(pmin(treat1, treat2), pmax(treat1, treat2), sep = " vs ")
    )
  
  direct_heterogeneity <- pw_diag %>%
    group_split(comparison) %>%
    map_dfr(function(df) {
      
      if (nrow(df) < 2) {
        return(data.frame(
          comparison = unique(df$comparison),
          k = nrow(df),
          tau2 = NA_real_,
          I2 = NA_real_,
          Q = NA_real_,
          Q_p = NA_real_
        ))
      }
      
      m <- meta::metagen(
        TE = TE,
        seTE = seTE,
        studlab = studlab,
        data = df,
        sm = "OR",
        common = FALSE,
        random = TRUE,
        method.tau = "REML"
      )
      
      data.frame(
        comparison = unique(df$comparison),
        k = m$k,
        tau2 = m$tau2,
        I2 = m$I2,
        Q = m$Q,
        Q_p = m$pval.Q
      )
    })
  
  write.csv(
    direct_heterogeneity,
    file.path(outdir, "Diagnostics_Direct_Homogeneity.csv"),
    row.names = FALSE
  )
  
  diagnostics$direct_heterogeneity <- direct_heterogeneity
  
  # ------------------------------------------------
  # E. Similarity / transitivity using effect modifiers
  # ------------------------------------------------
  message("Checking similarity/transitivity across treatments and comparisons")
  
  available_modifiers <- effect_modifiers[effect_modifiers %in% names(data)]
  
  if (length(available_modifiers) > 0) {
    
    transitivity_treatment <- dat_diag %>%
      transmute(
        study = as.character(ID_study),
        treatment = as.character(ID_abxDuration),
        across(all_of(available_modifiers), as.character)
      ) %>%
      distinct()
    
    transitivity_by_treatment <- transitivity_treatment %>%
      pivot_longer(
        cols = all_of(available_modifiers),
        names_to = "effect_modifier",
        values_to = "value"
      ) %>%
      group_by(effect_modifier, treatment, value) %>%
      summarise(n_arms = n(), .groups = "drop")
    
    
    study_modifiers <- dat_diag %>%
      transmute(
        study = as.character(ID_study),
        across(all_of(available_modifiers), as.character)
      ) %>%
      distinct()
    
    comparison_modifiers <- pw_diag %>%
      transmute(
        study = as.character(studlab),
        comparison = comparison
      ) %>%
      distinct() %>%
      left_join(study_modifiers, by = "study") %>%
      pivot_longer(
        cols = all_of(available_modifiers),
        names_to = "effect_modifier",
        values_to = "value"
      ) %>%
      group_by(comparison, effect_modifier, value) %>%
      summarise(n_studies = n_distinct(study), .groups = "drop")
    
    write.csv(
      comparison_modifiers,
      file.path(outdir, "Diagnostics_Transitivity_By_Comparison.csv"),
      row.names = FALSE
    )
    
    diagnostics$transitivity_by_treatment <- transitivity_by_treatment
    diagnostics$transitivity_by_comparison <- comparison_modifiers
  }
  
  # ------------------------------------------------
  # F. Consistency diagnostics: node-splitting
  # ------------------------------------------------
  message("Running node-splitting consistency diagnostics")
  
  nodesplit_res <- try(
    gemtc::mtc.nodesplit(
      network,
      linearModel = "random",
      likelihood = "binom",
      link = "logit",
      n.adapt = 5000,
      n.iter = 30000
    ),
    silent = TRUE
  )
  
  if (!inherits(nodesplit_res, "try-error")) {
    
    nodesplit_summary <- summary(nodesplit_res)
    
    capture.output(
      nodesplit_summary,
      file = file.path(outdir, "Diagnostics_NodeSplitting_Consistency.txt")
    )
    
    diagnostics$nodesplit <- nodesplit_summary
    
  } else {
    warning("Node-splitting failed. This can happen when there are too few closed loops or sparse direct evidence.")
    diagnostics$nodesplit <- NULL
  }
  
  # ------------------------------------------------
  # G. Consistency vs inconsistency model comparison
  # ------------------------------------------------
  message("Comparing consistency and inconsistency models")
  
  incons_model <- try(
    mtc.model(
      network,
      type = "inconsistency",
      linearModel = "random",
      likelihood = "binom",
      link = "logit",
      n.chain = 4
    ),
    silent = TRUE
  )
  
  if (!inherits(incons_model, "try-error")) {
    
    set.seed(1234)
    
    incons_res <- try(
      mtc.run(incons_model, n.adapt = 5000, n.iter = 30000),
      silent = TRUE
    )
    
    if (!inherits(incons_res, "try-error")) {
      
      consistency_summary <- try(summary(summary_res$results), silent = TRUE)
      inconsistency_summary <- try(summary(incons_res), silent = TRUE)
      
      capture.output(
        list(
          consistency_model = consistency_summary,
          inconsistency_model = inconsistency_summary
        ),
        file = file.path(outdir, "Diagnostics_Consistency_vs_Inconsistency_Model.txt")
      )
      
      diagnostics$inconsistency_model <- incons_res
    }
  }
  
  # ------------------------------------------------
  # H. Export combined diagnostics object
  # ------------------------------------------------
  saveRDS(
    diagnostics,
    file.path(outdir, "Diagnostics_All_Assumptions.rds")
  )
  
  
  #Export relative effects
  # ================================================================
  # EXPORT RELATIVE EFFECTS VS REFERENCE (EXCEL)
  # ================================================================
  
  relative_effects_xlsx <- summary_df %>%
    select(
      Reference,
      Comparator,
      mean,
      q2.5,
      q97.5,
      OR,
      OR_low,
      OR_high
    ) %>%
    arrange(OR)
  
  writexl::write_xlsx(
    relative_effects_xlsx,
    file.path(outdir, "Relative_Effects_vs_Reference.xlsx")
  )
  
  if (!is.null(sucra_uncertainty_df)) {
    writexl::write_xlsx(
      sucra_uncertainty_df,
      file.path(outdir, paste0("SUCRA_Uncertainty_", outcome, ".xlsx"))
    )
  }
  
  
  # ================================================================
  # 9. RETURN OBJECT
  # ================================================================
  return(list(
    summary_res   = summary_res,
    final_results = final_results,
    sucra_df      = sucra_df,
    rank_df       = rank_df,
    ni_df         = ni_df,
    diagnostics   = diagnostics,
    forest_plot   = forest_plot,
    sucra_plot    = sucra_plot,
    ni_plot       = ni_plot,
    network_plot  = network_plot,
    ref_treat     = ref_treat,
    posterior_best_plot = posterior_best_plot,
    sucra_uncertainty_df   = sucra_uncertainty_df,
    sucra_uncertainty_plot = sucra_uncertainty_plot,
    posterior_best_df = posterior_best_df,
    posterior_best_plot = posterior_best_plot,
    antibiotic_duration_summary = antibiotic_duration_summary,
    ni_reference = ni_reference
  ))
}



  
  