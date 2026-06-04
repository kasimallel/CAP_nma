run_full_bayesian_nma2 <- function(
    data,
    outcome,
    sample_size_var = "Overall participants ITT",
    study_var = "ID_study",
    treat_var = "ID_abxDuration",
    group_vars = NULL,
    effect_modifiers = NULL,
    route_path = getwd(),
    file_prefix = "Bayesian_NMA",
    reference = NULL,
    lower_is_better = TRUE,
    noninf_margin_relative = 0.10,
    n_adapt = 5000,
    n_iter = 30000,
    n_chain = 4,
    seed = 1234,
    absolute_risk_method = c("gemtc", "empirical_beta")
) {
  
  absolute_risk_method <- match.arg(absolute_risk_method)
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(meta)
    library(netmeta)
    library(gemtc)
    library(igraph)
    library(coda)
    library(ggplot2)
    library(openxlsx)
    library(purrr)
    library(stringr)
    library(tibble)
  })
  
  dir.create(route_path, recursive = TRUE, showWarnings = FALSE)
  
  q025 <- function(x) stats::quantile(x, 0.025, na.rm = TRUE)
  q975 <- function(x) stats::quantile(x, 0.975, na.rm = TRUE)
  
  safe_sheet <- function(x) {
    substr(gsub("[\\[\\]\\*\\?/\\\\:]", "_", x), 1, 31)
  }
  
  get_abs_risk_draws <- function(base_res, arm_data, treatments, n_draws = 10000) {
    
    if (absolute_risk_method == "gemtc") {
      if (!exists("absolute.effect", where = asNamespace("gemtc"), inherits = FALSE)) {
        stop(
          "gemtc::absolute.effect() is not available in this installation. ",
          "Use absolute_risk_method = 'empirical_beta', or estimate absolute risks externally."
        )
      }
      
      out <- list()
      
      for (trt in treatments) {
        draws <- tryCatch({
          as.numeric(as.matrix(gemtc::absolute.effect(base_res, t1 = trt)))
        }, error = function(e) {
          stop(
            "Could not extract model-based absolute risk for treatment ", trt, ". ",
            "Absolute-risk NI requires posterior absolute risk draws."
          )
        })
        
        out[[trt]] <- draws
      }
      
      return(out)
    }
    
    if (absolute_risk_method == "empirical_beta") {
      warning(
        "Using empirical_beta absolute risks. These are Bayesian arm-level posterior risks, ",
        "not network-model-implied absolute risks. Use only as sensitivity analysis."
      )
      
      out <- list()
      
      for (trt in treatments) {
        tmp <- arm_data %>%
          filter(treatment == trt) %>%
          summarise(
            events = sum(responders, na.rm = TRUE),
            total = sum(sampleSize, na.rm = TRUE),
            .groups = "drop"
          )
        
        out[[trt]] <- stats::rbeta(
          n_draws,
          shape1 = tmp$events + 0.5,
          shape2 = tmp$total - tmp$events + 0.5
        )
      }
      
      return(out)
    }
  }
  
  build_sucra <- function(ranks) {
    if ("note" %in% names(ranks)) return(ranks)
    
    rank_cols <- grep("^rank_", names(ranks), value = TRUE)
    
    rank_mat <- as.matrix(ranks[, rank_cols, drop = FALSE])
    storage.mode(rank_mat) <- "numeric"
    
    k <- ncol(rank_mat)
    
    cumulative <- t(apply(rank_mat, 1, cumsum))
    colnames(cumulative) <- paste0("cum_rank_", seq_len(k))
    
    cumulative_tbl <- as_tibble(cumulative)
    cumulative_tbl$treatment <- ranks$treatment
    cumulative_tbl <- cumulative_tbl %>%
      select(treatment, everything())
    
    sucra_tbl <- tibble(
      treatment = ranks$treatment,
      mean_rank = as.numeric(rank_mat %*% seq_len(k)),
      prob_best = rank_mat[, 1],
      prob_worst = rank_mat[, k],
      SUCRA = apply(rank_mat, 1, function(p) {
        sum(cumsum(p)[1:(k - 1)]) / (k - 1)
      })
    )
    
    list(
      cumulative_tbl = cumulative_tbl,
      sucra_tbl = sucra_tbl
    )
  }
  
  analyse_one <- function(dat0, label = "overall") {
    
    message("Running Bayesian NMA for: ", label)
    
    dat <- dat0 %>%
      mutate(
        outcome_val = as.numeric(.data[[outcome]]),
        N = as.numeric(.data[[sample_size_var]]),
        study = as.character(.data[[study_var]]),
        treatment = as.character(.data[[treat_var]])
      ) %>%
      filter(!is.na(outcome_val), !is.na(N), !is.na(study), !is.na(treatment)) %>%
      mutate(
        outcome_val = as.integer(outcome_val),
        N = as.integer(N)
      ) %>%
      filter(outcome_val >= 0, N > 0, outcome_val <= N) %>%
      distinct(study, treatment, .keep_all = TRUE)
    
    arm_audit <- dat %>%
      count(study, name = "n_arms") %>%
      mutate(flag = ifelse(n_arms < 2, "excluded_single_arm", "ok"))
    
    dat <- dat %>%
      semi_join(filter(arm_audit, n_arms >= 2), by = "study")
    
    if (nrow(dat) == 0) stop("No valid multi-arm study data for ", label)
    
    pw <- meta::pairwise(
      treat = treatment,
      event = outcome_val,
      n = N,
      studlab = study,
      data = dat,
      sm = "OR",
      method.incr = "only0",
      incr = 0.5,
      allstudies = TRUE
    )
    
    pw_clean <- as.data.frame(pw) %>%
      filter(!is.na(TE), !is.na(seTE), seTE > 0)
    
    if (nrow(pw_clean) == 0) stop("No valid pairwise contrasts for ", label)
    
    conn <- netmeta::netconnection(
      treat1 = pw_clean$treat1,
      treat2 = pw_clean$treat2,
      studlab = pw_clean$studlab
    )
    
    subnet_df <- data.frame(
      studlab = pw_clean$studlab,
      subnet = conn$subnet
    )
    
    largest_subnet <- subnet_df %>%
      count(subnet, sort = TRUE) %>%
      slice(1) %>%
      pull(subnet)
    
    keep_studies <- subnet_df %>%
      filter(subnet == largest_subnet) %>%
      pull(studlab)
    
    pw_conn <- pw_clean %>% filter(studlab %in% keep_studies)
    dat_conn <- dat %>% filter(study %in% keep_studies)
    
    network_edges <- pw_conn %>%
      count(treat1, treat2, name = "n_direct_studies")
    
    g <- igraph::graph_from_data_frame(
      network_edges[, c("treat1", "treat2")],
      directed = FALSE
    )
    
    connectivity <- tibble(
      connected = igraph::is_connected(g),
      n_treatments = igraph::vcount(g),
      n_direct_edges = igraph::ecount(g),
      n_studies = n_distinct(dat_conn$study),
      n_arms = nrow(dat_conn)
    )
    
    if (!connectivity$connected) {
      stop("Network remains disconnected for ", label)
    }
    
    ref_use <- ifelse(
      is.null(reference),
      sort(unique(dat_conn$treatment))[1],
      reference
    )
    
    nm <- netmeta::netmeta(
      TE = TE,
      seTE = seTE,
      treat1 = treat1,
      treat2 = treat2,
      studlab = studlab,
      data = pw_conn,
      sm = "OR",
      random = TRUE,
      common = TRUE,
      reference.group = ref_use
    )
    
    heterogeneity_tbl <- tibble(
      tau2_netmeta = nm$tau^2,
      tau_netmeta = nm$tau,
      I2_netmeta = nm$I2,
      Q_total = nm$Q,
      Q_df = nm$df.Q,
      Q_p = nm$pval.Q
    )
    
    netsplit_tbl <- tryCatch({
      ns <- netmeta::netsplit(nm)
      as.data.frame(ns$random) %>%
        tibble::rownames_to_column("comparison")
    }, error = function(e) {
      tibble(note = paste("netsplit failed:", e$message))
    })
    
    arm_data <- dat_conn %>%
      transmute(
        study = study,
        treatment = treatment,
        responders = outcome_val,
        sampleSize = N
      )
    
    network <- gemtc::mtc.network(data.ab = arm_data)
    
    model_re <- gemtc::mtc.model(
      network,
      type = "consistency",
      linearModel = "random",
      likelihood = "binom",
      link = "logit",
      n.chain = n_chain,
      dic = TRUE
    )
    
    model_fe <- gemtc::mtc.model(
      network,
      type = "consistency",
      linearModel = "fixed",
      likelihood = "binom",
      link = "logit",
      n.chain = n_chain,
      dic = TRUE
    )
    
    set.seed(seed)
    res_re <- gemtc::mtc.run(model_re, n.adapt = n_adapt, n.iter = n_iter)
    
    set.seed(seed)
    res_fe <- gemtc::mtc.run(model_fe, n.adapt = n_adapt, n.iter = n_iter)
    
    dic_re <- res_re$deviance
    dic_fe <- res_fe$deviance
    
    model_fit <- tibble(
      model = c("fixed_consistency", "random_consistency"),
      DIC = c(dic_fe$DIC, dic_re$DIC),
      pD = c(dic_fe$pD, dic_re$pD),
      residual_deviance = c(dic_fe$Dbar, dic_re$Dbar),
      n_data_points = nrow(arm_data)
    )
    
    if (is.finite(dic_re$DIC) && dic_re$DIC + 3 < dic_fe$DIC) {
      base_res <- res_re
      selected_model <- "random_consistency"
      selected_linear_model <- "random"
    } else {
      base_res <- res_fe
      selected_model <- "fixed_consistency"
      selected_linear_model <- "fixed"
    }
    
    convergence <- tryCatch({
      gd <- coda::gelman.diag(base_res, multivariate = FALSE)
      as.data.frame(gd$psrf) %>%
        tibble::rownames_to_column("parameter") %>%
        rename(
          PSRF = `Point est.`,
          upper_CI = `Upper C.I.`
        )
    }, error = function(e) {
      tibble(note = paste("Gelman diagnostic failed:", e$message))
    })
    
    treatments <- sort(unique(arm_data$treatment))
    comps <- t(combn(treatments, 2))
    
    relative_effects_logscale <- map_dfr(seq_len(nrow(comps)), function(i) {
      a <- comps[i, 1]
      b <- comps[i, 2]
      
      x <- tryCatch({
        as.numeric(as.matrix(gemtc::relative.effect(base_res, t1 = a, t2 = b)))
      }, error = function(e) {
        NA_real_
      })
      
      if (all(is.na(x))) {
        return(tibble(
          treat1 = a,
          treat2 = b,
          log_effect_median = NA_real_,
          log_effect_lCrI = NA_real_,
          log_effect_uCrI = NA_real_,
          prob_log_effect_less_0 = NA_real_,
          prob_log_effect_greater_0 = NA_real_,
          note = "relative.effect failed"
        ))
      }
      
      tibble(
        treat1 = a,
        treat2 = b,
        log_effect_median = median(x, na.rm = TRUE),
        log_effect_lCrI = q025(x),
        log_effect_uCrI = q975(x),
        prob_log_effect_less_0 = mean(x < 0, na.rm = TRUE),
        prob_log_effect_greater_0 = mean(x > 0, na.rm = TRUE),
        note = NA_character_
      )
    })
    
bayesian_league_long <- relative_effects_logscale %>%
  mutate(
    estimate = ifelse(
      all(c("log_effect_median", "log_effect_lCrI", "log_effect_uCrI") %in% names(.)),
      paste0(
        round(log_effect_median, 4),
        " (",
        round(log_effect_lCrI, 4),
        ", ",
        round(log_effect_uCrI, 4),
        ")"
      ),
      NA_character_
    )
  )
    
    bayesian_league_matrix <- {
      mat <- matrix(
        "",
        nrow = length(treatments),
        ncol = length(treatments),
        dimnames = list(treatments, treatments)
      )
      
      diag(mat) <- "Reference"
      
      for (i in seq_len(nrow(bayesian_league_long))) {
        a <- bayesian_league_long$treat1[i]
        b <- bayesian_league_long$treat2[i]
        est <- bayesian_league_long$estimate[i]
        
        mat[a, b] <- est
        mat[b, a] <- paste0("reverse of ", est)
      }
      
      as.data.frame(mat) %>%
        rownames_to_column("treatment")
    }
    
    abs_risk_draws <- get_abs_risk_draws(
      base_res = base_res,
      arm_data = arm_data,
      treatments = treatments,
      n_draws = n_iter * n_chain
    )
    
    ni_tbl <- map_dfr(seq_len(nrow(comps)), function(i) {
      a <- comps[i, 1]
      b <- comps[i, 2]
      
      p_a <- abs_risk_draws[[a]]
      p_b <- abs_risk_draws[[b]]
      
      m <- min(length(p_a), length(p_b))
      p_a <- p_a[seq_len(m)]
      p_b <- p_b[seq_len(m)]
      
      threshold_b <- pmin(1, p_b * (1 + noninf_margin_relative))
      threshold_a <- pmin(1, p_a * (1 + noninf_margin_relative))
      
      if (lower_is_better) {
        prob_a_noninf_vs_b <- mean(p_a <= threshold_b, na.rm = TRUE)
        prob_b_noninf_vs_a <- mean(p_b <= threshold_a, na.rm = TRUE)
        rule <- "lower risk better: P(risk_test <= min(1, risk_comparator * 1.10))"
      } else {
        prob_a_noninf_vs_b <- mean(p_a >= p_b * (1 - noninf_margin_relative), na.rm = TRUE)
        prob_b_noninf_vs_a <- mean(p_b >= p_a * (1 - noninf_margin_relative), na.rm = TRUE)
        rule <- "higher response better: P(risk_test >= risk_comparator * 0.90)"
      }
      
      tibble(
        treat1 = a,
        treat2 = b,
        margin_type = "relative_10_percent",
        decision_rule = rule,
        
        risk_treat1_median = median(p_a, na.rm = TRUE),
        risk_treat1_lCrI = q025(p_a),
        risk_treat1_uCrI = q975(p_a),
        
        risk_treat2_median = median(p_b, na.rm = TRUE),
        risk_treat2_lCrI = q025(p_b),
        risk_treat2_uCrI = q975(p_b),
        
        threshold_for_treat1_vs_treat2_median = median(threshold_b, na.rm = TRUE),
        threshold_for_treat2_vs_treat1_median = median(threshold_a, na.rm = TRUE),
        
        prob_treat1_noninferior_to_treat2 = prob_a_noninf_vs_b,
        prob_treat2_noninferior_to_treat1 = prob_b_noninf_vs_a
      )
    })
    
    ranks <- tryCatch({
      rp <- gemtc::rank.probability(
        base_res,
        preferredDirection = ifelse(lower_is_better, -1, 1)
      )
      
      rank_mat <- as.matrix(rp)
      
      rank_df <- as.data.frame(rank_mat)
      rank_df <- tibble::rownames_to_column(rank_df, "treatment")
      
      names(rank_df)[-1] <- paste0("rank_", seq_len(ncol(rank_df) - 1))
      
      rank_df <- rank_df %>%
        mutate(
          treatment = as.character(treatment),
          across(starts_with("rank_"), as.numeric)
        )
      
      rank_df
      
    }, error = function(e) {
      tibble(note = paste("rank.probability failed:", e$message))
    })
    
    sucra_objects <- tryCatch({
      build_sucra(ranks)
    }, error = function(e) {
      list(
        cumulative_tbl = tibble(note = paste("cumulative ranking failed:", e$message)),
        sucra_tbl = tibble(note = paste("SUCRA failed:", e$message))
      )
    })
    
    cumulative_ranking <- sucra_objects$cumulative_tbl
    ranking_summary <- sucra_objects$sucra_tbl
    
    nodesplit_summary <- tryCatch({
      nsmod <- gemtc::mtc.nodesplit(network)
      nsres <- gemtc::mtc.run(nsmod, n.adapt = n_adapt, n.iter = n_iter)
      as.data.frame(summary(nsres))
    }, error = function(e) {
      tibble(
        note = paste(
          "Node-splitting failed or no eligible independent closed loops:",
          e$message
        )
      )
    })
    
    ume_fit <- tryCatch({
      ume_model <- gemtc::mtc.model(
        network,
        type = "ume",
        linearModel = selected_linear_model,
        likelihood = "binom",
        link = "logit",
        n.chain = n_chain
      )
      
      set.seed(seed)
      ume_res <- gemtc::mtc.run(
        ume_model,
        n.adapt = n_adapt,
        n.iter = n_iter
      )
      
      ume_dic <- ume_res$deviance
      
      consistency_dic <- if (selected_linear_model == "random") dic_re else dic_fe
      
      tbl <- tibble(
        model = c("selected_consistency", "UME_inconsistency"),
        linear_model = selected_linear_model,
        DIC = c(consistency_dic$DIC, ume_dic$DIC),
        pD = c(consistency_dic$pD, ume_dic$pD),
        residual_deviance = c(
          consistency_dic$Dbar,
          ume_dic$Dbar),
        n_data_points = nrow(arm_data),
        delta_DIC_UME_minus_consistency = ume_dic$DIC - consistency_dic$DIC,
        interpretation = case_when(
          ume_dic$DIC + 3 < consistency_dic$DIC ~
            "UME fits meaningfully better: possible global inconsistency",
          consistency_dic$DIC + 3 < ume_dic$DIC ~
            "Consistency model fits meaningfully better",
          TRUE ~
            "No meaningful DIC difference; inspect residual deviance and local checks"
        )
      )
      
      list(
        table = tbl,
        model = ume_model,
        result = ume_res,
        dic = ume_dic
      )
    }, error = function(e) {
      list(
        table = tibble(
          note = paste(
            "UME model failed or is unsupported for this network/software setup:",
            e$message
          )
        ),
        model = NULL,
        result = NULL,
        dic = NULL
      )
    })
    
    global_inconsistency <- ume_fit$table
    
    transitivity_tbl <- if (!is.null(effect_modifiers)) {
      
      dat_conn %>%
        select(any_of(c("study", "treatment", effect_modifiers))) %>%
        mutate(across(any_of(effect_modifiers), as.character)) %>%
        pivot_longer(
          cols = any_of(effect_modifiers),
          names_to = "effect_modifier",
          values_to = "value"
        ) %>%
        group_by(treatment, effect_modifier) %>%
        summarise(
          n_nonmissing = sum(!is.na(value)),
          n_unique = n_distinct(value, na.rm = TRUE),
          values_observed = paste(sort(unique(na.omit(value))), collapse = "; "),
          .groups = "drop"
        )
      
    } else {
      tibble(
        note = "No effect_modifiers supplied. Transitivity/exchangeability cannot be statistically proven; assess clinically and descriptively."
      )
    }
    
    net_fig <- file.path(route_path, paste0(file_prefix, "_", label, "_network.png"))
    
    png(net_fig, width = 1800, height = 1400, res = 200)
    plot(nm, plastic = FALSE, thickness = "number.of.studies")
    dev.off()
    
    rank_fig <- file.path(route_path, paste0(file_prefix, "_", label, "_ranking_probabilities.png"))
    rankogram_fig <- file.path(route_path, paste0(file_prefix, "_", label, "_rankogram.png"))
    sucra_fig <- file.path(route_path, paste0(file_prefix, "_", label, "_SUCRA_curve.png"))
    
    if (!"note" %in% names(ranks)) {
      
      rank_long <- ranks %>%
        select(treatment, starts_with("rank_")) %>%
        pivot_longer(
          cols = starts_with("rank_"),
          names_to = "rank",
          values_to = "probability"
        ) %>%
        mutate(rank = as.integer(str_remove(rank, "rank_")))
      
      p_rank <- ggplot(rank_long, aes(x = rank, y = probability, group = treatment)) +
        geom_line() +
        geom_point() +
        facet_wrap(~ treatment) +
        theme_bw() +
        labs(
          x = "Rank",
          y = "Posterior probability",
          title = paste("Ranking probabilities:", label)
        )
      
      ggsave(rank_fig, p_rank, width = 10, height = 7, dpi = 300)
      
      p_rankogram <- ggplot(rank_long, aes(x = factor(rank), y = probability, fill = treatment)) +
        geom_col(position = "stack") +
        theme_bw() +
        labs(
          x = "Rank",
          y = "Posterior probability",
          fill = "Treatment",
          title = paste("Rankogram:", label)
        )
      
      ggsave(rankogram_fig, p_rankogram, width = 10, height = 7, dpi = 300)
    }
    
    if (!"note" %in% names(cumulative_ranking) &&
        nrow(cumulative_ranking) > 0 &&
        "treatment" %in% names(cumulative_ranking)) {
      
      sucra_long <- cumulative_ranking %>%
        mutate(treatment = as.character(treatment)) %>%
        pivot_longer(
          cols = starts_with("cum_rank_"),
          names_to = "rank",
          values_to = "cumulative_probability"
        ) %>%
        mutate(rank = as.integer(str_remove(rank, "cum_rank_")))
      
      p_sucra <- ggplot(
        sucra_long,
        aes(x = rank, y = cumulative_probability, group = treatment)
      ) +
        geom_line() +
        geom_point() +
        theme_bw() +
        labs(
          x = "Rank threshold",
          y = "Cumulative probability",
          title = paste("Cumulative ranking / SUCRA curve:", label)
        )
      
      ggsave(sucra_fig, p_sucra, width = 10, height = 7, dpi = 300)
    }

    
    list(
      label = label,
      arm_data = arm_data,
      arm_audit = arm_audit,
      connectivity = connectivity,
      network_edges = network_edges,
      heterogeneity = heterogeneity_tbl,
      model_fit = model_fit,
      selected_model = tibble(selected_model = selected_model),
      convergence = convergence,
      transitivity = transitivity_tbl,
      relative_effects_logscale = relative_effects_logscale,
      bayesian_league_long = bayesian_league_long,
      bayesian_league_matrix = bayesian_league_matrix,
      direct_indirect_mixed_netmeta = netsplit_tbl,
      noninferiority_10pct_relative = ni_tbl,
      ranking_probabilities = ranks,
      cumulative_ranking = cumulative_ranking,
      ranking_summary = ranking_summary,
      node_splitting = nodesplit_summary,
      global_inconsistency_UME = global_inconsistency,
      figures = tibble(
        network_plot = net_fig,
        ranking_plot = ifelse(file.exists(rank_fig), rank_fig, NA_character_),
        rankogram_plot = ifelse(file.exists(rankogram_fig), rankogram_fig, NA_character_),
        sucra_curve_plot = ifelse(file.exists(sucra_fig), sucra_fig, NA_character_)
      ),
      objects = list(
        network = network,
        model_re = model_re,
        model_fe = model_fe,
        result_re = res_re,
        result_fe = res_fe,
        selected_result = base_res,
        ume_model = ume_fit$model,
        ume_result = ume_fit$result,
        netmeta = nm
      )
    )
  }
  
  split_data <- if (is.null(group_vars)) {
    list(overall = data)
  } else {
    data %>%
      group_split(across(all_of(group_vars)), .keep = TRUE) %>%
      setNames(
        data %>%
          distinct(across(all_of(group_vars))) %>%
          tidyr::unite("grp", all_of(group_vars), sep = "_", remove = FALSE) %>%
          pull(grp)
      )
  }
  
  results <- purrr::imap(split_data, analyse_one)
  
  wb <- openxlsx::createWorkbook()
  
  add_sheet <- function(wb, sheet, df) {
    sheet <- safe_sheet(sheet)
    
    original <- sheet
    i <- 1
    
    while (sheet %in% names(wb)) {
      suffix <- paste0("_", i)
      sheet <- substr(paste0(substr(original, 1, 31 - nchar(suffix)), suffix), 1, 31)
      i <- i + 1
    }
    
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, df)
  }
  
  for (nm_i in names(results)) {
    r <- results[[nm_i]]
    prefix <- substr(gsub("[^A-Za-z0-9_]", "_", nm_i), 1, 10)
    
    add_sheet(wb, paste0(prefix, "_arm_data"), r$arm_data)
    add_sheet(wb, paste0(prefix, "_network"), r$network_edges)
    add_sheet(wb, paste0(prefix, "_connect"), r$connectivity)
    add_sheet(wb, paste0(prefix, "_heterog"), r$heterogeneity)
    add_sheet(wb, paste0(prefix, "_fit"), r$model_fit)
    add_sheet(wb, paste0(prefix, "_selected"), r$selected_model)
    add_sheet(wb, paste0(prefix, "_conv"), r$convergence)
    add_sheet(wb, paste0(prefix, "_transit"), r$transitivity)
    add_sheet(wb, paste0(prefix, "_log_effects"), r$relative_effects_logscale)
    add_sheet(wb, paste0(prefix, "_bayes_league_L"), r$bayesian_league_long)
    add_sheet(wb, paste0(prefix, "_bayes_league_M"), r$bayesian_league_matrix)
    add_sheet(wb, paste0(prefix, "_direct_ind"), r$direct_indirect_mixed_netmeta)
    add_sheet(wb, paste0(prefix, "_NI_10pct"), r$noninferiority_10pct_relative)
    add_sheet(wb, paste0(prefix, "_ranks"), r$ranking_probabilities)
    add_sheet(wb, paste0(prefix, "_cumrank"), r$cumulative_ranking)
    add_sheet(wb, paste0(prefix, "_SUCRA"), r$ranking_summary)
    add_sheet(wb, paste0(prefix, "_nodesplit"), r$node_splitting)
    add_sheet(wb, paste0(prefix, "_UME"), r$global_inconsistency_UME)
    add_sheet(wb, paste0(prefix, "_figures"), r$figures)
  }
  
  xlsx_path <- file.path(route_path, paste0(file_prefix, "_results.xlsx"))
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  
  message("Excel exported to: ", xlsx_path)
  
  invisible(list(
    excel_path = xlsx_path,
    results = results
  ))
}