library("stats")
library("rstatix")
library("dplyr")
library("tidyverse")
library("tidyr")
library("broom")
library("dunn.test")
library("patchwork")
library("vegan")
library("data.table")
library("ggplot2")
library("ggrepel")
library("ggforce")
library("concaveman")
library("gridExtra")
library("stringr")
library("shades")
library("ggpubr")
library("ggokabeito")
library("multcomp")
library("car")
library("asbio")
library("plotrix")
library("lmtest")
library("boot")
library("moments")
library("lme4")
library("DAAG")
library("performance")
library("emmeans")
library("purrr")
library("DT")
library("openxlsx")

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 8),
      panel.grid.major = element_line(linewidth = 0.1, color = "#cccccc"),
      panel.grid.minor = element_line(linewidth = 0.05, color = "#cccccc"),
      strip.background.y = ggh4x::element_part_rect(side = "l", fill = NA),
      strip.text.y = element_text(size = 6, face = "bold"),
      strip.background.x = ggh4x::element_part_rect(side = "b", fill = NA),
      strip.text.x = element_text(size = 6, face = "bold"),
      legend.position = "bottom",
      legend.title.position = "top",
      legend.key.spacing.x = unit(0, "mm"),
      legend.key.spacing.y = unit(0, "mm"),
      legend.key.size = unit(4, "mm")
    )
)



wcmd_and_adonis <- function(df, distance = "bray", exponent = 0.5) {
  asvs_to_remove <-
    df %>%
    dplyr::group_by(asv) %>%
    dplyr::summarise(
      total_reads = sum(reads)
    ) %>% 
    dplyr::filter(total_reads <= 5) %>%
    dplyr::pull(asv)
  
  df_wide <- 
    df %>%
    dplyr::filter(
      !asv %in% asvs_to_remove
    ) %>%
    dplyr::mutate(
      present = as.numeric(reads > 0)
    ) %>%
    dplyr::select(dplyr::any_of(c("asv", "sample", "exp", "day", "RA", "media", "chemical", "concentration_ug_L"))) %>%
    tidyr::pivot_wider(
      values_from = RA,
      names_from = asv
    ) %>%
    dplyr::arrange(concentration_ug_L)
  
  species <-
    df_wide %>%
    dplyr::select(!dplyr::any_of(c("sample", "exp", "media", "chemical", "concentration_ug_L", "day")))
  
  env <-
    df_wide %>%
    dplyr::select(dplyr::any_of(c("sample", "exp", "media", "chemical", "concentration_ug_L", "day")))
  
  wcmd <- 
    vegan::wcmdscale(
      vegdist(species),
      k = 2
    )
  
  colnames(wcmd) <- c("WCMD1", "WCMD2")
  
  remove_diuron <- env$chemical != "diuron"
  
  permanova <-
    adonis2(
      species[remove_diuron,]^exponent ~ concentration_ug_L,
      data = env[remove_diuron,],
      method = distance,
      by = "terms",
      permutations = 9999
    )
  
  list(
    adonis = permanova,
    wcmd = wcmd %>%
      dplyr::bind_cols(env)
  )
}

simper_summary <- function(
    df,
    resp = "RA",
    exponent = 0.5,
    tax_level = "phylum",
    n_taxa = 5
) {
  asvs_to_remove <-
    df %>%
    dplyr::filter(chemical != "diuron") %>%
    dplyr::group_by(asv) %>%
    dplyr::summarise(
      total_reads = sum(reads)
    ) %>% 
    dplyr::filter(total_reads <= 5) %>%
    dplyr::pull(asv)
  
  simper_df <-
    df %>%
    dplyr::filter(chemical != "diuron") %>%
    dplyr::filter(
      !asv %in% asvs_to_remove
    ) %>%
    replace(is.na(.), "Unassigned") %>%
    dplyr::arrange(concentration_ug_L)
  
  asvs_wide <-
    simper_df %>%
    dplyr::select(
      dplyr::any_of(c("asv", "exp", "media", "chemical", "concentration_ug_L", "day", resp))
    ) %>%
    tidyr::pivot_wider(names_from = asv, values_from = !!sym(resp))
  
  species <-
    asvs_wide %>%
    dplyr::select(-c(1:5))
  
  groups <-
    asvs_wide %>%
    dplyr::select(1:5)
  
  # see simper documentation for meaning of columns
  simper_object <-
    vegan::simper(
      species^exponent,
      groups$concentration_ug_L
    )
  
  subset_simper <-
    paste0(
      min(simper_df$concentration_ug_L),
      "_",
      max(simper_df$concentration_ug_L)
    )
  
  df_simper <-
    tibble::tibble(
      asv = names(simper_object[[subset_simper]]$average),
      avg = simper_object[[subset_simper]]$average
    ) 
  
  df_simper_summary <-
    df_simper %>%
    dplyr::left_join(
      simper_df %>%
        dplyr::select(asv:species) %>%
        dplyr::distinct()
    ) %>%
    dplyr::group_by(!!sym(tax_level)) %>%
    dplyr::summarise(
      simper_avg = sum(avg) * 100
    ) %>%
    dplyr::arrange(dplyr::desc(simper_avg))
  
  simper_names <- df_simper_summary %>% dplyr::pull(!!sym(tax_level))
  simper_avg <- df_simper_summary %>% dplyr::pull(simper_avg)
  
  plt_data <-
    simper_df %>%
    dplyr::filter(chemical != "diuron") %>%
    dplyr::left_join(df_simper, by = "asv") %>%
    dplyr::filter(
      !!sym(tax_level) %in% simper_names[n_taxa:1]
    ) %>%
    dplyr::group_by(asv) %>%
    dplyr::mutate(
      rel_change = RA - mean(RA[concentration_ug_L == 0], na.rm = TRUE),
      !!tax_level := factor(!!sym(tax_level), levels = simper_names[n_taxa:1]) %>%
        forcats::fct_relabel(
          ~ paste0(., " <span style = \"font-family: Monaco; font-size: 8px\">(", 
                   format(simper_avg[n_taxa:1], digits = 1), "%)</span>")
        )
    ) 
  
  plt_data_summary <-
    plt_data %>%
    dplyr::group_by(!!sym(tax_level)) %>%
    dplyr::mutate(rel_abund_ctrl = mean(RA[concentration_ug_L == 0])) %>%
    dplyr::group_by(!!sym(tax_level), concentration_ug_L) %>%
    dplyr::summarise(
      weighted_rel_change = 1 + sum(rel_change) / sum(rel_abund_ctrl),
      # Mean - Confidence interval
      weighted_rel_change_lower = weighted_rel_change - (sd(rel_change) / sqrt(n())) / mean(rel_abund_ctrl),
      # Confidence interval cannot be lower than -1
      weighted_rel_change_lower = ifelse(
        weighted_rel_change_lower < 0,
        0, 
        weighted_rel_change_lower
      ), 
      weighted_rel_change_upper = weighted_rel_change + (sd(rel_change) / sqrt(n())) / mean(rel_abund_ctrl)
    ) %>%
    dplyr::filter(concentration_ug_L != 0)
  
  plt_data_summary
}