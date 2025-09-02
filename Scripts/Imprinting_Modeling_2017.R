library(tidyr)
library(dplyr)
library(devtools)
library(imprinting)
library(tidyverse)
library(DescTools)


nrep <- 1000 # Number of bootstrap replicates

# retrieve imprinting probabilities
imprinting_usa_2017 <- get_imprinting_probabilities(observation_years = 2017, countries = "United States")
imprinting_usa_2017 <- subset(spread(imprinting_usa_2017, birth_year, imprinting_prob), select = -c(year, country))
imprinting_usa_2017 <- cbind(V0 = colnames(imprinting_usa_2017), t(imprinting_usa_2017))

write.csv(imprinting_df_2017, "Ort-et-al_N1-imprinting/Data/2017_Imprinting.csv")

# read in titer csv, with manually merged imprinting probabilities
titer_df <- read.csv("Ort-et-al_N1-imprinting/Data/2017_Titers_Imprinting.csv")

# filter to adults
titer_df_adults <- titer_df %>% filter(age >= 18)


# Given the data, and for a single antigen, fits separate linear models predicting
# log titers based on age, year of birth, group1 imprinting probability, and H1N1
# imprinting probability
fit_linear_models <- function(data, antigen){
  
  data <- data %>% mutate(log2_titer = log(titer, base = 2)) %>%
    filter(antigen == !! antigen)
  
  age_model <- lm(log2_titer ~ age, data = data)
  yob_model <- lm(log2_titer ~ yob, data = data)
  group1_ha_imprinting_model <- lm(log2_titer ~ imp_group1, data = data)
  group2_ha_imprinting_model <- lm(log2_titer ~ imp_group2, data = data)
  n1_imprinting_model <- lm(log2_titer ~ imp_n1, data = data)
  n2_imprinting_model <- lm(log2_titer ~ imp_n2, data = data)
  
  AIC_table <- AIC(age_model, yob_model, group1_ha_imprinting_model, group2_ha_imprinting_model,
                   n1_imprinting_model, n2_imprinting_model) %>%
    as_tibble(rownames = "model") %>%
    arrange(AIC) %>%
    mutate(loglik = sapply(model,
                           FUN = function(x){
                             logLik(get(x)) %>% as.numeric()
                           },
                           simplify = T, USE.NAMES = F)) %>%
    mutate(delta_AIC = AIC - min(AIC)) %>%
    mutate(antigen = antigen) %>%
    select(antigen, everything()) %>%
    mutate(model = str_remove(model, "_model"))
  
  return(AIC_table)
}

# Compute Spearman correlations between titers to each antigen and various predictors
compute_spearman_correlations <- function(data){
  
  candidate_correlates <- c("imp_group1", "imp_group2", "imp_n1", "imp_n2", "age", "yob")
  
  data %>%
    group_by(antigen) %>%
    summarise(across(candidate_correlates, function(x){cor.test(x, titer, method = 'spearman')$estimate},
                     .names = "coef_{.col}"),
              across(candidate_correlates, function(x){cor.test(x, titer, method = 'spearman')$p.value},
                     .names = "pvalue_{.col}"),
              n_obs = n()
    ) %>%
    pivot_longer(cols = !any_of(c("antigen","n_obs") )) %>%
    mutate(value_type = str_extract(name, "[^_]+")) %>%
    mutate(name = str_remove(name, paste0(value_type,"_"))) %>%
    rename(predictor = name) %>%
    pivot_wider(names_from = value_type, values_from = value) %>%
    rowwise() %>%
    mutate(lower = CorCI(coef, n_obs)['lwr.ci'],
           upper = CorCI(coef, n_obs)['upr.ci']) %>%
    ungroup()
}

# Takes difference in Spearman correlation between predictors 1 and 2,#
# tests their significance using bootstrapping
run_bootstrap_correlation_test <- function(data, predictor1, predictor2, nrep){
  
  # Internal function: given table of Spearman correlation results, 
  # retrieves comparison of interest, puts in wide format
  get_pw_comparison <- function(spearman_corrs, predictor1, predictor2){
    spearman_corrs %>%
      filter(predictor %in% c(predictor1, predictor2)) %>%
      select(antigen, predictor, coef) %>%
      mutate(predictor = case_match(predictor,
                                    predictor1 ~ 1,
                                    predictor2 ~ 2 )) %>%
      pivot_wider(names_from = predictor, values_from = coef, names_prefix = "cor_predictor_") %>%
      # We care about the absolute value of the correlation coefficient
      mutate(across(!matches("antigen"), abs)) %>%
      select(antigen, cor_predictor_1, cor_predictor_2) %>%
      mutate(cor_diff = cor_predictor_1 - cor_predictor_2)
  }
  
  obs_corr_diff <- compute_spearman_correlations(data) %>%
    get_pw_comparison(predictor1 = predictor1, predictor2 = predictor2) %>%
    rename_with(.fn = function(x){paste0(x,"_obs")}, .cols = !matches("antigen"))
  
  bootstrap_corr_diff <- c()
  
  for(i in 1:nrep){
    for(antigen in unique(data$antigen)){
      
      row_indices <- (1:nrow(data))[data$antigen == antigen]
      
      resampled_rows <- sample(row_indices, replace = T)
      
      resampled_data <- data[resampled_rows,]
      
      resampled_corr <- compute_spearman_correlations(resampled_data) %>%
        get_pw_comparison(predictor1 = predictor1, predictor2 = predictor2)
      
      bootstrap_corr_diff <- bootstrap_corr_diff %>%
        bind_rows(resampled_corr)
      
    }
    
  }
  
  
  bootstrap_results <- bootstrap_corr_diff %>% 
    left_join(obs_corr_diff) %>%
    group_by(antigen, cor_predictor_1_obs, cor_predictor_2_obs, cor_diff_obs) %>%
    summarise(bootstrap_diff_mean = mean(cor_diff),
              bootstrap_diff_lower = quantile(cor_diff, 0.025),
              bootstrap_diff_upper = quantile(cor_diff, 0.975)) %>%
    mutate(predictor1 = predictor1, predictor2 = predictor2) %>%
    select(antigen, predictor1, predictor2, everything())
  
  return(bootstrap_results)
  
}

label_antigens <- function(data){
  data %>%
    mutate(antigen = case_match(
      antigen,
      "ca09" ~ "CA09",
      "vn04" ~ "VN04",
      "dc24" ~ "DC24",
      "bc24" ~ "BC24"
    ))
  
}

label_predictor <- function(data){
  data %>%
    mutate(predictor = case_match(
      predictor,
      "imp_group1" ~ "Group 1 imprinting",
      "imp_group2" ~ "Group 2 imprinting",
      "imp_n1" ~ "N1 imprinting",
      "imp_n2" ~ "N2 imprinting",
      "yob" ~ "Birth year",
      "age" ~ "Age"
    ))
}

format_table <- function(spearman_cors_table){
  spearman_cors_table %>%
    filter(antigen %in% c("ca09", "vn04", "dc24", "bc24")) %>%
    label_antigens() %>%
    label_predictor() %>%
    select(-n_obs) %>%
    mutate(Correlation = paste(antigen, predictor, sep = ' titers - '),
           `Coefficient (95% CI)` = paste0(
             round(coef, 2), " (", round(lower,2), ",", round(upper,2), ")"),
           `P value` = pvalue) %>%
    select(Correlation, `Coefficient (95% CI)`, `P value`)
  
}

# For the 2017 cohort, compute Spearman correlations between titers to
# each antigen and imprinting probabilities and year of birth
spearman_cors_yob <- compute_spearman_correlations(
  titer_df_adults) %>%
  filter(predictor != 'age') %>%
  format_table()

write_csv(spearman_cors_yob, "Ort-et-al_N1-imprinting/Tables/Extended_Table1.csv")


# As a complement to correlations and linear models, run bootstrap tests

# Is N1 imprinting probability more strongly associated with titers than yob is?
bootstrap_n1_vs_yob <- run_bootstrap_correlation_test(titer_df_adults,
                                                          predictor1 = "imp_n1",
                                                          predictor2 = "yob",
                                                          nrep = nrep)

# Does group 1 imprinting beat YOB?
bootstrap_group1_vs_yob <- run_bootstrap_correlation_test(titer_df_adults,
                                                          predictor1 = "imp_group1",
                                                          predictor2 = "yob",
                                                          nrep = nrep)

# Does N1 imprinting beat group 1 imprinting?
bootstrap_n1_vs_group1 <- run_bootstrap_correlation_test(titer_df_adults,
                                                          predictor1 = "imp_n1",
                                                          predictor2 = "imp_group1",
                                                          nrep = nrep)

bootstrap_results <- 
  bootstrap_n1_vs_yob %>%
  bind_rows(bootstrap_group1_vs_yob) %>%
  bind_rows(bootstrap_n1_vs_group1)

bootstrap_results <- bootstrap_results %>%
  label_antigens() %>%
  rename(predictor = predictor1) %>%
  label_predictor() %>%
  rename(predictor1 = predictor, predictor = predictor2) %>%
  label_predictor() %>%
  rename(predictor2 = predictor) %>%
  arrange(antigen) %>%
  mutate(across(where(is.numeric), ~round(.x,3))) %>%
  mutate(`Bootstrap difference (95% CI)` = 
           paste0(
             bootstrap_diff_mean, 
             " (",
             bootstrap_diff_lower,
             ",",
             bootstrap_diff_upper,
             ")"
           )) %>%
  rename(Antigen = antigen,
         `Predictor 1` = predictor1,
         `Predictor 2` = predictor2,
         `Titer correlation with predictor 1` = cor_predictor_1_obs,
         `Titer correlation with predictor 2` = cor_predictor_2_obs,
         `Observed difference` = cor_diff_obs) %>%
  select(Antigen, `Predictor 1`, `Predictor 2`, `Titer correlation with predictor 1`,
         `Titer correlation with predictor 2`, `Observed difference`,
         `Bootstrap difference (95% CI)`)

write_csv(bootstrap_results,  "Ort-et-al_N1-imprinting/Tables/Extended_Table2.csv")


# Run the models for the 2017 cohort and export table with the results
model_comparison <- sapply(unique(titer_df_adults$antigen), FUN = fit_linear_models,
                           data = titer_df_adults, simplify = F) %>% bind_rows() %>%
  rename(k = df)

model_comparison <- model_comparison %>%
  label_antigens() %>%
  mutate(model = case_match(
    model,
    "yob" ~ "Birth year",
    "group1_ha_imprinting" ~ "Group 1 imprinting",
    "group2_ha_imprinting" ~ "Group 2 imprinting",
    "n1_imprinting" ~ "N1 imprinting",
    "n2_imprinting" ~ "N2 imprinting",
    "age" ~ "Age"
  )) %>%
  mutate(across(any_of(c('AIC','loglik','delta_AIC')), ~round(.x,2)))

write_csv(model_comparison, "Ort-et-al_N1-imprinting/Tables/Extended_Table3.csv")
