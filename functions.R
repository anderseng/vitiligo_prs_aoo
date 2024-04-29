# Function to recode specified SNPs so that risk alleles match weights
recode_snps <- function(df, snps_to_recode) {
  df %>%
    mutate(across(all_of(snps_to_recode), 
                  ~ ifelse(. == 0, 2,
                           ifelse(. == 2, 0,
                                  ifelse(. == 1, 1, NA)))))
}

# Function to perform matrix multiplication to generate the main PRS
calculate_main_prs <- function(risk_genos = risk_genos_recoded,
                               weights = weights,
                               iid_col = "IID",
                               weights_col = "Weights_risk_direction_only") {
  weight_mat <- weights %>%
    column_to_rownames(var = "SNP") %>%
    select(weights = !!sym(weights_col)) %>%
    as.matrix()
  
  geno_mat <- risk_genos %>%
    select(rownames(weight_mat)) %>%
    as.matrix()
  
  rownames(geno_mat) <- risk_genos %>% pull(!!sym(iid_col))
  
  res <- geno_mat %*% weight_mat
  res <- as.data.frame(res) %>%
    rownames_to_column(., var="tmp_id") %>%
    transmute(!!sym(iid_col) := tmp_id,
              score := weights)
  return(res)
}

# Function to calculate Class II risk score
calculate_class2_risk_score <- function(risk_genos) {
  high_risk_hap_ln_or <- 2.0541237
  generic_risk_ln_or <- 0.5721089
  
  risk_genos %>%
    mutate(mhc_class2_only_prs = case_when(
      RS145954018 >= 1 ~ high_risk_hap_ln_or,
      RS114448410 == 2 ~ 2 * generic_risk_ln_or,
      RS114448410 == 1 ~ generic_risk_ln_or,
      RS114448410 == 0 & RS145954018 == 0 ~ 0,
      TRUE ~ NA_real_
    )) %>%
    select(IID, mhc_class2_only_prs)
}

# write a function that scales the PRS based on controls
library(dplyr)

scale_prs <- function(data = analysis_df,
                      pheno_col = "vitiligo",
                      pheno_cond = 0,
                      scale_columns = c("mhc_class2_only_prs", "CONFIRMED_prs")) {
  for (i in scale_columns) {
    mean_val <- mean(data[data[[pheno_col]] == pheno_cond, i], na.rm = TRUE)
    std_dev_val <- sd(data[data[[pheno_col]] == pheno_cond, i], na.rm = TRUE)
    
    new_col_name <- paste0(i, "_scaled")
    percentile_col_name <- paste0(i, "_percentile")
    
    # Calculate percentile
    control_prs <- na.omit(data[data[[pheno_col]] == pheno_cond, i])
    ecdf_func <- ecdf(control_prs)
    percentiles <- ifelse(data$case_control == "control", ecdf_func(data[[i]]), NA)
    
    data <- data %>% 
      mutate(!!sym(new_col_name) := (data[[i]] - mean_val) / std_dev_val,
             !!sym(percentile_col_name) := ecdf_func(!!sym(i)))
  }
  return(data)
}


# Write a function that creates normal distributions for plotting
generate_plot_dist_data <- function(dist_df = dist_params,
                                    name_col = "AOO_category",
                                    mean_col = "FMM_mean",
                                    sd_col = "FMM_sd",
                                    size = 1000) {
  set.seed(123)
  data <- data.frame()
  for (i in seq_len(nrow(dist_df))) {
    dist_name <- dist_df %>% slice(i) %>% pull(!!sym(name_col))
    dist_mean <- dist_df %>% slice(i) %>% pull(!!sym(mean_col))
    dist_sd <- dist_df %>% slice(i) %>% pull(!!sym(sd_col))
    
    dist_data <- data.frame(
      x = rnorm(size, mean = dist_mean, sd = dist_sd),
      group = rep(dist_name, each = size)
    )
    
    data <- rbind(data, dist_data)
  }
  return(data)
}

# write a plotting function for the AOO histogram with overlayed distributions
plot_distributions <- function(histogram_data = case_only_analysis_df,
                               histogram_variable = "VITageonset",
                               dist_df = dist_params,
                               name_col = "AOO_category",
                               mean_col = "FMM_mean",
                               sd_col = "FMM_sd",
                               size = dim(case_only_analysis_df)[1],
                               binwidth = 2) {
  # Generate data using the function
  density_data <- generate_plot_dist_data(dist_df, name_col, mean_col, sd_col, size)
  
  # Assuming you have another data frame for the independent data
  independent_data <- tibble(x = histogram_data[[histogram_variable]])
  
  vlines <- purrr::map2(density_data[[mean_col]], density_data[[name_col]],
                        ~geom_vline(xintercept = .x,
                                    color = .y,
                                    linetype = "dashed",
                                    size = 1))
  # Create the plot
  ggplot() +
    geom_histogram(data = independent_data, aes(x = x, y = ..count..), 
                   binwidth = binwidth, fill = "gray", alpha = 0.7) +
    geom_density(data = density_data, aes(x = x, y = ..count.., fill=group),
                 alpha = 0.2, adjust = 2) +
    geom_vline(data = dist_df, aes(xintercept = !!sym(mean_col)),
               color = "black", linetype = "dashed", size = 0.5) +
    labs(title = "Early- and Late-Onset Groups Defined by FMM",
         x = "Vitiligo Age-of-Onset",
         y = "Count") +
    vlines +
    theme_bw()
}




# Function to compute logistic regression results
logistic_regression <- function(data = analysis_df,
                                pheno,
                                prs,
                                covariates = kCovariates) {
  formula <- as.formula(paste(pheno, "~",
                              paste(c(prs, covariates), collapse = " + ")))
  model <- glm(formula, data = data, family = binomial)
  
  # Create a tidy data frame
  tidy_output <- model %>%
    tidy() %>% 
    transmute(
    Phenotype = pheno,
    PRS = prs,
    term = term,
    pval = `p.value`,
    estimate,
    std_error = `std.error`,
    OR = exp(estimate),
    lower_ci = exp(estimate - qnorm(0.975) * std_error),
    upper_ci = exp(estimate + qnorm(0.975) * std_error)
  )
  
  return(tidy_output)
}


# Function to compute linear regression results
aoo_linear_regression <- function(data = analysis_df,
                                pheno,
                                prs,
                                covariates = kCovariates) {
  formula <- as.formula(paste(pheno, "~",
                              paste(c(prs, covariates), collapse = " + ")))
  model <- lm(formula, data = data)
  
  
  # Create a tidy data frame
  tidy_output <- model %>%
    tidy() %>% 
    transmute(
      Phenotype = pheno,
      PRS = prs,
      term,
      pval = `p.value`,
      estimate,
      std_error = `std.error`,
      lower_ci = estimate - qnorm(0.975) * std_error,
      upper_ci = estimate + qnorm(0.975) * std_error,
      n_obs = glance(model)$nobs,
      adj_r_sq = glance(model)$adj.r.squared
    )
  
  return(tidy_output)
}
