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
scale_prs <- function(data = analysis_df,
                      pheno_col = "vitiligo",
                      pheno_cond = 0,
                      scale_columns = c("mhc_class2_only_prs", "CONFIRMED_prs")) {
  for (i in scale_columns) {
    mean_val <- mean(data[data[[pheno_col]] == pheno_cond, i], na.rm = TRUE)
    std_dev_val <- sd(data[data[[pheno_col]] == pheno_cond, i], na.rm = TRUE)
    
    new_col_name <- paste0(i, "_scaled")
    data <- data %>% mutate(!!sym(new_col_name) := (data[[i]] - mean_val) / std_dev_val)
  }
  return(data)
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
  model <- glm(formula, data = data)
  
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
      upper_ci = estimate + qnorm(0.975) * std_error
    )
  
  return(tidy_output)
}


# Function to compute the interaction term
regression_interaction <- function(data = analysis_df,
                                       outcome,
                                       predictor,
                                       interactor,
                                       covariates = kCovariates) {
  
  #setup the formula
  formula <- as.formula(paste(outcome, "~",
                              paste(c(covariates,
                                      paste0(predictor, "*", interactor)),
                                    collapse = " + ")))
  
  # Check if the outcome variable is binary
  is_binary <- all(data[[outcome]] %in% c(0, 1))
  if (is_binary) {
    # Binary outcome variable, perform logistic regression
    model <- glm(formula, data = data, family = "binomial")
  } else {
    # Continuous outcome variable, perform linear regression
    model <- lm(formula, data = data)
  }
  
  # Create a tidy data frame
  tidy_output <- broom::tidy(model) %>% 
    transmute(
      outcome,
      predictor,
      term,
      pval = `p.value`,
      estimate,
      std_error = `std.error`,
      lower_ci = estimate - qnorm(0.975)*std_error,
      upper_ci = estimate + qnorm(0.975)*std_error
    )
  
  return(tidy_output)
}

