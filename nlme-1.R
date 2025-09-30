#-----------------------------------------------#
#-----------------------------------------------#
#               Read & prepare data             #
#-----------------------------------------------#
#-----------------------------------------------#

tb <- read_csv("data/data_ebola_ab.csv")

tb <- tb |>
  # clean the names
  janitor::clean_names() |>
  # transform chr into factors
  mutate(across(where(is.character), fct))

tb_clean <- tb |>
  filter(age >= 0)

tb_clean <- tb_clean |>
  mutate(age_updated = age + time / 365) |>
  relocate(age_updated, .after = age)

tb_after_7 <- tb_clean |>
  filter(time != 0)

# tb_clean_centered <- tb_clean |>
#   mutate(age_centered = age - mean(age)) |>
#   relocate(age_centered, .after = age) |>
#   select(-c(age, age_updated))
# 
# tb_after_7_centered <- tb_after_7 |>
#   mutate(age_centered = age - mean(age)) |>
#   relocate(age_centered, .after = age) |>
#   select(-c(age, age_updated))

tb_after_7 <- tb_after_7 |>
  mutate(
    id = parse_number(as.character(id))
  )
tb_after_7

#-----------------------------------------------#
#-----------------------------------------------#
#               No covariates model             #
#-----------------------------------------------#
#-----------------------------------------------#

saemix_data_nlme <- saemixData(
  name.data = tb_after_7,
  name.group = "id",
  name.predictors = "time",
  name.response = "abelisa",
  units = list(x = "Days", y = "EU/ml")
)
# saemix_data_nlme


function_approx <- function(psi, id, xidep){
  time <- xidep[, 1]
  y_min <- psi[id, 1]
  y_max <- psi[id, 2]
  rate <- psi[id, 3]
  fpred <- y_min + (y_max - y_min) * exp(-rate * time)
  return(fpred)
}
# function_approx


saemix_model_nlme <- saemixModel(
  model = function_approx,
  description = "Approx function", 
  psi0 = matrix(
    c(1.6, 4.9, 0.05),
    ncol = 3, 
    byrow = TRUE, 
    dimnames = list(NULL, c("y_min", "y_max", "rate"))),
  transform.par = c(0, 0, 1),   # log-transform rate only
  covariance.model = diag(1, 3),  # random effects on all
  error.model = "constant"
)
# saemix_model_nlme


saemix_options_nlme <- list(
  seed = 42,
  save = FALSE, 
  save.graph = FALSE,
  displayProgress = FALSE,
  print = FALSE
)
#saemix_options_nlme


saemix_fit_nlme <- saemix(
  saemix_model_nlme,
  saemix_data_nlme,
  saemix_options_nlme
)
# saemix_fit_nlme


saemix_fit_nlme@results


saemix_fit_nlme@results@fixed.effects


fit_with_preds <- saemix.predict(saemix_fit_nlme)
tb_data <- fit_with_preds@data@data
tb_predictions <- fit_with_preds@results@predictions

tb_complete <- cbind(tb_data, tb_predictions)

ggplot(subset(tb_complete, id == 1), aes(x = time)) +
  geom_point(aes(y = abelisa, color = "Observed"), size = 2) +
  geom_line(aes(y = ipred, color = "Individual prediction"), linewidth = 1.5) +
  geom_line(aes(y = ypred, color = "Population prediction"), linewidth = 1.5, linetype = "dashed") +
  scale_color_manual(
    name = NULL,
    values = c(
      "Observed" = "tomato",
      "Individual prediction" = "darkseagreen",
      "Population prediction" = "violet"
    )
  ) +
  labs(
    title = str_c("Subject", 1, ": Obs vs. Pred", sep = " "),
    x = "Time",
    y = "Antibody levels (EU/ml)"
  ) +
  theme_bw() +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
  )

ggplot(tb_complete, aes(x = time)) +
  geom_point(aes(y = abelisa), color = "tomato", size = 2) +
  geom_line(aes(y = ipred), color = "darkseagreen", linewidth = 1.5) +
  geom_line(aes(y = ypred), color = "violet", linewidth = 1.5, linetype = "dashed") +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  facet_wrap(~ id) +
  theme_bw()

psi <- psi(saemix_fit_nlme) |>
  as_tibble() |>
  rowid_to_column("iteration")

psi |>
  ggplot(aes(x = iteration, y = y_min)) +
  geom_line(color = "tomato") +
  labs(
    x = "Iteration",
    y = "y_min"    
  ) +
  theme_bw()

psi |>
  ggplot(aes(x = iteration, y = y_max)) +
  geom_line(color = "tomato") +
  labs(
    x = "Iteration",
    y = "y_max"    
  ) +
  theme_bw()

psi |>
  ggplot(aes(x = iteration, y = rate)) +
  geom_line(color = "tomato") +
  labs(
    x = "Iteration",
    y = "rate"    
  ) +
  theme_bw()






#-----------------------------------------------#
#-----------------------------------------------#
#               Covariates model                #
#-----------------------------------------------#
#-----------------------------------------------#

tb_after_7 <- tb_after_7 |>
  mutate(age_centered = age - mean(age, na.rm = TRUE))

function_approx <- function(psi, id, xidep) {
  time <- xidep[, 1]
  y_min <- psi[id, 1]
  y_max <- psi[id, 2]
  rate  <- psi[id, 3]
  fpred <- y_min + (y_max - y_min) * exp(-rate * time)
  return(fpred)
}

saemix_data_nlme <- saemixData(
  name.data = tb_after_7,
  name.group = "id",
  name.predictors = "time",
  name.response = "abelisa",
  name.covariates = "age_centered",
  units = list(x = "Days", y = "EU/ml")
)

psi0_matrix <- matrix(
  c(1.6, 4.9, 0.05,  # fixed effects: y_min, y_max, rate
    0.0, 0.0, 0.01), # covariate effect on rate only
  nrow = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("y_min", "y_max", "rate"))
)

covariate_model_matrix <- matrix(
  data = c(0, 0, 1),
  nrow = 1,
  ncol = 3,
  dimnames = list("age_centered", c("y_min", "y_max", "rate"))
)

saemix_model_nlme <- saemixModel(
  model = function_approx,
  description = "Antibody decay model with age covariate",
  psi0 = psi0_matrix,
  transform.par = c(0, 0, 1),           # log-transform rate
  covariance.model = diag(1, 3),
  covariate.model = covariate_model_matrix,
  error.model = "constant"
)

saemix_options_nlme <- list(
  seed = 42,
  save = FALSE,
  save.graph = FALSE,
  displayProgress = TRUE,
  print = FALSE
)

saemix_fit_nlme <- saemix(
  saemix_model_nlme,
  saemix_data_nlme,
  saemix_options_nlme
)

# Fixed effects and covariate effects
saemix_fit_nlme@results@fixed.effects


fit_with_preds <- saemix.predict(saemix_fit_nlme)

tb_data <- fit_with_preds@data@data
tb_predictions <- fit_with_preds@results@predictions
tb_complete <- cbind(tb_data, tb_predictions)

ggplot(subset(tb_complete, id == 1), aes(x = time)) +
  geom_point(aes(y = abelisa, color = "Observed"), size = 2) +
  geom_line(aes(y = ipred, color = "Individual prediction"), linewidth = 1.5) +
  geom_line(aes(y = ypred, color = "Population prediction"), linewidth = 1.5, linetype = "dashed") +
  scale_color_manual(
    name = NULL,
    values = c("Observed" = "tomato", "Individual prediction" = "darkseagreen", "Population prediction" = "violet")
  ) +
  labs(
    title = "Subject 1: Observed vs. Predicted",
    x = "Time (Days)",
    y = "Antibody Level (EU/ml)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggplot(tb_complete, aes(x = time)) +
  geom_point(aes(y = abelisa), color = "tomato", size = 2) +
  geom_line(aes(y = ipred), color = "darkseagreen", linewidth = 1.5) +
  geom_line(aes(y = ypred), color = "violet", linewidth = 1.5, linetype = "dashed") +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  facet_wrap(~ id) +
  theme_bw()

psi <- psi(saemix_fit_nlme) |>
  as_tibble() |>
  rowid_to_column("iteration")

psi |>
  ggplot(aes(x = iteration, y = y_min)) +
  geom_line(color = "tomato") +
  labs(
    x = "Iteration",
    y = "y_min"    
  ) +
  theme_bw()

psi |>
  ggplot(aes(x = iteration, y = y_max)) +
  geom_line(color = "tomato") +
  labs(
    x = "Iteration",
    y = "y_max"    
  ) +
  theme_bw()

psi |>
  ggplot(aes(x = iteration, y = rate)) +
  geom_line(color = "tomato") +
  labs(
    x = "Iteration",
    y = "rate"    
  ) +
  theme_bw()



#-----------------------------------------------#
#-----------------------------------------------#
#                Model Selection                #
#-----------------------------------------------#
#-----------------------------------------------#


# Define the function to fit and evaluate models
fit_model <- function(covariate_matrix) {
  # Define the saemix model with appropriate covariates
  saemix_model <- saemixModel(
    model = function_approx,
    description = "Antibody decay model with covariates",
    psi0 = psi0_matrix,
    transform.par = c(0, 0, 1),  # log-transform rate
    covariance.model = diag(1, 3),
    covariate.model = covariate_matrix,  # Covariate matrix based on input
    error.model = "constant"
  )
  
  # Fit the model
  saemix_fit <- saemix(
    saemix_model,
    saemix_data_nlme,
    saemix_options_nlme
  )
  
  # Get the summary of the model fit to extract AIC, BIC, logLik
  fit_summary <- summary(saemix_fit)
  
  # Extract AIC, BIC, and logLik from the summary
  aic <- fit_summary$AIC
  bic <- fit_summary$BIC
  logLik <- fit_summary$logLik
  
  return(list(aic = aic, bic = bic, logLik = logLik, fit = saemix_fit))
}

# Define all covariate configurations
covariate_configs <- list(
  # Single covariates
  list(age_centered = 1),  # Only age
  list(bmi = 1),           # Only BMI
  list(country = 1),       # Only country
  list(sex = 1),           # Only sex
  
  # Two covariates
  list(age_centered = 1, bmi = 1),                # Age + BMI
  list(age_centered = 1, country = 1),            # Age + Country
  list(age_centered = 1, sex = 1),                # Age + Sex
  list(bmi = 1, country = 1),                     # BMI + Country
  list(bmi = 1, sex = 1),                         # BMI + Sex
  list(country = 1, sex = 1),                     # Country + Sex
  
  # Three covariates
  list(age_centered = 1, bmi = 1, country = 1),   # Age + BMI + Country
  list(age_centered = 1, bmi = 1, sex = 1),       # Age + BMI + Sex
  list(age_centered = 1, country = 1, sex = 1),   # Age + Country + Sex
  list(bmi = 1, country = 1, sex = 1),            # BMI + Country + Sex
  
  # All four covariates
  list(age_centered = 1, bmi = 1, country = 1, sex = 1)  # All four
)

# Loop through all covariate configurations and fit the models
results <- lapply(covariate_configs, function(config) {
  covariate_matrix <- matrix(
    data = c(ifelse("age_centered" %in% names(config), 1, 0),
             ifelse("bmi" %in% names(config), 1, 0),
             ifelse("country" %in% names(config), 1, 0),
             ifelse("sex" %in% names(config), 1, 0)),
    nrow = 1,
    ncol = 3,
    dimnames = list(NULL, c("y_min", "y_max", "rate"))
  )
  # Fit the model and return fit criteria
  fit_model(covariate_matrix)
})

# Initialize an empty list to hold results
aic_values <- c()
bic_values <- c()
logLik_values <- c()

# Iterate over the results and extract AIC, BIC, and logLik (if available)
for (res in results) {
  aic_values <- c(aic_values, res$aic)
  bic_values <- c(bic_values, res$bic)
  logLik_values <- c(logLik_values, ifelse(is.null(res$logLik), NA, res$logLik))
}

# Convert to a tibble for comparison
model_comparison <- tibble(
  model = names(covariate_configs),
  aic = aic_values,
  bic = bic_values,
  logLik = logLik_values
)

# Print the model_comparison to check its structure
print(model_comparison)

# Sort models based on AIC (or BIC) to find the best model
model_comparison_sorted <- model_comparison %>%
  arrange(aic)  # Use `bic` for BIC-based selection, or AIC for AIC-based selection

# Display the sorted models
print(model_comparison_sorted)
