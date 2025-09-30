# Load libraries
library(saemix)
library(dplyr)
library(ggplot2)


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


# -------------------------
# Step 1: Prepare the data
# -------------------------

# Center the age covariate (if not already done)
tb_after_7 <- tb_after_7 %>%
  mutate(age_centered = age - mean(age, na.rm = TRUE))

# Create saemixData object (for saemix 3.4+)
saemix_data_nlme <- saemixData(
  name.data = tb_after_7,
  name.group = "id",
  name.predictors = c("time", "age_centered"),  # time and covariate
  name.response = "abelisa",
  units = list(x = "Days", y = "EU/ml"),
  verbose = FALSE
)

# -----------------------------------------------
# Step 2: Define initial values and covariate model
# -----------------------------------------------

# psi0: two rows (intercept, covariate effect), 3 columns (parameters)
psi0_matrix <- matrix(
  c(
    1.6, 4.9, 0.05,  # Intercepts (base values for y_min, y_max, rate)
    0.0, 0.0, 0.1    # Covariate effects (only rate affected)
  ),
  nrow = 2,
  byrow = TRUE
)
colnames(psi0_matrix) <- c("y_min", "y_max", "rate")

# covariate.model: 1 row (covariate), 3 columns (parameters)
covariate_model_matrix <- matrix(
  c(0, 0, 1),
  nrow = 1,
  ncol = 3,
  byrow = TRUE,
  dimnames = list(
    c("age_centered"),
    c("y_min", "y_max", "rate")
  )
)
covariate_model_matrix <- as(covariate_model_matrix, "matrix")

# --------------------------------------
# Step 3: Define the structural model
# --------------------------------------

function_approx <- function(psi, id, xidep) {
  time <- xidep[, "time"]
  y_min <- psi[id, "y_min"]
  y_max <- psi[id, "y_max"]
  rate  <- psi[id, "rate"]
  
  y_min + (y_max - y_min) * exp(-rate * time)
}

# --------------------------------------
# Step 4: Build the saemixModel
# --------------------------------------

saemix_model_nlme <- saemixModel(
  model = function_approx,
  description = "Antibody decay with age effect on rate",
  psi0 = psi0_matrix,
  transform.par = c(0, 0, 1),         # log-transform rate
  fixed.estim = c(1, 1, 1),           # estimate all parameters
  covariance.model = diag(1, 3),      # random effects on all
  covariate.model = covariate_model_matrix,
  error.model = "constant"
)

# --------------------------------------
# Step 5: Set estimation options
# --------------------------------------

saemix_options_nlme <- list(
  seed = 42,
  save = FALSE,
  save.graph = FALSE,
  displayProgress = TRUE,
  print = TRUE
)

# --------------------------------------
# Step 6: Fit the model
# --------------------------------------

saemix_fit_nlme <- saemix(
  saemix_model_nlme,
  saemix_data_nlme,
  saemix_options_nlme
)

# --------------------------------------
# Step 7: Inspect results
# --------------------------------------

print(saemix_fit_nlme)
print(saemix_fit_nlme@results@fixed.effects)

# --------------------------------------
# Step 8: Plot observed vs predicted (example: subject 1)
# --------------------------------------

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
    title = "Subject 1: Observed vs Predicted",
    x = "Time (days)",
    y = "Antibody levels (EU/ml)"
  ) +
  theme_bw() +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom"
  )
