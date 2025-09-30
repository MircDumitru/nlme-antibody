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


tb_clean_centered <- tb_clean |>
  mutate(age_centered = age - mean(age)) |>
  relocate(age_centered, .after = age) |>
  select(-c(age, age_updated))

tb_after_7_centered <- tb_after_7 |>
  mutate(age_centered = age - mean(age)) |>
  relocate(age_centered, .after = age) |>
  select(-c(age, age_updated))


# consntat corresponding to intial values
c_const <- 100

# Nonlinear function with log-parameters and fixed c
f_model_log <- function(time, log_phi_S, log_phi_L, 
                        log_delta_A, log_delta_S, log_delta_L) {
  phi_S <- exp(log_phi_S)
  phi_L <- exp(log_phi_L)
  delta_A <- exp(log_delta_A)
  delta_S <- exp(log_delta_S)
  delta_L <- exp(log_delta_L)
  
  c_const * exp(-delta_A * time) +
    (phi_S / (delta_A - delta_S)) * exp(-delta_S * time) +
    (phi_L / (delta_A - delta_L)) * exp(-delta_L * time)
}

# Starting values (log-scale for parameters)
start_vals <- c(
  log_phi_S = log(50),
  log_phi_L = log(30),
  log_delta_A = log(0.0005),
  log_delta_S = log(0.005),
  log_delta_L = log(0.001)
)

# Fit nlme model:
nlme_model <- nlme(
  abelisa ~ f_model_log(time, log_phi_S, log_phi_L, 
                        log_delta_A, log_delta_S, log_delta_L),
  
  # Data: 
  data = tb_after_7_centered,
  
  # Fixed effects:
  fixed = list(
    log_phi_S ~ sex + country,
    log_phi_L ~ sex + country,
    log_delta_A ~ sex + country,
    log_delta_S + log_delta_L ~ 1  # no covariates on delta S and L params
  ),
  
  # Random effects on log_phi_S, log_phi_L, per subject
  random = log_phi_S + log_phi_L ~ 1| id,
  # Starting values
  start = start_vals,
  
  control = nlmeControl(
    pnlsTol = 0.1, 
    msMaxIter = 100, 
    msMaxEval = 1000
  )
)

summary(nlme_model)










nlme_model <- nlme(
  abelisa ~ f_model_log(time, log_phi_S, log_phi_L, 
                        log_delta_A, log_delta_S, log_delta_L),
  
  data = tb_after_7_centered,
  
  fixed = log_phi_S + log_phi_L + log_delta_A + log_delta_S + log_delta_L ~ 1,
  
  random = ~ log_phi_S + log_phi_L | id,
  
  start = start_vals,
  
  control = nlmeControl(
    pnlsTol = 0.1,
    msMaxIter = 100,
    msMaxEval = 1000
  )
)

