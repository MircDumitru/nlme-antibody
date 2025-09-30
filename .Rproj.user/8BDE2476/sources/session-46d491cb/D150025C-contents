library(tidyverse)  # tidyverse 
library(nlme)       # non-linear mixed effect models



## Exploratory data analysis

# read the csv file
tb <- read_csv("data_r_codes/data_ebola_ab.csv")


# checking the var types
tb |>
  summarise(across(everything(), typeof)) |>
  knitr::kable(align = "ccccccc")

# checking for possible missing vals
tb |>
  summarise(across(everything(), \(x) sum(is.na(x)))) |>
  knitr::kable(align = "ccccccc")


tb <- tb |>
  # clean the names
  janitor::clean_names() |>
  # transform chr into factors
  mutate(across(where(is.character), fct))

# check the transformed tibble
# tb |>
#   glimpse()

# tb |>
#   summarise(across(where(is.factor), \(x) length(unique(x))))

tb |>
  summarise(across(where(is.factor), n_distinct)) |>
  knitr::kable(align = "ccc")

## The factor variables

### The `id` variable

tb_id_levels_number <- tb |>
  summarise(
    n = n(),
    .by = id
  ) |>
  mutate(
    prop = n / sum(n),
  )

tb_id_levels_number |>
  head() |>
  knitr::kable(align = "ccc")


tb_id_levels_number |>
  ggplot(aes(x = id, y = n)) +
  geom_point(size = 3, color = "tomato") +
  geom_segment(aes(x = id, xend = id, y = 0, yend = n), color = "tomato") +
  labs(
    title = "Number of observations of the ID levels",
    subtitle = "Each subject appears with 8 measurements in the data",
    x = "Subject's ID",
    y = "Number of observations"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))


tb_id_levels_number |>
  select(n) |>
  n_distinct()

n_obs_per_subject <- tb_id_levels_number |>
  pull(n) |>
  unique()
# n_obs_per_subject


tb_sex_levels_number <- tb |>
  summarise(
    n_obs = n(),
    .by = sex
  ) |>
  mutate(
    n_subjects = n_obs / n_obs_per_subject,
    prop = n_subjects / sum(n_subjects),
  )

tb_sex_levels_number |>
  knitr::kable(align = "cccc")


tb |>
  ggplot(aes(x = sex, y = after_stat(prop), group = 1)) +
  geom_bar(fill = "tomato") +
  labs(
    title = "Proportion of the sex levels in the data",
    x = "Subject's sex",
    y = "Proportion"
  ) +
  scale_y_continuous(labels = scales::percent)


### The `country` variable 

tb_country_levels_number <- tb |>
  summarise(
    n_obs = n(),
    .by = country
  ) |>
  mutate(
    n_subjects = n_obs / n_obs_per_subject,
    prop = n_subjects / sum(n_subjects),
  )

tb_country_levels_number |>
  knitr::kable(align = "ccc")


tb |>
  ggplot(aes(x = country,  y = after_stat(prop), group = 1)) +
  geom_bar(fill = "tomato") +
  labs(
    title = "Proportion of the country levels in the data",
    x = "Subjects's country",
    y = "Proportion"
  ) +
  scale_y_continuous(labels = scales::percent)


tb |>
  summarise(
    n_subjects = n() / n_obs_per_subject,
    .by = c(country, sex)
  ) |>
  arrange(country, sex) |>
  group_by(country) |>
  mutate(prop_across_regions = n_subjects / sum(n_subjects)) |>
  ungroup() |>
  mutate(prop = n_subjects / sum(n_subjects)) |>
  knitr::kable(align = "cccccc")


tb |>
  ggplot(aes(x = country, fill = sex)) +
  geom_bar(position = "fill") +
  labs(
    title = "Proportion of sex for each region",
    x = "Patient's country",
    y = "Counts", 
    fill = NULL
  ) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
  ) +
  scale_y_continuous(labels = scales::percent)


tb |>
  ggplot(aes(x = id, y = log(time))) +
  geom_point(color = "tomato") +
  labs(
    title = "Measurment times",
    x = "Patient's ID",
    y = "Log Time (days)", 
  ) +
  scale_y_continuous(breaks = seq(-2, 7, by = 2)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))


tb |>
  select(id, time) |>
  # pivot the tibble wider and get the time measurements 
  pivot_wider(
    names_from = time,
    names_prefix = "t_",
    values_from = time
    ) |>
  summarise(
    across(starts_with("t_"), n_distinct)
  )


tb |>
  ggplot(aes(x = age)) +
  geom_density(color = "tomato", fill = "tomato", alpha = .5) +
  labs(
    title = "Age density",
    x = "Age (years)",
    y = "Density"
  )


tb_negative_age <- tb |>
  filter(age < 0) |>
  distinct(id, age)

tb_negative_age

tb_clean <- tb |>
  filter(age >= 0)


tb_clean |>
  ggplot(aes(x = age, color = sex, fill = sex)) +
  geom_density(alpha = .5) +
  labs(
    title = "Age densities for sex levels",
    x = "Age (year)",
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom"
  )

tb_clean |>
  ggplot(aes(x = age, color = country, fill = country)) +
  geom_density(alpha = .5) +
  labs(
    title = "Age densities for region levels",
    x = "Age (years)",
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom"
  )


tb_clean |>
  summarise(
    min = min(age),
    median = median(age),
    IQR = IQR(age), 
    max = max(age),
    .by = sex
  ) |>
  knitr::kable(align = "ccccc")


tb_clean |>
  summarise(
    min = min(age),    
    median = median(age),
    IQR = IQR(age), 
    max = max(age),
    .by = country
  )  |>
  knitr::kable(align = "ccccc")


tb_clean |>
  ggplot(aes(y = age, x = sex, fill = sex)) +
  geom_boxplot(alpha = .5, show.legend = FALSE) +
  labs(
    title = "Age for sex levels",
    x = "Sex",
    y = "Age (years)"
  ) 

tb_clean |>
  ggplot(aes(y = age, x = country, fill = country)) +
  geom_boxplot(alpha = .5, show.legend = FALSE) +
  labs(
    title = "Age for region levels",
    x = "Region",
    y = "Age (years)"
  ) 


tb_clean |>
  nrow() / n_obs_per_subject == tb_clean |>
  distinct(id, age) |>
  nrow()


tb_clean <- tb_clean |>
  mutate(age_updated = age + time / 365) |>
  relocate(age_updated, .after = age)

tb_clean


tb_clean |>
  nrow() / n_obs_per_subject == tb_clean |>
  distinct(id, bmi) |>
  nrow()


tb_clean |>
  ggplot(aes(x = bmi)) +
  geom_density(color = "tomato", fill = "tomato", alpha = .5) +
  labs(
    title = "BMI density",
    x = "BMI (kg/m²)",
    y = "Density"
  )


tb_clean |>
  ggplot(aes(x = bmi, color = sex, fill = sex)) +
  geom_density(alpha = .5) +
  labs(
    title = "BMI densities for sex levels",
    x = "BMI (kg/m²)",
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom"
  )

tb_clean |>
  ggplot(aes(x = bmi, color = country, fill = country)) +
  geom_density(alpha = .5) +
  labs(
    title = "BMI densities for region levels",
    x = "BMI (kg/m²)",
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom"
  )


tb_clean |>
  ggplot(aes(y = bmi, x = sex, fill = sex)) +
  geom_boxplot(alpha = .5, show.legend = FALSE) +
  labs(
    title = "BMI for sex levels",
    x = "Sex",
    y = "BMI (kg/m²)"
  ) 

tb_clean |>
  ggplot(aes(y = bmi, x = country, fill = country)) +
  geom_boxplot(alpha = .5, show.legend = FALSE) +
  labs(
    title = "BMI for region levels",
    x = "Region",
    y = "BMI (kg/m²)"
  ) 


tb_clean |>
  ggplot(aes(x = abelisa)) +
  geom_density(color = "tomato", fill = "tomato", alpha = .5) +
  labs(
    title = "Antibody concentration density",
    x = "Antibody concentration (EU/mL)",
    y = "Density"
  ) 


tb_clean |>
  ggplot(aes(x = abelisa, color = sex, fill = sex)) +
  geom_density(alpha = .5) +
  labs(
    title = "Antibody conc. densities for sex levels",
    x = "Antibody conc. (EU/mL)",
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom"
  ) 

tb_clean |>
  ggplot(aes(x = abelisa, color = country, fill = country)) +
  geom_density(alpha = .5) +
  labs(
    title = "Antibody conc. densities for region levels",
    x = "Antibody conc. (EU/mL)",
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom"
  )


tb_clean |>
  ggplot(aes(y = abelisa, x = sex, fill = sex)) +
  geom_boxplot(alpha = .5, show.legend = FALSE) +
  labs(
    title = "Antibody concentration for sex levels",
    x = "Sex",
    y = "Antibody conc. (EU/mL)"
  ) 

tb_clean |>
  ggplot(aes(y = abelisa, x = country, fill = country)) +
  geom_boxplot(alpha = .5, show.legend = FALSE) +
  labs(
    title = "Antibody concentration for region levels",
    x = "Region",
    y = "Antibody conc. (EU/mL)"
  ) 


tb_clean |>
  ggplot(aes(x = time, y = abelisa, color = bmi, group = id)) +
  geom_line() +
  scale_color_distiller(palette = "RdPu") +
  labs(
    title = "Antibody concentration during time",
    x = "Time (days)",
    y = "Antibody conc. (EU/mL)",
    color = "BMI"
  ) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom"
  )


tb_clean |>
  ggplot(aes(x = time, y = abelisa, color = bmi, group = id)) +
  geom_line() +
  facet_grid(sex ~ country) +
  scale_color_distiller(palette = "RdPu") +
  labs(
    title = "Antibody concentration during time",
    x = "Time (days)",
    y = "Antibody conc. (EU/mL)",
    color = "BMI"
  ) +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom"
  )

## Non-linear model with mixed effects.

tb_clean |>
  filter(id == "ID1") |>
  ggplot(aes(x = time, y = abelisa)) +
  geom_point(color = "tomato") +
  geom_line(color = "tomato") +
  labs(
    title = "Antibody concentration during time for subject ID1",
    x = "Time (days)",
    y = "Antibody conc. (EU/mL)"
  ) 


compartment_model <- function(psi, t){
  k_a <- psi[1] 
  k_e <- psi[2]
  v  <- psi[3]
  return (k_a / (v * (k_a - k_e)) * (exp(-k_e * t) - exp(-k_a * t)))
}


# Parameters 
psi = c(10, 0.01, 50)
# Generate time 
t <- seq(0, 730, by = 1)

compartment_data <- compartment_model(psi, t)
tb_compartment <- tibble(
  t = t,
  f = compartment_data
)
tb_compartment |>
  ggplot(aes(x = t, y = f)) +
  geom_line()

# nls_one_sib <- nls(
#   abelisa ~ compartment_model(psi, time),
#   data = tb_clean,
#   start = list(psi=c(k_a = 1, k_e = 0.1, v = 20))
#   )
# nls_one_sib




# Parameters for Exponential Decay
# A_exp <- 100000  # Initial concentration
# k_exp <- 0.002 # Decay rate


# Generate time points
# time <- seq(0, 730, by = 1)
# concentration_exp <- exp_decay(time, A_exp, k_exp)


# Plot Exponential Decay
# plot(time, concentration_exp, type = "l", col = "blue", lwd = 2, 
     # xlab = "Time (Days)", ylab = "Antibody Concentration", 
     # main = "Exponential Decay Model")
# grid()








#----------------------
#----------------------

# install.packages("nlme")
library(tidyverse)  # tidyverse 
library(nlme)       # non-linear mixed effect models
library(saemix)


tb_theo <- tibble(Theoph)

tb_theo <- tb_theo |>
  janitor::clean_names() |>
  rename(
    id = subject,
    weight = wt,
    concentration = conc
  ) |>
  relocate(time, .after = id) |>
  relocate(concentration, .after = time) |>
  select(-dose) |>
  filter(time != 0) |>
  mutate(
    id = fct_relevel(id, as.character(1:12))
  )

tb_theo 

tb_theo |>
  ggplot(aes(x = time, y = concentration, group = id)) +
  geom_point(color = "tomato") +
  geom_line(color = "tomato")

tb_theo |>
  ggplot(aes(x = time, y = concentration, group = id)) +
  geom_point(color = "tomato") +
  geom_line(color = "tomato") +
  facet_wrap( ~ id)

#----------------------
#----------------------

subject_1 <- tb_theo |>
  filter(id == 1) 

subject_1 |>
  ggplot(aes(x = time, y = concentration)) +
  geom_point(color = "tomato") +
  geom_line(color = "tomato")

f1 <- function(psi, t){
  k_a <- psi[1]
  k_e <- psi[2]
  v  <- psi[3]
  f  <- k_a / (v * (k_a - k_e)) * (exp(-k_e * t) - exp(-k_a * t)) 
  f
}

model_subject_1 <- nls(
  concentration ~ f1(psi, time), 
  start = list(psi = c(k_a = 1, k_e = 0.1, v = 0.125)), 
  data = subject_1
)
model_subject_1
coef(model_subject_1)


subject_1_pred <- tibble(time = seq(0, 40, 0.2)) 
subject_1_pred <- subject_1_pred |>
  mutate(
    concentration_pred = predict(
      model_subject_1, 
      newdata = subject_1_pred
    )
  )
subject_1_pred


subject_1 |>
  ggplot() +
  geom_point(
    aes(x = time, y = concentration),
    color = "tomato",
    size = 2
  ) +
  geom_line(
    data = subject_1_pred, 
    aes(x = time, y = concentration_pred), 
    color = "darkseagreen",
    linewidth = 1
  )

#----------------------
#----------------------

model_all <- nls(
  concentration ~ f1(psi, time), 
  start = list(psi = c(k_a = 1, k_e = 0.1, v = 0.125)), 
  data = tb_theo
)
model_all
coef(model_all)


all_pred <- tibble(time = seq(0, 40, 0.2)) 
all_pred <- all_pred|>
  mutate(
    concentration_pred = predict(
      model_all, 
      newdata = all_pred
    )
  )
all_pred


tb_theo |>
  ggplot() +
  geom_point(
    aes(x = time, y = concentration, group = id),
    color = "tomato",
    size = 2
  ) +
  geom_line(
    data = all_pred, 
    aes(x = time, y = concentration_pred), 
    color = "darkseagreen",
    linewidth = 1
  )

tb_theo |>
  ggplot() +
  geom_point(
    aes(x = time, y = concentration, group = id),
    color = "tomato",
    size = 2
  ) +
  geom_line(
    data = all_pred, 
    aes(x = time, y = concentration_pred), 
    color = "darkseagreen",
    linewidth = 1
  ) +
  facet_wrap( ~ id)

#----------------------
#----------------------

subject_individual <- tibble()

for(i in 1: 12){
#subject_i <- tb_theo |>
#  filter(id == i) 

model_subject_i <- nls(
  concentration ~ f1(psi, time), 
  start = list(psi = c(k_a = 1, k_e = 0.1, v = 0.125)), 
  data = tb_theo |>
    filter(id == i) 
)

subject_i_pred <- tibble(
  time = seq(0, 40, 0.2),
  id = rep(i, length(time))
  ) 
subject_i_pred <- subject_i_pred |>
  mutate(
    concentration_pred = predict(
      model_subject_i, 
      newdata = subject_i_pred
    )
  )

subject_individual <- bind_rows(subject_individual, subject_i_pred)
}

subject_individual <- subject_individual |>
  relocate(id, .before = time) |>
  mutate(
    id = factor(id, levels = as.character(1:12), ordered = TRUE)
  ) 

subject_individual

tb_theo |>
  ggplot() +
  geom_point(
    aes(x = time, y = concentration, group = id),
    color = "tomato",
    size = 2
  ) +
  geom_line(
    data = subject_individual, 
    aes(x = time, y = concentration_pred, group = id), 
    color = "darkseagreen",
    linewidth = 1
  ) +
  facet_wrap( ~ id)

model_subject_9 <- nls(
  concentration ~ f1(psi, time), 
  start = list(psi = c(k_a = 1, k_e = 0.1, v = 0.125)), 
  data = tb_theo |>
    filter(id == 9) 
)
model_subject_9

#----------------------
#----------------------
saemix_data <- saemixData(name.data = tb_theo,
                          name.group = "id",
                          name.predictors = "time",
                          name.response = "concentration")

model_1_nlme <- function(psi, id, x){
  t   <- x[, 1]
  k_a <- psi[id, 1]
  k_e <- psi[id, 2]
  v  <- psi[id, 3]
  fpred  <- k_a / (v * (k_a - k_e)) * (exp(-k_e * t) - exp(-k_a * t)) 
  fpred
}

saemix_model <- saemixModel(
  model = model_1_nlme,
  psi0  = c(k_a = 1, k_e = 0.1, v = 0.125)
  )

saemix_options <- list(
  map = TRUE, 
  fim = TRUE, 
  ll.is = FALSE, 
  displayProgress = FALSE, 
  save = FALSE, 
  seed = 42
  )

saemix_fit1 <- saemix(
  saemix_model, 
  saemix_data, 
  saemix_options
  )

saemix_fit1@results

psi <- psi(saemix_fit1)
psi

saemix_fit <- saemix.predict(saemix_fit1)

saemix.plot.fits(saemix_fit1)

saemix.plot.obsvspred(saemix_fit1,level=1)

saemix.plot.scatterresiduals(saemix_fit1, level=1)












