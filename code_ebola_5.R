#----------------------
#----------------------

library(tidyverse)  # tidyverse 
library(nlme)       # non-linear mixed effect models
library(saemix)

#----------------------
#----------------------

theophylline <- tibble(Theoph)

theophylline <- theophylline |>
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

theophylline 

theo_plot <- theophylline %>% 
  ggplot() + aes(x = time, y = concentration) + geom_point(color="#993399", size=2) +
  xlab("time (h)") + ylab("concentration (mg/l)")

theo_plot + geom_line(color="#993399", aes(group = id))

theo_plot + geom_line() + facet_wrap( ~ id)

#----------------------
#----------------------

subject1 <- theophylline %>% 
  filter(id == 1) %>% select("time","concentration")

subject1_plot <- subject1 %>% 
  ggplot() + aes(x = time, y = concentration) + geom_point( color="#993399", size=3) + 
  xlab("time (h)") + ylab("concentration (mg/l)") + ylim(c(0,11))

subject1_plot + geom_line(color="#993399")

f1 <- function(psi, t){
  D  <- 320; ka <- psi[1]; V  <- psi[2]; ke <- psi[3]
  f  <- D*ka/V/(ka-ke)*(exp(-ke*t)-exp(-ka*t)) 
  f
}

model_1 <- nls(concentration ~ f1(psi, time), start = list(psi=c(ka=1, V=40, ke=0.1)), data=subject1)

coef(model_1)

dplot <- data.frame(time = seq(0, 40, by=0.2))

dplot$pred_1 <- predict(model_1, newdata = dplot)

subject1_plot + geom_line(data = dplot, aes(x = time, y = pred_1), colour = "#339900", linewidth=1)
#----------------------
#----------------------

model_all <- nls(concentration ~ f1(psi, time), start = list(psi=c(ka=1, V=40, ke=0.1)), data=theophylline)
coef(model_all)

dplot$pred_all <- predict(model_all, newdata = dplot)
theo_plot + geom_line(data = dplot, aes(x = time, y = pred_all), colour="#339900", size=1)

theo_plot +  
  geom_line(data = dplot, aes(x = time, y=pred_all), colour="#339900", linewidth=1) + 
  facet_wrap(~ id)

#----------------------
#----------------------

res <- split(theophylline, theophylline$id) %>% 
  map(~{
    model_i <- nls(concentration ~ f1(psi, time), 
                   start = list(psi=c(ka=1, V=40,k=0.08)), 
                   data = .x)
    list(psi = coef(model_i),
         y_hat = predict(model_i, newdata = dplot),
         id = unique(.x$id))
  })
psi_hat <- map_df(res, "psi") %>% 
  setNames(c("ka","V","ke")) %>% 
  add_column(id = factor(map_dbl(res, "id"))) 

theo_pred <-
  map_df(res, "y_hat") %>% 
  pivot_longer(everything(), names_to = "id", values_to = "concentration") %>% 
  add_column(time = rep(dplot$time, each = length(res)))

theo_plot + geom_line(data = theo_pred, aes(x=time,y=concentration), colour="#339900", size=0.75) + facet_wrap(~id)

model_9 <- nls(concentration ~ f1(psi, time),  start = list(psi=c(ka=1, V=40,k=0.08)), 
               data = filter(theophylline, id == 9))
model_9

#----------------------
#----------------------

saemix_data <- saemixData(name.data       = theophylline,
                          name.group      = "id",
                          name.predictors = "time",
                          name.response   = "concentration")

model1_nlme <- function(psi,id,x) {
  D   <- 320
  t   <- x[,1]
  ka  <- psi[id,1]
  V   <- psi[id,2]
  ke  <- psi[id,3]
  fpred <- D*ka/(V*(ka-ke))*(exp(-ke*t)-exp(-ka*t))
  fpred
}

saemix_model <- saemixModel(model = model1_nlme,
                            psi0  = c(ka=1,V=20,ke=0.5))

saemix_options <- list(map=TRUE, fim=TRUE, ll.is=FALSE, displayProgress=FALSE, save=FALSE, seed=632545)

saemix_fit1    <- saemix(saemix_model, saemix_data, saemix_options)

saemix_options <- list(map=TRUE, fim=TRUE, ll.is=FALSE, displayProgress=FALSE, save=FALSE, seed=632545)

saemix_fit1    <- saemix(saemix_model, saemix_data, saemix_options)

saemix_fit1@results
