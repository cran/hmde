## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
# remotes::install_github("traitecoevo/hmde")
# install.packages(c("dplyr", "ggplot2"))

library(hmde)
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
hmde_model("vb_single_ind")
#prior_pars_ind_max_size_sd_only is the argument name for the prior parameter

## -----------------------------------------------------------------------------
#Analytic solution in function form
solution <- function(t, pars = list(y_0, beta, S_max)){
  return(
    pars$S_max + (y_0 - pars$S_max)*exp(-t * pars$beta)
  )
}

#Parameters
beta <- 0.35 #Growth rate
y_0 <- 1 #Starting size
S_max <- 20 #Asymptotic max size
time <- c(0,30) 
pars_list <- list(y_0 = y_0,
                  beta = beta,
                  S_max = S_max)
y_final <- solution(time[2], pars_list)

#Plot of growth function
ggplot() +
  xlim(y_0, y_final) +
  ylim(0, beta*(S_max-y_0)*1.1) +
  labs(x = "Y(t)", y = "f", title = "von Berralanffy growth") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")) +
  geom_function(fun=hmde_model_des("vb_single_ind"), 
                args=list(pars = list(S_max, beta)),
                colour="green4", linewidth=1,
                xlim=c(y_0, y_final))

#Size over time
ggplot() +
  geom_function(fun=solution, 
                args=list(pars = pars_list),
                colour="green4", linewidth=1,
                xlim=c(time)) +
  xlim(time) +
  ylim(0, y_final*1.05) +
  labs(x = "Time", y = "Y(t)", title = "von Bertalanffy growth") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"))


## -----------------------------------------------------------------------------
lizard_vb_fit <- hmde_model("vb_multi_ind") |>
  hmde_assign_data(data = Lizard_Size_Data)  |>
  hmde_run(chains = 4, cores = 1, iter = 2000)

lizard_estimates <- hmde_extract_estimates(fit = lizard_vb_fit,
                                           input_measurement_data = Lizard_Size_Data)

## -----------------------------------------------------------------------------
measurement_data_transformed <- lizard_estimates$measurement_data %>%
  group_by(ind_id) %>%
  mutate(
    delta_y_obs = y_obs - lag(y_obs),
    obs_interval = time - lag(time),
    obs_growth_rate = delta_y_obs/obs_interval,
    delta_y_est = y_hat - lag(y_hat),
    est_growth_rate = delta_y_est/obs_interval
  ) %>%
  ungroup()


#Distributions of estimated growth and size
hist(measurement_data_transformed$y_hat, 
     main = "Estimated size distribution",
     xlab = "Size (cm)")
hist(measurement_data_transformed$delta_y_est, 
     main = "Estimated growth increments",
     xlab = "Growth increment (cm)")
hist(measurement_data_transformed$est_growth_rate, 
     main = "Estimated annualised growth rate distribution",
     xlab = "Growth rate (cm/yr)")

#Quantitative R^2
cor(measurement_data_transformed$y_obs, measurement_data_transformed$y_hat)^2
r_sq_est <- cor(lizard_estimates$measurement_data$y_obs,
                lizard_estimates$measurement_data$y_hat)^2
r_sq <- paste0("R^2 = ", 
               signif(r_sq_est,
                      digits = 3))

obs_scatter <- ggplot(data = lizard_estimates$measurement_data, 
       aes(x = y_obs, y = y_hat)) +
  geom_point(shape = 16, size = 1, colour = "green4") +
  xlab("Y obs.") +
  ylab("Y est.") +
  geom_abline(slope = 1, linetype = "dashed") +
  annotate("text", x = 25, y = 18, 
           label = r_sq) +
  theme_classic()

#Plots of size over time for a sample of 5 individuals
obs_est_ind <- hmde_plot_obs_est_inds(n_ind_to_plot = 5,
                       measurement_data = lizard_estimates$measurement_data) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.2))

## -----------------------------------------------------------------------------
#1-dimensional parameter distributions
s_max_hist <- ggplot(lizard_estimates$individual_data, 
       aes(ind_max_size_mean)) +
  geom_histogram(bins = 10,
                 colour = "black",
                 fill = "lightblue") +
  labs(x="S_max estimate") +
  theme_classic()

beta_hist <- ggplot(lizard_estimates$individual_data, 
       aes(ind_growth_rate_mean)) +
  geom_histogram(bins = 10,
                 colour = "black",
                 fill = "lightblue") +
  labs(x="beta estimate") +
  theme_classic()

#2-dimensional parameter distribution
par_scatter <- ggplot(data = lizard_estimates$individual_data, 
       aes(x = ind_max_size_mean, y = ind_growth_rate_mean)) +
  geom_point(shape = 16, size = 1, colour = "green4") +
  xlab("Individual max sizes (mm)") +
  ylab("Individual growth rate parameters") +
  theme_classic()

#Correlation of parameters
cor(lizard_estimates$individual_data$ind_max_size_mean,
    lizard_estimates$individual_data$ind_growth_rate_mean,
    method = "spearman")

#Plot function pieces over estimated sizes.
de_pieces <- hmde_plot_de_pieces(lizard_estimates)

## -----------------------------------------------------------------------------
pars_CI_names <- c(
  "mean log max size",
  "mean max size in mm",
  "log max size standard deviation",
  "mean log growth par",
  "mean growth par mm/yr",
  "log growth par standard deviation"
)

#Vector that picks out which pars to be exponentiated
exp_vec <- c(FALSE, TRUE, FALSE, 
             FALSE, TRUE, FALSE)

#Print mean estimates and CIs
for(i in 1:nrow(lizard_estimates$population_data)){
  if(!exp_vec[i]){
    lizard_estimates$population_data$mean[i] 
    print(paste0("95% CI for ", 
                 pars_CI_names[i],
                 ": (",
                 lizard_estimates$population_data$CI_lower[i],
                 ", ",
                 lizard_estimates$population_data$CI_upper[i],
                 ")"))
  } else {
    exp(lizard_estimates$population_data$mean[i])
    print(paste0("95% CI for ",
                 pars_CI_names[i], 
                 ": (",
                 exp(lizard_estimates$population_data$CI_lower[i]),
                 ", ",
                 exp(lizard_estimates$population_data$CI_upper[i]),
                 ")"))
  }
}

