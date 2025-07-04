## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = 'center',
  message = FALSE, 
  warning = FALSE
)

library(patchwork)

options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
# remotes::install_github("traitecoevo/hmde")
# install.packages(c("dplyr", "ggplot2"))

library(hmde)
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
hmde_model("constant_single_ind")
#prior_pars_ind_beta is the argument name for the prior parameters

## -----------------------------------------------------------------------------
beta <- 2 #Annual growth rate
y_0 <- 1 #Starting size
time <- 0:20 
sizes_over_time <- tibble(Y_t = 1 + beta*time, #Linear sizes over time
                          t = time)

sizes_over_time

## ----echo=FALSE, fig.width=6, fig.height=3------------------------------------
#Plot of growth function
constant_growth_function <- ggplot() +
  xlim(y_0, max(sizes_over_time$Y_t)) +
  ylim(0, beta*2) +
  labs(x = "Y(t)", y = "f", title = "Constant growth") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")) +
  geom_function(fun=hmde_model_des("constant_single_ind"), 
                args=list(pars = list(beta)),
                colour="green4", linewidth=1,
                xlim=c(y_0, max(sizes_over_time)))

#Sizes over time
sizes_over_time_plot <- ggplot(data = sizes_over_time, aes(x=t, y = Y_t)) +
  geom_line(colour="green4", linewidth=1) +
  xlim(0, max(sizes_over_time$t)) +
  ylim(0, max(sizes_over_time$Y_t)*1.05) +
  labs(x = "Time", y = "Y(t)", title = "Size over time") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"))

constant_growth_function + sizes_over_time_plot

## ----eval=FALSE, fold---------------------------------------------------------
# #Plot of growth function
# ggplot() +
#   geom_function(fun = hmde_model_des("constant_single_ind"), # Visualising the growth function
#                 args = list(pars = list(beta)),
#                 colour = "green4", linewidth=1,
#                 xlim = c(y_0, max(sizes_over_time))) +
#   xlim(y_0, max(sizes_over_time$Y_t)) + # Creating the x axis
#   ylim(0, beta*2) + # Creating the y axis
#   labs(x = "Y(t)", y = "f", title = "Constant growth") + # Axe labels and plot title
#   theme_classic() + # Theme for the plot
#   theme(axis.text=element_text(size=16), # Plot customising
#         axis.title=element_text(size=18,face="bold"))
# 
# #Sizes over time
# ggplot(data = sizes_over_time, aes(x=t, y = Y_t)) +
#   geom_line(colour="green4", linewidth=1) + # Line graph of sizes_over_time
#   xlim(0, max(sizes_over_time$t)) +
#   ylim(0, max(sizes_over_time$Y_t)*1.05) +
#   labs(x = "Time", y = "Y(t)", title = "Constant growth") +
#   theme_classic() +
#   theme(axis.text=element_text(size=16),
#         axis.title=element_text(size=18,face="bold"))

## -----------------------------------------------------------------------------
Trout_Size_Data

## -----------------------------------------------------------------------------
Trout_Size_Data_transformed <- Trout_Size_Data %>%
  group_by(ind_id) %>%
  mutate(
    delta_y_obs = y_obs - lag(y_obs),
    obs_interval = time - lag(time),
    obs_growth_rate = delta_y_obs/obs_interval
  ) %>%
  ungroup()

Trout_Size_Data_transformed

## ----echo = FALSE, fig.width=6, fig.height=6----------------------------------
histogram_func <- function(data, var, main, xlab, ...){
  ggplot2::ggplot(data = data, aes(x = {{var}})) + 
  geom_histogram(colour = "black", fill = "lightblue", ...) + 
  labs(title = main,
       x= xlab) + 
    theme_classic() 
}

hist_a <- histogram_func(Trout_Size_Data, 
                         y_obs, 
                         "Observed size distribution",
                         "Size (cm)",
                         binwidth = 5)

hist_b <- histogram_func(Trout_Size_Data_transformed, 
                         obs_interval, 
                         "Observed interval distribution",
                         "Time (yr)", 
                         binwidth = 0.55)

hist_c <- histogram_func(Trout_Size_Data_transformed,
                         delta_y_obs, 
                         "Observed growth increments",
                         "Growth increment (cm)", 
                         binwidth = 5.5)

hist_d <- histogram_func(Trout_Size_Data_transformed, 
                         obs_growth_rate, 
                         "Observed annualised growth \n rate distribution", 
                         "Growth rate (cm/yr)", 
                         binwidth = 80)

hist_a + hist_b + hist_c + hist_d + plot_layout(ncol = 2)

## -----------------------------------------------------------------------------
hmde_model("constant_multi_ind") |> 
  names()

## -----------------------------------------------------------------------------
trout_constant_fit <- hmde_model("constant_multi_ind") |>
  hmde_assign_data(data = Trout_Size_Data)  |>
  hmde_run(chains = 4, cores = 1, iter = 2000)

## -----------------------------------------------------------------------------
trout_constant_estimates <- hmde_extract_estimates(                                 
                                 fit = trout_constant_fit,
                                 input_measurement_data = Trout_Size_Data)

## -----------------------------------------------------------------------------
measurement_data_transformed <- trout_constant_estimates$measurement_data %>%
  group_by(ind_id) %>%
  mutate(
    delta_y_obs = y_obs - lag(y_obs),
    obs_interval = time - lag(time),
    obs_growth_rate = delta_y_obs/obs_interval,
    delta_y_est = y_hat - lag(y_hat),
    est_growth_rate = delta_y_est/obs_interval
  ) %>%
  ungroup()

## ----fig.width=5, fig.height=6------------------------------------------------
est_hist_y_hat <- histogram_func(measurement_data_transformed, y_hat, 
               "Estimated size distribution",
               xlab = "Size (cm)",
               bins = 5)

est_hist_delta_y_est <-  histogram_func(measurement_data_transformed, delta_y_est, 
               "Estimated growth  \n increments",
               xlab = "Growth increment (cm)",
               bins = 5)

est_hist_growth_rate <- histogram_func(measurement_data_transformed, est_growth_rate, 
               "Estimated annualised growth rate distribution", xlab = "Growth rate (cm/yr)",
               bins = 5)

(est_hist_y_hat + est_hist_delta_y_est) / est_hist_growth_rate 

## ----fig.width=4, fig.height=4------------------------------------------------
#Quantitative R^2
cor(measurement_data_transformed$y_obs, measurement_data_transformed$y_hat)^2
r_sq_est <- cor(trout_constant_estimates$measurement_data$y_obs,
                          trout_constant_estimates$measurement_data$y_hat)^2
r_sq <- paste0("R^2 = ", 
               signif(r_sq_est,
                      digits = 3))

obs_est_size_plot <- ggplot(data = trout_constant_estimates$measurement_data, 
       aes(x = y_obs, y = y_hat)) +
  geom_point(shape = 16, size = 1, colour = "green4") +
  xlab("Y obs.") +
  ylab("Y est.") +
  geom_abline(slope = 1, linetype = "dashed") +
  annotate("text", x = 45, y = 80, 
           label = r_sq) +
  theme_classic()
obs_est_size_plot

#Plots of size over time for a sample of 5 individuals
size_over_time_plot <- hmde_plot_obs_est_inds(n_ind_to_plot = 5,
                       measurement_data = trout_constant_estimates$measurement_data)
size_over_time_plot


## ----echo=FALSE, fig.width=6.5, fig.height=5----------------------------------
ind_hist_beta <- histogram_func(trout_constant_estimates$individual_data,
               ind_beta_mean,
               main = "Individual beta parameters", 
               xlab = "beta estimate")

de_pieces <- hmde_plot_de_pieces(trout_constant_estimates)


ind_hist_beta + de_pieces 

## -----------------------------------------------------------------------------
#Mean of normal distribution
trout_constant_estimates$population_data$mean[1] #Raw value
print(paste0("95% CI for mean log growth: (", 
             trout_constant_estimates$population_data$CI_lower[1], " , ",
             trout_constant_estimates$population_data$CI_upper[1], ")")) #Raw CI

exp(trout_constant_estimates$population_data$mean[1]) #In cm/yr units
print(paste0("95% CI for mean growth in cm/yr: (", 
             exp(trout_constant_estimates$population_data$CI_lower[1]), " , ",
             exp(trout_constant_estimates$population_data$CI_upper[1]), ")"))

#Standard deviation of underlying normal distribution
trout_constant_estimates$population_data$mean[2]
print(paste0("95% CI for log growth standard deviation: (", 
             trout_constant_estimates$population_data$CI_lower[2], " , ",
             trout_constant_estimates$population_data$CI_upper[2], ")")) #Raw CI

