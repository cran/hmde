## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE
)

options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
library(hmde)
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
g_max <- 1 #Max growth rate
S_max <- 10 #Size at which the maximum growth occurs
k <- 0.75
y_0 <- 1 #Starting size
y_final <- 40

#Plot of growth function
ggplot() +
  xlim(y_0, y_final) +
  labs(x = "Y(t)", y = "f", title = "Canham growth") +
  theme_classic() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold")) +
  geom_function(fun=hmde_model_des("canham_single_ind"),
                args=list(pars = list(g_max, S_max, k)),
                colour="green4", linewidth=1,
                xlim=c(y_0, y_final))

## -----------------------------------------------------------------------------
hmde_model("canham_single_ind")
#prior_pars_ind_max_growth is the argument name for the prior parameters

## -----------------------------------------------------------------------------
hist(Tree_Size_Data$y_obs,
     xlab = "Size (cm)", main ="")

Tree_Size_Data_Transformed <- Tree_Size_Data %>%
  group_by(ind_id) %>%
  mutate(Delta_y_obs = y_obs - lag(y_obs)) %>%
  ungroup() %>%
  arrange(ind_id, time) %>%
  filter(!is.na(Delta_y_obs))

hist(Tree_Size_Data_Transformed$Delta_y_obs,
     xlab = "Growth increment (cm)", main="")

## ----eval=FALSE---------------------------------------------------------------
# tree_canham_fit <- hmde_model("canham_multi_ind") |>
#   hmde_assign_data(data = Tree_Size_Data)  |>
#   hmde_run(chains = 4, cores = 4, iter = 2000)
# 
# tree_canham_estimates <- hmde_extract_estimates(fit = tree_canham_fit,
#                                  input_measurement_data = Tree_Size_Data)
# 
# #saveRDS(tree_canham_fit, "tree_canham_fit.rds")
# #saveRDS(tree_canham_estimates, "tree_canham_estimates.rds")

## -----------------------------------------------------------------------------
measurement_data_transformed <- Tree_Size_Ests$measurement_data %>%
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
ggplot(measurement_data_transformed, 
       aes(y_hat)) +
  geom_histogram(bins = 10,
                 colour = "black",
                 fill = "lightblue") +
  labs(x="Size (cm)",
       title = "Estimated size distribution") +
  theme_classic()

ggplot(measurement_data_transformed, 
       aes(delta_y_est)) +
  geom_histogram(bins = 10,
                 colour = "black",
                 fill = "lightblue") +
  labs(x = "Growth increment (cm)",
       title = "Estimated growth increments") +
  theme_classic()

ggplot(measurement_data_transformed, 
       aes(est_growth_rate)) +
  geom_histogram(bins = 10,
                 colour = "black",
                 fill = "lightblue") +
  labs(x = "Growth rate (cm/yr)",
       title = "Estimated annualised growth rate distribution") +
  theme_classic()

#Quantitative R^2
cor(measurement_data_transformed$y_obs, measurement_data_transformed$y_hat)^2

r_sq_est <- cor(Tree_Size_Ests$measurement_data$y_obs,
                Tree_Size_Ests$measurement_data$y_hat)^2
r_sq <- paste0("R^2 = ", 
               signif(r_sq_est,
                      digits = 3))

obs_est_scatter <- ggplot(data = Tree_Size_Ests$measurement_data, 
       aes(x = y_obs, y = y_hat)) +
  geom_point(shape = 16, size = 1, colour = "green4") +
  xlab("Y obs.") +
  ylab("Y est.") +
  geom_abline(slope = 1, linetype = "dashed") +
  annotate("text", x = 7, y = 22, 
           label = r_sq) +
  theme_classic()
obs_est_scatter

#Plots of size over time for a sample of 5 individuals
obs_est_size_plot <- hmde_plot_obs_est_inds(n_ind_to_plot = 5,
                       measurement_data = Tree_Size_Ests$measurement_data) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.1, 0.85)) +
  guides(colour=guide_legend(nrow=2,byrow=TRUE))
obs_est_size_plot


## -----------------------------------------------------------------------------
#1-dimensional parameter distributions
ggplot(Tree_Size_Ests$individual_data, 
       aes(ind_max_growth_mean)) +
  geom_histogram(bins = 10,
                 colour = "black",
                 fill = "lightblue") +
  labs(x="g_max estimate") +
  theme_classic()

ggplot(Tree_Size_Ests$individual_data, 
       aes(ind_size_at_max_growth_mean)) +
  geom_histogram(bins = 10,
                 colour = "black",
                 fill = "lightblue") +
  labs(x="S_max estimate") +
  theme_classic()

ggplot(Tree_Size_Ests$individual_data, 
       aes(ind_k_mean)) +
  geom_histogram(bins = 10,
                 colour = "black",
                 fill = "lightblue") +
  labs(x="k estimate") +
  theme_classic()

#2-dimensional parameter distributions
pairplot1 <- ggplot(data = Tree_Size_Ests$individual_data,
       aes(x = ind_max_growth_mean, y = ind_size_at_max_growth_mean)) +
  geom_point(shape = 16, size = 1, colour = "green4") +
  xlab("Ind. max growth (cm/yr)") +
  ylab("Ind. size at max growth (cm)") +
  theme_classic()

pairplot2 <- ggplot(data = Tree_Size_Ests$individual_data,
       aes(x = ind_max_growth_mean, y = ind_k_mean)) +
  geom_point(shape = 16, size = 1, colour = "green4") +
  xlab("Ind. max growth (cm/yr)") +
  ylab("Ind. spread par.") +
  theme_classic()

pairplot3 <- ggplot(data = Tree_Size_Ests$individual_data,
       aes(x = ind_k_mean, y = ind_size_at_max_growth_mean)) +
  geom_point(shape = 16, size = 1, colour = "green4") +
  xlab("Ind. spread par.") +
  ylab("Ind. size at max growth (cm)") +
  theme_classic()

pairplot1
pairplot2
pairplot3

#monotonic correlation of parameters
cor(Tree_Size_Ests$individual_data[,c(2,6,10)], method="spearman")

#Plot function pieces over estimated sizes.
est_de_plot <- hmde_plot_de_pieces(Tree_Size_Ests)
est_de_plot

## -----------------------------------------------------------------------------
pars_CI_names <- c(
  "mean log max growth rate",
  "mean max growth rate in cm/yr",
  "log max growth rate standard deviation",
  "mean log size at max growth rate",
  "mean max size at max growth rate in cm",
  "log max size at growth rate standard deviation",
  "mean mean log spread parameter",
  "mean mean spread parameter",
  "log spread parameter standard deviation"
)

#Vector that picks out which pars to be exponentiated
exp_vec <- c(FALSE, TRUE, FALSE, 
             FALSE, TRUE, FALSE,
             FALSE, TRUE, FALSE)

#Print mean estimates and CIs
for(i in 1:nrow(Tree_Size_Ests$population_data)){
  if(!exp_vec[i]){
    print(paste0(pars_CI_names[i],
                 " estimate: ",
                 Tree_Size_Ests$population_data$mean[i] ))
    print(paste0("95% CI for ", 
                 pars_CI_names[i],
                 ": (",
                 Tree_Size_Ests$population_data$CI_lower[i],
                 ", ",
                 Tree_Size_Ests$population_data$CI_upper[i],
                 ")"))
  } else {
    print(paste0(pars_CI_names[i],
                 " estimate: ",
                 exp(Tree_Size_Ests$population_data$mean[i])))
    print(paste0("95% CI for ",
                 pars_CI_names[i], 
                 ": (",
                 exp(Tree_Size_Ests$population_data$CI_lower[i]),
                 ", ",
                 exp(Tree_Size_Ests$population_data$CI_upper[i]),
                 ")"))
  }
}

