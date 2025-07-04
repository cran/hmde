## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE
)
# job::job({
# knitr::knit("vignettes/here-be-dragons.Rmd.orig", "vignettes/here-be-dragons.Rmd")
# })

options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
#remotes::install_github("traitecoevo/hmde")

library(hmde)
library(dplyr)
library(ggplot2)
library(deSolve)
library(mixtools)
library(MASS)

## -----------------------------------------------------------------------------
#Change these values to change the model parameters. Must be positive values.
beta_0 <- 10
beta_1 <- 1

#True initial condition
true_y_0 <- 1
max_time <- 9
time <- 0:max_time

#Analytic solution
analytic_solution <- function(x = NULL, pars = NULL){ #Pars is list of beta_0, beta_1, y_0
  return(
    (pars[[1]]/pars[[2]]) + (pars[[3]] - (pars[[1]]/pars[[2]])) * exp(-pars[[2]] * x)
  )
}
true_pars <- list(
  beta_0 = beta_0,
  beta_1 = beta_1,
  true_y_0 = true_y_0
)
true_args_list <- list(pars = c(beta_0,
                             beta_1,
                             true_y_0))

y_true <- analytic_solution(time, true_pars)

## -----------------------------------------------------------------------------
#Produce observations with error and limited precision
y_obs <- round(y_true + rnorm(length(y_true), 0, 0.1), digits = 1)

#Unrounded data if needed
#y_obs <- y_true + rnorm(length(y_true), 0, 0.1)

#Observed data frame
obs_data_frame <- tibble(
  time = time,
  y_obs = y_obs,
  obs_index = 1:length(y_obs)
)

#Have a look at the true and 'observed' data
plot_data <- tibble(
  x = c(time, time),
  y = c(y_true,y_obs),
  Data = c(rep("True sizes", times = length(y_true)),
             rep("Obs. sizes", times = length(y_obs)))
)

sizes_over_time <- ggplot(plot_data, aes(x = x, y = y, group = Data)) +
  geom_point(aes(colour = Data, shape = Data), size = 2) +
  geom_line(aes(colour = Data, linetype = Data), linewidth = 1) +
  scale_linetype_manual(values = c("dashed", "dotted")) +
  ylim(0, 10.5) +
  labs(x = "Time", y = "Y") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.3))
sizes_over_time

#Have a look at the observations against the analytic solution
analytic_observed <- ggplot(obs_data_frame, aes(x = time, y = y_obs)) +
  geom_function(fun=analytic_solution, args = true_args_list,
                linewidth = 1, colour = "black") +
  geom_point(colour = "darkorchid", size = 3) +
  geom_line(colour = "darkorchid", linewidth = 1,
            linetype = "dashed") +
  labs(x = "Time", y = "Y",
       title = "Data simulation") +
  ylim(0, 10.5) +
  theme_classic()
analytic_observed

## ----results='hide'-----------------------------------------------------------
runs <- 100
step_size = 0.5
par_est_tibble <- tibble(run = c(),
                         step_size = c(),
                         beta_0 = c(),
                         beta_1 = c())

for(i in 1:runs){
  #Run the model
  suppressWarnings(
    fit <- hmde_model("affine_single_ind") |>
    hmde_assign_data(n_obs = nrow(obs_data_frame),
                     y_obs = obs_data_frame$y_obs,
                     obs_index = obs_data_frame$obs_index,
                     time = obs_data_frame$time,
                     y_bar = mean(obs_data_frame$y_obs),
                     step_size = step_size,
                     int_method = 1)  |>  #RK4
    hmde_run(chains = 1, cores = 1, iter = 2000)
  )
  
  #Extract parameter estimates
  ests <- hmde_extract_estimates(fit = fit,
                                 input_measurement_data = obs_data_frame)
  
  temp <- tibble(
    run = i,
    step_size = step_size,
    beta_0 = ests$individual_data$ind_beta_0_mean,
    beta_1 = ests$individual_data$ind_beta_1_mean
  )
  
  par_est_tibble <- rbind(par_est_tibble, temp)
} 

## -----------------------------------------------------------------------------
step_size_mix_models <- list()
step_size_mix_model_plots <- list()
step_size_mix_models_par_ests <- tibble(
  good_beta_0 = c(),
  good_beta_1 = c(),
  error_beta_0 = c(),
  error_beta_1 = c(),
  step_size = c(),
  error_fraction = c(),
  dist = c()
)

for(i in 1:length(unique(par_est_tibble$step_size))){
  #Get data for single step size
  step_size_selected <- unique(par_est_tibble$step_size)[i]
  
  analysis_data <- par_est_tibble %>%
    filter(step_size == step_size_selected)
  
  #Get some extreme estimates
  possible_error <- analysis_data %>%
    filter(beta_0 > mean(analysis_data$beta_0))
  
  #To speed up the iterative algorithm we provide some initial conditions
  mu <- list( #Means from true parameters and extreme estimates
    true = c(beta_0, beta_1),
    error = c(mean(possible_error$beta_0), 
              mean(possible_error$beta_1))
  )
  
  #Fit multivariate normal finite mixture model to the estimates
  step_size_mix_models[[i]] <- mvnormalmixEM(x = analysis_data[,c(3,4)], mu = mu)
  
  print(paste0("Summary of mixture model for step size ", step_size_selected))
  print(summary(step_size_mix_models[[i]]))
  
  step_size_mix_model_plots[[i]] <- plot(step_size_mix_models[[i]], 
                               whichplots = 2, 
                               xlab2 = "Beta 0", 
                               ylab2 = "Beta 1")
  
  dist_table <- tibble( #Data to calculate distance
    b_0 = c(step_size_mix_models[[i]]$mu[[2]][1], 
            step_size_mix_models[[i]]$mu[[1]][1]),
    b_1 = c(step_size_mix_models[[i]]$mu[[2]][2], 
            step_size_mix_models[[i]]$mu[[1]][2])
  )
  
  #Extract values
  step_size_mix_models_par_ests_temp <- tibble(
    good_beta_0 = step_size_mix_models[[i]]$mu[[1]][1],
    good_beta_1 = step_size_mix_models[[i]]$mu[[1]][2],
    error_beta_0 = step_size_mix_models[[i]]$mu[[2]][1],
    error_beta_1 = step_size_mix_models[[i]]$mu[[2]][2],
    step_size = step_size_selected,
    error_prob = step_size_mix_models[[i]]$lambda[2],
    dist = dist(dist_table)
  )
  
  step_size_mix_models_par_ests <- rbind(step_size_mix_models_par_ests, 
                                         step_size_mix_models_par_ests_temp)
}

#Have a look at the estimates
step_size_mix_models_par_ests

## -----------------------------------------------------------------------------
legend_spec <- tibble(
  step_size_name = c("0.5", "0.25", "0.125"),
  step_size = c(0.5, 0.25, 0.125),
  x = c(-10, -10, -10),
  y = c(-10, -10, -10),
  colours = c("#f8766d", "#00ba38", "#609cff"),
  linetypes = c("longdash", "dashed", "dotted"),
  shapes = c(19, 17, 15)
) 

legend_spec_with_true <- tibble(
  step_size_name = c("0.5", "0.25", "0.125", "True pars"),
  step_size = c(0.5, 0.25, 0.125, NA),
  x = c(-10, -10, -10, -10),
  y = c(-10, -10, -10, -10),
  colours = c("#f8766d", "#00ba38", "#609cff", "black"),
  linetypes = c("longdash", "dashed", "dotted", "solid"),
  shapes = c(19, 17, 15, 3)
)

for(i in 1:nrow(legend_spec)){
  fancy_name_no_step_size <- 
  paste0("Beta_0 = ",
         signif(step_size_mix_models_par_ests$error_beta_0[i],
                digits = 3),
         ",\n Beta_1 = ",
         signif(step_size_mix_models_par_ests$error_beta_1[i], 
                             digits = 3))
  legend_spec$fancy_name_no_step_size[i] <- fancy_name_no_step_size
  legend_spec_with_true$fancy_name_no_step_size[i] <- fancy_name_no_step_size
  
  fancy_name <- paste0("Step size ", step_size_mix_models_par_ests$step_size[i], 
                       "\n", fancy_name_no_step_size)
  
  legend_spec$fancy_name[i] <- fancy_name
  legend_spec_with_true$fancy_name[i] <- fancy_name
}

legend_spec_with_true$fancy_name_no_step_size[4] <-
  paste0("Beta_0 = ",
         beta_0,
         ",\n Beta_1 = ",
         beta_1)

legend_spec_with_true$fancy_name[4] <-
  paste0("True values\n Beta_0 = ",
         beta_0,
         ",\n Beta_1 = ",
         beta_1)

## -----------------------------------------------------------------------------
scatterplot_errors_only <- list()
scatterplot_good_only <- list()

for(i in 1:length(unique(par_est_tibble$step_size))){
  step_size_select <- unique(par_est_tibble$step_size)[i]
  
  plot_data <- par_est_tibble %>%
    filter(step_size == step_size_select)
  
  #Get classification from mixture model
  plot_data$good_est <- step_size_mix_models[[i]][["posterior"]][,1]
  
  error_ests_scatter <- plot_data %>%
    filter(!as.logical(good_est))
  good_ests_scatter <- plot_data %>%
    filter(as.logical(good_est))
  
  #Scatter plot of erroneous parameters
  xpos <- (min(error_ests_scatter$beta_0) + 
             0.2*(max(error_ests_scatter$beta_0) - 
                    min(error_ests_scatter$beta_0)))
  ypos <- (max(error_ests_scatter$beta_1) - 
             0.1*(max(error_ests_scatter$beta_1) - 
                    min(error_ests_scatter$beta_1)))
  norm_data <- as.data.frame(mvrnorm(n = 10000,
                       mu = step_size_mix_models[[i]][["mu"]][[2]],
                       Sigma = step_size_mix_models[[i]][["sigma"]][[2]]))
  names(norm_data) <- c("beta_0", "beta_1")
    
  scatterplot_errors_only[[i]] <- ggplot(data = error_ests_scatter, 
                                         aes(x = beta_0, y = beta_1)) +
    geom_point(colour = legend_spec$colours[i], 
               shape = legend_spec$shapes[i], 
               alpha = 0.5,
               size = 2) +
    geom_density_2d(data = norm_data, colour = "black") +
    labs(x = "beta_0 est.",
           y = "beta_1 est.",
           title = "Second cluster") +
    annotate("text", x = xpos, y = ypos, 
             label = paste0("Probability: \n",
                            step_size_mix_models_par_ests$error_prob[i])) +
    theme_classic()
  
  #Scatter plot of good parameter estimates
  xpos <- (min(good_ests_scatter$beta_0) + 
             0.2*(max(good_ests_scatter$beta_0) - 
                    min(good_ests_scatter$beta_0)))
  ypos <- (max(good_ests_scatter$beta_1) - 
             0.1*(max(good_ests_scatter$beta_1) - 
                    min(good_ests_scatter$beta_1)))
  norm_data <- as.data.frame(mvrnorm(n = 10000,
                       mu = step_size_mix_models[[i]][["mu"]][[1]],
                       Sigma = step_size_mix_models[[i]][["sigma"]][[1]]))
  names(norm_data) <- c("beta_0", "beta_1")
  
  scatterplot_good_only[[i]] <- ggplot(data = good_ests_scatter, 
                                         aes(x = beta_0, y = beta_1)) +
    geom_point(colour = legend_spec$colours[i], 
               shape = legend_spec$shapes[i], 
               alpha = 0.5,
               size = 2) +
    geom_density_2d(data = norm_data, colour = "black") +
    labs(x = "beta_0 est.",
           y = "beta_1 est.",
           title = "First cluster") +
    annotate("text", x = xpos, y = ypos, 
             label = paste0("Probability: \n",
                            (1-step_size_mix_models_par_ests$error_prob[i]))) +
      theme_classic()
}

## -----------------------------------------------------------------------------
#install.packages("deSolve")
library(deSolve)

#Create DE function
DE <- function(Time, State, Pars) { #Implementation of DE
  with(as.list(c(State, Pars)), {
    dY <- beta_0 - beta_1 * Y

    return(list(c(dY)))
  })
}

## -----------------------------------------------------------------------------
yini  <- c(Y = true_y_0) #Initial condition
y_over_time <- tibble(model="True Sizes",
                      y_analytic = y_true,
                      y_numeric = y_true,
                      time = 0:max_time,
                      beta_0_par = beta_0,
                      beta_1_par = beta_1
                      )

#Generate Y(t) with RK4
for(i in 1:nrow(step_size_mix_models_par_ests)){
  pars_combo <- c(beta_0 = step_size_mix_models_par_ests$error_beta_0[i],
                    beta_1 = step_size_mix_models_par_ests$error_beta_1[i])
  times <- seq(0, max_time, by = step_size_mix_models_par_ests$step_size[i])
  
  solution_pars <- c(pars_combo, true_y_0)
  y_true_temp <- analytic_solution(times, solution_pars)
    
  numerical_output <- ode(yini, times, DE, pars_combo, method = "rk4")
  
  y_over_time_temp <- tibble(
    model = step_size_mix_models_par_ests$step_size[i],
    y_analytic = y_true_temp,
    y_numeric = numerical_output[,2],
    time = times,
    beta_0_par = step_size_mix_models_par_ests$error_beta_0[i],
    beta_1_par = step_size_mix_models_par_ests$error_beta_1[i]
  )
  
  y_over_time <- rbind(y_over_time, y_over_time_temp)
}

## -----------------------------------------------------------------------------
y_over_time_filtered <- y_over_time %>%
  filter(time %in% 0:max_time)
  
#Plot sizes over time for all models
compare_sizes_over_time <- ggplot(y_over_time_filtered, 
                                  aes(x=time, y=y_numeric, group_by = as.factor(model))) +
  geom_point(aes(colour = as.factor(model),
             shape = as.factor(model)),
             alpha=0.5, size = 2, stroke = 1.5) +
  geom_line(aes(colour = as.factor(model)), alpha=0.5, linewidth = 1) +
  scale_colour_manual(values = legend_spec_with_true$colours) +
  scale_shape_manual(values = legend_spec_with_true$shapes) +
  labs(x = "Time", y = "Y(t)", title = "Estimated Y(t) with bad parameters",
       colour = "Step size", shape = "Step size") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.3))

compare_sizes_over_time

## -----------------------------------------------------------------------------
#Generate y(t) with RK4 given the parameter estimates
y_over_time_smallstep <- tibble(model=legend_spec_with_true$fancy_name[4],
                      y_hat = y_true,
                      time = 0:max_time
                      )

for(i in 1:nrow(step_size_mix_models_par_ests)){
  pars_combo <- c(beta_0 = step_size_mix_models_par_ests$error_beta_0[i],
                    beta_1 = step_size_mix_models_par_ests$error_beta_1[i])
    times <- seq(0, max_time, by = 0.001)
    
    numerical_output <- ode(yini, times, DE, pars_combo, method = "rk4")
    
    y_over_time_temp <- tibble(
      model = legend_spec_with_true$fancy_name_no_step_size[i],
      y_hat = numerical_output[,2],
      time = times
    )
    
    y_over_time_smallstep <- rbind(y_over_time_smallstep, y_over_time_temp)
}

point_data <- y_over_time_smallstep %>%
  filter(time %in% 0:max_time)

#Plot sizes over time
compare_sizes_over_time_smallstep <- ggplot(y_over_time_smallstep, 
                                  aes(x=time, y=y_hat, grouping = as.factor(model))) +
  geom_line(aes(colour = as.factor(model),
                linetype = as.factor(model)), alpha=0.5, linewidth = 1) +
  geom_point(data = point_data,
             aes(colour = as.factor(model),
             shape = as.factor(model)),
             alpha=0.5, size = 2, stroke = 1.5) +
  geom_function(fun=analytic_solution, args=true_args_list,
                colour="black",
                linetype = "solid",
                linewidth=1) +
  scale_shape_manual(values = legend_spec_with_true$shapes) +
  scale_colour_manual(values = legend_spec_with_true$colours) +
  scale_linetype_manual(values = c(legend_spec$linetypes, NA)) +
  labs(x = "Time", y = "Y(t)", title = "Small step size",
       colour = "Parameters", 
       shape = "Parameters",
       linetype = "Parameters") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.4),
        legend.key.spacing.y = unit(2, 'mm')) +
  guides(colour = guide_legend(byrow = TRUE))

compare_sizes_over_time_smallstep

## -----------------------------------------------------------------------------
#Get asymptotic size
step_size_mix_models_par_ests$Y_max <- step_size_mix_models_par_ests$error_beta_0/step_size_mix_models_par_ests$error_beta_1

#Build points for start and end of lines
plot_data <- tibble(
  x = c(0, (beta_0/beta_1), 
        rep(0, times = nrow(step_size_mix_models_par_ests)), 
        step_size_mix_models_par_ests$Y_max),
  y = c(beta_0, 0, 
        step_size_mix_models_par_ests$error_beta_0, 
        rep(0, times = nrow(step_size_mix_models_par_ests))),
  step_size = c("True pars", "True pars", 
                step_size_mix_models_par_ests$step_size, 
                step_size_mix_models_par_ests$step_size)
)

#Plot DEs 
error_de_plot <- ggplot(data = plot_data, aes(x,y)) +
  geom_line(aes(colour = as.factor(step_size),
                linetype = as.factor(step_size)),
            linewidth = 1) +
  scale_colour_manual(values = c(legend_spec$colours[3:1], "black")) +
  scale_linetype_manual(values = c(legend_spec$linetypes[3:1], "solid")) +
  labs(title = "ODEs",
       x = "Y(t)", y = "f", 
       colour = "Step size", 
       linetype = "Step size") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.7))
error_de_plot

#Plot analytic solutions
error_solution_plot <- ggplot() +
  geom_function(fun=analytic_solution, args=true_args_list,
                  colour="black", 
                  linetype = "solid",
                  linewidth=1) 

for(i in 1:nrow(step_size_mix_models_par_ests)){ #Add the analytic solutions
  args_list <- list(pars=c(step_size_mix_models_par_ests$error_beta_0[i],
                           step_size_mix_models_par_ests$error_beta_1[i],
                           true_y_0))
  error_solution_plot <- error_solution_plot +
    geom_function(fun=analytic_solution, args=args_list,
                  colour=legend_spec$colours[i], 
                  linetype = legend_spec$linetypes[i],
                  linewidth=1)
}

error_solution_plot <- error_solution_plot +
  geom_line(data = legend_spec_with_true,
            linewidth=1,
            aes(colour = fancy_name_no_step_size,
                linetype = fancy_name_no_step_size,
                x = x, y = y)) +
  scale_colour_manual(values = c(legend_spec_with_true$colours[4], 
                                 legend_spec$colours[c(2,3,1)])) +
  scale_linetype_manual(values = c(legend_spec_with_true$linetypes[4], 
                                   legend_spec$linetypes[c(2,3,1)])) +
  xlim(0, max_time) +
  ylim(true_y_0, (beta_0/beta_1+0.5)) +
  labs(x = "Time", y = "Y(t)", 
       title = "Analytic solutions",
       colour = "Parameters", 
       linetype = "Parameters") +
  theme_classic() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.4),
        legend.key.spacing.y = unit(2, 'mm')) +
  guides(colour = guide_legend(byrow = TRUE))

error_solution_plot

