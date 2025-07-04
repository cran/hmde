## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# # install.packages("remotes")
# remotes::install_github("traitecoevo/hmde")

## -----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(hmde)

## -----------------------------------------------------------------------------
hmde_model("canham_multi_ind")

## ----eval=FALSE---------------------------------------------------------------
# # Build fit and extract estimates
# canham_multi_ind_fit <- hmde_model("canham_multi_ind") |>
#   hmde_assign_data(data = Tree_Size_Data)  |>
#   hmde_run(chains = 1, cores = 1, iter = 1000)
# 
# Tree_Size_Ests <- hmde_extract_estimates(fit = canham_multi_ind_fit,
#                                  input_measurement_data = Tree_Size_Data)

## -----------------------------------------------------------------------------
summary(Tree_Size_Ests)

# Plot fitted growth function pieces
plot_par_individual_data <- Tree_Size_Ests$individual_data[,c(1, 2, 6, 10)] #Pull out estimates only
hmde_plot_de_pieces(Tree_Size_Ests)

## -----------------------------------------------------------------------------
#Plots of size over time for a sample of 5 individuals
sample_ids <- sample(1:nrow(Tree_Size_Ests$individual_data), size=5) %>%
  sort()
plot_data <- Tree_Size_Ests$measurement_data %>%
  filter(ind_id %in% sample_ids)

ind_size_lims <- Tree_Size_Ests$measurement_data %>%
  filter(ind_id %in% sample_ids)%>%
  group_by(ind_id) %>%
  summarise(y_0 = min(y_hat),
         y_final = max(y_hat))

ggplot(data=plot_data, aes(group = ind_id)) +
  geom_point(aes(x = time, y=y_obs, colour = as.factor(ind_id)),
             shape = 1) +
  geom_line(aes(x = time, y=y_obs, colour = as.factor(ind_id)),
            linetype = "dashed") +
  geom_point(aes(x = time, y=y_hat, colour = as.factor(ind_id)),
             shape = 2) +
  geom_line(aes(x = time, y=y_hat, colour = as.factor(ind_id)),
            linetype = "solid") +
  labs(x="Time (years)", y="DBH (cm)", colour="Ind. ID") +
  theme_classic()

#Load DE for Canham
DE_function = hmde_model_des("canham_multi_ind")

#Produce plot with focus inds
function_plot <- hmde_plot_de_pieces(Tree_Size_Ests,
                    alpha = 0.2) 
for(i in 1:length(sample_ids)){
  args_list <- list(pars=Tree_Size_Ests$individual_data[sample_ids[i],c(2, 6, 10)])
  
  function_plot <- function_plot +
      geom_function(fun=DE_function, args=args_list,
                    colour="black", linewidth=1, alpha = 1,
                    xlim=c(ind_size_lims$y_0[i], ind_size_lims$y_final[i]))
    
}

function_plot

