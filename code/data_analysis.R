#####################################################
## package includes and options
#####################################################

# library for convenience functions (e.g. plotting)
library(tidyverse)

# library for Bayesian regression modeling
library(brms)
# option for Bayesian regression models: use all available cores for parallel computing
options(mc.cores = parallel::detectCores())

# package for function 'std.error' to obtain standard errors
library(plotrix)

#####################################################
## read and massage the data
#####################################################

politedata = read.csv('politeness_data.csv') 
head(politedata)

#####################################################
## plot means for each group
#####################################################

politedata %>% 
  group_by(gender, attitude) %>% 
  summarize(mean_frequency = mean(freq),
            standard_error = std.error(freq)) %>% 
  ggplot((aes(x = gender, y = mean_frequency, fill = attitude))) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = mean_frequency - standard_error,
                    ymax = mean_frequency + standard_error), position = "dodge")

ggsave(filename = "../text/pics/basic_data_plot.pdf",
       plot = last_plot(),
       width = 6, height = 4)

#####################################################
## run model with different random-effects structures
#####################################################

# model with only fixed effects (non-hierarchical)
model_FE = brm(formula = frequency ~ gender * attitude, data = politeness_data)

# hierarchical model with random intercepts
model_interceptOnly = brm(formula = frequency ~ gender * attitude +
                            (1 | scenario + subject), 
                          data = politeness_data)

# hierarchical model with the maximial RE structure licensed by the design
# (notice that factor 'gender' does not vary for a given value of variable 'subject')
model_MaxRE = brm(formula = frequency ~ gender * attitude +
                    (1 + gender * attitude | scenario) +
                    (1 + attitude | subject), 
                  data = politeness_data)

#####################################################
## convenience function to extract predictor values &
## compute the probability of the relevant hypotheses
#####################################################

extract_comparisons = function(model) {
  # get posterior samples
  post_samples = posterior_samples(model)
  # mnemonic names for reconstructed predictor values for all cells in the design matrix
  F_inf = post_samples$b_Intercept
  F_pol = post_samples$b_Intercept + 
    post_samples$b_attitudepol
  M_inf = post_samples$b_Intercept + 
    post_samples$b_genderM
  M_pol = post_samples$b_Intercept + 
    post_samples$b_genderM + 
    post_samples$b_attitudepol + 
    post_samples$`b_genderM:attitudepol` 
  tibble(hypothesis = c("Female-polite < Female-informal", 
                        "Male-polite < Male-informal",
                        "Male-informal < Female-polite"),
         probability = c(
           mean(F_pol < F_inf),
           mean(M_pol < M_inf),
           mean(M_inf < F_pol)
         ))
}

#####################################################
## access probability of the relevant hypotheses
## under different models
#####################################################

extract_comparisons(model_FE)
extract_comparisons(model_interceptOnly)
extract_comparisons(model_MaxRE)



