## Code to follow the tutorial: "Bayesian regression modeling (for factorial designs): A tutorial"
## Authors: Michael Franke & Timo Roettger
## Date last modified: 21/06/19
## Contact: mchfranke@gmail.com; timo.b.roettger@gmail.com


####################
## install packages
#####################

# package for convenience functions (e.g. plotting)
library(tidyverse)

# package for Bayesian regression modeling
library(brms)

# option for Bayesian regression models:
# use all available cores for parallel computing
options(mc.cores = parallel::detectCores())

# package for credible interval computation
library(HDInterval)

# set the random seed in order to make sure
# you can reproduce the same results
set.seed(1702)

##################
## load the data
##################

# load the data into variable "politedata"
politedata = read_csv("https://raw.githubusercontent.com/michael-franke/bayes_mixed_regression_tutorial/master/code/politeness_data.csv")

# inspect head of data
head(politedata)

#####################################################
## run & inspect model with only fixed effects
#####################################################

# formula for fixed effect regression model
formula_FE = pitch ~ gender * context

# run regression model in brms
model_FE = brm(
  formula = formula_FE,
  data = politedata,
  seed = 1702
)

# print out model_FE summary
model_FE

# extract posterior samples
post_samples_FE = posterior_samples(model_FE)
head(post_samples_FE %>% round(1))

# proportion of negative samples for parameter b_contextpol
# this number approximates P(b_contextpol < 0 | model, data)
mean(post_samples_FE$b_contextpol < 0)

# proportion of samples for which the mean of cell 2 was larger
# than that of cell 3
# this number approximates the quantity:
# P(b_contextpol > b_genderM | model, data)
mean(post_samples_FE$b_contextpol > post_samples_FE$b_genderM)


###########################################
## showcasing the faintr package
## (still hoping for the best)
###########################################

# load package to allow installation from GitHub
library(devtools)

# install package with convenience function for Bayesian regression
# models for factorial designs from GitHub
install_github(
  repo = "michael-franke/bayes_mixed_regression_tutorial",
  subdir = "faintr")

# load the just installed package
library(faintr)

# extract posterior cell means
extract_posterior_cell_means(model_FE)$predictor_values

# extract cell means and plot them
posterior_cell_means = extract_posterior_cell_means(model_FE)$predictor_values %>% 
  gather(key = "parameter", value = "posterior") 

# compare cell means with each other
compare_groups(
  model = model_FE,
  lower = list(gender = "M", context = "inf"),
  higher = list(gender = "F", context = "pol")
)


get_posterior_beliefs_about_hypotheses = function(model) {
  # insert the comparisons you are interested in as strings 
  tibble(
    hypothesis = c("Female-polite < Female-informal", 
                   "Male-polite < Male-informal",
                   "Male-informal < Female-polite"),
    probability = c(
      # insert the comparisons you are interested in referring to the extracted samples
      compare_groups(
        model = model, 
        lower = list(gender = "F", context = "pol"), 
        higher = list(gender = "F", context = "inf")
      )$probability,
      compare_groups(
        model = model, 
        lower = list(gender = "M", context = "pol"), 
        higher = list(gender = "M", context = "inf")
      )$probability,
      compare_groups(
        model = model, 
        lower = list(gender = "M", context = "inf"),
        higher = list(gender = "F", context = "pol")
      )$probability 
    )
  )
}

get_posterior_beliefs_about_hypotheses(model_FE)

#########################
## add prior information
#########################

# see the priors of your fitted model
prior_summary(model_FE)

# get all possible priors for your model before fitting it
# (here prior_summary and get_prior give identical outputs 
# because you have not specified any priors beyond the default priors
get_prior(formula = pitch ~ gender * context,
          data = politedata)

# define priors
priorFE <- c(
  # define a skeptical prior for the effect of context on female speakers
  prior(normal(0, 10), coef = contextpol)
)

# let's run our models with our specified priors
model_FE_prior = brm(formula = pitch ~ gender * context,
                     prior = priorFE, # add prior
                     data = politedata,
                     control = list(adapt_delta = 0.99),
                     seed = 1702)

# Extract posterior beliefs about our hypotheses
get_posterior_beliefs_about_hypotheses(model_FE_prior)

###############
## model check
###############

# run model without considering gender
model_FE_noGender = brm(formula =  pitch ~ context,
                        data = politedata,
                        control = list(adapt_delta = 0.99),
                        seed = 1702)

# perform posterior predictive check fo both models
pp_check(model_FE_noGender, nsample = 100)
pp_check(model_FE, nsample = 100) 

###############################################
## models with additional random effects
###############################################

# hierarchical model with the maximial RE structure licensed by the design
# (notice that factor 'gender' does not vary for a given value of variable 'subject')
model_MaxRE = brm(formula = pitch ~ gender * context +
                    (1 + gender * context | sentence) +
                    (1 + context | subject),
                  data = politedata,
                  control = list(adapt_delta = 0.99),
                  seed = 1702)

# Extract posterior beliefs about our hypotheses
get_posterior_beliefs_about_hypotheses(model_MaxRE)

##################################
## comparing selected hypotheses
##################################

get_posterior_beliefs_about_hypotheses(model_FE)
get_posterior_beliefs_about_hypotheses(model_interceptOnly)
get_posterior_beliefs_about_hypotheses(model_MaxRE)

#################################
## posteriors of cell differences
## for final data report
#################################

compare_groups(
  model = model_MaxRE, 
  lower = list(gender = "F", context = "pol"), 
  higher = list(gender = "F", context = "inf")
)

compare_groups(
  model = model_MaxRE, 
  lower = list(gender = "M", context = "pol"), 
  higher = list(gender = "M", context = "inf")
)

compare_groups(
  model = model_MaxRE, 
  lower = list(gender = "M", context = "inf"),
  higher = list(gender = "F", context = "pol")
)