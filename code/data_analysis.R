#####################################################
## package includes and options
#####################################################

# package for convenience functions (e.g. plotting)
library(tidyverse)

# package for Bayesian regression modeling
library(brms)
# option for Bayesian regression models: 
# use all available cores for parallel computing
options(mc.cores = parallel::detectCores())

# package for credible interval computation
library(HDInterval)

# set seed
set.seed(1702)

#####################################################
## read and massage the data
#####################################################

politedata = read_csv('https://raw.githubusercontent.com/michael-franke/bayes_mixed_regression_tutorial/master/code/politeness_data.csv') 
head(politedata)

#####################################################
## plot means for each group
#####################################################

politedata.agg <- 
  politedata %>% 
    group_by(gender, context, sentence) %>% 
    summarize(mean_frequency = mean(pitch))

politedata.agg2 <- 
  politedata %>%
  group_by(gender, context) %>% 
  summarize(mean_frequency = round(mean(pitch), 0))

ggplot(data = politedata.agg, 
       aes(x = gender, 
           y = mean_frequency, 
           colour = context)) + 
  geom_point(position = position_dodge(0.5), 
             alpha = 0.3, 
             size = 3) +
  geom_point(data = politedata.agg2, 
             aes(x = gender, 
                 y = mean_frequency, 
                 #colour = context,
                 fill = context),
             position = position_dodge(0.5), 
             pch = 21, 
             colour = "black",
             size = 5) +
  scale_x_discrete(breaks = c("F", "M"),
                  labels = c("female", "male")) +
  scale_y_continuous(expand = c(0, 0), breaks = (c(50,100,150,200,250,300)), limits = c(50,300)) +
  scale_colour_manual(breaks = c("inf", "pol"),
                      labels = c("informal", "polite"),
                      values = c("#f1a340", "#998ec3")) +
  scale_fill_manual(breaks = c("inf", "pol"),
                      labels = c("informal", "polite"),
                      values = c("#f1a340", "#998ec3")) +
  ylab("pitch in Hz\n") +
  xlab("\ngender") +
  theme_classic() +
  theme(legend.position = "right",
        legend.key.height = unit(2,"line"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, face = "bold"),
        axis.line.x = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"),
        plot.margin = unit(c(0.2,0.1,0.2,0.1),"cm"))
  
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

ggsave(filename = "../text/pics/basic_data_plot.pdf",
       plot = last_plot(),
       width = 6, height = 4)

#####################################################
## run & inspect model with only fixed effects
#####################################################

formulaFE = pitch ~ gender * context

# model with only fixed effects (non-hierarchical)
modelFE = brm(formula = formulaFE, data = politedata)

# extract posterior samples 
post_samples_FE = posterior_samples(modelFE)
head(post_samples_FE)

# plotting the posterior distributions
plot_posterior_density_FE = 
  modelFE %>% as.tibble() %>% 
  select(- lp__, - sigma) %>% 
  gather(key = "parameter", value = "posterior") %>% 
  mutate(parameter = case_when(parameter == "b_Intercept" ~ "Intercept",
                               parameter == "b_contextpol" ~ "context:pol",
                               parameter == "b_genderM" ~ "gender:M",
                               parameter == "b_genderM.contextpol" ~ "gender:M__context:pol")) %>% 
  mutate(parameter = as.factor(parameter)) %>% 
  mutate(parameter = factor(parameter, levels = c("Intercept", "context:pol", "gender:M", "gender:M__context:pol"))) %>% 
  ggplot(aes(x = posterior)) + 
    geom_density(fill = "grey") +
    facet_wrap(~ parameter, scales = "free") +
  ylab("density\n") +
  xlab("\nparameter value") +
  theme_classic() +
  theme(legend.position = "right",
        legend.key.height = unit(2,"line"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, face = "bold"),
        axis.line = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"),
        plot.margin = unit(c(0.2,0.1,0.2,0.1),"cm")) +
  geom_segment(
    mapping = aes(y = 0, yend = 0, x = low, xend = high),
    color = "firebrick",
    size = 3,
    #alpha = 0.7,
    data = fixef(modelFE) %>% as.tibble() %>%
      mutate(parameter = c("Intercept", "gender:M", "context:pol", "gender:M__context:pol")) %>%
      mutate(parameter = as.factor(parameter)) %>% 
      mutate(parameter = factor(parameter, levels = c("Intercept", "context:pol", "gender:M", "gender:M__context:pol"))) %>% 
      rename(low = Q2.5, high = Q97.5)
      )

# save the plotted figure
ggsave(plot = last_plot(), filename = "../text/pics/posterior_density_FE.pdf",
       width = 9, height = 6)

# proportion of negative samples for parameter p_contextpol
# this number approximates P(beta_pol < 0 | model, data)
mean(post_samples_FE$b_contextpol < 0)

# proportion of samples where the mean for cell 2 was bigger 
# than that of cell 3 
# this number approximates P(beta_pol > beta_male | model, data)
mean(post_samples_FE$b_contextpol > post_samples_FE$b_genderM)

###########################################
## showcasing the faintr package
## (still hoping for the best)
###########################################

# package to allow installation from github
library(devtools)

# package with convenience function for Bayesian regression models for factorial designs
install_github(
  repo = 'michael-franke/bayes_mixed_regression_tutorial', 
  subdir = 'faintr',
  build_vignettes = TRUE) # install from GitHub
library(faintr)

# extract cell means and plot them
posterior_cell_means = extract_posterior_cell_means(modelFE)$predictor_values %>% 
  select(- n_sample) %>% 
  gather(key = "parameter", value = "posterior") 

posterior_cell_means_HDIs = posterior_cell_means %>% 
  group_by(parameter) %>% 
  summarize(low = hdi(posterior)[1],
            high = hdi(posterior)[2])

posterior_cell_means_plot = posterior_cell_means %>% 
  # mutate(parameter = as.factor(parameter)) %>% 
  # mutate(parameter = factor(parameter, levels = c("Intercept", "context:pol", "gender:M", "gender:M__context:pol"))) %>% 
  ggplot(aes(x = posterior)) + 
  geom_density(fill = "grey") +
  facet_wrap(~ parameter, scales = "free") +
  ylab("density\n") +
  xlab("\nparameter value") +
  theme_classic() +
  theme(legend.position = "right",
        legend.key.height = unit(2,"line"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16),
        legend.background = element_rect(fill = "transparent"),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, face = "bold"),
        axis.line = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 18, face = "bold"),
        plot.margin = unit(c(0.2,0.1,0.2,0.1),"cm")) +
  geom_segment(
    mapping = aes(y = 0, yend = 0, x = low, xend = high),
    color = "firebrick",
    size = 3,
    #alpha = 0.7,
    data = posterior_cell_means_HDIs
  )
  
# save the plotted figure
ggsave(plot = last_plot(), filename = "../text/pics/posterior_density_cell_means.pdf",
       width = 9, height = 6)


compare_groups(
  model = modelFE, 
  lower = list(gender = "M", context = "inf"),
  higher = list(gender = "F", context = "pol")
)

get_posterior_beliefs_about_hypotheses_new = function(model) {
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

get_posterior_beliefs_about_hypotheses_new(modelFE)


############################
## add prior information
############################

# get all possible priors for your model
get_prior(formula = pitch ~ gender * context,
          data = politedata)

# define priors
priorFE <- c(
  # define a regularizing prior for the intercept within the range of possible pitch values
  prior(normal(194, 66), class = Intercept),
  # define a skeptical prior for the relevant coefficients
  prior(normal(0, 50), class = b)
)


# let's run our models with our specified priors
modelFE_prior = brm(formula = pitch ~ gender * context,
                        prior = priorFE,
                        data = politedata,
                        control = list(adapt_delta = 0.99))

get_posterior_beliefs_about_hypotheses_new(modelFE_prior)

####################
## model check
####################

# run model without considering gender
modelFE_noGender = brm(formula =  pitch ~ context,
                        data = politedata,
                        prior = priorFE,
                        control = list(adapt_delta = 0.99))

pp_check(modelFE_noGender, nsample = 100)

# save the plotted figure
ggsave(plot = last_plot(), filename = "../text/pics/pp_check_FE_noGender.pdf",
       width = 6, height = 6)

# model with gender
pp_check(modelFE, nsample = 100)
ggsave(plot = last_plot(), filename = "../text/pics/pp_check_FE.pdf",
       width = 6, height = 6)


###############################################
## models with additional random effects
###############################################

# hierarchical model with random intercepts

# define priors
prior_interceptOnly <- c(
  # define a regularizing prior for the intercept within the range of possible pitch values
  prior(normal(194, 66), class = Intercept),
  # define a skeptical prior for the relevant coefficients
  prior(normal(0, 50), class = b)
)

# model
model_interceptOnly = brm(formula = pitch ~ gender * context +
                            (1 | sentence + subject),
                          data = politedata,
                          prior = prior_interceptOnly,
                          control = list(adapt_delta = 0.99))

# hierarchical model with the maximial RE structure licensed by the design
# (notice that factor 'gender' does not vary for a given value of variable 'subject')


# define priors
prior_MaxRE <- c(
  # define a regularizing prior for the intercept within the range of possible pitch values
  prior(normal(194, 66), class = Intercept),
  # define a skeptical prior for the relevant coefficients
  prior(normal(0, 50), class = b)
)

# model
model_MaxRE = brm(formula = pitch ~ gender * context +
                    (1 + gender * context | sentence) +
                    (1 + context | subject),
                  prior = prior_MaxRE,
                  data = politedata,
                  control = list(adapt_delta = 0.99))

##################################
## comparing selected hypotheses
##################################

get_posterior_beliefs_about_hypotheses_new(modelFE)
get_posterior_beliefs_about_hypotheses_new(model_interceptOnly)
get_posterior_beliefs_about_hypotheses_new(model_MaxRE)

