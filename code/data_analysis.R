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

formula_FE = pitch ~ gender * context

# model with only fixed effects (non-hierarchical)
model_FE = brm(
  formula = formula_FE, 
  data = politedata, 
  seed = 1702)

# extract posterior samples 
post_samples_FE = posterior_samples(model_FE)
head(post_samples_FE %>% round(1))

# plotting the posterior distributions
plot_posterior_density_FE = 
  model_FE %>% as_tibble() %>% 
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
    data = fixef(model_FE) %>% as_tibble() %>%
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
# install_github(
#   repo = 'michael-franke/bayes_mixed_regression_tutorial',
#   subdir = 'faintr'
# ) # install from GitHub

library(faintr)

# extract cell means and plot them
posterior_cell_means = extract_posterior_cell_means(model_FE)$predictor_values %>% 
  gather(key = "parameter", value = "posterior") 

posterior_cell_means_HDIs = posterior_cell_means %>% 
  group_by(parameter) %>% 
  summarize(low = hdi(posterior)[1],
            high = hdi(posterior)[2])

posterior_cell_means_plot = posterior_cell_means %>% 
  mutate(parameter = as.factor(parameter)) %>% 
  mutate(parameter = factor(parameter, labels = c("female - informal", "female - polite", "male - informal", "male - polite"))) %>% 
  ggplot(aes(x = posterior)) + 
  geom_density(fill = "grey") +
  facet_wrap(~ parameter, scales = "free") +
  ylab("density\n") +
  xlab("\nparameter value") +
  scale_x_continuous(expand = c(0, 0), breaks = (c(100,200,300)), limits = c(80,300)) + 
  scale_y_continuous(expand = c(0, 0), breaks = (c(0,0.02,0.04,0.06)), limits = c(0,0.06)) + 
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
        plot.margin = unit(c(0.2,0.4,0.2,0.1),"cm")) +
  geom_segment(
    mapping = aes(y = 0, yend = 0, x = low, xend = high),
    color = "firebrick",
    size = 3,
    #alpha = 0.7,
    data = posterior_cell_means_HDIs %>% 
      mutate(parameter = as.factor(parameter)) %>% 
      mutate(parameter = factor(parameter, labels = c("female - informal", "female - polite", "male - informal", "male - polite"))) 
  )
  
# save the plotted figure
ggsave(plot = last_plot(), filename = "../text/pics/posterior_density_cell_means.pdf",
       width = 9, height = 6)

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

############################
## add prior information
############################

# get all possible priors for your model
get_prior(formula = pitch ~ gender * context,
          data = politedata)

# define priors
priorFE <- c(
  # define a skeptical prior for the relevant coefficients
  prior(normal(0, 10), coef = contextpol)
)

# let's run our models with our specified priors
model_FE_prior = brm(formula = pitch ~ gender * context,
                    prior = priorFE,
                    data = politedata,
                    control = list(adapt_delta = 0.99),
                    seed = 1702)                        

get_posterior_beliefs_about_hypotheses(model_FE_prior)

####################
## model check
####################

# run model without considering gender
model_FE_noGender = brm(formula =  pitch ~ context,
                       data = politedata,
                       control = list(adapt_delta = 0.99),
                       seed = 1702)

pp_check1 <- pp_check(model_FE_noGender, nsample = 100) +
  scale_color_manual(values = c("#f1a340", "lightgrey")) +
  scale_x_continuous(breaks = (c(0,100,200,300,400)), limits = c(-50,400)) +
  scale_y_continuous(limits = c(0,0.01)) +
  labs(title = "Model without gender") +
  theme(legend.position = "none",
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
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.4),"cm"))

# model with gender
pp_check2 <- pp_check(model_FE, nsample = 100) +
  scale_color_manual(values = c("#f1a340", "lightgrey")) +
  scale_x_continuous(breaks = (c(0,100,200,300,400)), limits = c(-50,400)) +
  scale_y_continuous(limits = c(0,0.01)) +
  labs(title = "Model including gender") +
  theme(legend.position = "none",
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
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.margin = unit(c(0.2,0.4,0.2,0.1),"cm"))

# combine plots
library(ggpubr)

pp_checks_plot <- 
  ggarrange(pp_check1, pp_check2,
            heights = c(1,1),
            widths = c(1,1),
            #labels = c("A - Model without gender", "B - Model including gender"), 
            font.label = list(size = 20), legend = "none",
            align = "h",
            ncol = 2, nrow = 1)

ggsave(plot = pp_checks_plot, filename = "../text/pics/pp_checks_plot.pdf",
       width = 8, height = 4)

###############################################
## models with additional random effects
###############################################

# hierarchical model with random intercepts

# model
model_interceptOnly = brm(formula = pitch ~ gender * context +
                            (1 | sentence + subject),
                          data = politedata,
                          control = list(adapt_delta = 0.99),
                          seed = 1702)

# hierarchical model with the maximial RE structure licensed by the design
# (notice that factor 'gender' does not vary for a given value of variable 'subject')

# model
model_MaxRE = brm(formula = pitch ~ gender * context +
                    (1 + gender * context | sentence) +
                    (1 + context | subject),
                  data = politedata,
                  control = list(adapt_delta = 0.99),
                  seed = 1702)

# extract cell means and 95% CIs
posterior_cell_means = extract_posterior_cell_means(model_MaxRE)$predictor_values %>% 
  gather(key = "parameter", value = "posterior") %>% 
  group_by(parameter) %>% 
  summarize(mean = mean(posterior),
              low = hdi(posterior)[1],
              high = hdi(posterior)[2])

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
