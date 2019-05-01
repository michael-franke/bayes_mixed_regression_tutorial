#####################################################
## package includes and options
#####################################################

# package for convenience functions (e.g. plotting)
library(tidyverse)

# package for Bayesian regression modeling
library(brms)
# option for Bayesian regression models: use all available cores for parallel computing
options(mc.cores = parallel::detectCores())

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
                      values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(breaks = c("inf", "pol"),
                      labels = c("informal", "polite"),
                      values = c("#0072B2", "#D55E00")) +
  ylab("pitch in Hz\n") +
  xlab("\nGender") +
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
plot_posterior_density_FE = modelFE %>% as.tibble() %>% 
  select(- lp__, - sigma) %>% 
  gather(key = "parameter", value = "posterior") %>% 
  ggplot(aes(x = posterior)) + geom_density() +
  facet_wrap(~ parameter, scales = "free")

# save the plotted figure
ggsave(plot = plot_posterior_density_FE, filename = "../text/pics/posterior_density_FE.pdf",
       width = 6, height = 4)

# proportion of negative samples for parameter p_contextpol
# this number approximates P(beta_pol < 0 | model, data)
mean(post_samples_FE$b_contextpol < 0)

# proportion of samples where the mean for cell 2 was bigger 
# than that of cell 3 
# this number approximates P(beta_pol > beta_male | model, data)
mean(post_samples_FE$b_contextpol - post_samples_FE$b_genderM > 0)

###########################################
## showcasing the faintr package
## (still hoping for the best)
###########################################

# package to allow installation from github
library(devtools)

# package with convenience function for Bayesian regression models for factorial designs
install_github('michael-franke/bayes_mixed_regression_tutorial/faintr') # install from GitHub
library(faintr)

extract_posterior_cell_means(modelFE)

get_posterior_beliefs_about_hypotheses_new = function(model) {
  # insert the comparisons you are interested in as strings 
  tibble(
    hypothesis = c("Female-polite < Female-informal", 
                   "Male-polite < Male-informal",
                   "Male-informal < Female-polite"),
    probability = c(
      # insert the comparisons you are interested in referring to the extracted samples
      get_cell_comparison(
        model = model, 
        cell_low = list(gender = "F", context = "pol"), 
        cell_high = list(gender = "F", context = "inf")
      ),
      get_cell_comparison(
        model = model, 
        cell_low = list(gender = "M", context = "pol"), 
        cell_high = list(gender = "M", context = "inf")
      ),
      get_cell_comparison(
        model = model, 
        cell_low = list(gender = "M", context = "inf"),
        cell_high = list(gender = "F", context = "pol")
      ) 
    )
  )
}

get_posterior_beliefs_about_hypotheses_new(model_FE)

###############################################
## models with additional random effects
###############################################

# hierarchical model with random intercepts
model_interceptOnly = brm(formula = pitch ~ gender * context +
                            (1 | sentence + subject),
                          data = politedata,
                          control = list(adapt_delta = 0.99))

# hierarchical model with the maximial RE structure licensed by the design
# (notice that factor 'gender' does not vary for a given value of variable 'subject')
model_MaxRE = brm(formula = pitch ~ gender * context +
                    (1 + gender * context | sentence) +
                    (1 + context | subject),
                  data = politedata,
                  control = list(adapt_delta = 0.9))

##################################
## comparing selected hypotheses
##################################

get_posterior_beliefs_about_hypotheses_new(model_FE)
get_posterior_beliefs_about_hypotheses_new(model_interceptOnly)
get_posterior_beliefs_about_hypotheses_new(model_MaxRE)

