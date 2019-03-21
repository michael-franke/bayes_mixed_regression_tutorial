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

# package to navigate to your source folder
library(rstudioapi)

## Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

#####################################################
## read and massage the data
#####################################################

politedata = read.csv('politeness_data.csv') 
head(politedata)

#####################################################
## plot means for each group
#####################################################

politedata.agg <- 
  politedata %>% 
    group_by(gender, attitude, scenario) %>% 
    summarize(mean_frequency = mean(freq))

politedata.agg2 <- 
  politedata %>%
  group_by(gender, attitude) %>% 
  summarize(mean_frequency = round(mean(freq), 0))

ggplot(data = politedata.agg, 
       aes(x = gender, 
           y = mean_frequency, 
           colour = attitude)) + 
  geom_point(position = position_dodge(0.5), 
             alpha = 0.3, 
             size = 3) +
  geom_point(data = politedata.agg2, 
             aes(x = gender, 
                 y = mean_frequency, 
                 #colour = attitude,
                 fill = attitude),
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
## run model with different random-effects structures
#####################################################

# model with only fixed effects (non-hierarchical)
model_FE = brm(formula = freq ~ gender * attitude, data = politedata)

# model with only a fixed effect for gender (non-hierarchical)
model_gender = brm(formula = freq ~ gender, data = politedata)

# hierarchical model with random intercepts
model_interceptOnly = brm(formula = freq ~ gender * attitude +
                            (1 | scenario + subject),
                          data = politedata)

# hierarchical model with the maximial RE structure licensed by the design
# (notice that factor 'gender' does not vary for a given value of variable 'subject')
model_MaxRE = brm(formula = freq ~ gender * attitude +
                    (1 + gender * attitude | scenario) +
                    (1 + attitude | subject),
                  data = politedata)

#####################################################
## convenience function to extract predictor values &
## compute the probability of the relevant hypotheses
#####################################################

extract_posterior_cell_means = function(model) {
  
  # extract information about dependent and independent variables from formula
  ## TODO :: check extensively, especially for mixed effects models whether this stuff here works
  dependent_variable = as.character(formula(model)[[1]])[[2]]
  independent_variables = strsplit(x = gsub(pattern = "\\(.*\\|.*\\)", "", as.character(formula(model)[[1]])[[3]]),
                                   split =  "(\\*|\\+)",
                                   fixed = FALSE)[[1]] %>% trimws()
  independent_variables = independent_variables[which(independent_variables != "")]
  
  
  # stop this if there are not at least two factors
  if (length(independent_variables) <= 1) {
    stop("Oeps! There do not seem to be at least two factors. If you have no more than one factor, computing cell means is not necessary. Use the estimated coffeficients instead.")
  }
  
  # construct three helpful representations of factors and their levels for the following
  ## factors :: a list with all factors and their levels
  factors = list()
  n_levels = c()
  for (iv in independent_variables) {
    new_levels = list(levels(as.factor(model.frame(model) %>% pull(iv))))
    factors = append(factors, new_levels)
    n_levels = c(n_levels, length(new_levels[[1]]))
  }
  names(factors) = independent_variables
  if (min(n_levels) <= 1){
    stop("Oeps! There seems to be a factor with less than 2 levels. Please check and possibly exclude that factor.")
  }
  ## ref_levels_list :: a list with all factors and their refernce levels
  ref_levels_list = factors
  for (iv in independent_variables) {
    ref_levels_list[[iv]] = ref_levels_list[[iv]][1]
  }
  ## ref_levels :: a string representation of each factor and its reference level
  ref_levels = c()
  for (iv in independent_variables) {
    ref_levels = c(ref_levels, paste0(iv, ref_levels_list[[iv]][1]))
  }
  
  # get the posterior samples for all regression coefficients
  post_samples = posterior_samples(model) %>% select(starts_with("b_"))
  
  # get a table of cells (factor-level combinations) in an ugly format (for internal use)
  cells = expand.grid(factors) 
  for (j in 1:ncol(cells)) {
    levels(cells[,j]) = paste0(colnames(cells)[j], levels(cells[,j]))
  }
  
  # get a table of cells (factor-level combinations) with more readable labels (for final output)
  cells_readable = expand.grid(factors) 
  for (j in 1:ncol(cells_readable)) {
    levels(cells_readable[,j]) = paste0(colnames(cells_readable)[j], ":", levels(cells_readable[,j]))
  }

  # get the names of all estimated coefficients
  coefficient_names = names(post_samples)
  # add the reference levels to the coefficient names (where it is missing/implicit)
  for (c in 1:length(coefficient_names)) {
    for (f in names(factors)) {
      if (!grepl(f, coefficient_names[c])) {
        coefficient_names[c] = paste0(coefficient_names[c], "_", f, ref_levels_list[[f]])
      }
    }
  }
  names(post_samples) = coefficient_names
  
  # two convenience functions to get all coefficients that belong to a design cell
  replace_with_ref_level_recursion = function(cell) {
    which_fcts_are_at_ref_level = map_lgl(cell, function(fl) fl %in% ref_levels)
    if (all(which_fcts_are_at_ref_level)) {
      return(list(cell))
    } 
    else {
      factors_to_replace_with_ref_level = which(which_fcts_are_at_ref_level == F)
      output = list(cell)
      for (f in factors_to_replace_with_ref_level) {
        replaced_cell = cell
        replaced_cell[f] = ref_levels[f]
        output = append(output, replace_with_ref_level(replaced_cell))
      }
      output
    }
  }
  replace_with_ref_level = function(cell) {
    unique(replace_with_ref_level_recursion(cell))
  }
  
  # get samples for the predictor values for each design cell
  predictor_values = map_df(
    1:nrow(cells), 
    function(i) {
      cell = cells[i,1:(ncol(cells))]
      coefficients_to_check = replace_with_ref_level(cell)
      out = 0
      column_indices = map_dbl(1:length(coefficients_to_check), 
              function(j){
                which(
                  map_lgl(coefficient_names, function(coefficient_in_question) {
                    all(map_lgl(coefficients_to_check[[j]], function(c) grepl(c, coefficient_in_question)))
                  }) == T  
                )
              }
            )
      tibble(
        cell = paste(map_chr(1:ncol(cell), function(j) as.character(cells_readable[i,j])), collapse = "__"),
        predictor_value = rowSums(post_samples[column_indices]),
        n_sample = 1:length(predictor_value)
      )
    }
  ) 
  predictor_values %>% spread(key = cell, value = predictor_value)
}

extract_posterior_cell_means(model_gender)
extract_posterior_cell_means(model_FE)
extract_posterior_cell_means(model_interceptOnly)
extract_posterior_cell_means(model_MaxRE)

get_posterior_beliefs_about_hypotheses = function(model) {
  posterior_cell_means = extract_posterior_cell_means(model)
  # insert the comparisons you are interested in as strings 
  tibble(hypothesis = c("Female-polite < Female-informal", 
                        "Male-polite < Male-informal",
                        "Male-informal < Female-polite"),
         probability = c(
           # insert the comparisons you are interested in referring to the extracted samples
           mean(posterior_cell_means$`gender:F__attitude:pol` < posterior_cell_means$`gender:F__attitude:inf`),
           mean(posterior_cell_means$`gender:M__attitude:pol` < posterior_cell_means$`gender:M__attitude:inf`),
           mean(posterior_cell_means$`gender:M__attitude:inf` < posterior_cell_means$`gender:F__attitude:pol`)
         ))
}


extract_comparisons = function(model) {
  # get posterior samples
  post_samples = posterior_samples(model) %>% as.tibble()
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


# new function that relies on generic extractions
get_posterior_beliefs_about_hypotheses(model_FE)
get_posterior_beliefs_about_hypotheses(model_interceptOnly)
get_posterior_beliefs_about_hypotheses(model_MaxRE)

# old function that only works for this particular case
extract_comparisons(model_FE)
extract_comparisons(model_interceptOnly)
extract_comparisons(model_MaxRE)



