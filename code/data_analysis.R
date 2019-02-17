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
    ggplot((aes(x = gender, 
                y = mean_frequency, 
                fill = attitude))) + 
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_frequency - standard_error,
                      ymax = mean_frequency + standard_error), 
                  position = "dodge")

ggsave(filename = "../text/pics/basic_data_plot.pdf",
       plot = last_plot(),
       width = 6, height = 4)

#####################################################
## run model with different random-effects structures
#####################################################

# model with only fixed effects (non-hierarchical)
model_FE = brm(formula = freq ~ gender * attitude, data = politedata)

# hierarchical model with random intercepts
model_interceptOnly = brm(formula = frequency ~ gender * attitude +
                            (1 | scenario + subject), 
                          data = politedata)

# hierarchical model with the maximial RE structure licensed by the design
# (notice that factor 'gender' does not vary for a given value of variable 'subject')
model_MaxRE = brm(formula = frequency ~ gender * attitude +
                    (1 + gender * attitude | scenario) +
                    (1 + attitude | subject), 
                  data = politedata)

#####################################################
## convenience function to extract predictor values &
## compute the probability of the relevant hypotheses
#####################################################

factors = list("gender" = c("F", "M"),
               "attitude" = c("inf", "pol"))

ref_levels_list = list("gender" = "F",
                  "attitude" = "inf")

extract_comparisons_generic = function(model, factors) {
  
  ref_levels = c("genderF", "attitudeinf")
  
  post_samples = posterior_samples(model) %>% select(starts_with("b_"))
  
  # get a table of cells (factor-level combinations)
  cells = expand.grid(factors) 
  for (j in 1:ncol(cells)) {
    levels(cells[,j]) = paste0(colnames(cells)[j], levels(cells[,j]))
  }
  
  coefficient_names = names(post_samples)
  
  # add the reference levels to the coefficient names
  for (c in 1:length(coefficient_names)) {
    for (f in names(factors)) {
      if (!grepl(f, coefficient_names[c])) {
        coefficient_names[c] = paste0(coefficient_names[c], "_", f, ref_levels_list[[f]])
      }
    }
  }
  names(post_samples) = coefficient_names
  
  # for each row in cells gather the columns in post_samples that contain the factor levels
  # and all other columns that contain any replacement of the current levels with the reference levels
  
  ## helper function that returns all coefficients to be added
  
  replace_with_ref_level_recursion = function(cell) {
    which_fcts_are_at_ref_level = map_lgl(cell, function(fl) fl %in% ref_levels)
    if (all(which_fcts_are_at_ref_level)) {
      return (list(cell))
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
  
  cell = c("genderF", "attitudeinf")
  
  replace_with_ref_level(cell)
  
  cells$out = "bla"
  
  for (i in 1:nrow(cells)) {
    coefficients_to_check = replace_with_ref_level(cells[i,1:(ncol(cells)-1)])
    # out = 0
    out = ""
    for (j in 1:length(coefficients_to_check)) {
      which_column = which(
        map_lgl(coefficient_names, function(coefficient_in_question) {
          all(map_lgl(coefficients_to_check[[j]], function(c) grepl(c, coefficient_in_question)))
        }) == T  
      )
      if(length(which_column) > 1) {stop("something went wrong; there have been several columns in the posterior samples that match where only one should")}
      # out = out + post_samples[coefficient_names[which_column]]
      out = paste0(out, "+", coefficient_names[which_column])
    }
    cells$out[i] = out
  }
  
}


extract_comparisons = function(model) {
  # get posterior samples
  post_samples = posterior_samples(model) %>% as.tibble()
  post_samples %>% select(everything(contains("genderM", "attitudepol")))
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



