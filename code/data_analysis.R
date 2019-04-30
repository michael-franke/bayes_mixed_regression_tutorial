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

stop()

#####################################################
## run model with different random-effects structures
#####################################################

# model with only fixed effects (non-hierarchical)
model_FE = brm(formula = pitch ~ gender * context, data = politedata)

# model with only a fixed effect for gender (non-hierarchical)
model_gender = brm(formula = pitch ~ gender, data = politedata)

# hierarchical model with random intercepts
model_interceptOnly = brm(formula = pitch ~ gender * context +
                            (1 | sentence + subject),
                          data = politedata)

# hierarchical model with the maximial RE structure licensed by the design
# (notice that factor 'gender' does not vary for a given value of variable 'subject')
model_MaxRE = brm(formula = pitch ~ gender * context +
                    (1 + gender * context | sentence) +
                    (1 + context | subject),
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

  predictor_values = predictor_values %>% spread(key = cell, value = predictor_value)
  
  ## towards an alternative (more versatile) output which compares all cells

  cells = expand.grid(factors)
  for (j in 1:ncol(cells_readable)) {cells_readable[,j] = as.character(cells_readable[,j])}  
  cells$cell_name = map_chr(
    1:nrow(cells_readable), 
    function(i) {paste(as.character(cells_readable[i,]), collapse = "__")}
  )
  for (j in 1:ncol(cells)) {cells[,j] = as.character(cells[,j])}
  
  cells_high = cells
  cells_low = cells
  names(cells_high) = map_chr(names(cells), function(c) {paste0(c, "_high")})
  names(cells_low)  = map_chr(names(cells), function(c) {paste0(c, "_low")})
  
  # from here: https://stackoverflow.com/questions/11693599/alternative-to-expand-grid-for-data-frames
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
  
  all_cells_compared = expand.grid.df(cells_high, cells_low)
  for (j in 1:ncol(all_cells_compared)) {all_cells_compared[,j] = as.character(all_cells_compared[,j])}  
  all_cells_compared = all_cells_compared %>% filter(cell_name_high != cell_name_low)
  
  all_cells_compared$posterior = map_dbl(
    1:nrow(all_cells_compared), 
    function(i) {
      mean(predictor_values[all_cells_compared$cell_name_high[i]] > 
        predictor_values[all_cells_compared$cell_name_low[i]])
    }  
  )

  ## output
  
  
  
  return(list(
    predictor_values = predictor_values, 
    all_cells_compared = all_cells_compared))
  
}

get_cell_comparison = function(model, cell_low, cell_high) {
  ### TODO: currently relying on correct order of factors in 'cell_low' and 'cell_high'
  cell_high = paste(names(cell_high), unlist(cell_high), sep = ":", collapse = "__")
  cell_low = paste(names(cell_low), unlist(cell_low), sep = ":", collapse = "__")
  all_cells_compared = extract_posterior_cell_means(model)$all_cells_compared
  all_cells_compared %>% filter(cell_name_high == cell_high, cell_name_low == cell_low) %>% pull(posterior)
}

##################################
## testing selected hypotheses
#### two possibilities
##################################

## version 1 (old)
#### uses `extract_posterior_cell_means(model)$predictor_values`

get_posterior_beliefs_about_hypotheses = function(model) {
  posterior_cell_means = extract_posterior_cell_means(model)$predictor_values
  # insert the comparisons you are interested in as strings 
  tibble(hypothesis = c("Female-polite < Female-informal", 
                        "Male-polite < Male-informal",
                        "Male-informal < Female-polite"),
         probability = c(
           # insert the comparisons you are interested in referring to the extracted samples
           mean(posterior_cell_means$`gender:F__context:pol` < posterior_cell_means$`gender:F__context:inf`),
           mean(posterior_cell_means$`gender:M__context:pol` < posterior_cell_means$`gender:M__context:inf`),
           mean(posterior_cell_means$`gender:M__context:inf` < posterior_cell_means$`gender:F__context:pol`)
         ))
}

get_posterior_beliefs_about_hypotheses(model_FE)
get_posterior_beliefs_about_hypotheses(model_interceptOnly)
get_posterior_beliefs_about_hypotheses(model_MaxRE)

## version 2 (new)
### uses new convenience function `get_cell_comparison`

get_posterior_beliefs_about_hypotheses_new = function(model) {
  # insert the comparisons you are interested in as strings 
  tibble(hypothesis = c("Female-polite < Female-informal", 
                        "Male-polite < Male-informal",
                        "Male-informal < Female-polite"),
         probability = c(
           # insert the comparisons you are interested in referring to the extracted samples
           get_cell_comparison(model, list(gender = "F", context = "pol"), list(gender = "F", context = "inf")),
           get_cell_comparison(model, list(gender = "M", context = "pol"), list(gender = "M", context = "inf")),
           get_cell_comparison(model, list(gender = "M", context = "inf"), list(gender = "F", context = "pol")) 
         ))
}

get_posterior_beliefs_about_hypotheses_new(model_FE)
get_posterior_beliefs_about_hypotheses_new(model_interceptOnly)
get_posterior_beliefs_about_hypotheses_new(model_MaxRE)

