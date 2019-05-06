#' Obtaining information about factors in regression model
#'
#' For a model for a factorial design, fitted with brms, this function returns information about the factors used, their levels, and the reference levels.
#' @param model Model fit from brms package.
#' @keywords regression, factorial design, brms
#' @export
#' @return list with names of factors and their levels, including the reference levels (in dummy coding)
#' @examples
#' library(brms)
#' m = brm(yield ~ N * P * K, npk)
#' get_factor_information(m)
get_factor_information = function(model) {
  
  # extract information about dependent and independent variables from formula
  ## TODO :: check extensively, especially for mixed effects models whether this stuff here works
  dependent_variable = as.character(formula(model)[[1]])[[2]]
  independent_variables = strsplit(x = gsub(pattern = "\\(.*\\|.*\\)", "", as.character(formula(model)[[1]])[[3]]),
                                   split =  "(\\*|\\+)",
                                   fixed = FALSE)[[1]] %>% trimws()
  independent_variables = independent_variables[which(independent_variables != "")]
  
  # stop this if there are not at least two factors
  if (length(independent_variables) < 1) {
    stop("Oeps! There do not seem to be any factors!")
  }
  
  # stop if any factor name or factor level contains a character that brms might not handle correctly
  check_permitted_characters = function(chr_vec) {
    str_which(chr_vec, "[^([:alnum:]|\\.|_)]")
  }
  # dependent variable
  if (length(check_permitted_characters(dependent_variable))>0) {
    stop("All variables, factor names, and factor levels must not contain any character except alpha-numeric characters (letters and numbers), dots '.', or underscores '_'.
The dependent variable '", dependent_variable, "' does not satisfy this constraint.")
  }
  # independent variables
  if (length(check_permitted_characters(independent_variables))>0) {
    stop("All variables, factor names, and factor levels must not contain any character except alpha-numeric characters (letters and numbers), dots '.', or underscores '_'.
The dependent variable(s) '", paste0(independent_variables[check_permitted_characters(independent_variables)], collapse = ", ") , "' does (do) not satisfy this constraint.")
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
  
  # check naming conventions in factor levels
  factor_levels = unlist(factors)
  if (length(check_permitted_characters(factor_levels))>0) {
    stop("All variables, factor names, and factor levels must not contain any character except alpha-numeric characters (letters and numbers), dots '.', or underscores '_'.
The factor level(s) '", paste0(factor_levels[check_permitted_characters(factor_levels)], collapse = ", ") , "' does (do) not satisfy this constraint.")
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
  
  return(list(
    dependent_variable = dependent_variable,
    independent_variables = independent_variables,
    factors = factors,
    ref_levels_list = ref_levels_list,
    ref_levels = ref_levels
  ))
}


#' Extracting cell means
#'
#' This function takes a brms model fit for a factorial design and outputs a comparison of all factor levels against each other, and posterior samples of all cell means. 
#' @param model Model fit from brms package.
#' @keywords regression, factorial design, brms
#' @export
#' @return list with (i) samples of estimated means of all cells and (ii) pairwise comparison of each cell (whether one has credibly a higher inferred mean than the other)
#' @examples
#' #' library(brms)
#' m = brm(yield ~ N * P * K, npk)
#' extract_posterior_cell_means(m)
extract_posterior_cell_means = function(model) {
  
  # get information about factors
  factor_info = get_factor_information(model)
  dependent_variable = factor_info[["dependent_variable"]]
  independent_variables = factor_info[["independent_variables"]]
  factors = factor_info[["factors"]]
  ref_levels_list = factor_info[["ref_levels_list"]]
  ref_levels = factor_info[["ref_levels"]]
  
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
        cell = paste(map_chr(1:ifelse(is.null(ncol(cell)), 1, ncol(cell)), 
                             function(j) as.character(cells_readable[i,j])), collapse = "__"),
        predictor_value = rowSums(post_samples[column_indices]),
        n_sample = 1:length(predictor_value)
      )
    }
  ) 
  
  predictor_values = predictor_values %>% spread(key = cell, value = predictor_value)
  
  ## an alternative (more versatile) output which compares all cells
  
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
  
  # borrowed from here: https://stackoverflow.com/questions/11693599/alternative-to-expand-grid-for-data-frames
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


#' Compare means of two subsets of factorial design cells
#'
#' This function takes a brms model fit for a factorial design and a specification of two groups (subsets of design cells) to compare. 
#' A group is specified as a named list, specifiying the factors and their levels which to include in the group.
#' It outputs the posterior mean of the 'higher' minus the 'lower' subset of cells, its 95 percent credible interval and the posterior probability that the 'higher' group has a higher mean than the the 'lower' group.
#' @param model Model fit from brms package.
#' @keywords regression, factorial design, brms
#' @importFrom HDInterval hdi
#' @export
#' @return list with posterior samples for each group, and the posterior probability that group 'higher' has a higher estimated coefficient in the posterior samples than the group 'lower'
#' @examples
#' library(brms)
#' m = brm(yield ~ N * P * K, npk)
#' # this compares two single cells in the factorial design
#' compare_groups(
#'  model = m, 
#'  higher = list("N" = "1", "P" = "1",  "K" = "1"), 
#'  lower  = list("N" = "0", "P" = "0",  "K" = "1")
#' )
#' # this compares the average of N=1 cells to the grand mean
#' # like in deviance conding
#' compare_groups(
#'  model = m, 
#'  higher = list("N" = "1"), 
#'  lower  = list()
#' )
#' 
compare_groups = function(model, higher, lower) {
  
  # get information about factors
  factor_info = get_factor_information(model)
  dependent_variable = factor_info[["dependent_variable"]]
  independent_variables = factor_info[["independent_variables"]]
  factors = factor_info[["factors"]]
  ref_levels_list = factor_info[["ref_levels_list"]]
  ref_levels = factor_info[["ref_levels"]]
  
  # check the input groups
  input_combined = c(higher, lower)
  ## check if all factor names specified are actually in the model
  input_factor_names = unique(names(input_combined))
  known_factor_names = map_lgl(input_factor_names, function(f) {f %in% independent_variables})
  if (sum(known_factor_names) < length(known_factor_names)) {
    stop("The following factor names specified in the groups to be compared do not match any independent variable in the specified model: 
  ", paste(input_factor_names[known_factor_names == F], collapse = ", "))
  }
  ## check if all factor levels specified are actually in the model
  for (i in length(input_combined)) {
    if (! input_combined[[i]] %in% factors[[names(input_combined)[i]]]) {
      stop("The level '", input_combined[[i]], "' is not part of the factor '", names(input_combined)[i], "' in the given model.")
    }
  }
  
  # get posterior samples for all cell means
  post_cell_samples = extract_posterior_cell_means(model)$predictor_values %>% 
    select(- n_sample)
  
  ## helper function :: recursive extraction of cell names
  collect_cell_names = function(remaining_names, remaining_factors) {
    if (length(remaining_factors) == 1) {
      return ( str_subset(remaining_names, remaining_factors) )
    } else {
      remaining_names = str_subset(remaining_names, remaining_factors[1]) 
      remaining_factors = remaining_factors[2:length(remaining_factors)]
      return(collect_cell_names(remaining_names, remaining_factors))
    }
  }
  
  ## helper function :: get names for cells
  get_group_names = function(group){
    if (length(group) == 0) {
      return("grand mean")
    }
    map_chr(1:length(group), 
            function(c) {paste(names(group[c]), unlist(group[c]), sep = ":")})  
  }
  
  ## helper function :: get means for each cell
  extract_group_samples = function(group) {
    if (length(group) == 0) {
      return(apply(as.matrix(post_cell_samples), 1, mean))
    }
    factors_group = get_group_names(group)
    cells_group = collect_cell_names(names(post_cell_samples), factors_group)
    apply(as.matrix(post_cell_samples %>% select(cells_group)), 1, mean)
  }
  
  post_samples_higher = extract_group_samples(higher)
  post_samples_lower  = extract_group_samples(lower)
  

  
  outlist = list(
    post_samples_higher = post_samples_higher,
    post_samples_lower = post_samples_lower,
    higher = get_group_names(higher),
    lower = get_group_names(lower),
    mean_diff = mean(post_samples_higher - post_samples_lower),
    l95_ci = as.vector(hdi(post_samples_higher - post_samples_lower)[1]),
    u95_ci = as.vector(hdi(post_samples_higher - post_samples_lower)[2]),
    probability = mean(post_samples_higher > post_samples_lower)
  )
  class(outlist) = "faintCompare"
  return(outlist)
}

#' Print comparison object between factor groups
#' @param model Model fit from brms package.
#' @keywords regression, factorial design, brms
#' @export
#' @return 
#' @examples
print.faintCompare = function(obj) {
  cat("Outcome of comparing groups:\n")
  cat(" * higher: ", obj$higher, "\n")
  cat(" * lower:  ", obj$lower, "\n")
  cat("Mean 'higher - lower': ", signif(obj$mean_diff, 4), "\n")
  cat("95% CI: [", signif(obj$l95_ci, 4), ";", signif(obj$u95_ci,4), "]\n")
  cat("P('higher - lower' > 0): ", signif(obj$probability,4), "\n")
}

# knit_print.faintCompare = function(obj) {
#   cat("Outcome of comparing groups:\n")
#   cat(" * higher: ", obj$higher, "\n")
#   cat(" * lower:  ", obj$lower, "\n")
#   cat("Mean 'higher - lower': ", signif(obj$mean_diff, 4), "\n")
#   cat("95% CI: [", signif(obj$l95_ci, 4), ";", signif(obj$u95_ci,4), "]\n")
#   cat("P('higher - lower' > 0): ", signif(obj$probability,4), "\n")
# }
