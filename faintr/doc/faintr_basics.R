## ----setup, include=FALSE, echo = FALSE, message = FALSE-----------------
knitr::opts_chunk$set(echo=TRUE, 
                      warning=FALSE,
                      message=FALSE, 
                      collapse = TRUE,
                      cache = TRUE,
                      dev.args = list(bg = 'transparent'), 
                      fig.align='center', 
                      fig.height = 3, 
                      fig.widht=4)
library(tidyverse)
theme_set(theme_bw() + theme(plot.background=element_blank()) )

## ---- eval = F-----------------------------------------------------------
#  devtools::install_github('michael-franke/bayes_mixed_regression_tutorial/faintr',
#                           build_vignettes = TRUE)
#  library(faintr)

## ---- echo = F-----------------------------------------------------------
library(faintr)

## ---- error=FALSE, warning=FALSE, message=FALSE--------------------------
library(tidyverse)
politedata = read_csv('https://raw.githubusercontent.com/michael-franke/bayes_mixed_regression_tutorial/master/code/politeness_data.csv') 
head(politedata)

## ------------------------------------------------------------------------
politedata %>% group_by(gender, context) %>% summarize(mean_pitch = mean(pitch))

## ---- error=FALSE, warning=FALSE, message=FALSE, results="hide"----------
library(brms)
m_dummy = brm(pitch ~ gender * context + (1 | subject + sentence), politedata)

## ------------------------------------------------------------------------
fixef(m_dummy)

## ------------------------------------------------------------------------
compare_groups(
  model = m_dummy, 
  higher = list(gender = "F", context = "pol"),
  lower = list(gender = "M", context = "inf")
)

## ------------------------------------------------------------------------
compare_groups(
  model = m_dummy, 
  higher = list(gender = "F"),
  lower = list()
)

## ------------------------------------------------------------------------
extract_posterior_cell_means(m_dummy)$all_cells_compared

## ------------------------------------------------------------------------
extract_posterior_cell_means(m_dummy)$cell_summary

