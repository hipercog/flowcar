# ---
# title: "Flow Short Scale pipeline"
# author: "Ben Cowley"
# ---
# attach packages
library(magrittr)
library(tidyverse)
library(data.table)
library(stringr)
library(purrr)
library(broom)
library(DT)
library(here)

source(file.path(here(), 'R', 'utils.R'))
# set folder with data files
# 
indir <- file.path(here(), '..', 'data')



### read Flow data ----
fss <-
  fread(file.path(indir, "fss_data_2019.csv")) %>%
  mutate(participant = participant + 9) %>%
  rbind(fread(file.path(indir, "fss_data_2017.csv")), .)
head(fss) # Flow + importance items Q1-Q13, skill-demand items A1-A3 (not used here)

# indicate fss questions for flow factors: fluency, absorption, flow, perceived importance
fluency_items <- c('Q2', 'Q4', 'Q5', 'Q7', 'Q8', 'Q9')
abs_items <- c('Q1', 'Q3', 'Q6', 'Q10')
flow_items <- c(fluency_items, abs_items)
pi_items <- c('Q11', 'Q12', 'Q13')
#pi_items <- c('Q11', 'Q12')


# compute flow factors per run
fss_items <- fss %>%
  mutate(fluency = rowMeans(dplyr::select(., all_of(fluency_items))),
         absorption = rowMeans(dplyr::select(., all_of(abs_items))),
         flow = rowMeans(dplyr::select(., all_of(flow_items))),
         pi_total = rowMeans(dplyr::select(., all_of(pi_items)))) %>%
  # add also total perceived importance (mean)
  rename(pi1 = Q11,
         pi2 = Q12,
         pi3 = Q13) %>%
  dplyr::select(participant:run, fluency, absorption, flow, pi1:pi3, pi_total) %>%
  # create ID = integer factor for participants
  mutate(ID = as.factor(as.integer(participant))) %>%
  # mutate participant variable into the right format (1 -> "01")
  mutate(participant = formatC(participant, width = 2, format = "d", flag = "0"))

head(fss_items)


