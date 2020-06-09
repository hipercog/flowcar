# ---
# title: "Flow data pipeline"
# author: "Ben Cowley (from T.Tammi)"
# ---
# attach packages
library(magrittr)
library(tidyverse)
library(data.table)
library(stringr)
library(purrr)
library(broom)
library(jsonlite)
library(DT)
library(here)

source(file.path(here(), 'R', 'utils.R'))
# set folder with data files
indir <- file.path(here(), '..', 'data')
# Also do something with blobs and steps (large tables)
BIGDATA <- FALSE
if (BIGDATA){
  t <- c("blob", "run", "step")
}else{
  t <- "run"
}


### Behavioral (CogCarSim) data ----
# import the cogcarsim2 databases
tb17 <- readCCSdb(file.path(indir, "cogcarsim2_2017.db"), t, 17)
tb19 <- readCCSdb(file.path(indir, "cogcarsim2_2019.db"), t, 19)

tbl_runs <- rbind(tb17$tbl_run, tb19$tbl_run)
head(tbl_runs)

# remove extra (test) runs
tbl_runs %<>% # %<>% operator overwrites dataframe
  filter(startsWith(participant, '0'))

if (BIGDATA){
  tbl_blob <- rbind(tb17$tbl_blob, tb19$tbl_blob)
  head(tbl_blob)
  tbl_step <- rbind(tb17$tbl_step, tb19$tbl_step)
  head(tbl_step)

  # filter other tables with updated runs
  tbl_blob %<>%
    filter(run %in% tbl_runs$run)
  
  tbl_step %<>%
    filter(run %in% tbl_runs$run)
  
}
rm(tb17, tb19)


# clean up / modify behavioural data
game_data <- tbl_runs %>%
  dplyr::select(-run) %>%
  # extract old participant variable into new participant, session, run variables
  tidyr::extract(participant, c('participant', 'session', 'run'),
                 "([[:alnum:]]+)-([[:alnum:]]+)-([[:alnum:]]+)") %>%
  mutate_at(vars(participant,session,run), as.numeric) %>%
  # this is a pure hack to change 2019 participant numbers to be 10:18
  mutate(participant = participant + ((provenance - 17) * 4.5)) %>%
  # arrange rows to create cumulative run variable
  arrange(participant, session, run) %>%
  # mutate participant variable into the right format (1 -> "01")
  mutate(participant = formatC(participant, width = 2, format = "d", flag = "0")) %>%
  group_by(participant) %>%
  # also log(duration) and log(cumrun) variables for plotting
  mutate(cumrun = row_number(),
         ln.duration = log(duration),
         ln.cumrun = log(cumrun)) %>% 
  ungroup # remember to ungroup!

head(game_data)



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
  # mutate participant variable into the right format (1 -> "01")
  mutate(participant = formatC(participant, width = 2, format = "d", flag = "0"))

head(fss_items)


# join game data and fss data into a single dataframe
fss_game <- game_data %>%
  dplyr::select(participant:run, collisions, duration, ln.duration, distance, cumrun, ln.cumrun) %>%
  left_join(fss_items, by = c('participant', 'session', 'run')) # if no vars given, use the ones with identical names

# see if everything is as it should be (e.g. no leftover rows of data, sane amount of sessions etc.)
summary(fss_game)



### Extracting learning curve coefficients and predicted learning curves ----
# Fit a log-log regression curve for each participant separately: 
#     log(duration) ~ log(cumulative run)
# to get slope and intercept coefficients
fss_learning <- fss_game %>%
  group_by(participant) %>%
  nest() %>% 
  # fit a linear model using log-log formula on durations and cumruns per participant
  mutate(fit = map(data, ~lm(ln.duration ~ ln.cumrun, data = .)), # make models and...
       coef = map(fit, tidy), # ...get coefficients and...
       data_aug = map(fit, augment)) %>% #...get predictions
  unnest(coef) %>%
  # predicted curve slope and interception
  select(participant, term, estimate, data_aug) %>%
  spread(term, estimate) %>%
  rename(slope = ln.cumrun,
       intercept = `(Intercept)`) %>%
  unnest(data_aug) %>%
  ungroup %>%
  # join new ones to fss_game -> fss_learning
  select(participant:.fitted, intercept, slope) %>%
  rename(learning_curve = .fitted) %>%
  left_join(fss_game, by = c('participant', 'ln.duration', 'ln.cumrun')) %>%
  select(participant, session, run, everything()) %>%
  # acquire  deviation from predicted curve
  mutate(deviation = ln.duration - learning_curve) %>%
  group_by(participant) %>%
  # acquire flow z-scores
  mutate(z.flow = ((flow-mean(flow)) / sd(flow))) %>%
  ungroup()

# spit out newly created data
head(fss_learning)

