# ---
# title: "Flow data pipeline"
# author: "Ben Cowley"
# ---
# attach packages
library(RSQLite)
library(data.table)
library(magrittr)
library(tidyverse)
library(stringr)
library(purrr)
library(broom)
library(jsonlite)
library(DT)
library(lme4)
library(viridis)
library(here)
library(gghalves)


source(file.path(here(), 'R', 'utils.R'))
# set folder with data files
indir <- file.path(here(), '..', 'data')
# create an output dir for figures
FIGOUT <- FALSE
if (FIGOUT){
  odir <- file.path(here(), 'figures')
  dir.create(odir, showWarnings = FALSE)
}
# Also do something with blobs and steps (large tables)
BIGDATA <- FALSE
if (BIGDATA){
  t <- c("blob", "run", "step")
}else{
  t <- "run"
}


### Behavioral (CogCarSim) + FSS data ----
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
  
  rm(tb17, tb19)
}


# assign table to modify data
game_data <- tbl_runs %>%
  dplyr::select(-run) %>%
  # extract old participant variable into new participant, session, run variables
  tidyr::extract(participant, c('participant', 'session', 'run'),
                 "([[:alnum:]]+)-([[:alnum:]]+)-([[:alnum:]]+)") %>%
  mutate_at(vars(participant,session,run), as.numeric) %>%
  mutate(tbl_runs$provenance == 19, participant = participant + ((provenance - 17) * 4.5)) %>%
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

  # group_by(participant) %>%
  # mutate(mean_flow = mean(fss_items$flow)) %>%
  # ungroup() %>%
  # arrange(mean_flow)

head(fss_items)


# Basic violin plot
fss_items %>%
  mutate(participant = fct_reorder(participant, desc(flow), .fun='median')) %>%
  ggplot( aes(x=participant, y=flow)) + 
  geom_half_violin(trim=FALSE, fill="gray", side = "r") +
  geom_boxplot(width=0.1, outlier.shape = NA) +
  geom_hline(yintercept = median(fss_items$flow), color = "red") +
  geom_text(aes(length(unique(fss_items$participant)), median(fss_items$flow), label = "med.\nFlow", vjust = -1)) +
  labs(title="Flow scores per participant", x="Participant", y = "Flow") +
  ylim(2, 7) +
  theme_classic()
if (FIGOUT)
  ggsave(file.path(odir, "FlowXsubj.svg"))



# join game data and fss data into a single dataframe
fss_game <- game_data %>%
  dplyr::select(participant:run, collisions, duration, ln.duration, distance, cumrun, ln.cumrun) %>%
  left_join(fss_items, by = c('participant', 'session', 'run')) # if no vars given, use the ones with identical names

# display the datatable created just before
# datatable(fss_game, filter='top', options = list(pageLength = 5, scrollX=T) )

# see if everything is as it should be (e.g. no leftover rows of data, sane amount of sessions etc.)
summary(fss_game)



### Extracting learning curve coefficients and predicted learning curves
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
  select(participant, session, run, everything()) 




# spit out newly created data
head(fss_learning)

# plot linear performance
ggplot(fss_learning, aes(cumrun, duration)) +
  geom_point(alpha=.6, size=2) +
  geom_smooth(method = "lm", se = FALSE, linetype = 1, size = 0.5, color="red") +
  facet_wrap(~participant) +
  theme_bw(base_size = 14)
if (FIGOUT)
  ggsave(file.path(odir, "PerfXsubj.svg"))

# plot linear performance with power law fit
ggplot(fss_learning, aes(cumrun, duration)) +
  geom_point(alpha=.6, size=2) +
  stat_smooth(method = 'nls', formula = 'y~a*x^b', size = 0.5, se=FALSE, color ="red") +
  facet_wrap(~participant) +
  theme_bw(base_size = 14)
if (FIGOUT)  ggsave(file.path(odir, "PerfXsubj_powerlaw.svg"))

# acquire flow z-scores
fss_learning %<>% group_by(participant) %>% 
  mutate(z.flow = ((flow-mean(flow)) / sd(flow)))

# plot linear performance with power law fit and flow z-scores coloring
ggplot(fss_learning, aes(cumrun, duration, color = z.flow)) +
  geom_point(alpha=.6, size=2) +
  stat_smooth(method = 'nls', formula = 'y~a*x^b', size = 0.5, se=FALSE, color ="red") +
  facet_wrap(~participant) +
  scale_color_viridis() +
  theme_bw(base_size = 14)
if (FIGOUT)
  ggsave(file.path(odir, "PerfXsubj_powlxFlow.svg"))

# plot log-log with flow z-scores coloring
ggplot(fss_learning, aes(ln.cumrun, ln.duration, color = z.flow)) +
  geom_point(alpha=.6, size=2) +
  geom_smooth(method = "lm", se = FALSE, linetype = 1, size = 0.5, color="red") +
  facet_wrap(~participant) +
  scale_color_viridis() +
  theme_bw(base_size = 14)
if (FIGOUT)
  ggsave(file.path(odir, "PerfXsubj_powlxFlow_loglog.svg"))

# plot deviation from predicted curve
fss_learning <- mutate(fss_learning, deviation = ln.duration - learning_curve)
ggplot(fss_learning, aes(deviation, flow)) +
  geom_point(alpha=.6, size=2) +
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~participant) +
  theme_bw(base_size = 14)
if (FIGOUT)
  ggsave(file.path(odir, "FlowXdevXsubj.svg"))

# statistical test (linear mixed model)
fss_learning_lmer <- lmer(flow ~ deviation + (deviation|participant), data=fss_learning)
summary(fss_learning_lmer) # model summary
plot(fss_learning_lmer) # model diagnostics
qqnorm(residuals(fss_learning_lmer)) #qq-plot
qqline(residuals(fss_learning_lmer)) #line of "perfect normality"

