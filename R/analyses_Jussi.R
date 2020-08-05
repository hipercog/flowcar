#Cogcarsim #2 flow data analyses (c) Jussi Palomäki

#Load libraries
library(ggplot2)
library(gridBase)
library(colorRamps)  
library(grDevices)
library(lme4)
library(effects)
library(multcomp)
library(tidyverse)
library(lmerTest)
library(ggpubr)
library(MuMIn)
library(lattice)
library(robustlmm)
library(grid)
library(sjPlot)
library(dynlm)
library(viridis)
library(psych)
library(lavaan)
library(GPArotation)
library(semTools)
library(lavaanPlot)
library(lavaan)
library(jtools)
library(interactions)
library(XML)
library(semPlot)
library(GGally)
library(ggbiplot)
library(gghalves)
library(bmlm) # see https://link.springer.com/article/10.3758/s13428-017-0980-9#Sec10


#Remove first row (cumrun == 1) for each participant. NOT USED IN THE INTERACTION ANALYSIS.
fss_learning_sub <- subset(fss_learning, cumrun!=1)
anova(lmer(flow ~ deviation + (deviation|participant), data=fss_learning_sub))
anova(lmer(flow ~ deviation + (deviation|participant), data=fss_learning))

#RQ2: Does the effect of Deviation on Flow (main finding from previous paper) depend on level cumrun (how many runs the participant has played)

#Fit model (note: could also try random slope/intercept for deviation*cumrun):
fss_learning2_lmer <- lmer(flow ~ deviation*cumrun + (deviation*cumrun|participant), data=fss_learning) #consider dropping random slope entirely; it won't really affect the results at all.

#Simple slopes:
sim_slopes(fss_learning2_lmer, pred = deviation, modx = cumrun, johnson_neyman = TRUE)

#Model R^2:
r.squaredGLMM(fss_learning2_lmer)
summ(fss_learning2_lmer)

#Test model assumptions:
plot(fss_learning2_lmer)
plot(resid(fss_learning2_lmer), fss_learning$deviation)
plot(resid(fss_learning2_lmer), fss_learning$cumrun)
qqnorm(residuals(fss_learning2_lmer))
qqline(residuals(fss_learning2_lmer))
qqmath(ranef(fss_learning2_lmer, condVar = TRUE), strip = TRUE)

#Model assumptions suggest the model is "not great, not terrible" (i.e. "3.6").
#Fit model using robust lmm:
fss_learning2_rlmer <- rlmer(flow ~ deviation*cumrun + (deviation|participant), data=fss_learning)
summary(fss_learning2_rlmer)
plot(fss_learning2_rlmer)

#Plot interaction:
plot_model(fss_learning2_lmer, type = "pred", terms = c("deviation", "cumrun"), mdrt.values = "meansd")

#Double check by plotting interaction by hand
#Calculate -1, mean, and +1 SD for cumrun
cumrun_minusSD <- mean(fss_learning$cumrun) - sd(fss_learning$cumrun)
cumrun_meanSD <- mean(fss_learning$cumrun)
cumrun_plusSD <- mean(fss_learning$cumrun) + sd(fss_learning$cumrun)

flow_interaction <- effect(c("deviation*cumrun"), fss_learning2_lmer,
                   xlevels=list(deviation=c(-0.2, 0, 0.2), #x-axis range
                                cumrun=c(cumrun_minusSD, cumrun_meanSD, cumrun_plusSD),
                                se=TRUE, confidence.level=.95, typical=mean))

flow_interaction <- as.data.frame(flow_interaction)
flow_interaction$cumrun <- as.factor(flow_interaction$cumrun)
ggplot(flow_interaction, aes(x=deviation, y=fit, colour=cumrun)) + geom_line()

#For participant-wise visualization, calculate slopes and intercepts of deviation at various levels of moderator (cumrun) for each participant.
#Here we're essentially just fitting LMs to each participant separately.
slope_intercept_dataframes = list()
fss_models <- list()
for (i in 1:length(unique(fss_learning$participant))) {
  temp_model = lm(flow ~ deviation*cumrun, data=subset(fss_learning, ID==i))
  flow_interaction <- effect(c("deviation*cumrun"), temp_model, 
                             xlevels=list(deviation=c(-0.2, 0, 0.2), #arbitrary levels, 0 to get intercept
                                          cumrun=c(cumrun_minusSD, cumrun_meanSD, cumrun_plusSD), #values are -1SD and +1SD
                                          se=TRUE, confidence.level=.95, typical=mean))
  flow_interaction <- as.data.frame(flow_interaction)
  flow_interaction$ID <- as.factor(i)
  flow_interaction$cumrun <- as.factor(flow_interaction$cumrun)
  #flow_interaction$ID <- as.factor(flow_interaction$ID)
  flow_interaction <- flow_interaction %>% mutate(slope_minusSD = fit[3]-fit[1], #slope is the difference in FIT between 1 unit increment in deviation
                                                  slope_MEAN = fit[6]-fit[4],
                                                  slope_plusSD = fit[9]-fit[7],
                                                  intercept_minusSD = fit[2], #intercept is the value of fit when deviation = 0
                                                  intercept_MEAN = fit[5],
                                                  intercept_plusSD = fit[8])
  slope_intercept_dataframes = rbind(slope_intercept_dataframes, flow_interaction)
  fss_models[[i]] = temp_model
}

#Grab p-values, sort from lowest to highest, and adjust using Bonferroni-Holm
fss_Pvalues <- list()
for (i in 1:18) {
  fss_Pvalues[[i]] <- anova(fss_models[[i]])$"Pr(>F)"[3]
}
fss_Pvalues <- as.data.frame(fss_Pvalues) 
colnames(fss_Pvalues) <- c(1:18)
fss_Pvalues <- sort(fss_Pvalues)
fss_Pvalues_holm <- p.adjust(fss_Pvalues, "holm")

#Grab order of participants for plotting
part_levels <- as.integer(names(fss_Pvalues))

#Reorder levels in data; now according to p-value, from low (significant) to high (nonsignificant):
fss_learning$ID <- factor(fss_learning$ID, levels=part_levels)

#Start plotting results. First assign initial plot as flow_fig
flow_fig <- ggplot(fss_learning, aes(deviation, flow)) +
  geom_point(alpha=.2, size=1.5) +
  #geom_smooth(method = "lm", se = FALSE, fullrange=TRUE) + 
  facet_wrap(~ID) +
  geom_abline(data = slope_intercept_dataframes, aes(intercept=intercept_minusSD, slope=slope_minusSD, linetype="-1 SD (= 8.81 runs)"), size=0.8) + 
  geom_abline(data = slope_intercept_dataframes, aes(intercept=intercept_plusSD, slope=slope_plusSD, linetype="+1 SD (= 32.19 runs)"), size=0.8) + 
  labs(colour = NULL, linetype = "Cumulative runs", title = NULL) + xlab("Deviation score") + ylab("Flow score") +
  theme_bw(base_size=14) + 
  theme(legend.position = c(x=.8, y=.1),
        legend.background = element_rect(fill="white"),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(size=9),
        panel.background = element_rect(fill = "white",
                                        colour = "lightgrey",
                                        size = 0.5, linetype = "solid"),
        panel.border = element_blank(), 
        panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey"),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="lightgray"))


#Function for plotting. (Custom colors for facets couldn't be done using ggplot() directly)
#COLOR THEMES, please change your desired alpha level (below which results are "statistically significant")
#INTERPRETATION: 
#Light green panels = interaction deviation*cumrum is statistically significant after bonferroni-holm correction for 18 multiple comparisons
#Dark green panels = interaction deviation*cumrum is statistically significant without correcting for multiple comparisons
#Grey panels = deviation*cumrum is not statistically significant

panel_colors <- function(alpha_level) {
  
  holm_adjusted_values <- rep("yellowgreen", sum(fss_Pvalues_holm < alpha_level))
  non_adjusted_values <- rep("yellow4", sum(fss_Pvalues < alpha_level)-sum(fss_Pvalues_holm < alpha_level)) #overlap between holm_adjusted removed!
  non_significant_values <- rep("grey69", 18-sum(fss_Pvalues < alpha_level)+sum(fss_Pvalues_holm < alpha_level)) #overlap between holm_adjusted and non_adjusted removed!
  all_values <- c(holm_adjusted_values,non_adjusted_values,non_significant_values)
  
  g <- ggplot_gtable(ggplot_build(flow_fig))
  stripr <- which(grepl('strip-t', g$layout$name))
  stripr <- rev(stripr[-c(4:5)])
  stripr <- c(rev(stripr[1:5]), rev(stripr[6:10]), rev(stripr[11:15]), rev(stripr[16:18]))
  fills <- all_values
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  
  #Draw figure
  grid.draw(g)
    
}


panel_colors(.05) #these functions may take a bit of time to load so be patient!
panel_colors(.1)


#SLIDING WINDOW (ACROSS PARTICIPANTS) window, calculates a sliding window of models fitted for N cumruns
#Now also calculates partially pooled results for participant-wise visualisation

#An example would be d(y) ~ L(y, 2), where d(x, k) is diff(x, lag = k) and L(x, k) is lag(x, lag = -k), note the difference in sign. The default for k is in both cases 1.


sliding_window_across <- function(window_width = 10, result=3) {
  estimates <- list()
  estimates_IDwise <- list()
  for (i in 1:(length(unique(fss_learning$cumrun))-window_width)) {
    fss_learning_temp <- subset(fss_learning, cumrun >= i & cumrun < i+window_width)
    
    ###calculate deviation for each segment
    fss_learning_temp <- fss_learning_temp %>% group_by(participant) %>%
      mutate(learning_curve = predict(lm(ln.duration ~ ln.cumrun)), deviation = ln.duration - learning_curve)
    ###
  
    fss_lmer_temp <- lmer(flow ~ deviation + (deviation|participant), 
                          data=fss_learning_temp)
    estimates[[i]] <- coef(summary(fss_lmer_temp))[2,1] #fixed effects slope, i.e. b-value for whole model
    estimates_IDwise[[i]] <-ranef(fss_lmer_temp)[["participant"]][2] #slopes for participants
  }
  estimates <- as.data.frame(unlist(estimates))
  estimates_IDwise <- as.data.frame(unlist(estimates_IDwise))
  estimates <- estimates %>% mutate(sliding_window=1:(length(unique(fss_learning$cumrun))-window_width))
  estimates_IDwise <- estimates_IDwise %>% mutate(ID=rep(1:18, length(unique(fss_learning$cumrun))-window_width),
                                                  sliding_window=rep(1:(length(unique(fss_learning$cumrun))-window_width), 
                                                                     each=18))
  names(estimates)[1] <- "estimate"
  names(estimates_IDwise)[1] <- "estimate"
  permuted_estimates <- transform(estimates, estimate = sample(estimate))
  
  #Choose which result to return (1=partially pooled figures, participant-wise, 2=figure across participants, 3=permuted estimates, 4=observed estimates):

  #FIGURE: participant-wise, partially pooled
  if (result == 1) {
    result <- ggplot(estimates_IDwise, aes(x=sliding_window, y=estimate)) + geom_line() + 
      geom_smooth(method="lm", size=0.3) + facet_wrap("ID") + labs(title = "Partial pooling, LMER model")
  }
  
  #FIGURE: across participants
  else if (result == 2) {
    result <- ggplot(estimates, aes(x=sliding_window, y=estimate)) + geom_line() + 
      geom_smooth(method="lm", size=0.3)
    
  #DATA VECTOR: permuted estimates
  }
  else if (result == 3) {
    result <- permuted_estimates
  }
  
  #DATA VECTOR: observed estimates
  else {
    result <- estimates
  }
  return(result)
}

#SLIDING WINDOW (PARTICIPANT-WISE)
sliding_window_within <- function(window_width = 10, result=2) {
  estimates_IDwise <- list()
  counter = 0
  for (x in 1:length(unique(fss_learning$participant))) {
    for (i in 1:(length(unique(fss_learning$cumrun))-window_width)) {
      fss_learning_temp <- subset(fss_learning, ID==x & cumrun >= i & cumrun < i+window_width)
      
      ###calculate deviation again for each segment
      fss_learning_temp <- fss_learning_temp %>% group_by(participant) %>%
        mutate(learning_curve = predict(lm(ln.duration ~ ln.cumrun)), deviation = ln.duration - learning_curve)
      ###
      
      temp_model <- lm(flow ~ deviation, data=fss_learning_temp)
      counter = counter + 1
      estimates_IDwise[[counter]] <- coef(summary(temp_model))[2,1]
    }
  }
  
  estimates_IDwise <- as.data.frame(unlist(estimates_IDwise)) 
  estimates_IDwise <- estimates_IDwise %>% mutate(sliding_window=rep(c(1:(length(unique(fss_learning$cumrun))-window_width)), 18), 
                                                  ID=rep(1:18, each=length(unique(fss_learning$cumrun))-window_width)) #18 is number of participants
  names(estimates_IDwise)[1] <- "estimate"
  permuted_estimates <- estimates_IDwise %>% group_by(ID) %>% 
    transform(estimate = sample(estimate))
  
  #choose which result to return (1=Participant-wise figure, 2=DATA VECTOR: permuted estimates, 3=DATA VECTOR: observed estimates):
  #participant-wise
  if (result == 1) {
    result <- ggplot(estimates_IDwise, aes(x=sliding_window, y=estimate)) + geom_line() + 
      ylab("Slope estimate") + xlab("Sliding window") +
      theme_bw() +
      facet_wrap("ID") + geom_smooth(method="lm", se=FALSE, size=0.3) + labs(title = NULL)
  }
  
  #permuted estimates:
  else if (result == 2) {
    result <- permuted_estimates
  }
  
  #observed estimates:
  else {
    result <- estimates_IDwise
  }
  return(result)
}

#PERMUTATIONS
permutations <- lapply(rep(10,30), sliding_window_within)
temp_list <- list()
for (i in 1:length(permutations)) {
  temp_list[[i]] <- permutations[[i]][1]
}
temp_frame <- temp_list %>% as.data.frame() #%>% rowMeans() #rowmeans being the mean of all permutations old.

#temp_frame <- temp_frame %>% mutate(ID = rep(1:18, each=30))


#new attempt:subtract each permutation from observed, then store
difference_list <- list()
for (i in 1:(ncol(temp_frame))) {
  difference_list[[i]] <- sliding_window_within(10,3)[1] - temp_frame[i]
}
difference_list <- difference_list %>% as.data.frame()
#transpose(difference_list[c(1:30),]) #this is one participant! (30 columns, rows = permutations) feed this to findcurves for each participant!

#OLD
# temp_frame <- temp_frame %>% as.data.frame() %>% mutate(mean_z = scale(.), ID = rep(c(1:18), each=30)) #need to be manually changed! 1:18 = participans, 30 = number of segments
# difference <- sliding_window_within(10,3)[1] - temp_frame[1] #consider scaling the sliding_window
# difference <- difference %>% mutate(ID = rep(c(1:18), each=30), sliding_window = rep(c(1:30), 18)) #these need to be manually changed! see above
# #names(difference)[1] <- "estimate" #not necessary if scale() not used for sliding_window_within() above
# ggplot(difference, aes(x=sliding_window, y=estimate)) + geom_line() + geom_smooth(method="lm", size=0.3) + facet_wrap("ID")

#MWEs
#difference measure
# difference_wide <- difference %>% pivot_wider(names_from = sliding_window, values_from = estimate) %>%
#   dplyr::select(-ID)
curve <- findcurves(transpose(difference_list[c(1:30),])) #this is participant 1
idx = c(1:30)
rng <- range(curve[,c("mean0","lo","up")])
plot(c(idx[1],idx[length(idx)]), rng
     , type="n", bty="n", xlab="sliding window", ylab="estimate", main = "Just a test")
plotlines(curve, idx, col ="blue")

#MWE of just the actual observations
estimates_IDwise_wide <- sliding_window_within(10,3) %>% pivot_wider(names_from = sliding_window, values_from = estimate) %>% dplyr::select(-ID)
curve <- findcurves(estimates_IDwise_wide)
idx = c(1:30)
rng <- range(curve[,c("mean0","lo","up")])
plot(c(idx[1],idx[length(idx)]), rng-0.5
     , type="n", bty="n", xlab="Sliding window", ylab="Slope estimate", main = NULL)
plotlines(curve, idx, col ="blue")


#OLD
# #random permutations
# random_permutations <- temp_frame[1]
# random_permutations <- random_permutations %>% mutate(ID = rep(c(1:18), each=30), sliding_window = rep(c(1:30), 18))
# names(random_permutations)[1] <- "estimate"
# permutations_wide <- random_permutations %>% pivot_wider(names_from = sliding_window, values_from = estimate) %>% 
#   dplyr::select(-ID)
# curve2 <- findcurves(permutations_wide)
# plotlines(curve2, idx, col = "red")
# 
#actual observations (same as above so redundant!)
# actual_observations <- sliding_window_within(10,3)[1]
# actual_observations <- actual_observations %>% mutate(ID = rep(c(1:18), each=30), sliding_window = rep(c(1:30), 18))
# #names(actual_observations)[1] <- "estimate"
# observations_wide <- actual_observations %>% pivot_wider(names_from = sliding_window, values_from = estimate) %>%
#   dplyr::select(-ID)
# curve3 <- findcurves(observations_wide)
# idx = c(1:30)
# rng <- range(curve3[,c("mean0","lo","up")])
# plot(c(idx[1],idx[length(idx)]), rng
#      , type="n", bty="n", xlab="sliding window", ylab="estimate", main = "Just a test")
# plotlines(curve3, idx, col = "blue")

#SLIDING WINDOW, ACROSS PARTICIPANTS, fixed increments of 5 (session length), 8 segments (1...5, 6...10, 11...15, etc.)
sliding_segment_across <- function(segment=5, result=2) {
  estimates <- list()
  estimates_IDwise <- list()
  for (i in seq(1, length(unique(fss_learning$cumrun)), segment)) {
    fss_learning_temp <- subset(fss_learning, cumrun >= i & cumrun < i+segment)
  
  ###calculate deviation again for each segment
    fss_learning_temp <- fss_learning_temp %>% group_by(participant) %>%
      mutate(learning_curve = predict(lm(ln.duration ~ ln.cumrun)), deviation = ln.duration - learning_curve)
  ###
  
    fss_lmer_temp <- lmer(flow ~ deviation + (deviation|participant), 
                          data=fss_learning_temp)
    estimates[[i]] <- coef(summary(fss_lmer_temp))[2,1] #slope, i.e. b-value (this results in list with NULL values due to "i" (1, 6, etc.)
    estimates_IDwise[[i]] <-ranef(fss_lmer_temp)[["participant"]][2] #slopes for participants
  
  }
  estimates <- as.data.frame(unlist(estimates)) #deal with NULL values
  estimates_IDwise <- as.data.frame(unlist(estimates_IDwise))
  estimates <- estimates %>% mutate(segment=1:(length(unique(fss_learning$cumrun))/segment))
  estimates_IDwise <- estimates_IDwise %>% mutate(ID=rep(1:18, length(unique(fss_learning$cumrun))/segment),
                                                  segment=rep(1:(length(unique(fss_learning$cumrun))/segment), each=18))

  names(estimates)[1] <- "estimate"
  names(estimates_IDwise)[1] <- "estimate"
  permuted_estimates <- transform(estimates, estimate = sample(estimate))

  #FIGURE, across participants
  if (result == 1) {
    result <- ggplot(estimates, aes(x=segment, y=estimate)) + geom_line() + geom_smooth(method="lm", size=0.3)
  }
  
  #DATA VECTOR: permuted estimates
  else if (result == 2) {
    result <- permuted_estimates
  }
  
  #DATA VECTOR: observed estimates
  else {
    result <- estimates
  }
  return(result)

}

#OLD
# #PERMUTATIONS
# permutations <- lapply(rep(5,20), sliding_segment_across)
# temp_list <- list()
# for (i in 1:length(permutations)) {
#   temp_list[[i]] <- permutations[[i]][1]
#   temp_frame <- temp_list %>% as.data.frame() %>% rowMeans()
#   temp_frame <- temp_frame %>% as.data.frame() %>% mutate(mean_z = scale(.))
# }
# 
# difference <- sliding_segment_across(5,3)[1] - temp_frame[1]
# difference <- difference %>% mutate(sliding_window = seq(1:8))
# ggplot(difference, aes(x=sliding_window, y=estimate)) + geom_line() + geom_smooth(method="lm", size=0.3)


#SLIDING WINDOW, PARTICIPANT-WISE, fixed increments of 5 (session length), 8 segments (1...5, 6...10, 11...15, etc.)
sliding_segment_within <- function(segment = 5, result=2) {
  estimates_IDwise <- list()
  counter = 0
  for (x in 1:length(unique(fss_learning$participant))) {
    for (i in seq(1, length(unique(fss_learning$cumrun)), 5)) {
      fss_learning_temp <- subset(fss_learning, ID==x & cumrun >= i & cumrun < i+5)
      
      ###calculate deviation again for each segment
      fss_learning_temp <- fss_learning_temp %>% group_by(participant) %>%
        mutate(learning_curve = predict(lm(ln.duration ~ ln.cumrun)), deviation = ln.duration - learning_curve)
      ###
      
      temp_model = lm(flow ~ deviation, data=fss_learning_temp)
      counter = counter + 1
      estimates_IDwise[[counter]] <- coef(summary(temp_model))[2,1]
    }
  }
  estimates_IDwise <- as.data.frame(unlist(estimates_IDwise)) 
  estimates_IDwise <- estimates_IDwise %>% mutate(segment=rep(c(1:8), 18), 
                                                  ID=rep(1:18, each=8))
  names(estimates_IDwise)[1] <- "estimate"

  permuted_estimates <- estimates_IDwise %>% group_by(ID) %>% 
    transform(estimate = sample(estimate))

  if (result == 1) {
    result <- ggplot(estimates_IDwise, aes(x=segment, y=estimate)) + geom_line() + 
      facet_wrap("ID") + geom_smooth(method="lm", size=0.3) + labs(title = "Individual level LMs")
  }
  
  #permuted estimates:
  else if (result == 2) {
    result <- permuted_estimates
  }
  
  #observed estimates:
  else {
    result <- estimates_IDwise
  }
  return(result)
}

#OLD
# #PERMUTATIONS
# permutations <- lapply(rep(5,20), sliding_segment_within) #figure out how to input function arguments directly here, now default = 2
# temp_list <- list()
# for (i in 1:length(permutations)) {
#   temp_list[[i]] <- permutations[[i]][1]
#   temp_frame <- temp_list %>% as.data.frame() %>% rowMeans()
#   temp_frame <- temp_frame %>% as.data.frame() %>% mutate(mean_z = scale(.), ID = rep(c(1:18), each=8)) #need to be manually hanged! 1:18 = participans, 8 = number of segments
# }
# 
# difference <- sliding_segment_within(5,3)[1] - temp_frame[1] #estimates vs permuted estimates
# difference <- difference %>% mutate(ID = rep(c(1:18), each=8), sliding_window = rep(c(1:8), 18)) #these need to be manually changed! see above
# ggplot(difference, aes(x=sliding_window, y=estimate)) + geom_line() + geom_smooth(method="lm", size=0.3) + facet_wrap("ID")



###OLDER VERSION:
#Across participants (calculates b-value for models with 4,5,6...40 cumruns):
estimates <- list()
for (i in 4:length(unique(fss_learning$cumrun))) {
  
  fss_learning_temp <- subset(fss_learning, cumrun<=i)
  fss_lmer_temp <- lmer(flow ~ deviation + (deviation|participant), 
                        data=fss_learning_temp)
  estimates[[i]] <- coef(summary(fss_lmer_temp))[2,1]
  
}

estimates <- as.data.frame(unlist(estimates)) 
estimates <- estimates %>% mutate(cumrun=4:40) 
names(estimates)[1] <- "estimate"
across_participants <-ggplot(estimates, aes(x=cumrun, y=estimate)) + geom_line() +
  geom_smooth(method="lm", size=0.3)


#Participant-wise (calculates b-value for models with 4,5,6...40 cumruns):
estimates_IDwise <- list()
counter = 0
for (x in 1:length(unique(fss_learning$participant))) {
  for (i in 4:length(unique(fss_learning$cumrun))) {
    temp_model = lm(flow ~ deviation, data=subset(fss_learning, ID==x & cumrun<=i))
    counter = counter + 1
    estimates_IDwise[[counter]] <- coef(summary(temp_model))[2,1]
  }
}

estimates_IDwise <- as.data.frame(unlist(estimates_IDwise)) 
estimates_IDwise <- estimates_IDwise %>% mutate(cumrun=rep(c(4:40), 18), 
                                                ID=rep(1:18, each=37))
names(estimates_IDwise)[1] <- "estimate"
participant_wise <- ggplot(estimates_IDwise, aes(x=cumrun, y=estimate)) + geom_line() + 
  facet_wrap("ID") + geom_smooth(method="lm", size=0.3)

ggarrange(across_participants, participant_wise, 
          font.label = list(size = 10, color = "black", face = "bold"), 
          common.legend = TRUE, legend="bottom", hjust = -0.25, ncol = 2, nrow = 1)


######
#LAGGED REGRESSIONS:

lag_fss_learning <- fss_learning %>% group_by(participant) %>% mutate(lag_flow1 = lag(flow, 1), lag_flow2 = lag(flow, 2), lag_flow3 = lag(flow, 3), 
                                                                      lag_flow4 = lag(flow, 4), lag_flow5 = lag(flow, 5),
                                                                      lag_deviation1 = lag(deviation, 1), lag_deviation2 = lag(deviation, 2), lag_deviation3 = lag(deviation, 3),
                                                                      lag_deviation4 = lag(deviation, 4), lag_deviation5 = lag(deviation, 5),
                                                                      lag_duration1 = lag(ln.duration, 1), lag_duration2 = lag(ln.duration, 2),
                                                                      lag_duration3 = lag(ln.duration, 3), lag_duration4 = lag(ln.duration, 4),
                                                                      duration_dif = lag(ln.duration, 1) - ln.duration)


#effect of trial N duration - trial N-1 duration on current trial flow. So: SDT predicts that feedback influences flow, but why is there an intereaction with cumrun? That is, feedback 
#influences flow MORE the more players have played. Alternative explanation: flow is in fact predicted by deviation from expected skill/performance!

summary(lmer(flow ~ duration_dif*ln.cumrun + (duration_dif|participant), data=lag_fss_learning)) #Here you can also add *deviation to see which of the two interactions, cumrun*deviation or cumrun*duration_dif is more important
#summary(lmer(flow ~ duration_dif + deviation + (deviation+duration_dif|participant), data=lag_fss_learning)) #double check syntax for variability in TWO slopes
summary(lmer(flow ~ lag_flow1 + lag_flow2 + lag_flow3 + lag_flow4 + (1|participant), data=lag_fss_learning))

lag_fss_learning_longer <- lag_fss_learning %>% dplyr::select(ID, flow, lag_flow1, lag_flow2, lag_flow3, lag_flow4) %>%
  gather(key=lag_amount, value=lagged_flow_score, lag_flow1, lag_flow2, lag_flow3, lag_flow4)

#You need to use identical name and labels values for both shape and colour scale.
ggplot(lag_fss_learning_longer, aes(lagged_flow_score, flow, size=lag_amount, linetype=lag_amount)) + 
  geom_smooth(method = "lm", fullrange=T, alpha=.15) + 
  geom_point(alpha=.2) +
  facet_wrap("participant") +
  scale_size_discrete(range=c(0.8, 0.2), labels=c("Flow lag = 1", "Flow lag = 2", "Flow lag = 3", "Flow lag = 4")) + 
  scale_linetype_discrete(labels=c("Flow lag = 1", "Flow lag = 2", "Flow lag = 3", "Flow lag = 4")) +
  theme_bw() + xlab("Lagged flow score") + ylab("Flow score") +
  guides(size=guide_legend(title="Lag amount"), linetype=guide_legend(title="Lag amount"))


#ALTERNATIVE #1:
ggplot(lag_fss_learning_longer, aes(lagged_flow_score, flow)) + 
  geom_smooth(method = "lm", fullrange=T, alpha=.15) + 
  geom_point(alpha=.2) +
  facet_grid(lag_amount~participant) +
  scale_linetype_discrete(labels=c("Flow lag = 1", "Flow lag = 2", "Flow lag = 3", "Flow lag = 4")) +
  theme_bw() + xlab("Lagged flow score") + ylab("Flow score")


#ALTERNATIVE #2:
ggplot(lag_fss_learning_longer, aes(lagged_flow_score, flow, color=lag_amount, size=lag_amount)) + 
  geom_smooth(method = "lm", fullrange=T, se=F) + 
  geom_point(alpha=.15) +
  facet_wrap("participant") +
  scale_colour_viridis_d(labels=c("Flow lag = 1", "Flow lag = 2", "Flow lag = 3", "Flow lag = 4"),
                         direction = -1) +
  scale_size_discrete(range=c(1, 0.2), labels=c("Flow lag = 1", "Flow lag = 2", "Flow lag = 3", "Flow lag = 4")) +
  theme_bw(base_size=13) + xlab("Lagged flow score") + ylab("Flow score") +
  guides(size=guide_legend(title=NULL), color=guide_legend(title=NULL)) +
  theme(legend.position = c(x=.8, y=.1),
        legend.background = element_rect(fill="white"),
        axis.text.x = element_text(size=8),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="grey"))

  
summary(lmer(flow ~ deviation + lag_deviation1 + lag_deviation2 + lag_deviation3 + lag_deviation4 + (deviation|participant), data=lag_fss_learning))


#probably wrong (have to create lag-variables grouped by participant...):
#summary(lmer(flow ~ lag(flow, 2) + lag(flow, 3) + lag(flow, 4) + lag(flow, 5) + (1|participant), data=fss_learning))


# fss_models = list()
# for (i in 1:length(unique(fss_learning$participant))) {
#   fss_models[[i]] = lm(flow ~ deviation*cumrun, data=subset(fss_learning, ID==i))
# }
# 
# fss_plots = list()
# for (i in 1:length(fss_models)) {
#   fss_plots[[i]] = plot_model(fss_models[[i]], type = "int", terms = c("deviation", "cumrun"), 
#                               mdrt.values = "meansd", axis.title="", title="", axis.labels="") + 
#     xlab(NULL)
# }
# 
# 
# fss_all_plots <- ggarrange(fss_plots[[1]], fss_plots[[2]], fss_plots[[3]], fss_plots[[4]], fss_plots[[5]],
#                            fss_plots[[6]], fss_plots[[7]], fss_plots[[8]], fss_plots[[9]], fss_plots[[10]],
#                            fss_plots[[11]], fss_plots[[12]], fss_plots[[13]], fss_plots[[14]], fss_plots[[15]],
#                            fss_plots[[16]], fss_plots[[17]], fss_plots[[18]], 
#                            font.label = list(size = 10, color = "black", face = "bold"), common.legend = TRUE, legend="right", hjust = -0.25, ncol = 4, nrow = 5)
# 
# fss_all_plots <- annotate_figure(fss_all_plots,
#                                  left = text_grob("Flow scores", color = "black", size=14, rot = 90),
#                                  bottom= text_grob("Deviation score", color ="black", size=14),
#                                  fig.lab = NULL, fig.lab.face = NULL)

#shuffle elements simultaneously:
#df[sample(nrow(df)),]

#Background variables
background19 <- read.csv("/Users/jpalomak/Downloads/background_2019.csv") #does not work out of the box! Fix later.
names(background19)[2] <- "gender"
names(background19)[3] <- "age"
names(background19)[4] <- "license"
names(background19)[6] <- "driving_experience"
names(background19)[7] <- "gaming_experience" #unclear how this is coded in 2017 so not analysed yet!
background19 <- background19 %>% dplyr::select(ID, gender, age, driving_experience, gaming_experience)
background19$gender <- factor(background19$gender, labels=c("m", "f"))
background19$driving_experience <- as.numeric(as.character((factor(background19$driving_experience, labels=c("2", "3", "5", "4", "1")))))
background19$gaming_experience <- as.numeric(as.character((factor(background19$gaming_experience, labels=c("3", "4", "1", "2")))))

background17 <- data.frame(ID = seq(1:9),
                           gender = c("m", "m", "f", "f", "m", "m", "f", "m", "m"),
                           age = c(27,23,23,28,27,22,31,38,25), 
                           driving_experience = c(2,3,1,1,4,4,3,5,3),
                           gaming_experience = c(4,4,3,4,4,3,1,3,4))

background_flowcar <- rbind(background17, background19)
background_flowcar <- background_flowcar %>% mutate(ID = seq(1:18))

#expand each row into 40, except for participant 8 (who has 39 rows)
background_flowcar <- background_flowcar[rep(row.names(background_flowcar), each=40), 1:5]

#remove one row, can be any row, of participant 8
background_flowcar <- background_flowcar[-317,] %>% dplyr::select(-ID)

fss_learning_full <- cbind(fss_learning, background_flowcar)
fss_learning_test <- lmer(flow ~ deviation*cumrun + age + gender + driving_experience + gaming_experience + (1|participant), data=fss_learning_full)

#gaming- and driving experience predict duration, but not after controlling for gender
fss_learning_test2 <- lmer(ln.duration ~ ln.cumrun + age + driving_experience + gaming_experience + (1|participant), data=fss_learning_full)
fss_learning_test3 <- lmer(ln.duration ~ ln.cumrun + age + gender + driving_experience + gaming_experience + (1|participant), data=fss_learning_full)


#Driving experience coding:
# 1 = 0-1000
# 2 = 1000-10000
# 3 = 10000-30000
# 4 = 30000-100000
# 5 = 100000-300000

#Gaming experience coding:
# 1 = Little / none
# 2 = A few times per year
# 3 = 1-3 hours per month
# 4 = At least one hour per week


#age: summary(lmer(ln.duration~age*ln.cumrun + (1|participant), data=fss_learning_full)) (no effect)

#gender summary(lmer(ln.duration~gender*ln.cumrun + (1|participant), data=fss_learning_full))

#driving experience summary(lmer(ln.duration~driving_experience*ln.cumrun + (1|participant), data=fss_learning_full))

apx1 <- ggplot(fss_learning_full, aes(cumrun, duration, colour=gender)) +
  geom_point(alpha=.1, size=2) +
  geom_smooth(method = "nls", se = FALSE, formula = 'y~a*x^b') +
  theme_bw(base_size = 14) +
  theme(legend.position = c(x=.7, y=.7),
        legend.background = element_rect(fill="white"),
        legend.text = element_text(size=14),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="grey")) +
  ylab("Run duration in seconds") + xlab(NULL) +
  labs(title=NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_colour_manual(name=NULL, labels=c("Female", "Male"), values=c("red3", "blue3"))
  # annotate(geom="text", x=1.2, y=6, label="Gender*ln(Cumulative runs) interaction:\nF(1, 16.006) = 18.6, p <.001", size=3.5,
  #        color="black") +
  # annotate(geom="text", x=1.2, y=5.87, label="Gender main effect:\nF(1, 16.005) = 27.35, p <.001", size=3.5,
  #          color="black")


apx2 <- ggplot(fss_learning_full, aes(cumrun, duration, colour=factor(driving_experience))) +
  geom_point(alpha=.1, size=2) +
  geom_smooth(method = "nls", se = FALSE, formula = 'y~a*x^b') +
  theme_bw(base_size = 14) +
  theme(legend.position = c(x=.6, y=.7),
        legend.background = element_rect(fill="white"),
        #legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="grey")) +
  ylab(NULL) + xlab(NULL) +
  labs(title=NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_colour_discrete(name=NULL, labels=c("No lifetime driving experience", "1000-10000 km", "10000-30000 km", "30000-100000 km", "100000-300000 km"))


apx3 <- ggplot(fss_learning_full, aes(cumrun, duration, colour=factor(gaming_experience))) +
  geom_point(alpha=.1, size=2) +
  geom_smooth(method = "nls", se = FALSE, formula = 'y~a*x^b') +
  theme_bw(base_size = 14) +
  theme(legend.position = c(x=.6, y=.7),
        legend.background = element_rect(fill="white"),
        #legend.title = element_text(size=12),
        legend.text = element_text(size=10),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="grey")) +
  ylab(NULL) + xlab(NULL) +
  labs(title=NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_colour_discrete(name=NULL, labels=c("Little/no gaming experience", "A few times per year", "1-3 hours per month", "At least one hour per week"))


##multilevel EFA and CFA

#fa(r = cor(fss[,4:13]), nfactors = 2, fm="ml", rotate="oblimin") ## this is not accounting for multilevel structure!

twolevel <- '
level: 1
  f1 =~ Q2 + Q4 + Q5 + Q7 + Q8 + Q9
  f2 =~ Q1 + Q3 + Q6 + Q10
level: 2
  fb1 =~ Q2 + Q4 + Q5 + Q7 + Q8 + Q9
  fb2 =~ Q1 + Q3 + Q6 + Q10
'

fss2 <- fss %>% mutate(session = rep(rep(1:8, each = 5), 18))
results <- cfa(twolevel, data = fss2[,c(1,4:13)], cluster = 'participant')
summary(results, fit.measures = T, standardized = T)

#measEq.syntax(results, data = fss2, group = "session") #ei toimi, selvittele joskus?


#perceived competence (item: "Osaamiseni taso on... matala/korkea")
#Perceived competence increases over time (sessions):
summary(lmer(comp2 ~ session + (1|participant), data=fss_learning))
ggplot(fss_learning, aes(session, comp2)) + geom_point() + facet_wrap("participant") + geom_smooth(method="lm", se=FALSE) +
  ylab("Perceived competence")

#controlling for duration, session is not significant (which means that perceived competence is linked primarily to success in trials):
summary(lmer(comp2 ~ ln.duration + session + (1|participant), data=fss_learning))
ggplot(fss_learning, aes(ln.duration, comp2)) + geom_point() + facet_wrap("participant") + geom_smooth(method="lm", se=FALSE) +
  ylab("Perceived competence") + coord_cartesian(xlim=c(5.07,5.5))

#Interesting interaction, not sure what to make of it:
fss_learning_comp1 <- lmer(flow ~ comp2*deviation + (deviation|participant), data=fss_learning)
plot_model(fss_learning_comp1, type = "pred", terms = c("comp2", "deviation"), mdrt.values = "meansd")


#Miscellaneous:

fss_learning_misc <- lmer(flow ~ cumrun + (cumrun|participant / session), data=fss_learning)

#partially pooled lines + lms within session:
ggplot(fss_learning, aes(cumrun, flow, group=session)) +
  geom_point(alpha=.3, size=1) +
  geom_line(y=predict(fss_learning_misc), size=0.6, color="blue") +
  geom_smooth(method="lm", size=0.6, color="red", se=FALSE) +
  facet_wrap(~participant) +
  theme_bw(base_size = 14)


#Check is effect of duration_dif on flow is mediated by deviation
mediation <- ' # direct effect
             flow ~ c*duration_dif
           # mediator
             deviation ~ a*duration_dif
             flow ~ b*deviation
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '

#Multilevel version, not sure if correct...
mediation <- ' 
level: 1
deviation ~ a*duration_dif
flow ~ b*deviation
flow ~ duration_dif

level: 2
flow ~ 1

#indirect effect
ab := a*b
'

fit <- sem(mediation, data = lag_fss_learning, cluster="participant") #se = "bootstrap" for bootstrapped CI
summary(fit, ci = T)

labels = list(duration_dif = "Feedback", deviation = "Deviation score", flow = "Flow")
labels2 = c("Flow", "Deviation", "Feedback")

# significant standardized paths only
lavaanPlot(model = fit, labels = labels, node_options = list(shape = "box", fontname = "Helvetica"), 
           edge_options = list(color = "grey"), coefs = TRUE, stars = "regress", stand = TRUE)

#semPaths(fit,rotation=1,sizeMan=15,style="lisrel",
 #        curvePivot=TRUE, what="par", nodeLabels=labels2)

#Note semPaths doesn't seem to easily handle multilevel SEMs.. google!
semPaths(fit, "est", "std", nodeLabels=labels2, sizeMan=15, 
         edge.label.cex=1.5, style="lisrel", intercepts=FALSE, residuals=FALSE, fade=TRUE)


#come back to this later!!!:::
# Multilevel models assume that the observed variables have at least two potential levels of variation. 
# Because temporal lag was experimentally manipulated within subjects, it does not vary between subjects. 
# On the other hand, hit rates (HR) vary both between and within participants: At the lower (within-person) level, 
# HR varies from trial to trial. At the upper level, we may also expect that HR varies, on average, 
# between participants. We are most interested in the within-person process, and therefore it is useful to 
# transform the variables such that these two levels are explicitly separated from each other (Bolger & Laurenceau, 2013). 
# Notice that this transformation is not strictly required, but this reasoning suggests that it is often useful and meaningful in 
# data sets where the predictor values vary both between and within subjects. 
# We first averaged the grand-mean-centered trial-level HR for each person to create a between-person component of HR. 
# We subtracted these means from the raw HR to create within-subject trial-by-trial deviations from the subject-means that 
# represent an entirely within-person version of HR. Isolating the within-person process from variables can 
# be done by using bmlm’s isolate() function:

#using mlm() from https://link.springer.com/article/10.3758/s13428-017-0980-9#Tab1
#https://mvuorre.github.io/bmlm/articles/bmlm-blch9/bmlm-blch9.html
#need first to remove "NA" value from duration_dif (first value)
lag_fss_learning2 <- lag_fss_learning %>% dplyr::select(duration_dif, deviation, flow)
lag_fss_learning2$duration_dif[1] <- 0
lag_fss_learning2 <- isolate(lag_fss_learning2, by = "participant",
                             value = c("duration_dif", "deviation", "flow")) #separate between subjects variability
lag_fss_learning2 <- as.data.frame(lag_fss_learning2) #change into dataframe -- cannot be tibble for mlm()!

fit <- mlm(d = lag_fss_learning2, 
           id = "participant",
           x = "duration_dif_cw",
           m = "deviation_cw",
           y = "flow_cw",
           iter = 1000, 
           cores = 1)

mlm_summary(fit)

mlm_path_plot(fit, level = .95, text = T,
              xlab = "Feedback",
              mlab = "Deviation score",
              ylab = "Self-reported flow", digits = 2)

mlm_pars_plot(fit, pars = "u_me", type = "coef", level = .90) #mediated effect for all participants, 80% credible intervals


#simple LM fit per session, and lm and loess across sessions
# ggplot(fss_learning, aes(cumrun, flow, group=session)) +
#   geom_point(alpha=.3, size=1) +
#   geom_smooth(method="lm", size=0.5, color="blue", se=FALSE) +
#   geom_smooth(method="lm", size=0.5, color="red", se=FALSE, inherit.aes=FALSE, aes(x=cumrun, y=flow)) +
#   geom_smooth(method="loess", size=0.8, color="black", se=FALSE, linetype=2, inherit.aes=FALSE, aes(x=cumrun, y=flow)) +
#   facet_wrap(~participant) +
#   theme_bw(base_size = 14)


#fluency and absorption mean, median, range
fss_learning %>% group_by(participant) %>% dplyr::summarize(fluency_m=mean(fluency), fluency_med=median(fluency), fluency_range=paste(as.character(min(fluency)), as.character(max(fluency)), sep=" - "), 
                                                     absorption_m=mean(absorption), absorption_med=median(absorption), absorption_range=paste(as.character(min(absorption)), as.character(max(absorption)), sep=" - "))

#summary stats, age, gaming and driving experience:
fss_learning_full %>% group_by(participant) %>% dplyr::summarize(age = mean(age), gender = mean(as.numeric(gender)), driving_experience = mean(driving_experience), gaming_experience = mean(gaming_experience))

# Basic violin plot, flow
flow1 <- fss_items %>%
  mutate(participant = fct_reorder(participant, desc(flow), .fun='median')) %>%
  ggplot( aes(x=participant, y=flow)) + 
  geom_violin() + geom_hline(aes(yintercept=median(flow)), color="firebrick") +
  geom_jitter(color="black", size=0.6, alpha=0.2, width=0.07) +
  # geom_half_violin(trim=FALSE, fill="gray", side = "r") +
  # geom_boxplot(width=0.1, outlier.shape = NA) +
  geom_hline(yintercept = median(fss_items$flow), color = "firebrick") +
  #geom_text(aes(length(unique(fss_items$participant)), median(fss_items$flow), label = "med.\nFlow", vjust = -1)) +
  labs(title="A", x="Participant", y = "Flow") +
  ylim(2, 7) +
  stat_summary(fun = median, geom = "point",
               size = 2, color="firebrick") +
  theme_classic(base_size=14)

#Session-wise flow scores, new data
flow2 <- fss_learning_full %>% 
  mutate(session = fct_reorder(as.factor(session), desc(flow), .fun='median')) %>% 
  ggplot(aes(as.factor(session), flow)) + 
  geom_violin() + geom_hline(aes(yintercept=median(flow)), color="firebrick") +
  geom_jitter(color="black", size=0.6, alpha=0.2, width=0.07) +
  ylab(NULL) + xlab("Session") +
  ylim(2, 7) +
  labs(title="B") +
  stat_summary(fun = median, geom = "point",
               size = 2.5, color="firebrick") +
  theme_classic(base_size=14)


#Combine flow1 and flow2 figures
ggarrange(flow1, flow2, font.label = list(size = 10, color = "black", face = "bold"), common.legend = FALSE, hjust = -0.25, ncol = 2, nrow = 1)


# Basic violin plot, perceived importance items
apx4 <- fss_items %>%
  mutate(participant = fct_reorder(participant, desc(pi_total), .fun='median')) %>%
  ggplot( aes(x=participant, y=pi_total)) + 
  geom_violin() + geom_hline(aes(yintercept=median(pi_total)), color="firebrick") +
  geom_jitter(color="black", size=0.6, alpha=0.2, width=0.07) +
  # geom_half_violin(trim=FALSE, fill="gray", side = "r") +
  # geom_boxplot(width=0.1, outlier.shape = NA) +
  geom_hline(yintercept = median(fss_items$pi_total), color = "firebrick") +
  #geom_text(aes(length(unique(fss_items$participant)), median(fss_items$pi_total), label = "med.\nFlow", vjust = -1)) +
  labs(title="Perceived importance", x=NULL, y = NULL) +
  ylim(2, 7) +
  stat_summary(fun = median, geom = "point",
               size = 2, color="firebrick") +
  theme_classic(base_size=12)

apx5 <- fss_items %>%
  mutate(participant = fct_reorder(participant, desc(fluency), .fun='median')) %>%
  ggplot( aes(x=participant, y=fluency)) + 
  geom_violin() + geom_hline(aes(yintercept=median(fluency)), color="firebrick") +
  geom_jitter(color="black", size=0.6, alpha=0.2, width=0.07) +
  # geom_half_violin(trim=FALSE, fill="gray", side = "r") +
  # geom_boxplot(width=0.1, outlier.shape = NA) +
  geom_hline(yintercept = median(fss_items$fluency), color = "firebrick") +
  #geom_text(aes(length(unique(fss_items$participant)), median(fss_items$fluency), label = "med.\nFlow", vjust = -1)) +
  labs(title="Fluency of performance", x=NULL, y = NULL) +
  ylim(2, 7) +
  stat_summary(fun = median, geom = "point",
               size = 2, color="firebrick") +
  theme_classic(base_size=12)

apx6 <- fss_items %>%
  mutate(participant = fct_reorder(participant, desc(absorption), .fun='median')) %>%
  ggplot( aes(x=participant, y=absorption)) + 
  geom_violin() + geom_hline(aes(yintercept=median(absorption)), color="firebrick") +
  geom_jitter(color="black", size=0.6, alpha=0.2, width=0.07) +
  # geom_half_violin(trim=FALSE, fill="gray", side = "r") +
  # geom_boxplot(width=0.1, outlier.shape = NA) +
  geom_hline(yintercept = median(fss_items$absorption), color = "firebrick") +
  #geom_text(aes(length(unique(fss_items$participant)), median(fss_items$absorption), label = "med.\nFlow", vjust = -1)) +
  labs(title="Absorption by activity", x=NULL, y = NULL) +
  ylim(2, 7) +
  stat_summary(fun = median, geom = "point",
               size = 2, color="firebrick") +
  theme_classic(base_size=12)


#combine apx1-3 (gender, driving- and gaming experience)
background_apx <- ggarrange(apx1, apx2, apx3, font.label = list(size = 10, color = "black", face = "bold"), 
          common.legend = FALSE, hjust = -0.25, ncol = 3, nrow = 1)

background_apx <- annotate_figure(background_apx, 
                bottom = text_grob("Cumulative runs", color = "black", size=14),
                fig.lab = NULL, fig.lab.face = NULL)

#combine apx4-6 (PI items, fluency and absorption)
scales_apx <- ggarrange(apx4, apx5, apx6, font.label = list(size = 10, color = "black", face = "bold"), 
          common.legend = FALSE, hjust = -0.25, ncol = 3, nrow = 1)

scales_apx <- annotate_figure(scales_apx, 
                                  bottom = text_grob("Participant", color = "black", size=14),
                                  fig.lab = NULL, fig.lab.face = NULL)


#scatterplot matrix of game variables
my_fn <- function(data, mapping, method="lm", ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(alpha=.2, size=0.5) + 
    geom_smooth(method=method, size=0.5)
  p
}
game_data %>%
  dplyr::select(cumrun, min_velocity, max_velocity, end_velocity, avg_velocity, collisions, speed_drops, duration) %>%
  ggpairs(columnLabels = c("Cumulative\nruns", "Min\nvelocity", "Max\nvelocity", "End\nvelocity", "Average\nvelocity", "Collisions", "Speed\ndrops", "Run\nduration"),
          lower = list(continuous = my_fn)) +
#          upper = list(continuous = wrap("cor", size = 3.5))) + 
  ggtitle("Scatterplot matrix of game variables")


#Or using psych-library:
game_data %>%
  dplyr::select(cumrun, min_velocity, max_velocity, end_velocity, avg_velocity, collisions, speed_drops, duration) %>% 
  pairs.panels(method = "spearman", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE,
             smoother = TRUE) # show correlation ellipses


#Biplot of principle components of data
PCA_biplotdata <- cbind(game_data, fss_learning_full)
PCA_biplotdata <- PCA_biplotdata[,!duplicated(names(PCA_biplotdata))] #remove duplicate columns

components <- PCA_biplotdata %>%
  dplyr::select(session, duration, run, min_velocity, max_velocity, end_velocity, avg_velocity, collisions, speed_drops,
               ln.duration, cumrun, collisions, absorption, fluency, flow, pi_total, age, driving_experience, gaming_experience) %>% prcomp()

#biplot(components) FIX THIS, ugly
ggbiplot(components, alpha=.1) + theme_minimal() #+ coord_cartesian(xlim=c(-2.5,1.5), ylim=c(-1, 1))



#Plot base R and ggplot in the same grid (sliding windows -analyses):
within_windows <- sliding_window_within(10,1)

par(mfrow = c(1,2), mar=c(4.2,6,2,0), oma=c(0,0,0,4))

# leave top-left quadrant empty!
plot.new()

# plot regular R graphic in top-right quadrant
plot(c(idx[1],idx[length(idx)]), c(rng[1], rng[2]+3) ,type="n", bty="n", xlab="", ylab="", main = NULL)
plotlines(curve, idx, col ="blue")

# now add the first ggplot element which is going to take
# up the left two quadrants
vp <- viewport(height = unit(1,"npc"), width=unit(0.5, "npc"), 
               just = c("left","top"), y = 1, x = 0)
               
print(within_windows, vp = vp)

#CUSTOM Cronbach's alpha for each trial. 1= full scale, 2 = fluency, 3 = absorption
alpha_calculator <- function(subscale=1) {
  cronbachs_alphas <- cbind(fss2[-317,], fss_learning_full) #remove one row from participant 8
  cronbachs_alphas$run_session <- as.integer(interaction(fss_learning$run, fss_learning$session)) #expand run integer into 1:40
  alphas_list <- list()
  alphas <- list()
  counter = 0
  if (subscale == 1) {
    for (i in 1:length(unique(cronbachs_alphas$run_session))) { #calculate cronbach's alpha separately for each run
      counter = counter + 1
      alphas_list[[counter]] <- subset(cronbachs_alphas, run_session==i) %>% dplyr::select(Q2:Q10) %>% psych::alpha()
      alphas[[counter]] <- summary(alphas_list[[counter]])$std.alpha
    }
    alphas <- as.data.frame(unlist(alphas))
    names(alphas)[1] <- "alpha"
    
    return(ggplot(alphas, aes(x=alpha)) + geom_histogram(alpha=.7, color="grey", fill="firebrick") +
             scale_x_continuous(name = NULL, breaks=c(0.55, 0.65, 0.75, 0.85, 0.95), limits=c(0.50, 1)) +
             ylab("Count") + labs(title="Full FSS scale") + scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9), limits=c(0,9)) +
             theme_classic(base_size=12))
  }
  else if (subscale == 2) {
    for (i in 1:length(unique(cronbachs_alphas$run_session))) { #calculate cronbach's alpha separately for each run
      counter = counter + 1
      alphas_list[[counter]] <- subset(cronbachs_alphas, run_session==i) %>% dplyr::select(Q2, Q7, Q8, Q9) %>% psych::alpha()
      alphas[[counter]] <- summary(alphas_list[[counter]])$std.alpha
    }
    alphas <- as.data.frame(unlist(alphas))
    names(alphas)[1] <- "alpha"
    
    return(ggplot(alphas, aes(x=alpha)) + geom_histogram(alpha=.7, color="grey", fill="firebrick") +
             scale_x_continuous(name = NULL, breaks=c(0.55, 0.65, 0.75, 0.85, 0.95), limits=c(0.50, 1)) +
             ylab(NULL) + labs(title="Fluency of performance") + scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9), limits=c(0,9)) +
             theme_classic(base_size=12))
  }
  else if (subscale == 3) {
    for (i in 1:length(unique(cronbachs_alphas$run_session))) { #calculate cronbach's alpha separately for each run
      counter = counter + 1
      alphas_list[[counter]] <- subset(cronbachs_alphas, run_session==i) %>% dplyr::select(Q3, Q4, Q5, Q6, Q10) %>% psych::alpha() 
      alphas[[counter]] <- summary(alphas_list[[counter]])$std.alpha
    }
    alphas <- as.data.frame(unlist(alphas))
    names(alphas)[1] <- "alpha"
    
    return(ggplot(alphas, aes(x=alpha)) + geom_histogram(alpha=.7, color="grey", fill="firebrick") +
             scale_x_continuous(name = NULL, breaks=c(0.55, 0.65, 0.75, 0.85, 0.95), limits=c(0.50, 1)) +
             ylab(NULL) + labs(title="Absorption by activity") + scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9), limits=c(0,9)) +
             theme_classic(base_size=12))
  }
}

#all in same figure
alpha_reliability <- ggarrange(alpha_calculator(1), alpha_calculator(2), alpha_calculator(3), 
          font.label = list(size = 10, color = "black", face = "bold"), 
          common.legend = TRUE, legend="bottom", hjust = -0.25, ncol = 3, nrow = 1)

alpha_reliability <- annotate_figure(alpha_reliability, 
                                     bottom = text_grob("Cronbach's alpha", color = "black", size=14),
                                     fig.lab = NULL, fig.lab.face = NULL)

#multilevel.reliability
mlr(cronbachs_alphas, grp = "ID", Time = "run_session", items = c("Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10"),alpha=TRUE,icc=FALSE, aov=TRUE,
    lmer=FALSE,lme = TRUE, long=FALSE,values=NA,na.action="na.omit",plot=TRUE,
    main="Lattice Plot by subjects over time")


#effect of flow on learning (slope of model ln.duration~ln.cumrun); this is a bad analysis
# uusfreimi <- as.data.frame(cbind(ranef(fss_learning_lmer_subset2p)[["participant"]][2], fss_learning %>% group_by(participant) %>% dplyr::summarize(mean(flow))))
# names(uusfreimi)[3] <- "flow"
# names(uusfreimi)[1] <- "learning"
# ggplot(uusfreimi, aes(flow, learning)) + geom_point() + geom_smooth(method="lm")

#################################
#################################
#################################
###Below here stuff for safekeeping that I will probably use later.

#+theme(panel.background=element_rect(fill='white', colour='black'),
#      strip.background=element_rect(fill='white', colour='white'))


# flow_all_figs <- ggarrange(one, two, three, font.label = list(size = 10, color = "black", face = "bold"), common.legend = TRUE, legend="bottom", hjust = -0.25, ncol = 3, nrow = 1)
# 
# flow_all_figs <- annotate_figure(flow_all_figs,
#                                 bottom = text_grob("Deviation score", color = "black", size=14),
#                                 fig.lab = NULL, fig.lab.face = NULL)


#Flow subscales and perceived importance

#PI scale:
#1: Koin pelissä onnistumisen tärkeäksi
#2: Minusta tuntui siltä, etten saisi tehdä yhtäkään virhettä
#3: Pelkäsin epäonnistuvani

#2 and 3 are negatively framed, and #3 has an inverse relation with deviation! Thus, probably not wise to combine as singular scale

# psych::alpha(fss_learning[,16:18])
# fss_learning_lmer <- lmer(pi2 ~ deviation + (deviation|participant), data=fss_learning)
# summary(fss_learning_lmer)
# 
# ggplot(fss_learning, aes(deviation, pi2)) +
#   geom_point(alpha=.6, size=2) +
#   geom_smooth(method = "lm", se = FALSE, fullrange=TRUE) + 
#   facet_wrap(~participant) +
#   theme_bw(base_size = 14)


# jet.colors <- colorRampPalette(matlab.like(9))
# ggplot(gg, aes(x=age, y=hrs, z=charges))+
#   stat_contour(aes(color=..level..),binwidth=200, size=2)+
#   scale_color_gradientn(colours=jet.colors(8))

#maybe useful at some point..
# all(is.na(temp_data$flow))
# all(is.na(temp_data$deviation))
# all(is.na(temp_data$cumrun))


#However, as we discussed earlier, the impact of positive feedback may decrease over time, 
#as a participant becomes more familiar with the task, and develops an increasingly stable sense 
#of how good or bad he/she is at the task.


###
#' Manuel Koller (2016). robustlmm: An R Package for Robust Estimation of Linear Mixed-Effects Models. Journal of Statistical Software, 75(6), 1-24.<doi:10.18637/jss.v075.i06>
#'   
#'   Corresponding BibTeX entry:
#'   
#'   @Article{,
#'     title = {{robustlmm}: An {R} Package for Robust Estimation of
#'       Linear Mixed-Effects Models},
#'     author = {Manuel Koller},
#'     journal = {Journal of Statistical Software},
#'     year = {2016},
#'     volume = {75},
#'     number = {6},
#'     pages = {1--24},
#'     doi = {10.18637/jss.v075.i06},
#'   }
#'   
#'   
