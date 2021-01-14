# ---
# title: "Data analysis pipeline for the paper "The Link Between Flow and Performance is Moderated by Task Experience"
# authors: "Jussi Palomäki, Tuisku Tammi, Noora Lehtonen, Niina Peltonen, Sami Abuhamdeh, Otto Lappi, Benjamin Ultan Cowley"
# note: This script depends on data wrangled in combine_data.R; analyses done by JP, TT and BC
# ---

# Load libraries
library(easypackages)
libraries("ggplot2", "gridBase", "colorRamps", "grDevices", "lme4", "effects",
          "multcomp", "tidyverse", "lmerTest", "ggpubr", "MuMIn", "lattice",
          "robustlmm", "grid", "sjPlot", "dynlm", "viridis", "psych", "lavaan",
          "GPArotation", "lavaanPlot", "jtools", "interactions", "XML",
          "semPlot", "GGally", "ggbiplot", "gghalves", "bmlm") # see https://link.springer.com/article/10.3758/s13428-017-0980-9#Sec10


#Remove first row (cumrun == 1) for each participant. NOT USED IN THE INTERACTION ANALYSES
fss_learning_sub <- subset(fss_learning, cumrun!=1)
anova(lmer(flow ~ deviation + (deviation|participant), data=fss_learning_sub))
anova(lmer(flow ~ deviation + (deviation|participant), data=fss_learning))


############
#HYPOTHESIS 1A, ANALYSIS: Replication of results of Cowley et al. (2019), effect of learning, power-law curve
############

H1a_lmer <- lmer(ln.duration ~ ln.cumrun + (ln.cumrun|participant), data=subset(fss_learning, as.numeric(ID) > 9)) #NOTE that this way of subsetting no longer works correctly when the ID factor levels are reordered (as will be done later in the pipeline!)
summary(H1a_lmer) # model summary
plot(H1a_lmer) # model diagnosticsi
qqnorm(residuals(H1a_lmer)) #qq-plot
qqline(residuals(H1a_lmer)) #line of "perfect normality"


############
#HYPOTHESIS 1B, ANALYSIS: The f~d effect
############

H1b_lmer <- lmer(flow ~ deviation + (deviation|participant), 
                        data=filter(fss_learning, as.numeric(ID) > 9))
plot(H1b_lmer) # model diagnostics
qqnorm(residuals(H1b_lmer)) #qq-plot
qqline(residuals(H1b_lmer)) #line of "perfect normality"


############
#HYPOTHESIS 1, PLOTTING
############

#Plot log-log with flow z-scores coloring
#First sort participants according to f~d
fss_learning_fdsort <- fss_learning
fss_learning_fdsort$ID <- factor(fss_learning_fdsort$ID, levels=c(1:18))
estimates_slope <- list()
counter = 0
for (x in 1:length(unique(fss_learning_fdsort$ID))) { #iterate through each ID
  estimates_temp <- subset(fss_learning_fdsort, ID==x) #subset one ID at a time
  temp_model <- lm(flow ~ deviation, data=estimates_temp)  #fit model for subsetted ID
  counter = counter + 1 #counter iteration for indexing
  estimates_slope[[counter]] <- coef(summary(temp_model))[2,1] #add values to estimates_slope list
}
estimates_slope <- as.data.frame(rbind(unlist(estimates_slope))) #make list into dataframe
colnames(estimates_slope) <- c(1:18) #rename columns (from V1...V18 to 1...18)
estimates_slope <- sort(estimates_slope, decreasing=F) #sort according to values
part_levels <- as.integer(names(estimates_slope)) #turn into integer values
fss_learning_fdsort$ID <- factor(fss_learning_fdsort$ID, levels=part_levels) #reorder ID factor levels in data according to sorted integers

# plot hypothesis 1a
rq1a <- ggplot(fss_learning_fdsort, aes(ln.cumrun, ln.duration, color = z.flow)) +
  geom_point(alpha=.6, size=1) +
  geom_smooth(method = "lm", se = FALSE, linetype = 1, size = 0.5, color="red") +
  facet_wrap(~ID) +
  scale_color_viridis(name="Z_Flow", guide = guide_colorbar(title.position = "top")) +
  xlab("ln(Cumulative runs)") + ylab("ln(Duration)") +
  labs(title = "A") +
  theme(legend.position = c(x=.8, y=.1),
        legend.background = element_rect(fill="white"),
        legend.box.background = element_rect(colour = "black"),
        legend.direction='horizontal',
        panel.background = element_rect(fill = "white",
                                        colour = "lightgrey",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey"),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="lightgrey", color="black"))
#if (FIGOUT)  ggsave(file.path(odir, "PerfXsubj_powlxFlow_loglog.svg"))

# plot hypothesis 1b
rq1b <- ggplot(fss_learning_fdsort, aes(deviation, flow)) +
  geom_point(alpha=.4, size=1) +
  geom_vline(xintercept=0, linetype=2, alpha=.5) + #line to separate better-than-expected and worse-than-expected performance
  geom_rect(fill="green",alpha=0.005, xmin=min(fss_learning_fdsort$deviation), xmax=0, ymin=0, ymax=max(fss_learning_fdsort$flow)) +
  geom_rect(fill="red",alpha=0.005, xmin=0, xmax=max(fss_learning_fdsort$deviation), ymin=0, ymax=max(fss_learning_fdsort$flow)) +
  geom_smooth(method = "lm", se = FALSE) + 
  xlab("Deviation score") + ylab("Flow score") +
  facet_wrap(~ID) +
  labs(title = "B") +
  theme(legend.position = c(x=.8, y=.1),
        legend.background = element_rect(fill="white"),
        legend.box.background = element_rect(colour = "black"),
        legend.direction='horizontal',
        panel.background = element_rect(fill = "white",
                                        colour = "lightgrey",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey"),
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="lightgrey", color="black"))
#if (FIGOUT)  ggsave(file.path(odir, "FlowXdevXsubj.svg"))

#Combine Rq1a and b into panel figure
RQ1 <- ggarrange(rq1a, rq1b, font.label = list(size = 10, color = "black", face = "bold"), common.legend = FALSE, hjust = -0.25, ncol = 2, nrow = 1)
#ggsave("figure4.pdf", width=12, height=6) #Uncomment if you want the figure saved as pdf.


############
#HYPOTHESIS 2, ANALYSES: Does the effect of Deviation on Flow (main finding from previous paper) depend on level cumrun (how many runs the participant has played)
############

#Fit model (note: could also try random slope/intercept for deviation*cumrun):
H2_lmer <- lmer(flow ~ deviation*cumrun + (deviation*cumrun|participant), data=fss_learning)
summ(H2_lmer)

#Random intercept only
H2_lmer_simple <- lmer(flow ~ deviation*cumrun + (1|participant), data=fss_learning)
summ(H2_lmer_simple)

#Simple slopes:
sim_slopes(H2_lmer_simple, pred = deviation, modx = cumrun, johnson_neyman = TRUE)

#Test model assumptions:
plot(H2_lmer_simple)
plot(resid(H2_lmer_simple), fss_learning$deviation)
plot(resid(H2_lmer_simple), fss_learning$cumrun)
qqnorm(residuals(H2_lmer_simple))
qqline(residuals(H2_lmer_simple))
qqmath(ranef(H2_lmer_simple, condVar = TRUE), strip = TRUE)

#Fit model using robust lmm:
H2_rlmer_simple <- rlmer(flow ~ deviation*cumrun + (1|participant), data=fss_learning)
summary(H2_rlmer_simple)
plot(H2_rlmer_simple)

#Plot interaction from model (figure not shown in paper!):
plot_model(H2_lmer_simple, type = "pred", terms = c("deviation", "cumrun"), mdrt.values = "meansd")


############
#HYPOTHESIS 2, PLOTTING
############

#Participant-wise visualizations of deviation*cumrun interaction effect
#Calculate slopes and intercepts of deviation at various levels of moderator (cumrun) for each participant.
cumrun_minusSD <- mean(fss_learning$cumrun) - sd(fss_learning$cumrun)
cumrun_meanSD <- mean(fss_learning$cumrun)
cumrun_plusSD <- mean(fss_learning$cumrun) + sd(fss_learning$cumrun)
slope_intercept_dataframes = list()
fss_models <- list()
for (i in 1:length(unique(fss_learning$participant))) {
  temp_model = lm(flow ~ deviation*cumrun, data=subset(fss_learning, ID==i)) #Here we're essentially just fitting LMs to each participant separately.
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
fss_learning$ID <- factor(fss_learning$ID, levels=part_levels) #NOTE 10.1. this causes issue when subsetting using ">" on ID

#Plot
flow_fig6x3 <- ggplot(fss_learning, aes(deviation, flow)) +
  geom_point(alpha=.2, size=1.5) +
  #geom_smooth(method = "lm", se = FALSE, fullrange=TRUE) + 
  facet_wrap(~ID, nrow=3, ncol=6) + 
  geom_abline(data = slope_intercept_dataframes, aes(intercept=intercept_minusSD, slope=slope_minusSD, linetype="-1 SD\n(8.81 runs)"), size=0.8) + #-1SD is at 8.81 runs
  geom_abline(data = slope_intercept_dataframes, aes(intercept=intercept_plusSD, slope=slope_plusSD, linetype="+1 SD\n(32.19 runs)"), size=0.8) +  #+1SD is at 32.19 runs
  geom_vline(xintercept=0, linetype=2, alpha=.5) + #line to separate better-than-expected and worse-than-expected performance
  geom_rect(fill="green",alpha=0.005, xmin=min(fss_learning$deviation), xmax=0, ymin=0, ymax=max(fss_learning$flow)) +
  geom_rect(fill="red",alpha=0.005, xmin=0, xmax=max(fss_learning$deviation), ymin=0, ymax=max(fss_learning$flow)) +
  labs(colour = NULL, linetype = "Cumulative runs", title = "A") + xlab("Deviation score") + ylab("Flow score") +
  #scale_x_discrete(breaks = c()) +
  theme_bw(base_size=13) + 
  theme(legend.justification = "top",
        legend.background = element_rect(fill="white"),
        legend.box.background = element_rect(colour = "black"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9),
        axis.text.x = element_text(size=9),
        panel.background = element_rect(fill = "white",
                                        colour = "lightgrey",
                                        size = 0.5, linetype = "solid"),
        panel.border = element_blank(), 
        panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey"),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="lightgray"))


#Participant-wise Sliding-window regression, custom function
sliding_window_within <- function(window_width = 10, result=2) {
  estimates_IDwise <- list()
  counter = 0
  for (x in 1:length(unique(fss_learning$participant))) { #iterate through all participants
    for (i in 1:(length(unique(fss_learning$cumrun))+1-window_width)) {
      fss_learning_temp <- subset(fss_learning, ID==x & cumrun >= i & cumrun < i+window_width) #subset based on selected window width
      
      ###calculate deviation score again for each segment
      fss_learning_temp <- fss_learning_temp %>% group_by(participant) %>%
        mutate(learning_curve = predict(lm(ln.duration ~ ln.cumrun)), deviation = ln.duration - learning_curve)
      ###
      temp_model <- lm(flow ~ deviation, data=fss_learning_temp)
      counter = counter + 1
      estimates_IDwise[[counter]] <- coef(summary(temp_model))[2,1]
    }
  }
  
  estimates_IDwise <- as.data.frame(unlist(estimates_IDwise)) 
  estimates_IDwise <- estimates_IDwise %>% mutate(sliding_window=rep(c(1:(length(unique(fss_learning$cumrun))+1-window_width)), 18),
                                                  ID=rep(1:18, each=length(unique(fss_learning$cumrun))+1-window_width)) #18 is number of participants
 
  names(estimates_IDwise)[1] <- "estimate"
  permuted_estimates <- estimates_IDwise %>% group_by(ID) %>% 
    transform(estimate = sample(estimate))
  
  #code for sorting facets based on slope BEGINS HERE
  #First, fit linear model separately for each participant with sliding window as predictor and slope estimate as DV
  #then grab the slope estimate of that model, and sort according to the values (this is essentially the same as above for p-values)
  estimates_slope <- list()
  counter = 0
  for (x in 1:length(unique(estimates_IDwise$ID))) {
    estimates_temp <- subset(estimates_IDwise, ID==x)
    temp_model <- lm(estimate ~ sliding_window, data=estimates_temp) 
    counter = counter + 1
    estimates_slope[[counter]] <- coef(summary(temp_model))[2,1]

  }
  
  estimates_slope <- as.data.frame(rbind(unlist(estimates_slope)))
  colnames(estimates_slope) <- c(1:18)
  estimates_slope <- sort(estimates_slope)
  part_levels <- as.integer(names(estimates_slope))
  estimates_IDwise$ID <- factor(estimates_IDwise$ID, levels=part_levels)
  #code for sorting facets based on slope ENDS HERE
  
  #choose which result to return (1=Participant-wise figure, 2=DATA VECTOR: permuted estimates, 3=DATA VECTOR: observed estimates):
  #participant-wise
  if (result == 1) {
    result <- ggplot(estimates_IDwise, aes(x=sliding_window, y=estimate)) + geom_line() + 
      ylab("Slope estimate") + xlab(NULL) + #xlab should be "sliding window" but changed now for cogcarsim #2 paper
      theme_bw(base_size=13) + scale_x_continuous(breaks = c(1, 10, 20, 31)) +
      theme(axis.text.x = element_text(size=8)) +
      facet_wrap("ID", ncol=6, nrow=3) + geom_smooth(method="lm", se=FALSE, size=0.3) + labs(title = "B") #note that title and ncol, nrow for facet_wrap are now specific to cogcarsim #2 paper!
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

#Combine the simple slopes figure (flow_fig6x3) with sliding windows regression figure.
#Use Minimum Width Envelope of the actual observations (sliding windows of slopes within participants)
mwe_fig <- ggplot(cbind(curve, idx), aes(x = idx, y = mean0)) +
  geom_ribbon(aes(ymin = lo, ymax = up, fill = "red"), alpha = 0.2, linetype = 0) +
  geom_ribbon(aes(ymin = lo0, ymax = up0, fill = "blue"), alpha = 0.2, linetype = 0) +
  ylim(-20,0) + scale_x_continuous(breaks=c(1, 10, 20, 31)) +
  geom_line(size = 0.75) +
  labs(title = "", x = NULL, y = NULL) + #title = "" to force empty space near top
  theme_bw() +
  theme(legend.position = "none")

combined_plot1 <- ggarrange(sliding_window_within(10,1), mwe_fig, ncol=2, nrow=1, widths = c(2,1.5))
combined_plot1 <- annotate_figure(combined_plot1, 
                                  bottom = text_grob("Sliding window", color = "black", size=14),
                                  fig.lab = NULL, fig.lab.face = NULL)

ggarrange(flow_fig6x3, combined_plot1, ncol=1, nrow=2)

#ggsave("combined_interaction_sliding.pdf", width=10, height=12) #Uncomment to save figure as pdf


############
#Plot basic stats
############

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
flow2 <- fss_learning %>% 
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


############
#Visualize effect of deviation vs. skill
############

testilmer <- lmer(flow ~ ln.duration + deviation + (1|ID), data=fss_learning)

nd <- with(fss_learning,
           rbind(data.frame(expand.grid(ID=levels(ID),
                                        ln.duration=seq(min(ln.duration),max(ln.duration),length=51)),
                            deviation=mean(deviation),focus="ln.duration"), #to predict for values of duration while holding deviation constant at its mean
                 data.frame(expand.grid(ID=levels(ID),
                                        deviation=seq(min(deviation),max(deviation),length=51)),
                            ln.duration=mean(ln.duration),focus="deviation"))) #to predict for values of deviation while holding duration constant at its mean

nd$val <- with(nd, c(ln.duration[focus=="ln.duration"], deviation[focus=="deviation"]))
nd$ID <- factor(nd$ID, levels=1:18)
nd$focus <- factor(nd$focus, labels = c("log(Duration)", "Deviation score"))

pframe <- data.frame(nd,resp=predict(testilmer,newdata=nd))

ggplot(pframe,aes(x=val,y=resp,colour=ID)) + geom_line()+
  facet_wrap(~focus,scale="free_x") + xlab("Value") + ylab("Flow score") + theme_bw(base_size=12) +
  ylim(2.5, 6.5) + scale_colour_discrete(name = "Participant") + theme(legend.text = element_text(size=9),
                                                                       legend.title = element_text(size=12))

#ggsave("partial_effects_lmer.pdf", width=6, height=5) #uncomment to save figure as pdf


############
############
############

#APPENDIX MATERIALS BEGIN HERE#

############
############
############

############
#Background variables
############

#background19 <- read.csv("/Users/jpalomak/Downloads/background_2019.csv") #already done in pipeline
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

#remove one row, can be any row, of participant 8 (row number 306; session 6, run 1, is missing in fss_learning; here doesn't matter which row since background is measured once per participant)
background_flowcar <- background_flowcar[-306,] %>% dplyr::select(-ID)
fss_learning_full <- cbind(fss_learning, background_flowcar)

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

#combine apx1-3 (gender, driving- and gaming experience)
background_apx <- ggarrange(apx1, apx2, apx3, font.label = list(size = 10, color = "black", face = "bold"), 
                            common.legend = FALSE, hjust = -0.25, ncol = 3, nrow = 1)

background_apx <- annotate_figure(background_apx, 
                                  bottom = text_grob("Cumulative runs", color = "black", size=14),
                                  fig.lab = NULL, fig.lab.face = NULL)

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

#combine apx4-6 (PI items, fluency and absorption)
scales_apx <- ggarrange(apx4, apx5, apx6, font.label = list(size = 10, color = "black", face = "bold"), 
                        common.legend = FALSE, hjust = -0.25, ncol = 3, nrow = 1)

scales_apx <- annotate_figure(scales_apx, 
                              bottom = text_grob("Participant", color = "black", size=14),
                              fig.lab = NULL, fig.lab.face = NULL)


############
#Lagged flow regressions
############

#NOTE: dplyr::mutate works, but plyr:mutate() doesn't work correctly with group_by and lag()
lag_fss_learning <- fss_learning %>% group_by(participant) %>% dplyr::mutate(lag_flow1 = lag(flow, 1), lag_flow2 = lag(flow, 2), lag_flow3 = lag(flow, 3), 
                                                                      lag_flow4 = lag(flow, 4), lag_flow5 = lag(flow, 5),
                                                                      lag_deviation1 = lag(deviation, 1), lag_deviation2 = lag(deviation, 2), lag_deviation3 = lag(deviation, 3),
                                                                      lag_deviation4 = lag(deviation, 4), lag_deviation5 = lag(deviation, 5),
                                                                      lag_duration1 = lag(ln.duration, 1), lag_duration2 = lag(ln.duration, 2),
                                                                      lag_duration3 = lag(ln.duration, 3), lag_duration4 = lag(ln.duration, 4),
                                                                      duration_dif = lag(ln.duration, 1) - ln.duration)

lag_fss_learning_longer <- lag_fss_learning %>% dplyr::select(ID, flow, lag_flow1, lag_flow2, lag_flow3, lag_flow4) %>%
  gather(key=lag_amount, value=lagged_flow_score, lag_flow1, lag_flow2, lag_flow3, lag_flow4)

#Plotting
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

#ggsave("lagged_flow.pdf", width=7, height=6) #Uncomment to save figure as pdf


############
#Multilevel Confirmatory Factor Analysis (CFA)
############

twolevel <- '
level: 1
  f1 =~ Q2 + Q4 + Q5 + Q7 + Q8 + Q9
  f2 =~ Q1 + Q3 + Q6 + Q10
level: 2
  fb1 =~ Q2 + Q4 + Q5 + Q7 + Q8 + Q9
  fb2 =~ Q1 + Q3 + Q6 + Q10
 
 #alternatively
  # Q1 ~~ Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10
  # Q2 ~~ Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10
  # Q3 ~~ Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10
  # Q4 ~~ Q4+Q5+Q6+Q7+Q8+Q9+Q10
  # Q5 ~~ Q5+Q6+Q7+Q8+Q9+Q10
  # Q6 ~~ Q6+Q7+Q8+Q9+Q10
  # Q7 ~~ Q7+Q8+Q9+Q10
  # Q8 ~~ Q8+Q9+Q10
  # Q9 ~~ Q9+Q10
  # Q10 ~~ Q10
'

fss2 <- fss %>% mutate(session = rep(rep(1:8, each = 5), 18))
fss2 <- fss2[-which(fss2$session==6 & fss2$run==1 & fss2$participant==8),]
results <- cfa(twolevel, data = fss2[,c(1,4:13)], cluster = 'participant')
summary(results, fit.measures = T, standardized = T)

#Factor scores
factor_scores <- as.data.frame(lavPredict(results))


############
#Scatterplot matrix of game variables:
############

game_data %>%
  dplyr::select(cumrun, min_velocity, max_velocity, end_velocity, avg_velocity, collisions, speed_drops, duration) %>% 
  pairs.panels(method = "spearman", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE,
             smoother = TRUE) # show correlation ellipses


############
#Cronbach's alpha calculator for each trial. 1= full scale, 2 = fluency, 3 = absorption
############

alpha_calculator <- function(subscale=1) {
  cronbachs_alphas <- cbind(fss2, fss_learning_full)
  #cronbachs_alphas <- cbind(fss2[-306,], fss_learning_full) #remove one row from participant 8 (session 6 run 1), this row is missing in fss_learning (game data). NOT NECESSARY, WAS DONE EARLIER!
  #coud also use: fss2[-which(fss2$session==6 & fss2$run==1 & fss2$participant==8),] to double check row is correct
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
            scale_x_continuous(name = NULL, breaks=c(0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)) +
            ylab("Count") + labs(title="Full FSS scale") + scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9), limits=c(0,9)) +
            #coord_cartesian(xlim = c(0.55, max(alphas))) +
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
             scale_x_continuous(name = NULL, breaks=c(0.55, 0.65, 0.75, 0.85, 0.95)) +
             ylab(NULL) + labs(title="Fluency of performance") + scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9), limits=c(0,9)) +
             #coord_cartesian(xlim = c(0.55, max(alphas))) +
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
             scale_x_continuous(name = NULL, breaks=c(0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)) +
             ylab(NULL) + labs(title="Absorption by activity") + scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9), limits=c(0,9)) +
             #coord_cartesian(xlim = c(0.55, max(alphas))) +
             theme_classic(base_size=12))
  }
}

#all in the same figure
alpha_reliability <- ggarrange(alpha_calculator(1), alpha_calculator(2), alpha_calculator(3), 
          font.label = list(size = 10, color = "black", face = "bold"), 
          common.legend = TRUE, legend="bottom", hjust = -0.25, ncol = 3, nrow = 1)

alpha_reliability <- annotate_figure(alpha_reliability, 
                                     bottom = text_grob("Cronbach's alpha", color = "black", size=14),
                                     fig.lab = NULL, fig.lab.face = NULL)

#ggsave("histogram_alpha.pdf", width=10, height=5) #uncomment to save figure as pdf


############
#Multilevel reliability (not reported in paper)
############

mlr(cronbachs_alphas, grp = "ID", Time = "run_session", items = c("Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9", "Q10"),alpha=TRUE,icc=FALSE, aov=TRUE,
    lmer=FALSE,lme = TRUE, long=FALSE,values=NA,na.action="na.omit",plot=TRUE,
    main="Lattice Plot by subjects over time")


############
#Visualize spread of measurement sessions in time
############

game_data2 <- game_data %>% group_by(participant, session) %>% dplyr::summarize(session = mean(session), date=mean.Date(as.Date(date))) %>%
  mutate(year=as.factor(ifelse(as.integer(participant)>9, 2019, 2017)))

#Tricky way to make the plotting work: here we substitute year 2017 with year 2019 so that free-scale faceting works better
game_data2$date <- sub("2017", "2019", game_data2$date)
game_data2$date <- as.Date(game_data2$date)

game_data2 <- game_data2 %>% mutate(physiology = as.factor(ifelse(session == 2 | session == 3 | session == 4, 0, 1)))

#Create separate dataframes to call for geom_text (enables drawing two different texts for facetet grids)
data.year2017 <- data.frame(participant=6,date=as.Date("2020-1-6"),session=NA,year=as.factor(2017))
data.year2019 <- data.frame(participant=6,date=as.Date("2020-1-6"),session=NA,year=as.factor(2019))

#Add gap between measurement dates
game_data2 <- game_data2 %>% group_by(participant) %>% dplyr::mutate(gap = as.integer(date - lag(date, 1)), mean_gap = mean(gap, na.rm=TRUE)) #note that with lag(), must use dplyr::mutate()

ggplot(game_data2, aes(participant, date, group=factor(session))) + geom_point(aes(color=physiology), size=2.5, alpha=.8) + coord_flip() + 
  facet_wrap("year") + #use scales="free_x" if year 2017 is NOT replaced by 2019
  scale_y_date(date_breaks = "weeks" , date_labels = "%b-%d") +
  scale_color_discrete(name=NULL, labels=c("No physiologial measures", "Physiological measures")) +
  theme(axis.text.x = element_text(angle = -60, hjust = 0),
        axis.text.y = element_blank(), #hide participant numbers for anonymity
        panel.background = element_rect(fill = "white",
                                        colour = "lightgrey",
                                        size = 0.5),
        legend.key=element_blank(),
        legend.position="bottom",
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        legend.text = element_text(size=12),
        panel.border = element_blank(),
        panel.grid.major = element_line(size = 0.25, colour = "grey"),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="lightgray"),
        strip.text.x = element_text(size=12)) + 
  geom_hline(yintercept = as.Date("2020-01-1"), linetype="dashed", color ="blue", size=0.5) +
  geom_text(aes(label=session), hjust=0, vjust=0, size=3) +
  geom_text(data=data.year2017, aes(x=participant, y=date), inherit.aes=FALSE, label="2018", 
            color="blue", size=3) +
  geom_text(data=data.year2019, aes(x=participant, y=date), inherit.aes=FALSE, label="2020", 
            color="blue", size=3) +
  ylab(NULL) + xlab("Participants")

#ggsave("measurement_dates.pdf", width=10, height=5) #uncomment to save figure as pdf


############
#Check if physiological measures affect results
############

fss_learning_full_phys <- fss_learning_full %>% mutate(physiology_dich = as.factor(ifelse(session == 2 | session == 3 | session == 4, 0, 1)))

#replication, using only 2017 data, controlling was physiology vs not:
summary(lmer(flow ~ deviation + physiology_dich + (1|participant), data=subset(fss_learning_full_phys, as.numeric(ID)>9)))

#main finding, using combined datasets:
summary(lmer(flow ~ deviation*cumrun + physiology_dich + (1|participant), data=fss_learning_full_phys))

#association with flow:
physiology_model <- lmer(flow ~ physiology_dich + (1|participant), data=fss_learning_full_phys)

#plots (ugly ones):
ggplot(fss_learning_full_phys, aes(deviation, flow, colour=physiology_dich)) + geom_smooth(method="lm", se=FALSE, size=0.4) + geom_point(alpha=.4, size=0.5) + facet_wrap("ID") 

ggplot(fss_learning_full_phys, aes(physiology_dich, flow)) + geom_point(alpha=.3) + 
  stat_summary(fun = median, geom = "point",
               size = 2, color="firebrick")

############
#Effect of perceived competence
############

#Perceived competence increases session by session for participants 2-6, 10, 13-14, 16-17, i.e. 10/18, 55% of the participants
#Overall trend and effect is positive.
fss_learning_tempo <- fss_learning %>% dplyr::select(session, ID, comp2) %>% na.omit()
fss_learning_tempo$ID <- factor(fss_learning_tempo$ID, levels=c(1:18))

#Sort participants according to slope values
estimates_slope <- list()
counter = 0
for (x in 1:length(unique(fss_learning_tempo$ID))) { #iterate through each ID
  estimates_temp <- subset(fss_learning_tempo, ID==x) #subset one ID at a time
  temp_model <- lm(comp2 ~ session, data=estimates_temp)  #fit model for subsetted ID
  counter = counter + 1 #counter iteration for indexing
  estimates_slope[[counter]] <- coef(summary(temp_model))[2,1] #add values to estimates_slope list
}
estimates_slope <- as.data.frame(rbind(unlist(estimates_slope))) #make list into dataframe
colnames(estimates_slope) <- c(1:18) #rename columns (from V1...V18 to 1...18)
estimates_slope <- sort(estimates_slope, decreasing=T) #sort according to values
part_levels <- as.integer(names(estimates_slope)) #turn into integer values
fss_learning_tempo$ID <- factor(fss_learning_tempo$ID, levels=part_levels) #reorder ID factor levels in data according to sorted integers

perc_comp_lmer1 <- lmer(comp2 ~ session + (1|ID), data=fss_learning_tempo)
perc_comp_lmer2 <- lmer(comp2 ~ session + (session|ID), data=fss_learning_tempo)

ggplot(fss_learning_tempo, aes(session, comp2)) + geom_point(alpha=.45) +
  geom_line(aes(y=predict(perc_comp_lmer2), color="Random slope and intercept"), size=0.6) +
  geom_line(aes(y=predict(perc_comp_lmer1), color="Random intercept"), size=0.6) +
  labs(color=NULL) +
  facet_wrap("ID") + ylab("Perceived competence") + xlab("Session") + theme_bw() +
  theme(legend.position = c(x=.8, y=.1),
        legend.background = element_rect(fill="white"),
        legend.box.background = element_rect(colour = "black"),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=9),
        panel.background = element_rect(fill = "white",
                                        colour = "lightgrey",
                                        size = 0.5, linetype = "solid"),
        panel.border = element_blank(), 
        panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey"),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="lightgray"))
  
#ggsave("perceived_importance_session.pdf", width=7.5, height=6) #uncomment to save figure as pdf


############
#Participant-wise partial regression plots
############

###Effect of duration on flow controlling for deviation##
skill_deviation <- fss_learning %>% dplyr::select(session, ID, deviation, ln.duration, flow) %>% na.omit()
skill_deviation$ID <- factor(skill_deviation$ID, levels=c(1:18))

#Remove effect of deviation on flow
residuals_list <- list()
counter = 0
for (x in 1:length(unique(skill_deviation$ID))) { #iterate through each ID
  residuals_temp <- subset(skill_deviation, ID==x) #subset one ID at a time
  residuals_temp <- resid(lm(flow ~ deviation, data=residuals_temp))  #remove effect of deviation on flow
  counter = counter + 1 #counter iteration for indexing
  residuals_list[[counter]] <- residuals_temp #add values to residuals_list list
}

residuals_list <- as.data.frame(unlist(residuals_list)) #make list into dataframe
names(residuals_list)[1] <- "resid_flow_dev"

#Remove effect of deviation on ln.duration
residuals_list2 <- list()
counter = 0
for (x in 1:length(unique(skill_deviation$ID))) { #iterate through each ID
  residuals_temp <- subset(skill_deviation, ID==x) #subset one ID at a time
  residuals_temp <- resid(lm(ln.duration ~ deviation, data=residuals_temp))  #remove effect of deviation on flow
  counter = counter + 1 #counter iteration for indexing
  residuals_list2[[counter]] <- residuals_temp #add values to residuals_list list
}

residuals_list2 <- as.data.frame(unlist(residuals_list2)) #make list into dataframe
names(residuals_list2)[1] <- "resid_dur_dev"


###Effect of deviation on flow controlling for duration###
residuals_list3 <- list()
counter = 0
for (x in 1:length(unique(skill_deviation$ID))) { #iterate through each ID
  residuals_temp <- subset(skill_deviation, ID==x) #subset one ID at a time
  residuals_temp <- resid(lm(flow ~ ln.duration, data=residuals_temp))  #remove effect of ln.duration on flow
  counter = counter + 1 #counter iteration for indexing
  residuals_list3[[counter]] <- residuals_temp #add values to residuals_list list
}

residuals_list3 <- as.data.frame(unlist(residuals_list3)) #make list into dataframe
names(residuals_list3)[1] <- "resid_flow_dur"

#Remove effect of deviation on ln.duration
residuals_list4 <- list()
counter = 0
for (x in 1:length(unique(skill_deviation$ID))) { #iterate through each ID
  residuals_temp <- subset(skill_deviation, ID==x) #subset one ID at a time
  residuals_temp <- resid(lm(deviation ~ ln.duration, data=residuals_temp))  #remove effect of deviation on flow
  counter = counter + 1 #counter iteration for indexing
  residuals_list4[[counter]] <- residuals_temp #add values to residuals_list list
}

residuals_list4 <- as.data.frame(unlist(residuals_list4)) #make list into dataframe
names(residuals_list4)[1] <- "resid_dev_dur"

#Combine everything into a dataframe
fss_residualtest <- cbind(fss_learning, residuals_list, residuals_list2, residuals_list3, residuals_list4)

#Plot the effect of duration on flow controlling for deviation
partial_plot1 <- ggplot(fss_residualtest, aes(resid_dur_dev, resid_flow_dev)) + geom_point(alpha=.35) + geom_smooth(method="lm", se=FALSE) + 
  facet_wrap("participant") + theme_bw() + ylab("Flow score") + xlab("log(Duration)") +
  labs(title="Controlling for Deviation scores")

#Plot the effect of deviation on flow controlling for deviation
partial_plot2 <- ggplot(fss_residualtest, aes(resid_dev_dur, resid_flow_dur)) + geom_point(alpha=.35) + geom_smooth(method="lm", se=FALSE) + 
  facet_wrap("participant") + theme_bw() + ylab(NULL) + xlab("Deviation scores") +
  labs(title="Controlling for log(Duration)")

ggarrange(partial_plot1, partial_plot2, 
          font.label = list(size = 10, color = "black", face = "bold"), 
          common.legend = TRUE, legend="bottom", hjust = -0.25, ncol = 2, nrow = 1)

#ggsave("partial_regression.pdf", width=12, height=7) #uncomment to save figure as pdf


############
#Main analyses controlling for covariates
############

fss_learning_covariates <- cbind(fss_learning, background_flowcar)
fss_learning_covariates$ID <- factor(fss_learning_covariates$ID, levels=1:18) #relevel factor to enable using ">" on subsetting

fss_learning_covariates2 <- subset(fss_learning_covariates, cumrun!=1)
fss_learning_covariates2$ID <- factor(fss_learning_covariates2$ID, levels=1:18)


#f~d effect (only new data; hypothesis 1)
fss_learning_lmer_h1a <- lmer(flow ~ deviation + age + gender + driving_experience + gaming_experience + (deviation|participant), 
                                data=subset(fss_learning_covariates, as.numeric(ID) > 9))

#learning effect (only new data; hypothesis 1)
fss_learning_lmer_h1b <- lmer(ln.duration ~ ln.cumrun + age + gender + driving_experience + gaming_experience + (ln.cumrun|participant), 
                                      data=subset(fss_learning_covariates, as.numeric(ID) > 9))

#moderating effect of experience (hypothesis 2)
fss_learning_lmer_h2 <- lmer(flow ~ deviation*cumrun + age + gender + driving_experience + gaming_experience + (1|participant), data=fss_learning_covariates)

#run duration vs. deviation
fss_learning_lmer_add1 <- lmer(flow ~ ln.duration + age + gender + driving_experience + gaming_experience + (1|participant), data=fss_learning_covariates)
fss_learning_lmer_add2 <- lmer(flow ~ deviation + ln.duration + age + gender + driving_experience + gaming_experience + (1|participant), data=fss_learning_covariates)


############
#Miscellaneous
############

#fluency and absorption mean, median, range (for table)
fss_learning %>% group_by(participant) %>% dplyr::summarize(fluency_m=mean(fluency), fluency_med=median(fluency), fluency_range=paste(as.character(min(fluency)), as.character(max(fluency)), sep=" - "), 
                                                            absorption_m=mean(absorption), absorption_med=median(absorption), absorption_range=paste(as.character(min(absorption)), as.character(max(absorption)), sep=" - "))

#summary stats, age, gaming and driving experience (for table)
fss_learning_full %>% group_by(participant) %>% dplyr::summarize(age = mean(age), gender = mean(as.numeric(gender)), driving_experience = mean(driving_experience), gaming_experience = mean(gaming_experience))

