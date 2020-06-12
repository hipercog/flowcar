#Preliminary flow data analyses (c) Jussi Palomäki

#Load libraries
library(ggplot2)
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

#Create integer variable for participants (ID)
fss_learning <- fss_learning %>% mutate(ID=as.integer(participant))
fss_learning$ID <- as.factor(fss_learning$ID)

#RQ: Does the effect of Deviation on Flow (main finding from previous paper) depend on level cumrun (how many runs the participant has played)

#Remove first row (cumrun == 1) for each participant. NOT USED IN THE INTERACTION ANALYSIS.
fss_learning_sub <- subset(fss_learning, cumrun!=1)
anova(lmer(flow ~ deviation + (deviation|participant), data=fss_learning_sub))
anova(lmer(flow ~ deviation + (deviation|participant), data=fss_learning))

#Fit model (note: could also try random slope/intercept for deviation*cumrun):
fss_learning2_lmer <- lmer(flow ~ deviation*cumrun + (deviation|participant), data=fss_learning)

#Model R^2:
r.squaredGLMM(fss_learning2_lmer)

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

#Grab p-values, sort from lowst to highest, and adjust using Bonferroni-Holm
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
  geom_abline(data = slope_intercept_dataframes, aes(intercept=intercept_minusSD, slope=slope_minusSD, linetype="-1 SD (at 8.81 runs)"), size=0.8) + 
  geom_abline(data = slope_intercept_dataframes, aes(intercept=intercept_plusSD, slope=slope_plusSD, linetype="+1 SD (at 32.19 runs)"), size=0.8) + 
  labs(colour = NULL, linetype = "Cumulative runs", title = "Participant-wise models") + xlab("Deviation score") + ylab("Flow score") +
  theme_bw(base_size=14) + 
  theme(legend.position = c(x=.8, y=.1),
        legend.background = element_rect(fill="white"),
        legend.box.background = element_rect(colour = "black"),
        axis.text.x = element_text(size=8),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        strip.background=element_rect(fill="grey"))

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
panel_colors(.09)
panel_colors(.005)
panel_colors(.1)




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
