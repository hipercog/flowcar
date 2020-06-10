# ---
# title: "Flow data pipeline"
# author: "Ben Cowley (from T.Tammi)"
# note: This script depends on data wrangled in combine_data.R
# ---
# attach packages
library(lme4)
library(viridis)
library(here)
library(gghalves)
library(tidyverse)

source(file.path(here(), 'R', 'utils.R'))
# create an output dir for figures
FIGOUT <- FALSE
if (FIGOUT){
  odir <- file.path(here(), 'figures')
  dir.create(odir, showWarnings = FALSE)
}

# Basic violin plot
fss_items %>%
  mutate(participant = fct_reorder(participant, desc(flow), .fun='var')) %>%
  ggplot( aes(x=participant, y=flow)) + 
  geom_half_violin(trim=FALSE, fill="gray", side = "r") +
  geom_boxplot(width=0.1, outlier.shape = NA) +
  geom_hline(yintercept = median(fss_items$flow), color = "red") +
  geom_text(aes(length(unique(fss_items$participant)), median(fss_items$flow), label = "med.\nFlow", vjust = -1)) +
  labs(title="Flow scores per participant", x="Participant", y = "Flow") +
  ylim(2, 7) +
  theme_classic()
if (FIGOUT) ggsave(file.path(odir, "FlowXsubj.svg"))

# plot linear performance
ggplot(fss_learning, aes(cumrun, duration)) +
  geom_point(alpha=.6, size=2) +
  geom_smooth(method = "lm", se = FALSE, linetype = 1, size = 0.5, color="red") +
  facet_wrap(~participant) +
  theme_bw(base_size = 14)
if (FIGOUT)  ggsave(file.path(odir, "PerfXsubj.svg"))

# plot linear performance with power law fit
ggplot(fss_learning, aes(cumrun, duration)) +
  geom_point(alpha=.6, size=2) +
  stat_smooth(method = 'nls', formula = 'y~a*x^b', size = 0.5, se=FALSE, color ="red") +
  facet_wrap(~participant) +
  theme_bw(base_size = 14)
if (FIGOUT)  ggsave(file.path(odir, "PerfXsubj_powerlaw.svg"))


# plot linear performance with power law fit and flow z-scores coloring
ggplot(fss_learning, aes(cumrun, duration, color = z.flow)) +
  geom_point(alpha=.6, size=2) +
  stat_smooth(method = 'nls', formula = 'y~a*x^b', size = 0.5, se=FALSE, color ="red") +
  facet_wrap(~participant) +
  scale_color_viridis() +
  theme_bw(base_size = 14)
if (FIGOUT)  ggsave(file.path(odir, "PerfXsubj_powlxFlow.svg"))

# plot log-log with flow z-scores coloring
ggplot(fss_learning, aes(ln.cumrun, ln.duration, color = z.flow)) +
  geom_point(alpha=.6, size=2) +
  geom_smooth(method = "lm", se = FALSE, linetype = 1, size = 0.5, color="red") +
  facet_wrap(~participant) +
  scale_color_viridis() +
  theme_bw(base_size = 14)
if (FIGOUT)  ggsave(file.path(odir, "PerfXsubj_powlxFlow_loglog.svg"))

# plot deviation from predicted curve
ggplot(fss_learning, aes(deviation, flow)) +
  geom_point(alpha=.6, size=2) +
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~participant) +
  theme_bw(base_size = 14)
if (FIGOUT)  ggsave(file.path(odir, "FlowXdevXsubj.svg"))


### statistical test (linear mixed model) ----
fss_learning_lmer <- lmer(flow ~ deviation + (deviation|participant), data=fss_learning)
summary(fss_learning_lmer) # model summary
plot(fss_learning_lmer) # model diagnostics
qqnorm(residuals(fss_learning_lmer)) #qq-plot
qqline(residuals(fss_learning_lmer)) #line of "perfect normality"



### exploration with MWE confidence bands as documented in Korpela et al. (2014) ----
# Plots of all signals with mean sample means (thick solid line), their 95% MWEs (dashed lines), 
# and naive 95% quantiles (dotted thin lines) are shown below. The naive quantiles are computed per 
# time instance, without taking other time instances into account. In a way, a naive quantile 
# corresponds to an area where the “unadjusted p-values” are at least 0.05, while the MWE correspond 
# to area where the “adjusted p-values” (here adjusted taking the multiplicity correction due to 
# autocorrelation structure into account) are at least 0.05. The width of the naive quantile 
# typically provides a lower bound for the width of the respective MWE.
dat <- game_data %>%
  select(participant, cumrun, duration) %>%
  pivot_wider(names_from = cumrun, values_from = duration) %>%
  select(-participant)
idx <- which(apply(dat,2,function(x) all(!is.na(x))))

# Commented here is canonical example for finding curves of multiple datasets, e.g. subgroups
# curves <- lapply(data, function(x) findcurves(x))
# rng <- range(sapply(curves,function(x) range(x[,c("mean0","lo","up")])))
curve <- findcurves(as.matrix(dat[,idx]))
rng <- range(curve[,c("mean0","lo","up")])
plot(c(idx[1],idx[length(idx)]), rng
     , type="n", bty="n", xlab="cumrun", ylab="ms", main = "Group performance 95% CB")
plotlines(curve, idx)

curve17 <- findcurves(as.matrix(dat[1:9,idx]))
curve19 <- findcurves(as.matrix(dat[10:18,idx]))
plot(c(idx[1],idx[length(idx)]), rng
     , type="n", bty="n", xlab="cumrun", ylab="ms", main = "Group performance 95% CB")
plotlines_comp(curve17, idx, col="red")
plotlines_comp(curve19, idx, col="blue")
