### Code was influenced by https://www.statmethods.net/stats/power.html and source code from https://cran.r-project.org/web/packages/wmwpow/wmwpow.pdf
### Wilcox rank sum test simulations of a unimodal distribution against a bimodal distribution

### reset R workspace
rm(list = ls())
if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
'%!in%' <- function(x, y)
  ! ('%in%'(x, y))
set.seed(1)

setwd("./") ### set directory
### check if packages are installed
list.of.packages <- c("tidyverse", "ggpubr","cowplot")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

library(tidyverse) ### Data manipulation
library(ggpubr) ### used to put stat test (wilcox) on plots
library(cowplot) ### for graph purposes

### Used for visualization purposes of the problem
### efficiency = 50%
### d = 0.5
### number of cells in each sample = 10,000
### cohen's d, effect size d = (u1 - u2)/s, assumption is that s1 = s2 in bimodal sample

n <- 10000
eff <- 0.50
d <- 0.5
sd <- 1 / d

y1 <- rnorm(n, 0, sd) ### normal distribution around 0
y2 <- rnorm(n, 1, sd) ### normal distribution around 1
w <-
  rbinom(n, 1, eff) ### number of successes given an efficiency, picks approx eff
dist1 <-
  rnorm(n, 0, sd) ### cell population 1, wild-type, same distribution as y1
dist2 <-
  c(y1[w == 0], y2[w == 1]) ### cell population 2, ko, mixed distribution of y1 & y2

### visualizing knockout cell population
wt_sub <- y1[w == 0]
ko_sub <- y2[w == 1]

temp <- data.frame(wt_sub)
temp2 <- data.frame(ko_sub)
colnames(temp)[1] <- "Meas"
colnames(temp2)[1] <- "Meas"
temp$Type <- "WT_Phenotype"
temp2$Type <- "KO_Phenotype"
temp <- rbind(temp, temp2)
title <-
  paste(
    "Knock-Out Cell Population Comparison \nCells(n) = ",
    n,
    " Edit Eff = ",
    eff,
    " Effect Size = ",
    d
  )

pdf("./Density_Knockout_Population.pdf",
    height = 10,
    width = 12)
temp %>% ggplot(aes(x = Meas, fill = Type)) +
  geom_density(alpha = 0.4) +
  xlab("Phenotype Measurement") +
  ylab("Density") +
  ggtitle(title)
dev.off()

pdf("./Density_Knockout_Population_Combined.pdf",
    height = 10,
    width = 10)
temp %>% ggplot(aes(x = Meas)) +
  geom_density(alpha = 0.4) +
  xlab("Phenotype Measurement") +
  ylab("Density") +
  ggtitle(title)
dev.off()

### visualizing wilcox comparison

dat_WT <- data.frame(dist1)
dat_KO <- data.frame(dist2)
dat_WT$Type <- "Wild-Type"
dat_KO$Type <- "Knock-Out"
colnames(dat_WT)[1] <- "Meas"
colnames(dat_KO)[1] <- "Meas"
dat <- rbind(dat_WT, dat_KO)

title <-
  paste("Cell Population Comparison \nCells(n) = ",
        n,
        " Edit Eff = ",
        eff,
        " Effect Size = ",
        d)
pdf("./Violin_Wilcox_Comparison.pdf",
    height = 10,
    width = 12)
dat %>% ggplot(aes(x = Type, y = Meas, fill = Type)) +
  geom_violin(alpha = 0.4) +
  stat_compare_means(method = "wilcox") +
  ylab("Phenotype Measurement") +
  xlab("") +
  ggtitle(title)
dev.off()

###simulation function, feed in the efficiency, number of samples, and standard deviation calculated from effect size
wilcox_sim <- function(eff, n, sd) {
  y1 <- rnorm(n, 0, sd) ### normal distribution around 0
  y2 <- rnorm(n, 1, sd) ### normal distribution around 1
  w <-
    rbinom(n, 1, eff) ### select randomly from population based on cutting efficiency
  dist1 <- rnorm(n, 0, sd) ### cell population 1 (wild-type)
  dist2 <-
    c(y1[w == 0], y2[w == 1]) ### cell population 2 (KO), mixed distribution
  dist1 <-
    round(dist1, 5) ### rounding the values due to computation processing
  dist2 <- round(dist2, 5)
  ### getting the p-value from the wilcox test
  p <-
    wilcox.test(
      dist1,
      dist2,
      paired = FALSE,
      correct = FALSE,
      alternative = "two.sided",
      exact = FALSE
    )$p.value
  return(p) ### return p-value
}

edit_effs <- c(.20, .50, .80, .90) ### tested editing efficiencies
effect_size <-
  seq(0.05, 1.0, 0.05) ### tested effect sizes of bimodal distribution from KO population

### empty vectors, used to store tested efficiencies (effs),
### effect sizes (ds), empircal power (emp_ps),
### and overlap percentages - overlapping coefficient (overlap)

effs <- c()
ds <- c()
emp_ps <- c()
overlap <- c()
n <- 10000

for (eff in edit_effs) {
  ### for loop through efficiencies
  for (d in effect_size) {
    ### for loop through effect sizes
    sd <- 1 / d
    
    ###This is an alternative for loop, if replicate is too computationally expensive to run
    # pvals <- c()
    # for (i in 1:1000){
    #   print(paste(eff,d,i))
    #   p <- wilcox_sim(eff,n,sd)
    #   pvals <- c(pvals,p)
    # }
    
    print(paste(eff, d))
    ### only running once for example purposes
    pvals <- replicate(1, wilcox_sim(eff, n, sd))
    
    ### ran on a server 1000 times.
    #pvals <- replicate(1000,wilcox_sim(eff,n,sd))
    
    emp_power <- round(sum(pvals < 0.05) / length(pvals), 3)
    ### calculating emp_power
    ### power = probability to reject H0 given H1
    
    ### overlap coefficient, percentage of overlapping populations (alternative way to think of effect size)
    ### identified from: https://rpsychologist.com/calculating-the-overlap-of-two-normal-distributions-using-monte-carlo-integration
    ovl <- 2 * pnorm(-abs(d) / 2)
    
    ### storing the simulations
    overlap <- c(overlap, ovl)
    effs <- c(effs, eff)
    ds <- c(ds, d)
    emp_ps <- c(emp_ps, emp_power)
  }
}
### Saving the simulations into a data frame
#results <- data.frame(effs,ds,emp_ps,overlap)
#write_delim(results,"./wilcox_mix_dist_simulations.txt")

power_data <- read.csv("./wilcox_mix_dist_simulations.txt", sep = "")
power_data$Editing_Efficiency <- as.character(power_data$effs) ###Setting as a character for plot
pdf("./Power_Wilcox_Mix_Dist.pdf",height = 10,width = 12)
power_data %>% ggplot(aes(x = ds, y = emp_ps, color = Editing_Efficiency)) +
  geom_line(alpha = 0.8,size = 1.2) +
  xlab("Effect Size in KO Sample") +
  ylab("Empirical Power")
dev.off()

