library(ggplot2)
library(dplyr)
setwd(getSrcDirectory(function(){})[1])
source("DraftModelRework.R")

#parameters----

# K is implemented using the Beverton-Holt difference equation (Bohner and Chieochan 2013)
K <- 50000

#YY fish are stocked annually at age 1 during the summer field season
numMyy <- 3000
numFyy <- 500

#Movement: how many fish leave/enter the area each year, and what size must they be below?
movingFish <- 500
cutoffSize <- 300

# YY survival is calculated as a proportion of wild pikeminnow survival
#   (1=equivalent survival rate to wild fish)
yyRelSurvival <- 1

# Suppression is size-selective based on WNRD 2022 efforts
# Level is a multiplier of the selectivity function (2022 estimate: 0.18)
#   A level of 1 means the most-selected lengths have a 100% chance of being suppressed
# Stocked fish cannot be suppressed
suppressionLevel <- 0.5

# Choose how many simulations will be run and how many will be plotted
# On the plots, black line = total population, red line = wild-type females,
#   verticals dashed lines = management actions start/stop
numSimulations <- 1
numPlots <- 1

#simulation----

#FUNCTION CALL: simulate(Carrying Capacity, Number Myy, Number Fyy, YY Relative Survival Rate,
#                        Suppression Level, Number of Simulations, Number of Plots)

#results <- data.frame(matrix(ncol=8,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "Eliminated", "Years", "MinFemales"))))

print(simulate(K,numMyy,numFyy,yyRelSurvival,movingFish,suppressionLevel,numSimulations,numPlots))