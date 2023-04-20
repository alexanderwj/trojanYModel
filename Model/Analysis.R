library(ggplot2)
library(dplyr)
source("DraftModelRework.R")

#parameters----

# K is implemented using the Beverton-Holt difference equation ()
K <- 10000

#YY fish are stocked annually at age 1 during the summer field season
numMyy <- 350
numFyy <- 350

# YY survival is calculated as a proportion of wild pikeminnow survival
#   (1=equivalent survival rate to wild fish)
yyRelSurvival <- 1

# Suppression is size-selective based on WNRD 2022 efforts
# Level is a multiplier of the selectivity function 
#   (Level of 1 means the most-selected lengths have a 100% chance of being suppressed)
# Stocked fish cannot be suppressed
suppressionLevel <- 1

# Choose how many simulations will be run and how many will be plotted
# On the plots, black line = total population, red line = wild-type females,
#   verticals red lines = management actions start/stop
numSimulations <- 5
numPlots <- 5

#simulation----

#FUNCTION CALL: simulate(Carrying Capacity, Number Myy, Number Fyy, YY Relative Survival Rate,
#                        Suppression Level, Number of Simulations, Number of Plots)

#results <- data.frame(matrix(ncol=8,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "Eliminated", "Years", "MinFemales"))))

print(simulate(K,numMyy,numFyy,yyRelSurvival,suppressionLevel,numSimulations,numPlots))