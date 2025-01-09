library(ggplot2)
library(dplyr)
#library(rowr)

setwd(getSrcDirectory(function(){})[1])
source("DraftModelRework.R")


#parameters----

# K is implemented using the Beverton-Holt difference equation (Bohner and Chieochan 2013)
K <- 6643

#YY fish are stocked annually at age 1 during the summer field season
numMyy <- 0
numFyy <- 0

#Movement: how many fish leave/enter the area each year, and what size must they be below?
movingFish <- 0
cutoffSize <- 999

# YY survival is calculated as a proportion of wild pikeminnow survival
#   (1=equivalent survival rate to wild fish)
yyRelSurvival <- 1

#Age of stocked YY fish
suppressionLevel <- 0

# Suppression is size-selective based on WNRD 2022 efforts
# Level is a multiplier of the selectivity function (2022 estimate: 0.18)
#   A level of 1 means the most-selected lengths have a 100% chance of being suppressed
# Stocked fish cannot be suppressed
stockedAge <- 0

# Choose how many simulations will be run and how many will be plotted
# On the plots, black line = total population, red line = wild-type females,
#   verticals dashed lines = management actions start/stop
numSimulations <- 5
numPlots <- 5

#simulation----

#FUNCTION CALL: simulate(K,Myy,Fyy,survival,movers,suppression,stockAge,simulations,plots) 

results <- (simulate(K,numMyy,numFyy,yyRelSurvival,movingFish,suppressionLevel,stockedAge,numSimulations,numPlots))

#print(results <- (simulate(K,numMyy,numFyy,yyRelSurvival,movingFish,suppressionLevel,numSimulations,numPlots)))
#results <- data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "Eliminated", "Years", "MinFemales","EndPop","EndBiomass","MovingFish"))))

