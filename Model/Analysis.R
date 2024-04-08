library(ggplot2)
library(dplyr)
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

# Suppression is size-selective based on WNRD 2022 efforts
# Level is a multiplier of the selectivity function (2022 estimate: 0.18)
#   A level of 1 means the most-selected lengths have a 100% chance of being suppressed
# Stocked fish cannot be suppressed
suppressionLevel <- 0

# Choose how many simulations will be run and how many will be plotted
# On the plots, black line = total population, red line = wild-type females,
#   verticals dashed lines = management actions start/stop
numSimulations <- 10
numPlots <- 10

#simulation----

#FUNCTION CALL: simulate(Carrying Capacity, Number Myy, Number Fyy, YY Relative Survival Rate,
#                        Suppression Level, Number of Simulations, Number of Plots)


#print(results <- (simulate(K,numMyy,numFyy,yyRelSurvival,movingFish,suppressionLevel,numSimulations,numPlots)))

#results <- data.frame(matrix(ncol=8,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "Eliminated", "Years", "MinFemales"))))

inds <- (simulate(K,numMyy,numFyy,yyRelSurvival,movingFish,suppressionLevel,numSimulations,numPlots))

# vec1 <- c(0.7,0.75,0.8,0.9,1)
# vec2 <- c(0.5,0.6,0.7,0.75,0.8,0.9,1)
# vec3 <- c(0.6,0.7,0.8,0.9,1)
# vec4 <- 6643*seq(0.0125,1,0.0125)
# vec5 <- c(seq(0,1,0.0125))
# 
# for (i in 1:length(vec4)) {
#   results <- (simulate(6643,round(6643*vec4[i],0),0,1,0,0,1000,50))
#   write.csv(results,paste("base_hs_lp_",as.integer(vec4[i]*100),"k.csv",sep=''))
#   print(vec4[i])
# }
# 
# results <- data.frame()
# for (i in 1:length(vec4)) {
#   for (j in 1:length(vec5)) {
#     tresults <- (simulate(6643,round(vec4[i]*vec5[j],0),round(vec4[i]*(1-vec5[j]),0),1,0,0,10,1))
#     results <- rbind(results,tresults)
#     print(c(vec4[i],round(vec4[i]*vec5[j],0),round(vec4[i]*(1-vec5[j]))))
#   }
# }
# write.csv(results,paste("base_hs_hp_k_sr_low.csv",sep=''))
# 
# for (i in 1:length(vec5)) {
#   tresults <- (simulate(6643,round(6643*vec5[i],0),0,1,0,0,10,1))
#   results <- rbind(results,tresults)
#   print(c(round(6643*vec5[i],0)))
# }
# write.csv(results,paste("base_hs_hp_k_sr_5.csv",sep=''))