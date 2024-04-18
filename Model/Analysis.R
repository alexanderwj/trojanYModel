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
numSimulations <- 5
numPlots <- 5

#simulation----

#FUNCTION CALL: simulate(Carrying Capacity, Number Myy, Number Fyy, YY Relative Survival Rate,
#                        Suppression Level, Number of Simulations, Number of Plots)


#print(results <- (simulate(K,numMyy,numFyy,yyRelSurvival,movingFish,suppressionLevel,numSimulations,numPlots)))

results <- data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "Eliminated", "Years", "MinFemales","EndPop","EndBiomass","MovingFish"))))

#results <- (simulate(K,numMyy,numFyy,yyRelSurvival,movingFish,suppressionLevel,numSimulations,numPlots))

inds <- (simulate(K,0,0,1,0,0,1,1))

vec1 <- c(0.4,0.6,0.8,1)
vec2 <- seq(0.1,1,0.1)


results <- (simulate(66430,13286,0,1,0,0,200,50))
write.csv(results,"base_lp_0.2k_200.csv")
#write.csv(results,"base_lp.csv")
results <- data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "Eliminated", "Years", "MinFemales","EndPop","EndBiomass","MovingFish"))))
for (i in 1:length(vec1)) {
  print(round(6643*vec1[i],0))
  tresults <- (simulate(6643,round(6643*vec1[i],0),0,1,0,0,1000,50))
  results <- rbind(results,tresults)
}
write.csv(results,paste("base_hs_lp_sr2.csv",sep=''))



# for (i in 1:length(vec1)) {
#   for (j in 1:length(vec2)) {
#     tresults <- (simulate(6643,round(6643*vec2[j],0),0,1,round(6643*vec1[i],0),0,5,1))
#     results <- rbind(results,tresults)
#     print(paste("Stocked: ",round(vec2[j]*6643,0),", Moving: ",round(vec1[i]*6643,0),sep=''))
#   }
# }
# write.csv(results,paste("base_hs_lp_pctK_movement_allM3.csv",sep=''))
# results <- data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "Eliminated", "Years", "MinFemales","EndPop","EndBiomass","MovingFish"))))
# for (i in 1:length(vec1)) {
#   for (j in 1:length(vec2)) {
#     tresults <- (simulate(6643,round(6643*vec2[j]/2,0),round(6643*vec2[j]/2,0),1,round(6643*vec1[i],0),0,5,1))
#     results <- rbind(results,tresults)
#     print(paste("Stocked: ",round(vec2[j]*6643,0),", Moving: ",round(vec1[i]*6643,0),sep=''))
#   }
# }
# write.csv(results,paste("base_hs_lp_pctK_movement_half3.csv",sep=''))

# 
# results <- data.frame()

# for (i in 1:length(vec5)) {
#   tresults <- (simulate(6643,round(6643*vec5[i],0),0,1,0,0,10,1))
#   results <- rbind(results,tresults)
#   print(c(round(6643*vec5[i],0)))
# }
# write.csv(results,paste("base_hs_hp_k_sr_5.csv",sep=''))