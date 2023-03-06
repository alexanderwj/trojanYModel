library(ggplot2)
source("DraftModelRework.R")

#parameters----

# K is implemented using the Beverton-Holt recruitment model (K=(R0-1)*M))
K <- 50000

#YY fish are stocked at age 1 during the summer field season
numFyy <- 0
numMyy <- 0

# YY survival is calculated as a proportion of wild pikeminnow survival
#   (1=equivalent survival rate to wild fish)
yyRelSurvival <- 1

# Suppression is size-selective based on WNRD 2022 efforts 
# Suppression level = relative probability of a fish being suppressed at length l
#   (for any length l, double the level = double the probability of suppression)
# Suppression probability scales inverse to population 
#   (1/2 the population = 2x the probability for each fish at the same level)
# A level of 1 roughly corresponds to WNRD 2022 efforts (~500 removed)
# Stocked fish cannot be suppressed
suppressionLevel <- 1

# Choose how many simulations will be run and how many will be plotted
# On the plots, black line = total population, red line = wild-type females
# You will get a report summarizing the results of all simulations
numSimulations <- 1
numPlots <- 0

#simulation----

results <- data.frame(matrix(ncol=8,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "Eliminated", "Years", "MinFemales"))))
suppressedLevels <- seq(0,1,0.2)
startTime <- Sys.time()

for(i in seq(1,5)) {
  for (suppressLevel in suppressedLevels) {
    results <- rbind(results,simulate(K,10000,2500,yyRelSurvival,suppressLevel,numSimulations,1))
  }
}

endTime <- Sys.time()
timeTaken <- round(difftime(endTime, startTime, units='secs'), digits=2)

print(results)
cat("\nTotal Execution Time: ", timeTaken, " seconds", sep='')
print(ggplot(results,aes(x=SuppressionLevel,y=Years,group=SuppressionLevel))+xlab("Suppression Level")+ylab("Years to Extirpation")+geom_boxplot())