source("DraftModelRework.R")

k <- 10000
numFyy <- 0
numMyy <- 5000
yySurvival <- 1
suppressionLevel <- 1
numSimulations <- 5
numPlots <- 5

results <- simulate(k,numMyy,numFyy,yySurvival,suppressionLevel,numSimulations,numPlots)

print(results)