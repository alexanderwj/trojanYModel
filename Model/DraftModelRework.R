#parameters----

# Simulation begins with age-2 fish, sex is randomly chosen
# K is implemented using the Beverton-Holt recruitment model (K=(R0-1)*M))
startingFish <- 500
k <- 80000

# Treatment years are when YY fish are stocked and suppression is applied, if applicable
burnInYears <- 25
treatmentYears <- 50
afterYears <- 25

# YY survival is calculated as a proportion of wild pikeminnow survival
#   (1=equivalent survival rate to wild fish)
wildMortality <- 0.57873786895
yySurvival <- 1

# YY fish are stocked at age 1
numFyy <- 2000
numMyy <- 0

# Suppression is size-selective based on WNRD 2022 efforts
# Suppression level = relative probability of a fish being suppressed at length l
#   (for any length l, double the level = double the probability of suppression)
# A level of 1 roughly corresponds to WNRD 2022 efforts (~500 removed at a pop. of ~80,000)
# Stocked fish cannot be suppressed
suppressionLevel <- 5

# Choose how many simulations will be run and how many will be plotted
# On the plots, black line = total population, red line = wild-type females
# You will get a report summarizing the results of all simulations
numSimulations <- 5
numPlots <- 5

#functions----

maturity <- function(inds) {
  lengths <- inds$length
  numFish <- length(lengths)
  probMature <- 1/(1+exp(log(19)*(lengths-250)/(250-308)))
  inds$mature<-rbinom(numFish, 1, probMature)
  inds
}

growth <- function(inds) {
  inds$age <- inds$age+1
  numFish <- length(inds$age)
  inds$length <- ifelse(inds$age > 1, 726.7*(1-exp(-.1437*(inds$age+.04422)))+rnorm(numFish, 0, 40), 726.7*(1-exp(-.1437*(inds$age+.04422)))+rnorm(numFish, 0, 14))
  inds
}

death <- function(inds) {
  numFish <- length(inds$age)
  inds$dead <- ifelse(inds$yy == 0, rbinom(numFish, 1, wildMortality), rbinom(numFish, 1, (1-((1-0.5701418464)*yySurvival))))
  inds <- subset(inds, dead==0)
  inds
}

birth <- function(inds) {
  matureFxx <- nrow(subset(inds, mature == 1 & sex == 1 & yy == 0))
  matureMxy <- nrow(subset(inds, mature == 1 & sex == 0 & yy == 0))
  matureFyy <- nrow(subset(inds, mature == 1 & sex == 1 & yy == 1))
  matureMyy <- nrow(subset(inds, mature == 1 & sex == 0 & yy == 1))
  
  totalPairs <- (min(c((matureFxx+matureFyy), (matureMxy+matureMyy))))
  if (totalPairs == 0) {
    return(inds)
  }
  spawners <- totalPairs*2
  
  newFish <- (((sqrt(k))*spawners)/(1+spawners/(sqrt(k))))/2+rnorm(1,0,k/10)
  if (newFish <= 0) {
     return(inds)
  }
  
  percMyy <- matureMyy/(matureMyy+matureMxy)
  percFyy <- matureFyy/(matureFyy+matureFxx)
  percBothyy <- percMyy*percFyy
  
  pairsBothyy <- rbinom(1,totalPairs,percBothyy)
  pairsFxxMyy <- rbinom(1,(totalPairs-pairsBothyy),percMyy)
  pairsFyyMxy <- rbinom(1,(totalPairs-pairsBothyy-pairsFxxMyy),percFyy)
  pairsWild <- (totalPairs-pairsBothyy-pairsFxxMyy-pairsFyyMxy)

  recruitsBothyy <- round((newFish*(pairsBothyy/totalPairs)), 0)
  recruitsFxxMyy <- round((newFish*(pairsFxxMyy/totalPairs)), 0)
  recruitsFyyMxy <- round((newFish*(pairsFyyMxy/totalPairs)), 0)
  recruitsWild <- (newFish-recruitsBothyy-recruitsFxxMyy-recruitsFyyMxy)

  if(recruitsBothyy > 0) {
    newFish_Bothyy <- data.frame(age=rep(0, recruitsBothyy), sex=rep(0, recruitsBothyy), length=rep(0, recruitsBothyy), mature=rep(0,recruitsBothyy), dead=rep(0, recruitsBothyy), yy=rep(1, recruitsBothyy), stocked=rep(0, recruitsBothyy))
    inds <- rbind(inds, newFish_Bothyy);
  }

  if(recruitsFxxMyy > 0) {
    newFish_FxxMyy <- data.frame(age=rep(0, recruitsFxxMyy), sex=rep(0, recruitsFxxMyy), length=rep(0, recruitsFxxMyy), mature=rep(0, recruitsFxxMyy), dead=rep(0, recruitsFxxMyy), yy=rep(0, recruitsFxxMyy), stocked=rep(0, recruitsFxxMyy))
    inds <- rbind(inds, newFish_FxxMyy);
  }

  if(recruitsFyyMxy > 0) {
    newFish_FyyMxy <- data.frame(age=rep(0, recruitsFyyMxy), sex=rep(0, recruitsFyyMxy), length=rep(0, recruitsFyyMxy), mature=rep(0, recruitsFyyMxy), dead=rep(0, recruitsFyyMxy), yy=rbinom(recruitsFyyMxy, 1, 0.5), stocked=rep(0, recruitsFyyMxy))
    inds <- rbind(inds, newFish_FyyMxy);
  }

  if(recruitsWild > 0) {
    newFish_Wild<- data.frame(age=rep(0, recruitsWild), sex=rbinom(recruitsWild, 1, 0.5), length=rep(0, recruitsWild), mature=rep(0, recruitsWild), dead=rep(0, recruitsWild), yy=rep(0, recruitsWild), stocked=rep(0, recruitsWild))
    inds <- rbind(inds, newFish_Wild);
  }
  inds
}

stockYY <- function(inds){
  if (numMyy > 0) {new_Myy <- data.frame(age=rep(1, numMyy), sex=rep(0, numMyy), length=rep(0, numMyy), mature=rep(0, numMyy), dead=rep(0, numMyy), yy=rep(1, numMyy), stocked=rep(1,numMyy)); inds <- rbind(inds, new_Myy)}
  if (numFyy > 0) {new_Fyy <- data.frame(age=rep(1, numFyy), sex=rep(1, numFyy), length=rep(0, numFyy), mature=rep(0, numFyy), dead=rep(0, numFyy), yy=rep(1, numFyy), stocked=rep(1,numFyy)); inds <- rbind(inds, new_Fyy)}
  inds
}

suppress <-function(inds) {
  lengths <- inds$length
  numFish <- length(lengths)
  suppressProb <-  (0.00008526*1.025^(.8753*(lengths))+.01753)*suppressionLevel*0.1877
  suppressProb <- ifelse(suppressProb>1, 1, suppressProb)
  inds$dead <- rbinom(numFish, 1, suppressProb)
  inds$dead <- ifelse(inds$stocked==1, 0, inds$dead)
  inds <- subset(inds, dead==0)
  inds
}

#simulation----

plotYears <- sample.int(numSimulations, numPlots)
cat("\nSimulation runs to be plotted:", sort(plotYears), "\n")

yearsToElim <- c()
startTime <- Sys.time()

for (y in 1:numSimulations) {
  inds <- data.frame(age=rep(2, startingFish), sex=rbinom(startingFish,1,0.5), length=rep(0,startingFish), mature=rep(0,startingFish), dead=rep(0,startingFish), yy=rep(0,startingFish), stocked=rep(0,startingFish))
  eliminationYear <- 0
  Population <- c(startingFish)
  numFxx <- nrow(subset(inds, sex == 1 & yy == 0))
  
  for (year in 1:(burnInYears+treatmentYears+afterYears)) {
    if (year > burnInYears && year <= (burnInYears+treatmentYears)) {
      if (suppressionLevel>0) {inds <- suppress(inds)}
      inds <- stockYY(inds)
    }
    
    if (nrow(subset(inds, sex == 1 & yy == 0)) == 0) {
      Population <- append(Population, nrow(inds))
      numFxx <- append(numFxx, 0)
      eliminationYear <- year
      break
    }
    
    inds <- death(inds)
    inds <- growth(inds)
    inds <- maturity(inds)
    inds <- birth(inds)
    
    Population <- append(Population, nrow(inds))
    numFxx <- append(numFxx, nrow(subset(inds, sex == 1 & yy == 0)))
  }
  
  Year <- 0:(burnInYears+treatmentYears+afterYears)
  Year <- (head(Year,(length(Population))))
  numFxx <- (head(numFxx,(length(Population))))
  if(eliminationYear != 0) {yearsToElim <- append(yearsToElim, eliminationYear-burnInYears)}
  
  if (is.element(y, plotYears)) {
    plot(Year, Population, type='l', main=(ifelse(eliminationYear == 0, (paste("Run", y, "- Not Extirpated. Min. Females:", min(tail(numFxx,(treatmentYears+afterYears))))), (paste("Run", y, "-", (eliminationYear-burnInYears), "years to extirpation")))), ylim=c(0,max(Population)))
    lines(Year, numFxx, type='l', col="red")
    abline(v=burnInYears, col="red")
    abline(v=burnInYears+treatmentYears, col="red")
  }
}
endTime <- Sys.time()
timeTaken <- round(difftime(endTime, startTime, units='secs'), digits=2)

#analysis----

cat("\nTotal Execution Time: ", timeTaken, " seconds (", round(timeTaken/numSimulations, digits=2), " seconds per run)", sep='')
cat("\nExtirpation Percentage: ", (((length(yearsToElim))/numSimulations)*100), "% (", length(yearsToElim), " runs out of ", numSimulations, ")\n", sep='')
if (is.null(yearsToElim) == FALSE) {cat("Of Extirpating Runs, Mean: ", (mean(yearsToElim)), " years, SD: ", (sd(yearsToElim)), ", Range: ", (range(yearsToElim)[1]), "-", (range(yearsToElim)[2]),  sep='')}