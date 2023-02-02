#parameters----

#maturity parameters
lm50 <- 250
lm95 <- 308

#Von Bert parameters
lInf <- 726.7
rate <- .1437
tZero <- -.04422

#annual mortality rate A
wildMortality <- 0.57873786895

# Simulation begins with age-2 fish, sex is randomly chosen
startingFish <- 500

# Treatment years are when YY fish are stocked and suppression is applied, if applicable
burnInYears <- 25
treatmentYears <- 100
afterYears <- 25

#component functions----

maturity <- function(inds) {
  lengths <- inds$length
  numFish <- length(lengths)
  probMature <- 1/(1+exp(log(19)*(lengths-lm50)/(lm50-lm95)))
  inds$mature<-rbinom(numFish, 1, probMature)
  inds
}

growth <- function(inds) {
  inds$age <- inds$age+1
  numFish <- length(inds$age)
  inds$length <- ifelse(inds$age > 1, lInf*(1-exp(-rate*(inds$age-tZero)))+rnorm(numFish, 0, 40), lInf*(1-exp(-rate*(inds$age-tZero)))+rnorm(numFish, 0, 14))
  inds
}

death <- function(inds, survival) {
  numFish <- length(inds$age)
  inds$dead <- ifelse(inds$yy == 0, rbinom(numFish, 1, wildMortality), rbinom(numFish, 1, (1-((1-wildMortality)*survival))))
  inds <- subset(inds, dead==0)
  inds
}

birth <- function(inds) {
  matureFxx <- nrow(subset(inds, mature == 1 & sex == 1 & yy == 0))
  matureMxy <- nrow(subset(inds, mature == 1 & sex == 0 & yy == 0))
  matureFyy <- nrow(subset(inds, mature == 1 & sex == 1 & yy == 1))
  matureMyy <- nrow(subset(inds, mature == 1 & sex == 0 & yy == 1))
  
  totalPairs <- matureFxx+matureFyy
  if (totalPairs == 0) {
    return(inds)
  }
  spawners <- totalPairs*2
  
  newFish <- (((sqrt(k)+1)*spawners)/(1+spawners/(sqrt(k))))/2+rnorm(1,0,k/10)
  if (newFish <= 0) {
     return(inds)
  }
  
  percMyy <- ifelse((matureMyy+matureMxy)==0,0,matureMyy/(matureMyy+matureMxy))
  percFyy <- ifelse((matureFyy+matureFxx)==0,0,matureFyy/(matureFyy+matureFxx))
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

stockYY <- function(inds,Myy,Fyy){
  if (Myy > 0) {new_Myy <- data.frame(age=rep(1, Myy), sex=rep(0, Myy), length=rep(0, Myy), mature=rep(0, Myy), dead=rep(0, Myy), yy=rep(1, Myy), stocked=rep(1,Myy)); inds <- rbind(inds, new_Myy)}
  if (Fyy > 0) {new_Fyy <- data.frame(age=rep(1, Fyy), sex=rep(1, Fyy), length=rep(0, Fyy), mature=rep(0, Fyy), dead=rep(0, Fyy), yy=rep(1, Fyy), stocked=rep(1,Fyy)); inds <- rbind(inds, new_Fyy)}
  inds
}

suppress <-function(inds,suppression) {
  lengths <- inds$length
  numFish <- length(lengths)
  suppressProb <-  (0.00008526*1.025^(.8753*(lengths))+.01753)*suppression*0.1877*(1/(nrow(inds)/80000))
  suppressProb <- ifelse(suppressProb>1, 1, suppressProb)
  inds$dead <- ifelse(inds$stocked==1, 0, rbinom(numFish, 1, suppressProb))
  inds <- subset(inds, dead==0)
  inds
}

#simulation function----

simulate <- function(k,Myy,Fyy,survival,suppression,simulations,plots) {
  plotYears <- sample.int(numSimulations, plots)
  #cat("\nSimulation runs to be plotted:", sort(plotYears), "\n")
  results <- data.frame(matrix(ncol=8,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "Eliminated", "Years", "MinFemales"))))
  
  for (y in 1:simulations) {
    inds <- data.frame(age=rep(2, startingFish), sex=rbinom(startingFish,1,0.5), length=rep(0,startingFish), mature=rep(0,startingFish), dead=rep(0,startingFish), yy=rep(0,startingFish), stocked=rep(0,startingFish))
    eliminationYear <- 0
    Population <- c(startingFish)
    numFxx <- nrow(subset(inds, sex == 1 & yy == 0))
    
    for (year in 1:(burnInYears+treatmentYears+afterYears)) {
      if (year > burnInYears && year <= (burnInYears+treatmentYears)) {
        if (suppressionLevel>0) {inds <- suppress(inds,suppression)}
        inds <- stockYY(inds,Myy,Fyy)
      }
      
      if (nrow(subset(inds, sex == 1 & yy == 0)) == 0) {
        yearResults <- data.frame(K=k,Myy=Myy,Fyy=Fyy,YYSurvival=survival,SuppressionLevel=suppression,Eliminated=1,Years=year-burnInYears,MinFemales=0)
        results <- rbind(results, yearResults)
        Population <- append(Population, nrow(inds))
        numFxx <- append(numFxx, 0)
        eliminationYear <- year
        break
      }
      
      inds <- death(inds,survival)
      inds <- growth(inds)
      inds <- maturity(inds)
      inds <- birth(inds)
      
      Population <- append(Population, nrow(inds))
      numFxx <- append(numFxx, nrow(subset(inds, sex == 1 & yy == 0)))
    }
    
    Year <- 0:(burnInYears+treatmentYears+afterYears)
    Year <- (head(Year,(length(Population))))
    numFxx <- (head(numFxx,(length(Population))))
    if(eliminationYear == 0) {yearResults <- data.frame(K=k,Myy=Myy,Fyy=Fyy,YYSurvival=survival,SuppressionLevel=suppression,Eliminated=0,Years=NA,MinFemales=min(tail(numFxx,(treatmentYears+afterYears)))); results <- rbind(results, yearResults)}
    
    if (is.element(y, plotYears)) {
      plot(Year, Population, type='l', main=(ifelse(eliminationYear == 0, (paste("Run", y, "- Not Extirpated. Min. Females:", min(tail(numFxx,(treatmentYears+afterYears))))), (paste("Run", y, "-", (eliminationYear-burnInYears), "years to extirpation")))), ylim=c(0,max(Population)))
      lines(Year, numFxx, type='l', col="red")
      abline(v=burnInYears, col="red")
      abline(v=burnInYears+treatmentYears, col="red")
    }
  }
  results
}