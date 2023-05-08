#parameters----

#maturity parameters: these are from McCormick et al. 2021 for common carp
lm50 <- 250
lm95 <- 308

#Von Bert parameters: these are from aged scales collected by WNRD in 2019 and 2022
lInf <- 726.7
rate <- .1437
tZero <- -.04422

#annual mortality rate A for wild fish: calculated using weighted catch-curve analysis of WNRD snorkel surveys
#wildMortality <- 0.57873786895
wildMortality <- .605741101

# Simulation begins with age-2 fish, sex is randomly chosen
startingFish <- 1000

# Treatment years are when YY fish are stocked and suppression is applied, if applicable
burnInYears <- 50
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

birth <- function(inds,K) {
  matureFxx <- nrow(subset(inds, mature == 1 & sex == 1 & yy == 0))
  matureMxy <- nrow(subset(inds, mature == 1 & sex == 0 & yy == 0))
  matureFyy <- nrow(subset(inds, mature == 1 & sex == 1 & yy == 1))
  matureMyy <- nrow(subset(inds, mature == 1 & sex == 0 & yy == 1))
  
  totalPairs <- matureFxx+matureFyy
  if (totalPairs == 0) {
    return(inds)
  }
  spawners <- totalPairs*2
  newFish <- ((20*K*spawners)/(K+(19*spawners)))*exp(rnorm(1,0,0.2))
  if (newFish <= 0) {
     return(inds)
  }
  
  percMyy <- ifelse((matureMyy+matureMxy)==0,0,matureMyy/(matureMyy+matureMxy))
  percFyy <- ifelse((matureFyy+matureFxx)==0,0,matureFyy/(matureFyy+matureFxx))
  
  pairsBothyy <- rbinom(1,totalPairs,percMyy*percFyy)
  pairsWild <- rbinom(1,totalPairs,(1-percFyy)*(1-percMyy))
  pairsFxxMyy <- rbinom(1,(totalPairs-pairsBothyy-pairsWild),(percMyy))
  pairsFyyMxy <- (totalPairs-pairsBothyy-pairsFxxMyy-pairsWild)

  recruitsBothyy <- round((newFish*(pairsBothyy/totalPairs)), 0)
  recruitsFxxMyy <- round((newFish*(pairsFxxMyy/totalPairs)), 0)
  recruitsFyyMxy <- round((newFish*(pairsFyyMxy/totalPairs)), 0)
  recruitsWild <- round((newFish-recruitsBothyy-recruitsFxxMyy-recruitsFyyMxy), 0)

  if(recruitsBothyy > 0) {
    newFish_Bothyy <- data.frame(age=rep(0, recruitsBothyy), sex=0, length=0, mature=0, dead=0, yy=1, stocked=0)
    inds <- rbind(inds, newFish_Bothyy);
  }

  if(recruitsFxxMyy > 0) {
    newFish_FxxMyy <- data.frame(age=rep(0, recruitsFxxMyy), sex=0, length=0, mature=0, dead=0, yy=0, stocked=0)
    inds <- rbind(inds, newFish_FxxMyy);
  }

  if(recruitsFyyMxy > 0) {
    newFish_FyyMxy <- data.frame(age=rep(0, recruitsFyyMxy), sex=0, length=0, mature=0, dead=0, yy=rbinom(recruitsFyyMxy, 1, 0.5), stocked=0)
    inds <- rbind(inds, newFish_FyyMxy);
  }

  if(recruitsWild > 0) {
    newFish_Wild<- data.frame(age=rep(0, recruitsWild), sex=rbinom(recruitsWild, 1, 0.5), length=0, mature=0, dead=0, yy=0, stocked=0)
    inds <- rbind(inds, newFish_Wild);
  }
  inds
}

stockYY <- function(inds,Myy,Fyy){
  if (Myy > 0) {new_Myy <- data.frame(age=rep(1, Myy), sex=0, length=0, mature=0, dead=0, yy=1, stocked=1); inds <- rbind(inds, new_Myy)}
  if (Fyy > 0) {new_Fyy <- data.frame(age=rep(1, Fyy), sex=1, length=0, mature=0, dead=0, yy=1, stocked=1); inds <- rbind(inds, new_Fyy)}
  inds
}

suppress <-function(inds,rate) {
  lengths <- inds$length
  numFish <- length(lengths)
  suppressProb <-  (0.00008526*1.025^(.8753*(lengths))+.01753)
  suppressProb <- ifelse(suppressProb>1, 1, suppressProb)
  suppressProb <- suppressProb*rate
  inds$dead <- ifelse(inds$stocked==1, 0, rbinom(numFish, 1, suppressProb))
  inds <- subset(inds, dead==0)
  inds
}

#simulation function----

simulate <- function(K,Myy,Fyy,survival,suppression,simulations,plots) {
  plotYears <- sample.int(simulations, plots)
  results <- data.frame(matrix(ncol=8,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "Eliminated", "Years", "MinFemales"))))
  
  for (y in 1:simulations) {
    inds <- data.frame(age=rep(2, startingFish), sex=rbinom(startingFish,1,0.5), length=0, mature=0, dead=0, yy=0, stocked=0)
    eliminationYear <- 0
    Population <- c(startingFish)
    numFxx <- nrow(subset(inds, sex == 1 & yy == 0))
    
    for (year in 1:(burnInYears+treatmentYears+afterYears)) {
      
      if (nrow(subset(inds, sex == 1 & yy == 0)) == 0) {
        yearResults <- data.frame(K=K,Myy=Myy,Fyy=Fyy,YYSurvival=survival,SuppressionLevel=suppression,Eliminated=1,Years=year-burnInYears,MinFemales=0)
        results <- rbind(results, yearResults)
        Population <- append(Population, nrow(inds))
        numFxx <- append(numFxx, 0)
        eliminationYear <- year
        break
      }
      
      if (year > burnInYears && year <= (burnInYears+treatmentYears)) {
        if (suppression>0) {inds <- suppress(inds,suppression)}
        inds <- stockYY(inds,Myy,Fyy)
      }
  
      inds <- death(inds,survival)
      inds <- growth(inds)
      inds <- maturity(inds)
      inds <- birth(inds,K)
      
      Population <- append(Population, nrow(inds))
      numFxx <- append(numFxx, nrow(subset(inds, sex == 1 & yy == 0)))
    }
    
    Year <- 0:(burnInYears+treatmentYears+afterYears)
    Year <- (head(Year,(length(Population))))
    numFxx <- (head(numFxx,(length(Population))))
    if(eliminationYear == 0) {yearResults <- data.frame(K=K,Myy=Myy,Fyy=Fyy,YYSurvival=survival,SuppressionLevel=suppression,Eliminated=0,Years=NA,MinFemales=min(tail(numFxx,(treatmentYears+afterYears)))); results <- rbind(results, yearResults)}
    
    if (is.element(y, plotYears)) {
      plot(Year, Population, type='l', xaxt="none", xlab = "",ylab="",main=(ifelse(eliminationYear == 0, (paste("Run", y, "- Not Extirpated. Min. Females:", min(tail(numFxx,(treatmentYears+afterYears))))), (paste("Run", y, "-", (eliminationYear-burnInYears), "years to extirpation")))), ylim=c(0,max(Population)))
      title(ylab = "Population", mgp = c(3, 1, 0))  
      title(xlab = "Year", mgp = c(2.5, 2, 0))
      axis(1,mgp = c(4, 1, 0))
      lines(Year, numFxx, type='l', col="red")
      abline(v=burnInYears, col="black",lty=2,lwd=2)
      abline(v=burnInYears+treatmentYears, col="black",lty=2,lwd=2)
    }
  }
  results
}