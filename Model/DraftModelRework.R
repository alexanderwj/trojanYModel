#parameters----

#maturity parameters: from prespawn/spawning fish collected from the SFER 5/4/23-7/12/23
b0.F <- -17.47307
b1.F <- 0.04949
b0.M <- -8.14029
b1.M <- 0.03789

#Von Bert parameters: these are from aged scales collected by WNRD in 2019 and 2022
lInf.F <- 597.13281136
rate.F <- 0.23430774
tZero.F <- -0.09127991
sigma.F <- 50.82
lInf.M <- 340.79797831
rate.M <- 0.47953492
tZero.M <- -0.03931952
sigma.M <- 30.42

#recruitment parameters: inherent growth rate and maximum recruitment as a multiplier of k.
r <- 50
k.mult <- 1.025638

#annual mortality rate A for wild fish: calculated using weighted catch-curve analysis of 2023 electrofishing catch
Z <- 0.9339988
wildMortality <- 1-exp(-Z)

#suppression selectivity parameters
s50 <- 417.8510
s95 <- 462.6664

# Simulation begins with age-2 fish, sex is randomly chosen
startingFish <- 1000
stockedAge <- 1

# Treatment years are when YY fish are stocked and suppression is applied, if applicable
burnInYears <- 50
treatmentYears <- 100
afterYears <- 50

#component functions----

maturity <- function(inds) {
  lengths <- inds$length
  numFish <- length(lengths)
  probMature <- ifelse(inds$sex==1,(exp(b0.F+(b1.F*lengths))/(1+exp(b0.F+(b1.F*lengths)))),(exp(b0.M+(b1.M*lengths))/(1+exp(b0.M+(b1.M*lengths)))))
  inds$mature<-rbinom(numFish, 1, probMature)
  inds
}

growth <- function(inds) {
  inds$age <- inds$age+1
  numFish <- length(inds$age)
  inds$length <- ifelse(inds$sex==0, lInf.M*(1-exp(-rate.M*(inds$age-tZero.M)))+rnorm(numFish, 0, sigma.M), lInf.F*(1-exp(-rate.F*(inds$age-tZero.F)))+rnorm(numFish, 0, sigma.F))
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
  spawners <- totalPairs
  newFish <- ifelse(spawners<79.63529,50*spawners,50*79.63529)*exp(rnorm(1,0,0.4))
  if (newFish == 0) {
    return(inds)
  }
  
  percMyy <- ifelse((matureMyy+matureMxy)==0,0,matureMyy/(matureMyy+matureMxy))
  percFyy <- ifelse((matureFyy+matureFxx)==0,0,matureFyy/(matureFyy+matureFxx))
  pBothyy <- percMyy*percFyy
  pWild <- (1-percFyy)*(1-percMyy)
  pFxxMyy <- (1-percFyy)*percMyy
  pFyyMxy <- (1-percMyy)*percFyy
  pGen <- runif(totalPairs)
  
  pairsBothyy <- sum(pGen<pBothyy)
  pairsWild <- sum(pGen>pBothyy & pGen<pBothyy+pWild)
  pairsFxxMyy <- sum(pGen>pBothyy+pWild & pGen<pBothyy+pWild+pFxxMyy)
  pairsFyyMxy <- sum(pGen>pBothyy+pWild+pFxxMyy & pGen<=1)
  
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
  if (Myy > 0 & stockedAge < 2) {new_Myy <- data.frame(age=rep(stockedAge, Myy), sex=0, length=0, mature=0, dead=0, yy=1, stocked=1); inds <- rbind(inds, new_Myy)}
  if (Fyy > 0 & stockedAge < 2) {new_Fyy <- data.frame(age=rep(stockedAge, Fyy), sex=1, length=0, mature=0, dead=0, yy=1, stocked=1); inds <- rbind(inds, new_Fyy)}
  if (Myy > 0 & stockedAge >= 2) {new_Myy <- data.frame(age=rep(stockedAge, Myy), sex=0, length=(lInf.M*(1-exp(-rate.M*(stockedAge-tZero.M)))+rnorm(Myy, 0, sigma.M)), mature=(1/(1+exp(log(19)*((lInf.M*(1-exp(-rate.M*(stockedAge-tZero.M))))-lm50.M)/(lm50.M-lm95.M)))), dead=0, yy=1, stocked=1); inds <- rbind(inds, new_Myy)}
  if (Fyy > 0 & stockedAge >= 2) {new_Fyy <- data.frame(age=rep(stockedAge, Fyy), sex=1, length=(lInf.F*(1-exp(-rate.F*(stockedAge-tZero.F)))+rnorm(Fyy, 0, sigma.F)), mature=(1/(1+exp(log(19)*((lInf.F*(1-exp(-rate.F*(stockedAge-tZero.F))))-lm50.F)/(lm50.F-lm95.F)))), dead=0, yy=1, stocked=1); inds <- rbind(inds, new_Fyy)}
  inds
}

suppress <-function(inds,rate) {
  lengths <- inds$length
  numFish <- length(lengths)
  suppressProb <- ifelse(lengths>240,(((1/(1+exp(log(19)*(lengths-s50)/(s50-s95))))*(1-0.1006185))+0.1006185),0)
  suppressProb <- suppressProb*rate
  inds$dead <- ifelse(inds$stocked==1, 0, rbinom(numFish, 1, suppressProb))
  inds <- subset(inds, dead==0)
  inds
}

immigration <- function(inds,num,size) {
  if (size>lInf.F) {
    vecM <- rexp(round(num/2),rate=Z)
    vecF <- rexp(round(num/2),rate=Z)
    lengthsM <- lInf.M*(1-exp(-rate.M*(vecM-tZero.M)))
    lengthsF <- lInf.F*(1-exp(-rate.F*(vecF-tZero.F)))
    newFish_M <- data.frame(age=vecM, sex=rep(0,round(num/2)), length=lengthsM, mature=rbinom(round(num/2),1,exp(b0.M+(b1.M*lengthsM))/(1+exp(b0.M+(b1.M*lengthsM)))), dead=0, yy=0, stocked=0)
    newFish_F <- data.frame(age=vecF, sex=rep(1,round(num/2)), length=lengthsF, mature=rbinom(round(num/2),1,exp(b0.F+(b1.F*lengthsF))/(1+exp(b0.F+(b1.F*lengthsF)))), dead=0, yy=0, stocked=0)
  } else {
    ifelse(size>lInf.M,ageM <- 100,ageM <- tZero.M-(log(1-size/lInf.M)/rate.M))
    ageF <- tZero.F-(log(1-size/lInf.F)/rate.F)
    mortM <- 10000-(10000*exp(-Z*ageM))
    mortF <- 10000-(10000*exp(-Z*ageF))
    numM <- round(num*(mortM/(mortM+mortF)),0)
    numF <- round(num*(mortF/(mortM+mortF)),0)
    
    vecM <- c()
    vecF <- c()
    while(length(vecM)<numM) {
      tvec <- rexp(numM*10,rate=Z)
      tvec <- tvec[tvec<=ageM]
      vecM <- append(vecM,tvec)
    }
    while(length(vecF)<numF) {
      tvec <- rexp(numF*10,rate=Z)
      tvec <- tvec[tvec<=ageF]
      vecF <- append(vecF,tvec)
    }
    
    lengthsM <- lInf.M*(1-exp(-rate.M*(vecM[0:numM]-tZero.M)))
    lengthsF <- lInf.F*(1-exp(-rate.F*(vecF[0:numF]-tZero.F)))
    
    newFish_M <- data.frame(age=vecM[0:numM], sex=rep(0,numM), length=lengthsM, mature=rbinom(round(num/2),1,exp(b0.M+(b1.M*lengthsM))/(1+exp(b0.M+(b1.M*lengthsM)))), dead=0, yy=0, stocked=0)
    newFish_F <- data.frame(age=vecF[0:numF], sex=rep(1,numF), length=lengthsF, mature=rbinom(round(num/2),1,exp(b0.F+(b1.F*lengthsF))/(1+exp(b0.F+(b1.F*lengthsF)))), dead=0, yy=0, stocked=0)
  }

  inds <- rbind(inds, newFish_M);
  inds <- rbind(inds, newFish_F);
  inds
}

emigration <- function(inds,num,size) {
  fish <- inds[sample(nrow(inds),num),]
  fish <- subset(fish,length<size)
  
  while (nrow(fish) < num) {
    fishes <- inds[sample(nrow(inds),num),]
    fishes <- subset(fishes,length<size)
    fish <- unique(rbind(fish,fishes))
  }
  
  fish <- head(fish,num)
  inds <- inds[!(row.names(inds) %in% row.names(fish)),]
  inds
}

#simulation function----

simulate <- function(K,Myy,Fyy,survival,movers,suppression,simulations,plots) {
  plotYears <- sample.int(simulations, plots)
  results <- data.frame(matrix(ncol=10,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "Eliminated", "Years", "MinFemales","EndPop", "MovingFish"))))
  
  for (y in 1:simulations) {
    inds <- data.frame(age=rep(3, startingFish), sex=rbinom(startingFish,1,0.5), length=0, mature=0, dead=0, yy=0, stocked=0)
    eliminationYear <- 0
    Population <- c(startingFish)
    numFxx <- nrow(subset(inds, sex == 1 & yy == 0))
    
    for (year in 1:(burnInYears+treatmentYears+afterYears)) {
      
      if (nrow(subset(inds, sex == 1 & yy == 0)) == 0) {
        yearResults <- data.frame(K=K,Myy=Myy,Fyy=Fyy,YYSurvival=survival,SuppressionLevel=suppression,Eliminated=1,Years=year-burnInYears,MinFemales=0,EndPop=NA, MovingFish=movers)
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
      
      if (movers>0) {
        inds <- immigration(inds,movers,999)
        inds <- emigration(inds,movers,999)
      }
      
      Population <- append(Population, nrow(inds))
      numFxx <- append(numFxx, nrow(subset(inds, sex == 1 & yy == 0)))
    }
    
    Year <- 0:(burnInYears+treatmentYears+afterYears)
    Year <- (head(Year,(length(Population))))
    numFxx <- (head(numFxx,(length(Population))))
    if(eliminationYear == 0) {yearResults <- data.frame(K=K,Myy=Myy,Fyy=Fyy,YYSurvival=survival,SuppressionLevel=suppression,Eliminated=0,Years=NA,MinFemales=min(tail(numFxx,(treatmentYears+afterYears))),EndPop=Population[burnInYears+treatmentYears],MovingFish=movers); results <- rbind(results, yearResults)}

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
  inds
}