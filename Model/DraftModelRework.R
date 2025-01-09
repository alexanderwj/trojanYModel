#parameters----

#maturity parameters: from prespawn/spawning fish collected from the SFER 5/4/23-7/12/23
b0.F <- -17.767447
b1.F <- 0.049971
b0.M <- -14.58834
b1.M <- 0.06888

#Von Bert parameters: these are from aged scales collected by WNRD in 2019 and 2022
lInf.F <- 598.02281
lInf.F.se <- 36.85475
rate.F <- 0.23366
rate.F.se <- 0.03481
tZero.F <- -0.08817
tZero.F.se <- 0.16027
sigma.F <- 50.71
lInf.M <- 339.84491
lInf.M.se <- 27.95761
rate.M <- 0.48474
rate.M.se <- 0.13155
tZero.M <- -0.03253
tZero.M.se <- 0.20391
sigma.M <- 30.68

#length mean and sd for ages 0-20: calculated using predictNLS()
meanLen.M <- c(127.5151,205.9272,254.0768,284.1627,303.2775,315.6077,323.6663,328.9904,332.5377,334.9159,336.5173,337.5985,338.3295,338.8238,339.1579,339.3836,339.5357,339.6381,339.7069,339.7530,339.7838,339.8044,339.8181,339.8272,339.8332)
sdLen.M <- c(32.46580,32.92120,33.03912,33.08162,33.76177,34.96708,36.33854,37.61338,38.67365,39.49747,40.10942,40.54986,40.85959,41.07357,41.21937,41.31763,41.38326,41.42677,41.45545,41.47425,41.48653,41.49451,41.49969,41.50303,41.50519)
meanLen.F <- c(131.7307,227.8030,303.7316,363.8131,411.4131,449.1706,479.1569,503.0002,521.9812,537.1092,549.1800,558.8222,566.5327,572.7050,577.6510,581.6182,584.8032,587.3625,589.4209,591.0777,592.4123,593.4880,594.3558,595.0562,595.6219,596.0790,596.4486,596.7475,596.9894,597.1852,597.3438,597.4722,597.5763,597.6607,597.7291)
sdLen.F <- c(51.81010,51.29143,51.31864,51.32377,51.30150,51.36909,51.61522,52.06812,52.70691,53.48364,54.34240,55.23149,56.10934,56.94592,57.72184,58.42640,59.05546,59.60953,60.09212,60.50858,60.86519,61.16854,61.42513,61.64114,61.82220,61.97341,62.09928,62.20377,62.29027,62.36173,62.42063,62.46910,62.50891,62.54157,62.56831)

#recruitment parameters: productivity and critical spawner level (crit.lev pre scaled to produce carrying capacity of 3142 age-1+ fish)
#prod <- 25; crit.lev <- 137.5
prod <- 41.6666666; crit.lev <- 82.5
#prod <- 125; crit.lev <- 27.5

#annual mortality rate A for wild fish: calculated using Chapman-Robsin method on electrofishing catch
wildMortality.M <- .5952374
wildMortality.F <- .4494907
Z.M <- -log(1-wildMortality.M)
Z.F <- -log(1-wildMortality.F)

#suppression selectivity parameters
g0 <- -6.41966
g1 <- 0.01769

# Simulation begins with age-1 fish, sex is randomly chosen
startingFish <- 1000

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
  inds$length <- ifelse(inds$sex==0,qnorm(inds$lenPct,mean=meanLen.M[inds$age],sd=sdLen.M[inds$age]),qnorm(inds$lenPct,mean=meanLen.F[inds$age],sd=sdLen.F[inds$age]))
  inds
}

death <- function(inds, survival) {
  numFish <- length(inds$age)
  inds$dead <- ifelse(inds$yy == 0, ifelse(inds$sex==0,rbinom(numFish, 1, wildMortality.M),rbinom(numFish, 1, wildMortality.F)), ifelse(inds$sex==0,rbinom(numFish, 1, 1-(1-wildMortality.M)*survival),rbinom(numFish, 1, 1-(1-wildMortality.F)*survival)))
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
  newFish <- ifelse(spawners<crit.lev,prod*spawners,prod*crit.lev)*exp(rnorm(1,0,0.4))
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
    newFish_Bothyy <- data.frame(age=rep(0, recruitsBothyy), sex=0, length=0, mature=0, dead=0, yy=1, stocked=0, lenPct=runif(recruitsBothyy))
    inds <- rbind(inds, newFish_Bothyy);
  }
  
  if(recruitsFxxMyy > 0) {
    newFish_FxxMyy <- data.frame(age=rep(0, recruitsFxxMyy), sex=0, length=0, mature=0, dead=0, yy=0, stocked=0, lenPct=runif(recruitsFxxMyy))
    inds <- rbind(inds, newFish_FxxMyy);
  }
  
  if(recruitsFyyMxy > 0) {
    newFish_FyyMxy <- data.frame(age=rep(0, recruitsFyyMxy), sex=0, length=0, mature=0, dead=0, yy=rbinom(recruitsFyyMxy, 1, 0.5), stocked=0, lenPct=runif(recruitsFyyMxy))
    inds <- rbind(inds, newFish_FyyMxy);
  }
  
  if(recruitsWild > 0) {
    recM <- rbinom(1,recruitsWild,0.5)
    recF <- recruitsWild-recM
    if (recM>0) {newFish_WildM <- data.frame(age=rep(0, recM), sex=0, length=0, mature=0, dead=0, yy=0, stocked=0, lenPct=runif(recM)); inds <- rbind(inds, newFish_WildM)}
    if (recF>0) {newFish_WildF <- data.frame(age=rep(0, recF), sex=1, length=0, mature=0, dead=0, yy=0, stocked=0, lenPct=runif(recF)); inds <- rbind(inds, newFish_WildF)}
  }
  inds
}

stockYY <- function(inds,Myy,Fyy,sAge){
  if (Myy>0) {new_Myy <- data.frame(age=rep(sAge, Myy), sex=0, length=0, mature=0, dead=0, yy=1, stocked=1, lenPct=runif(Myy)); inds <- rbind(inds, new_Myy)}
  if (Fyy>0) {new_Fyy <- data.frame(age=rep(sAge, Fyy), sex=1, length=0, mature=0, dead=0, yy=1, stocked=1, lenPct=runif(Fyy)); inds <- rbind(inds, new_Fyy)}
  inds
}

suppress <-function(inds,rate) {
  lengths <- inds$length
  numFish <- length(lengths)
  suppressProb <- (exp(g0+(g1*lengths))/(1+exp(g0+(g1*lengths))))
  suppressProb <- suppressProb*rate
  inds$dead <- ifelse(inds$stocked==1, 0, rbinom(numFish, 1, suppressProb))
  inds <- subset(inds, dead==0)
  inds
}

immigration <- function(inds,num,size) {
  if (size>lInf.F) {
    vecM <- rexp(round(num/2),rate=Z.M)
    vecF <- rexp(round(num/2),rate=Z.F)
    lengthsM <- lInf.M*(1-exp(-rate.M*(vecM-tZero.M)))
    lengthsF <- lInf.F*(1-exp(-rate.F*(vecF-tZero.F)))
    newFish_M <- data.frame(age=vecM, sex=rep(0,round(num/2)), length=lengthsM, mature=rbinom(round(num/2),1,exp(b0.M+(b1.M*lengthsM))/(1+exp(b0.M+(b1.M*lengthsM)))), dead=0, yy=0, stocked=0, lenPct=runif(round(num/2)))
    newFish_F <- data.frame(age=vecF, sex=rep(1,round(num/2)), length=lengthsF, mature=rbinom(round(num/2),1,exp(b0.F+(b1.F*lengthsF))/(1+exp(b0.F+(b1.F*lengthsF)))), dead=0, yy=0, stocked=0, lenPct=runif(round(num/2)))
  } else {
    ifelse(size>lInf.M,ageM <- 100,ageM <- tZero.M-(log(1-size/lInf.M)/rate.M))
    ageF <- tZero.F-(log(1-size/lInf.F)/rate.F)
    mortM <- 10000-(10000*exp(-Z.M*ageM))
    mortF <- 10000-(10000*exp(-Z.F*ageF))
    numM <- round(num*(mortM/(mortM+mortF)),0)
    numF <- round(num*(mortF/(mortM+mortF)),0)
    
    vecM <- c()
    vecF <- c()
    while(length(vecM)<numM) {
      tvec <- rexp(numM*10,rate=Z.M)
      tvec <- tvec[tvec<=ageM]
      vecM <- append(vecM,tvec)
    }
    while(length(vecF)<numF) {
      tvec <- rexp(numF*10,rate=Z.F)
      tvec <- tvec[tvec<=ageF]
      vecF <- append(vecF,tvec)
    }
    
    lengthsM <- qnorm(inds$lenPct,mean=meanLen.M[vecM[0:numM]],sd=sdLen.M[vecM[0:numM]])
    lengthsF <- qnorm(inds$lenPct,mean=meanLen.F[vecF[0:numF]],sd=sdLen.F[vecF[0:numF]])
    
    newFish_M <- data.frame(age=vecM[0:numM], sex=rep(0,numM), length=lengthsM, mature=rbinom(round(num/2),1,exp(b0.M+(b1.M*lengthsM))/(1+exp(b0.M+(b1.M*lengthsM)))), dead=0, yy=0, stocked=0, lenPct=runif(round(num/2)))
    newFish_F <- data.frame(age=vecF[0:numF], sex=rep(1,numF), length=lengthsF, mature=rbinom(round(num/2),1,exp(b0.F+(b1.F*lengthsF))/(1+exp(b0.F+(b1.F*lengthsF)))), dead=0, yy=0, stocked=0, lenPct=runif(round(num/2)))
  }

  inds <- rbind(inds, newFish_M);
  inds <- rbind(inds, newFish_F);
  inds
}

emigration <- function(inds,num,size) {
  if (num <= nrow(inds)) {
    fish <- inds[sample(nrow(inds),num),]
    fish <- subset(fish,length<size)
    while (nrow(fish) < num) {
      fishes <- inds[sample(nrow(inds),num),]
      fishes <- subset(fishes,length<size)
      fish <- unique(rbind(fish,fishes))
    }
  }
  else {fish <- subset(inds,length<size)}
  
  inds <- inds[!(row.names(inds) %in% row.names(fish)),]
  inds
}

#simulation function----

simulate <- function(K,Myy,Fyy,survival,movers,suppression,stockAge,simulations,plots) {
  plotYears <- sample.int(simulations, plots)
  results <- data.frame(matrix(ncol=11,nrow=0, dimnames=list(NULL, c("K", "Myy", "Fyy", "YYSurvival", "SuppressionLevel", "StockedAge", "Eliminated", "Years", "MinFemales","EndPop", "MovingFish"))))
  
  for (y in 1:simulations) {
    inds <- data.frame(age=rep(1, startingFish), sex=rbinom(startingFish,1,0.5), length=0, mature=0, dead=0, yy=0, stocked=0,lenPct=runif(startingFish))
    eliminationYear <- 0
    Population <- c(startingFish)
    numFxx <- nrow(subset(inds, sex == 1 & yy == 0))
    
    for (year in 1:(burnInYears+treatmentYears+afterYears)) {

      if (nrow(subset(inds, sex == 1 & yy == 0)) == 0) {
        yearResults <- data.frame(K=K,Myy=Myy,Fyy=Fyy,YYSurvival=survival,SuppressionLevel=suppression,StockedAge=stockAge,Eliminated=1,Years=year-burnInYears,MinFemales=0,EndPop=NA,MovingFish=movers)
        results <- rbind(results, yearResults)
        eliminationYear <- year
        break
      }
      
      if (year > burnInYears && year <= (burnInYears+treatmentYears) && (nrow(subset(inds, sex == 1 & yy == 0)) != 0)) {
        if (suppression>0) {inds <- suppress(inds,suppression)}
        if (Myy+Fyy>0) {inds <- stockYY(inds,Myy,Fyy,stockAge)}
      }
      
      inds <- death(inds,survival)
      inds <- growth(inds)
      inds <- maturity(inds)
      inds <- birth(inds,K)
      
      if (movers>0) {
        inds <- immigration(inds,movers,999)
        inds <- emigration(inds,movers,999)
      }
      
      #Biomass Calculations, Not Used Anymore
      # if (year==(burnInYears+treatmentYears)) {
      #   endBM <- sum(((-6.0913+0.9387*inds$length^2.842)*0.00003)/1000,na.rm=T)
      # }
      
      Population <- append(Population, nrow(subset(inds,length>100)))
      numFxx <- append(numFxx, nrow(subset(inds, sex == 1 & yy == 0 & length>100)))
    }
    
    Year <- 0:(burnInYears+treatmentYears+afterYears)
    Year <- (head(Year,(length(Population))))
    numFxx <- (head(numFxx,(length(Population))))
    if(eliminationYear == 0) {yearResults <- data.frame(K=K,Myy=Myy,Fyy=Fyy,YYSurvival=survival,SuppressionLevel=suppression,StockedAge=stockAge,Eliminated=0,Years=NA,MinFemales=min(tail(numFxx,(treatmentYears+afterYears))),EndPop=Population[burnInYears+treatmentYears+1],MovingFish=movers); results <- rbind(results, yearResults)}

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