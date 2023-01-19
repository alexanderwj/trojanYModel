library('omics')

#parameters----
numFish <- 500
burnInYears <- 20
stockYears <- 50
afterYears <- 20
numFyy <- 1000
numMyy <- 0
yyFitness <- 1
suppressionLevel <- 0 
numSimulations <- 5
numSets <- 1
numPlots <- 5
k <- 10000

maturity <- function(inds) {
  if (is.matrix(inds)==FALSE) {return (inds)}
  for (i in 1:dim(inds)[1]) {
    if (inds[i,4] == 0) {
      if (inds[i,1] == 3) {
        inds[i,4] <- rbinom(1,1,0.4)
      }
      if (inds[i,1] > 3) {
        inds[i,4] <- rbinom(1,1,0.95)
      }
    }
  }
  return(inds)
}

growth <- function(inds) {
  if (is.matrix(inds)==FALSE) {return (inds)}
  if(nrow(inds) < 1) {return ("No Fish:Growth")}
  inds[,1] <- inds[,1] + 1
  for (i in 1:dim(inds)[1]) {
    if (inds[i,1] == 1) {
      inds[i,2] <- (680.7*(1-exp(-.1587*(inds[i,1]-.3714)))) + rnorm(1,0,14)
    }
    else {
      inds[i,2] <- (680.7*(1-exp(-.1587*(inds[i,1]-.3714)))) + rnorm(1,0,40)
    }
  }
  return(inds)
}

death <- function(inds, yyFit) {
  if (is.matrix(inds)==FALSE) {
    if (inds[6] == 1) {
      dead <- rbinom(1,1,yyFit)
    }
    if (inds[6] == 0) {
      dead <- rbinom(1,1,0.5701418464)
    }
    if (dead == 1) {
      return ("No Fish:Death")
    }
    return (inds)
  }
  
  if (nrow(inds) < 1) {return ("No Fish:Death")}
  for (i in 1:dim(inds)[1]) {
    if (inds[i,6] == 0) {
      inds[i,3] <- rbinom(1,1,0.5701418464)
    }
    if (inds[i,6] == 1) {
      inds[i,3] <- rbinom(1,1,yyFit)
    }
  }
  inds <- inds[inds[,3] == 0,]
  return(inds)
}

birth <- function(inds, k) {
  if (is.matrix(inds)==FALSE) {return (inds)}
  if (nrow(inds) < 1) {return ("No Fish:Birth")}
  mature <- ((inds[inds[,4]==1,]))
  if (is.matrix(mature)==FALSE || nrow(mature) < 1) {return (inds)}
  
  matureFxx <- 0
  matureMxy <- 0
  matureFyy <- 0
  matureMyy <- 0
  
  matureF <- (((mature[mature[,5]==1,])))
  matureM <- (((mature[mature[,5]==0,])))

  if((is.null((dim(matureM)))) || (is.null((dim(matureF))))) {
    if(is.null((dim(matureF)))) {
      if(matureF[6]==1) {
        matureFyy <- 1
        matureFxx <- 0
      }
      else {
        matureFyy <- 0
        matureFxx <- 1
      }
      matureF <- 1
    }
    if(is.null((dim(matureM)))) {
      if(matureM[6]==1) {
        matureMyy <- 1
        matureMxy <- 0
      }
      else {
        matureMyy <- 0
        matureMxy <- 1
      }
      matureM <- 1
    }
  }
  else {
    matureFyy <- (nrow(((matureF[matureF[,6]==1,]))))
    matureMyy <- (nrow(((matureM[matureM[,6]==1,]))))
    matureFxx <- (nrow(((matureF[matureF[,6]==0,]))))
    matureMxy <- (nrow(((matureM[matureM[,6]==0,]))))
    if (is.null(matureFxx)) {matureFxx <- 0}
    if (is.null(matureMxy)) {matureMxy <- 0}
    if (is.null(matureMyy)) {matureMyy <- 0}
    if (is.null(matureFyy)) {matureFyy <- 0}
  }
  
  #print(c((matureFxx+matureFyy), (matureMxy+matureMyy)))
  adults <- matureFyy+matureFxx+matureMyy+matureMxy
  totalPairs <- (min(c((matureFxx+matureFyy), (matureMxy+matureMyy))))
  if(totalPairs == 0) {return(inds)}
  #print(totalPairs)
  
  percMyy <- matureMyy/(matureMyy+matureMxy)
  percFyy <- matureFyy/(matureFyy+matureFxx)
  percBothyy <- percMyy*percFyy
  
  spawners <- totalPairs*2
  newFish = (((sqrt(k)+1)*spawners)/(1+spawners/(sqrt(k))))+rnorm(1,0,k/5)
  #newFish = (18.004*spawners*(exp(-.005*spawners)))+rnorm(1,0,200/4)
  if (newFish <= 0) {
    return(inds)
  }
  newFish <- round(newFish,0)
  
  pairsBothyy <- rbinom(1,totalPairs,percBothyy)
  pairsFxxMyy <- rbinom(1,(totalPairs-pairsBothyy),percMyy)
  pairsFyyMxy <- rbinom(1,(totalPairs-pairsBothyy-pairsFxxMyy),percFyy)
  pairsWild <- (totalPairs-pairsBothyy-pairsFxxMyy-pairsFyyMxy)

  recruitsBothyy <- round((newFish*(pairsBothyy/totalPairs)), 0)
  recruitsFxxMyy <- round((newFish*(pairsFxxMyy/totalPairs)), 0)
  recruitsFyyMxy <- round((newFish*(pairsFyyMxy/totalPairs)), 0)
  recruitsWild <- (newFish-recruitsBothyy-recruitsFxxMyy-recruitsFyyMxy)
  
  # print("Mature: ")
  # print(c(matureMxy, matureFxx, matureMyy, matureFyy))
  # cat("Pairs: ", totalPairs, "\n")
  # print(c(pairsWild, pairsFxxMyy, pairsFyyMxy, pairsBothyy))
  # print("Recruits: ")
  # print(c(recruitsBothyy, recruitsFxxMyy, recruitsFyyMxy, recruitsWild))

  if(recruitsBothyy > 0) {
    newFish_Bothyy <- array(data = 0, dim = c(recruitsBothyy, 7))
    for (i in 1:recruitsBothyy) {
      newFish_Bothyy[i,5] <- 0
      newFish_Bothyy[i,6] <- 1
    }
    inds <- rbind(inds, newFish_Bothyy);
  }
  
  if(recruitsFxxMyy > 0) {
    newFish_FxxMyy <- array(data = 0, dim = c(recruitsFxxMyy, 7))
    for (i in 1:recruitsFxxMyy) {
      newFish_FxxMyy[i,5] <- 0
      newFish_FxxMyy[i,6] <- 0
    }
    inds <- rbind(inds, newFish_FxxMyy);
  }
  
  if(recruitsFyyMxy > 0) {
    newFish_FyyMxy <- array(data = 0, dim = c(recruitsFyyMxy, 7))
    newFish_FyyMxy[,6] <- rbinom(recruitsFyyMxy,1,0.5)
    newFish_FyyMxy[,5] <- 0
    inds <- rbind(inds, newFish_FyyMxy);
  }
  
  if(recruitsWild > 0) {
    newFish_Wild<- array(data = 0, dim = c(recruitsWild, 7))
    for (i in 1:recruitsWild) {
      newFish_Wild[i,5] <- rbinom(1,1,0.5)
      newFish_Wild[i,6] <- 0
    }
    inds <- rbind(inds, newFish_Wild);
  }
  
  return(inds)
}

stockMyy <- function(inds, numFish){
  if (numFish == 0) {return(inds)}
  new_inds <- array(data = 0, dim = c(numFish, 7))
  new_inds[,1] <- 1
  new_inds[,5] <- 0
  new_inds[,6] <- 1
  inds <- rbind(inds, new_inds);
  return(inds)
}

stockFyy <- function(inds, numFish){
  if (numFish == 0) {return(inds)}
  new_inds <- array(data = 0, dim = c(numFish, 7))
  new_inds[,1] <- 1
  new_inds[,5] <- 1
  new_inds[,6] <- 1
  inds <- rbind(inds, new_inds);
  return(inds)
}

suppress <-function(inds, level) {
  if (level == 0) {return(inds)}
  if (is.null(numFish)) {return(inds)}
  if (is.matrix(inds)==FALSE) {return(inds)}
  numFish <- nrow(inds)
  lengthZero <- nrow(inds[inds[,2]==0,])
  currentRemoved <- 0
  numsRemoved <- c()
  for (g in 1:50000) {
    number <- sample(1:numFish, 1)
    if (is.element(number,numsRemoved)) {next}
    length <- (inds[number,2])
    if (length == 0) {next}
    suppressProb <- 0.00008526*1.025^(.8753*(length))+.01753
    if (suppressProb > 1) {suppressProb <- 1}
    suppressed <- rbinom(1,1,suppressProb)
    if (suppressed == 1) {
      inds[number,7] <- 1
      numsRemoved <- append(numsRemoved,number)
      currentRemoved <- currentRemoved + 1
    }
    if (currentRemoved == level || currentRemoved == (numFish-lengthZero)) {currentRemoved <- 0; numsRemoved <- c(); break}
  }
  #notRemoved <- inds[inds[,7]==0,]
  #if (is.matrix(notRemoved) == FALSE) {currentRemoved <- 0; numsRemoved <- c(); return(inds)}
  #if (nrow(notRemoved[notRemoved[,2]==0,]) == nrow(notRemoved)) {currentRemoved <- 0; numsRemoved <- c(); return(inds[inds[,7] == 0,])}
  inds <- inds[inds[,7] == 0,]
  return(inds)
}

#totalYearsToElim <- matrix(0, ncol=31, nrow=31)
#yyFitnessByYear <- seq(0.4, 1, 0.02)
#FyyRatioByYear <- seq(0, 1.01, 0.033)

for (j in 1:numSets) {
  yearsToElim <- c()
  #numFyyByYear <- c()
  #numMyyByYear <- c()
  #yyFitnessByYear <- seq(0.4, 1, 0.02)
  #FyyRatioByYear <- c()
  #suppressionLevelByYear <- c()
  #yyFitness <- 0.38
  #numFyy <- -33
  plotYears <- sample.int(numSimulations, numPlots)
  print(plotYears)
  yyFit <- (1-((1-0.5701418464)*yyFitness))
  #row <- 0
  
  start_time <- Sys.time()
  for (y in  1:numSimulations) {
    cat((numSimulations*(j-1)+y),'')
    inds  <- array(data = 0, dim = c(numFish, 7))
    colnames(inds) <- c("age", "length", "death", "mature", "sex", "yy", "suppressed")
    rownames(inds) <- NULL
    for (i in 1:numFish) {
      inds[i,5] <- rbinom(1,1,0.5)
      inds[i,1] <- 3
      inds[i,4] <- 1
    }
    currentYear <- -1
    breakCond <- 0
    Population <- c(numFish)
    Females <- (inds[inds[,5] == 1,])
    numFxx <- c((nrow(Females[Females[,6] == 0,])))
    
    # row <- row+1
    # numFyy <- numFyy+33
    # if (y%%31 == 1) {
    #   yyFitness <- yyFitness + 0.02
    #   numFyy <- 0
    #   row <- 1
    # }
    #suppressionLevel <- 500
    #suppressionLevel <- sample(1:500,1)
    #numFyy <- sample(1:1000, 1)
    #numMyy <- 1000-numFyy
    #FyyRatio <- numFyy/1000
    
    for (year in 1:(stockYears+burnInYears+afterYears)) {
      #cat("Year", year,"\n")
      
      inds <- growth(inds=inds)
      #print("Growth")
      if (is.character(inds)) {Population <- append(Population, 0); numFxx <- append(numFxx, 0); breakCond <- 1; break}
      inds <- maturity(inds=inds)
      #print("Maturity")
      if (is.character(inds)) {Population <- append(Population, 0); numFxx <- append(numFxx, 0); breakCond <- 1; break}
      inds <- birth(inds=inds, k=k)
      #print("Birth")
      if (is.character(inds)) {Population <- append(Population, 0); numFxx <- append(numFxx, 0); breakCond <- 1; break}
      inds <- death(inds=inds, yyFit=yyFit)
      #print("Death")
      if (is.character(inds)) {Population <- append(Population, 0); numFxx <- append(numFxx, 0); breakCond <- 1; break}
      
      if ((year>burnInYears) && year<=(burnInYears+stockYears)) {
        inds <- stockMyy(inds,numMyy)
        inds <- stockFyy(inds,numFyy)
        inds <- suppress(inds,suppressionLevel)
      }
      #print("Stock/Suppress")
      if (is.matrix(inds) == FALSE) {
        if (inds[5] == 0) {Population <- append(Population, 0); numFxx <- append(numFxx, 0); breakCond <- 1; break}
        if (inds[5] == 1) {
          if (inds[6] == 0) {
            Population <- append(Population, 0); numFxx <- append(numFxx, 1); next
          }
          else {Population <- append(Population, 0); numFxx <- append(numFxx, 0); breakCond <- 1; break}
        }
      }
      
      Population <- append(Population, nrow(inds))
      Females <- (inds[inds[,5] == 1,])
      if (is.matrix(Females) == FALSE) {
        if (Females[6] == 0) {numFxx <- append(numFxx, 1); next}
        else if (Females[6] == 1) {numFxx <- append(numFxx, 0); breakCond <- 1; break}
      }
      if (nrow(Females) == 0) {numFxx <- append(numFxx, 0); breakCond <- 1; break}
      
      Fxx <- (Females[Females[,6] == 0,])
      if (is.matrix(Fxx) == FALSE) {numFxx <- append(numFxx, 1); breakCond <- 1; next}
      if (nrow(Fxx) == 0) {numFxx <- append(numFxx, 0); breakCond <- 1; break}
      else {numFxx <- append(numFxx, (nrow(Fxx))); next}
      
      numFxx <- append(numFxx, Fxx)
    }
    
    Year <- 0:(burnInYears+stockYears+afterYears)
    Year <- (head(Year,(length(Population))))
    numFxx <- (head(numFxx,(length(Population))))
    
    eliminatedYears <- (which(numFxx==0))
    if (identical(eliminatedYears, integer(0)) == FALSE) {
      eliminatedYear <- eliminatedYears[1]
      yearsToElim <- append(yearsToElim, (eliminatedYear-burnInYears-1))
      #numFyyByYear <- append(numFyyByYear, numFyy)
      #numMyyByYear <- append(numMyyByYear, numMyy)
      #yyFitnessByYear <- append(yyFitnessByYear, yyFitness)
      #FyyRatioByYear <- append(FyyRatioByYear, FyyRatio)
      #print(FyyRatio)
      #print(yyFitness)
      #print(row)
      #print(ceiling(y/31))
      #totalYearsToElim[row, ceiling(y/31)] <- totalYearsToElim[row, ceiling(y/31)]+(eliminatedYear-burnInYears-1)
      #print(totalYearsToElim)
      #suppressionLevelByYear <- append(suppressionLevelByYear, suppressionLevel)
      
      if (is.element(y, plotYears)) {
        # print(inds)
        # print(numFxx)
        # print(Population)
        plot(Year, Population, type='l', main=(paste("Run", y, "-", (eliminatedYear-burnInYears-1), "years to extirpation")), ylim=c(0,max(Population)))
        lines(Year, numFxx, type='l', col="red")
        abline(v=burnInYears, col="red")
        abline(v=burnInYears+stockYears, col="red")
      }
    }
    else {
      #yearsToElim <- append(yearsToElim, 120)
      #FyyRatioByYear <- append(FyyRatioByYear, FyyRatio)
      #print(FyyRatio)
      #print(yyFitness)
      #print(row)
      #print(ceiling(y/31))
      #totalYearsToElim[row, ceiling(y/31)] <- totalYearsToElim[row, ceiling(y/31)]+120
      #print(totalYearsToElim)
      if (is.element(y, plotYears)) {
        plot(Year, Population, type = "l", main=(paste("Run", y, "- Not Extirpated. Min. Females:", min(tail(numFxx,(stockYears+afterYears))))), ylim=c(0,max(Population)))
        lines(Year,numFxx,type="l",col="red")
        abline(v=burnInYears, col="red")
        abline(v=burnInYears+stockYears, col="red")
      }
    }
    
    #print(length(Year))
    #print(length(Population))
    #print(length(numFxx))
    #print(Females)
    #print(numFxx)
    
    remove(Population)
    remove(numFxx)
    remove(inds)
  }
  #filled.contour(x = FyyRatioByYear, y = yyFitnessByYear, z = totalYearsToElim, xlab = 'Fyy:Myy Ratio', ylab = 'yy Relative Survival', plot.axes = {contour(FyyRatioByYear,yyFitnessByYear,totalYearsToElim, add=TRUE, lwd=0.5); axis(1); axis(2)})
}
end_time <- Sys.time()

#avgTotalYearsToElim <- totalYearsToElim/5
#totalYearsToElim <- (head(yearsToElim,(length(yyFitnessByYear))))
#logYearsToElim <- log(totalYearsToElim)
#plot(FyyRatioByYear, yearsToElim, main=("Fyy:Myy Ratio vs. Years to Extirpation"))
#plot(numFyyByYear, logYearsToElim, main=("log(Num Fyy Stocked out of 2000) vs. Years to Extirpation"))
#abline((lm(yearsToElim~numFyyByYear)),col='green')

#filled.contour(x = FyyRatioByYear, y = yyFitnessByYear, z = avgTotalYearsToElim, xlab = 'Fyy:Myy Ratio', ylab = 'yy Relative Survival', plot.axes = {contour(FyyRatioByYear,yyFitnessByYear,avgTotalYearsToElim, add=TRUE, lwd=0.5); axis(1); axis(2)})
#contour(xgrid, ygrid, mtrx3d, col="black", add=TRUE)

cat("\nTotal Execution Time:", end_time-start_time, "seconds")
cat("\nExtirpation Percentage: ", (((length(yearsToElim))/numSimulations)*100), "%\n", sep='')
if (is.null(yearsToElim) == FALSE) {
  cat("Of Extirpating Runs, Mean: ", (mean(yearsToElim)), ", SD: ", (sd(yearsToElim)), ", Range: ", (range(yearsToElim)[1]), "-", (range(yearsToElim)[2]),  sep='')
}