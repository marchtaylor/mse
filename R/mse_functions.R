
# iter of FLIndex or FLIndices --------------------------------------------
iterIndicies <- function(indicies, iter=NULL){
  if(is.null(iter)) iter <- 1 # if no iter defined, use 1
  if(! class(indicies) %in% c("FLIndices", "FLIndex")) stop("object must be of class 'FLIndices' or 'FLIndex'")
  indicies <- FLIndices(indicies) # convert to FLIndices for level loop below
  for(lev in seq(length(idxErr))){
    indicies[[lev]] <- iter(indicies[[lev]], iter)
  }
  if(length(indicies)==1) indicies <- indicies[[1]]
  return(indicies)
}


# Stock assessment -------------------------------------------------------
assessStock <- function(
stock, method="XSA", 
iter=NULL, minyear=NULL, maxyear=NULL,
index = NULL, control = NULL,
fratio = "missing", fit.plusgroup = TRUE, desc = NULL,
stockReal=NULL, assessErr=NULL, alpha=1,
...){
  if(is.null(iter)) iter <- 1 # if no iter defined, use 1
  if(is.null(minyear)) minyear <- range(stock)["minyear"] # if no minyear, assume full min year extent
  if(is.null(maxyear)) maxyear <- range(stock)["maxyear"] # if no minyear, assume full max year extent
  
  # perform assessment on stockSub
  if(method == "XSA"){
    # create subset of stock and indices
    stockSub <- window(iter(stock, iter), start=minyear, end=maxyear)
    indexSub <- window(iterIndicies(index, iter), start=minyear, end=maxyear)
    # run XLS
    res <- FLXSA(stockSub, indexSub, control)
    stockSub <- stockSub + res
  }
  
  if(method == "VPA"){
    # create subset of stock and indices
    stockSub <- window(iter(stock, iter), start=minyear, end=maxyear)
    
    # take terminal F from previous years average
    DIMS <- dim(stockSub@stock.n)
    range(stockSub)
    stockSub@harvest[,ac(maxyear)] <- apply(stockSub@harvest[,ac((maxyear-3):(maxyear-1))], 1, mean)    
  
    # run VPA
    res <- FLAssess::VPA(stockSub) #, fratio = fratio, fit.plusgroup = fit.plusgroup, desc = desc)
    stockSub <- stockSub + res
  }
  
  if(method == "brute"){ # Restrepo, V.R., Legault, C.M., 1995. Approximations for solving the catch equation when it involves a Ilplus group"·. Fishery bulletin 93, 308–314.
    # create subset of stock and error
    stockSub <- window(iter(stock, iter), start=minyear, end=maxyear)
    stockRealSub <- window(iter(stockReal, iter), start=minyear, end=maxyear)
    
    # project stock.n at t+1
    Ns <- window(stockRealSub@stock.n, end=maxyear+1)
    Cs <- window(stockRealSub@catch.n, end=maxyear+1)
    DIMS <- dim(Ns)
    Zs <-  (stockRealSub@m + stockRealSub@harvest) # natural death numbers
    for(i in seq(DIMS[1])){
      if(i < (DIMS[1])){ Ns[i+1,ac(maxyear+1)] <- (Ns[i,ac(maxyear)] * exp(-Zs[i,ac(maxyear)])) }
      if(i == DIMS[1]){ Ns[i,ac(maxyear+1)] <- Ns[i,ac(maxyear+1)] + (Ns[i,ac(maxyear)] * exp(-Zs[i,ac(maxyear)])) } # terminal year
    }
    Fs <- stockSub@harvest
    Ms <- stockSub@m
    
    # add error to estimated numbers
    errorSub <- window(iter(assessErr, iter), start=minyear, end=maxyear+1)
    NaNs <- which(is.na(errorSub))
    if(length(NaNs)>0){
      tmp <- c(errorSub)
      tmp[NaNs] <- 1
      errorSub <- FLQuant(tmp, dim=dim(errorSub), dimnames=dimnames(errorSub))
    }
    Ns <- Ns * errorSub
   
    # calc Fs
    non.plus.ages <- seq(DIMS)[-c(DIMS[1]-1, DIMS[1])]
    plus.ages <- c(DIMS[1]-1, DIMS[1])
    
    # stkReal_nn@harvest[seq(ageDim)[ageDim],] / stkReal_nn@harvest[seq(ageDim)[ageDim-1],]
    Fs[non.plus.ages,]  <- log( (Cs[-DIMS[1],-DIMS[2]] / Ns[-1,-1]) * exp(-Ms[-DIMS[1],-DIMS[2]]/2) + 1 )[non.plus.ages,]
    Fs[plus.ages,]      <- log( ((Cs[DIMS[1],-DIMS[2]]/alpha + Cs[DIMS[1]-1,-DIMS[2]]) / Ns[DIMS[1],-1]) * exp(-Ms[DIMS[1]-1,-DIMS[2]]/2) + 1)
    
    # empirical correction?
    Fs <- Fs/(0.9970 + 0.0808*Ms)
    
    # Replace values in stockSub
    stockSub@stock.n <- window(Ns, end=maxyear)
    stockSub@harvest <- Fs
  }
  
  # return results to stock
  stock[,ac(minyear:maxyear),,,,iter] <- stockSub[,ac(minyear:maxyear)]
  return(stock)
}


# TAC calculation ------------------------------------------------------------
hcrSform <- function(
  stock, tac, # observed stock and tac objects
  yr_ass = NULL, # year of assessment
  iter=NULL, # iteration number
  Fmin=0.1, Fmax=0.3, # min and max F associated with min and max SSB
  SSBmin=250000, SSBmax=500000, # min and max SSB limits for min and max F
  relcatchLim=0.15 # restriction of relative allowed yearly changes
){
  
  if(is.null(iter)) iter <- 1
  if(is.null(yr_ass)) yr_ass <- range(stock)["maxyear"]
  
  # create subset of stock and tac
  stockSub <- window(iter(stock, iter), end=yr_ass)
  tacSub <- iter(tac, iter)
  
  SSB <- c(ssb(stockSub[,ac(yr_ass)]))
  
  # if SSB below rebuilding size
  if(SSB < SSBmin){
    ctrl_target <- rbind(
      f_df <- data.frame(
        year = yr_ass+2,
        rel.year = NA,
        quantity = "f",
        val = Fmin,
        max = NA,
        min = NA),
      catch2_df <- data.frame(
        year = yr_ass+2,
        rel.year = yr_ass+1,
        quantity = "catch",
        val = NA,
        max = 1 + relcatchLim,
        min = 1 - relcatchLim)
    )
  }
  # if SSB in rebuilding size
  if(SSB >= SSBmin & SSB < SSBmax){
    ctrl_target <- rbind(
      f_df <- data.frame(
        year = yr_ass+2,
        rel.year = NA,
        quantity = "f",
        val = Fmin + (SSB-SSBmin) * ((Fmax-Fmin)/(SSBmax-SSBmin)),
        max = NA,
        min = NA),
      catch2_df <- data.frame(
        year = yr_ass+2,
        rel.year = yr_ass+1,
        quantity = "catch",
        val = NA,
        max = 1 + relcatchLim,
        min = 1 - relcatchLim)
    )
  }
  # if SSB above biological target
  if(SSB >= SSBmax){
    ctrl_target <- rbind(
      f_df <- data.frame(
        year = yr_ass+2,
        rel.year = NA,
        quantity = "f",
        val = Fmax,
        max = NA,
        min = NA),
      catch2_df <- data.frame(
        year = yr_ass+2,
        rel.year = yr_ass+1,
        quantity = "catch",
        val = NA,
        max = 1 + relcatchLim,
        min = 1 - relcatchLim)
    )
  }
  
  # if no TAC defined for yr_ass+1, then maintain F from last year
  if( is.na(tacSub[,ac(yr_ass+1)]) ){
    ctrl_target <- rbind(
      ctrl_target,
      catch1_df <- data.frame(
        year = yr_ass+1,
        rel.year = yr_ass,
        quantity = "f",
        val = 1,
        max = NA,
        min = NA
      )
    )
  } else { # else, grab previosly defined TAC for yr_ass+1
    ctrl_target <- rbind(
      ctrl_target,
      catch1_df <- data.frame(
        year = yr_ass+1,
        rel.year = NA,
        quantity = "catch",
        val = c(tacSub[,ac(yr_ass+1)]),
        max = NA,
        min = NA
      )
    )
  }
  
  # Order by year
  ctrl_target <- ctrl_target[order(ctrl_target$year),]
  
  # Make the control object
  ctrl_obj <- fwdControl(ctrl_target)
  
  # determine average recruitment from last 3 years
  stockSub_sr <- as.FLSR(stockSub, model="geomean")
  mean_rec <- mean(rec(stockSub)[,ac((yr_ass-3):(yr_ass-1))])
  params(stockSub_sr)['a',] <- mean_rec
  
  # translate control object F value to TAC (catch)
  stockSub <- stf(stockSub, 2)
  stockSub <- fwd(stockSub, ctrl = ctrl_obj, sr = stockSub_sr) 
  
  # record TAC
  tacSub[,ac(yr_ass:(yr_ass+2))] <- stockSub@catch[,ac(yr_ass:(yr_ass+2))]
  
  # return results to tac
  iter(tac, iter) <- tacSub
  return(tac)
  
}

hcrConstF <- function(
  stock, tac, # observed stock and tac objects
  yr_ass = NULL, # year of assessment
  iter=NULL, # iteration number
  harvest=0.3 # constant F
){
  
  if(is.null(iter)) iter <- 1
  if(is.null(yr_ass)) yr_ass <- range(stock)["maxyear"]
  
  # create subset of stock and tac
  stockSub <- window(iter(stock, iter), end=yr_ass)
  tacSub <- iter(tac, iter)
  
  ctrl_target <- data.frame(
    year = yr_ass+2,
    rel.year = NA,
    quantity = "f",
    val = harvest,
    max = NA,
    min = NA
  )
  
  # if no TAC defined for yr_ass+1, then maintain F from last year
  if( is.na(tacSub[,ac(yr_ass+1)]) ){
    ctrl_target <- rbind(
      ctrl_target,
      catch1_df <- data.frame(
        year = yr_ass+1,
        rel.year = yr_ass,
        quantity = "f",
        val = 1,
        max = NA,
        min = NA
      )
    )
  } else { # else, grab previosly defined TAC for yr_ass+1
    ctrl_target <- rbind(
      ctrl_target,
      catch1_df <- data.frame(
        year = yr_ass+1,
        rel.year = NA,
        quantity = "catch",
        val = c(tacSub[,ac(yr_ass+1)]),
        max = NA,
        min = NA
      )
    )
  }
  
  # Order by year
  ctrl_target <- ctrl_target[order(ctrl_target$year),]
  
  # Make the control object
  ctrl_obj <- fwdControl(ctrl_target)
  
  # determine average recruitment from last 3 years
  stockSub_sr <- as.FLSR(stockSub, model="geomean")
  mean_rec <- mean(rec(stockSub)[,ac((yr_ass-3):(yr_ass-1))])
  params(stockSub_sr)['a'] <- mean_rec
  
  # translate control object F value to TAC (catch)
  stockSub <- stf(stockSub, 2)
  stockSub <- fwd(stockSub, ctrl = ctrl_obj, sr = stockSub_sr) 
  
  # record TAC
  tacSub[,ac(yr_ass:(yr_ass+2))] <- stockSub@catch[,ac(yr_ass:(yr_ass+2))]
  
  # return results to tac
  iter(tac, iter) <- tacSub
  return(tac)
  
}

calcTac <- function(stock, tac, hcrFun, ...){
  tac <- hcrFun(stock, tac, ...)
  return(tac)
}


# Advance stock (operational model) ---------------------------------
advanceStock <- function(stock, sr, tac, iter=NULL, minyear=NULL, maxyear=NULL, ...){
  if(is.null(minyear)) minyear <- range(stock)["minyear"]
  if(is.null(maxyear)) maxyear <- range(stock)["maxyear"]
  if(is.null(iter)) iter <- 1
  
  # create subset of stock, tac
  stockSub <- window(iter(stock, iter), start=minyear, end=maxyear+1)
  tacSub <- window(iter(tac, iter), start=minyear, end=maxyear+1)
  
  ctrl_target <- data.frame(
    year = maxyear+1,
    quantity = "catch",
    val = c(tacSub[,ac(maxyear+1)])
  )
  ctrl_obj <- fwdControl(ctrl_target)
  
  # project forward with real biology
  # stockSub <- stf(stockSub, 1)
  stockSub <- fwd(stockSub, ctrl = ctrl_obj, sr = sr, ...)
  
  # return results to stock
  stock[,ac(minyear:(maxyear+1)),,,,iter] <- stockSub
  return(stock)
}


# Observe catches (with possible error) -----------------------------------
observeCatch <- function(stockReal, stockObs, obsErr=NULL, iter=NULL, minyear=NULL, maxyear=NULL){
  
  if(is.null(obsErr)) obsErr <- replace(stockReal@catch.n, is.na(stockReal@catch.n), 0)*0 + 1 # If missing, then no observation error
  if(is.null(minyear)) minyear <- as.numeric(dims(obsErr)["minyear"])
  if(is.null(maxyear)) maxyear <- as.numeric(dims(obsErr)["maxyear"])
  if(is.null(iter)) iter <- 1
  
  # create subset of stocks and error
  stockRealSub <- window(iter(stockReal, iter), start=minyear, end=maxyear)
  stockObsSub <- window(iter(stockObs, iter), start=minyear, end=maxyear)
  obsErrSub <- window(iter(obsErr, iter), start=minyear, end=maxyear)
  
  # Observe real stock catches and pass to observed with error
  stockObsSub@catch.n <- stockRealSub@catch.n * obsErrSub
  stockObsSub@discards.n <- stockRealSub@discards.n * obsErrSub
  stockObsSub@landings.n <- stockRealSub@landings.n * obsErrSub
  stockObsSub@catch <- computeCatch(stockObsSub)
  stockObsSub@discards <- computeDiscards(stockObsSub)
  stockObsSub@landings <- computeLandings(stockObsSub)
  
  # return results to stock
  stockObs[,ac(minyear:maxyear),,,,iter] <- stockObsSub
  return(stockObs)
}


# Observe indices (with possible error) -----------------------------------
observeIndex <- function(stock, index, q.hat, obsErr=NULL, iter=NULL, minyear=NULL, maxyear=NULL){
  
  if(is.null(iter)) iter <- 1
  index <- FLIndices(index) # convert to FLIndices for level loop below
  if(!is.null(obsErr)){
    obsErr <- FLIndices(obsErr) # convert to FLIndices for level loop below
  }else{
    obsErr <- index
    for(lev in seq(obsErr)){
      # If missing error object, then no observation error
      obsErr[[lev]]@index <- replace(obsErr[[lev]]@index, is.na(obsErr[[lev]]@index), 0)*0 + 1 
    }
  }
  
  # for each level of index
  for(lev in seq(index)){
    # create subset of stock, index and index error
    if(is.null(minyear)) minyear <- as.numeric(dims(obsErr[[lev]])["minyear"])
    if(is.null(maxyear)) maxyear <- as.numeric(dims(obsErr[[lev]])["maxyear"])
    stockSub <- window(iter(stock, iter), start=minyear, end=maxyear)
    indexSub <- window(iter(index[[lev]], iter), start=minyear, end=maxyear)
    obsErrSub <- window(iter(obsErr[[lev]], iter), start=minyear, end=maxyear)
    
    # identify extent of index
    indexAgeMin <- ifelse(is.na(range(indexSub)["min"]), NA, range(indexSub)["min"])
    indexAgeMax <- ifelse(is.na(range(indexSub)["max"]), NA, range(indexSub)["max"])
    indexAges <- indexAgeMin:indexAgeMax
    alpha <- ifelse(is.na(range(indexSub)["startf"]), 0, range(indexSub)["startf"])  # year fraction start of survey
    beta <- ifelse(is.na(range(indexSub)["endf"]), 0.01, range(indexSub)["endf"]) # year fraction end of survey
    Zay <- c(stockSub@m[ac(indexAges),] + stockSub@harvest[ac(indexAges),]) # Population at age (a) at beginning of year (y)
    Pay <- c(stockSub@stock.n[ac(indexAges),]) # Total mortality at age
    qa <- c(q.hat[[lev]]) # catchability at age
    Pay.plus <- exp( log(Pay) - alpha*Zay + log( (1-exp(-(beta-alpha)*Zay)) / ((beta-alpha)*Zay) ) ) # Population at age (a) at projected forward in year (y) to between alpha and beta
    CPUE <- exp(log(qa)+log(Pay.plus))
    
    # record indices of CPUE
    indexSub@index <- CPUE * obsErrSub@index
    
    # record catch
    indexSub@catch.n <- indexSub@index * 1
    indexSub@effort <- FLQuant(1, dimnames = list(year=minyear:maxyear, iter=1))
    
    # return results to index
    index[[lev]][,ac(minyear:maxyear),,,,iter] <- indexSub
  }
  
  if(length(index)==1) index <- index[[1]]
  return(index)
}

