
### Setup --------------------------------------------------------------
# setwd("/home/")
results_path <- "results/msy_brute/"
dir.create(results_path)


# General settings
idxCV <- 0.38
nyear <- 10
niter <- 10
wts.nyears <- 5
fbar.nyears <- 5
loop_multi <- expand.grid(F=seq(0.1,0.7,0.1), stringsAsFactors = FALSE)
loop_multi
save(loop_multi, file="loop_multi.Rdata")
ncores <- 2 #as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS"))
unlink("output.txt")


# Required packages -------------------------------------------------------
library(FLCore) # (version: ‘2.5.20160107’)
library(FLash)  # (version: ‘2.5.2’)
library(FLAssess)  # (version: ‘2.5.20130716’)
library(FLXSA) # (version: ‘2.5.20140808’)
library(ggplotFL) # (version: ‘2.5.20160119’)
library(snow) # for multicore use (version: ‘0.4.1’)


# Required functions ------------------------------------------------------
source("R/mse_functions.R")


# Required data -----------------------------------------------------------
# stock object
data(ple4)
stkReal <- ple4

# sequence of forcast years
year_seq <- seq(range(stkReal)["maxyear"]+1, length.out = nyear)


# Read survey indicies
data("ple4.indices")
idxReal <- ple4.indices


# Conduct initial assessment ----------------------------------------------
# XSA control
xsaCtrl <- FLXSA.control(
  tol = 1e-09, maxit = 50, min.nse = 0.3, fse = 1,
  rage = 0, qage = 10, shk.n = TRUE, shk.f = TRUE,
  shk.yrs = 5, shk.ages= 5, window = 100, tsrange = 99,
  tspower = 1
)

# XSA and update stkReal
xsaReal <- FLXSA(stkReal, idxReal, xsaCtrl)
stkReal <- stkReal + xsaReal

# record survey catchability
q.hat <- xsaReal@q.hat 



# SRR ---------------------------------------------------------------------
# Assessment-assumed SRR
stkRealSr <- as.FLSR(stkReal)
range(stkRealSr)
model(stkRealSr) <- "ricker"
stkRealSr <- fmle(stkRealSr) 
plot(stkRealSr)


# Expand objects and duplicate observed objects ---------------------------
# expand years and carrying forward weight at age, selectivity (fbar), etc.
stkReal <- stf(stkReal, nyear=nyear, wts.nyears=wts.nyears, fbar.nyears=fbar.nyears)
idxReal <- window(idxReal, end=tail(year_seq,1))

# expand to niter
stkReal <- propagate(stkReal, niter)
for(i in seq(idxReal)){
  idxReal[[i]] <- propagate(idxReal[[i]], niter)
}

# duplicate observed stock and indices
stkObs <- stkReal
idxObs <- idxReal


# Create TAC object -------------------------------------------------------
TAC <- catch(stkReal) # TAC
dims(TAC)


# Noise and error generation --------------------------------------------------------
# recruitment variability
set.seed(1)
recrErr <- FLQuant(exp(sample(c(residuals(stkRealSr)), nyear*niter, replace=TRUE)), dimnames = list(year=year_seq, iter=1:niter))

# Indices error from surveys
idxErr <- idxReal
for(i in seq(idxReal)){
  idxErr[[i]]@index <- FLQuant(rep(1, prod(dim(idxReal[[i]]@index))), dimnames = dimnames(idxReal[[i]]@index))
  idxErrDim <- dim(idxReal[[i]]@index[,ac(year_seq)])
  idxErr[[i]]@index[,ac(year_seq)] <- rlnorm(prod(idxErrDim), meanlog = 0, sdlog = idxCV)
}
sd(log(c(idxErr[[1]]@index[3,ac(year_seq),,,,]))); idxCV # a check

# Assessment error (for use with assessment method = "Restrepo")
assCV <- 0.1
assErr <- FLQuant(rep(1, prod(dim(stkReal@stock.n))), dimnames = dimnames(stkReal@stock.n))
assErrDim <- dim(stkReal@stock.n[,ac(year_seq)])
assErr[,ac(year_seq)] <- rlnorm(prod(assErrDim), meanlog = 0, sdlog = assCV)
sd(log(c(assErr[3,ac(year_seq),,,,]))); assCV # a check
# plot(assErr)

# Catch numbers error
catchCV <- 0
catchErr <- FLQuant(rep(1, prod(dim(stkReal@catch.n))), dimnames = dimnames(stkReal@catch.n))
catchErrDim <- dim(stkReal@catch.n[,ac(year_seq)])
catchErr[,ac(year_seq)] <- rlnorm(prod(catchErrDim), meanlog = 0, sdlog = catchCV)
sd(log(c(catchErr[3,ac(year_seq),,,,]))); catchCV # a check
# plot(catchErr)


### MSE ---------------------------------------------------------------------

for (rn in seq(nrow(loop_multi))){
 
  # make duplicates of objects
  stkReal_dup  <- stkReal
  stkObs_dup   <- stkObs
  idxObs_dup   <- idxObs
  idxReal_dup  <- idxReal
  idxErr_dup   <- idxErr
  TAC_dup      <- TAC

  # MSE function to run for each iter ====
  MSE_fun <- function(x){
    # load packages and functions for each cluster
    library(FLCore)
    library(FLash) 
    library(FLAssess) 
    library(FLXSA)
    
    # create iter versions of objects
    stkObs_nn <- iter(stkObs, x)
    stkReal_nn <- iter(stkReal, x)
    idxObs_nn <- iterIndicies(idxObs, x)
    idxReal_nn <- iterIndicies(idxReal, x)
    TAC_nn <- iter(TAC, x)
    idxErr_nn <- iterIndicies(idxErr, x)
    assErr_nn <- iter(assErr, x)
    catchErr_nn <- iter(catchErr, x)

    for(ii in seq(year_seq)){
      
      # Calculate TAC if year ii is not the final simulation year
      if(ii < nyear){
        TAC_nn <- calcTac(stock=stkObs_nn, tac=TAC_nn, hcrFun=hcrConstF, yr_ass = year_seq[ii]-1, harvest=loop_multi$F[rn])
      }
      
      # advance operational model by one year using TAC
      stkReal_nn <- advanceStock(
        stock = stkReal_nn, tac = TAC_nn, sr = stkRealSr, 
        maxyear = year_seq[ii]-1,
        sr.residuals = iter(trim(recrErr, year=year_seq[ii]), x),
        sr.residuals.mult = TRUE
      )
      
      # record observed catches
      stkObs_nn <- observeCatch(stockReal=stkReal_nn, stockObs=stkObs_nn, obsErr = catchErr)
      
#       # record real and observed indices (only necessary for XSA)
#       idxReal_nn <- observeIndex(stock=stkReal_nn, index=idxReal_nn, q.hat=q.hat, obsErr=NULL, minyear = year_seq[1])
#       idxObs_nn <- observeIndex(stock=stkReal_nn, index=idxObs_nn, q.hat=q.hat, obsErr=idxErr_nn, minyear = year_seq[1])
#       
#       # Assess stock from catch and indices (only XSA)
#       stkObs_nn <- assessStock(
#         stock = stkObs_nn, index = idxObs_nn,
#         control = xsaCtrl,
#         maxyear = year_seq[ii],
#         method = "XSA"
#       )
      
      # Assess stock
      stkObs_nn <- assessStock(
        stock = stkObs_nn,
        minyear = year_seq[ii]-1,
        maxyear = year_seq[ii],
        method = "brute",
        stockReal=stkReal_nn, assessErr=assErr_nn, alpha=1
      )
      
      print(paste("iter =", x, "; year =", year_seq[ii]))
    }

    RES <- list(
      stkObs_nn = stkObs_nn,
      stkReal_nn = stkReal_nn,
      idxObs_nn = idxObs_nn,
      idxReal_nn = idxReal_nn,
      TAC_nn = TAC_nn
    )
    
    sink(file="output.txt", append = TRUE)
    print(paste("run", rn , "iter", x, "completed @", Sys.time()))
    sink()
    
    return(RES)
  }
  
  # Perform multicore MSE ====
  nn <- split(1:niter, 1:niter)
  
  # unlink("output.txt")
  cl <- makeCluster(ncores, type="SOCK", outfile="tmp.txt")
  clusterExport(cl, list=ls(), envir = .GlobalEnv)
  t1 <- Sys.time()
  res <- parLapply(cl, nn, MSE_fun )
  stopCluster(cl)
  Sys.time()- t1
  
  # Rebuild objects with multicore results (by iter) ====
  for(n in nn){
    iter(stkReal, n) <- res[[n]]$stkReal_nn
    iter(stkObs, n) <- res[[n]]$stkObs_nn
    iter(TAC, n) <- res[[n]]$TAC_nn
    for(lev in seq(res[[n]]$idxReal_nn)){
      iter(idxReal[[lev]], n) <- res[[n]]$idxReal_nn[[lev]]
      iter(idxObs[[lev]], n) <- res[[n]]$idxObs_nn[[lev]]
    }
  }
  
  
  
  ### Output / Analysis -------------------------------------------------------
  
  # Save / Load Results ====
  
  SAVE <- TRUE
  if(SAVE){
    rn_pars <- colnames(loop_multi)
    save( stkObs, file=paste0(results_path, "stkObs_", paste(rn_pars, loop_multi[rn,], sep="_", collapse="_"), ".Rdata") )
    save( stkReal, file=paste0(results_path, "stkReal_", paste(rn_pars, loop_multi[rn,], sep="_", collapse="_"), ".Rdata") )
    save( idxObs, file=paste0(results_path, "idxObs_", paste(rn_pars, loop_multi[rn,], sep="_", collapse="_"), ".Rdata") )
    save( idxReal, file=paste0(results_path, "idxReal_", paste(rn_pars, loop_multi[rn,], sep="_", collapse="_"), ".Rdata") )
    save( TAC, file=paste0(results_path, "TAC_", paste(rn_pars, loop_multi[rn,], sep="_", collapse="_"), ".Rdata") )
  }
  
  # restore objects
  stkReal  <- stkReal_dup
  stkObs   <- stkObs_dup
  idxObs   <- idxObs_dup
  idxReal  <- idxReal_dup
  idxErr   <- idxErr_dup
  TAC      <- TAC_dup
  
}
