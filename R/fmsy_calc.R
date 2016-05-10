# Required packages ====
library(FLCore) # (version: ‘2.5.20160107’)
library(ggplotFL) # (version: ‘2.5.20160119’)

results_path <- "results/msy_xsa/"

# MSY last 5 years
load(file="loop_multi.Rdata")
rn_pars <- colnames(loop_multi)
loop_multi
res <- loop_multi
for(rn in seq(nrow(dat5))){
  load(file=paste0(results_path, "stkReal_", paste(rn_pars, loop_multi[rn,], sep="_", collapse="_"), ".Rdata") )
  nyear <- dims(stkReal)$year
  res$Y[rn] <- median(stkReal@catch[,(nyear-5):nyear])
  res$Fbar[rn] <- median(fbar(stkReal)[,(nyear-5):nyear])
  rm(stkReal)
}


# plot Y ~ F and determin MSY
plot(Y~Fbar, res)
grid()
newdat <- data.frame(Fbar=seq(min(res$Fbar), max(res$Fbar), len=1000))
spl <- smooth.spline(x=res$Fbar, y=res$Y, spar = 0.2)
newdat$Y <- predict(spl,  x=newdat$Fbar)$y
head(newdat)
lines(Y ~ Fbar, newdat)
Fmsy <- newdat$Fbar[which.max(newdat$Y)]
Ymsy <- newdat$Y[which.max(newdat$Y)]
Fmsy; Ymsy
usr <- par()$usr
segments(
  x0=c(Fmsy, Fmsy),
  x1=c(Fmsy,usr[1]),
  y0=c(usr[3],Ymsy),
  y1=c(Ymsy,Ymsy),
  lty=2,
  col=1
)
legend("topright", legend = paste(c("MSY", "FMSY"), "=", c(round(Ymsy), round(Fmsy,3))), bty="n")


# plot of F target versus Fbar
plot(Fbar ~ F, res); abline(0,1, lty=2)


# plot of stock at F=0.4
load(file=paste0(results_path, "stkReal_", paste(rn_pars, loop_multi[4,], sep="_", collapse="_"), ".Rdata") )
plot(stkReal)  


