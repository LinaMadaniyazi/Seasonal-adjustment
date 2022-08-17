################################################################################################################################################
# R codes  for analysis in:
# Should we adjust for season in time-series studies of the short-term association between temperature and mortality?
################################################################################################################################################

################################################################################################################################################
#ASSESS TEMPERATURE-MORTALITY ASSOCIATION
################################################################################################################################################

# LOAD PACKAGES
library(splines); library(tsModel); library(dlnm); library(mgcv); library(scales)


# FUNCTION FOR COMPUTING THE QAIC
QAIC <- function(model) {
  phi <- summary(model)$dispersion
  loglik <- sum(dpois(model$y, model$fitted.values, log=TRUE))
  return(-2*loglik + 2*summary(model)$df[3]*phi)
}


# LOAD THE DATASET
example<- readRDS("example.RDS")


# FORMAT DATE
example$date <- as.Date(example$date)


# GENERATING TIME VARIABLE
example$time <- seq(nrow(example))
numyears<- length(unique(example$year))


# DEFINE CROSSBASIS FOR TEMPERATURE
cb.temp <- crossbasis(example$tmean, lag=21,
                       argvar=list(fun="bs", knots=quantile(example$tmean, c(.10,.75,.90), na.rm=T),degree=2),
                       arglag=list(fun="ns", knots=logknots(21,3)))


# RUN THE TIME SERIES REGRESSIONS, OBTAIN RESIDUALS, AND EXTRACT REDUCED PRED
# MODEL WITH TEMPERATURE ONLY
model0 <- glm(nonext ~ cb.temp, example, family=quasipoisson,na.action = "na.exclude")
# MODEL WITH TEMPERATURE AND TIME TREND
model1<- glm(nonext ~ cb.temp+ns(time,2), example, family=quasipoisson,na.action = "na.exclude")
# MODEL WITH TEMPERATURE, TIME TREND, AND SEASON
model2<- glm(nonext ~ cb.temp+ns(time,8*numyears), example, family=quasipoisson,na.action = "na.exclude")


# OBTAIN QAIC
QAIC(model0); QAIC(model1); QAIC(model2)


# EXTRACT REDUCED PRED AND RE-CENTERED PRED
cp0<- crosspred(cb.temp,model0,by=0.1)
cen0<- cp0$predvar[which.min(cp0$allRRfit)]
pred0<- crosspred(cb.temp, model0, cen=cen0, by=1)

cp1<- crosspred(cb.temp,model1,by=0.1)
cen1<- cp1$predvar[which.min(cp1$allRRfit)]
pred1<- crosspred(cb.temp, model1, cen=cen1, by=1)

cp2<- crosspred(cb.temp,model2,by=0.1)
cen2<- cp2$predvar[which.min(cp2$allRRfit)]
pred2<- crosspred(cb.temp, model2, cen=cen2, by=1)


# PLOT PRED
xlab <- expression(paste("Temperature (",degree,"C)"))
plot(pred0,"overall",col="black",lty=3,ylim=c(0.95,1.8),axes=T,lab=c(6,5,7),xlab=xlab,ylab="Relative Risk",ci="n",lwd=3)
lines(pred1,"overall",col="black",ylim=c(0.95,1.8),axes=T,lab=c(6,5,7),xlab=xlab,lwd=3,lty=5,ylab="Relative Risk",ci="n")
lines(pred2,"overall",col="black",ylim=c(0.95,1.8),axes=T,lab=c(6,5,7),xlab=xlab,lwd=3,ylab="Relative Risk",ci="n")

# PLOT PACF
pacf(resid(model0, type="pearson"), lag.max=30,main = "",ylim=c(-0.1,0.2),na.action = na.pass)
pacf(resid(model1, type="pearson"), lag.max=30,main = "",ylim=c(-0.1,0.2),na.action = na.pass)
pacf(resid(model2, type="pearson"), lag.max=30,main = "",ylim=c(-0.1,0.2),na.action = na.pass)

################################################################################################################################################
# ASSESS SEASONALITY BEFORE AND AFTER TEMPERATURE ADJUSTMENT TO GENERATE FIGURE 4
# Adapted from our previsou tutorial paper:https://doi.org/10.1093/ije/dyac115
# R codes: https://github.com/LinaMadaniyazi/Assessing-seasonality-and-the-role-of-its-potential-drivers-in-environmental-epidemiology-a-tutorial
################################################################################################################################################


