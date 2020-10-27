############################################
##### Econnomics 480 -- Forecasting ########
##### 'Forecasting Coffee Index Prices' ####
############################################


#############################
########### Setup ###########
#############################

library(sandwich)

projecttext <- read.table ("projecttext.txt", header=TRUE)
attach(projecttext)

### Variables ###
#coffee_cpi
#price_arabica
#price_tea
#price_milk

##Log Variables##
Lcoffee <- log(coffee_cpi)
Larabica <- log(price_arabica)
Ltea <- log(price_tea)
Lmilk <- log(price_milk)

## First Differences ##
Tcoffee <- length(Lcoffee) 
DLcoffee <- Lcoffee[2:Tcoffee]-Lcoffee[1:(Tcoffee-1)]
T2coffee <- length(DLcoffee)

DLarabica <- Larabica[2:Tcoffee]-Larabica[1:(Tcoffee-1)]

DLtea <- Ltea[2:Tcoffee]-Ltea[1:(Tcoffee-1)]

DLmilk <- Lmilk[2:Tcoffee]-Lmilk[1:(Tcoffee-1)]


##Dickey-Fuller##
maxlag <- 10

df_dep <- DLcoffee[(maxlag+1):T2coffee]
df_indep <- Lcoffee[(maxlag+1):T2coffee]
constant <- matrix(1,length(df_dep),1)
timetrend <- seq(1,length(df_dep),by=1)
xmat <- cbind(constant,timetrend,df_indep)

for (i in 1:maxlag) {
	xmat <- cbind(xmat,DLcoffee[(maxlag+1-i):(T2coffee-i)])
	}

aic_mat <- matrix(0,(maxlag+1),1)
df_mat <- matrix(0,(maxlag+1),1)

for (i in 1:(maxlag+1)) {
	xreg <- xmat[,1:(2+i)]
	reg <- lm(df_dep~xreg+0)
	vcov <- vcovHC(reg,type=c("HC"))
	se <- sqrt(diag(vcov))
	aic_mat[i] <- AIC(reg)
	df_mat[i] <- reg$coef[3]/se[3]
	}

df_statistic <- df_mat[which(aic_mat==min(aic_mat))]
df_statistic

######################################
############ Rolling #################
######################################

start_window <- 50
forecast_num <- T2coffee-start_window-1

forecasts_AR1 <- matrix(0,forecast_num,1)
forecasts_DL1 <- matrix(0,forecast_num,1)
forecasts_ADL1 <- matrix(0,forecast_num,1)
forecasts_MP1 <- matrix(0,forecast_num,1)

forecast_errors_AR1 <- matrix(0,forecast_num,1)
forecast_errors_DL1 <- matrix(0,forecast_num,1)
forecast_errors_ADL1 <- matrix(0,forecast_num,1)
forecast_errors_MP1 <- matrix(0,forecast_num,1)

for (i in 1:forecast_num) {
	yf <- DLcoffee[(1+i):(start_window+i)]
	xf <- DLcoffee[(i):(start_window+i-1)]
	af <- DLarabica[(i):(start_window+i-1)]
	tf <- DLtea[(i):(start_window+i-1)]
	mf <- DLmilk[(i):(start_window+i-1)]
	AR1 <- lm(yf~xf)
	DL1 <- lm(yf~af)
	ADL1 <- lm(yf~xf+af)
	MP1 <- lm(yf~xf+af+tf+mf)

	coefsAR <- AR1$coef
	coefsDL <- DL1$coef
	coefsADL <- ADL1$coef
	coefsMP <- MP1$coef

	forecasts_AR1[i] <- coefsAR[1]+coefsAR[2]*yf[length(yf)]
	forecasts_DL1[i] <- coefsDL[1]+coefsDL[2]*DLarabica[(start_window+i)]
	forecasts_ADL1[i] <- coefsADL[1]+coefsADL[2]*yf[length(yf)]+coefsADL[3]*DLarabica[(start_window+i)]
	forecasts_MP1[i] <- coefsMP[1]+coefsMP[2]*yf[length(yf)]+coefsMP[3]*DLarabica[(start_window+i)]+coefsMP[4]*DLtea[(start_window+i)]+coefsMP[5]*DLmilk[(start_window+i)]

	forecast_errors_AR1[i] <- DLcoffee[(start_window+i+1)]-forecasts_AR1[i]
	forecast_errors_DL1[i] <- DLcoffee[(start_window+i+1)]-forecasts_DL1[i]
	forecast_errors_ADL1[i] <- DLcoffee[(start_window+i+1)]-forecasts_ADL1[i]
	forecast_errors_MP1[i] <- DLcoffee[(start_window+i+1)]-forecasts_MP1[i]
	}

sq_forecast_errors_AR1 <- forecast_errors_AR1*forecast_errors_AR1
sq_forecast_errors_DL1 <- forecast_errors_DL1*forecast_errors_DL1
sq_forecast_errors_ADL1 <- forecast_errors_ADL1*forecast_errors_ADL1
sq_forecast_errors_MP1 <- forecast_errors_MP1*forecast_errors_MP1

sq_mean_AR1 <-sqrt(mean(sq_forecast_errors_AR1))
sq_mean_DL1 <-sqrt(mean(sq_forecast_errors_DL1))
sq_mean_ADL1 <-sqrt(mean(sq_forecast_errors_ADL1))
sq_mean_MP1 <-sqrt(mean(sq_forecast_errors_MP1))

##Run##
sq_mean_AR1 
sq_mean_DL1 
sq_mean_ADL1 
sq_mean_MP1 


#######################################
########### Diebold-Mariano ###########
#######################################

dm_quad_ARvsDL <- sq_forecast_errors_AR1-sq_forecast_errors_DL1
dm_quad_ARvsADL <- sq_forecast_errors_AR1-sq_forecast_errors_ADL1
dm_quad_ARvsMP <- sq_forecast_errors_AR1-sq_forecast_errors_MP1
dm_quad_DlvsADL <- sq_forecast_errors_DL1-sq_forecast_errors_ADL1
dm_quad_DLvsMP <- sq_forecast_errors_DL1-sq_forecast_errors_MP1
dm_quad_ADLvsMP <- sq_forecast_errors_ADL1-sq_forecast_errors_MP1

dm1 <- lm(dm_quad_ARvsDL~1)
dm2 <- lm(dm_quad_ARvsADL~1)
dm3 <- lm(dm_quad_ARvsMP~1)
dm4 <- lm(dm_quad_DlvsADL~1)
dm5 <- lm(dm_quad_DLvsMP~1)
dm6 <- lm(dm_quad_ADLvsMP~1)

vcov1<- vcovHC(dm1,type=c("HC"))
vcov2<- vcovHC(dm2,type=c("HC"))
vcov3<- vcovHC(dm3,type=c("HC"))
vcov4<- vcovHC(dm4,type=c("HC"))
vcov5<- vcovHC(dm5,type=c("HC"))
vcov6<- vcovHC(dm6,type=c("HC"))

Std_Error_1 <- sqrt(diag(vcov1))
Std_Error_2 <- sqrt(diag(vcov2))
Std_Error_3 <- sqrt(diag(vcov3))
Std_Error_4 <- sqrt(diag(vcov4))
Std_Error_5 <- sqrt(diag(vcov5))
Std_Error_6 <- sqrt(diag(vcov6))

##Run##
Std_Error_1 
Std_Error_2 
Std_Error_3 
Std_Error_4 
Std_Error_5 
Std_Error_6 

#------------#

t_stat_1<- dm1$coef[1]/Std_Error_1[1]
t_stat_2<- dm2$coef[1]/Std_Error_2[1]
t_stat_3<- dm3$coef[1]/Std_Error_3[1]
t_stat_4<- dm4$coef[1]/Std_Error_4[1]
t_stat_5<- dm5$coef[1]/Std_Error_5[1]
t_stat_6<- dm6$coef[1]/Std_Error_6[1]

##Run##

t_stat_1 
t_stat_2
t_stat_3
t_stat_4
t_stat_5
t_stat_6

### Forecast ###
forecasts_final <- coefsDL[1]+coefsDL[2]*DLarabica[(449)]
forecasts_final

#200.233 Actual
#200.824 Last
#200.872 Predicted

summary(AR1)
summary(DL1)
summary(ADL1)
summary(MP1)