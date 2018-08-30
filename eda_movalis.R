# 7_IMS data.skv
sales <- read.delim("Movalis/7_IMS data2.skv",header=TRUE,sep=";",dec=",",
                    colClasses=c("character","factor","factor","factor","numeric","numeric",
                                 "factor","character","factor","factor","numeric","factor"
                                 ,"factor","factor","factor","factor","factor","factor","factor"
                                 ,"factor","factor","factor","factor","factor","factor","factor"))
print(dim(sales))
print(names(sales))
sales$Time.AsOfDate <- as.Date(sales$Time.AsOfDate,format="%d.%m.%Y")

# remove single-valued/null fields
sales$ATC.ATC_3 <- NULL # M01A ANTIRHEUMATIC NON-STEROID
sales$Corporation.Producer_Name <- NULL # missing
sales$Market_Name.Market_Name <- NULL # Movalis
sales$OTC_Name.OTC_Name <- NULL # PRESC

# MOVALIS sales by month
table(sales$Brands.TradeName[sales$Brands.Brand=="MOVALIS"])
table(sales$Brands.Brand[sales$Brands.TradeName=="MOVALIS"])

sales_movalis <- aggregate(sales[sales$Brands.Brand=="MOVALIS",c("Units.month","Sales.LC.MTH")],
                           list(sales$Time.AsOfDate[sales$Brands.Brand=="MOVALIS"]),sum)
names(sales_movalis) <- c("Month","Units","Value")

par(mar = c(5, 4, 4, 4))
plot(sales_movalis$Month,sales_movalis$Units/1e5,type="l",col="blue",xlab="",ylab="Units (00,000s)",
     main="Movalis Sales Russia",ylim=c(0,5))
par(new = TRUE)
plot(sales_movalis$Month, sales_movalis$Value/1e6, type="l", col="black", axes=FALSE, 
     xlab="", ylab="", ylim=c(0,200))
axis(side=4, at = pretty(c(0,200)))
mtext("Value (000,000s)", side=4, line=3)

# decomposition
plot(decompose(ts(sales_movalis$Units/1e6,start=c(2011,1),frequency=12)))

# ARIMA forecasting
n <- length(sales_movalis$Units)
print(n)
ntrain <- 36
ntest <- 12
ntrials <- n - (ntrain + ntest) + 1
rmse <- numeric(ntrials)
for(i in 1:ntrials) {
  train <- sales_movalis$Units[i:(i+ntrain-1)]
  test <- sales_movalis$Units[(i+ntrain):(i+ntrain+ntest-1)]
  arima.mod <- try(arima(train,order=c(1,1,1),seasonal=list(order=c(1,1,1),period=12)),
                   silent=TRUE)
  if (class(arima.mod)=="try-error")
    rmse[i] <- NA
  else
    rmse[i] <- sqrt(mean((predict(arima.mod,ntest)$pred[1:ntest]-test)^2))/mean(test)
}
print(mean(rmse,na.rm=TRUE))
print(sd(rmse,na.rm=TRUE)/sqrt(ntrials))

# share of market
sales_total <- aggregate(sales[,c("Units.month","Sales.LC.MTH")],
                           list(sales$Time.AsOfDate),sum)
names(sales_total) <- c("Month","Units","Value")
sales_total <- merge(sales_total,sales_movalis,by="Month")
sales_total$som_units <- sales_total$Units.y/sales_total$Units.x
sales_total$som_value <- sales_total$Value.y/sales_total$Value.x

par(mar = c(5, 4, 4, 4))
plot(sales_total$Month,sales_total$som_units,type="l",col="blue",xlab="",ylab="SOM (Units)",
     main="Movalis Share of Market Russia",ylim=c(0,0.3))
par(new = TRUE)
plot(sales_total$Month, sales_total$som_value, type="l", col="black", axes=FALSE, 
     xlab="", ylab="", ylim=c(0,0.3))
axis(side=4, at = pretty(c(0,0.3)))
mtext("SOM (Value)", side=4, line=3)
legend("topright",c("Units","Value"),col=c("blue","black"),lty=1)

# main competitors by product (TradeName)
setdiff(levels(sales$Brands.TradeName),levels(sales$Brands.Brand))
setdiff(levels(sales$Brands.Brand),levels(sales$Brands.TradeName))
tab <- xtabs(Units.month ~ Brands.TradeName,sales)
stab <- sort(tab,decreasing=TRUE)
majorprods <- names(stab)[cumsum(stab)/sum(tab)<0.9]
sales$Product <- sales$Brands.TradeName
levels(sales$Product) <- c(levels(sales$Brands.TradeName),"Other")
sales$Product[!(sales$Brands.TradeName %in% majorprods)] <- "Other"
sales <- droplevels(sales)
tab1 <- xtabs(Units.month ~ Product,sales)
stab1 <- sort(tab1,decreasing=TRUE)/sum(tab1)
barplot(stab1,las=3,main="Market Share (Units) 01.01.2011 to 01.08.2017")

# competitive VAR (Units)
products <- c("NISE","DICLOFENAC","ORTOPHEN","KETONAL","MOVALIS","MELOXICAM")
sales_all <- aggregate(sales$Units.month,
                          list(sales$Time.AsOfDate,sales$Product),sum)
names(sales_all) <- c("Month","Product","Units")
dat <- sales_all[sales_all$Product=="MOVALIS",c("Month","Units")]
names(dat) <- c("Month","MOVALIS")
for (name in setdiff(products,"MOVALIS")) {
  newdat <- sales_all[sales_all$Product==name,c("Month","Units")]
  validname <- sub(" ",".",name)
  validname <- sub("-",".",name)
  names(newdat) <- c("Month",validname)
  dat <- merge(dat,newdat,by="Month",all.x=TRUE)
}
dat[is.na(dat)] <- 0
dat <- merge(dat,sales_total[,c("Month","Units.x")],by="Month")
for (name in setdiff(names(dat),c("Month","Units.x")))
  dat[,name] <- dat[,name]/dat$Units.x
dat$Units.x <- NULL

library(vars)
rmse <- numeric(ntrials)
for(i in 1:ntrials) {
  train <- dat[i:(i+ntrain-1),2:(length(products))]
  test <- dat[(i+ntrain):(i+ntrain+ntest-1),2:(length(products))]
  var.mod <- VAR(train,p=1,season=12)
  rmse[i] <- sqrt(mean((predict(var.mod,n.ahead=12)$fcst$MOVALIS[,"fcst"]-test$MOVALIS)^2))/mean(test$MOVALIS)
}
print(mean(rmse,na.rm=TRUE))
print(sd(rmse,na.rm=TRUE)/sqrt(ntrials))

# marima
library(marima)
ar <- c(1,2,4)
ma <- c(1)
Mod <- define.model(kvar=length(products),ar=ar,ma=ma)
arp <- Mod$ar.pattern
map <- Mod$ma.pattern

train <- dat[i:(i+ntrain-1),2:(length(products)+1)]
test <- dat[(i+ntrain):(i+ntrain+ntest-1),2:(length(products)+1)]

train.dif <- define.dif(t(train),difference=c(1,1,2,1,3,1,4,1,5,1))$y.dif
test.dif <- define.dif(t(rbind(train[ntrain,],test)),difference=c(1,1,2,1,3,1,4,1,5,1))$y.dif

marima1 <-marima(t(train.dif),ar.pattern=arp,ma.pattern=map,weight=0.66,max.iter=200,penalty=2)
summary(marima1)
plot(2:ntrain,train.dif[5,],type="l",col="blue")
lines(2:ntrain,marima1$fitted[5,],col="red")

p <- arma.forecast(cbind(train.dif,test.dif),nstart=ntrain-1,nstep=ntest,marima=marima1)
tmp <- cbind(c(train[ntrain,5],test[1:(ntest-1),5]),p$forecasts[5,ntrain:(ntrain+ntest-1)])
fitted <- apply(tmp,1,sum)
plot(1:ntest,test[,5],type="l",col="blue")
lines(1:ntest, fitted,col="red")

rmse <- numeric(ntrials)
for(i in 1:ntrials) {
  train <- dat[i:(i+ntrain-1),2:(length(products)+1)]
  test <- dat[(i+ntrain):(i+ntrain+ntest-1),2:(length(products)+1)]
  
  train.dif <- define.dif(t(train),difference=c(1,1,2,1,3,1,4,1,5,1))$y.dif
  test.dif <- define.dif(t(rbind(train[ntrain,],test)),difference=c(1,1,2,1,3,1,4,1,5,1))$y.dif
  
  marima1 <-marima(t(train.dif),ar.pattern=arp,ma.pattern=map,weight=0.66,max.iter=200,penalty=2)
  p <- arma.forecast(cbind(train.dif,test.dif),nstart=ntrain-1,nstep=ntest,marima=marima1)
  tmp <- cbind(c(train[ntrain,5],test[1:(ntest-1),5]),p$forecasts[5,ntrain:(ntrain+ntest-1)])
  fitted <- apply(tmp,1,sum)
  
  rmse[i] <- sqrt(mean((fitted-test$MOVALIS)^2))/mean(test$MOVALIS)
}
print(mean(rmse,na.rm=TRUE))
print(sd(rmse,na.rm=TRUE)/sqrt(ntrials))
print(median(rmse,na.rm=TRUE))
print(median(abs(rmse-median(rmse)))/sqrt(ntrials))
idx <- log(rmse)<0
ntrials1 <- sum(idx)
print(mean(rmse[idx],na.rm=TRUE))
print(sd(rmse[idx],na.rm=TRUE)/sqrt(ntrials1))
min(rmse)

# keep data from 2014-10 as we only have call data from 2014-10
sales <- sales[sales$Time.AsOfDate >= as.Date("01.10.2014",format="%d.%m.%Y"),]

# sales of Movalis by brick by month
sales_mb <- aggregate(sales[sales$Brands.Brand=="MOVALIS","Units.month"],
                      sales[sales$Brands.Brand=="MOVALIS",
                            c("Time.AsOfDate","IMS_Region.IMS_Brick")],sum)
names(sales_mb) <- c("Month","Brick","Units")
boxplot(Units ~ Month, sales_mb, outline=FALSE, main="Distribution of Sales Units over Bricks",
        sub="Movalis, Russia")
print(length(levels(sales_mb$Brick)))



# Calls
#options(java.parameters="-XX:-UseGCOverheadLimit")
#options(java.parameters="-XmX4G")
#library(xlsx)
#calls <- read.xlsx("Movalis/1_Veeva statistics.xlsx",sheetName="Calls",
#                   colClasses=c("factor","factor","factor","character","factor","character",
#                                "numeric"))
#library(gdata)
#calls <- read.xls("Movalis/1_Veeva statistics.xlsx",sheet="Calls",colClasses=
#                    c("factor","factor","factor","character","factor","character","numeric"))

calls <- read.delim("R/Movalis/1_Veeva statistics_calls.txt",header=TRUE,sep="\t",dec=",",
                    colClasses=c("character","factor","factor","character","factor","character",
                                 "numeric","numeric","numeric","numeric","numeric"
                                 ,"numeric","numeric","numeric","numeric","numeric","numeric"
                                 ,"numeric","numeric","numeric","numeric","numeric","numeric"
                                 ,"numeric","numeric","numeric","numeric","numeric","numeric"
                                 ,"numeric","numeric","numeric","numeric","numeric","numeric"
                                 ,"numeric","numeric","numeric","numeric","numeric","numeric"
                                 ,"numeric"),fileEncoding="utf-16")
print(dim(calls))
calls[,7:42][is.na(calls[,7:42])] <- 0
totcalls <- data.frame(data=apply(calls[,7:42],2,sum),row.names=names(calls)[7:42])
period <- as.Date(paste(substring(rownames(totcalls),2),"01",sep="."),format="%Y.%m.%d")
plot(period,totcalls$data,type="l",main="Total calls Movalis RU",ylim=c(0,35000),
     xlab="",ylab="")

# total sales vs total calls
x <- totcalls$data[period %in% dat$Month]
y <- dat$MOVALIS[dat$Month %in% period]
plot(x,y,type="p",xlab="Total calls",
     ylab="Total sales",main="Sales vs Calls Movalis RU")
abline(lm(y~x),col="red")

# adstock 
adstock<-function(x,rate=0){
  return(as.numeric(filter(x=x,filter=rate,method="recursive")))
}

xy<-merge(data.frame(period=period,calls=totcalls$data),sales_movalis,
          by.x="period",by.y="Month")
n<-length(xy$period)
tt<-as.numeric(xy$period-xy$period[1])
#library(minpack.lm)
modFit<-nls(Units~a+b*adstock(calls,phi)+bt*tt,data=xy,
              start=list(a=min(xy$Units),
                         b=mean(xy$Units)/mean(xy$calls),
                         bt=(xy$Units[n]-xy$Units[1])/(tt[n]-tt[1]),
                         phi=0.5),
              lower=c(a=0,
                         b=0,
                         bt=-Inf,
                         phi=0),
              upper=c(a=Inf,
                      b=Inf,
                      bt=Inf,
                      phi=1),
            algorithm = "port")

summary(modFit)$coefficients
#phip<-summary(modFit)$coefficients["phi","Estimate"]
# no adstock

lm2<-lm(Units~calls+tt,data=xy)
summary(lm2)
plot(lm2)

