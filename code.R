#######################
#data import
allergy <- read.table("~/project/allergy.txt", header=TRUE, quote="\"")

#######################
#get lsmooth function
lsmooth <- function(x, y, poly=1, h=1.06*sd(x)*length(x)^(-0.2),
                    adjust=1, kernel=function(u) dnorm(u), npoints=100, 
                    xmin=min(x, na.rm=T), xmax=max(x, na.rm=T), varband=T, obsfit=F)
{
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  nn <- length(x)
  h1 <- h*adjust
  mse <- df <- 0
  if(varband)
  {
    #Compute MSE
    Z <- matrix(x, ncol=nn, nrow=nn, byrow=F)
    X <- t(Z)
    W <- kernel((X-Z)/h1)
    U1 <- U2 <- U3 <- U4 <- U5 <- X-X
    U0 <- U1+1
    if (poly>0) U1 <- (X-Z)
    if (poly>1) U2 <- (X-Z)^2
    if (poly>2) U3 <- (X-Z)^3
    if (poly>3) U4 <- (X-Z)^4
    if (poly>4) U5 <- (X-Z)^5
    yhat <- Wmat <- NULL
    for (i in 1:nn)
    {
      Xmat <- cbind(U0[i,],U1[i,],U2[i,],U3[i,],U4[i,],U5[i,])[,1:(poly+1)]
      Vmat <- diag(W[i,])
      Pmat <- solve(t(Xmat)%*%Vmat%*%Xmat)%*%t(Xmat)%*%Vmat
      betahat <- Pmat%*%y
      yhat <- c(yhat, betahat[1])
      Wmat <- rbind(Wmat, Pmat[1,]) 
    }
    sse <- sum((y-yhat)^2)
    df <- nn - sum(diag(Wmat))
    mse <- sse/df
  }
  
  #Compute fitted curve and (possibly) variability bands
  if (!obsfit) z <- seq(xmin, xmax, length=npoints)
  if(obsfit) {z <- x; npoints=nn}
  Z <- matrix(z, ncol=length(x), nrow=npoints, byrow=F)
  X <- matrix(x, ncol=length(x), nrow=npoints, byrow=T)
  U1 <- U2 <- U3 <- U4 <- U5 <- X-X
  U0 <- U1+1
  if (poly>0) U1 <- (X-Z)
  if (poly>1) U2 <- (X-Z)^2
  if (poly>2) U3 <- (X-Z)^3
  if (poly>3) U4 <- (X-Z)^4
  if (poly>4) U5 <- (X-Z)^5
  W <- kernel((X-Z)/h1)
  ssw <- yhat <- NULL
  for (i in 1:npoints)
  {
    Xmat <- cbind(U0[i,],U1[i,],U2[i,],U3[i,],U4[i,],U5[i,])[,1:(poly+1)]
    Vmat <- diag(W[i,])
    Pmat <- solve(t(Xmat)%*%Vmat%*%Xmat)%*%t(Xmat)%*%Vmat
    betahat <- Pmat%*%y
    yhat <- c(yhat, betahat[1])
    ssw <- c(ssw, sum(Pmat[1,]^2)) 
  }
  fitres <- cbind(z, yhat)
  if (varband) fitres <- cbind(fitres, cbind(yhat-1.96*sqrt(mse*ssw), yhat+1.96*sqrt(mse*ssw))) 
  list(data=cbind(x, y), fit=fitres, mse=mse, df=df, bandwidth=h1, kernel=kernel)
}
#######################
#1
cor.test(allergy$age,allergy$bmi,method="kendall")
plot(allergy$age, allergy$bmi,xlab="Age",ylab="BMI",main="BMI-Age Scatterplot")
allergy$agecat[allergy$age<30]<-"30-"
allergy$agecat[allergy$age<40 & allergy$age>=30]<-"30-40"
allergy$agecat[allergy$age<50 & allergy$age>=40]<-"40-50"
allergy$agecat[allergy$age>=50]<-"50+"
agecat.bmi.mean<-aggregate(allergy$bmi, by=list(allergy$agecat), FUN=mean, na.rm=TRUE)
lines(x=c(25,35,45,55),y=agecat.bmi.mean[,2])
legend(20,40,c("BMI Mean by Decade"),lty=1)
colnames(agecat.bmi.mean)<-c("Age Group","BMI Mean")
agecat.bmi.mean

#######################
#2
kruskal.test(day~treat,data=allergy)

#######################
#3
plot(x=allergy$day, y=allergy$logpollen,xlab="Day",ylab="Log-transformed Pollen Count",
     main="Pattern of log-transformed pollen levels over time")
lines(lsmooth(x=allergy$day, y=allergy$logpollen,npoints=100,poly=2)$fit[,1],
      lsmooth(x=allergy$day, y=allergy$logpollen,npoints=100,poly=2)$fit[,2])

#######################
#4
fit<-lsmooth(x=allergy$day, y=allergy$logpollen,obsfit=T,poly=2)$fit[,2]
residual<-as.matrix(allergy$logpollen)-fit
n<-length(residual)
fit.star<-NULL
for (i in 1:300){
  fit.star<-cbind(fit.star,as.matrix(lsmooth(x=allergy$day, 
                                             y=fit+residual[sample(1:n,replace=T),],
                                             obsfit=T,poly=2)$fit[,2] ) )
}
yhatlb<-apply(fit.star,1,quantile,probs=0.025)
yhatup<-apply(fit.star,1,quantile,probs=0.975)
plot<-cbind(as.matrix(allergy$day),as.matrix(yhatlb),as.matrix(yhatup))
colnames(plot) <- c("day", "yhatlb", "yhatup")
plot<-as.data.frame(plot)
plot<-plot[order(plot$day),]
lines(x=plot$day,y=plot$yhatlb,lty=2)
lines(x=plot$day,y=plot$yhatup,lty=2)

data<-allergy[c("day","logpollen")]
fitbt<-NULL
n<-length(data[,1])
for (i in 1:300) {
  tempdata<-data[sample(1:n,replace=T),]
  fitbt<-cbind(fitbt,lsmooth(x=tempdata$day,y=tempdata$logpollen,poly=2,npoints=100)$fit[,2])
}
yhatlb<-apply(fitbt,1,quantile,probs=0.025)
yhatup<-apply(fitbt,1,quantile,probs=0.975)
x0<-lsmooth(x=tempdata$day,y=tempdata$logpollen,poly=2,npoints=100)$fit[,1]
lines(x=x0,y=yhatlb,lty=3)
lines(x=x0,y=yhatup,lty=3)

legend(220,4, c("Original Data","Residual"),lty=c(3,2),cex=0.7)

#######################
#5
allergy$score<-allergy$itchy+allergy$sneezy+allergy$runny+allergy$stuffy
par(mfrow=c(2,2))
for (i in 1:4){
  plot(x=allergy[which(allergy$treat==i),]$logpollen,y=allergy[which(allergy$treat==i),]$score,
       xlab="Logpollen",ylab="Score",main=paste("Treatment",i,"Score Trend",sep=" "))
lines(lsmooth(x=allergy[which(allergy$treat==i),]$logpollen,y=allergy[which(allergy$treat==i),]$score
              ,poly=2,npoints=300)$fit[,1],
      lsmooth(x=allergy[which(allergy$treat==i),]$logpollen,y=allergy[which(allergy$treat==i),]$score
              , poly=2,npoints=300)$fit[,2])
}
par(mfrow=c(1,1))

#######################
#6
pvalue<-NULL
for (i in 1:4){
  mse0<-lsmooth(x=allergy[which(allergy$treat==i),]$logpollen,
                y=allergy[which(allergy$treat==i),]$score
                ,poly=0,npoints=100,h=100)$mse
  df0<-lsmooth(x=allergy[which(allergy$treat==i),]$logpollen,
               y=allergy[which(allergy$treat==i),]$score
               ,poly=0,npoints=100,h=100)$df
  mse1<-lsmooth(x=allergy[which(allergy$treat==i),]$logpollen,
                y=allergy[which(allergy$treat==i),]$score
                ,poly=2,npoints=100)$mse
  df1<-lsmooth(x=allergy[which(allergy$treat==i),]$logpollen,
               y=allergy[which(allergy$treat==i),]$score
               ,poly=2,npoints=100)$df
  q<-((mse0*df0-mse1*df1)/(df0-df1))/mse1
  pvalue<-cbind(pvalue,pf(q, df0-df1, df1,lower.tail = F))
}
pvalue

#######################
#7 
allergy$dayrank<-(rank(allergy$day))
allergy$scorerank<-(rank(allergy$score))
allergy$logpollenrank<-(rank(allergy$logpollen))
allergy$bmirank<-(rank(allergy$bmi))
allergy$agerank<-(rank(allergy$age))
fit<-lm(scorerank~female+dayrank+logpollenrank+bmirank+agerank+treat,data=allergy)
summary(fit)

allergy$pollencat<-cut(allergy$logpollen,breaks=quantile(allergy$logpollen, probs = seq(0, 1, 1/5)))
allergy$bmicat<-cut(allergy$bmi,breaks=quantile(allergy$bmi, probs = seq(0, 1, 1/5)))
allergy$agecat<-cut(allergy$age,breaks=quantile(allergy$age, probs = seq(0, 1, 1/3)))
allergy$daycat<-cut(allergy$day,breaks=quantile(allergy$day, probs = seq(0, 1, 1/5)))
pollencat<-levels(allergy$pollencat)
bmicat<-levels(allergy$bmicat)
agecat<-levels(allergy$agecat)
daycat<-levels(allergy$daycat)
ppollen<-pbmi<-page<-pday<-ptreat<-NULL
for (i in 1:5){
  ppollen<-cbind(ppollen,kruskal.test(allergy[which(allergy$pollencat==pollencat[i]),]$score~
                                        allergy[which(allergy$pollencat==pollencat[i]),]$female)$p.value)
  pbmi<-cbind(pbmi,kruskal.test(allergy[which(allergy$bmicat==bmicat[i]),]$score~
                                  allergy[which(allergy$bmicat==bmicat[i]),]$female)$p.value)
  pday<-cbind(pday,kruskal.test(allergy[which(allergy$daycat==daycat[i]),]$score~
                                  allergy[which(allergy$daycat==daycat[i]),]$female)$p.value)
}
for (i in 1:4){
  ptreat<-cbind(ptreat,kruskal.test(allergy[which(allergy$treat==i),]$score~
                                  allergy[which(allergy$treat==i),]$female)$p.value)
}
for (i in 1:3){
  page<-cbind(page,kruskal.test(allergy[which(allergy$agecat==agecat[i]),]$score~
                                  allergy[which(allergy$agecat==agecat[i]),]$female)$p.value)
}
ppollen
pbmi
pday
page
ptreat

boxplot(allergy[which(allergy$pollencat==pollencat[5]),]$score~
                       allergy[which(allergy$pollencat==pollencat[5]),]$female, xlab="Gender",
        ylab="Score",main="Score Boxplot for Subjects with Logpollen in (3.53,4.72]")
kruskal.test(allergy$score~allergy$female)