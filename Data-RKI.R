if (plot_out == 1) pdf("Plots-Data-RKI.pdf")

#####################################################################################
#
# Data
#   
#####################################################################################

rki_data_raw <- read.delim("https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv", sep=",", header=T)

#rki_data <- as.data.frame(rki_data_raw[rki_data$Bundesland=="Sachsen-Anhalt",])
rki_data <- rki_data_raw

start_date <- min(strptime(rki_data$Refdatum, "%Y/%m/%d %H:%M:%S", tz="CET"))
start_date_2 <- strptime("2020/3/1 0:00:00", "%Y/%m/%d %H:%M:%S", tz="CET")

dates <- strptime(levels(rki_data$Meldedatum), "%Y/%m/%d %H:%M:%S", tz="CET")
days <- as.numeric(dates - start_date, units="days")

ages <- levels(rki_data$Altersgruppe)
ages[7] <- "unknown"
levels(rki_data$Altersgruppe) <- ages


res <- aggregate(data.frame(ConfirmedCases=rki_data$AnzahlFall, RecoverdCases=rki_data$AnzahlGenesen, DeathCases=rki_data$AnzahlTodesfall, IsBeginning=rki_data$IstErkrankungsbeginn*rki_data$AnzahlFall),
                 by=data.frame(Time=rki_data$Meldedatum, Age=rki_data$Altersgruppe),
                 FUN=sum)


res$Time <- strptime(res$Time, "%Y/%m/%d %H:%M:%S", tz="CET")
res$Time  <- as.numeric(res$Time - start_date, units="days")

t <- as.numeric(strptime(rki_data$Meldedatum, "%Y/%m/%d %H:%M:%S", tz="CET")-start_date, units="days")
tmp <- rki_data[ t> 130 & t<140,]
sum(tmp$AnzahlFall)

#####################################################################################
#
# Constants I
#   
#####################################################################################

# introduction of dexamethasone
t_0 = as.numeric(strptime("2020/3/1 00:00:00", "%Y/%m/%d %H:%M:%S", tz="CET") - start_date, units="days")
#tm = 107 + 7
tm = strptime("2020/06/16 13:00:00", "%Y/%m/%d %H:%M:%S", tz="CET") - start_date + 7
phi_m = 0.65*0.41+0.8*0.25+(1-0.41-0.25)*1
phi_m = 0.83

# confirmed to death (determined by overlapping the two curves)
td <-10.5
# assume incubation time from https://www.ncbi.nlm.nih.gov/pubmed/32150748
te = 5.1
ti = 14 # determined by overlapping the two curves
th = 8 
thi = 10 
thd = td - th 
thii = 10

#####################################################################################
#
# Sitation Analysis
#   
#####################################################################################
################################
# Confirmed Cases
################################
if (plot_out == 2) png("Confirmed_Cases-1.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)

plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),")"), ylab="Confirmed Cases", xlim = c(0,max(res$Time)), ylim = c(0,max(res$ConfirmedCases)))
for(i in 1:7)
{
  sdata<-subset(res,Age==ages[i])
  points(sdata$Time, sdata$ConfirmedCases, col=i)
  lines(lowess(sdata$Time, sdata$ConfirmedCases,f=0.03), col=i)
}

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topright", legend=ages, col=c(1:7), bty="n", pch=1)
box()

if (plot_out != 2) title("Confirmed Cases COVID-19 in Germany", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()

################################
# Death Cases
################################

if (plot_out == 2) png("Death_Cases-1.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)

plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),")"), ylab="Death Cases", xlim = c(0,max(res$Time)), ylim = c(0,250))
for(i in 1:7)
{
  sdata<-subset(res,Age==ages[i])
  points(sdata$Time, sdata$DeathCases, col=i)
  lines(lowess(sdata$Time, sdata$DeathCases,f=0.05), col=i)
}

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topright", legend=ages, col=c(1:7), bty="n", pch=1)
box()

if (plot_out != 2) title("Death Cases COVID-19 in Germany", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()

################################
# Time from Confirmed to Death
#######################-#########

if (plot_out == 2) png("Time_from_Confirmed_to_Death.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)

if(length(tcd) <= 1)
{
  tcd_i <- -4 -1 * c(0:300)/100
  qs <- array(NA,dim=c(7,length(tcd_i)))
  for(i in c(1,3:5))
  {
    sdata<-subset(res,Age==ages[i])
    l <- length(sdata$Time)
    f <- c(3:l)-1
    
    sd_tcd <- c(
      ((sdata$Time[1] - sdata$Time[2])^2)^0.5,
      (((sdata$Time[f-1] - sdata$Time[f])^2 + (sdata$Time[f] - sdata$Time[f+1])^2)/2)^0.5,
      ((sdata$Time[l-1] - sdata$Time[l])^2)^0.5
    )
    
    for (iii in  1:length(tcd_i))
    {
      tcd <- tcd_i[iii]
      t=c(1:l)
      r=t
      for(ii in t)
      {
        t[ii] <- sdata$Time[ii]-tcd
        s  <- sum(exp(-0.5*((sdata$Time-t[ii])/sd_tcd)^2) )
        ss <- sum(exp(-0.5*((sdata$Time-t[ii])/sd_tcd)^2) * sdata$DeathCases ) 
        r[ii] <- ss/s
      }
      rr <- (sdata$ConfirmedCases/quantile(sdata$ConfirmedCases,probs=0.9) - r/quantile(r,probs=0.9))^2
      qs[i,iii] <- sum(rr[sdata$Time < 200 & sdata$Time > 50])
    }
  }
  tcd <-  1:7
  for(i in tcd) tcd[i] <- tcd_i[min(qs[i,])==qs[i,]][1]
  
  tcd[c(2,6,7)]<- c(mean(tcd[c(1,3)]),0,mean(tcd[c(1,3:5)]))
} 


plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),"). Deaths Cases are shifted"), ylab="Confirmed Cases vs Normalized Deaths", xlim = c(0,max(res$Time)), ylim = c(0,max(res$ConfirmedCases)))
for(i in 1:6)
{
  sdata<-subset(res,Age==ages[i])
  l <- length(sdata$Time)
  f <- c(3:l)-1
  
  sd_tcd <- c(
    ((sdata$Time[1] - sdata$Time[2])^2)^0.5,
    (((sdata$Time[f-1] - sdata$Time[f])^2 + (sdata$Time[f] - sdata$Time[f+1])^2)/2)^0.5,
    ((sdata$Time[l-1] - sdata$Time[l])^2)^0.5
  )
  
  if(!is.na(tcd[i]))
  {
    t=c(1:l)
    r=t
    for(ii in t)
    {
      t[ii] <- sdata$Time[ii]-tcd[i]
      s  <- sum(exp(-0.5*((sdata$Time-t[ii])/sd_tcd)^2) )
      ss <- sum(exp(-0.5*((sdata$Time-t[ii])/sd_tcd)^2) * sdata$DeathCases ) 
      r[ii] <- ss/s
    }
    points(sdata$Time, r/quantile(r,probs=0.9)*quantile(sdata$ConfirmedCases,probs=0.9), col=i)
  }
  lines(lowess(sdata$Time, sdata$ConfirmedCases,f=0.03), col=i)
}

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topright", legend=ages, col=c(1:7), bty="n", pch=1)
box()

if (plot_out != 2) title("Time from Confirmed to Death COVID-19 in Germany", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()

################################
# Death Cases per Day
################################

if (plot_out == 2) png("Death_Cases_per_Day-1.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)

rr <- list()
plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),")"), ylab="Death Cases (%)", xlim = c(0,max(res$Time)), ylim = c(0,50))
for(i in 1:7)
{
  sdata<-subset(res,Age==ages[i])
  l <- length(sdata$Time)
  f <- c(3:l)-1
  
  sd_tcd <- c(
    ((sdata$Time[1] - sdata$Time[2])^2)^0.5,
    (((sdata$Time[f-1] - sdata$Time[f])^2 + (sdata$Time[f] - sdata$Time[f+1])^2)/2)^0.5,
    ((sdata$Time[l-1] - sdata$Time[l])^2)^0.5
  )
  
  if(!is.na(tcd[i]))
  {
    t=c(1:l)
    r2=t2=r1=r=t
    for(ii in t)
    {
      t[ii] <- sdata$Time[ii]
      t2[ii] <- t[ii]+tcd[i]
      s  <- sum(exp(-0.5*((sdata$Time-t[ii])/sd_tcd)^2) )
      ss <- sum(exp(-0.5*((sdata$Time-t[ii])/sd_tcd)^2) * sdata$DeathCases ) 
      r[ii] <- ss/s
      ss <- sum(exp(-0.5*((sdata$Time-t[ii])/sd_tcd)^2) * sdata$ConfirmedCases ) 
      r1[ii] <- ss/s
      s  <- sum(exp(-0.5*((sdata$Time-t2[ii])/sd_tcd)^2) )
      ss <- sum(exp(-0.5*((sdata$Time-t2[ii])/sd_tcd)^2) * sdata$ConfirmedCases ) 
      r2[ii] <- ss/s
    }
    r3 <- (r+1e-12)/r2
    m <- array(1,length(sdata$Time))
    m[sdata$Time > tm] <- phi_m
    
    s <- (r3>0) 
    ss <- c(10:17)+t_0
    r4 <- r3/quantile((r3/m)[ss], probs=0.5, na.rm=T)/m
    #ss <- (r3>0) & (r > .1)
    #r4 <- r3/quantile((r3/m)[ss], probs=0.25)/m
    d <- (quantile(r4[s], na.rm=T, probs=3/4)-quantile(r4[s], na.rm=T, probs=1/4))/(quantile(t[s], na.rm=T, probs=3/4)-quantile(t[s], na.rm=T, probs=1/4))
    r4 <- r4 + d * t * m
    rr[[i]] <- list("Time"=t,"DeathCases" = r,"ConfirmedCases" = r1,"TimeShifted"=t2,"ConfirmedCasesShifted" = r2, "DeathFraction" = r3, "DeathFractionNormalized" = r4, "DeathFractionNormalized" = r4, "Drift" = d)
    points(rr[[i]]$Time, rr[[i]]$DeathFraction*100, col=i)
    lines(lowess(rr[[i]]$Time, rr[[i]]$DeathFraction*100, f=0.1), col=i)
  }
}



for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topright", legend=ages, col=c(1:7), bty="n", pch=1)

box()

phi_d_drift <- phi_d_median <- phi_d_mean <- c(1:8)
s1 <- 0
s2 <- 0
ss <- 0
for(i in 1:7) 
{
  m <- array(1,length(rr[[i]]$Time))
  m[rr[[i]]$Time > tm] <- phi_m
  phi_d_mean[i] <- mean(rr[[i]]$DeathFraction/m)
  phi_d_median[i] <- median(rr[[i]]$DeathFraction/m)
  s1 <- s1 + sum(rr[[i]]$ConfirmedCasesShifted)*phi_d_mean[i]
  s2 <- s2 + sum(rr[[i]]$ConfirmedCasesShifted)*phi_d_median[i]
  ss <- ss + sum(rr[[i]]$ConfirmedCasesShifted)
  phi_d_mean[i] <- s1/ss
  phi_d_median[i] <- s2/ss
  phi_d_drift[i] <- rr[[i]]$Drift
  lines(c(-1e9,tm,tm,1e9), c(phi_d_mean[i],phi_d_mean[i],phi_d_mean[i]*phi_m,phi_d_mean[i]*phi_m)*100,lty=3,col=i)
}
box()
if (plot_out != 2) title("Death Cases per Day COVID-19 in Germany", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()


plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),")"), ylab="Death Cases Normalized", xlim = c(0,max(res$Time)), ylim = c(0,15))
for(i in 1:7) 
{
  points(rr[[i]]$Time, rr[[i]]$DeathFractionNormalized, col=i)
  lines(lowess(rr[[i]]$Time, rr[[i]]$DeathFractionNormalized, f=0.08), col=i)
}

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topright", legend=ages, col=c(1:7), bty="n", pch=1)

box()

################################
# Death Cases per Day (II)
################################
if (plot_out == 2) png("Death_Cases_per_Day-2.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)

plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),")"), ylab="Death Cases (%)", xlim = c(t_0,max(res$Time)), ylim = c(0,1))
for(i in 1:7)
{
  points(rr[[i]]$Time, rr[[i]]$DeathFraction*100, col=i)
  lines(lowess(rr[[i]]$Time, rr[[i]]$DeathFraction*100, f=0.1), col=i)
}

phi_d <- c(1:6)
for(i in phi_d) 
{
  m <- array(1,length(rr[[i]]$Time))
  m[rr[[i]]$Time > tm+t_0] <- phi_m
  phi_d[i] <- mean(rr[[i]]$DeathFraction/m)
  lines(c(-1e9,tm,tm,1e9), c(phi_d[i],phi_d[i],phi_d[i]*phi_m,phi_d[i]*phi_m)*100,lty=3,col=i)
}

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topright", legend=ages, col=c(1:7), bty="n", pch=1)
box()

if (plot_out != 2) title("Death Cases per Day COVID-19 in Germany", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()

################################
# Death Cases per Day (III)
################################
if (plot_out == 2) png("Death_Cases_per_Day-3.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)

plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),")"), ylab="Death Cases (%)", xlim = c(0,max(res$Time)), ylim = c(0,10))

phi_d = t = days

for(ii in 1:length(dates))
{
  s=0
  ss=0
  for(i in 1:7)
  {

    l <- length(rr[[i]]$Time)
    f <- c(3:l)-1
    sd_tcd <- c(
      ((rr[[i]]$Time[1] - rr[[i]]$Time[2])^2)^0.5,
      (((rr[[i]]$Time[f-1] - rr[[i]]$Time[f])^2 + (rr[[i]]$Time[f] - rr[[i]]$Time[f+1])^2)/2)^0.5,
      ((rr[[i]]$Time[l-1] - rr[[i]]$Time[l])^2)^0.5
    )

    s  = s + sum(exp(-0.5*((rr[[i]]$Time-t[ii])/sd_tcd)^2) * rr[[i]]$ConfirmedCasesShifted)
    ss = ss + sum(exp(-0.5*((rr[[i]]$Time-t[ii])/sd_tcd)^2) * rr[[i]]$DeathFraction * rr[[i]]$ConfirmedCasesShifted)
  }
  phi_d[ii] <- ss/s
}

m <- array(1,length(t))
m[t > tm] <- phi_m
phi_d_mean[8] <- mean(phi_d/m)
phi_d_median[8] <- median(phi_d/m)

s <- (phi_d>0) 
ss <- c(10:17)+t_0
r4 <- phi_d/quantile((phi_d/m)[ss], probs=0.5)/m
phi_d_drift[8] <- (quantile(r4[s], na.rm=T, probs=3/4)-quantile(r4[s], na.rm=T, probs=1/4))/(quantile(t[s], na.rm=T, probs=3/4)-quantile(t[s], na.rm=T, probs=1/4))

plot(days,phi_d)
points(days,phi_d_mean[8] * ( phi_d_drift[8]*(max(days)-days))*m, col=2)

#phi_d[length(t)-c(14:0)] = pmax(phi_d[length(t)-c(14:0)], mean(phi_d[length(t)-c(21:14)]))

data <- data.frame(x=t, y=phi_d*100)
mean <- lowess(data,f=0.1)
sd <- data.frame(x=data$x,y=(mean$y - data$y)^2)
sd$y[2:(length(sd$x)-1)] <- (sd$y[1:(length(sd$x)-2)] + sd$y[2:(length(sd$x)-1)] + sd$y[3:(length(sd$x))])/2
sd$y[c(1,length(sd$x))]<- sd$y[c(1,length(sd$x))] *3/2
sd <- lowess(sd, f=0.15)
sd$y <- abs(sd$y)^0.5*qt(0.975,length(sd$y))

polygon(c(mean$x,rev(mean$x)),c(mean$y-sd$y,rev(mean$y+sd$y)),lty=0,col="gray95")
lines(mean$x,mean$y-sd$y,lty=2)
lines(mean$x,mean$y+sd$y,lty=2)

points(data)
lines(mean,lty=2,col=2)

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

box()

if (plot_out != 2) title("Death Cases per Day COVID-19 in Germany", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()

################################
# Hidden Cases per Day
################################

if (plot_out == 2) png("Hidden_Cases_per_Day.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)

plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),")"), ylab="Hidden Cases Per Day (%)", xlim = c(min(days),max(res$Time)), ylim = c(0,500))

phi_hidden = t = days

for(ii in 1:length(dates))
{
  s=0
  ss=0
  
  for(i in 4:6)
  {
    l <- length(rr[[i]]$Time)
    f <- c(3:l)-1
    sd_tcd <- c(
      ((rr[[i]]$Time[1] - rr[[i]]$Time[2])^2)^0.5,
      (((rr[[i]]$Time[f-1] - rr[[i]]$Time[f])^2 + (rr[[i]]$Time[f] - rr[[i]]$Time[f+1])^2)/2)^0.5,
      ((rr[[i]]$Time[l-1] - rr[[i]]$Time[l])^2)^0.5
    )
    
    #s  = s + sum(exp(-0.5*((rr[[i]]$Time-t[ii])/sd_tcd)^2) * (rr[[i]]$DeathCases + 1))
    #ss = ss + sum(exp(-0.5*((rr[[i]]$Time-t[ii])/sd_tcd)^2) * (rr[[i]]$DeathFractionNormalized * rr[[i]]$DeathCases + 1))
    s  = s + sum(exp(-0.5*((rr[[i]]$Time-t[ii])/sd_tcd)^2) * (rr[[i]]$DeathCases))
    ss = ss + sum(exp(-0.5*((rr[[i]]$Time-t[ii])/sd_tcd)^2) * (rr[[i]]$DeathFractionNormalized * rr[[i]]$DeathCases))
  }
  phi_hidden[ii] <- ss/s
}

data <- data.frame(x=t, y=phi_hidden*100)
mean <- lowess(data,f=0.1)
sd <- data.frame(x=data$x,y=(mean$y - data$y)^2)
sd$y[2:(length(sd$x)-1)] <- (sd$y[1:(length(sd$x)-2)] + sd$y[2:(length(sd$x)-1)] + sd$y[3:(length(sd$x))])/2
sd$y[c(1,length(sd$x))]<- sd$y[c(1,length(sd$x))] *3/2
sd <- lowess(sd, f=0.15)
sd$y <- abs(sd$y)^0.5*qt(0.975,length(sd$y))

phi_no_IsBeginning <- aggregate(res$ConfirmedCases, by=list(res$Time), FUN=sum)$x/aggregate(res$IsBeginning, by=list(res$Time), FUN=sum)$x
phi_no_IsBeginning[is.infinite(phi_no_IsBeginning)] <- 1

data2 <- data.frame(x=t, y=phi_no_IsBeginning*100)
mean2 <- lowess(data2,f=0.1)
sd2 <- data.frame(x=data2$x,y=(mean2$y - data2$y)^2)
sd2$y[2:(length(sd2$x)-1)] <- (sd2$y[1:(length(sd2$x)-2)] + sd2$y[2:(length(sd2$x)-1)] + sd2$y[3:(length(sd2$x))])/2
sd2$y[c(1,length(sd2$x))]<- sd2$y[c(1,length(sd2$x))] *3/2
sd2 <- lowess(sd2, f=0.15)
sd2$y <- abs(sd2$y)^0.5*qt(0.975,length(sd2$y))

polygon(c(mean$x,rev(mean$x)),c(mean$y-sd$y,rev(mean$y+sd$y)),lty=0,col="gray95")
polygon(c(mean2$x,rev(mean2$x)),c(mean2$y-sd2$y,rev(mean2$y+sd2$y)),lty=0,col="gray95")

lines(mean$x,mean$y-sd$y,lty=2, col="red")
lines(mean$x,mean$y+sd$y,lty=2, col="red")
lines(mean2$x,mean2$y-sd2$y,lty=2, col="blue")
lines(mean2$x,mean2$y+sd2$y,lty=2, col="blue")

points(data, col="rosybrown1")
points(data2, col="deepskyblue")

lines(mean,lwd=1.5,lty=2,col="red")
lines(mean2,lwd=1.5,lty=2,col="blue")
# Bad Feilnbach
lines(c(min(data$x),max(data$x)), c(260,260),lwd=1.5,lty=2,col="darkgreen")
# Kupferzell
lines(c(min(data$x),max(data$x)), c(390,390),lwd=1.5,lty=2,col="orange")

lines(c(-1e9,1e9), c(100,100),lty=3,col="darkgrey")

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topleft", legend=c("by death cases", "by untracted cases", "Corona-Monitoring lokal: Bad Feilnbach", "Corona-Monitoring lokal: Kupferzell "),  col=c("red","blue", "darkgreen","orange"), lwd=2) 


box()

if (plot_out != 2) title("Hidden Cases per Day COVID-19 in Germany", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()

################################
# Cases Cumulative
################################

if (plot_out == 2) png("Cases_Cumulative-1.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)

plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),")"), ylab="Cumulative Cases", xlim = c(t_0,max(res$Time)), ylim = c(0,3e5))

cumulative = list()

for(i in 1:8)
{
  t <- days
  s = array(0,dim=c(3,length(t)))
  
  for(ii in 1:length(t))
  {
    c = d = 0
    if(i < 8)
    {
      c = d = 0 
      f <- rr[[i]]$Time == t[ii]
      if(any(f))
      {
        c <- sum(rr[[i]]$ConfirmedCases[f])
        d <- sum(rr[[i]]$DeathCases[f])
      }
    } 
    
    if(i == 8)
    {
      c = d = 0 
      for(iii in 1:7)
      {
        f <- rr[[iii]]$Time == t[ii]
        if(any(f))
        {
          c <- sum(c, rr[[iii]]$ConfirmedCases[f])
          d <- sum(d, rr[[iii]]$DeathCases[f])
        }
      }
      print(paste(ii, c))
    }
    
    if(ii == 1) 
    { 
      s[1,ii] = c 
    } else 
    {
      s[1,ii] = s[1,ii-1] + s[1,ii] + c
    }
    
    if(ii == 1) 
    {
      s[3,ii] = d
    } else 
    {
      s[3,ii] = s[3,ii-1] + s[3,ii] + d
    }
    
    if(ii + ti <= length(t))
    {
      s[2, ii+14] = s[2,ii+14] + c
    }
    
    if(ii == 1)
    {
      s[2, ii] = s[2,ii] - d
    } else 
    {
      s[2, ii] = s[2,ii-1] + s[2,ii] - d
    }
    
    if (s[2, ii] <= 0) 
    {
      s[2, ii] = 0
      s[2, ii] = s[2, ii+1] + s[2, ii]
    }
  }
  
  cumulative[[i]] <- list("Time"=t, "ConfirmedCases"=s[1,], "RecoverdCases"=s[2,], "DeathCases"=s[3,])
  for(ii in 1:3) lines(t, s[ii,], col=i, lty=ii)
}


legend("topleft", legend=c(ages,"Confirmed cases (RKI)", "Recovered cases", "Death cases (RKI)"), col=c(1:7,rep(1,3)), bty="n", pch=-1, lty=c(rep(1,7),1:3) )

for(i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

box()

if (plot_out != 2) title("Cases Cumulative COVID-19 in Germany", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()

################################
# Cases Cumulative Extrapolated
################################

if (plot_out == 2) png("Cases_Cumulative-2.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)

plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),")"), ylab="Extrapolated Total Cumulative Cases", xlim = c(t_0,max(res$Time)), ylim = c(0,0.8e6))

cumulative_extrapolated = list()

for(i in 1:8)
{
  t <- days
  s = array(0,dim=c(3,length(t)))
  
  for(ii in 1:length(t))
  {
    
    if(i < 8)
    {
      c = d = 0 
      f <- rr[[i]]$Time == t[ii]
      if(any(f))
      {
        c <- sum(rr[[i]]$ConfirmedCases[f] * phi_hidden[ii])
        d <- sum(rr[[i]]$DeathCases[f])
      }
    } 
    
    if(i == 8)
    {
      c = d = 0 
      for(iii in 1:7)
      {
        f <- rr[[iii]]$Time == t[ii]
        if(any(f))
        {
          c <- sum(c, rr[[iii]]$ConfirmedCases[f] * phi_hidden[ii]) 
          d <- sum(d, rr[[iii]]$DeathCases[f])
        }
      }
    }
    
    if(ii == 1) 
    { 
      s[1,ii] = c 
    } else 
    {
      s[1,ii] = s[1,ii-1] + s[1,ii] + c
    }
    
    if(ii == 1) 
    {
      s[3,ii] = d
    } else 
    {
      s[3,ii] = s[3,ii-1] + s[3,ii] + d
    }
    
    if(ii + ti <= length(t))
    {
      s[2, ii+14] = s[2,ii+14] + c
    }
    
    if(ii == 1)
    {
      s[2, ii] = s[2,ii] - d
    } else 
    {
      s[2, ii] = s[2,ii-1] + s[2,ii] - d
    }
    
    if (s[2, ii] <= 0) 
    {
      s[2, ii] = 0
      s[2, ii] = s[2, ii+1] + s[2, ii]
    }
  }
  
  cumulative_extrapolated[[i]] <- list("Time"=t, "ConfirmedCases"=s[1,], "RecoverdCases"=s[2,], "DeathCases"=s[3,])
  for(ii in 1:3) lines(t, s[ii,], col=i, lty=ii)
}


legend("topleft", legend=c(ages,"all ages","confirmed cases (extrapolated)", "recovered and hidden recovered cases (extrapolated)", "death cases (RKI)"), col=c(1:8,rep(1,3)), bty="n", pch=-1, lty=c(rep(1,8),1:3) )

for(i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

box()

if (plot_out != 2) title("Cases Cumulative Extrapolated COVID-19 in Germany", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()

################################
# Basic Reproductive Number R0
################################

if (plot_out == 2) png("Basic_Reproductive_Number_R0-1.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)

n_R0 <- 3
R0 <- list()

for(ii in 1:8)
{
  t <- cumulative[[ii]]$Time[cumulative[[ii]]$Time-te>0]
  r=0
  t <- t[-1]
  c1 <- approx(cumulative[[ii]]$Time,cumulative[[ii]]$ConfirmedCases, t-te)$y
  c2 <- approx(cumulative[[ii]]$Time,cumulative[[ii]]$ConfirmedCases, t)$y
  t <- t[1:(length(t)-n_R0)]
  
  for(i in 1:(length(t)))
  {
    r[i] = sum(c2[i+c(1:n_R0)]-c2[i+c(1:n_R0)-1]) / sum(c1[i+c(1:n_R0)]-c1[i+c(1:n_R0)-1]) 
  }
  f <- !(r==Inf | is.na(r))
  R0[[ii]] <- list("Time" = t[f], "R0" = r[f])
}


plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),")"), ylab="Basic Reproductive Number R0", xlim = c(t_0,max(res$Time)), ylim = c(0,6))

for(ii in 1:8)
{
  data <- data.frame(x=R0[[ii]]$Time, y=R0[[ii]]$R0)
  mean <- lowess(data,f=0.15)
  sd <- data.frame(x=data$x,y=(mean$y - data$y)^2)
  sd$y[2:(length(sd$x)-1)] <- (sd$y[1:(length(sd$x)-2)] + sd$y[2:(length(sd$x)-1)] + sd$y[3:(length(sd$x))])/2
  sd$y[c(1,length(sd$x))]<- sd$y[c(1,length(sd$x))] *3/2
  sd <- lowess(sd, f=0.15)
  sd$y <- abs(sd$y)^0.5*qt(0.975,length(sd$y))
  
  #polygon(c(mean$x,rev(mean$x)),c(mean$y-sd$y,rev(mean$y+sd$y)),lty=0,col="gray95")
  #lines(mean$x,mean$y-sd$y,lty=2)
  #lines(mean$x,mean$y+sd$y,lty=2)
  
  points(data, col=ii)
  lines(mean,lty=2,col=ii)
  
}
lines(c(-1e9,1e9), c(1,1),lty=3,col=3)

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topright", legend=c(ages,"all ages","R0 (3 days)", "trend", "steady state"), col=c(1:8,1,1,3), bty="n", pch=c(rep(-1,8),1,-1,-1), lty=c(rep(1,8),-1,2,3) )

box()

if (plot_out != 2) title("Basic Reproductive Number R0 COVID-19 in Germany", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()

################################
# Extrapolated Basic Reproductive Number R0
################################

if (plot_out == 2) png("Extrapolated_Basic_Reproductive_Number_R0.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)

n_R0 <- 3
R0_Extrapolated <- list()

for(ii in 1:8)
{
  t <- cumulative_extrapolated[[ii]]$Time[cumulative_extrapolated[[ii]]$Time-te>0]
  r=0
  t <- t[-1]
  c1 <- approx(cumulative_extrapolated[[ii]]$Time,cumulative_extrapolated[[ii]]$ConfirmedCases, t-te)$y
  c2 <- approx(cumulative_extrapolated[[ii]]$Time,cumulative_extrapolated[[ii]]$ConfirmedCases, t)$y
  t <- t[1:(length(t)-n_R0)]
  
  for(i in 1:(length(t)))
  {
    r[i] = sum(c2[i+c(1:n_R0)]-c2[i+c(1:n_R0)-1]) / sum(c1[i+c(1:n_R0)]-c1[i+c(1:n_R0)-1]) 
  }
  f <- !(r==Inf | is.na(r))
  R0_Extrapolated[[ii]] <- list("Time" = t[f], "R0" = r[f])
}


plot(NA,NA, xlab = paste0("Days after (",format(start_date, "%m/%d/%Y"),")"), ylab="Extrapolated Basic Reproductive Number R0", xlim = c(t_0,max(res$Time)), ylim = c(0,6))

for(ii in 1:8)
{
  data <- data.frame(x=R0_Extrapolated[[ii]]$Time, y=R0_Extrapolated[[ii]]$R0)
  mean <- lowess(data,f=0.15)
  sd <- data.frame(x=data$x,y=(mean$y - data$y)^2)
  sd$y[2:(length(sd$x)-1)] <- (sd$y[1:(length(sd$x)-2)] + sd$y[2:(length(sd$x)-1)] + sd$y[3:(length(sd$x))])/2
  sd$y[c(1,length(sd$x))]<- sd$y[c(1,length(sd$x))] *3/2
  sd <- lowess(sd, f=0.15)
  sd$y <- abs(sd$y)^0.5*qt(0.975,length(sd$y))
  
  #polygon(c(mean$x,rev(mean$x)),c(mean$y-sd$y,rev(mean$y+sd$y)),lty=0,col="gray95")
  #lines(mean$x,mean$y-sd$y,lty=2)
  #lines(mean$x,mean$y+sd$y,lty=2)
  
  points(data, col=ii)
  lines(mean,lty=2,col=ii)
  
}
lines(c(-1e9,1e9), c(1,1),lty=3,col=3)

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topright", legend=c(ages,"all ages","R0 (3 days)", "trend", "steady state"), col=c(1:8,1,1,3), bty="n", pch=c(rep(-1,8),1,-1,-1), lty=c(rep(1,8),-1,2,3) )

box()

if (plot_out != 2) title("Extrapolated Basic Reproductive Number R0 COVID-19 in Germany", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()

median(R0[[8]]$R0[(length(R0[[8]]$R0)-107):(length(R0[[8]]$R0)-7)])
median(R0_Extrapolated[[8]]$R0[(length(R0[[8]]$R0)-107):(length(R0[[8]]$R0)-7)])
m1 = median(R0_Extrapolated[[8]]$R0[(length(R0[[8]]$R0)-107):(length(R0[[8]]$R0)-7)])/median(R0[[8]]$R0[(length(R0[[8]]$R0)-107):(length(R0[[8]]$R0)-7)])
mean(R0[[8]]$R0[(length(R0[[8]]$R0)-107):(length(R0[[8]]$R0)-7)])
mean(R0_Extrapolated[[8]]$R0[(length(R0[[8]]$R0)-107):(length(R0[[8]]$R0)-7)])
m2 = mean(R0_Extrapolated[[8]]$R0[(length(R0[[8]]$R0)-107):(length(R0[[8]]$R0)-7)])/mean(R0[[8]]$R0[(length(R0[[8]]$R0)-107):(length(R0[[8]]$R0)-7)])


R0_1 <- R0[[8]]$R0#[order(R0[[8]]$R0)]
R0_2 <- R0_Extrapolated[[8]]$R0#[order(R0_Extrapolated[[8]]$R0)]

fit = lm(R0_2~R0_1+0, weights = 1/(R0_1*R0_2))
#plot((R0_1-mean(R0_1))/sd(R0_1), (R0_2-R0_1)/sd(R0_2-R0_1),ylim = c(-1,1),xlim=c(-1,1))
#plot(R0_1, R0_2, ylim=c(0,5),xlim=c(0,5))
#x<-c(-1,10)
#lines(x,x*fit$coefficients[1], col=1)
#lines(x,x*m1, col=2)
#lines(x,x*m2, col=3)
m1
m2
fit
summary(fit)
summary.aov(fit)

#####################################################################################
#
# Constants II
#   
#####################################################################################

phi_d <- phi_d_mean[8]
phi_dd <- phi_d_drift[8]
# assume median from the first three german lock down weeks (good)
phi_d_1 <- phi_d_mean[8]

# assume crtical from Zunyou and McGoogan (2020)
phi_h <- 2087/44415*phi_d_mean[8]/(1023/44672)

# assume phi_d_2 = phi_h * 95%
phi_d_2 <- phi_h * 0.95

# http://www.gbe-bund.de/gbe10/I?I=838:37792217D 
n_h_max <- 28031 

# assume population of germany from https://de.wikipedia.org/wiki/Deutschland (03/13/20)
n_max <- 83.019213e6


if (plot_out == 1) dev.off()