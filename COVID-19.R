# set current working directory
setwd("~/SimCOVID-19")

library(JuliaCall)
julia <- julia_setup()
diffeqr::diffeq_setup()

# Data 
# https://github.com/CSSEGISandData/COVID-19/
#
# download:
# git clone https://github.com/CSSEGISandData/COVID-19/ COVID-19
#
# up date:
# cd ./COVID-19
# git pull
#
library(readr)
csse_covid_19_time_series_confirmed <- read_csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv")
csse_covid_19_time_series_recovered <- read_csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Recovered.csv")
csse_covid_19_time_series_deaths    <- read_csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Deaths.csv")

ger_data_confirmed <- t(subset(csse_covid_19_time_series_confirmed, `Country/Region`=="Germany", select =-c(1:4)))
ger_data_recovered <- t(subset(csse_covid_19_time_series_recovered, `Country/Region`=="Germany", select =-c(1:4)))
ger_data_deaths    <- t(subset(csse_covid_19_time_series_deaths,    `Country/Region`=="Germany", select =-c(1:4)))

rnames <- rownames(ger_data_confirmed)
len_ger_data <- length(ger_data_confirmed)

#View(data.frame(day=1:len_ger_data, confirmed=ger_data_confirmed, recovered=ger_data_recovered))

plot(ger_data_confirmed , type="l", col = 1, lty = 1)
lines(ger_data_recovered , type="l", col = 2, lty = 2)


lm_data_confirmed <- data.frame(t = 1:len_ger_data,ln_x=log(ger_data_confirmed))
lm_data_recovered <- data.frame(t = 1:len_ger_data,ln_x=log(ger_data_recovered))

plot(lm_data_confirmed , type="p", col = 1, pch = 1)
lines(lm_data_recovered , type="p", col = 2, pch = 2)

t_0 <- 40

lm_res_1 = lm(ln_x~t, lm_data_confirmed[15:35,])
lines(15:35,lm_res_1$fitted.values, col = 1, lty = 1)
lm_res_2 = lm(ln_x~t, lm_data_confirmed[t_0:len_ger_data,])
lines(t_0:len_ger_data,lm_res_2$fitted.values, col = 1, lty = 2)
lm_res_1
lm_res_2
dt_r = 14
lm_res_3 = lm(ln_x~t, lm_data_recovered[(14+dt_r):(35+dt_r),])
lines((14+dt_r):(35+dt_r),lm_res_3$fitted.values, col = 2, lty = 1)
lm_res_4 = lm(ln_x~t, lm_data_recovered[(t_0+dt_r):len_ger_data,])
lines((t_0+dt_r):len_ger_data,lm_res_4$fitted.values, col = 2, lty = 2)


# assume incubation time from https://www.ncbi.nlm.nih.gov/pubmed/32150748
ti = 5.1
# assume k from ln(x)-plot and correct it by the influece of ti
k <- lm_res_2$coefficients[2] + lm_res_2$coefficients[2] * exp(-lm_res_2$coefficients[2] * ti)
# assume 0.2% deaths from https://www.lungenaerzte-im-netz.de/krankheiten/covid-19/symptome-krankheitsverlauf/
kd <- 0.002
# assume population of germany from https://de.wikipedia.org/wiki/Deutschland (03/13/20)
x_max <- 83.019213e6

#dx/dt = k * x
#1/x * dx = k * dt
#integrate
# ln(x) - ln(x_0) = k * t
# x / x_0 = exp(k * t)
# x = x_0 exp(k * t)


f <- function(x, p, t) 
{
  # - infection
  dx_1 = - k * x[1] / x_max * x[2]
  
  # + infektion - (recovered + deaths)
  dx_2 = + k * x[1] / x_max * x[2] - k * exp(-k * ti) * x[2]
  
  # + recovered
  dx_3 = + (1 - kd) * k * exp(-k * ti) * x[2]
  
  # + deaths
  dx_4 = + kd * k * exp(-k * ti) * x[2]
  
  return(c(dx_1, dx_2, dx_3, dx_4))
}


t_0 <- 40
x_0 <- c(x_max - ger_data_confirmed[t_0], ger_data_confirmed[t_0] - ger_data_recovered[t_0] - ger_data_deaths[t_0], ger_data_recovered[t_0], ger_data_deaths[t_0])
tspan <- list(0,250)
t = 0:10000/10000*(250)

sol = diffeqr::ode.solve(f, x_0, tspan, saveat = t)

plot(c(t_0:50)-t_0, ger_data_confirmed[c(t_0:50),1], 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Confirmed cases")
lines(sol$t, sol$u[,2]+sol$u[,3]+sol$u[,4], lty=2)
legend("topleft", legend <- c("Confirmed cases", "Modeled cases"), 
       col=1,
       bg="white",
       pch=c(1,-1),
       lwd=1,
       lty=c(-1,2),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
title("Situation COVID-19 in Germany",
      sub="Created by Sören Thiering 3/13/20. Email: soeren.thiering@hs-anhalt.de")


plot(c(t_0:50)-t_0, ger_data_confirmed[c(t_0:50),1]/x_max*100, 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases [%]")
lines(sol$t, (sol$u[,2]+sol$u[,3]+sol$u[,4])/x_max*100, lty=2)
legend("topleft", legend <- c("Confirmed cases", "Modeled cases"), 
       col=1,
       bg="white",
       pch=c(1,-1),
       lwd=1,
       lty=c(-1,2),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
title("Situation COVID-19 in Germany",
      sub="Created by Sören Thiering 3/13/20. Email: soeren.thiering@hs-anhalt.de")


t_0 <- 51
x_0 <- c(x_max - ger_data_confirmed[t_0], ger_data_confirmed[t_0] - ger_data_recovered[t_0] - ger_data_deaths[t_0], ger_data_recovered[t_0], ger_data_deaths[t_0])
tspan <- list(0,250)
t = 0:10000/10000*(250)

sol = diffeqr::ode.solve(f, x_0, tspan, saveat = t)

plot(sol$t-ti,sol$u[,1],
     type="l", 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases", 
     ylim=c(0,x_max), 
     xlim=c(0,100), 
     lty=1, col=1)
lines(sol$t-ti, sol$u[,2], lty=2, col=2)
lines(sol$t-ti, sol$u[,3], lty=3, col=3)
lines(sol$t-ti, sol$u[,4], lty=4, col=4)
lines(sol$t-ti, (sol$u[,2]+sol$u[,3]+sol$u[,4]), lty=5, col=16)
legend("topright", legend <- c("noninfected","incubation","recovered","deaths","total infected"), 
       col=c(1:4,16),
       bg="white",
       lwd=1,
       lty=c(1:5),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
title("Forecast COVID-19 in Germany", 
      sub="Created by Sören Thiering 3/13/20. Email: soeren.thiering@hs-anhalt.de")


plot(sol$t-ti,sol$u[,1]/x_max*100,
     type="l", 
     xlab=paste("Days after", rnames[t_0]),
     ylab="Cases [%]",
     ylim=c(0,100),
     xlim=c(0,100), 
     lty=1, col=1)
lines(sol$t-ti, sol$u[,2]/x_max*100, lty=2, col=2)
lines(sol$t-ti, sol$u[,3]/x_max*100, lty=3, col=3)
lines(sol$t-ti, sol$u[,4]/x_max*100, lty=4, col=4)
lines(sol$t-ti, (sol$u[,2]+sol$u[,3]+sol$u[,4])/x_max*100, lty=5, col=16)
legend("topright", legend <- c("noninfected","incubation","recovered","deaths","total infected"), 
       col=c(1:4,16),
       bg="white",
       lwd=1,
       lty=c(1:5),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
title("Forecast COVID-19 in Germany", 
      sub="Created by Sören Thiering 3/13/20. Email: soeren.thiering@hs-anhalt.de")
 
