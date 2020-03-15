# set current working directory
setwd("~/SimCOVID-19")

# output
#
# grapfical: 0
# pdf: 1
# png: 2
plot_out <- 2


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

if (plot_out == 1) pdf("Plots.pdf")

if (plot_out == 2) png("Situation-1.png", width = 640, height = 480)
plot(ger_data_confirmed , type="l", col = 1, lty = 1)
lines(ger_data_recovered , type="l", col = 2, lty = 2)
if (plot_out > 1) dev.off()

lm_data_confirmed <- data.frame(t = 1:len_ger_data,ln_x=log(ger_data_confirmed))
lm_data_recovered <- data.frame(t = 1:len_ger_data,ln_x=log(ger_data_recovered))

t_0 <- 42
dt_r = 14
lm_res_1 = lm(ln_x~t, lm_data_confirmed[15:35,])
lm_res_2 = lm(ln_x~t, lm_data_confirmed[t_0:len_ger_data,])
lm_res_3 = lm(ln_x~t, lm_data_recovered[(14+dt_r):(35+dt_r),])
lm_res_4 = lm(ln_x~t, lm_data_recovered[(t_0+dt_r):len_ger_data,])
lm_res_1
lm_res_2
lm_res_3
lm_res_4


if (plot_out == 2) png("Situation-2.png", width = 640, height = 480)
plot(lm_data_confirmed , type="p", col = 1, pch = 1)
lines(lm_data_recovered , type="p", col = 2, pch = 2)
lines(15:35,lm_res_1$fitted.values, col = 1, lty = 1)
lines(t_0:len_ger_data,lm_res_2$fitted.values, col = 1, lty = 2)
lines((14+dt_r):(35+dt_r),lm_res_3$fitted.values, col = 2, lty = 1)
lines((t_0+dt_r):len_ger_data,lm_res_4$fitted.values, col = 2, lty = 2)
if (plot_out > 1) dev.off()


# assume incubation time from https://www.ncbi.nlm.nih.gov/pubmed/32150748
ti = 5.1

# assume k_raw from ln(x)-plot
k_raw = lm_res_2$coefficients[2]

# example from Kentaro Iwata and Chisato Miyakoshi
# https://www.preprints.org/manuscript/202002.0179/v1
# https://gist.github.com/skoba/abc760104be559881ab7269372bb03ea#file-covid19-py
# ti = 10
# R0 = 3
# k_raw = log(R0 + 1)/ti

# dx/dt = k * x
# 1/x * dx = k * dt
# integrate
# ln(x) - ln(x_0) = k * t
# x / x_0 = exp(k * t)
# x = x_0 exp(k * t)
# x_0 = 1
# R0 = exp(k * ti) -1

R0 = as.numeric(exp(k_raw * ti) - 1)
R0

# correct k_raw by the influece of ti
k <- k_raw + k_raw * exp(-k_raw * ti)

# assume 0.2% deaths from https://www.lungenaerzte-im-netz.de/krankheiten/covid-19/symptome-krankheitsverlauf/
kd_1 <- 0.002

# assume Hubei / China
kd_2 <- as.numeric(
        subset(csse_covid_19_time_series_deaths, `Province/State` == "Hubei", select = length(csse_covid_19_time_series_confirmed[1,])) 
        /
        (
          subset(csse_covid_19_time_series_deaths, `Province/State` == "Hubei", select = length(csse_covid_19_time_series_confirmed[1,]))
          +
          subset(csse_covid_19_time_series_recovered, `Province/State` == "Hubei", select = length(csse_covid_19_time_series_confirmed[1,]))
        )
      )

# assume k_bed from max kd_2 (Hubei / China)
k_bed <- kd_2

# http://www.gbe-bund.de/gbe10/I?I=838:37792217D 
n_bed_max <- 28031 

# assume population of germany from https://de.wikipedia.org/wiki/Deutschland (03/13/20)
x_max <- 83.019213e6

#      1  2      3   4          5      6     7
p <- c(k, x_max, ti, n_bed_max, k_bed, kd_1, kd_2)

f = JuliaCall::julia_eval("function f(du, u, h, p, t)
  # parameter
  k         = p[1]
  x_max     = p[2]
  ti        = p[3]
  n_bed_max = p[4]
  k_bed     = p[5]
  kd_1      = p[6]
  kd_2      = p[7]
  
  # - infection
  du[1] = - k * u[1] / x_max * u[2]
  
  # + infektion - (recovered + deaths)
  du[2] = + k * u[1] / x_max * u[2] - k * h(p, t - ti)[1] / x_max * h(p, t - ti)[2] 
  
  # kd
  if ((u[2] * k_bed) > n_bed_max) 
    kd = kd_1 * n_bed_max / (u[2] * k_bed) + kd_2 * ((u[2] * k_bed) - n_bed_max) / (u[2] * k_bed)
  else 
    kd = kd_1
  end

  # + recovered
  # approximation of history data
  du[3] = + (1 - kd) * k * h(p, t - ti)[1] / x_max * h(p, t - ti)[2]
  
  # + deaths
  du[4] = + kd * k * h(p, t - ti)[1] / x_max * h(p, t - ti)[2]

end")


# dx/dt = k * x
# 1/x * dx = k * dt
# integrate
# ln(x) - ln(x_0) = k * t
# x / x_0 = exp(k * t)
# x = x_0 exp(k * t)
# x_0 = 1
h_0 <- function(p, t)
{
  return(u_0 * exp(p[1] * t))
}

t_0 <- 40
u_0 <- c(x_max - ger_data_confirmed[t_0], ger_data_confirmed[t_0] - ger_data_recovered[t_0] - ger_data_deaths[t_0], ger_data_recovered[t_0], ger_data_deaths[t_0])
tspan <- list(0,250)
constant_lags = c(ti)
t = 0:10000/10000*(250)

sol = diffeqr::dde.solve('f', u_0, h_0, p=p, tspan, saveat = t, constant_lags=constant_lags)

if (plot_out == 2) png("Modell_vs_Situation-1.png", width = 640, height = 480)
plot(c(t_0:len_ger_data)-t_0, ger_data_confirmed[c(t_0:len_ger_data),1], 
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
      sub="Created by Sören Thiering 3/15/20. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


if (plot_out == 2) png("Modell_vs_Situation-2.png", width = 640, height = 480)
plot(c(t_0:len_ger_data)-t_0, ger_data_confirmed[c(t_0:len_ger_data),1]/x_max*100, 
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
      sub="Created by Sören Thiering 3/15/20. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()

t_0 <- 53
u_0 <- c(x_max - ger_data_confirmed[t_0], ger_data_confirmed[t_0] - ger_data_recovered[t_0] - ger_data_deaths[t_0], ger_data_recovered[t_0], ger_data_deaths[t_0])
tspan <- list(0,250)
constant_lags = c(ti)
t = 0:10000/10000*(250)


sol = diffeqr::dde.solve('f', u_0, h_0, p=p, tspan, saveat = t, constant_lags=constant_lags)

if (plot_out == 2) png("Forecast-1.png", width = 640, height = 480)
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
      sub="Created by Sören Thiering 3/15/20. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()

if (plot_out == 2) png("Forecast-2.png", width = 640, height = 480)
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
      sub="Created by Sören Thiering 3/15/20. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()

if (plot_out == 1) dev.off() 

