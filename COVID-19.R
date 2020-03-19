# set current working directory
setwd("~/SimCOVID-19")

# output
#
# grapfical: 0
# pdf: 1
# png: 2
plot_out <- 1


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
plot(ger_data_confirmed, 
     type="l", 
     col = 1, 
     lty = 1,
     xlab=paste("Days after", rnames[1]), 
     ylab="Cases")
lines(ger_data_recovered , type="l", col = 2, lty = 2)
legend("topleft", legend <- c("Confirmed cases", "Recovered cases"), 
       col=c(1,2),
       bg="white",
       pch=c(-1,-1),
       lwd=1,
       lty=c(1,2),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation COVID-19 in Germany",
      sub="Created by Sören Thiering 03/18/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()

ln_data_confirmed <- data.frame(t = 1:len_ger_data,ln_x=log(ger_data_confirmed))
ln_data_recovered <- data.frame(t = 1:len_ger_data,ln_x=log(ger_data_recovered))

tc_0 <- 42
tr_0 <- 35
dt_r = 14

ln_res_1 = lm(ln_x~t, ln_data_confirmed[15:35,])
ln_res_2 = lm(ln_x~t, ln_data_confirmed[tc_0:len_ger_data,])
ln_res_3 = lm(ln_x~t, ln_data_recovered[(14+dt_r):(35+dt_r),])
ln_res_4 = lm(ln_x~t, ln_data_recovered[(tr_0+dt_r):len_ger_data,])
ln_res_1
ln_res_2
ln_res_3
ln_res_4


if (plot_out == 2) png("Situation-2.png", width = 640, height = 480)
plot(ln_data_confirmed,
     type="p", 
     col = 1, 
     pch = 1,
     xlab=paste("Days after", rnames[1]), 
     ylab="ln( cases )")
lines(ln_data_recovered , type="p", col = 2, pch = 2)
lines(15:35,ln_res_1$fitted.values, col = 1, lty = 1)
lines(tc_0:len_ger_data,ln_res_2$fitted.values, col = 1, lty = 2)
lines((14+dt_r):(35+dt_r),ln_res_3$fitted.values, col = 2, lty = 1)
lines((tr_0+dt_r):len_ger_data,ln_res_4$fitted.values, col = 2, lty = 2)
legend("topleft", legend <- c("Confirmed cases", "Recovered cases", "Regression 1st phase", "Regression 2nd phase (used for modeling)"), 
       col=c(1,2,16,16),
       bg="white",
       pch=c(1,2,-1,-1),
       lwd=1,
       lty=c(-1,-1,1,2),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation COVID-19 in Germany",
      sub="Created by Sören Thiering 03/18/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


# assume incubation time from https://www.ncbi.nlm.nih.gov/pubmed/32150748
te = 5.1
ti = 14
th = 8
thi = 14

# assume k_raw from ln(x)-plot
k_raw = ln_res_2$coefficients[2]

# example from Kentaro Iwata and Chisato Miyakoshi
# https://www.preprints.org/manuscript/202002.0179/v1
# https://gist.github.com/skoba/abc760104be559881ab7269372bb03ea#file-covid19-py
# ti = 10
# R0 = 3
# k_raw = log(R0 + 1)/te

# dx/dt = k * x
# 1/x * dx = k * dt
# integrate
# ln(x) - ln(x_0) = k * t
# x / x_0 = exp(k * t)
# x = x_0 exp(k * t)
# x_0 = 1
# R0 = exp(k * te) -1

if (plot_out == 2) png("Situation-3.png", width = 640, height = 480)
R0_plot <- data.frame(t = c(tc_0:(len_ger_data)), R0 = 0)
for (i in R0_plot$t) R0_plot$R0[i-R0_plot$t[1]-1] <- exp(lm(ln_x~t, ln_data_confirmed[(i-10):i,])$coefficients[2] * te) - 1
R0_plot$t <- R0_plot$t - R0_plot$t[1] 
plot(R0_plot$t, R0_plot$R0, 
     type = "p", 
     ylim=c(0,6),
     xlab=paste("Days after", rnames[40]),
     ylab="Basic Reproductive Number R0")
lines(lowess(R0_plot),lty=2,col=2)
lines(c(0,length(R0_plot$R0)-1), c(1, 1),lty=3,col=3)
legend("topright", legend <- c("Determined from a period of 10 days.", "Trend of Estimate", "Steady State"), 
       col=c(1,2,3),
       bg="white",
       pch=c(1,-1,-1),
       lwd=1,
       lty=c(-1,2,3),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation COVID-19 in Germany",
                         sub="Created by Sören Thiering 03/18/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


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
kh <- kd_2

# http://www.gbe-bund.de/gbe10/I?I=838:37792217D 
nh_max <- 28031 

# assume population of germany from https://de.wikipedia.org/wiki/Deutschland (03/13/20)
x_max <- 83.019213e6

# correct k_raw by the influece of te, ti, th, thi and kh
k <- k_raw + k_raw * exp(-k_raw * (ti * (1-kh) + (th+thi) * kh - te)) 

# parameter vector
#      1  2      3   4   5   6    7       8   9     10
p <- c(k, x_max, te, ti, th, thi, nh_max, kh, kd_1, kd_2)

f = JuliaCall::julia_eval("function f(du, u, h, p, t)
  # parameter
  k       = p[1]
  x_max   = p[2]
  te      = p[3]
  ti      = p[4]
  th      = p[5]
  thi     = p[6]
  nh_max  = p[7]
  kh      = p[8]
  kd_1    = p[9]
  kd_2    = p[10]
  
  # Suspected
  du[1] = - k * u[1] / x_max * u[2]
  
  # Exposed / Incubating
  du[2] = + k * u[1] / x_max * u[2] - k * h(p, t - te)[1] / x_max * h(p, t - te)[2] 
  
  # kd
  if (u[4] > nh_max) 
    kd = kd_1 * nh_max / u[4] + kd_2 * (u[4] - nh_max) / u[4]
  else 
    kd = kd_1
  end

  # Infected
  du[3] = + k * h(p, t - te)[1] / x_max * h(p, t - te)[2] - k * h(p, t - te - ti)[1] / x_max * h(p, t - te - ti)[2] * (1-kh)  - k * h(p, t - te - th)[1] / x_max * h(p, t - te - th)[2] * kh

  # Hostspitalisation 
  du[4] = + k * h(p, t - te - th)[1] / x_max * h(p, t - te - th)[2] * kh - k * h(p, t - te - th - thi)[1] / x_max * h(p, t - te - th - thi)[2] * kh
  
  # Recovered
  du[5] = + k * h(p, t - te - ti)[1] / x_max * h(p, t - te - ti)[2] * (1-kh) + k * h(p, t - te - th - thi)[1] / x_max * h(p, t - te - th - thi)[2] * (kh - kd)

  # Deaths
  du[6] = + k * h(p, t - te - th - thi)[1] / x_max * h(p, t - te - th - thi)[2] * kd


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
S <-x_max - ger_data_confirmed[t_0]
I <- ger_data_confirmed[t_0] - ger_data_recovered[t_0] - ger_data_deaths[t_0]
E <- I * exp(k*5.1)
I <- I * 0.95
H <- I * 0.05
R <- ger_data_recovered[t_0]
D <- ger_data_deaths[t_0]
u_0 <- c(S, E, I, H, R, D)
tspan <- list(0,250)
constant_lags = c(ti+8)
t = 0:10000/10000*(250)

sol = diffeqr::dde.solve('f', u_0, h_0, p=p, tspan, saveat = t, constant_lags=constant_lags)

if (plot_out == 2) png("Model_vs_Situation-1.png", width = 640, height = 480)
plot(c(t_0:len_ger_data)-t_0, ger_data_confirmed[c(t_0:len_ger_data),1], 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Confirmed cases")
lines(sol$t, sol$u[,3]+sol$u[,4]+sol$u[,5]+sol$u[,6], lty=2)
legend("topleft", legend <- c("Confirmed cases", "Modeled cases"), 
       col=1,
       bg="white",
       pch=c(1,-1),
       lwd=1,
       lty=c(-1,2),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation vs Model COVID-19 in Germany",
      sub="Created by Sören Thiering 03/18/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


if (plot_out == 2) png("Model_vs_Situation-2.png", width = 640, height = 480)
plot(c(t_0:len_ger_data)-t_0, ger_data_confirmed[c(t_0:len_ger_data),1]/x_max*100, 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases [%]")
lines(sol$t, (sol$u[,3]+sol$u[,4]+sol$u[,5]+sol$u[,6])/x_max*100, lty=2)
legend("topleft", legend <- c("Confirmed cases", "Modeled cases"), 
       col=1,
       bg="white",
       pch=c(1,-1),
       lwd=1,
       lty=c(-1,2),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation vs Model COVID-19 in Germany",
      sub="Created by Sören Thiering 03/18/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()

t_0 <- 55
S <-x_max - ger_data_confirmed[t_0]
I <- ger_data_confirmed[t_0] - ger_data_recovered[t_0] - ger_data_deaths[t_0]
E <- I * exp(k*5.1)
I <- I * 0.95
H <- I * 0.05
R <- ger_data_recovered[t_0]
D <- ger_data_deaths[t_0]
u_0 <- c(S, E, I, H, R, D)
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
lines(sol$t-ti, sol$u[,5], lty=5, col=5)
lines(sol$t-ti, sol$u[,6], lty=6, col=6)
lines(sol$t-ti, (sol$u[,2]+sol$u[,3]+sol$u[,4]+sol$u[,5]+sol$u[,6]), lty=7, col=16)
legend("topright", legend <- c("Susceptible (noninfected)", "Exposed (incubation time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Total Infected"),       col=c(1:6,16),
       bg="white",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany", 
      sub="Created by Sören Thiering 03/18/2020. Email: soeren.thiering@hs-anhalt.de")
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
lines(sol$t-ti, sol$u[,5]/x_max*100, lty=5, col=5)
lines(sol$t-ti, sol$u[,6]/x_max*100, lty=6, col=6)
lines(sol$t-ti, (sol$u[,2]+sol$u[,3]+sol$u[,4]+sol$u[,5]+sol$u[,6])/x_max*100, lty=7, col=16)
legend("topright", legend <- c("Susceptible (noninfected)", "Exposed (incubation time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Total Infected"), 
       col=c(1:6,16),
       bg="white",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany", 
      sub="Created by Sören Thiering 03/18/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


max(sol$u[,4])
max(sol$u[,4])/x_max*100
max(sol$u[,4])/nh_max


if (plot_out == 2) png("Forecast-ARDS-1.png", width = 640, height = 480)
plot(sol$t-ti, sol$u[,4], 
     type="l",
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases", 
     ylim=c(0,max(sol$u[,6])), 
     xlim=c(0,100), 
     lty=1, col=1)
lines(sol$t-ti, sol$u[,6], type="l", lty=2, col=2)
legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths"), 
       col=c(1:6,16),
       bg="white",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany", 
      sub="Created by Sören Thiering 03/18/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()

if (plot_out == 2) png("Forecast-ARDS-2.png", width = 640, height = 480)
plot(sol$t-ti, sol$u[,4]/x_max*100, 
     type="l",
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases [%]", 
     ylim=c(0,max(sol$u[,6])/x_max*100), 
     xlim=c(0,100), 
     lty=1, col=1)
lines(sol$t-ti, sol$u[,6]/x_max*100, type="l", lty=2, col=2)
legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths"), 
       col=c(1:6,16),
       bg="white",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany", 
      sub="Created by Sören Thiering 03/18/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()

if (plot_out == 1) dev.off() 
