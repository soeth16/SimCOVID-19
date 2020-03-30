# set current working directory
setwd("~/SimCOVID-19")
# setwd("C:/Users/soeth/Desktop/SimCOVID-19")

# git pull
setwd("COVID-19")
system("git pull")
setwd("..")

library(JuliaCall)
julia <- julia_setup()
diffeqr::diffeq_setup()


# output
#
# grapfical: 0
# pdf: 1
# png: 2
#plot_out <- 0
for (plot_out in c(2:0)) {

if (plot_out == 1) pdf("Plots.pdf")


#####################################################################################
#
# Data
#   
#####################################################################################
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
csse_covid_19_time_series_confirmed <- read_csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
csse_covid_19_time_series_recovered <- read_csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv")
csse_covid_19_time_series_deaths    <- read_csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")

ger_data_confirmed <- t(subset(csse_covid_19_time_series_confirmed, `Country/Region`=="Germany", select =-c(1:4)))
ger_data_recovered <- t(subset(csse_covid_19_time_series_recovered, `Country/Region`=="Germany", select =-c(1:4)))
ger_data_deaths    <- t(subset(csse_covid_19_time_series_deaths,    `Country/Region`=="Germany", select =-c(1:4)))

rnames <- rownames(ger_data_confirmed)
len_ger_data <- length(ger_data_confirmed)
#View(data.frame(day=1:len_ger_data, confirmed=ger_data_confirmed, recovered=ger_data_recovered))



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
       bty="n",
       pch=c(-1,-1),
       lwd=1,
       lty=c(1,2),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation COVID-19 in Germany",
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()

ln_data_confirmed <- data.frame(t = 1:len_ger_data,ln_x=log(ger_data_confirmed[1:len_ger_data]))
ln_data_recovered <- data.frame(t = 1:len_ger_data,ln_x=log(ger_data_recovered[1:len_ger_data]))

tc_0 <- 42
tc_max <- 59
tc2_0 <- 61
tc2_max <- len_ger_data
tr_0 <- 35
tr_max <- len_ger_data
dt_r = 14

ln_res_1a = lm(ln_x~t, ln_data_confirmed[15:35,])
ln_res_2a = lm(ln_x~t, ln_data_confirmed[tc_0:tc_max,])
ln_res_3a = lm(ln_x~t, ln_data_confirmed[tc2_0:tc2_max,])
ln_res_1b = lm(ln_x~t, ln_data_recovered[(14+dt_r):(35+dt_r),])
ln_res_2b = lm(ln_x~t, ln_data_recovered[(tr_0+dt_r):tr_max,])
ln_res_1a
ln_res_2a
ln_res_3a
ln_res_1b
ln_res_2b



if (plot_out == 2) png("Situation-2.png", width = 640, height = 480)
plot(ln_data_confirmed,
     type="p", 
     col = 1, 
     pch = 1,
     xlab=paste("Days after", rnames[1]), 
     ylab=expression(log['e']('Casses')))
lines(ln_data_recovered , type="p", col = 2, pch = 2)
lines(15:35,ln_res_1a$fitted.values, col = 1, lty = 1)
lines(tc_0:tc_max,ln_res_2a$fitted.values, col = 1, lty = 2)
lines(tc2_0:tc2_max,ln_res_3a$fitted.values, col = 1, lty = 3)
lines((14+dt_r):(35+dt_r),ln_res_1b$fitted.values, col = 2, lty = 1)
lines((tr_0+dt_r):tr_max,ln_res_2b$fitted.values, col = 2, lty = 2)
legend("topleft", legend <- c("Confirmed cases", "Recovered cases", "Regression 1st phase", "Regression 2nd phase (used for modeling)", "Regression 3nd phase (used for modeling)"), 
       col=c(1,2,16,16),
       bty="n",
       pch=c(1,2,-1,-1,-1),
       lwd=1,
       lty=c(-1,-1,1,2,3),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation COVID-19 in Germany",
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()



#####################################################################################
#
# Constants
#   
#####################################################################################

# assume incubation time from https://www.ncbi.nlm.nih.gov/pubmed/32150748
te = 5.1
ti = 14
th = 8
thi = 14

# assume k_raw from ln(x)-plot
k_raw = ln_res_2a$coefficients[2]
k_raw2 = ln_res_3a$coefficients[2]

# correct k_raw by the influece of te, ti, th, thi and kh
# x2'   = + k(t) x1(t) / x_max x2(t) - k(t - te) * x1 (t - te) / x_max * x2(t - te)
# x3:6' = + k(t - te) * x1 (t - te) / x_max * x2(t - te) # x3:6 as kumulativ infectet --> x3:6' = k_raw * x3:6 && x1 / x_max = 1 && k(t) == const.
#         + k * x1 / x_max * x2(t - te) = k_raw * x3:6(t) # xn(t0 + t) = xn(t0) * exp(k_raw * t) in exponetial phase
#       = + k * x2(t) * exp(-k_raw * te)
# x2(t) = x3:6(t) * k_raw / (k * exp(-k_raw * te))
# x2' = 

f_k <- function(k) k *(1-exp(-k * te)) - k_raw
k <- uniroot(f_k , c(0.01, 1))$root
f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - k_raw2
k2 <- uniroot(f_k2 , c(0.01, 1))$root


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

R0 = exp(k_raw* te) -1
R0
R0 = exp(k_raw*0.7 * te) -1
R0  
log(1+1)/te


if (plot_out == 2) png("Situation-3.png", width = 640, height = 480)
R0_plot <- data.frame(t = c(tc_0:(len_ger_data)), R0 = 0)
for (i in R0_plot$t) R0_plot$R0[i-R0_plot$t[1]+1] <- exp(lm(ln_x~t, ln_data_confirmed[(i-3):i,])$coefficients[2] * te) - 1
R0_plot$t <- R0_plot$t - R0_plot$t[1] 
plot(R0_plot$t, R0_plot$R0, 
     type = "p", 
     ylim=c(0,6),
     xlab=paste("Days after", rnames[tc_0]),
     ylab="Basic Reproductive Number R0")
lines(lowess(R0_plot),lty=2,col=2)
lines(c(0,length(R0_plot$R0)-1), c(1, 1),lty=3,col=3)
legend("topright", legend <- c("Determined from a period of 3 days.", "Trend of Estimate", "Steady State"), 
       col=c(1,2,3),
       bty="n",
       pch=c(1,-1,-1),
       lwd=1,
       lty=c(-1,2,3),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation COVID-19 in Germany",
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


# assume 0.2% deaths from https://www.lungenaerzte-im-netz.de/krankheiten/covid-19/symptome-krankheitsverlauf/
kd_1 <- 0.002

# assume Hubei / China
kd_2 <- as.numeric(
  subset(csse_covid_19_time_series_deaths, `Province/State` == "Hubei", select = length(csse_covid_19_time_series_confirmed[1,1:len_ger_data])) 
  /
    (
      subset(csse_covid_19_time_series_deaths, `Province/State` == "Hubei", select = length(csse_covid_19_time_series_confirmed[1,1:len_ger_data]))
      +
        subset(csse_covid_19_time_series_recovered, `Province/State` == "Hubei", select = length(csse_covid_19_time_series_confirmed[1,1:len_ger_data]))
    )
)

# assume k_bed from max kd_2 (Hubei / China)
kh <- kd_2

# http://www.gbe-bund.de/gbe10/I?I=838:37792217D 
nh_max <- 28031 

# assume population of germany from https://de.wikipedia.org/wiki/Deutschland (03/13/20)
x_max <- 83.019213e6

#        1  2      3   4   5   6    7       8   9     10    11    12
# p <- c(k, x_max, te, ti, th, thi, nh_max, kh, kd_1, kd_2, tk_2, k2)

f = JuliaCall::julia_eval("function f(du, u, h, p, t)
    # parameter
    
    x_max   = p[2]    # max inahbitants
    te      = p[3]    # exposed / incubation time 
    ti      = p[4]    # time of usual infection
    th      = p[5]    # time until hostspitalisation
    thi     = p[6]    # time for hostspitalisation with ARDS
    nh_max  = p[7]    # max intensive beds for treating ARDS
    kh      = p[8]    # propotion factor for hostspitalisation with ARDS
    kd_1    = p[9]    # propotion factor deaths with hostspitalisation
    kd_2    = p[10]   # propotion factor deaths without hostspitalisation
    tk_2    = p[11]   # time of change from k to k2
    k2      = p[12]   # k after change
    
    #k      = p[1]    # growth rate
    du[7] = 0
    if (t >= tk_2) 
      #u[7] = p[1]*0.65
      u[7] = k2
    else
      u[7] = p[1]
    end
    
    # kd
    if (u[4] > nh_max) 
      kd = kd_1 * nh_max / u[4] + kd_2 * (u[4] - nh_max) / u[4]
    else 
      kd = kd_1
    end
    
    # Suspected
    du[1] = - u[7] * u[1] / x_max * u[2]
    
    # Exposed / Incubating
    du[2] = + u[7] * u[1] / x_max * u[2] - h(p, t - te)[7] * h(p, t - te)[1] / x_max * h(p, t - te)[2] 
    
    # Infected
    du[3] = + h(p, t - te)[7] * h(p, t - te)[1] / x_max * h(p, t - te)[2] - h(p, t - te - ti)[7] * h(p, t - te - ti)[1] / x_max * h(p, t - te - ti)[2] * (1-kh)  - h(p, t - te - th)[7] * h(p, t - te - th)[1] / x_max * h(p, t - te - th)[2] * kh

    # Hostspitalisation 
    du[4] = + h(p, t - te - th)[7] * h(p, t - te - th)[1] / x_max * h(p, t - te - th)[2] * kh - h(p, t - te - th)[7] * h(p, t - te - th - thi)[1] / x_max * h(p, t - te - th - thi)[2] * kh
    
    # Recovered
    du[5] = + h(p, t - te - ti)[7] * h(p, t - te - ti)[1] / x_max * h(p, t - te - ti)[2] * (1-kh) + h(p, t - te - th - thi)[7] * h(p, t - te - th - thi)[1] / x_max * h(p, t - te - th - thi)[2] * (kh - kd)
    
    # Deaths
    du[6] = + h(p, t - te - th - thi)[7] * h(p, t - te - th - thi)[1] / x_max * h(p, t - te - th - thi)[2] * kd
  
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
  tmp <- u_0 * exp(p[1] * t)
  tmp[7] <- p[1]
  return(tmp)
}



#####################################################################################
#
# Modell vs Situation
#   
#####################################################################################

t_0 <- 40
tk_2 <- tc_max - t_0 - te

I <- ger_data_confirmed[t_0] - ger_data_recovered[t_0] - ger_data_deaths[t_0]
E <- I * k_raw/(k * exp(-k*te))
I <- I * 0.95
H <- I * 0.05
R <- ger_data_recovered[t_0]
D <- ger_data_deaths[t_0]
S <-x_max - I - E - H - R - D

u_0 <- c(S, E, I, H, R, D, k)
p <- c(k, x_max, te, ti, th, thi, nh_max, kh, kd_1, kd_2, tk_2, k2)
tspan <- list(0,250)
constant_lags = c(1,1,1,1,1,1,1)
t = 0:1000/1000*(250)

sol = diffeqr::dde.solve('f', u_0, h_0, p=p, tspan, saveat = t, constant_lags=constant_lags, abstol = 1e-8, reltol = 1e-8)

if (plot_out == 2) png("Model_vs_Situation-1.png", width = 640, height = 480)
plot(c(t_0:len_ger_data)-t_0, ger_data_confirmed[c(t_0:len_ger_data),1], 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Confirmed cases")
lines(sol$t, sol$u[,3]+sol$u[,4]+sol$u[,5]+sol$u[,6], lty=2)
legend("topleft", legend <- c("Confirmed cases", "Modeled cases"), 
       col=1,
       bty="n",
       pch=c(1,-1),
       lwd=1,
       lty=c(-1,2),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation vs Model COVID-19 in Germany",
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


if (plot_out == 2) png("Model_vs_Situation-2.png", width = 640, height = 480)
plot(c(t_0:len_ger_data)-t_0, ger_data_confirmed[c(t_0:len_ger_data),1]/x_max*100, 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases [%]")
lines(sol$t, (sol$u[,3]+sol$u[,4]+sol$u[,5]+sol$u[,6])/x_max*100, lty=2)
legend("topleft", legend <- c("Confirmed cases", "Modeled cases"), 
       col=1,
       bty="n",
       pch=c(1,-1),
       lwd=1,
       lty=c(-1,2),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation vs Model COVID-19 in Germany",
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


#####################################################################################
#
# Forecast
#   
#####################################################################################

t_0 <- 40
tk_2 <- tc_max - t_0 - te

I <- ger_data_confirmed[t_0] - ger_data_recovered[t_0] - ger_data_deaths[t_0]
E <- I * k_raw/(k * exp(-k*te))
I <- I * 0.95
H <- I * 0.05
R <- ger_data_recovered[t_0]
D <- ger_data_deaths[t_0]
S <- x_max - I - E - H - R - D

u_0 <- c(S, E, I, H, R, D, k)
p <- c(k, x_max, te, ti, th, thi, nh_max, kh, kd_1, kd_2, tk_2, k2)
tspan <- list(0,250)
constant_lags = c(1,1,1,1,1,1,1)
t = 0:1000/1000*(250)


sol = diffeqr::dde.solve('f', u_0, h_0, p=p, tspan, saveat = t, constant_lags=constant_lags, abstol = 1e-8, reltol = 1e-8)

if (plot_out == 2) png("Forecast-1.png", width = 640, height = 480)
plot(sol$t,sol$u[,1],
     type="l", 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases", 
     ylim=c(0,x_max), 
     xlim=c(0,250), 
     lty=1, col=1)
lines(sol$t, sol$u[,2], lty=2, col=2)
lines(sol$t, sol$u[,3], lty=3, col=3)
lines(sol$t, sol$u[,4], lty=4, col=4)
lines(sol$t, sol$u[,5], lty=5, col=5)
lines(sol$t, sol$u[,6], lty=6, col=6)
lines(sol$t, (sol$u[,2]+sol$u[,3]+sol$u[,4]+sol$u[,5]+sol$u[,6]), lty=7, col=16)
legend("topright", legend <- c("Susceptible (noninfected)", "Exposed (incubation time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Total Infected"),       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()





if (plot_out == 2) png("Forecast-2.png", width = 640, height = 480)
plot(sol$t,sol$u[,1]/x_max*100,
     type="l", 
     xlab=paste("Days after", rnames[t_0]),
     ylab="Cases [%]",
     ylim=c(0,100),
     xlim=c(0,250), 
     lty=1, col=1)
lines(sol$t, sol$u[,2]/x_max*100, lty=2, col=2)
lines(sol$t, sol$u[,3]/x_max*100, lty=3, col=3)
lines(sol$t, sol$u[,4]/x_max*100, lty=4, col=4)
lines(sol$t, sol$u[,5]/x_max*100, lty=5, col=5)
lines(sol$t, sol$u[,6]/x_max*100, lty=6, col=6)
lines(sol$t, (sol$u[,2]+sol$u[,3]+sol$u[,4]+sol$u[,5]+sol$u[,6])/x_max*100, lty=7, col=16)
legend("topright", legend <- c("Susceptible (noninfected)", "Exposed (incubation time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Total Infected"), 
       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


max(sol$u[,4])
max(sol$u[,4])/x_max*100
max(sol$u[,4])/nh_max



if (plot_out == 2) png("Forecast-ARDS-1.png", width = 640, height = 480)
plot(sol$t, sol$u[,4], 
     type="l",
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases", 
     ylim=c(0,max(sol$u[,6])), 
     xlim=c(0,100), 
     lty=1, col=1)
lines(sol$t, sol$u[,6], type="l", lty=2, col=2)
par(new=T)
plot(sol$t, sol$u[,4] / nh_max * 100, 
     type="l",
     xlab="", 
     ylab="",
     axes = F,
     xlim=c(0,100), 
     lty=3, col=3)
axis(4, pretty(sol$u[,4] / nh_max * 100,5))
mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths"), 
       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()

if (plot_out == 2) png("Forecast-ARDS-2.png", width = 640, height = 480)
plot(sol$t, sol$u[,4]/x_max*100, 
     type="l",
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases [%]", 
     ylim=c(0,max(sol$u[,6])/x_max*100), 
     xlim=c(0,100), 
     lty=1, col=1)
lines(sol$t, sol$u[,6]/x_max*100, type="l", lty=2, col=2)
par(new=T)
plot(sol$t, sol$u[,4] / nh_max * 100, 
     type="l",
     xlab="", 
     ylab="",
     axes = F,
     xlim=c(0,100), 
     lty=3, col=3)
axis(4, pretty(sol$u[,4] / nh_max * 100,5))
mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths"), 
       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()





#####################################################################################
#
# week bevore 03/23/2020
#   
#####################################################################################

t_0 <- 40
tk_2 <- tc_max - t_0 - te

I <- ger_data_confirmed[t_0] - ger_data_recovered[t_0] - ger_data_deaths[t_0]
E <- I * k_raw/(k * exp(-k*te))
I <- I * 0.95
H <- I * 0.05
R <- ger_data_recovered[t_0]
D <- ger_data_deaths[t_0]
S <- x_max - I - E - H - R - D

u_0 <- c(S, E, I, H, R, D, k)
p <- c(k, x_max, te, ti, th, thi, nh_max, kh, kd_1, kd_2, tk_2, k)
tspan <- list(0,250)
constant_lags = c(1,1,1,1,1,1,1)
t = 0:1000/1000*(250)


sol = diffeqr::dde.solve('f', u_0, h_0, p=p, tspan, saveat = t, constant_lags=constant_lags, abstol = 1e-8, reltol = 1e-8)

if (plot_out == 2) png("Forecast-1-old.png", width = 640, height = 480)
plot(sol$t,sol$u[,1],
     type="l", 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases", 
     ylim=c(0,x_max), 
     xlim=c(0,250), 
     lty=1, col=1)
lines(sol$t, sol$u[,2], lty=2, col=2)
lines(sol$t, sol$u[,3], lty=3, col=3)
lines(sol$t, sol$u[,4], lty=4, col=4)
lines(sol$t, sol$u[,5], lty=5, col=5)
lines(sol$t, sol$u[,6], lty=6, col=6)
lines(sol$t, (sol$u[,2]+sol$u[,3]+sol$u[,4]+sol$u[,5]+sol$u[,6]), lty=7, col=16)
legend("topright", legend <- c("Susceptible (noninfected)", "Exposed (incubation time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Total Infected"),       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany (week bevore 03/23/2020", 
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()





if (plot_out == 2) png("Forecast-2-old.png", width = 640, height = 480)
plot(sol$t,sol$u[,1]/x_max*100,
     type="l", 
     xlab=paste("Days after", rnames[t_0]),
     ylab="Cases [%]",
     ylim=c(0,100),
     xlim=c(0,250), 
     lty=1, col=1)
lines(sol$t, sol$u[,2]/x_max*100, lty=2, col=2)
lines(sol$t, sol$u[,3]/x_max*100, lty=3, col=3)
lines(sol$t, sol$u[,4]/x_max*100, lty=4, col=4)
lines(sol$t, sol$u[,5]/x_max*100, lty=5, col=5)
lines(sol$t, sol$u[,6]/x_max*100, lty=6, col=6)
lines(sol$t, (sol$u[,2]+sol$u[,3]+sol$u[,4]+sol$u[,5]+sol$u[,6])/x_max*100, lty=7, col=16)
legend("topright", legend <- c("Susceptible (noninfected)", "Exposed (incubation time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Total Infected"), 
       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany (week bevore 03/23/2020)", 
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


max(sol$u[,4])
max(sol$u[,4])/x_max*100
max(sol$u[,4])/nh_max


if (plot_out == 2) png("Forecast-ARDS-1-old.png", width = 640, height = 480)
plot(sol$t, sol$u[,4], 
     type="l",
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases", 
     ylim=c(0,max(sol$u[,6])), 
     xlim=c(0,100), 
     lty=1, col=1)
lines(sol$t, sol$u[,6], type="l", lty=2, col=2)
par(new=T)
plot(sol$t, sol$u[,4] / nh_max * 100, 
     type="l",
     xlab="", 
     ylab="",
     axes = F,
     xlim=c(0,100), 
     lty=3, col=3)
axis(4, pretty(sol$u[,4] / nh_max * 100,5))
mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths"), 
       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()

if (plot_out == 2) png("Forecast-ARDS-2-old.png", width = 640, height = 480)
plot(sol$t, sol$u[,4]/x_max*100, 
     type="l",
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases [%]", 
     ylim=c(0,max(sol$u[,6])/x_max*100), 
     xlim=c(0,100), 
     lty=1, col=1)
lines(sol$t, sol$u[,6]/x_max*100, type="l", lty=2, col=2)
par(new=T)
plot(sol$t, sol$u[,4] / nh_max * 100, 
     type="l",
     xlab="", 
     ylab="",
     axes = F,
     xlim=c(0,100), 
     lty=3, col=3)
axis(4, pretty(sol$u[,4] / nh_max * 100,5))
mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths"), 
       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


#####################################################################################
#
# variant calculation
#   
#####################################################################################

t_0 <- 40
tk_2 <- tc_max - t_0 - te

I <- ger_data_confirmed[t_0] - ger_data_recovered[t_0] - ger_data_deaths[t_0]
E <- I * k_raw/(k * exp(-k*te))
I <- I * 0.95
H <- I * 0.05
R <- ger_data_recovered[t_0]
D <- ger_data_deaths[t_0]
S <-x_max - I - E - H - R - D

u_0 <- c(S, E, I, H, R, D, k)
tspan <- list(0,250)
constant_lags = c(1,1,1,1,1,1,1)
# k2 = k to assume nothing is happend
p <- c(k, x_max, te, ti, th, thi, nh_max, kh, kd_1, kd_2, tk_2, k)
t = 0:1000/1000*(250)

sol1 = diffeqr::dde.solve('f', u_0, h_0, p=p, tspan, saveat = t, constant_lags=constant_lags, abstol = 1e-8, reltol = 1e-8)

I <- ger_data_confirmed[t_0] - ger_data_recovered[t_0] - ger_data_deaths[t_0]
E <- I * k_raw/(k * exp(-k*te))
I <- I * 0.95
H <- I * 0.05
R <- ger_data_recovered[t_0]
D <- ger_data_deaths[t_0]
S <-x_max - I - E - H - R - D

u2_0 <- c(S, E, I, H, R, D, k)
tspan <- list(0,250)
constant_lags = c(1,1,1,1,1,1,1)
t = 0:1000/1000*(250)

rbcol = rainbow(11)

if (plot_out == 2) png("Forecast-ARDS-3.png", width = 640, height = 480)


par(mar=c(5,6,6,5)+0.1)

plot(sol1$t, (sol1$u[,4]), 
     type="l",
     axes=F, xlab="", ylab="", 
     ylim=c(-0.2,1)*max((sol1$u[,6])),
     xlim=c(0,150), 
     lty=1, col=1)
lines(sol1$t, (sol1$u[,6]), type="l", lty=2, col=1)

axis(2, pretty(range(sol1$u[,6]),10), col="black",las=1)  ## las=1 makes horizontal labels
mtext("Cases",side=2,line=4.5,las=0)
axis(1,pretty(range(sol1$t),5))
mtext(paste("Time (Days after ", rnames[t_0],")", sep = ""),side=1,col="black",line=2.5)


for (i in c(0:10))
{
  tk_2 <- 7
  R0_2 <- 0.8 + i * 0.2
  f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - log(R0_2 + 1) / te
  k2 <- uniroot(f_k2 , c(0.01, 1))$root
  # parameter vector
  #       1  2      3   4   5   6    7       8   9     10    11    12
  p2 <- c(k, x_max, te, ti, th, thi, nh_max, kh, kd_1, kd_2, tk_2, k2)
  sol2 = diffeqr::dde.solve('f', u_0, h_0, p=p2, tspan, saveat = t, constant_lags=constant_lags, abstol = 1e-8, reltol = 1e-8)
  
  par(new=T)
  plot(sol2$t, (sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max((sol1$u[,6])), xlim=c(0,150))
  lines(sol2$t, (sol2$u[,6]), type="l", lty=2, col=rbcol[i])
  
  par(new=T)
  plot(sol2$t, sol2$u[,7], ylim = c(0.2,1.5), xlim=c(0,150), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
  if(i==0)
  {
    axis(4, pretty(range(sol2$u[,7]),10))
    mtext(expression('Growths Rate '*(D^-1)), side=4,line=3,las=0)
  }
}

par(new=T)
plot(c(0,250), log(c(nh_max, nh_max)), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(0,1)*max(log(sol1$u[,6])), xlim=c(0,150))
for (i in c(0:25)) 
  lines(c(i*7+2,i*7+2), (c(-10, 120)), type="l", lty = 5, col="light gray" )
lines(c(23,23), (c(-10, 25)), type="l", lty = 5, col="grey" )


legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", paste("R0 =", 0.8 + c(0:10)*.2+.8)), 
       col=c(1,1,rbcol[0:8]),
       bty="n",
       lwd=1,
       lty=c(1,2 ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
box()

if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()



rbcol = rainbow(11)

if (plot_out == 2) png("Forecast-ARDS-4.png", width = 640, height = 480)

par(mar=c(5,6,6,5)+0.1)

plot(sol1$t, (sol1$u[,4]), 
     type="l",
     axes=F, xlab="", ylab="", 
     ylim=c(-0.2,1)*max((sol1$u[,6])),
     xlim=c(0,150), 
     lty=1, col=1)
lines(sol1$t, (sol1$u[,6]), type="l", lty=2, col=1)

axis(2, pretty(range(sol1$u[,6]),10), col="black",las=1)  ## las=1 makes horizontal labels
mtext("Cases",side=2,line=4.5,las=0)
axis(1,pretty(range(sol1$t),5))
mtext(paste("Time (Days after ", rnames[t_0],")", sep = ""),side=1,col="black",line=2.5)


for (i in c(0:10))
{
  tk_2 <- 17 + i * 3
  R0_2 <- 1
  f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - log(R0_2 + 1) / te
  k2 <- uniroot(f_k2 , c(0.01, 1))$root
  # parameter vector
  #       1  2      3   4   5   6    7       8   9     10    11    12
  p2 <- c(k, x_max, te, ti, th, thi, nh_max, kh, kd_1, kd_2, tk_2, k2)
  sol2 = diffeqr::dde.solve('f', u_0, h_0, p=p2, tspan, saveat = t, constant_lags=constant_lags, abstol = 1e-8, reltol = 1e-8)
  
  par(new=T)
  plot(sol2$t, (sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max((sol1$u[,6])), xlim=c(0,150))
  lines(sol2$t, (sol2$u[,6]), type="l", lty=2, col=rbcol[i])
  
  par(new=T)
  plot(sol2$t, sol2$u[,7], ylim = c(0.2,1.5), xlim=c(0,150), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
  if(i==0)
  {
    axis(4, pretty(range(sol2$u[,7]),10))
    mtext(expression('Growths Rate '*(D^-1)), side=4,line=3,las=0)
  }
}

par(new=T)
plot(c(0,250), (c(nh_max, nh_max)), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(0,1)*max(log(sol1$u[,6])), xlim=c(0,150))
for (i in c(0:25)) 
  lines(c(i*7+2,i*7+2), (c(-10, 25)), type="l", lty = 5, col="light gray" )
lines(c(23,23), (c(-10, 25)), type="l", lty = 5, col="grey" )

legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", paste("t =", c(0:10)*3+17)), 
       col=c(1,1,rbcol[0:17]),
       bty="n",
       lwd=1,
       lty=c(1,2 ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
box()


if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()




rbcol = rainbow(11)

if (plot_out == 2) png("Forecast-ARDS-5.png", width = 640, height = 480)

par(mar=c(5,6,6,5)+0.1)

sol1$u[log(sol1$u[,4])<=0,4] <- 1 
sol1$u[log(sol1$u[,6])<=0,6] <- 1 

plot(sol1$t, log(sol1$u[,4]), 
     type="l",
     axes=F, xlab="", ylab="", 
     ylim=c(-0.2,1)*max(log(sol1$u[,4])),
     xlim=c(0,150), 
     lty=1, col=1)
lines(sol1$t, log(sol1$u[,6]), type="l", lty=2, col=1)

axis(2, pretty(range(log(sol1$u[,4])),10), col="black",las=1)  ## las=1 makes horizontal labels
mtext(expression(log['e']('Cases')),side=2,line=4,las=0)
axis(1,pretty(range(sol1$t),5))
mtext(paste("Time (Days after ", rnames[t_0],")", sep = ""),side=1,col="black",line=2.5)


for (i in c(0:10))
{
  tk_2 <- 17 + i * 3
  R0_2 <- 1
  f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - log(R0_2 + 1) / te
  k2 <- uniroot(f_k2 , c(0.01, 1))$root
  # parameter vector
  #       1  2      3   4   5   6    7       8   9     10    11    12
  p2 <- c(k, x_max, te, ti, th, thi, nh_max, kh, kd_1, kd_2, tk_2, k2)
  sol2 = diffeqr::dde.solve('f', u_0, h_0, p=p2, tspan, saveat = t, constant_lags=constant_lags, abstol = 1e-8, reltol = 1e-8)
  
  par(new=T)
  sol2$u[log(sol2$u[,4])<=0,4] <- 1 
  sol2$u[log(sol2$u[,6])<=0,6] <- 1 
  plot(sol2$t, log(sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max(log(sol1$u[,4])),xlim=c(0,150))
  lines(sol2$t, log(sol2$u[,6]), type="l", lty=2, col=rbcol[i])
  
  par(new=T)
  plot(sol2$t, sol2$u[,7], ylim = c(0.2,1.5), xlim=c(0,150), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
  if(i==0)
  {
    axis(4, pretty(range(sol2$u[,7]),10))
    mtext(expression('Growths Rate '*(D^-1)), side=4,line=2.5,las=0)
  }
}

par(new=T)
plot(c(0,250), log(c(nh_max, nh_max)), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(0,1)*max(log(sol1$u[,6])), xlim=c(0,150))
for (i in c(0:25)) 
  lines(c(i*7+2,i*7+2), (c(-10, 25)), type="l", lty = 5, col="light gray" )
lines(c(23,23), (c(-10, 25)), type="l", lty = 5, col="grey" )

legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", paste("t =", c(0:10)*3+17)), 
       col=c(1,1,rbcol[0:17]),
       bty="n",
       lwd=1,
       lty=c(1,2 ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
box()


if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()




min(sol1$t[(sol1$u[,4] > nh_max)]) + 16


# the end 
if (plot_out == 1) dev.off()
}

