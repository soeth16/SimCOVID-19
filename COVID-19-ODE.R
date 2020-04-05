# set current working directory
setwd("~/SimCOVID-19")
# setwd("C:/Users/soeth/Desktop/SimCOVID-19")

# git pull
setwd("COVID-19")
system("git pull")
setwd("..")

library(JuliaCall)
JuliaCall::julia_install_package_if_needed("DifferentialEquations")
JuliaCall::julia_library("DifferentialEquations")
JuliaCall::julia_install_package_if_needed("DiffEqParamEstim")
JuliaCall::julia_library("DiffEqParamEstim")
JuliaCall::julia_install_package_if_needed("Optim")
JuliaCall::julia_library("Optim")
JuliaCall::julia_install_package_if_needed("Distributions")
JuliaCall::julia_library("Distributions")
JuliaCall::julia_install_package_if_needed("BlackBoxOptim")
JuliaCall::julia_library("BlackBoxOptim")


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

#####################################################################################
#
# Modell
#   
#####################################################################################
#        1  2     3
# p <- c(k, tk_2, k2)

JuliaCall::julia_assign("x_max",x_max) 
JuliaCall::julia_assign("te", te) 
JuliaCall::julia_assign("ti", ti) 
JuliaCall::julia_assign("th", th)
JuliaCall::julia_assign("thi", thi)
JuliaCall::julia_assign("nh_max", nh_max)
JuliaCall::julia_assign("kh", kh)
JuliaCall::julia_assign("kd_1", kd_1)
JuliaCall::julia_assign("kd_2", kd_2)

f = JuliaCall::julia_eval("function f(du, u, p, t)
    # parameter
    # k     = p[1]  # growth rate
    # tk_2  = p[2]  # time of change from k to k2
    # k2    = p[3]  # k after change
    du[8] = 0
    if (t >= p[2]) 
      k = p[3]
    else
      k = p[1]
    end
    
    # kd
    if (u[4] > nh_max) 
      kd = kd_1 * nh_max / u[4] + kd_2 * (u[4] - nh_max) / u[4]
    else 
      kd = kd_1
    end
    
    # Susceptibles
    du[7] = - k * u[7] / x_max * u[6]
    
    # Exposed / Incubating
    du[6] = + k * u[7] / x_max * u[6] - 1/te * u[6] 
    
    # Infected
    du[5] = + 1/te * u[6] - 1/ti * u[5] * (1-kh)  - 1/th * u[5] * kh

    # Hostspitalisation 
    du[4] = + 1/th * u[5] * kh - 1/thi * u[4]
    
    # Recovered
    du[3] = + 1/ti * u[5] * (1-kh) + 1/thi * u[4] * (kh - kd)
    
    # Deaths
    du[2] = + 1/thi * u[4] * kd
  
    # Confirmed
    du[1] = du[5] + du[4] + du[3] + du[2]
    
  end")



#####################################################################################
#
# Modell parameter fit
#   
#####################################################################################

t_0 <- 40
tk_2 <- tc_max - t_0 - te


t = c(t_0:len_ger_data)
C = ger_data_confirmed
I <- ger_data_confirmed - ger_data_recovered - ger_data_deaths
E <- I * k_raw/(k * exp(-k*te))
I <- I * 0.95
H <- I * 0.05
R <- ger_data_recovered
D <- ger_data_deaths
S <-x_max - I - E - H - R - D

data_df <- data.frame(C[t], D[t], R[t])

JuliaCall::julia_assign("u0", c(C[t_0], D[t_0], R[t_0], H[t_0], I[t_0], E[t_0], S[t_0], k))
JuliaCall::julia_assign("tspan", c(0,250))
JuliaCall::julia_assign("p", c(k, tk_2, k))
JuliaCall::julia_assign("saveat", c(0:1000/1000*250))
JuliaCall::julia_eval("prob = ODEProblem(f,u0,tspan,p)")

data <- array(dim=c(3,length(t),100))
for (i in 1:3)
  for (j in 1:length(t))
    data[i,j,] <- rnorm(100,data_df[j,i],1000)
JuliaCall::julia_assign("t", t-t_0)
JuliaCall::julia_assign("data", data)
JuliaCall::julia_eval("distributions = [fit_mle(Normal,data[i,j,:]) for i in 1:3, j in 1:length(t)]")
JuliaCall::julia_eval("obj = build_loss_objective(prob,Tsit5(),LogLikeLoss(t,distributions),maxiters=10000,verbose=false); nothing")

# L2Loss
#data <- array(dim=c(3,length(t),1))
#for ( i in 1:3)
#  for (j in 1:length(t))
#    data[i,j,] <- data_df[j,i]
#JuliaCall::julia_assign("t", t-t_0)
#JuliaCall::julia_assign("data", data)
#JuliaCall::julia_eval("obj = build_loss_objective(prob,Tsit5(),L2Loss(t,data), maxiters=10000,verbose=false); nothing")


JuliaCall::julia_eval("bound1 = Tuple{Float64,Float64}[(0.1,1),(13,15),(0.1,0.5)]")
JuliaCall::julia_eval("res = bboptimize(obj;SearchRange = bound1, MaxSteps = 1e4)")

p2 <- JuliaCall::julia_eval("p = best_candidate(res)")

#####################################################################################
#
# Output
#   
#####################################################################################
# grapfical: 0
# pdf: 1
# png: 2
plot_out <- 0

#####################################################################################
#
# Modell vs Situation
#   
#####################################################################################

JuliaCall::julia_assign("p", p2)
JuliaCall::julia_eval("prob = ODEProblem(f,u0,tspan,p)")
JuliaCall::julia_eval("sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-9, maxiters = 100000, saveat=saveat); nothing")
sol = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))


if (plot_out == 2) png("Model_vs_Situation-1.png", width = 640, height = 480)
plot(c(t_0:len_ger_data)-t_0, ger_data_confirmed[c(t_0:len_ger_data),1], 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Confirmed cases")
lines(sol$t, sol$u[,1], lty=2)
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
     ylab="Cases (%)")
lines(sol$t, (sol$u[,1])/x_max*100, lty=2)
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

if (plot_out == 2) png("Forecast-1.png", width = 640, height = 480)

plot(sol$t,sol$u[,7],
     type="l", 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases", 
     ylim=c(0,x_max), 
     xlim=c(0,250), 
     lty=1, col=1)
lines(sol$t, sol$u[,6], lty=2, col=2)
lines(sol$t, sol$u[,5], lty=3, col=3)
lines(sol$t, sol$u[,4], lty=4, col=4)
lines(sol$t, sol$u[,3], lty=5, col=5)
lines(sol$t, sol$u[,2], lty=6, col=6)
lines(sol$t, sol$u[,1], lty=7, col=16)
legend("topright", legend <- c("Susceptibles (noninfected)", "Exposed (incubation time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Total Infected"),       col=c(1:6,16),
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
plot(sol$t,sol$u[,7]/x_max*100,
     type="l", 
     xlab=paste("Days after", rnames[t_0]),
     ylab="Cases (%)",
     ylim=c(0,100),
     xlim=c(0,250), 
     lty=1, col=1)
lines(sol$t, sol$u[,6]/x_max*100, lty=2, col=2)
lines(sol$t, sol$u[,5]/x_max*100, lty=3, col=3)
lines(sol$t, sol$u[,4]/x_max*100, lty=4, col=4)
lines(sol$t, sol$u[,3]/x_max*100, lty=5, col=5)
lines(sol$t, sol$u[,2]/x_max*100, lty=6, col=6)
lines(sol$t, sol$u[,1]/x_max*100, lty=7, col=16)
legend("topright", legend <- c("Susceptibles (noninfected)", "Exposed (incubation time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Total Infected"), 
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

par(mar=c(5,6,6,5)+0.1)

plot(sol$t, sol$u[,4], 
     type="l",
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases", 
     ylim=c(0,max(sol$u[,2],sol$u[,4])), 
     xlim=c(0,250), 
     lty=1, col=1)
lines(sol$t, sol$u[,2], type="l", lty=2, col=2)
par(new=T)
plot(sol$t, sol$u[,4] / nh_max * 100, 
     type="l",
     xlab="", 
     ylab="",
     axes = F,
     xlim=c(0,250), 
     lty=3, col=3)
axis(4, pretty(sol$u[,4] / nh_max * 100,5))
mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hostpital Workload"), 
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

par(mar=c(5,6,6,5)+0.1)

plot(sol$t, sol$u[,4]/x_max*100, 
     type="l",
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases (%)", 
     ylim=c(0,max(sol$u[,2],sol$u[,4])/x_max*100), 
     xlim=c(0,250), 
     lty=1, col=1)
lines(sol$t, sol$u[,2]/x_max*100, type="l", lty=2, col=2)
par(new=T)
plot(sol$t, sol$u[,4] / nh_max * 100, 
     type="l",
     xlab="", 
     ylab="",
     axes = F,
     xlim=c(0,250), 
     lty=3, col=3)
axis(4, pretty(sol$u[,4] / nh_max * 100,5))
mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hostpital Workload"), 
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



