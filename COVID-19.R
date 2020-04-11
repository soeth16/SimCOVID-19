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
JuliaCall::julia_library("Distributed")
JuliaCall::julia_eval("addprocs(10)")


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

f = JuliaCall::julia_eval("function f(du, u, h, p, t)
    # parameter
    # k     = p[1]  # growth rate
    # tk_2  = p[2]  # time of change from k to k2
    # k2    = p[3]  # k after change
    du[8] = 0
    if (t >= p[4]) 
      u[8] = p[5]
    elseif (t >= p[2]) 
      u[8] = p[3]
    elseif (t < 1) 
      u[8] = p[1]
    end
    
    # kd
    if (u[4] > nh_max) 
      kd = kd_1 * nh_max / u[4] + kd_2 * (u[4] - nh_max) / u[4]
    else 
      kd = kd_1
    end
    
    # Susceptibles
    du[7] = - u[8] * u[7] / x_max * u[6]
    
    # Exposed / Incubating
    du[6] = + u[8] * u[7] / x_max * u[6] - h(p, t - te)[8] * h(p, t - te)[7] / x_max * h(p, t - te)[6] 
    
    # Infected
    du[5] = + h(p, t - te)[8] * h(p, t - te)[7] / x_max * h(p, t - te)[6] - h(p, t - te - ti)[8] * h(p, t - te - ti)[7] / x_max * h(p, t - te - ti)[6] * (1-kh)  - h(p, t - te - th)[8] * h(p, t - te - th)[7] / x_max * h(p, t - te - th)[6] * kh

    # Hostspitalisation 
    du[4] = + h(p, t - te - th)[8] * h(p, t - te - th)[7] / x_max * h(p, t - te - th)[6] * kh - h(p, t - te - th)[8] * h(p, t - te - th - thi)[7] / x_max * h(p, t - te - th - thi)[6] * kh
    
    # Recovered
    du[3] = + h(p, t - te - ti)[8] * h(p, t - te - ti)[7] / x_max * h(p, t - te - ti)[6] * (1-kh) + h(p, t - te - th - thi)[8] * h(p, t - te - th - thi)[7] / x_max * h(p, t - te - th - thi)[6] * (kh - kd)
    
    # Deaths
    du[2] = + h(p, t - te - th - thi)[8] * h(p, t - te - th - thi)[7] / x_max * h(p, t - te - th - thi)[6] * kd
  
    # Confirmed
    du[1] = du[5] + du[4] + du[3] + du[2]
    
  end")

JuliaCall::julia_assign("lags", c(te, te + ti, te + th,  te + th + thi))


JuliaCall::julia_eval("cb = ContinuousCallback((u,t,integrator)->(t-6),(integrator)->(integrator.u[8] = 0.05))")


# dx/dt = k * x
# 1/x * dx = k * dt
# integrate
# ln(x) - ln(x_0) = k * t
# x / x_0 = exp(k * t)
# x = x_0 exp(k * t)
# x_0 = 1

h = JuliaCall::julia_eval("function h(p, t)
    u = u0 * exp(p[1] * t)
    u[8] = p[1]
    u
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
JuliaCall::julia_assign("p", c(k, tk_2, k2, tk_2+7, k2))
JuliaCall::julia_assign("saveat", c(0:1000/1000*250))
#JuliaCall::julia_assign("lags", lags)
JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")

data <- array(dim=c(3,length(t),1000))
for (j in 1:length(t))
{
  data[1,j,] <- rnorm(1000,data_df[j,1],10)
  data[2,j,] <- rnorm(1000,data_df[j,2],100+data_df[j,2]/1e5)
  data[3,j,] <- rnorm(1000,data_df[j,3],100+data_df[j,3]/1e5)
}
JuliaCall::julia_assign("t", t-t_0)
JuliaCall::julia_assign("data", data)
JuliaCall::julia_eval("distributions = [fit_mle(Normal,data[i,j,:]) for i in 1:3, j in 1:length(t)]")
JuliaCall::julia_eval("obj = build_loss_objective(prob,MethodOfSteps(Rodas5()),LogLikeLoss(t,distributions),maxiters=10000,verbose=false); nothing")

# L2Loss
#data <- array(dim=c(1,length(t),1))
#for ( i in 1:1)
#  for (j in 1:length(t))
#    data[i,j,] <- data_df[j,i]
#JuliaCall::julia_assign("t", t-t_0)
#JuliaCall::julia_assign("data", data)
#JuliaCall::julia_eval("obj = build_loss_objective(prob,Tsit5(),L2Loss(t,data), maxiters=10000,verbose=false); nothing")


JuliaCall::julia_eval("bound1 = Tuple{Float64,Float64}[(0.25,0.4),(10,15),(0.2,0.3),(20,25),(0.15,0.25)]")
JuliaCall::julia_eval("res = bboptimize(obj;SearchRange = bound1, MaxSteps = 1e4)")

p2 <- JuliaCall::julia_eval("p = best_candidate(res)")

rnames[p2[2]+t_0]
rnames[p2[4]+t_0]

#####################################################################################
#
# Output
#   
#####################################################################################
# grapfical: 0
# pdf: 1
# png: 2
#plot_out <- 0
for (plot_out in c(2:0)) {
  
  if (plot_out == 1) pdf("Plots.pdf")
  
#####################################################################################
#
# Modell vs Situation
#   
#####################################################################################

JuliaCall::julia_assign("p", p2)
JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
JuliaCall::julia_eval("sol = solve(prob, MethodOfSteps(Rodas5()), reltol=1e-5, abstol=1e-6, maxiters = 10000, saveat=saveat); nothing")
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
# Situation
#   
#####################################################################################

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
for (i in c(0:25)) 
  lines(c(i*7-22,i*7-22), (c(-10, 25)), type="l", lty = 5, col="light gray" )

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
for (i in c(0:25)) 
  lines(c(i*7-22-tc_0,i*7-22-tc_0), (c(-10, 25)), type="l", lty = 5, col="light gray" )
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
for (i in c(0:100)) 
  lines(c(i*7-2-tc_0,i*7-2-tc_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
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
for (i in c(0:100)) 
  lines(c(i*7-2-tc_0,i*7-2-tc_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
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
for (i in c(0:100)) 
  lines(c(i*7-2-tc_0,i*7-2-tc_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
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
for (i in c(0:100)) 
  lines(c(i*7-2-tc_0,i*7-2-tc_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
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




#####################################################################################
#
# week bevore 03/23/2020
#   
#####################################################################################

JuliaCall::julia_assign("p", c(p2[1], p2[2], p2[1], p2[4], p2[1]))
JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
JuliaCall::julia_eval("sol = solve(prob, Tsit5(), reltol=1e-5, abstol=1e-6, maxiters = 10000, saveat=saveat); nothing")
sol = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))


if (plot_out == 2) png("Forecast-1-old.png", width = 640, height = 480)
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
for (i in c(0:100)) 
  lines(c(i*7-2-tc_0,i*7-2-tc_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
legend("topright", legend <- c("Susceptibles (noninfected)", "Exposed (incubation time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Total Infected"),       col=c(1:6,16),
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
for (i in c(0:100)) 
  lines(c(i*7-2-tc_0,i*7-2-tc_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
legend("topright", legend <- c("Susceptibles (noninfected)", "Exposed (incubation time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Total Infected"), 
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
for (i in c(0:100)) 
  lines(c(i*7-2-tc_0,i*7-2-tc_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
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

if (plot_out == 2) png("Forecast-ARDS-2-old.png", width = 640, height = 480)

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
for (i in c(0:100)) 
  lines(c(i*7-2-tc_0,i*7-2-tc_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
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


#####################################################################################
#
# variant calculation
#   
#####################################################################################

JuliaCall::julia_assign("p", c(p2[1], p2[2], p2[1], p2[3], p2[1]))
JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
JuliaCall::julia_eval("sol = solve(prob, Tsit5(), reltol=1e-5, abstol=1e-6, maxiters = 10000, saveat=saveat); nothing")
sol1 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))

rbcol = rainbow(11)

if (plot_out == 2) png("Forecast-ARDS-3.png", width = 640, height = 480)

par(mar=c(5,6,6,5)+0.1)

plot(sol1$t, (sol1$u[,4]), 
     type="l",
     axes=F, xlab="", ylab="", 
     ylim=c(-0.2,1)*max((sol1$u[,2])),
     xlim=c(0,150), 
     lty=1, col=1)
lines(sol1$t, (sol1$u[,2]), type="l", lty=2, col=1)

axis(2, pretty(range(sol1$u[,2]),10), col="black",las=1)  ## las=1 makes horizontal labels
mtext("Cases",side=2,line=4.5,las=0)
axis(1,pretty(range(sol1$t),5))
mtext(paste("Time (Days after ", rnames[t_0],")", sep = ""),side=1,col="black",line=2.5)


for (i in c(0:10))
{
  tk_2 <- 7
  R0_2 <- 0.4 + i * 0.2
  f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - log(R0_2 + 1) / te
  k2 <- uniroot(f_k2 , c(0.01, 1))$root
  
  JuliaCall::julia_assign("p", c(p2[1], tk_2, k2, tk_2+1, k2))
  JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
  JuliaCall::julia_eval("sol = solve(prob, Tsit5(), reltol=1e-5, abstol=1e-6, maxiters = 10000, saveat=saveat); nothing")
  sol2 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
  
  par(new=T)
  plot(sol2$t, (sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max((sol1$u[,2])), xlim=c(0,150))
  lines(sol2$t, (sol2$u[,2]), type="l", lty=2, col=rbcol[i])
  
  par(new=T)
  plot(sol2$t, sol2$u[,8], ylim = c(0.2,1.5), xlim=c(0,150), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
  if(i==0)
  {
    axis(4, pretty(range(sol2$u[,8]),10))
    mtext(expression('Growths Rate '*(D^-1)), side=4,line=3,las=0)
  }
}

par(new=T)
plot(c(0,250), log(c(nh_max, nh_max)), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(0,1)*max(log(sol1$u[,2])), xlim=c(0,150))
for (i in c(0:25)) 
  lines(c(i*7+2,i*7+2), (c(-10, 120)), type="l", lty = 5, col="light gray" )
lines(c(23,23), (c(-10, 25)), type="l", lty = 5, col="grey" )


legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", paste("R0 =", 0.4 + c(0:10)*.2)), 
       col=c(1,1,rbcol[0:11]),
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
     ylim=c(-0.2,1)*max((sol1$u[,2])),
     xlim=c(0,150), 
     lty=1, col=1)
lines(sol1$t, (sol1$u[,2]), type="l", lty=2, col=1)

axis(2, pretty(range(sol1$u[,2]),10), col="black",las=1)  ## las=1 makes horizontal labels
mtext("Cases",side=2,line=4.5,las=0)
axis(1,pretty(range(sol1$t),5))
mtext(paste("Time (Days after ", rnames[t_0],")", sep = ""),side=1,col="black",line=2.5)


for (i in c(0:10))
{
  tk_2 <- 17 + i * 3
  R0_2 <- 1
  f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - log(R0_2 + 1) / te
  k2 <- uniroot(f_k2 , c(0.01, 1))$root
  
  JuliaCall::julia_assign("p", c(p2[1], tk_2, k2, tk_2+1, k2))
  JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
  JuliaCall::julia_eval("sol = solve(prob, Tsit5(), reltol=1e-5, abstol=1e-6, maxiters = 10000, saveat=saveat); nothing")
  sol2 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
  
  par(new=T)
  plot(sol2$t, (sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max((sol1$u[,2])), xlim=c(0,150))
  lines(sol2$t, (sol2$u[,2]), type="l", lty=2, col=rbcol[i])
  
  par(new=T)
  plot(sol2$t, sol2$u[,8], ylim = c(0.2,1.5), xlim=c(0,150), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
  if(i==0)
  {
    axis(4, pretty(range(sol2$u[,8]),10))
    mtext(expression('Growths Rate '*(D^-1)), side=4,line=3,las=0)
  }
}

par(new=T)
plot(c(0,250), (c(nh_max, nh_max)), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(0,1)*max(log(sol1$u[,2])), xlim=c(0,150))
for (i in c(0:25)) 
  lines(c(i*7+2,i*7+2), (c(-10, 25)), type="l", lty = 5, col="light gray" )
lines(c(23,23), (c(-10, 25)), type="l", lty = 5, col="grey" )

legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", paste("t =", c(0:10)*3+17)), 
       col=c(1,1,rbcol[0:11]),
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
sol1$u[log(sol1$u[,2])<=0,6] <- 1 

plot(sol1$t, log(sol1$u[,4]), 
     type="l",
     axes=F, xlab="", ylab="", 
     ylim=c(-0.2,1)*max(log(sol1$u[,4])),
     xlim=c(0,150), 
     lty=1, col=1)
lines(sol1$t, log(sol1$u[,2]), type="l", lty=2, col=1)

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
  
  JuliaCall::julia_assign("p", c(p2[1], tk_2, k2, tk_2+1, k2))
  JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
  JuliaCall::julia_eval("sol = solve(prob, Tsit5(), reltol=1e-5, abstol=1e-6, maxiters = 10000, saveat=saveat); nothing")
  sol2 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
  
  par(new=T)
  sol2$u[log(sol2$u[,4])<=0,4] <- 1 
  sol2$u[log(sol2$u[,2])<=0,6] <- 1 
  plot(sol2$t, log(sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max(log(sol1$u[,4])),xlim=c(0,150))
  lines(sol2$t, log(sol2$u[,2]), type="l", lty=2, col=rbcol[i])
  
  par(new=T)
  plot(sol2$t, sol2$u[,8], ylim = c(0.2,1.5), xlim=c(0,150), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
  if(i==0)
  {
    axis(4, pretty(range(sol2$u[,8]),10))
    mtext(expression('Growths Rate '*(D^-1)), side=4,line=2.5,las=0)
  }
}

par(new=T)
plot(c(0,250), log(c(nh_max, nh_max)), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(0,1)*max(log(sol1$u[,2])), xlim=c(0,150))
for (i in c(0:25)) 
  lines(c(i*7+2,i*7+2), (c(-10, 25)), type="l", lty = 5, col="light gray" )
lines(c(23,23), (c(-10, 25)), type="l", lty = 5, col="grey" )

legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", paste("t =", c(0:10)*3+17)), 
       col=c(1,1,rbcol[0:11]),
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

