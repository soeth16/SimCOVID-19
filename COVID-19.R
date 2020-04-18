# set current working directory
setwd("~/SimCOVID-19")
# setwd("C:/Users/soeth/Desktop/SimCOVID-19")

# git pull
setwd("COVID-19")
system("git pull")
setwd("..")

library(JuliaCall)
#JuliaCall::julia_install_package_if_needed("DifferentialEquations")
#JuliaCall::julia_install_package_if_needed("DiffEqParamEstim")
#JuliaCall::julia_install_package_if_needed("Optim")
#JuliaCall::julia_install_package_if_needed("Distributions")
#JuliaCall::julia_install_package_if_needed("BlackBoxOptim")
JuliaCall::julia_library("Distributed")
JuliaCall::julia_eval("addprocs(8 - nprocs(), topology=:all_to_all, lazy = true)")
JuliaCall::julia_library("DifferentialEquations")
JuliaCall::julia_library("DiffEqParamEstim")
JuliaCall::julia_library("Optim")
JuliaCall::julia_library("Distributions")
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
ti = 12
th = 8
thi = 14
thd = 2


# assume k_raw from ln(x)-plot
k_raw = ln_res_2a$coefficients[2]
k_raw2 = ln_res_3a$coefficients[2]

# correct k_raw by the influece of te, ti, th, thi and phi_h
# x2'   = + k(t) x1(t) / n_max x2(t) - k(t - te) * x1 (t - te) / n_max * x2(t - te)
# x3:6' = + k(t - te) * x1 (t - te) / n_max * x2(t - te) # x3:6 as kumulativ infectet --> x3:6' = k_raw * x3:6 && x1 / n_max = 1 && k(t) == const.
#         + k * x1 / n_max * x2(t - te) = k_raw * x3:6(t) # xn(t0 + t) = xn(t0) * exp(k_raw * t) in exponetial phase
#       = + k * x2(t) * exp(-k_raw * te)
# x2(t) = x3:6(t) * k_raw / (k * exp(-k_raw * te))
# x2' = 

f_k <- function(k) k *(1-exp(-k * te)) - k_raw
k1 <- uniroot(f_k , c(0.01, 1))$root
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

# assume median from 40:80 in germany (good)
phi_d_1 <- as.numeric(quantile(
  ger_data_deaths[50:60]/(ger_data_confirmed[40:50]), 
  probs = 0.1
))

# assume Hubei / China (bad)
phi_d_2 <- as.numeric(quantile(
  subset(csse_covid_19_time_series_deaths, `Province/State` == "Hubei", select = c(20:80)) 
  /
    (
      subset(csse_covid_19_time_series_confirmed, `Province/State` == "Hubei", select = c(10:70))
    ),
  probs = 0.9
))


# assume phi_h from max phi_d_2 (Hubei / China)
phi_h <- 0.14

# http://www.gbe-bund.de/gbe10/I?I=838:37792217D 
n_h_max <- 28031 

# assume population of germany from https://de.wikipedia.org/wiki/Deutschland (03/13/20)
n_max <- 83.019213e6

#####################################################################################
#
# Modell
#   
#####################################################################################
#        1  2     3
# p <- c(k, tk_2, k2)
# 
# k     = p[1]  # 1st growth rate
# tk_2  = p[2]  # time of change to 2nd growth rate
# k2    = p[3]  # 2nd growth rate
# tk_3  = p[2]  # time of change to 3rd  growth rate
# k3    = p[3]  # 3rd growth rate



JuliaCall::julia_assign("n_max",n_max) 
JuliaCall::julia_assign("te", te) 
JuliaCall::julia_assign("ti", ti) 
JuliaCall::julia_assign("th", th)
JuliaCall::julia_assign("thi", thi)
JuliaCall::julia_assign("thd", thd)
JuliaCall::julia_assign("n_h_max", n_h_max)
JuliaCall::julia_assign("phi_h", phi_h)
JuliaCall::julia_assign("phi_d_1", phi_d_1)
JuliaCall::julia_assign("phi_d_2", phi_d_2)


JuliaCall::julia_eval("@everywhere n_max = $n_max")
JuliaCall::julia_eval("@everywhere te = $te")
JuliaCall::julia_eval("@everywhere ti = $ti") 
JuliaCall::julia_eval("@everywhere th = $th")
JuliaCall::julia_eval("@everywhere thi = $thi")
JuliaCall::julia_eval("@everywhere thd = $thd")
JuliaCall::julia_eval("@everywhere n_h_max = $n_h_max")
JuliaCall::julia_eval("@everywhere phi_h = $phi_h")
JuliaCall::julia_eval("@everywhere phi_d_1 = $phi_d_1")
JuliaCall::julia_eval("@everywhere phi_d_2 = $phi_d_2")

k = JuliaCall::julia_eval("@everywhere function k(p, t)
    if (t >= p[6]) 
      p[7]
    elseif (t >= p[4] && t < p[6]) 
      p[5]
    #if (t >= p[4]) 
    #  p[5]
    elseif (t >= p[2] && t < p[4]) 
      p[3]
    elseif (t < p[2]) 
      p[1]
    end
  end")

k <- function(p, t) {
  tmp <- 0
  for (tt in 1:length(t))
  {
    if (t[tt] >= p[6]) tmp[tt] <- p[7]
    else if (t[tt] >= p[4] && t[tt] < p[6]) tmp[tt] <- p[5]
    else if (t[tt] >= p[2] && t[tt] < p[4]) tmp[tt] <- p[3]
    else if (t[tt] < p[2]) tmp[tt] <- p[1]
  }
  return (tmp)
}


phi_d = JuliaCall::julia_eval("@everywhere function phi_d(u, p)
    # phi_d for distributed hospital usage
    #if (u[4] > n_h_max) 
    #  phi_d_1 * n_h_max / u[4] + phi_d_2 * (u[4] - n_h_max) / u[4]
    #else 
    #  phi_d_1
    #end
    
    # phi_d for hotspots without distrubution
    if (u[4] / n_h_max < 1) 
      phi_d_1 * (1 - u[4] / n_h_max) + phi_d_2 * u[4] / n_h_max
    else
      phi_d_2
    end
  end")

f = JuliaCall::julia_eval("@everywhere function f(du, u, h, p, t)

    # Susceptibles
    du[7] = (
      - k(p, t) * u[7] / n_max * u[6]
    )
    
    # Exposed / Incubating
    du[6] = (
      + k(p, t) * u[7] / n_max * u[6] 
      - k(p, t - te) * h(p, t - te)[7] / n_max * h(p, t - te)[6] 
    )
    
    # Infected
    du[5] = (
      + k(p, t - te) * h(p, t - te)[7] / n_max * h(p, t - te)[6] 
      - k(p, t - te - ti) * h(p, t - te - ti)[7] / n_max * h(p, t - te - ti)[6] * (1-phi_h)
      - k(p, t - te - th) * h(p, t - te - th)[7] / n_max * h(p, t - te - th)[6] * phi_h
    )

    # Hostspitalisation 
    du[4] = (
      + k(p, t - te - th) * h(p, t - te - th)[7] / n_max * h(p, t - te - th)[6] * phi_h 
      - k(p, t - te - th - thi) * h(p, t - te - th - thi)[7] / n_max * h(p, t - te - th - thi)[6] * (phi_h - phi_d(h(p, t - thi + thd), p))
      - k(p, t - te - th - thd) * h(p, t - te - th - thd)[7] / n_max * h(p, t - te - th - thd)[6] * phi_d(u, p)
    )
    
    # Recovered
    du[3] = ( 
      + k(p, t - te - ti) * h(p, t - te - ti)[7] / n_max * h(p, t - te - ti)[6] * (1-phi_h) 
      + k(p, t - te - th - thi) * h(p, t - te - th - thi)[7] / n_max * h(p, t - te - th - thi)[6] * (phi_h - phi_d(h(p, t - thi + thd), p))
    )
    
    # Deaths
    du[2] = ( 
      + k(p, t - te - th - thd) * h(p, t - te - th - thd)[7] / n_max * h(p, t - te - th - thd)[6] * phi_d(u, p)
    )
    
    # Confirmed
    du[1] = du[5] + du[4] + du[3] + du[2]
    
  end")

JuliaCall::julia_assign("lags", c(te, te + ti, te + th,  te + th + thi, te + th + thd, thi - thd))


JuliaCall::julia_eval("cb = ContinuousCallback((u,t,integrator)->(t-6),(integrator)->(integrator.u[8] = 0.05))")


# dx/dt = k * x
# 1/x * dx = k * dt
# integrate
# ln(x) - ln(x_0) = k * t
# x / x_0 = exp(k * t)
# x = x_0 exp(k * t)
# x_0 = 1

h = JuliaCall::julia_eval("@everywhere function h(p, t)
    u = u0 * exp(p[1] * t)
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
E <- I * k_raw/(k1 * exp(-k1*te))
I <- I * 0.95
H <- I * 0.05
R <- ger_data_recovered
D <- ger_data_deaths
S <-n_max - I - E - H - R - D

data_df <- data.frame(C[t], D[t], R[t])

JuliaCall::julia_assign("u0", c(C[t_0], D[t_0], R[t_0], H[t_0], I[t_0], E[t_0], S[t_0]))
JuliaCall::julia_assign("tspan", c(0,250))
JuliaCall::julia_assign("p", c(k1, tk_2, k2, tk_2+7, k2,1000,k1))
JuliaCall::julia_assign("saveat", c(0:1000/1000*250))
JuliaCall::julia_eval("@everywhere u0 = $u0");
JuliaCall::julia_eval("@everywhere tspan = $tspan");
JuliaCall::julia_eval("@everywhere p = $p");
JuliaCall::julia_eval("@everywhere saveat = $saveat");
JuliaCall::julia_eval("@everywhere lags = $lags");
JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")

data <- array(dim=c(3,length(t),1000))
for (j in 1:length(t))
{
  data[1,j,] <- rnorm(1000,data_df[j,1],max(data_df[,1])*0.01)
  data[2,j,] <- rnorm(1000,data_df[j,2],(max(data_df[,2])+data_df[j,2])*0.01)
  data[3,j,] <- rnorm(1000,data_df[j,3],(max(data_df[,3])+data_df[j,3])*0.01)
}
JuliaCall::julia_assign("t", t-t_0)
JuliaCall::julia_assign("data", data)
JuliaCall::julia_eval("@everywhere t = $t")
JuliaCall::julia_eval("@everywhere data = $data")
JuliaCall::julia_eval("distributions = [fit_mle(Normal,data[i,j,:]) for i in 1:3, j in 1:length(t)]")
JuliaCall::julia_eval("obj = build_loss_objective(prob, Rodas5(), reltol=1e-6, abstol=1e-9, maxiters = 10000, LogLikeLoss(t,distributions), verbose=true); nothing")



JuliaCall::julia_eval("bound1 = Tuple{Float64,Float64}[(0.25,0.4),(10,15),(0.2,0.3),(20,23),(0.15,0.25),(1000,1000+1e-9),(0.15,0.15+1e-9)]")
#JuliaCall::julia_eval("res1 = bboptimize(obj;SearchRange = bound1, MaxSteps = 1e3, NumDimensions = 5,
#    Workers = workers(),
#    Method = :dxnes)")

#p2 <- JuliaCall::julia_eval("p = best_candidate(res1)")
p2 <- c(0.317511, 11.7132, 0.251697, 21.9461, 0.179116, 1000, 0.1)
p2
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
  JuliaCall::julia_eval("sol = solve(prob, Rodas5(), reltol=1e-7, abstol=1e-9, maxiters = 10000, saveat=saveat); nothing")
  sol = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
  
  if (plot_out == 2) png("Model_vs_Situation-1.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  plot(c(t_0:len_ger_data)-t_0, ger_data_confirmed[c(t_0:len_ger_data),1], 
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Confirmed cases",
       col=1,pch=1)
  points(c(t_0:len_ger_data)-t_0, ger_data_recovered[c(t_0:len_ger_data),1],col=3,pch=2)
  points(c(t_0:len_ger_data)-t_0, ger_data_deaths[c(t_0:len_ger_data),1],col=2,pch=3)
  lines(sol$t, sol$u[,1], lty=2, col=1)
  lines(sol$t, sol$u[,2], lty=2, col=2)
  lines(sol$t, sol$u[,3], lty=2, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Confirmed cases", "Recovered cases", "Deaths", "Modelled cases"),
         col=c(1,3,2,1),
         bty="n",
         pch=c(1,2,3,-1),
         lwd=1,
         lty=c(-1,-1,-1,2),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Situation vs Model COVID-19 in Germany",
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  if (plot_out == 2) png("Model_vs_Situation-2.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  plot(c(t_0:len_ger_data)-t_0, ger_data_confirmed[c(t_0:len_ger_data),1]/n_max*100, 
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases (%)",
       pch=1)
  points(c(t_0:len_ger_data)-t_0, ger_data_recovered[c(t_0:len_ger_data),1]/n_max*100,col=3,pch=2)
  points(c(t_0:len_ger_data)-t_0, ger_data_deaths[c(t_0:len_ger_data),1]/n_max*100,col=2,pch=3)
  lines(sol$t, sol$u[,1]/n_max*100, lty=2, col=1)
  lines(sol$t, sol$u[,2]/n_max*100, lty=2, col=2)
  lines(sol$t, sol$u[,3]/n_max*100, lty=2, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Confirmed cases", "Recovered cases", "Deaths", "Modelled cases"),
         col=c(1,3,2,1),
         bty="n",
         pch=c(1,2,3,-1),
         lwd=1,
         lty=c(-1,-1,-1,2),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Situation vs Model COVID-19 in Germany",
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  #####################################################################################
  #
  # Situation
  #   
  #####################################################################################
  
  if (plot_out == 2) png("Situation-1.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  plot(ger_data_confirmed, 
       type="l", 
       col = 1, 
       lty = 1,
       xlab=paste("Days after", rnames[1]), 
       ylab="Cases")
  lines(ger_data_recovered , type="l", col = 2, lty = 2)
  for (i in c(0:100)) 
    lines(c(i*7-22,i*7-22), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
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
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  if (plot_out == 2) png("Situation-2.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
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
  for (i in c(0:100)) 
    lines(c(i*7-22,i*7-22), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
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
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  if (plot_out == 2) png("Situation-3.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
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
  for (i in c(0:100)) 
    lines(c(i*7-22,i*7-22), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
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
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  #####################################################################################
  #
  # Forecast
  #   
  #####################################################################################
  
  if (plot_out == 2) png("Forecast-1.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  plot(sol$t,sol$u[,7],
       type="l", 
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases", 
       ylim=c(0,1e5), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,6], lty=2, col=2)
  lines(sol$t, sol$u[,5], lty=3, col=3)
  lines(sol$t, sol$u[,4], lty=4, col=4)
  lines(sol$t, sol$u[,3], lty=5, col=5)
  lines(sol$t, sol$u[,2], lty=6, col=6)
  lines(sol$t, sol$u[,1], lty=7, col=16)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topright", legend <- c("Susceptibles (Noninfected)", "Exposed (Incubation Time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Confirmed Cases (Total Infected)"),       col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  if (plot_out == 2) png("Forecast-2.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  plot(sol$t,sol$u[,7]/n_max*100,
       type="l", 
       xlab=paste("Days after", rnames[t_0]),
       ylab="Cases (%)",
       ylim=c(0,100),
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,6]/n_max*100, lty=2, col=2)
  lines(sol$t, sol$u[,5]/n_max*100, lty=3, col=3)
  lines(sol$t, sol$u[,4]/n_max*100, lty=4, col=4)
  lines(sol$t, sol$u[,3]/n_max*100, lty=5, col=5)
  lines(sol$t, sol$u[,2]/n_max*100, lty=6, col=6)
  lines(sol$t, sol$u[,1]/n_max*100, lty=7, col=16)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topright", legend <- c("Susceptibles (Noninfected)", "Exposed (Incubation Time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Confirmed Cases (Total Infected)"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  max(sol$u[,4])
  max(sol$u[,4])/n_max*100
  max(sol$u[,4])/n_h_max
  
  
  
  if (plot_out == 2) png("Forecast-ARDS-1.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  
  plot(sol$t, sol$u[,4], 
       type="l",
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases", 
       ylim=c(0,max(sol$u[,2],sol$u[,4])), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,2], type="l", lty=2, col=2)
  par(new=T)
  plot(sol$t, sol$u[,4] / n_h_max * 100, 
       type="l",
       xlab="", 
       ylab="",
       axes = F,
       xlim=c(0,250), 
       lty=3, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  axis(4, pretty(sol$u[,4] / n_h_max * 100,5))
  mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hostpital Workload"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  if (plot_out == 2) png("Forecast-ARDS-2.png", width = 640, height = 480)
  
  par(mar=c(5,6,7,5)+0.1)
  
  plot(sol$t, sol$u[,4]/n_max*100, 
       type="l",
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases (%)", 
       ylim=c(0,max(sol$u[,2],sol$u[,4])/n_max*100), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,2]/n_max*100, type="l", lty=2, col=2)
  par(new=T)
  plot(sol$t, sol$u[,4] / n_h_max * 100, 
       type="l",
       xlab="", 
       ylab="",
       axes = F,
       xlim=c(0,250), 
       lty=3, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  axis(4, pretty(sol$u[,4] / n_h_max * 100,5))
  mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hostpital Workload"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  
  
  #####################################################################################
  #
  # Scenario 1 no social distancing and no shutdown
  #   
  #####################################################################################
  
  JuliaCall::julia_assign("p", c(p2[1], 1000, p2[1], 1001, p2[1], 1002, p2[1]))
  JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
  JuliaCall::julia_eval("sol = solve(prob, Rodas5(), reltol=1e-7, abstol=1e-9, maxiters = 100000, saveat=saveat); nothing")
  sol = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
  
  
  if (plot_out == 2) png("Forecast-1-scenario-1.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  plot(sol$t,sol$u[,7],
       type="l", 
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases", 
       ylim=c(0,n_max), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,6], lty=2, col=2)
  lines(sol$t, sol$u[,5], lty=3, col=3)
  lines(sol$t, sol$u[,4], lty=4, col=4)
  lines(sol$t, sol$u[,3], lty=5, col=5)
  lines(sol$t, sol$u[,2], lty=6, col=6)
  lines(sol$t, sol$u[,1], lty=7, col=16)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topright", legend <- c("Susceptibles (Noninfected)", "Exposed (Incubation Time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Confirmed Cases (Total Infected)"),
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 1 no social distancing and no shutdown", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  
  if (plot_out == 2) png("Forecast-2-scenario-1.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  plot(sol$t,sol$u[,7]/n_max*100,
       type="l", 
       xlab=paste("Days after", rnames[t_0]),
       ylab="Cases (%)",
       ylim=c(0,100),
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,6]/n_max*100, lty=2, col=2)
  lines(sol$t, sol$u[,5]/n_max*100, lty=3, col=3)
  lines(sol$t, sol$u[,4]/n_max*100, lty=4, col=4)
  lines(sol$t, sol$u[,3]/n_max*100, lty=5, col=5)
  lines(sol$t, sol$u[,2]/n_max*100, lty=6, col=6)
  lines(sol$t, sol$u[,1]/n_max*100, lty=7, col=16)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topright", legend <- c("Susceptibles (Noninfected)", "Exposed (Incubation Time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Confirmed Cases (Total Infected)"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Senario 1 no social distenining / no shutdown ", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  max(sol$u[,4])
  max(sol$u[,4])/n_max*100
  max(sol$u[,4])/n_h_max
  
  
  if (plot_out == 2) png("Forecast-ARDS-1-scenario-1.png", width = 640, height = 480)
  
  par(mar=c(5,6,7,5)+0.1)
  
  plot(sol$t, sol$u[,4], 
       type="l",
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases", 
       ylim=c(0,max(sol$u[,2],sol$u[,4])), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,2], type="l", lty=2, col=2)
  par(new=T)
  plot(sol$t, sol$u[,4] / n_h_max * 100, 
       type="l",
       xlab="", 
       ylab="",
       axes = F,
       xlim=c(0,250), 
       lty=3, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  axis(4, pretty(sol$u[,4] / n_h_max * 100,5))
  mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hostpital Workload"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 1 no social distancing and no shutdown", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  if (plot_out == 2) png("Forecast-ARDS-2-scenario-1.png", width = 640, height = 480)
  
  par(mar=c(5,6,7,5)+0.1)
  
  plot(sol$t, sol$u[,4]/n_max*100, 
       type="l",
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases (%)", 
       ylim=c(0,max(sol$u[,2],sol$u[,4])/n_max*100), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,2]/n_max*100, type="l", lty=2, col=2)
  par(new=T)
  plot(sol$t, sol$u[,4] / n_h_max * 100, 
       type="l",
       xlab="", 
       ylab="",
       axes = F,
       xlim=c(0,250), 
       lty=3, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  axis(4, pretty(sol$u[,4] / n_h_max * 100,5))
  mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hostpital Workload"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 1 no social distancing and no shutdown", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  #####################################################################################
  #
  # Scenario 2 exit after day 50 without social distancing (worest case)
  #   
  #####################################################################################
  
  JuliaCall::julia_assign("p", c(p2[1], p2[2], p2[3], p2[4], p2[5], 50, p2[1]))
  JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
  JuliaCall::julia_eval("sol = solve(prob, Rodas5(), reltol=1e-7, abstol=1e-9, maxiters = 100000, saveat=saveat); nothing")
  sol = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
  
  
  if (plot_out == 2) png("Forecast-1-scenario-2.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  plot(sol$t,sol$u[,7],
       type="l", 
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases", 
       ylim=c(0,n_max), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,6], lty=2, col=2)
  lines(sol$t, sol$u[,5], lty=3, col=3)
  lines(sol$t, sol$u[,4], lty=4, col=4)
  lines(sol$t, sol$u[,3], lty=5, col=5)
  lines(sol$t, sol$u[,2], lty=6, col=6)
  lines(sol$t, sol$u[,1], lty=7, col=16)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topright", legend <- c("Susceptibles (Noninfected)", "Exposed (Incubation Time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Confirmed Cases (Total Infected)"),       col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 2 exit after day 50 without social distancing ", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  
  if (plot_out == 2) png("Forecast-2-scenario-2.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  plot(sol$t,sol$u[,7]/n_max*100,
       type="l", 
       xlab=paste("Days after", rnames[t_0]),
       ylab="Cases (%)",
       ylim=c(0,100),
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,6]/n_max*100, lty=2, col=2)
  lines(sol$t, sol$u[,5]/n_max*100, lty=3, col=3)
  lines(sol$t, sol$u[,4]/n_max*100, lty=4, col=4)
  lines(sol$t, sol$u[,3]/n_max*100, lty=5, col=5)
  lines(sol$t, sol$u[,2]/n_max*100, lty=6, col=6)
  lines(sol$t, sol$u[,1]/n_max*100, lty=7, col=16)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topright", legend <- c("Susceptibles (Noninfected)", "Exposed (Incubation Time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Confirmed Cases (Total Infected)"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 2 exit after day 50 without social distancing ", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  max(sol$u[,4])
  max(sol$u[,4])/n_max*100
  max(sol$u[,4])/n_h_max
  
  
  if (plot_out == 2) png("Forecast-ARDS-1-scenario-2.png", width = 640, height = 480)
  
  par(mar=c(5,6,7,5)+0.1)
  
  plot(sol$t, sol$u[,4], 
       type="l",
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases", 
       ylim=c(0,max(sol$u[,2],sol$u[,4])), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,2], type="l", lty=2, col=2)
  par(new=T)
  plot(sol$t, sol$u[,4] / n_h_max * 100, 
       type="l",
       xlab="", 
       ylab="",
       axes = F,
       xlim=c(0,250), 
       lty=3, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  axis(4, pretty(sol$u[,4] / n_h_max * 100,5))
  mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hostpital Workload"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 2 exit after day 50 without social distancing ", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  if (plot_out == 2) png("Forecast-ARDS-2-scenario-2.png", width = 640, height = 480)
  
  par(mar=c(5,6,7,5)+0.1)
  
  plot(sol$t, sol$u[,4]/n_max*100, 
       type="l",
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases (%)", 
       ylim=c(0,max(sol$u[,2],sol$u[,4])/n_max*100), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,2]/n_max*100, type="l", lty=2, col=2)
  par(new=T)
  plot(sol$t, sol$u[,4] / n_h_max * 100, 
       type="l",
       xlab="", 
       ylab="",
       axes = F,
       xlim=c(0,250), 
       lty=3, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  axis(4, pretty(sol$u[,4] / n_h_max * 100,5))
  mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hostpital Workload"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 2 exit after day 50 without social distancing ", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  #####################################################################################
  #
  # Scenario 3 exit after day 50 with social distancing (better case)
  #   
  #####################################################################################
  
  JuliaCall::julia_assign("p", c(p2[1], p2[2], p2[3], p2[4], p2[5], 50, p2[3]))
  JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
  JuliaCall::julia_eval("sol = solve(prob, Rodas5(), reltol=1e-7, abstol=1e-9, maxiters = 100000, saveat=saveat); nothing")
  sol = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
  
  
  if (plot_out == 2) png("Forecast-1-scenario-3.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  plot(sol$t,sol$u[,7],
       type="l", 
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases", 
       ylim=c(0,n_max), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,6], lty=2, col=2)
  lines(sol$t, sol$u[,5], lty=3, col=3)
  lines(sol$t, sol$u[,4], lty=4, col=4)
  lines(sol$t, sol$u[,3], lty=5, col=5)
  lines(sol$t, sol$u[,2], lty=6, col=6)
  lines(sol$t, sol$u[,1], lty=7, col=16)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topright", legend <- c("Susceptibles (Noninfected)", "Exposed (Incubation Time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Confirmed Cases (Total Infected)"),       col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 3 exit after day 50 with social distancing", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  
  if (plot_out == 2) png("Forecast-2-scenario-3.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  plot(sol$t,sol$u[,7]/n_max*100,
       type="l", 
       xlab=paste("Days after", rnames[t_0]),
       ylab="Cases (%)",
       ylim=c(0,100),
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,6]/n_max*100, lty=2, col=2)
  lines(sol$t, sol$u[,5]/n_max*100, lty=3, col=3)
  lines(sol$t, sol$u[,4]/n_max*100, lty=4, col=4)
  lines(sol$t, sol$u[,3]/n_max*100, lty=5, col=5)
  lines(sol$t, sol$u[,2]/n_max*100, lty=6, col=6)
  lines(sol$t, sol$u[,1]/n_max*100, lty=7, col=16)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topright", legend <- c("Susceptibles (Noninfected)", "Exposed (Incubation Time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Confirmed Cases (Total Infected)"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 3 exit after day 50 with social distancing", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  max(sol$u[,4])
  max(sol$u[,4])/n_max*100
  max(sol$u[,4])/n_h_max
  
  
  if (plot_out == 2) png("Forecast-ARDS-1-scenario-3.png", width = 640, height = 480)
  
  par(mar=c(5,6,7,5)+0.1)
  
  plot(sol$t, sol$u[,4], 
       type="l",
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases", 
       ylim=c(0,max(sol$u[,2],sol$u[,4])), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,2], type="l", lty=2, col=2)
  par(new=T)
  plot(sol$t, sol$u[,4] / n_h_max * 100, 
       type="l",
       xlab="", 
       ylab="",
       axes = F,
       xlim=c(0,250), 
       lty=3, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  axis(4, pretty(sol$u[,4] / n_h_max * 100,5))
  mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hostpital Workload"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 3 exit after day 50 with social distancing", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  if (plot_out == 2) png("Forecast-ARDS-2-scenario-3.png", width = 640, height = 480)
  
  par(mar=c(5,6,7,5)+0.1)
  
  plot(sol$t, sol$u[,4]/n_max*100, 
       type="l",
       xlab=paste("Days after", rnames[t_0]), 
       ylab="Cases (%)", 
       ylim=c(0,max(sol$u[,2],sol$u[,4])/n_max*100), 
       xlim=c(0,250), 
       lty=1, col=1)
  lines(sol$t, sol$u[,2]/n_max*100, type="l", lty=2, col=2)
  par(new=T)
  plot(sol$t, sol$u[,4] / n_h_max * 100, 
       type="l",
       xlab="", 
       ylab="",
       axes = F,
       xlim=c(0,250), 
       lty=3, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  axis(4, pretty(sol$u[,4] / n_h_max * 100,5))
  mtext("Hostpital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hostpital Workload"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 3 exit after day 50 with social distancing", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  #####################################################################################
  #
  # variant calculation
  #   
  #####################################################################################
  
  JuliaCall::julia_assign("p", c(p2[1], p2[2], p2[1], p2[3], p2[1], p2[3], p2[1]))
  JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
  JuliaCall::julia_eval("sol = solve(prob, Rodas5(), reltol=1e-9, abstol=1e-12, maxiters = 10000, saveat=saveat); nothing")
  sol1 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
  
  rbcol = rainbow(11)
  
  if (plot_out == 2) png("Forecast-ARDS-3.png", width = 640, height = 480)
  
  par(mar=c(5,6,7,5)+0.1)
  
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
    tk_2 <- 22
    R0_2 <- 0.4 + i * 0.2
    f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - log(R0_2 + 1) / te
    k2 <- uniroot(f_k2 , c(0.01, 1))$root
    
    p_var <- c(p2[1], tk_2, k2, tk_2+1, k2, tk_2+2, k2)
    JuliaCall::julia_assign("p", p_var)
    JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
    JuliaCall::julia_eval("sol = solve(prob, Rodas5(), reltol=1e-7, abstol=1e-9, maxiters = 10000, saveat=saveat); nothing")
    sol2 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
    
    par(new=T)
    plot(sol2$t, (sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max((sol1$u[,2])), xlim=c(0,150))
    lines(sol2$t, (sol2$u[,2]), type="l", lty=2, col=rbcol[i])
    
    par(new=T)
    plot(sol2$t, k(p_var, sol2$t), ylim = c(0.2,1.5), xlim=c(0,150), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
    if(i==0)
    {
      axis(4, pretty(range(k(p_var, sol2$t)),10))
      mtext(expression('Growths Rate '*(D^-1)), side=4,line=3,las=0)
    }
  }
  
  par(new=T)
  plot(c(0,250), log(c(n_h_max, n_h_max)), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(0,1)*max(log(sol1$u[,2])), xlim=c(0,150))
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  axis(4, pretty(sol1$u[,4] / n_h_max * 100,5))
  lines(c(22,22), (c(-10, 25)), type="l", lty = 5, col="grey" )
  
  mtext("Week number (2020)",side=3,line=2,las=0)
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
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  
  rbcol = rainbow(11)
  
  if (plot_out == 2) png("Forecast-ARDS-4.png", width = 640, height = 480)
  
  par(mar=c(5,6,7,5)+0.1)
  
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
    tk_2 <- 16 + i * 3
    R0_2 <- 1
    f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - log(R0_2 + 1) / te
    k2 <- uniroot(f_k2 , c(0.01, 1))$root
    
    p_var <- c(p2[1], tk_2, k2, tk_2+1, k2, tk_2+2, k2)
    JuliaCall::julia_assign("p", p_var)
    JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
    JuliaCall::julia_eval("sol = solve(prob, Rodas5(), reltol=1e-9, abstol=1e-12, maxiters = 10000, saveat=saveat); nothing")
    sol2 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
    
    par(new=T)
    plot(sol2$t, (sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max((sol1$u[,2])), xlim=c(0,150))
    lines(sol2$t, (sol2$u[,2]), type="l", lty=2, col=rbcol[i])
    
    par(new=T)
    plot(sol2$t, k(p_var, sol2$t), ylim = c(0.2,1.5), xlim=c(0,150), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
    if(i==0)
    {
      axis(4, pretty(range(k(p_var, sol2$t)),10))
      mtext(expression('Growths Rate '*(D^-1)), side=4,line=3,las=0)
    }
  }
  
  par(new=T)
  plot(c(0,250), (c(n_h_max, n_h_max)), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(0,1)*max(log(sol1$u[,2])), xlim=c(0,150))
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  lines(c(22,22), (c(-10, 25)), type="l", lty = 5, col="grey" )
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", paste("t =", c(0:10)*3+19)), 
         col=c(1,1,rbcol[0:11]),
         bty="n",
         lwd=1,
         lty=c(1,2 ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  box()
  
  
  if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  
  
  rbcol = rainbow(11)
  
  if (plot_out == 2) png("Forecast-ARDS-5.png", width = 640, height = 480)
  
  par(mar=c(5,6,7,5)+0.1)
  
  sol1$u[log(sol1$u[,4])<=0,4] <- 1 
  sol1$u[log(sol1$u[,2])<=0,6] <- 1 
  
  plot(sol1$t, log(sol1$u[,4]), 
       type="l",
       axes=F, xlab="", ylab="", 
       ylim=c(-0.2,1)*max(log(sol1$u[,6])),
       xlim=c(0,150), 
       lty=1, col=1)
  lines(sol1$t, log(sol1$u[,2]), type="l", lty=2, col=1)
  
  axis(2, pretty(range(log(sol1$u[,4])),10), col="black",las=1)  ## las=1 makes horizontal labels
  mtext(expression(log['e']('Cases')),side=2,line=4,las=0)
  axis(1,pretty(range(sol1$t),5))
  mtext(paste("Time (Days after ", rnames[t_0],")", sep = ""),side=1,col="black",line=2.5)
  
  
  for (i in c(0:10))
  {
    tk_2 <- 16 + i * 3
    R0_2 <- 1
    f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - log(R0_2 + 1) / te
    k2 <- uniroot(f_k2 , c(0.01, 1))$root
    
    
    p_var <- c(p2[1], tk_2, k2, tk_2+1, k2, tk_2+2, k2)
    JuliaCall::julia_assign("p", p_var)
    JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
    JuliaCall::julia_eval("sol = solve(prob, Rodas5(), reltol=1e-9, abstol=1e-12, maxiters = 10000, saveat=saveat); nothing")
    sol2 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
    
    par(new=T)
    sol2$u[log(sol2$u[,4])<exp(1),4] <- 1 
    sol2$u[log(sol2$u[,2])<exp(1),6] <- 1 
    plot(sol2$t, log(sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max(log(sol1$u[,6])),xlim=c(0,150))
    lines(sol2$t, log(sol2$u[,2]), type="l", lty=2, col=rbcol[i])
    
    par(new=T)
    plot(sol2$t, k(p_var, sol2$t), ylim = c(0.2,1.5), xlim=c(0,150), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
    if(i==0)
    {
      axis(4, pretty(range(k(p_var, sol2$t)),10))
      mtext(expression('Growths Rate '*(D^-1)), side=4,line=2.5,las=0)
    }
  }
  
  par(new=T)
  plot(c(0,250), log(c(n_h_max, n_h_max)), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(0,1)*max(log(sol1$u[,2])), xlim=c(0,150))
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  lines(c(22,22), (c(-10, 25)), type="l", lty = 5, col="grey" )
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", paste("t =", c(0:10)*3+19)), 
         col=c(1,1,rbcol[0:11]),
         bty="n",
         lwd=1,
         lty=c(1,2 ,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  box()
  
  
  if (plot_out != 2) title("Forecast COVID-19 in Germany", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  
  
  min(sol1$t[(sol1$u[,4] > n_h_max)]) + 16
  
  
  # the end 
  if (plot_out == 1) dev.off()
}

