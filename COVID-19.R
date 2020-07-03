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
JuliaCall::julia_eval("addprocs(8 - nprocs(), topology=:all_to_all)")
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


#####################################################################################
#
# Constants
#   
#####################################################################################

# confirmed to death (determined by overlapping the two curves)
td <-10.5
tcd <- trunc(td)
# assume incubation time from https://www.ncbi.nlm.nih.gov/pubmed/32150748
te = 5.1
ti = 14 # determined by overlapping the two curves
th = 8 
thi = 10 
thd = td - th 
thii = 10 # == 4 weeks


# ln plots
t_0 <- 40
tc0_0 <- 15
tc0_max <- 35
tc1_0 <- 36
tc1_max <- t_0+11
tc2_0 <- t_0+12
tc2_max <- t_0+21
tc3_0 <- t_0+22
tc3_max <- t_0+49
tc4_0 <- t_0+50
tc4_max <- len_ger_data

dt_r = 14

# approx hidden cases
phi_d <- ger_data_deaths[t_0:len_ger_data]/ger_data_confirmed[(t_0-tcd):(len_ger_data-tcd)]
#phi_d_median <- quantile(phi_d[23:45],probs=0.5)
phi_d_median <- mean(phi_d[23:45])
phi_c <- phi_d / phi_d_median
phi_c[1:22] <- 1
plot(phi_c,col=1)
phi_r <-c(1:ti, phi_c[1:(length(phi_c)-ti)])
phi_r[1:ti] <- 1
points(phi_r,col=2)
phi_c <- lowess(1:length(phi_c),phi_c,f=0.3)$y
phi_r <- lowess(1:length(phi_r),phi_r,f=0.3)$y
lines(phi_c,col=1,lty=2)
lines(phi_r,col=2,lty=2)

ln_data_confirmed <- data.frame(t = 1:len_ger_data,ln_x=log(ger_data_confirmed[1:len_ger_data]))
ln_data_confirmed$`ln_x`[t_0:len_ger_data] <- log(ger_data_confirmed[t_0:len_ger_data]*phi_c)
ln_data_recovered <- data.frame(t = 1:len_ger_data,ln_x=log(ger_data_recovered[1:len_ger_data]))
ln_data_recovered$`ln_x`[t_0:len_ger_data] <- log(ger_data_recovered[t_0:len_ger_data]*phi_r)

ln_res_1a = lm(ln_x~t, ln_data_confirmed[tc0_0:tc0_max,])
ln_res_2a = lm(ln_x~t, ln_data_confirmed[tc1_0:tc1_max,])
ln_res_3a = lm(ln_x~t, ln_data_confirmed[tc2_0:tc2_max,])
ln_res_4a = lm(ln_x~t, ln_data_confirmed[tc3_0:tc3_max,])
ln_res_5a = lm(ln_x~t, ln_data_confirmed[tc4_0:tc4_max,])

ln_res_1b = lm(ln_x~t, ln_data_recovered[(tc0_0:tc0_max)+dt_r,])
ln_res_2b = lm(ln_x~t, ln_data_recovered[(tc1_0:tc1_max)+dt_r,])
ln_res_3b = lm(ln_x~t, ln_data_recovered[(tc2_0:tc2_max)+dt_r,])
ln_res_4b = lm(ln_x~t, ln_data_recovered[(tc3_0+dt_r):min(tc3_max+dt_r,len_ger_data),])

ln_res_1a
ln_res_2a
ln_res_3a
ln_res_4a
ln_res_5a
ln_res_1b
ln_res_2b
ln_res_3b
ln_res_4b


# assume k_raw from ln(x)-plot
k_raw = ln_res_2a$coefficients[2]
k_raw2 = ln_res_3a$coefficients[2]

k_raw
f_k <- function(k) k *(1-exp(-k * te)) - k_raw
k1 <- uniroot(f_k , c(0.01, 1))$root
f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - k_raw2
k2 <- uniroot(f_k2 , c(0.01, 1))$root

# assume median from the first three german lock down weeks (good)
phi_d_1 <- phi_d_median

# assume Hubei / China (bad)
phi_d_2 <- as.numeric(quantile(
  subset(csse_covid_19_time_series_deaths, `Province/State` == "Hubei", select = c(21:81)) /subset(csse_covid_19_time_series_confirmed, `Province/State` == "Hubei", select = c(10:70)),
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
# Model
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
JuliaCall::julia_assign("thii", thii)
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
JuliaCall::julia_eval("@everywhere thii = $thii")
JuliaCall::julia_eval("@everywhere thd = $thd")
JuliaCall::julia_eval("@everywhere n_h_max = $n_h_max")
JuliaCall::julia_eval("@everywhere phi_h = $phi_h")
JuliaCall::julia_eval("@everywhere phi_d_1 = $phi_d_1")
JuliaCall::julia_eval("@everywhere phi_d_2 = $phi_d_2")

JuliaCall::julia_eval("@everywhere function k(p, t)
    if (t >= p[12]) 
      p[13]
    elseif (t >= p[10] && t < p[12]) 
      p[11]
    elseif (t >= p[8] && t < p[10]) 
      p[9]
    elseif (t >= p[6] && t < p[8]) 
      p[7]
    elseif (t >= p[4] && t < p[6]) 
      p[5]
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
    if (t[tt] >= p[12]) tmp[tt] <- p[13]
    else if (t[tt] >= p[10] && t[tt] < p[12]) tmp[tt] <- p[11]
    else if (t[tt] >= p[8] && t[tt] < p[10]) tmp[tt] <- p[9]
    else if (t[tt] >= p[6] && t[tt] < p[8])  tmp[tt] <- p[7]
    else if (t[tt] >= p[4] && t[tt] < p[6])  tmp[tt] <- p[5]
    else if (t[tt] >= p[2] && t[tt] < p[4])  tmp[tt] <- p[3]
    else if (t[tt] < p[2]) tmp[tt] <- p[1]
  }
  return (tmp)
}


JuliaCall::julia_eval("@everywhere function phi_d(u, p)
    # phi_d for distributed hospital usage
    if (u > n_h_max) 
      phi_d_1 * n_h_max / u + phi_d_2 * (u - n_h_max) / u
    else 
      phi_d_1
    end
    
    # phi_d for hot spots without distribution
    #if (u / n_h_max < 1) 
    #  phi_d_1 * (1 - (u / n_h_max)) + phi_d_2 * u / n_h_max
    #else
    #  phi_d_2
    #end
  end")

f = JuliaCall::julia_eval("@everywhere function f(du, u, h, p, t)

    # Susceptibles
    du[7] = (
      - k(p, t) * u[7] / n_max * u[6]
    )
    
    # Exposed / Incubating
    du[6] = (
      - du[7]
      +h(p, t - te, Val{1}; idxs = 7)
    )
    
    # Infected
    du[5] = (
      -h(p, t - te, Val{1}; idxs = 7)
      +h(p, t - te - ti, Val{1}; idxs = 7) * (1-phi_h)
      +h(p, t - te - th, Val{1}; idxs = 7) * phi_h
      -h(p, t - te - th - thi, Val{1}; idxs = 7) * (phi_h - phi_d(h(p, t - thi + thd; idxs = 4), p))
      +h(p, t - te - th - thi - thii, Val{1}; idxs = 7) * (phi_h - phi_d(h(p, t - thi + thd - thii; idxs = 4), p))
    )

    # Hospitalization 
    du[4] = (
      - h(p, t - te - th, Val{1}; idxs = 7) *  phi_h 
      + h(p, t - te - th - thi, Val{1}; idxs = 7) * (phi_h - phi_d(h(p, t - thi + thd; idxs = 4), p))
      + h(p, t - te - th - thd, Val{1}; idxs = 7) * phi_d(u[4], p)
    )
    
    # Recovered
    du[3] = ( 
      - h(p, t - te - ti, Val{1}; idxs = 7) *  (1-phi_h) 
      - h(p, t - te - th - thi - thii, Val{1}; idxs = 7) * (phi_h - phi_d(h(p, t - thi + thd - thii; idxs = 4), p))
    )
    
    # Deaths
    du[2] = ( 
      - h(p, t - te - th - thd, Val{1}; idxs = 7) *  phi_d(u[4], p)
    )
    
    # Confirmed
    du[1] = du[5] + du[4] + du[3] + du[2]
    
  end")

JuliaCall::julia_eval("@everywhere lags = [te, te + ti, te + th, te + th + thi, te + th + thd, thi - thd, te + th + thi + thii]")



# dx/dt = k * x
# 1/x * dx = k * dt
# integrate
# ln(x) - ln(x_0) = k * t
# x / x_0 = exp(k * t)
# x = x_0 exp(k * t)
# x_0 = 1

h = JuliaCall::julia_eval("@everywhere function h(p, t; idxs::Union{Nothing,Int} = nothing)
    if t > -te
      #if idxs === nothing
      #  u0 * exp(p[1] * t)
      #else
        u0[idxs] * exp(p[1] * t)
      #end
    else
      #if idxs === nothing
      #  [0, 0, 0, 0, 0, 0, n_max]
      #else
        0
      #end
    end
  end")

h = JuliaCall::julia_eval("@everywhere function h(p, t, deriv::Type{Val{1}}; idxs::Union{Nothing,Int} = nothing)
    if t > -te
      - k(p, t) * u0[7] / n_max * u0[6] * (exp(k(p, t - te) * t))
    else
        0
    end
  end")





#####################################################################################
#
# Model parameter fit
#   
#####################################################################################

t_0 <- 40
tk_2 <- tc0_max - t_0 - te


t = c(t_0:len_ger_data)
C = ger_data_confirmed
R <- ger_data_recovered
D <- ger_data_deaths
C[t] <- C[t]*phi_c
R[t] <- R[t]*phi_r
I <- ger_data_confirmed - ger_data_recovered - ger_data_deaths
E <- I * k_raw/(k1 * exp(-k1*te))
I <- I * 0.95
H <- I * 0.05
S <-n_max - I - E - H - R - D

data_df <- data.frame(C = C[t], D = D[t], R = R[t])
data_sd <- abs(data.frame(C = lowess(t,C[t],f=0.1)$y - data_df$C, D = lowess(t,D[t],f=0.3)$y - data_df$D, R = lowess(t,R[t],f=0.3)$y - data_df$R))

JuliaCall::julia_assign("u0", c(C[t_0], D[t_0], R[t_0], H[t_0], I[t_0], E[t_0], S[t_0]))
JuliaCall::julia_assign("tspan", c(0,len_ger_data-t_0+1))
JuliaCall::julia_assign("p", c(k1, tk_2, k2, tk_2+7, k2,50,k2))
JuliaCall::julia_assign("saveat", c(0,t-t_0+1))
JuliaCall::julia_eval("@everywhere u0 = $u0");
JuliaCall::julia_eval("@everywhere tspan = $tspan");
JuliaCall::julia_eval("@everywhere p = $p");
JuliaCall::julia_eval("@everywhere saveat = $saveat");
JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")

data <- array(dim=c(3,length(t),1000))
for (j in 1:length(t))
{
  data[1,j,] <- rnorm(1000,data_df[j,1], data_sd[j,1]+50)
  data[2,j,] <- rnorm(1000,data_df[j,2], data_sd[j,2]+50)
  data[3,j,] <- rnorm(1000,data_df[j,3], data_sd[j,3]+50)
}
JuliaCall::julia_assign("t", t-t_0)
JuliaCall::julia_assign("data", data)
JuliaCall::julia_eval("@everywhere t = $t")
JuliaCall::julia_eval("@everywhere data = $data")
JuliaCall::julia_eval("distributions = [fit_mle(Normal,data[i,j,:]) for i in 1:3, j in 1:length(t)]")
JuliaCall::julia_eval("obj = build_loss_objective(prob, Rodas5(), reltol=1e-4, abstol=1e-7, maxiters = 1e5, LogLikeLoss(t,distributions), verbose=true); nothing")





JuliaCall::julia_eval("bound1 = Tuple{Float64,Float64}[(0.25,0.35),(11.6,11.8),(0.2,0.3),(21.8,22),(0.15,0.2),(35,40),(0.15,0.25),(52,54),(0.15,0.25),(60,65),(0.15,0.25),(95,105),(0.15,0.25)]")
JuliaCall::julia_eval("res1 = bboptimize(obj;SearchRange = bound1, MaxSteps = 11e3, NumDimensions = 15,
    Workers = workers(),
    TraceMode = :compact,
    Method = :adaptive_de_rand_1_bin_radiuslimited)")





p2 <- JuliaCall::julia_eval("p = best_candidate(res1)")
#p2 <- c(0.3165994, 11.7403526, 0.2604021, 21.9828423, 0.1865795, 36.7995045, 0.1781261, 50.0007786, 0.1765100, 61.8112355, 0.1733849, 107.6686829, 0.1827476)
p2
rnames[p2[2]+t_0]
rnames[p2[4]+t_0]
rnames[p2[6]+t_0]
rnames[p2[8]+t_0]
rnames[p2[10]+t_0]
rnames[p2[12]+t_0]

JuliaCall::julia_assign("saveat", c(0:1000/1000*250))
JuliaCall::julia_assign("tspan", c(0,250))
JuliaCall::julia_eval("@everywhere tspan = $tspan")
JuliaCall::julia_eval("@everywhere saveat = $saveat")

#####################################################################################
#
# Output
#   
#####################################################################################
# graphical: 0
# pdf: 1
# png: 2
 plot_out <- 0 
#for (plot_out in c(2:0)) {
  
  if (plot_out == 1) pdf("Plots.pdf")
  
  #####################################################################################
  #
  # Model vs Situation
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
       ylim=c(0,max( phi_c*ger_data_confirmed[c(t_0:len_ger_data),1])),
       col="dark gray",pch=1)
  points(c(t_0:len_ger_data)-t_0,  ger_data_recovered[c(t_0:len_ger_data),1],col="light green",pch=2)
  
  points(c(t_0:len_ger_data)-t_0, phi_c*ger_data_confirmed[c(t_0:len_ger_data),1],col=1,pch=1)
  points(c(t_0:len_ger_data)-t_0, phi_r*ger_data_recovered[c(t_0:len_ger_data),1],col=3,pch=2)
  points(c(t_0:len_ger_data)-t_0, ger_data_deaths[c(t_0:len_ger_data),1],col=2,pch=3)
  lines(sol$t, sol$u[,1], lty=2, col=1)
  lines(sol$t, sol$u[,2], lty=2, col=2)
  lines(sol$t, sol$u[,3], lty=2, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Confirmed cases", "Recovered cases", "Deaths", "Modeled cases"),
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
       ylim=c(0,max( phi_c*ger_data_confirmed[c(t_0:len_ger_data),1]/n_max*100)),
       col="dark gray",pch=1)
  points(c(t_0:len_ger_data)-t_0, ger_data_recovered[c(t_0:len_ger_data),1]/n_max*100,col="light green",pch=2)
  points(c(t_0:len_ger_data)-t_0, phi_c*ger_data_confirmed[c(t_0:len_ger_data),1]/n_max*100,col=1,pch=1)
  points(c(t_0:len_ger_data)-t_0, phi_r*ger_data_recovered[c(t_0:len_ger_data),1]/n_max*100,col=3,pch=2)
  points(c(t_0:len_ger_data)-t_0, ger_data_deaths[c(t_0:len_ger_data),1]/n_max*100,col=2,pch=3)
  
  lines(sol$t, sol$u[,1]/n_max*100, lty=2, col=1)
  lines(sol$t, sol$u[,2]/n_max*100, lty=2, col=2)
  lines(sol$t, sol$u[,3]/n_max*100, lty=2, col=3)
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Confirmed cases", "Recovered cases", "Deaths", "Modeled cases"),
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
       col = "dark gray", 
       lty = 1,
       xlab=paste("Days after", rnames[1]), 
       ylab="Cases",
       ylim=c(0,max(ger_data_confirmed[t_0:len_ger_data] * phi_c)*1.1))
  lines(c(ger_data_confirmed[1:(t_0-1)],ger_data_confirmed[t_0:len_ger_data] * phi_c) , type="l", col = "black", lty = 1)
  lines(ger_data_recovered , type="l", col = "light green", lty = 2)
  lines(c(ger_data_recovered[1:(t_0-1)],ger_data_recovered[t_0:len_ger_data] * phi_c) , type="l", col = "green", lty = 2)
  lines(ger_data_deaths, type="l", col = "red", lty = 3)
  
  for (i in c(0:100)) 
    lines(c(i*7-22,i*7-22), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Confirmed and hidden cases","Confirmed cases (JHU CSSE)", "Recovered and hidden recovered cases", "Recovered cases (JHU CSSE)", "Death cases"), 
         col=c("black","dark gray","green","light green","red"),
         bty="n",
         pch=c(-1,-1),
         lwd=1,
         lty=c(1,1,2,2,3),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  box()
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
       ylab=expression(log['e']('Cases')))
  lines(ln_data_recovered , type="p", col = 2, pch = 2)
  lines(tc0_0:tc0_max,ln_res_1a$fitted.values, col = 1, lty = 1)
  lines(tc1_0:tc1_max,ln_res_2a$fitted.values, col = 1, lty = 2)
  lines(tc2_0:tc2_max,ln_res_3a$fitted.values, col = 1, lty = 3)
  lines(tc3_0:tc3_max,ln_res_4a$fitted.values, col = 1, lty = 4)
  lines(tc4_0:tc4_max,ln_res_5a$fitted.values, col = 1, lty = 5)
  
  lines(c(tc0_0:tc0_max)+dt_r,ln_res_1b$fitted, col = 2, lty = 1)
  lines(c(tc1_0:tc1_max)+dt_r,ln_res_2b$fitted, col = 2, lty = 2)
  lines(c(tc2_0:tc2_max)+dt_r,ln_res_3b$fitted, col = 2, lty = 3)
  lines((tc3_0+dt_r):min(tc3_max+dt_r,len_ger_data),ln_res_4b$fitted, col = 2, lty = 4)

  for (i in c(1:100)) 
    lines(c(i*7-22,i*7-22), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Confirmed cases", "Recovered cases", "1st phase (diffusion / no growth)", "2nd phase (uncontrolled growth)", "3rd phase (social distancing)", "4th phase (lockdown)", "5th phase (exit from lockdown)"), 
         col=c(1,2,16,16,16,16,16,16),
         bty="n",
         pch=c(1,2,-1,-1,-1,-1,-1,-1),
         lwd=1,
         lty=c(-1,-1,1,2,3,4,5,6),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  box()
  if (plot_out != 2) title("Situation COVID-19 in Germany",
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  if (plot_out == 2) png("Situation-3.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  
  R0_plot <- data.frame(t = c(tc0_0:(len_ger_data)), R0 = NA)
  for (i in 1:nrow(R0_plot)) R0_plot$R0[i] <- sum(ger_data_confirmed[tc0_0+i-(1:3),]-ger_data_confirmed[tc0_0+i-(1:3)-1,])/sum(ger_data_confirmed[tc0_0+i-(1:3)-te,]-ger_data_confirmed[tc0_0+i-(1:3)-te-1,])
  R0_plot$t <- R0_plot$t - R0_plot$t[1]
  R0_plot$R0[R0_plot$R0==Inf] = 10
  
  
  plot(c(-1e9,1e9),c(1, 1),lty=3,col=3,
       type = "l", 
       ylim=c(0,6),
       xlim=c(t_0,len_ger_data-tc0_0),
       xlab=paste("Days after", rnames[tc0_0]),
       ylab="Basic Reproductive Number R0")

  data <- data.frame(x=R0_plot$t, y=R0_plot$R0)
  mean <- lowess(data,f=0.2)
  sd <- data.frame(x=data$x,y=(mean$y - data$y)^2)
  sd$y[2:(length(sd$x)-1)] <- (sd$y[1:(length(sd$x)-2)] + sd$y[2:(length(sd$x)-1)] + sd$y[3:(length(sd$x))])/2
  sd$y[c(1,length(sd$x))]<- sd$y[c(1,length(sd$x))] *3/2
  sd <- lowess(sd, f=0.15)
  sd$y <- abs(sd$y)^0.5*qt(0.975,length(sd$y))
  
  polygon(c(mean$x,rev(mean$x)),c(mean$y-sd$y,rev(mean$y+sd$y)),lty=0,col="gray95")
  lines(mean$x,mean$y-sd$y,lty=2)
  lines(mean$x,mean$y+sd$y,lty=2)
  
  lines(c(-1e9,1e9), c(1,1),lty=3,col=3)
  
  points(data)
  lines(mean,lty=2,col=2)
  
  for (i in c(0:100)) 
    lines(c(i*7-22-tc0_0,i*7-22-tc0_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-tc0_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
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
  box()
  if (plot_out != 2) title("Situation COVID-19 in Germany",
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  
  if (plot_out == 2) png("Situation-4.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  R0_plot <- data.frame(t = c(tc0_0:(len_ger_data)), R0 = 0)
  for (i in R0_plot$t) R0_plot$R0[i-R0_plot$t[1]+1] <- lm(ln_x~t, ln_data_confirmed[(i-3):i,])$coefficients[2] + 1/te
  R0_plot$t <- R0_plot$t - R0_plot$t[1] 
  plot(c(-1e9,1e9), c(1/te, 1/te), lty=3,col=3,
       type = "l", 
       ylim=c(0,0.5),
       xlim=c(t_0,len_ger_data-tc0_0),
       xlab=paste("Days after", rnames[tc0_0]),
       ylab=expression('Growth Rate '*(D^-1)))
  
  data <- data.frame(x=R0_plot$t, y=R0_plot$R0)
  mean <- lowess(data,f=0.2)
  sd <- data.frame(x=data$x,y=(mean$y - data$y)^2)
  sd$y[2:(length(sd$x)-1)] <- (sd$y[1:(length(sd$x)-2)] + sd$y[2:(length(sd$x)-1)] + sd$y[3:(length(sd$x))])/2
  sd$y[c(1,length(sd$x))]<- sd$y[c(1,length(sd$x))] *3/2
  sd <- lowess(sd, f=0.15)
  sd$y <- abs(sd$y)^0.5*qt(0.975,length(sd$y))
  
  polygon(c(mean$x,rev(mean$x)),c(mean$y-sd$y,rev(mean$y+sd$y)),lty=0,col="gray95")
  lines(mean$x,mean$y-sd$y,lty=2)
  lines(mean$x,mean$y+sd$y,lty=2)
  
  lines(c(-1e9,1e9), c(1/te,1/te),lty=3,col=3)
  
  points(data)
  lines(mean,lty=2,col=2)
  
  for (i in c(0:100)) 
    lines(c(i*7-22-tc0_0,i*7-22-tc0_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-tc0_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
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
  box()
  if (plot_out != 2) title("Situation COVID-19 in Germany",
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()

  
  
  if (plot_out == 2) png("Situation-5.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)

  ger_data_deaths_day <- lowess(ger_data_deaths[t_0:len_ger_data]-ger_data_deaths[(t_0-1):(len_ger_data-1)],f=0.1)$y
  ger_data_confirmed_day <- lowess(ger_data_confirmed[t_0:len_ger_data]-ger_data_confirmed[(t_0-1):(len_ger_data-1)],f=0.1)$y
  
  ger_data_phi_d <- ger_data_deaths_day[(tcd+1):length(ger_data_deaths_day)]/ger_data_confirmed_day[1:(length(ger_data_deaths_day)-tcd)]
  
  #plot(1:length(ger_data_deaths_day),ger_data_confirmed_day)
  #points((1:length(ger_data_deaths_day))-tcd,ger_data_deaths_day/0.043,col=2)
  
  plot(c(-1e9,1e9), c(phi_d_median,phi_d_median)*100,
       lty=3,col=3,
       type = "l", 
       ylim=c(2,8),
       xlim=c((t_0+tcd),len_ger_data)-t_0,
       xlab=paste("Days after", rnames[t_0]),
       ylab="Death Cases (%)")
  data <- data.frame(x=c((t_0+tcd):len_ger_data)-t_0, y=ger_data_phi_d*100)
  mean <- lowess(data,f=0.3)
  sd <- lowess(data$x,(mean$y - data$y)^2, f=0.15)
  sd$y <- abs(sd$y)^0.5*3/2
  
  polygon(c(mean$x,rev(mean$x)),c(mean$y-sd$y,rev(mean$y+sd$y)),lty=0,col="gray95")
  lines(mean$x,mean$y-sd$y,lty=2)
  lines(mean$x,mean$y+sd$y,lty=2)
  
  lines(c(-1e9,1e9), c(phi_d_median,phi_d_median)*100,lty=3,col=3)
  
  points(data)
  lines(mean,lty=2,col=2)
  
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  box()
  if (plot_out != 2) title("Situation COVID-19 in Germany",
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  
  
  if (plot_out == 2) png("Situation-6.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  
  ger_data_phi_h <- ger_data_phi_d / phi_d_median 

  plot(c(-1e9,1e9), c(1,1)*100,
       type="l",
       lty=3,col=3,
       ylim=c(0,250),
       xlim=c((t_0+tcd),len_ger_data)-t_0,
       xlab=paste("Days after", rnames[t_0]),
       ylab="Hidden Cases Per Day (%)")
  
  data <- data.frame(x=c((t_0+tcd):len_ger_data)-t_0, y=ger_data_phi_h*100)
  mean <- lowess(data,f=0.3)
  sd <- data.frame(x=data$x,y=(mean$y - data$y)^2)
  sd$y[2:(length(sd$x)-1)] <- (sd$y[1:(length(sd$x)-2)] + sd$y[2:(length(sd$x)-1)] + sd$y[3:(length(sd$x))])/2
  sd$y[c(1,length(sd$x))]<- sd$y[c(1,length(sd$x))] *3/2
  sd <- lowess(sd, f=0.25)
  sd$y <- abs(sd$y)^0.5*qt(0.975,length(sd$y))
  
  polygon(c(mean$x,rev(mean$x)),c(mean$y-sd$y,rev(mean$y+sd$y)),lty=0,col="gray95")
  lines(mean$x,mean$y-sd$y,lty=2)
  lines(mean$x,mean$y+sd$y,lty=2)
  
  lines(c(-1e9,1e9), c(1,1)*100,lty=3,col=3)
  
  points(data)
  lines(mean,lty=2,col=2)
  
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  box()
  if (plot_out != 2) title("Situation COVID-19 in Germany",
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  
  if (plot_out == 2) png("Situation-7.png", width = 640, height = 480)
  par(mar=c(5,6,7,5)+0.1)
  
  phi_c_plot <- phi_d / phi_d_median
  
  plot(c(-1e9,1e9), c(1,1)*100,
       type="l",
       lty=3,col=3,
       ylim=c(80,140),
       xlim=c(t_0+tcd,len_ger_data)-t_0,
       xlab=paste("Days after", rnames[t_0]),
       ylab="Hidden Cases Cumulative (%)")
  
  data <- data.frame(x=c((t_0):len_ger_data)-t_0, y=phi_c_plot*100)
  data <- data[1:length(data$x),]
  mean <- lowess(data,f=0.3)
  mean$y[1:22]<- mean$y[22]
  sd <- data.frame(x=data$x,y=(mean$y - data$y)^2)
  sd$y[2:(length(sd$x)-1)] <- (sd$y[1:(length(sd$x)-2)] + sd$y[2:(length(sd$x)-1)] + sd$y[3:(length(sd$x))])/2
  sd$y[c(1,length(sd$x))]<- sd$y[c(1,length(sd$x))] *3/2
  sd <- lowess(sd, f=0.25)
  sd$y <- abs(sd$y)^0.5*qt(0.975,length(sd$y))
  
  polygon(c(mean$x,rev(mean$x)),c(mean$y-sd$y,rev(mean$y+sd$y)),lty=0,col="gray95")
  lines(mean$x,mean$y-sd$y,lty=2)
  lines(mean$x,mean$y+sd$y,lty=2)
  
  lines(c(-1e9,1e9), c(1,1)*100,lty=3,col=3)
  
  points(data)
  mean$y[1:22]<-NA
  lines(mean,lty=2,col=2)
  
  for (i in c(1:100)) 
    lines(c(i*7-22-t_0,i*7-22-t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray" )
  axis(3, c(1:100)*7-22-t_0, c(1:100)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
  mtext("Week number (2020)",side=3,line=2,las=0)
  box()
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
       ylim=c(0,3e5), 
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
  mtext("Hospital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hospital Workload"), 
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
  mtext("Hospital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hospital Workload"), 
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
  # Scenario 1 no social distancing and no lockdown
  #   
  #####################################################################################
  
  JuliaCall::julia_assign("p", c(p2[1], 1000, p2[1], 1001, p2[1], 1002, p2[1], 1003, p2[1], 1004, p2[1], 1005, p2[1]))
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
  if (plot_out != 2) title("Scenario 1 no social distancing and no lockdown", 
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
  if (plot_out != 2) title("Scenario 1 no social distancing / no lockdown ", 
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
  mtext("Hospital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hospital Workload"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 1 no social distancing and no lockdown", 
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
  mtext("Hospital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hospital Workload"), 
         col=c(1:6,16),
         bty="n",
         lwd=1,
         lty=c(1:7),
         cex = 0.8,
         x.intersp = 2.5,
         ncol=1)
  if (plot_out != 2) title("Scenario 1 no social distancing and no lockdown", 
                           sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
  if (plot_out > 1) dev.off()
  
  #####################################################################################
  #
  # Scenario 2 exit after day 50 without social distancing (worst case)
  #   
  #####################################################################################
  
  JuliaCall::julia_assign("p", c(p2[1], p2[2], p2[3], p2[4], p2[5], p2[6], p2[7], 50, p2[1], 80, p2[1], 100, p2[1]))
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
  mtext("Hospital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hospital Workload"), 
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
  mtext("Hospital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hospital Workload"), 
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
  
  JuliaCall::julia_assign("p", c(p2[1], p2[2], p2[3], p2[4], p2[5], p2[6], p2[7], 50, p2[3], 80, p2[3], 100, p2[3]))
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
  mtext("Hospital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hospital Workload"), 
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
  mtext("Hospital Workload compared to 2017 (%)", side=4,line=3,las=0)
  mtext("Week number (2020)",side=3,line=2,las=0)
  legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hospital Workload"), 
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
  
  JuliaCall::julia_assign("p", c(p2[1], p2[2], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1]))
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
    
    p_var <- c(p2[1], tk_2, k2, tk_2+1, k2, tk_2+2, k2, tk_2+3, k2, tk_2+4, k2, tk_2+5, k2)
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
      mtext(expression('Growth Rate '*(D^-1)), side=4,line=3,las=0)
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
    
    p_var <- c(p2[1], tk_2, k2, tk_2+1, k2, tk_2+2, k2, tk_2+3, k2, tk_2+4, k2, tk_2+5, k2)
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
      mtext(expression('Growth Rate '*(D^-1)), side=4,line=3,las=0)
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
    
    
    p_var <- c(p2[1], tk_2, k2, tk_2+1, k2, tk_2+2, k2, tk_2+3, k2, tk_2+4, k2, tk_2+5, k2)
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
      mtext(expression('Growth Rate '*(D^-1)), side=4,line=2.5,las=0)
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
#}

system("git add *")
system(paste("git commit -m \"Update Data ", format(Sys.time(), "%m/%d/%Y %H:%M"),"\"", sep=""))
system("git push")

