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


# output
#
# grapfical: 0
# pdf: 1
# png: 2
plot_out <- 0
#for (plot_out in c(0)) {

#if (plot_out == 1) pdf("Plots.pdf")


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
tc_max <- 60
tc2_0 <- 59
tc2_max <- len_ger_data
tr_0 <- 35
tr_max <- len_ger_data
dt_r = 14

ln_res_1a = lm(ln_x~t, ln_data_confirmed[15:35,])
ln_res_2a = lm(ln_x~t, ln_data_confirmed[tc_0:tc_max,])
ln_res_3a = lm(ln_x~t, ln_data_confirmed[tc2_0:tc2_max,])
ln_res_1b = lm(ln_x~t, ln_data_recovered[(14+dt_r):(35+dt_r),])
ln_res_2b = lm(ln_x~t, ln_data_recovered[(tr_0+dt_r):tr_max,])



if (plot_out == 2) png("Situation-2.png", width = 640, height = 480)
plot(ln_data_confirmed,
     type="p", 
     col = 1, 
     pch = 1,
     xlab=paste("Days after", rnames[1]), 
     ylab="ln( cases )")
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




# assume incubation time from https://www.ncbi.nlm.nih.gov/pubmed/32150748
te = 5.1
ti = 14

# assume k_raw from ln(x)-plot
k_raw = ln_res_2a$coefficients[2]

# correct k_raw by the influece of te, ti, 

f_k <- function(k) k *(1-exp(-k * ti)) - k_raw
k <- uniroot(f_k , c(0.01, 1))$root

# assume population of germany from https://de.wikipedia.org/wiki/Deutschland (03/13/20)
x_max <- 83.019213e6


JuliaCall::julia_assign("x_max", x_max)
JuliaCall::julia_assign("te", te)
JuliaCall::julia_assign("ti", ti)
f = JuliaCall::julia_eval("function f(du, u, p, t)
    # parameter
    # k     = p[1]  # growth rate
    # tk_2  = p[2]  # time of change from k to k2
    # k2    = p[3]  # k after change
    if (t >= p[2]) 
      k = p[3]
    else
      k = p[1]
    end
    #ti     = p[4]  # infection time

    # Susceptibles
    du[5] = - k * u[5] / x_max * u[4]
    
    # Exposed
    du[4] = + k * u[5] / x_max * u[4] - 1/te * u[4]

    # Infected
    du[3] = + 1/te * u[4] - 1/ti * u[3]
    
    # Recovered
    du[2] = + 1/ti * u[3]

    # Confirmed
    du[1] = du[2] + du[3]

  end")



t_0 <- 40
tk_2 <- tc_max - t_0 - te
t = c(t_0:len_ger_data)
C = ger_data_confirmed
I = ger_data_confirmed - ger_data_recovered - ger_data_deaths
R = ger_data_recovered + ger_data_deaths
E <- I * k_raw/(k * exp(-k*te))
S = x_max - E - I - R

data_df <- data.frame(C[t], R[t], I[t])
JuliaCall::julia_assign("u0", c(C[t_0], R[t_0], I[t_0], E[t_0], S[t_0]))
JuliaCall::julia_assign("tspan", c(0,250))
JuliaCall::julia_assign("p", c(k,tk_2,k2))
JuliaCall::julia_assign("saveat", c(0:1000/1000*250))
JuliaCall::julia_eval("prob = ODEProblem(f,u0,tspan,p)")

# LogLikeLoss
data <- array(dim=c(3,length(t),100))
for (j in 1:length(t))
{
  data[1,j,] <- rnorm(100,data_df[j,1],100)
  data[2,j,] <- rnorm(100,data_df[j,2],100)
  data[3,j,] <- rnorm(100,data_df[j,3],100)
}
JuliaCall::julia_assign("t", t-t_0)
JuliaCall::julia_assign("data", data)
JuliaCall::julia_eval("distributions = [fit_mle(Normal,data[i,j,:]) for i in 1:3, j in 1:length(t)]")
JuliaCall::julia_eval("obj = build_loss_objective(prob,Tsit5(),LogLikeLoss(t,distributions),maxiters=10000,verbose=false); nothing")

# L2Loss
#data <- array(dim=c(3,length(t),1))
#for (j in 1:length(t))
#{
#  data[1,j,] <- data_df[j,1]
#  data[2,j,] <- data_df[j,2]
#  data[3,j,] <- data_df[j,3]
#}
#JuliaCall::julia_assign("t", t-t_0)
#JuliaCall::julia_assign("data", data)
#JuliaCall::julia_eval("obj = build_loss_objective(prob,Tsit5(),L2Loss(t,data), maxiters=10000,verbose=false); nothing")

JuliaCall::julia_eval("bound1 = Tuple{Float64, Float64}[(0.1, 2),(10, 25),(0.1, 2)]")
JuliaCall::julia_eval("res = bboptimize(obj;SearchRange = bound1, MaxSteps = 11e4)")
p2 <- JuliaCall::julia_eval("p = best_candidate(res)")

#JuliaCall::julia_eval("res = Optim.optimize(obj, [0.3,20,0.15,14], Optim.BFGS())")
#p2 <- JuliaCall::julia_eval("p = res.minimizer")

JuliaCall::julia_eval("prob = ODEProblem(f,u0,tspan,p)")
JuliaCall::julia_eval("sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-9, maxiters = 1000000, saveat=saveat); nothing")
sol = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))

if (plot_out == 2) png("Model_vs_Situation-1.png", width = 640, height = 480)
plot(t-t_0,C[t],
     pch=1,
     col=1,
     ylim=c(0,max(data[1,,])),
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases")
points(t-t_0, R[t], pch=2, col=2)
points(t-t_0, I[t], pch=3, col=3)
points(t-t_0, S[t], pch=4, col=4)

lines(sol$t, sol$u[,1], lty=1, col=1)
lines(sol$t, sol$u[,2], lty=2, col=2)
lines(sol$t, sol$u[,3], lty=3, col=3)
lines(sol$t, sol$u[,4], lty=4, col=4)
lines(sol$t, sol$u[,5], lty=4, col=5)

legend("topleft", legend <- c("Confirmed", "Recovered", "Infected", "Exposed", "Susceptibles"), 
       col=c(1:5),
       bty="n",
       pch=c(1:5),
       lwd=1,
       lty=c(1:5),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation vs Model COVID-19 in Germany",
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


if (plot_out == 2) png("Model_vs_Situation-2.png", width = 640, height = 480)

plot(t-t_0,C[t]/x_max*100,
     pch=1,
     col=1,
     ylim=c(0,max(data[1,,]/x_max*100)),
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases (%)")
points(t-t_0, R[t]/x_max*100, pch=2, col=2)
points(t-t_0, I[t]/x_max*100, pch=3, col=3)
points(t-t_0, S[t]/x_max*100, pch=4, col=4)

lines(sol$t, sol$u[,1]/x_max*100, lty=1, col=1)
lines(sol$t, sol$u[,2]/x_max*100, lty=2, col=2)
lines(sol$t, sol$u[,3]/x_max*100, lty=3, col=3)
lines(sol$t, sol$u[,5]/x_max*100, lty=4, col=4)
lines(sol$t, sol$u[,4]/x_max*100, lty=4, col=5)

legend("topleft", legend <- c("Confirmed", "Recovered", "Infected", "Susceptibles", "Exposed"), 
       col=c(1:4),
       bty="n",
       pch=c(1:4),
       lwd=1,
       lty=c(1:4),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)
if (plot_out != 2) title("Situation vs Model COVID-19 in Germany",
                         sub="Created by Sören Thiering 03/27/2020. Email: soeren.thiering@hs-anhalt.de")
if (plot_out > 1) dev.off()


if (plot_out == 2) png("Forecast-1.png", width = 640, height = 480)

plot(sol$t,sol$u[,5],
     type="l", 
     xlab=paste("Days after", rnames[t_0]), 
     ylab="Cases", 
     ylim=c(0,x_max), 
     xlim=c(0,250), 
     lty=1, col=1)
lines(sol$t, sol$u[,4], lty=2, col=2)
lines(sol$t, sol$u[,3], lty=3, col=3)
lines(sol$t, sol$u[,2], lty=4, col=4)
lines(sol$t, sol$u[,1], lty=5, col=5)
legend("topright", legend <- c("Susceptibles (noninfected)", "Exposed (incubation time)", "Infected", "Recovered", "Confirmed"),
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


if (plot_out == 2) png("Forecast-2.png", width = 640, height = 480)
plot(sol$t,sol$u[,5]/x_max*100,
     type="l", 
     xlab=paste("Days after", rnames[t_0]),
     ylab="Cases (%)",
     ylim=c(0,100),
     xlim=c(0,250), 
     lty=1, col=1)
lines(sol$t, sol$u[,4]/x_max*100, lty=2, col=2)
lines(sol$t, sol$u[,3]/x_max*100, lty=3, col=3)
lines(sol$t, sol$u[,2]/x_max*100, lty=4, col=4)
lines(sol$t, sol$u[,1]/x_max*100, lty=5, col=5)
legend("topright", legend <- c("Susceptibles (noninfected)", "Exposed (incubation time)", "Infected", "Recovered", "Confirmed"),
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



# the end 
#if (plot_out == 1) dev.off()
#}

