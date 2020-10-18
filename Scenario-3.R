if (plot_out == 1) pdf("Plots-Scenario-3.pdf")

#####################################################################################
#
# Scenario 3 exit after day 50 with social distancing (better case)
#   
#####################################################################################

JuliaCall::julia_assign("p", c(p2[1], p2[2], p2[3], p2[4], p2[5], p2[6], p2[7], 50, p2[3], 80, p2[3], 100, p2[3], 130, p2[3], 150, p2[3], 150, p2[3], 150, p2[3], 150, p2[3], 150, p2[3]))
JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
JuliaCall::julia_eval("sol = solve(prob,  MethodOfSteps(Tsit5()), reltol=1e-4, abstol=1e-9, maxiters = 100000, saveat=saveat); nothing")
sol = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))


if (plot_out == 2) png("Forecast-1-scenario-3.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)
plot(sol$t+t_0,sol$u[,7],
     type="l", 
     xlab=paste("Days after", format(start_date, "%m/%d/%Y")), 
     ylab="Cases", 
     ylim=c(0,n_max), 
     xlim=c(t_0,400), 
     lty=1, col=1)
lines(sol$t+t_0, sol$u[,6], lty=2, col=2)
lines(sol$t+t_0, sol$u[,5], lty=3, col=3)
lines(sol$t+t_0, sol$u[,4], lty=4, col=4)
lines(sol$t+t_0, sol$u[,3], lty=5, col=5)
lines(sol$t+t_0, sol$u[,2], lty=6, col=6)
lines(sol$t+t_0, sol$u[,1], lty=7, col=16)

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:53,1:47)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topright", legend <- c("Susceptibles (Noninfected)", "Exposed (Incubation Time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Confirmed Cases (Total Infected)"),       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)

box()

if (plot_out != 2) title("Scenario 3 exit after day 50 with social distancing", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()



if (plot_out == 2) png("Forecast-2-scenario-3.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)
plot(sol$t+t_0,sol$u[,7]/n_max*100,
     type="l", 
     xlab=paste("Days after", format(start_date, "%m/%d/%Y")),
     ylab="Cases (%)",
     ylim=c(0,100),
     xlim=c(t_0,400), 
     lty=1, col=1)
lines(sol$t+t_0, sol$u[,6]/n_max*100, lty=2, col=2)
lines(sol$t+t_0, sol$u[,5]/n_max*100, lty=3, col=3)
lines(sol$t+t_0, sol$u[,4]/n_max*100, lty=4, col=4)
lines(sol$t+t_0, sol$u[,3]/n_max*100, lty=5, col=5)
lines(sol$t+t_0, sol$u[,2]/n_max*100, lty=6, col=6)
lines(sol$t+t_0, sol$u[,1]/n_max*100, lty=7, col=16)

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:53,1:47)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topright", legend <- c("Susceptibles (Noninfected)", "Exposed (Incubation Time)", "Infected", "Hospitalized (ARDS)", "Recovered", "Deaths","Confirmed Cases (Total Infected)"), 
       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)

box()

if (plot_out != 2) title("Scenario 3 exit after day 50 with social distancing", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()


max(sol$u[,4])
max(sol$u[,4])/n_max*100
max(sol$u[,4])/n_h_max


if (plot_out == 2) png("Forecast-ARDS-1-scenario-3.png", width = 640, height = 480)

par(mar=c(5,6,7,5)+0.1)

plot(sol$t+t_0, sol$u[,4], 
     type="l",
     xlab=paste("Days after", format(start_date, "%m/%d/%Y")), 
     ylab="Cases", 
     ylim=c(0,max(sol$u[,2],sol$u[,4])), 
     xlim=c(t_0,400), 
     lty=1, col=1)
lines(sol$t+t_0, sol$u[,2], type="l", lty=2, col=2)
par(new=1)
plot(sol$t+t_0, sol$u[,4] / n_h_max * 100, 
     type="l",
     xlab="", 
     ylab="",
     axes = F,
     xlim=c(t_0,400), 
     lty=3, col=3)

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:53,1:47)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

axis(4, pretty(sol$u[,4] / n_h_max * 100,5))
mtext("Hospital Workload compared to 2017 (%)", side=4,line=3,las=0)

legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hospital Workload"), 
       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)

box()

if (plot_out != 2) title("Scenario 3 exit after day 50 with social distancing", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()

if (plot_out == 2) png("Forecast-ARDS-2-scenario-3.png", width = 640, height = 480)

par(mar=c(5,6,7,5)+0.1)

plot(sol$t+t_0, sol$u[,4]/n_max*100, 
     type="l",
     xlab=paste("Days after", format(start_date, "%m/%d/%Y")), 
     ylab="Cases (%)", 
     ylim=c(0,max(sol$u[,2],sol$u[,4])/n_max*100), 
     xlim=c(t_0,400), 
     lty=1, col=1)
lines(sol$t+t_0, sol$u[,2]/n_max*100, type="l", lty=2, col=2)
par(new=1)
plot(sol$t+t_0, sol$u[,4] / n_h_max * 100, 
     type="l",
     xlab="", 
     ylab="",
     axes = F,
     xlim=c(t_0,400), 
     lty=3, col=3)

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:53,1:47)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

axis(4, pretty(sol$u[,4] / n_h_max * 100,5))
mtext("Hospital Workload compared to 2017 (%)", side=4,line=3,las=0)

legend("topleft", legend <- c("Hospitalized (ARDS)", "Deaths", "Hospital Workload"), 
       col=c(1:6,16),
       bty="n",
       lwd=1,
       lty=c(1:7),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)

box()

if (plot_out != 2) title("Scenario 3 exit after day 50 with social distancing", 
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()


if (plot_out == 1) dev.off()