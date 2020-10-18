if (plot_out == 1) pdf("Plots-Model_vs_Situation3.pdf")

#####################################################################################
#
# Model vs Situation
#   
#####################################################################################

JuliaCall::julia_assign("p", p2)
JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
JuliaCall::julia_eval("sol = solve(prob,  MethodOfSteps(Tsit5()), reltol=1e-4, abstol=1e-9, maxiters = 10000, saveat=saveat); nothing")
sol = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))

if (plot_out == 2) png("Model_vs_Situation-1.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)
plot(cumulative[[8]]$Time, cumulative[[8]]$ConfirmedCases, 
     xlab=paste("Days after",format(start_date, "%m/%d/%Y")), 
     ylab="Confirmed cases",
     ylim=c(0,max(cumulative_extrapolated[[8]]$ConfirmedCases)),
     col="dark gray",pch=1)
points(cumulative[[8]]$Time, cumulative[[8]]$RecoverdCases, col="light green",pch=2)
points(cumulative_extrapolated[[8]]$Time, cumulative_extrapolated[[8]]$ConfirmedCases, col=1, pch=1)
points(cumulative_extrapolated[[8]]$Time, cumulative_extrapolated[[8]]$RecoverdCases, col=3, pch=2)
points(cumulative_extrapolated[[8]]$Time, cumulative_extrapolated[[8]]$DeathCases, col=2, pch=3)
lines(sol$t+t_0, sol$u[,1], lty=2, col=1)
lines(sol$t+t_0, sol$u[,2], lty=2, col=2)
lines(sol$t+t_0, sol$u[,3], lty=2, col=3)

lk <- (length(p2)-1)/2
tk <- p2[(1:lk)*2]
for (i in tk) lines(c(i,i)+t_0, (c(-1e9, 1e9)), type="l", lty = 8, col="blue")

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:53,1:47)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topleft", legend <- c("Confirmed cases", "Recovered cases", "Deaths", "Modeled cases"),
       col=c(1,3,2,1),
       bty="n",
       pch=c(1,2,3,-1),
       lwd=1,
       lty=c(-1,-1,-1,2),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)

box()

if (plot_out != 2) title("Situation vs Model COVID-19 in Germany",
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()


if (plot_out == 2) png("Model_vs_Situation-2.png", width = 640, height = 480)
par(mar=c(5,6,7,5)+0.1)
plot(cumulative[[8]]$Time, cumulative[[8]]$ConfirmedCases/n_max*100, 
     xlab=paste("Days after", format(start_date, "%m/%d/%Y")), 
     ylab="Cases (%)",
     ylim=c(0,max(cumulative_extrapolated[[8]]$ConfirmedCases/n_max*100)),
     col="dark gray",pch=1)
points(cumulative[[8]]$Time, cumulative[[8]]$RecoverdCases/n_max*100, col="light green",pch=2)
points(cumulative_extrapolated[[8]]$Time, cumulative_extrapolated[[8]]$ConfirmedCases/n_max*100, col=1, pch=1)
points(cumulative_extrapolated[[8]]$Time, cumulative_extrapolated[[8]]$RecoverdCases/n_max*100, col=3, pch=2)
points(cumulative_extrapolated[[8]]$Time, cumulative_extrapolated[[8]]$DeathCases/n_max*100, col=2, pch=3)

lines(sol$t+t_0, sol$u[,1]/n_max*100, lty=2, col=1)
lines(sol$t+t_0, sol$u[,2]/n_max*100, lty=2, col=2)
lines(sol$t+t_0, sol$u[,3]/n_max*100, lty=2, col=3)

lk <- (length(p2)-1)/2
tk <- p2[(1:lk)*2]
for (i in tk) lines(c(i,i)+t_0, (c(-1e9, 1e9)), type="l", lty = 8, col="blue")

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:53,1:47)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

legend("topleft", legend <- c("Confirmed cases", "Recovered cases", "Deaths", "Modeled cases"),
       col=c(1,3,2,1),
       bty="n",
       pch=c(1,2,3,-1),
       lwd=1,
       lty=c(-1,-1,-1,2),
       cex = 0.8,
       x.intersp = 2.5,
       ncol=1)

box()

if (plot_out != 2) title("Situation vs Model COVID-19 in Germany",
                         sub=paste("Created by Sören Thiering (",format(Sys.Date(), "%m/%d/%Y"),"). Email: soeren.thiering@hs-anhalt.de",sep=""))
if (plot_out > 1) dev.off()


if (plot_out == 1) dev.off()