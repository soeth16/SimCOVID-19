if (plot_out == 1) pdf("Plots-Variants.pdf")

#####################################################################################
#
# variant calculation
#   
#####################################################################################

JuliaCall::julia_assign("p", c(p2[1], p2[2], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1], p2[3], p2[1]))
JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
JuliaCall::julia_eval("sol = solve(prob, MethodOfSteps(Tsit5()), reltol=1e-4, abstol=1e-9, maxiters = 10000, saveat=saveat); nothing")
sol1 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))

rbcol = rainbow(11)

if (plot_out == 2) png("Forecast-ARDS-3.png", width = 640, height = 480)

par(mar=c(5,6,7,5)+0.1)

plot(sol1$t+t_0, (sol1$u[,4]), 
     type="l",
     axes=F, xlab="", ylab="", 
     ylim=c(-0.2,1)*max((sol1$u[,2])),
     xlim=c(t_0,200), 
     lty=1, col=1)
lines(sol1$t+t_0, (sol1$u[,2]), type="l", lty=2, col=1)

axis(2, pretty(range(sol1$u[,2]),10), col="black",las=1)  ## las=1 makes horizontal labels
mtext("Cases",side=2,line=4.5,las=0)
axis(1,pretty(range(t_0:200),5))
mtext(paste("Days after", format(start_date, "%m/%d/%Y")),side=1,col="black",line=2.5)


for (i in c(0:10))
{
  tk_2 <- 22
  R0_2 <- 0.4 + i * 0.2
  f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - log(R0_2 + 1) / te
  k2 <- uniroot(f_k2 , c(0.01, 1))$root
  
  p_var <- c(p2[1], tk_2, k2, tk_2+1, k2, tk_2+2, k2, tk_2+3, k2, tk_2+4, k2, tk_2+5, k2, tk_2+6, k2, tk_2+7, k2,tk_2+7, k2, tk_2+7, k2, tk_2+7, k2, tk_2+7, k2)
  JuliaCall::julia_assign("p", p_var)
  JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
  JuliaCall::julia_eval("sol = solve(prob, MethodOfSteps(Tsit5()), reltol=1e-4, abstol=1e-9, maxiters = 10000, saveat=saveat); nothing")
  sol2 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
  
  par(new=1)
  plot(sol2$t+t_0, (sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max((sol1$u[,2])), xlim=c(t_0,200))
  lines(sol2$t+t_0, (sol2$u[,2]), type="l", lty=2, col=rbcol[i])
  
  par(new=1)
  plot(sol2$t+t_0, k(p_var, sol2$t), ylim = c(0.2,1.5), xlim=c(t_0,200), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
  if(i==0)
  {
    axis(4, pretty(range(k(p_var, sol2$t)),10))
    mtext(expression('Growth Rate '*(D^-1)), side=4,line=3,las=0,adj=0)
  }
}

par(new=1)
plot(c(-1e9,1e9), c(n_h_max, n_h_max), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max((sol1$u[,2])), xlim=c(t_0,200))

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:53,1:47)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

lines(c(22+t_0,22+t_0), (c(-1e9, 1e9)), type="l", lty = 5, col="grey" )

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

plot(sol1$t+t_0, (sol1$u[,4]), 
     type="l",
     axes=F, xlab="", ylab="", 
     ylim=c(-0.2,1)*max((sol1$u[,2])),
     xlim=c(t_0,200), 
     lty=1, col=1)
lines(sol1$t+t_0, (sol1$u[,2]), type="l", lty=2, col=1)

axis(2, pretty(range(sol1$u[,2]),10), col="black",las=1)  ## las=1 makes horizontal labels
mtext("Cases",side=2,line=4.5,las=0)
axis(1,pretty(range(t_0:200),5))
mtext(paste("Days after", format(start_date, "%m/%d/%Y")),side=1,col="black",line=2.5)


for (i in c(0:10))
{
  tk_2 <- 16 + i * 3
  R0_2 <- 1
  f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - log(R0_2 + 1) / te
  k2 <- uniroot(f_k2 , c(0.01, 1))$root
  
  p_var <- c(p2[1], tk_2, k2, tk_2+1, k2, tk_2+2, k2, tk_2+3, k2, tk_2+4, k2, tk_2+5, k2, tk_2+6, k2, tk_2+7, k2, tk_2+7, k2, tk_2+7, k2, tk_2+7, k2, tk_2+7, k2)
  JuliaCall::julia_assign("p", p_var)
  JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
  JuliaCall::julia_eval("sol = solve(prob, MethodOfSteps(Tsit5()), reltol=1e-4, abstol=1e-9, maxiters = 10000, saveat=saveat); nothing")
  sol2 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
  
  par(new=1)
  plot(sol2$t+t_0, (sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max((sol1$u[,2])), xlim=c(t_0,200))
  lines(sol2$t+t_0, (sol2$u[,2]), type="l", lty=2, col=rbcol[i])
  
  par(new=1)
  plot(sol2$t+t_0, k(p_var, sol2$t), ylim = c(0.2,1.5), xlim=c(t_0,200), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
  if(i==0)
  {
    axis(4, pretty(range(k(p_var, sol2$t)),10))
    mtext(expression('Growth Rate '*(D^-1)), side=4,line=3,las=0,adj=0)
  }
}

par(new=1)
plot(c(0,250), (c(n_h_max, n_h_max)), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(0,1)*max(log(sol1$u[,2])), xlim=c(t_0,200))

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:53,1:47)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

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

plot(sol1$t+t_0, log(sol1$u[,4]), 
     type="l",
     axes=F, xlab="", ylab="", 
     ylim=c(-0.2,1)*max(log(sol1$u[,6])),
     xlim=c(t_0,200), 
     lty=1, col=1)
lines(sol1$t+t_0, log(sol1$u[,2]), type="l", lty=2, col=1)

axis(2, pretty(range(log(sol1$u[,6])),10), col="black",las=1)  ## las=1 makes horizontal labels
mtext(expression(log['e']('Cases')),side=2,line=4,las=0)
axis(1,pretty(range(t_0:200),5))
mtext(paste("Days after", format(start_date, "%m/%d/%Y")),side=1,col="black",line=2.5)


for (i in c(0:10))
{
  tk_2 <- 16 + i * 3
  R0_2 <- 1
  f_k2 <- function(k2) k2 *(1-exp(-k2 * te)) - log(R0_2 + 1) / te
  k2 <- uniroot(f_k2 , c(0.01, 1))$root
  
  
  p_var <- c(p2[1], tk_2, k2, tk_2+1, k2, tk_2+2, k2, tk_2+3, k2, tk_2+4, k2, tk_2+5, k2, tk_2+6, k2, tk_2+6, k2, tk_2+7, k2, tk_2+7, k2, tk_2+7, k2, tk_2+7, k2)
  JuliaCall::julia_assign("p", p_var)
  JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")
  JuliaCall::julia_eval("sol = solve(prob, MethodOfSteps(Tsit5()), reltol=1e-8, abstol=1e-9, maxiters = 10000, saveat=saveat); nothing")
  sol2 = list(t=JuliaCall::julia_eval("sol.t"), u=JuliaCall::julia_eval("typeof(u0)<:Number ? Array(sol) : sol'"))
  
  par(new=1)
  sol2$u[log(sol2$u[,4])<exp(1),4] <- 1 
  sol2$u[log(sol2$u[,2])<exp(1),6] <- 1 
  plot(sol2$t+t_0, log(sol2$u[,4]), type="l", lty=1, col=rbcol[i], axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max(log(sol1$u[,6])),xlim=c(t_0,200))
  lines(sol2$t+t_0, log(sol2$u[,2]), type="l", lty=2, col=rbcol[i])
  
  par(new=1)
  plot(sol2$t+t_0, k(p_var, sol2$t), ylim = c(0.2,1.5), xlim=c(t_0,200), type="l", lty=3, col=rbcol[i], axes=F, xlab="", ylab="")
  if(i==0)
  {
    axis(4, pretty(range(k(p_var, sol2$t)),10))
    mtext(expression('Growth Rate '*(D^-1)), side=4,line=2.5,las=0,adj=0)
  }
}

par(new=1)
plot(c(0,250), log(c(n_h_max, n_h_max)), type="l", lty = 4, col="grey", axes=F, xlab="", ylab="", ylim=c(-0.2,1)*max(log(sol1$u[,6])), xlim=c(t_0,200))

for (i in c(0:100)) lines(c(i*7-2,i*7-2), (c(-1e9, 1e9)), type="l", lty = 5, col="light gray")
axis(3, c(0:100)*7-2, c(0:53,1:47)+1, col="light gray", las=0)  ## las=1 makes horizontal labels
mtext("Week number (2020/21)",side=3,line=2,las=0)

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


if (plot_out == 1) dev.off()