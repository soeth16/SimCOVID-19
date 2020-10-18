#####################################################################################
#
# Model parameter fit
#   
#####################################################################################
JuliaCall::julia_assign("tspan", c(t[1],t[length(t)]))
JuliaCall::julia_assign("saveat", t)
JuliaCall::julia_eval("prob = DDEProblem(f,u0,h,tspan,p,constant_lags=lags)")

data <- array(dim=c(3,length(t),1000))
for (j in 1:length(t))
{
  data[1,j,] <- rnorm(1000,data_df[j,1], (data_sd[j,1]+50))
  data[2,j,] <- rnorm(1000,data_df[j,2], (data_sd[j,2]+50)*1000)
  data[3,j,] <- rnorm(1000,data_df[j,3], (data_sd[j,3]+50)*1004)
}
JuliaCall::julia_assign("t", t)
JuliaCall::julia_assign("data", data)
JuliaCall::julia_eval("distributions = [fit_mle(Normal,data[i,j,:]) for i in 1:1, j in 1:length(t)]")
JuliaCall::julia_eval("obj = build_loss_objective(prob, MethodOfSteps(Tsit5()), reltol=1e-4, abstol=1e-9, maxiters = 1e4, LogLikeLoss(t,distributions), verbose=true); nothing")


JuliaCall::julia_eval("bound1 = Tuple{Float64,Float64}[(0.3,0.4),(11,12),(0.25,0.3),(21,22),(0.15,0.20),(52,54),(0.15,0.20),(61,63),(0.15,0.20),(90,99),(0.3,0.50),(103,107),(0.15,0.20),(125,128),(0.15,0.25),(133,136),(0.15,0.25),(160,165),(0.25,0.40),(165,170),(0.18,0.25),(180,220),(0.19,0.25)]")
JuliaCall::julia_eval("opt = bbsetup(obj;SearchRange = bound1, MaxSteps = 11e2, NumDimensions = 23,
    TraceMode = :compact,
    Method = :adaptive_de_rand_1_bin_radiuslimited)")


for (i in 1:10)
{
  JuliaCall::julia_eval("res1 = bboptimize(opt)")
  p2 <- JuliaCall::julia_eval("p = best_candidate(res1)")
  source('./Modell-Situation.R')
}


p2 <- JuliaCall::julia_eval("p = best_candidate(res1)")

#p2 <- c(0.30971, 11.9905, 0.3, 21.4314, 0.160225, 52.1844, 0.2, 61.2596, 0.182962, 98.9658, 0.346067, 103.072, 0.15, 125.009, 0.249945, 135.757, 0.171135, 162.013, 0.4, 167.06, 0.188391, 186.523, 0.198328)
p2
print(paste(p2[1], 0, dates[days== t[1]+t_0]))
print(paste(p2[3], p2[2], dates[days==t[p2[2]]+t_0]))
print(paste(p2[5], p2[4], dates[days==t[p2[4]]+t_0]))
print(paste(p2[7], p2[6], dates[days==t[p2[6]]+t_0]))
print(paste(p2[9], p2[8], dates[days==t[p2[8]]+t_0]))
print(paste(p2[11], p2[10], dates[days==t[p2[10]]+t_0]))
print(paste(p2[13], p2[12], dates[days==t[p2[12]]+t_0]))
print(paste(p2[15], p2[14], dates[days==t[p2[14]]+t_0]))
print(paste(p2[17], p2[16], dates[days==t[p2[16]]+t_0]))
print(paste(p2[19], p2[18], dates[days==t[p2[18]]+t_0]))
print(paste(p2[21], p2[20], dates[days==t[p2[20]]+t_0]))
print(paste(p2[23], p2[22], dates[days==t[p2[22]]+t_0]))


JuliaCall::julia_assign("saveat", c(0:1000/1000*400))
JuliaCall::julia_assign("tspan", c(0,400))

