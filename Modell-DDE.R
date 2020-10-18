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

################################
# Model Subfunctions
################################

JuliaCall::julia_eval("@inbounds function k(p, t)::Float64
    if (t < p[2]) 
      p[1]
    elseif (t >= p[2] && t < p[4]) 
      p[3]
    elseif (t >= p[4] && t < p[6]) 
      p[5]
    elseif (t >= p[6] && t < p[8]) 
      p[7]
    elseif (t >= p[8] && t < p[10]) 
      p[9]
    elseif (t >= p[10] && t < p[12]) 
      p[11]
    elseif (t >= p[12] && t < p[14]) 
      p[13]
    elseif (t >= p[14] && t < p[16]) 
      p[15]
    elseif (t >= p[16] && t < p[18]) 
      p[17]
    elseif (t >= p[18] && t < p[20]) 
      p[19]
    elseif (t >= p[20] && t < p[22]) 
      p[21]
    elseif (t >= p[22]) 
      p[23]
    end
  end")

k <- function(p, t) {
  tmp <- 0
  for (tt in 1:length(t))
  {
    if (t[tt] < p[2]) tmp[tt] <- p[1]
    else if (t[tt] >= p[2] && t[tt] < p[4])  tmp[tt] <- p[3]
    else if (t[tt] >= p[4] && t[tt] < p[6])  tmp[tt] <- p[5]
    else if (t[tt] >= p[6] && t[tt] < p[8])  tmp[tt] <- p[7]
    else if (t[tt] >= p[8] && t[tt] < p[10]) tmp[tt] <- p[9]
    else if (t[tt] >= p[10] && t[tt] < p[12]) tmp[tt] <- p[11]
    else if (t[tt] >= p[12] && t[tt] < p[14]) tmp[tt] <- p[13]
    else if (t[tt] >= p[14] && t[tt] < p[16]) tmp[tt] <- p[15]
    else if (t[tt] >= p[16] && t[tt] < p[18]) tmp[tt] <- p[17]
    else if (t[tt] >= p[18] && t[tt] < p[20]) tmp[tt] <- p[19]
    else if (t[tt] >= p[20] && t[tt] < p[22]) tmp[tt] <- p[21]
    else if (t[tt] >= p[22]) tmp[tt] <- p[23]
  }
  return (tmp)
}


# introduction of dextamethasone
JuliaCall::julia_eval(paste0("@fastmath phi_dextamethasone(t)::Float64 = (t > ", tm, " ? ", phi_m, " : 1)"))

JuliaCall::julia_eval(paste0("@fastmath phi_d_drift(t)::Float64 = (t > ", max(days) - t_0, " ? 1 : ", max(days) - t_0, "- t + 1) * ", phi_dd))

# phi_d for distributed hospital usage
JuliaCall::julia_eval(paste0("@fastmath phi_d(u, t)::Float64 = ((u > ",n_h_max," ?  (",phi_d_1 * n_h_max,"/ u + ",phi_d_2," * (u - ",n_h_max,") / u) : ",phi_d_1," ) * phi_d_drift(t)) * phi_dextamethasone(t)"))

# phi_d for hot spots without distribution
# JuliaCall::julia_eval(paste0("@fastmath phi_d(u, t)::Float64 = (u / ",n_h_max," < 1 ?  (",phi_d_1," * (1 - (u / ",n_h_max,")) + ",phi_d_2," * u / ",n_h_max,") : ",phi_d_2," ) * phi_dextamethasone(t)"))


################################
# Model Function
################################

f = JuliaCall::julia_eval(paste0("@fastmath function f(du, u, h, p, t)
    
    
    # Susceptibles
    du[7] = (
      - k(p, t) * u[7]::Float64 / ", n_max, " * u[6]::Float64
    )
    
    # Exposed / Incubating
    du[6] = (
      - du[7]
      +h(p, t ", - te, ", Val{1}; idxs = 7)::Float64
    )
    
    # Infected
    du[5] = (
      -h(p, t ", - te, ", Val{1}; idxs = 7)::Float64
      +h(p, t ", - te - ti, ", Val{1}; idxs = 7)::Float64* (1-",phi_h,")
      +h(p, t ", - te - th, ", Val{1}; idxs = 7)::Float64 * ",phi_h,"
      -h(p, t ", - te - th - thi, ", Val{1}; idxs = 7)::Float64 * (",phi_h," - phi_d(h(p, t ", - thi + thd, "; idxs = 4)::Float64, t ", - thi + thd, "))
      +h(p, t ", - te - th - thi - thii, ", Val{1}; idxs = 7)::Float64 * (",phi_h," - phi_d(h(p, t ", - thi + thd - thii, "; idxs = 4)::Float64, t ", - thi + thd - thii, "))
    )
    # Hospitalization 
    du[4] = (
      - h(p, t ", - te - th, ", Val{1}; idxs = 7)::Float64 *  ",phi_h," 
      + h(p, t ", - te - th - thi, ", Val{1}; idxs = 7)::Float64 * (",phi_h," - phi_d(h(p, t ", - thi + thd, "; idxs = 4)::Float64, t ", - thi + thd, "))
      + h(p, t ", - te - th - thd, ", Val{1}; idxs = 7)::Float64 * phi_d(u[4]::Float64, t)
    )
    
    # Recovered
    du[3] = ( 
      - h(p, t ", - te - ti, ", Val{1}; idxs = 7)::Float64 *  (1-",phi_h,") 
      - h(p, t ", - te - th - thi - thii, ", Val{1}; idxs = 7)::Float64 * (",phi_h," - phi_d(h(p, t ", - thi + thd - thii, "; idxs = 4)::Float64, t ", - thi + thd - thii, "))
    )
    
    # Deaths
    du[2] = ( 
      - h(p, t ", - te - th - thd, ", Val{1}; idxs = 7)::Float64 *  phi_d(u[4]::Float64, t)
    )
    
    # Confirmed
    du[1] = du[5] + du[4] + du[3] + du[2]
    
  end"))

JuliaCall::julia_eval(paste0("lags = [",te,", ",te + ti,", ",te + th,", ",te + th + thi,", ",te + th + thd,", ", thi - thd,", ",te + th + thi + thii,"]"))


################################
# Past
################################
# dx/dt = k * x
# 1/x * dx = k * dt
# integrate
# ln(x) - ln(x_0) = k * t
# x / x_0 = exp(k * t)
# x = x_0 exp(k * t)
# x_0 = 1

JuliaCall::julia_eval(paste0("h(p, t; idxs::Union{Nothing,Int} = nothing)::Float64 =  t > ", - te, " ? u0[idxs] * exp(p[1] * t) : 0"))

JuliaCall::julia_eval(paste0("h(p, t, deriv::Type{Val{1}}; idxs::Union{Nothing,Int} = nothing)::Float64 = t > ", - te, " ? (- k(p, t) * u0[7] / ", n_max, " * u0[6] * (exp(k(p, t ", - te, ") * t))) : 0"))


################################
# Parameter
################################
#t_0 <- 60
k_raw = 0.3626059 
k1 = 0.412884

t <- cumulative_extrapolated[[8]]$Time  >= t_0
C <- cumulative_extrapolated[[8]]$ConfirmedCases[t]
R <- cumulative_extrapolated[[8]]$RecoverdCases[t]
D <- cumulative_extrapolated[[8]]$DeathCases[t]
I <- C - R - D
E <- I * k_raw/(k1 * exp(-k1*te))
I <- I * 0.95
H <- I * 0.05
S <- n_max - I - E - H - R - D
t <- cumulative_extrapolated[[8]]$Time[t]-t_0
t
data_df <- data.frame(C, D, R)
data_sd <- abs(data.frame(C = lowess(t,C,f=0.1)$y - data_df$C, D = lowess(t,D,f=0.3)$y - data_df$D, R = lowess(t,R,f=0.3)$y - data_df$R))
data_df

p1 <- c(0.3488631, 11.3885073, 0.2633408, 21.5223237, 0.1649659, 52.5617909, 0.1733183, 61.8837740, 0.1553055, 95.7754461, 0.3427507, 103.2323733, 0.1701468, 126.6378092, 0.1844588, 133.9112688, 0.2032979, 164.5808251, 0.3328118, 166.3242064, 0.1967201, 202.4674021, 0.2190007)

JuliaCall::julia_assign("u0", c(C[1], D[1], R[1], H[1], I[1], E[1], S[1]))
JuliaCall::julia_assign("saveat", c(0:1000/1000*400))
JuliaCall::julia_assign("tspan", c(0,400))
JuliaCall::julia_assign("p", p1)