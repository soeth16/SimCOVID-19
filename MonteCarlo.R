#####################################################################################
#
# SimCOVID-19
# 
# Simulation of the SARS-CoV-2 pandemic virus outbreak which cause the corona virus 
# disease 2019 (COVID-19) by using a modified S. E. I. R. D. model. This repository 
# is currently in a preview state! It should be used with caution! Currently only 
# german data has been processed, other countries will follow. Further questions can
# answered by soeren.thiering@hs-anhalt.de
#
#####################################################################################

#######################################
# Setup
#######################################

# Working Directory
setwd("~/SimCOVID-19/")

# Output
# graphical: 0
# pdf: 1
# png: 2
plot_outs <- c(2:0)

#######################################
# Libraries
#######################################
library(JuliaCall)
#JuliaCall::julia_eval("import Pkg; Pkg.add(\"DifferentialEquations\")")
#JuliaCall::julia_eval("Pkg.add(\"DiffEqParamEstim\")")
#JuliaCall::julia_eval("Pkg.add(\"Optim\")")
#JuliaCall::julia_eval("Pkg.add(\"Distributions\")")
#JuliaCall::julia_eval("Pkg.add(\"BlackBoxOptim\")")
#JuliaCall::julia_eval("Pkg.add(\"DiffEqUncertainty\")")
#JuliaCall::julia_eval("Pkg.add(\"Distributed\")")
#JuliaCall::julia_eval("using Distributed")
#JuliaCall::julia_eval("addprocs(8 - nprocs(), topology=:all_to_all)")
JuliaCall::julia_eval("using DifferentialEquations")
JuliaCall::julia_eval("using DiffEqParamEstim")
JuliaCall::julia_eval("using Optim")
JuliaCall::julia_eval("using Distributions")
JuliaCall::julia_eval("using BlackBoxOptim")
JuliaCall::julia_eval("using DiffEqUncertainty")

#######################################
# Data & Constants & Sitation Analysis 
#######################################
tcd <- c(-5.710, -5.685, -5.660, -6.360, -5.830,  0.000, -5.890)
for (plot_out in plot_outs) source('./Data-RKI.R')

#######################################
# Modell
#######################################
source('./Modell-DDE.R')

#######################################
# Modell Fit
#######################################
source('./Fit-DDE.R')
p_monte <- array(NA, dim=c(6,length(p2)) )
p_monte[1,] <- p2
for (i_fit in 2:6)
{
  source('./Fit-DDE.R')
  p_monte[i_fit,] <- p2
}
print(p_monte)

