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

# Base Directory
setwd("~/SimCOVID-19")
if(!dir.exists("../SimCOVID-19_MonteCarlo")) dir.create("../SimCOVID-19_MonteCarlo", showWarnings = TRUE, recursive = FALSE, mode = "0775")
if(!file.exists("../SimCOVID-19_MonteCarlo/Data-RKI.R")) file.copy('Data-RKI.R', "../SimCOVID-19_MonteCarlo/Data-RKI.R")
if(!file.exists("../SimCOVID-19_MonteCarlo/Modell-DDE.R")) file.copy('Modell-DDE.R', "../SimCOVID-19_MonteCarlo/Modell-DDE.R")
if(!file.exists("../SimCOVID-19_MonteCarlo/Fit-DDE.R")) file.copy('Fit-DDE.R', "../SimCOVID-19_MonteCarlo/Fit-DDE.R")
setwd("../SimCOVID-19_MonteCarlo/")
dir()

# Working Directory
now_date <- format(Sys.time(), "%m-%d-%Y_%H:%M:%OS3")
dir.create(now_date, showWarnings = TRUE, recursive = FALSE, mode = "0775")
setwd(now_date)

# Output
# graphical: 0
# pdf: 1
# png: 2
plot_out <<- 1

#######################################
# Libraries
#######################################
library(JuliaCall)
#JuliaCall::julia_eval("import Pkg; Pkg.add(\"DifferentialEquations\"); nothing")
#JuliaCall::julia_eval("Pkg.add(\"DiffEqCallbacks\"); nothing")
#JuliaCall::julia_eval("Pkg.add(\"DiffEqParamEstim\"); nothing")
#JuliaCall::julia_eval("Pkg.add(\"NLSolversBase\"); nothing")
#JuliaCall::julia_eval("Pkg.add(\"Optim\"); nothing")
#JuliaCall::julia_eval("Pkg.add(\"Distributions\"); nothing")
#JuliaCall::julia_eval("Pkg.add(\"BlackBoxOptim\"); nothing")
#JuliaCall::julia_eval("Pkg.add(\"DiffEqUncertainty\"); nothing")
#JuliaCall::julia_eval("Pkg.add(\"Distributed\"); nothing")
#JuliaCall::julia_eval("using Distributed; nothing")
#JuliaCall::julia_eval("addprocs(8 - nprocs(), topology=:all_to_all); nothing")
JuliaCall::julia_eval("using NLSolversBase; nothing")
JuliaCall::julia_eval("using DifferentialEquations; nothing")
JuliaCall::julia_eval("using DiffEqCallbacks; nothing")
JuliaCall::julia_eval("using DiffEqParamEstim; nothing")
JuliaCall::julia_eval("using Optim; nothing")
JuliaCall::julia_eval("using Distributions; nothing")
JuliaCall::julia_eval("using BlackBoxOptim; nothing")
JuliaCall::julia_eval("using DiffEqUncertainty; nothing")

#######################################
# Data & Constants & Sitation Analysis 
#######################################
tcd <<- c(-5.710, -5.685, -5.660, -6.360, -5.830,  0.000, -5.890)
source('../Data-RKI.R')

#######################################
# Modell
#######################################
source('../Modell-DDE.R')

#######################################
# Modell Fit
#######################################

p2 <<- 0
source('../Fit-DDE.R')
p_monte <<- array(NA, dim=c(6,length(p2)) )
p_monte[1,] <<- p2
for (i_fit in 2:6)
{
  setwd("..")
  now_date <- format(Sys.time(), "%m-%d-%Y_%H:%M:%OS3")
  dir.create(now_date, showWarnings = TRUE, recursive = FALSE, mode = "0775")
  setwd(now_date)
  source('../Fit-DDE.R')
  p_monte[i_fit,] <- p2
}

setwd("..")
write_delim(p_monte, "result.csv", append=T, delim="\t")
print(p_monte)
