setwd('C:/Users/Timothée/Documents/GitHub/CarpentryExercise/')
#source('C:/Users/Timothée/Documents/GitHub/CarpentryExercise/S4imulations.R')
#source('C:/Users/Timothée/Documents/GitHub/CarpentryExercise/Application.R')
#package.skeleton(name = "S4imul",code_files = "S4imulations.R",force = TRUE)
#library("S4imul")

library(devtools)
library(Rcpp)
Rcpp.package.skeleton(name =  "S4imulC",code_files = "S4imulations.R",cpp_files = "sizeSurvivalC.cpp",force = TRUE)

load_all("S4imulC/")

library("S4imulC")
system.time(replicate(n = 100000, expr = bathtubC(0.2))) # twice as fast only... not very useful!
system.time(replicate(n = 100000, expr = bathtub(0.2)))

