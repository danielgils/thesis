
args <- commandArgs(TRUE)

dataT <- as.double(args[1])

#functions and data
source("settings_Rsims_4atts.R")

estim_T4(dataT)