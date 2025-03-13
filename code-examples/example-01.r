################################################################################
# File: alcohol_cr01.R
# Purpose: prep the alcohol dataset

# Notes: current plan:
#	* load samplesheet, chop into timepoints, load covariates in sourced scripts,
#	* merge, estimate
################################################################################
rm(list=ls())
gc()

packages <- c("reshape2", "gdata", "plyr")
lapply(packages, require, character.only=T)

################################################################################
# SET WD
setwd("~/alcohol_jan2017")	

load("mydata.rda")

## load new collaborator data
load("/home/user/data/WS_LIX/LIX_analyses/PdP_MC/new_Results/ReCent/FinalFINAl")

## 
wt <- read.csv("/home/user/consoritiumData/Alcohol analysis/2019updated-weights")





