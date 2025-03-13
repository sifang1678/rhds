################################################################################
packages <- c("tidyverse", "ewaff", "R.utils")
lapply(packages, require, character.only=T)

# dir <- "/project/alcohol-use/working" #use 4 Isambard
dir <- "/tiedye/pebble/alcohol-use/working" #4 pebble

dir.data <- file.path(dir, "data")
dir.figs <- file.path(dir, "figs")
dir.logs <- file.path(dir, "logs")



