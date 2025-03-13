args <- commandArgs(trailingOnly=TRUE)

# config paths
config.name <- "default"
if (length(args) > 0)
    config.name <- args[1]

source("config.r") #paths
paths <- paths[[config.name]]
print(paths)




