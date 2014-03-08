library(devtools)
library(roxygen2)
roxygenize(".")
dev_mode()
install(".")
library(EcoSimR)

extdat <- read.csv("~/Documents/macwarb.csv",header=T,stringsAsFactors=F)

test <- null_model_engine(macwarb)
