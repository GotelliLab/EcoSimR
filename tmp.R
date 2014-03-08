library(devtools)
library(roxygen2)
roxygenize(".")
dev_mode()
install(".")
library(EcoSimR)


test <- null_model_engine(macwarb)
summary(test)
plot(test)
