library(devtools)
dev_mode()
install(".")
library(EcoSimR)


test <- null_model_engine(macwarb)
summary(test)
plot(test)
