EcoSimR
=======

Test repository for EcoSimR, by Gotelli, N.J. and A.M. Ellison. 2013. EcoSimR 1.00.  http://www.uvm.edu/~ngotelli/EcoSim/EcoSim.html

QuickStart
=======

First install the dev branch
```coffescript
library(devtools)
## use dev_mode() if you want to just play with this in a sandbox
install_github("gotellilab/EcoSimR@dev")
```

Here's a quick example with the MacAurthur Warbler data.
 
```coffeescript
library(EcoSimR)

test <- null_model_engine(macwarb)
summary(test)
plot(test)
