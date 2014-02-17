## EcoSimR: R code for Null Model Analysis
##
##
##
##############################################
## EcoSimR Niche Overlap Shell
## Nicholas J. Gotelli & Aaron M. Ellison
##
##
##############################################
## Version 1.0
## 12 May 2013
#############################################
# 17 May default parameters function added by NJG
# 19 May Source files reorganized and simplified by NJG
# 25 May Graphic designation added to param list by NJG

#############################################
##
## 
source("EcoSimR - Algorithms Source.R")
source("EcoSimR - General Functions Source.R")
source("EcoSimR - Metrics Source.R")
source("EcoSimR - Graphics Source.R")
#############################################
Data.File <- "Macarthur Warblers.csv"
Output.File <-"Niche Overlap Output.txt"
Algorithm <- "RA3"  #choices are "RA1", "RA2", "RA3", "RA4"; default is "RA3"
Metric <- "Pianka"  #choices are "Pianka", "Czekanowski", 
#"Pianka.var", "Czekanowski.var",
# "Pianka.skew", "Czekanowski.skew"; default is Pianka
N.Reps <- 1000  # 1000 is the typical number of replicates, but any number > 2 will run
Random.Seed <- 0 ## If 0, uses random integer. User can replace 0 with your integer of choice e.g. 107
Plot.Output <- "screen" 	#choices are "file", "screen", "none"; default is "screen"
Print.Output <- "screen"	#choices are "file", "screen", "none"; default is "screen"
Display.About <- "none" # choices are "screen", "none"; default is "none"
Graphic <- "Niche.Overlap.Plot" # other choices will be added with other modules


Param.List <- Get.Params(Data.File,Output.File,Algorithm,Metric,
                         N.Reps,Random.Seed,Plot.Output,Print.Output,Display.About)
RandomInteger <- Set.The.Seed(Param.List)

