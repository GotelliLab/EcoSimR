##############################################
## Ant community analysis from warming experiment at Cedar Creek LTER, MN
## Adam Clark, 23 Dec 2013
## University of Minnesota, Tilman Lab
##############################################

rm(list=ls())   # remove all objects in memory


## load EcoSimR
setwd("../../base_functions/")
source("EcoSimR - Main Source.R")
setwd("../contributed/adam_ants/")
source("Dist_test.R")

## Create data matrix
bacdat<-read.csv("bacants.csv")

## Load location data
locdat<-read.csv("bacplot_locations.csv")

## Calculate distances between all plots
posmatrix<-cbind(locdat$Lat, locdat$Long)
distances<-as.matrix(dist(posmatrix))

## For preliminary analysis, we use only the 2011 data
bacdat2011<-subset(bacdat, Year==2011)

antdat<-table(bacdat2011$spCode, bacdat2011$PlotID)
write.csv(antdat, "antdat.csv")

## Load vegan package
require(vegan)

Data.File <- "antdat.csv"
Output.File <-"2011_ant_distance"
Algorithm <- "RA3"	#choices are "RA1", "RA2", "RA3", "RA4"; default is "RA3"
Metric <- "Dist_test"	#
N.Reps <- 1000	# 1000 is the typical number of replicates, but any number > 2 will run
Random.Seed <- 625 ## If 0, uses random integer. User can replace 0 with your integer of choice e.g. 107
Plot.Output <- "file" 	#choices are "file", "screen", "none"; default is "screen"
Print.Output <- "file"	#choices are "file", "screen", "none"; default is "screen"
Display.About <- "none" # choices are "screen", "none"; default is "none"
Graphic <- "Niche.Overlap.Plot" # other choices will be added with other modules
#############################################


##############################################

Param.List <- Get.Params(Data.File,Output.File,Algorithm,Metric,
                         N.Reps,Random.Seed,Plot.Output,Print.Output,Display.About,Graphic)
RandomInteger <- Set.The.Seed(Param.List)


Null.Result <- Null.Model.Engine(Param.List)

Output.Results(Param.List,Null.Result)