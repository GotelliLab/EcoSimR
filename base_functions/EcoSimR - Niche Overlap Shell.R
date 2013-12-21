##############################################
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
## Version 1.00
## 15 June 2013
#############################################
## Modified on 18 May 2013 by NJG to pass a single Param.List to all functions
#############################################
####For beginners - start here####

## clean the slate ##

rm(list=ls())   # remove all objects in memory


## load EcoSimR
 
source("EcoSimR - Main Source.R")

#############################################
## Model input parameters
## USER CAN MODIFY PARAMETERS IN THIS SECTION
Data.File <- "Macarthur Warblers.csv"
Output.File <-"Niche Overlap Output.txt"
Algorithm <- "RA3"	#choices are "RA1", "RA2", "RA3", "RA4"; default is "RA3"
Metric <- "Pianka"	#choices are "Pianka", "Czekanowski", 
                              #"Pianka.var", "Czekanowski.var",
                              # "Pianka.skew", "Czekanowski.skew"; default is Pianka
N.Reps <- 1000	# 1000 is the typical number of replicates, but any number > 2 will run
Random.Seed <- 625 ## If 0, uses random integer. User can replace 0 with your integer of choice e.g. 107
Plot.Output <- "screen" 	#choices are "file", "screen", "none"; default is "screen"
Print.Output <- "screen"	#choices are "file", "screen", "none"; default is "screen"
Display.About <- "none" # choices are "screen", "none"; default is "none"
Graphic <- "Niche.Overlap.Plot" # other choices will be added with other modules
#############################################


##############################################
## Execute analyses
## Beginning users should NOT modify this section
##
## First command initialized the parameter list from the user inputs
## Second command runs niche overlap analysis using Data.File, Algorithm, and Metric from user inputs
## Third command outputs graphics and statistics to devices specified from user inputs

Param.List <- Get.Params(Data.File,Output.File,Algorithm,Metric,
                         N.Reps,Random.Seed,Plot.Output,Print.Output,Display.About,Graphic)
RandomInteger <- Set.The.Seed(Param.List)


Null.Result <- Null.Model.Engine(Param.List)

Output.Results(Param.List,Null.Result)
