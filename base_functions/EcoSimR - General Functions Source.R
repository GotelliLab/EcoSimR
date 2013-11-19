##############################################
## EcoSimR - R Code For Null Model Analysis
##
##
##  
##############################################
## EcoSimR General Functions Source File
## Nicholas J. Gotelli & Aaron M. Ellison
##
##
##############################################
## Version 1.00
## 12 May 2013
#############################################
##################################
# Get.Params function
# Bundles null models parameters into a single list
Get.Params <- function (P1="Macarthur Warblers.csv",
                        P2="Niche Overlap Output.txt",
                        P3="RA3",
                        P4="Pianka",
                        P5=1000,
                        P6=0,
                        P7="screen",
                        P8="screen",
                        P9="none",
                        P10="Niche.Overlap.Plot")
{
  Params <- list(Data.File=P1,Output.File=P2,
                 Algorithm=P3,Metric=P4,
                 N.Reps=P5,Random.Seed=P6,
                 Plot.Output=P7,Print.Output=P8,Display.About=P9,
                 Graphic=P10)
  return(Params)
}
####################################
#############################################
## DataMatrix function

## Update May 17 2013 AME to error trap on NA

## Takes a .csv declared as "Data.File" and turns it into a matrix, retains header rows and row names

Data.Read <- function(Data.File="Macarthur Warblers.csv")
{
Data.Matrix <- as.matrix(read.csv(file=Data.File,header=TRUE,row.names=1))

	ifelse(sum(is.na(Data.Matrix)) >0, stop("Input Matrix has NA values", call. = FALSE), return(Data.Matrix))
}
##
#############################################

#####################################
## Random seed generator function
# The default "System" uses a random integer and stores it
# Otherwise, the user supplies a random integer that is stored


Set.The.Seed <- function(p=Param.List) {
  
  ifelse (p$Random.Seed==0, RandomInteger <- trunc(runif(1,-2000000000,2000000000)), RandomInteger <- p$Random.Seed)

  set.seed(RandomInteger)
  
  return(RandomInteger)
} 
##
#############################################
#############################################
## Null.Model.Engine function
## Input is a data matrix, number of replicates, an algorithm, and a metric
## Output is a list with the Observed metric, a vector of simulated metrics, Elapsed Time, and Time Stamp
##
##


Null.Model.Engine <- function(p=Param.List)
{
  Data.Matrix <- Data.Read(p$Data.File)
  N.Reps <- p$N.Reps
  Algorithm <- get(p$Algorithm)
  Metric <- get(p$Metric)
  
  
  Start.Time <- Sys.time()
  
  
  if (N.Reps < 2) N.Reps <- 2
  Sim <- replicate(N.Reps,Metric(Algorithm(Data.Matrix)))
  Obs <- Metric(Data.Matrix)
  
  End.Time <- Sys.time()
  Elapsed.Time <- format(End.Time-Start.Time,digits=2)
  Time.Stamp <- date()
  
  Null.Model.Out <- list(Obs=Obs,Sim=Sim, Elapsed.Time=Elapsed.Time, Time.Stamp=Time.Stamp)
  
  return(Null.Model.Out)
  
}

##

#############################################
## Null.Model.Summary function
## Generic function for calculating null model summary statistics
## Takes as input a list of Null.Model.Out, with Obs, Sim, Elapsed Time, and Time Stamp values
##
## Rewritten 24 January 2013 by AME to allow for nice screen or file output
## Also eliminates need for Print.Results function
##

Null.Model.Summary <- function(p=Param.List,Null.Model.Out=list(Obs=runif(1),Sim=runif(n=1000), 
								Elapsed.Time=NULL, Time.Stamp=NULL), Output.File=NULL)
{ 


Obs <- Null.Model.Out$Obs
Sim <- Null.Model.Out$Sim
Time.Stamp <- Null.Model.Out$Time.Stamp
Elapsed.Time <- Null.Model.Out$Elapsed.Time 


if (!is.null(Output.File)) outfile <- file(p$Output.File, "w") else outfile <-""

cat("Time Stamp: " , Time.Stamp, file=outfile,  "\n") 
cat("Data File: ", p$Data.File, file=outfile, "\n")
cat("Output File: ", p$Output.File, file=outfile, "\n") 
cat("Random Number Seed: ",RandomInteger, file=outfile, "\n")
cat("Number of Replications: ",p$N.Reps, file=outfile, "\n")
cat("Elapsed Time: ", Elapsed.Time, file=outfile, "\n")
cat("Metric: ", p$Metric, file=outfile, "\n")
cat("Algorithm: ", p$Algorithm, file=outfile, "\n") 
  
cat("Observed Index: ", format(Obs,digits=5), file=outfile, "\n")
cat("Mean Of Simulated Index: ",format(mean(Sim),digits=5), file=outfile, "\n")
cat("Variance Of Simulated Index: ",format(var(Sim),digits=5), file=outfile, "\n")
cat("Lower 95% (1-tail): ",format(quantile(Sim,0.05),digits=5), file=outfile, "\n")
cat("Upper 95% (1-tail): ",format(quantile(Sim,0.95),digits=5), file=outfile, "\n")
cat("Lower 95% (2-tail): ",format(quantile(Sim,0.025),digits=5), file=outfile, "\n")
cat("Upper 95% (2-tail): ",format(quantile(Sim,0.975),digits=5), file=outfile, "\n")

#P-values
 if (Obs > max(Sim)) {
	cat("P(Obs <= null) < ",1/length(Sim), file=outfile, "\n")
	cat("P(Obs >= null) > ",(length(Sim) - 1)/length(Sim), file=outfile, "\n")
 } else if (sum(Obs <= Sim)==0) {
	cat("P(Obs <= null) > ",(length(Sim) - 1)/length(Sim), file=outfile, "\n")
	cat("P(Obs >= null) < ", 1/length(Sim), file=outfile, "\n")
 } else {
	cat("P(Obs <= null) = ", format(sum(Obs >= Sim)/length(Sim),digits=5), file=outfile, "\n")
	cat("P(Obs >= null) = ", format(sum(Obs <= Sim)/length(Sim),digits=5), file=outfile, "\n")
 }

cat("P(Obs = null) = ",format(sum(Obs == Sim)/length(Sim),digits=5), file=outfile, "\n")
cat("Standardized Effect Size (SES): ", format((Obs - mean(Sim))/sd(Sim),digits=5), file=outfile, "\n")

if(!is.null(Output.File)) close(outfile)
}
##
#############################################

#############################################
## Output Results function
## Takes Null.Result, Data.File, user choices for Plot and Print outputs, and filename for Output.File
## Includes calls to plotting functions in Graphics Source
## Includes calls to textfile output in Null.Model.Summary

Output.Results <- function(p=Param.List,Null.Result=Null.Model.Engine())

{
  Graphic <-get(p$Graphic)
  if (p$Display.About == "screen")
  {
    cat("###############################################################", "\n")
    cat("\n")
    cat("EcoSimR - R Code For Null Model Analysis", "\n")
    cat("\n")
    cat("Nicholas J. Gotelli & Aaron M. Ellison", "\n")
    cat("\n")
    cat("website: http://www.uvm.edu/~ngotelli/EcoSim/EcoSim.html", "\n")
    cat("\n")
    cat("Version 1.00", "\n")
    cat("15 June 2013", "\n")
    cat("###############################################################", "\n")
    cat("\n")
    cat("Citation Format:", "\n")
    cat("Gotelli, N.J. and A.M. Ellison. 2013. EcoSimR. Version 1.00.", "\n")
    cat("http://www.uvm.edu/~ngotelli/EcoSim/EcoSim.html", "\n")
    cat("\n") 
    cat("###############################################################", "\n")
    cat("\n")
    cat("Contacting the authors:", "\n")
    cat("\n")
    cat("Nicholas J. Gotelli               Aaron M. Ellison", "\n")
    cat("ngotelli@uvm.edu                  aellison@fas.harvard.edu", "\n")
    cat("\n") 
    cat("Department of Biology             Harvard University", "\n")
    cat("University of Vermont             Harvard Forest", "\n")
    cat("Burlington, Vermont               Petersham, Massachusetts", "\n")
    cat("05405 USA                         01366 USA", "\n")
    cat("###############################################################", "\n")
    cat("\n") 
    cat("\n") 
    cat("\n") 
    cat("\n") 
  }
  
if (p$Plot.Output == "screen")
	{
		par(mfrow=c(3,1))
		Null.Model.Plot(Null.Result,Null.Result$Time.Stamp)
	  Graphic(Data.Read(p$Data.File),p$Algorithm,Null.Result$Time.Stamp, p$Plot.Output)

	} else if (p$Plot.Output == "file") 

	{
  		pdf("Null Distribution Plot.pdf") 
			Null.Model.Plot(Null.Result,Null.Result$Time.Stamp)
  		dev.off()

  		pdf("Data Visualization Plot.pdf")
			Graphic(Data.Read(p$Data.File),p$Algorithm,Null.Result$Time.Stamp, p$Plot.Output)
  		dev.off()

	} else {
		print("No graphics requested")
	}


if (p$Print.Output == "file") 
	{
		Null.Model.Summary(p,Null.Result, p$Output.File)
	} else if (p$Print.Output == "screen") 

	{
    
		Null.Model.Summary(p,Null.Result, Output.File=NULL) 
	} else {
		cat("No statistical output requested", "\n")
		cat("No output file created", "\n")
	} 


}

##
#############################################

