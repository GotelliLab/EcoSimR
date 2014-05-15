
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


##


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

