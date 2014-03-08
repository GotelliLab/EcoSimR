#'Run null model
#'@description Create a null model object <put a better description in here>
#'@param species_data a dataframe <put some guidelines in here>
#'@param algo the algorithm to use, must be "RA1", "RA2", "RA3", "RA4"
#'@param metric the metric used to caluclate the null model: choices are "Pianka", "Czekanowski", "Pianka.var", "Czekanowski.var", "Pianka.skew", "Czekanowski.skew"; default is Pianka
#'@param N.reps the number of replicates to run the null model for
#'@examples \dontrun{
#' Put an example in here
#'  
#'}
#'
#'@export

null_model_engine <- function(species_data, algo = "RA3", metric = "Pianka", N.Reps = 1000)
{


  algoF <- eval(parse(text = match.arg(algo,choices = c("RA1","RA2","RA3","RA4"))))
  metricF <- eval(parse(text = match.arg(metric,choices = c("Pianka", "Czekanowski", "Pianka.var", "Czekanowski.var", "Pianka.skew", "Czekanowski.skew"))))
  
  
  Start.Time <- Sys.time()
  
  
  if (N.Reps < 2) N.Reps <- 2
  Sim <- replicate(N.Reps,metricF(algoF(species_data)))
  Obs <- metricF(species_data)
  
  End.Time <- Sys.time()
  Elapsed.Time <- format(End.Time-Start.Time,digits=2)
  Time.Stamp <- date()
  
  Null.Model.Out <- list(Obs=Obs,Sim=Sim, Elapsed.Time=Elapsed.Time, Time.Stamp=Time.Stamp,Metric = metric, Algorithm = algo)
  class(Null.Model.Out) <- "nichenullmod"
  return(Null.Model.Out)
  
}

#############################################
## Null.Model.Summary function
## Generic function for calculating null model summary statistics
## Takes as input a list of Null.Model.Out, with Obs, Sim, Elapsed Time, and Time Stamp values
##
## Rewritten 24 January 2013 by AME to allow for nice screen or file output
## Also eliminates need for Print.Results function
##

summary.nichenullmod <- function(nullmodObj)
{ 
  

  
  
  #if (!is.null(Output.File)) outfile <- file(p$Output.File, "w") else outfile <-""
  
  cat("Time Stamp: " , nullmodObj$Time.Stamp,   "\n") 
 # cat("Data File: ", p$Data.File,  "\n")
#  cat("Output File: ", p$Output.File,  "\n") 
  cat("Random Number Seed: ",RandomInteger,  "\n")
  cat("Number of Replications: ",nullmodObj$N.Reps,  "\n")
  cat("Elapsed Time: ", nullmodObj$Elapsed.Time, "\n")
  cat("Metric: ", nullmodObj$Metric,  "\n")
  cat("Algorithm: ", nullmodObj$Algorithm,  "\n") 
  
  cat("Observed Index: ", format(nullmodObj$Obs,digits=5),  "\n")
  cat("Mean Of Simulated Index: ",format(mean(nullmodObj$Sim),digits=5),  "\n")
  cat("Variance Of Simulated Index: ",format(var(Sim),digits=5),  "\n")
  cat("Lower 95% (1-tail): ",format(quantile(nullmodObj$Sim,0.05),digits=5),  "\n")
  cat("Upper 95% (1-tail): ",format(quantile(nullmodObj$Sim,0.95),digits=5), "\n")
  cat("Lower 95% (2-tail): ",format(quantile(nullmodObj$Sim,0.025),digits=5), "\n")
  cat("Upper 95% (2-tail): ",format(quantile(nullmodObj$Sim,0.975),digits=5),  "\n")
  
  #P-values
  if (Obs > max(Sim)) {
    cat("P(Obs <= null) < ",1/length(nullmodObj$Sim),  "\n")
    cat("P(Obs >= null) > ",(length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim),  "\n")
  } else if (sum(nullmodObj$Obs <= nullmodObj$Sim)==0) {
    cat("P(Obs <= null) > ",(length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim), "\n")
    cat("P(Obs >= null) < ", 1/length(nullmodObj$Sim), "\n")
  } else {
    cat("P(Obs <= null) = ", format(sum(nullmodObj$Obs >= nullmodObj$Sim)/length(nullmodObj$Sim),digits=5),  "\n")
    cat("P(Obs >= null) = ", format(sum(nullmodObj$Obs <= nullmodObj$Sim)/length(nullmodObj$Sim),digits=5), "\n")
  }
  
  cat("P(Obs = null) = ",format(sum(nullmodObj$Obs == nullmodObj$Sim)/length(nullmodObj$Sim),digits=5),  "\n")
  cat("Standardized Effect Size (SES): ", format((nullmodObj$Obs - mean(nullmodObj$Sim))/sd(nullmodObj$Sim),digits=5), "\n")
  
  #if(!is.null(Output.File)) close(outfile)
}
