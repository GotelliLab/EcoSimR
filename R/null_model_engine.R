#'Run null model
#'@description Create a null model object <put a better description in here>
#'@param species_data a dataframe <put some guidelines in here>
#'@param output name of file to write results to
#'@param algo the algorithm to use, must be "RA1", "RA2", "RA3", "RA4"
#'@param metric the metric used to caluclate the null model: choices are "Pianka", "Czekanowski", "Pianka.var", "Czekanowski.var", "Pianka.skew", "Czekanowski.skew"; default is Pianka
#'@param N.reps the number of replicates to run the null model for
#'@examples \dontrun{
#' Put an example in here
#'  
#'}
#'
#'@export

null_model_engine <- function(species_data, output, algo = "RA3", metric = "Pianka", N.Reps = 1000)
{


  
  
  Start.Time <- Sys.time()
  
  
  if (N.Reps < 2) N.Reps <- 2
  Sim <- replicate(N.Reps,Metric(algo(species_data)))
  Obs <- metric(Data.Matrix)
  
  End.Time <- Sys.time()
  Elapsed.Time <- format(End.Time-Start.Time,digits=2)
  Time.Stamp <- date()
  
  Null.Model.Out <- list(Obs=Obs,Sim=Sim, Elapsed.Time=Elapsed.Time, Time.Stamp=Time.Stamp)
  
  return(Null.Model.Out)
  
}