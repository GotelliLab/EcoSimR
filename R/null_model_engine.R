#'Run null model
#'@description This drives the null models for all the different kinds of null models that can be run. It is the underlying engine.
#'@param species_data a dataframe <put some guidelines in here>
#'@param algo the algorithm to use, must be "RA1", "RA2", "RA3", "RA4"
#'@param metric the metric used to caluclate the null model: choices are "Pianka", "Czekanowski", "Pianka.var", "Czekanowski.var", "Pianka.skew", "Czekanowski.skew"; default is Pianka
#'@param n.reps the number of replicates to run the null model.
#'@param row.names Does your dataframe have row names? If yes, they are stripped, otherwise FALSE for data that has no row names
#'@param random.seed Choose a seed to start your random number.  0 will choose a random seed, otherwise set the seed with any integer.
#'@examples \dontrun{
#' ## Load MacAruthur warbler data
#' data(macwarb)
#' 
#' ## Run the null model
#' warbMod <- null_model_engine(macwarb)
#' ## Summary and plot info
#' summary(warbMod)
#' plot(warbMod)
#'  
#'}
#'
#'@export


null_model_engine <- function(species_data, algo, metric, n.reps = 1000, row.names = TRUE, random.seed = 0)
{
  ## Set the seed
  ifelse (random.seed==0, RandomInteger <- trunc(runif(1,-2000000000,2000000000)), RandomInteger <- random.seed)
  
  set.seed(RandomInteger)
  
  ### Strip out row names
  
  if(row.names){
    species_data <- species_data[,-1] 
  }
  ## Convert to matrix for type consistency
  if(!is.matrix(species_data)){ species_data <- as.matrix(species_data)}
  
  ### check if rownames accidentally was set to FALSE
  if(is.character(species_data)){stop("Did you forget to set row.names to TRUE?  Your data is non-numeric")}
  
  
  algoF <- eval(parse(text = algo))
  metricF <- eval(parse(text = metric))
  
  
  Start.Time <- Sys.time()
  
  
  if (n.reps < 2) n.reps <- 2
  Sim <- replicate(n.reps,metricF(algoF(species_data)))
  Obs <- metricF(species_data)
  
  End.Time <- Sys.time()
  Elapsed.Time <- format(End.Time-Start.Time,digits=2)
  Time.Stamp <- date()
  
  Null.Model.Out <- list(Obs=Obs,Sim=Sim, Elapsed.Time=Elapsed.Time, Time.Stamp=Time.Stamp,Metric = metric, Algorithm = algo, N.Reps = n.reps, RandomInteger = RandomInteger, Data = species_data)
  class(Null.Model.Out) <- "nullmod"
  return(Null.Model.Out)
  
}

