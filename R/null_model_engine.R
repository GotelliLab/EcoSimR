#'Run null model
#'@description This drives the null models for all the different kinds of null models that can be run. It is the underlying engine.
#'@param speciesData a dataframe <put some guidelines in here>
#'@param algo the algorithm to use, must be "RA1", "RA2", "RA3", "RA4"
#'@param metric the metric used to caluclate the null model: choices are "Pianka", "Czekanowski", "Pianka.var", "Czekanowski.var", "Pianka.skew", "Czekanowski.skew"; default is Pianka
#'@param nReps the number of replicates to run the null model.
#'@param rowNames Does your dataframe have row names? If yes, they are stripped, otherwise FALSE for data that has no row names
#'@param random.seed Choose a seed to start your random number.  0 will choose a random seed, otherwise set the seed with any integer.
#'
#'@export


null_model_engine <- function(speciesData, algo, metric, nReps = 1000, rowNames = TRUE, saveSeed = FALSE, algoOpts = list(), metricOpts = list())
{
  pb <- txtProgressBar(min = 0, max = nReps, style = 3)
  ## Set the seed
  if(saveSeed){
    randomSeed <- .Random.seed
  } else {
    randomSeed <- NULL
  }
  
  ### Strip out row names
  
  if(rowNames){
    speciesData <- speciesData[,-1] 
  }
  ## Convert to matrix for type consistency
  if(!is.matrix(speciesData)){ speciesData <- as.matrix(speciesData)}
  
  ### check if rownames accidentally was set to FALSE
  if(is.character(speciesData)){stop("Did you forget to set rowNames to TRUE?  Your data is non-numeric")}
  
  
  algoF <- get(algo)
  metricF <- get(metric)
  
  
  startTime <- Sys.time()
  
  
  if (nReps < 2) nReps <- 2
  
  sim <- rep(NA,nReps)

  algoOpts[["speciesData"]] <- speciesData
  metricOpts[["m"]] <- speciesData
  obs <- do.call(metricF,metricOpts)
  
  
  for(i in 1:nReps){
    m <- do.call(algoF,algoOpts)
    metricOpts[["m"]] <- m
    sim[i] <- do.call(metricF,metricOpts)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  ## Save final matrix for plotting
  finalRandomData <- m
  
  endTime <- Sys.time()
  elapsedTime <- format(endTime-startTime,digits=2)
  timeStamp <- date()
  
  ### Reverse engineers the naming for consistent output
  ### Be sure to update the code below if new algos and metrics are added
  

  
  nullModelOut <- list(Obs=obs,Sim=sim, Elapsed.Time=elapsedTime, Time.Stamp=timeStamp,Metric = metric, Algorithm = algo, n.reps = nReps, 
                       Reproducible = saveSeed,RandomSeed = randomSeed, Data = speciesData,Randomized.Data = finalRandomData)
  class(nullModelOut) <- "nullmod"
  return(nullModelOut)
  
}

