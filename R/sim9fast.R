#'sim9.fast
#'@description A special implementation of the sequential swap algorithm
#'@param speciesData a dataframe in which rows are species, columns are sites,
#' and the entries indicate the absence (0) or presence (1) of a species in a 
#' site. Empty rows and empty columns should not be included in the matrix.
#'@param algo the algorithm to use, must be "sim1", "sim2", "sim3", "sim4", "sim5", "sim6", "sim7", "sim8", "sim9", "sim10"
#'@param metric the metric used to caluclate the null model: choices are "species_combo", "checker", "c_score", "c_score_var", "c_score_skew", "v_ratio"; default is "c_score"
#'@param nReps the number of replicates to run the null model.
#'@param saveSeed TRUE or FALSE.  If TRUE the current seed is saved so the simulation can be repeated
#'@param burn_in The number of burn_in iterations to use with the simFast algorithm
#'@param algoOpts a list containing all the options for the specific algorithm you want to use.  Must match the algorithm given in the `algo` argument
#'@param metricOpts a list containing all the options for the specific metric you want to use.  Must match the metric given in the `metric` argument
#'@param suppressProg a parameter to suppress the progress bar. Mostly this is just used for creating documentation with knitr
#'@examples \dontrun{
#' 
#' ## Run the null model
#' finchMod <- cooc_null_model(dataWiFinches, algo="sim1",nReps=1000000,burn_in = 500)
#' ## Summary and plot info
#' summary(finchMod)
#' plot(finchMod,type="burn_in")
#' plot(finchMod,type="hist")
#' plot(finchMod,type="cooc")
#'}
#'
#'@export





sim9 <- function (speciesData,algo,metric, nReps = 1000 ,rowNames = TRUE, saveSeed = FALSE,burn_in = 0,suppressProg = TRUE)
{
  
  if(saveSeed){
    randomSeed <- .Random.seed
  } else {
    randomSeed <- NULL
  }

  ## Convert to matrix for type consistency
  if(!is.matrix(speciesData)){ speciesData <- as.matrix(speciesData)}
  
  ### Check for row names hidden in the data frame and automagically strip them.
  
  if(suppressWarnings(is.na(as.numeric(speciesData[2,1])))){
    speciesData <- speciesData[,-1] 
    class(speciesData) <- "numeric"
  }
  
  
  Start.Time <- Sys.time()
  metricF <- get(metric)
  
  Obs <- metricF(speciesData)
  #Trim the matrix to be just rowssums > 0 
  msim <- speciesData[rowSums(speciesData) > 0, ]
  ifelse(burn_in == 0, burn_in <- max(1000,10*nrow(msim)),burn_in <- burn_in)
  burn.in.metric <- vector(mode="numeric",length = burn_in)
  simulated.metric <- vector(mode="numeric",length = nReps)
  # run sequential swap for burn in series
  if(suppressProg){
    bi_pb <- txtProgressBar(min = 0, max = nReps, style = 3, file = stderr())
  } else{
    cat("Burn-in Progress \n")
    bi_pb <- txtProgressBar(min = 0, max = nReps, style = 3)
  }
  for (i in 1:burn_in)
  {
    msim <-sim9_single(msim)
    burn.in.metric[i] <- metricF(msim)
    setTxtProgressBar(bi_pb, i)
  }
  close(bi_pb)
  # run sequential swap for simulated series
  if(suppressProg){
    cat("Swap Progress \n")
    stat_pb <- txtProgressBar(min = 0, max = nReps, style = 3, file = stderr())
  } else{
    cat("Burn-in Progress \n")
    stat_pb <- txtProgressBar(min = 0, max = nReps, style = 3)
  }
  for (i in 1: nReps)
  {
    msim <-sim9_single(msim)
    simulated.metric[i] <- metricF(msim)    
    setTxtProgressBar(stat_pb, i)
  }
  close(stat_pb)
  
  Sim <- simulated.metric
  End.Time <- Sys.time()
  Elapsed.Time <- format(End.Time-Start.Time,digits=2)
  Time.Stamp <- date()

  sim9.fast.out <- list(Obs=Obs,Sim=Sim, Elapsed.Time=Elapsed.Time, Time.Stamp=Time.Stamp,Metric = metric, Algorithm = algo, N.Reps = nReps, SaveSeed = saveSeed, RandomSeed = randomSeed,Randomized.Data = msim , Data = speciesData,burn.in = burn_in,burn.in.metric= burn.in.metric)
  # plot to screen the trace function for the burn in
  
  class(sim9.fast.out) <- "nullmod"
  return(sim9.fast.out)
}



#' sim9_single 
#' @description Function for a single iteration of the sequential swap
#' @param speciesData binary presence-absence matrix
#' @export
sim9_single <- function (speciesData = matrix(rbinom(100, 1, 0.5), nrow = 10)) 
{
  # select two random rows and create submatrix
  ran.rows <- sample.int(nrow(speciesData), 2)
  m.pair <- speciesData[ran.rows, ]
  
  # find columns if any in pair for which colsum =1; these can be swapped
  Sum.Is.One = colSums(m.pair) == 1
  
  # Only swap if there are two or more columns to swap
  if(sum(Sum.Is.One) > 1){
    columns <- which(Sum.Is.One)
    
    # return swap entries in the two rows
    speciesData[ran.rows, columns] <- m.pair[, sample(columns)]
  }
  
  return(speciesData)
}


