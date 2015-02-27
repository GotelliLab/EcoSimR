#' Sim9.fast function
#' @description Special implementation of sequential swap algorithm
#' @export




sim9 <- function (speciesData,algo,metric, nReps = 1000 ,rowNames = TRUE, saveSeed = FALSE,burn_in = 0)
{
  
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
  
  Start.Time <- Sys.time()
  metricF <- get(metric)
  
  Obs <- metricF(speciesData)
  #Trim the matrix to be just rowssums > 0 
  msim <- speciesData[rowSums(speciesData) > 0, ]
  ifelse(burn_in == 0, burn_in <- max(1000,10*nrow(msim)),burn_in <- burn_in)
  burn.in.metric <- vector(mode="numeric",length = burn_in)
  simulated.metric <- vector(mode="numeric",length = nReps)
  # run sequential swap for burn in series
  cat("Burn-in Progress \n")
  
  bi_pb <- txtProgressBar(min = 0, max = burn_in, style = 3)
  for (i in 1:burn_in)
  {
    msim <-sim9_single(msim)
    burn.in.metric[i] <- metricF(msim)
    setTxtProgressBar(bi_pb, i)
  }
  close(bi_pb)
  # run sequential swap for simulated series
  cat("Swap Progress \n")
  
  stat_pb <- txtProgressBar(min = 0, max = nReps, style = 3)
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



#' sim9_single function
#' @description Function for a single iteration of the fast swap
#' @export
sim9_single <- function (m = matrix(rbinom(100, 1, 0.5), nrow = 10)) 
{
  # select two random rows and create submatrix
  ran.rows <- sample.int(nrow(m), 2)
  m.pair <- m[ran.rows, ]
  
  # find columns if any in pair for which colsum =1; these can be swapped
  Sum.Is.One = colSums(m.pair) == 1
  
  # Only swap if there are two or more columns to swap
  if(sum(Sum.Is.One) > 1){
    columns <- which(Sum.Is.One)
    
    # return swap entries in the two rows
    m[ran.rows, columns] <- m.pair[, sample(columns)]
  }
  
  return(m)
}


