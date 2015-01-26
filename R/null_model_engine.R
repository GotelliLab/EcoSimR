#'Run null model
#'@description This drives the null models for all the different kinds of null models that can be run. It is the underlying engine.
#'@param speciesData a dataframe <put some guidelines in here>
#'@param algo the algorithm to use, must be "RA1", "RA2", "RA3", "RA4"
#'@param metric the metric used to caluclate the null model: choices are "Pianka", "Czekanowski", "Pianka.var", "Czekanowski.var", "Pianka.skew", "Czekanowski.skew"; default is Pianka
#'@param nReps the number of replicates to run the null model.
#'@param rowNames Does your dataframe have row names? If yes, they are stripped, otherwise FALSE for data that has no row names
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


null_model_engine <- function(speciesData, algo, metric, nReps = 1000, rowNames = TRUE, randomSeed = 0,algoOpts = NULL, metricOpts = NULL)
{
  pboptions(type="txt",char="=")
  ## Set the seed
  ifelse (randomSeed==0, randomInteger <- trunc(runif(1,-2000000000,2000000000)), randomInteger <- randomSeed)
  
  set.seed(randomInteger)
  
  ### Strip out row names
  
  if(rowNames){
    speciesData <- speciesData[,-1] 
  }
  ## Convert to matrix for type consistency
  if(!is.matrix(speciesData)){ speciesData <- as.matrix(speciesData)}
  
  ### check if rownames accidentally was set to FALSE
  if(is.character(speciesData)){stop("Did you forget to set rowNames to TRUE?  Your data is non-numeric")}
  
  
  algoF <- eval(parse(text = algo))
  metricF <- eval(parse(text = metric))
  
  
  startTime <- Sys.time()
  
  
  if (nReps < 2) nReps <- 2
  
  Sim <- rep(NA,nReps)
  for(i in 1:nReps){
    algo.opts["speciesData"] <- 
    z <- do.call(algoF,)
  }
  
  sim <- pbreplicate(nReps,metricF(algoF(speciesData)))
  obs <- metricF(speciesData)
  
  endTime <- Sys.time()
  elapsedTime <- format(endTime-startTime,digits=2)
  timeStamp <- date()
  
  ### Reverse engineers the naming for consistent output
  ### Be sure to update the code below if new algos and metrics are added
  
  aChoice <- c("Uniform.Size", "Uniform.Size.User", "Source.Pool", "Gamma","ra1","ra2","ra3","ra4",
                paste("sim",1:10,sep=""),"simFast")
  aFunc <- c("uniform_size", "uniform_size_user", "source_pool_draw", "Gamma","ra1","ra2","ra3","ra4",
              paste("sim",1:10,sep=""),"simFast")
  
  mChoice <- c("Min.Diff", "Min.Ratio", "Var.Diff", "Var.Ratio","Pianka", "Czekanowski", "Pianka.var", "Czekanowski.var", "Pianka.skew", "Czekanowski.skew",
                "Species.Combo", "Checker", "C.Score", "C.Score.var", "C.Score.skew", "V.Ratio")
  mFunc <- c("min_diff", "min_ratio", "var_diff", "var_ratio","pianka", "czekanowski", "pianka_var", "czekanowski_var", "pianka_skew", "czekanowski_skew",
              "species_combo", "checker", "c_score", "c_score_var", "c_score_skew", "v_ratio")
  
  metric <- mChoice[which(mFunc==metric)]
  algo <- aChoice[which(aFunc==algo)]
  
  
  nullModelOut <- list(Obs=obs,Sim=sim, Elapsed.Time=elapsedTime, Time.Stamp=timeStamp,Metric = metric, Algorithm = algo, nReps = nReps, Random.Integer = randomInteger, Data = speciesData)
  class(nullModelOut) <- "nullmod"
  return(nullModelOut)
  
}

