#' Sim9.fast function
#' @description Special implementation of sequential swap algorithm
#' @export




sim9.fast <- function (speciesData,algo,metric, nReps = 1000 ,rowNames = TRUE, saveSeed = FALSE,burnin = 0)
{
  
  if(saveSeed){
    randomSeed <- .Random.seed
  } else {
    randomSeed <- NULL
  }

  
  
  ### Strip out row names
  
  if(row.names){
    speciesData <- speciesData[,-1] 
  }
  ## Convert to matrix for type consistency
  if(!is.matrix(speciesData)){ speciesData <- as.matrix(speciesData)}
  
  Start.Time <- Sys.time()
  Metric <- eval(parse(text = metric))
  
  Obs <- Metric(speciesData)
  msim <- speciesData
  ifelse(burnin ==0, burnin <- max(1000,10*nrow(msim)),burnin <- burnin)
  burn.in.metric <- vector(mode="numeric",length = burnin)
  simulated.metric <- vector(mode="numeric",length = n.reps)
  # run sequential swap for burn in series
  cat("Burn-in Progress \n")
  
  bi_pb <- txtProgressBar(min = 0, max = burnin, style = 3)
  for (i in 1:burnin)
  {
    msim <-sim9.single(msim)
    burn.in.metric[i] <- Metric(msim)
    setTxtProgressBar(bi_pb, i)
  }
  close(bi_pb)
  # run sequential swap for simulated series
  cat("Swap Progress \n")
  
  stat_pb <- txtProgressBar(min = 0, max = n.reps, style = 3)
  for (i in 1: n.reps)
  {
    msim <-sim9.single(msim)
    simulated.metric[i] <- Metric(msim)    
    setTxtProgressBar(stat_pb, i)
  }
  close(stat_pb)
  
  Sim <- simulated.metric
  End.Time <- Sys.time()
  Elapsed.Time <- format(End.Time-Start.Time,digits=2)
  Time.Stamp <- date()

  ### Reverse engineers the naming for consistent output
  ### Be sure to update the code below if new algos and metrics are added
  
  a.choice <- c("Uniform.Size", "Uniform.Size.User", "Source.Pool", "Gamma","ra1","ra2","ra3","ra4",
                paste("sim",1:10,sep=""),"simFast")
  a.func <- c("uniform_size", "uniform_size_user", "source_pool_draw", "Gamma","ra1","ra2","ra3","ra4",
              paste("sim",1:10,sep=""),"simFast")
  
  m.choice <- c("Min.Diff", "Min.Ratio", "Var.Diff", "Var.Ratio","Pianka", "Czekanowski", "Pianka.var", "Czekanowski.var", "Pianka.skew", "Czekanowski.skew",
                "Species.Combo", "Checker", "C.Score", "C.Score.var", "C.Score.skew", "V.Ratio")
  m.func <- c("min_diff", "min_ratio", "var_diff", "var_ratio","pianka", "czekanowski", "pianka_var", "czekanowski_var", "pianka_skew", "czekanowski_skew",
              "species_combo", "checker", "c_score", "c_score_var", "c_score_skew", "v_ratio")
  
  metric <- m.choice[which(m.func==metric)]
  algo <- a.choice[which(a.func==algo)]
  
  
  
  sim9.fast.out <- list(Obs=Obs,Sim=Sim, Elapsed.Time=Elapsed.Time, Time.Stamp=Time.Stamp,Metric = metric, Algorithm = algo, N.Reps = n.reps, SaveSeed = saveSeed, RandomSeed = randomSeed, Data = speciesData,burn.in = burnin,burn.in.metric=burn.in.metric)
  # plot to screen the trace function for the burn in
  
  class(sim9.fast.out) <- "nullmod"
  return(sim9.fast.out)
}



#' Sim9.Single function
#' @description Function for a single iteration of the fast swap
#' @export
sim9.single <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  
  # select two random rows and create submatrix
  ran.rows <- sample(nrow(m),2)
  
  m.pair <- m[ran.rows,]
  
  
  
  # find columns if any in pair for which colsum =1; these can be swapped
  
  Sum.Is.One <- which(colSums(m.pair)==1)
  
  if(length(Sum.Is.One)>1)
  {
    # swap them in m.pair.swapped
    m.pair.swapped <- m.pair
    m.pair.swapped[,Sum.Is.One] <- m.pair[,sample(Sum.Is.One)]
    
    # return the two swapped rows to the original m matrix
    
    m[ran.rows,] <- m.pair.swapped
    
  }
  
  return(m)
}

