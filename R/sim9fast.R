#' Sim9.fast function
#' @description Special implementation of sequential swap algorithm
#' @export




sim9.fast <- function (species_data,algo,metric, n.reps = 1000 ,row.names = TRUE, random.seed = 0,burnin = 0)
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
  
  Start.Time <- Sys.time()
  Metric <- eval(parse(text = metric))
  
  Obs <- Metric(species_data)
  msim <- species_data
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
  sim9.fast.out <- list(Obs=Obs,Sim=Sim, Elapsed.Time=Elapsed.Time, Time.Stamp=Time.Stamp,Metric = metric, Algorithm = algo, N.Reps = n.reps, RandomInteger = RandomInteger, Data = species_data,burn.in = burnin,burn.in.metric=burn.in.metric)
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


###############################################################
# Burn.In Plot Function
# Plots the trace of burn in values for sequential swap
Burn.In.Plot <- function(v=runif(1000),z=0.9)
{
  v <- c(z,v)
  plot(x=1:length(v),y=v,xlab="Iteration",ylab="Index",
       las=1,type="l",col="royalblue3")
  abline(h=z,col="red3")
  lines(lowess(1:length(v),v), col="gray",lwd=4) # lowess line (x,y) 
}
