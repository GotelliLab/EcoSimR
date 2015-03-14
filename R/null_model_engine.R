#'Run null model
#'@description This drives the null models for all the different kinds of null models that can be run. It is the underlying engine.
#'@param speciesData a dataframe of data that will work with metrics and algorithms.
#'@param algo the algorithm to use
#'@param metric the metric used to caluclate the null model
#'@param nReps the number of replicates to run the null model.
#'@param saveSeed Should the existing random seed be saved to make the model reproducible? 
#'@param algoOpts a list containing options for a supplied alogrithm
#'@param metricOpts a list containing options for a supplied metric
#'@param type The type of null model you are running,  if it's meant to integrate with existing null models the type should be "size","niche!","cooc", or leave it NULL if you it's not meant to be compatible with existing null models and a generic null model will be returned.
#'@param suppressProg a parameter to suppress the progress bar. Mostly this is just used for creating documentation with knitr
#'@examples \dontrun{
#' # User defined function
#' 
#'
#'}
#'@export


null_model_engine <- function(speciesData, algo, metric, nReps = 1000, saveSeed = FALSE, algoOpts = list(), metricOpts = list(),type=NULL,suppressProg=FALSE)
{
  if(suppressProg){
    pb <- txtProgressBar(min = 0, max = nReps, style = 3, file = stderr())
  } else{
  pb <- txtProgressBar(min = 0, max = nReps, style = 3)
  }
  ## Set the seed
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

  algoF <- get(algo)
  metricF <- get(metric)
  
  ### Error check for input functions
  if(!grepl("speciesData",names(formals(algoF)[1]))){
    stop("Please enter a valid algorithm with 'speciesData' as the first parameter")
  }
  if(!grepl("m",names(formals(metricF)[1])) && nchar(names(formals(metricF)[1])) == 1 ){
    stop("Please enter a valid metric with 'm' as the first parameter")
  }
    
  
  
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
  

  
  nullModelOut <- list(Obs=obs,Sim=sim, Elapsed.Time=elapsedTime, Time.Stamp=timeStamp,Metric = metric, Algorithm = algo, nReps = nReps, 
                       Reproducible = saveSeed,RandomSeed = randomSeed, Data = speciesData,Randomized.Data = finalRandomData)
  
  if(is.null(type)){
  class(nullModelOut) <- "nullmod"
} else if (type %in% c("niche","cooc","size")){
  class(nullModelOut) <- paste(type,"nullmod",sep="")
  
  
}
  
  return(nullModelOut)
  
}



#' Generic function for calculating null model summary statistics
#' @description Takes a null model object and prints a nice summary.
#' @param object the null model object to print a summary of. 
#' @param ... Extra parameters for summary
#' @export

summary.nullmod <- function(object,...)
{ 
  nullmodObj <- object
  #if (!is.null(Output.File)) outfile <- file(p$Output.File, "w") else outfile <-""
  
  cat("Time Stamp: " , nullmodObj$Time.Stamp,   "\n") 
  cat("Reproducible: ",nullmodObj$Reproducible,  "\n")
  cat("Number of Replications: ",nullmodObj$nReps,  "\n")
  cat("Elapsed Time: ", nullmodObj$Elapsed.Time, "\n")
  cat("Metric: ", nullmodObj$Metric,  "\n")
  cat("Algorithm: ", nullmodObj$Algorithm,  "\n") 
  
  cat("Observed Index: ", format(nullmodObj$Obs,digits=5),  "\n")
  cat("Mean Of Simulated Index: ",format(mean(nullmodObj$Sim),digits=5),  "\n")
  cat("Variance Of Simulated Index: ",format(var(nullmodObj$Sim),digits=5),  "\n")
  cat("Lower 95% (1-tail): ",format(quantile(nullmodObj$Sim,0.05),digits=5),  "\n")
  cat("Upper 95% (1-tail): ",format(quantile(nullmodObj$Sim,0.95),digits=5), "\n")
  cat("Lower 95% (2-tail): ",format(quantile(nullmodObj$Sim,0.025),digits=5), "\n")
  cat("Upper 95% (2-tail): ",format(quantile(nullmodObj$Sim,0.975),digits=5),  "\n")
  
  #P-values
  if (nullmodObj$Obs > max(nullmodObj$Sim)) {
    cat("Lower-tail P < ",(length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim),  "\n")
    cat("Upper-tail P > ",1/length(nullmodObj$Sim),  "\n")
  } else if(nullmodObj$Obs < min(nullmodObj$Sim)) {
    cat("Lower-tail P > ", 1/length(nullmodObj$Sim), "\n")
    cat("Upper-tail P < ",(length(nullmodObj$Sim) - 1)/length(nullmodObj$Sim), "\n")
  } else {
    cat("Lower-tail P = ", format(sum(nullmodObj$Obs >= nullmodObj$Sim)/length(nullmodObj$Sim),digits=5),  "\n")
    cat("Upper-tail P = ", format(sum(nullmodObj$Obs <= nullmodObj$Sim)/length(nullmodObj$Sim),digits=5), "\n")
  }
  
  cat(paste("Observed metric > ",sum(nullmodObj$Obs > nullmodObj$Sim)," simulated metrics",sep="") , "\n")
  cat(paste("Observed metric < ",sum(nullmodObj$Obs < nullmodObj$Sim)," simulated metrics",sep="")  ,"\n")
  cat(paste("Observed metric = ",sum(nullmodObj$Obs == nullmodObj$Sim)," simulated metrics",sep="") , "\n")
  cat("Standardized Effect Size (SES): ", format((nullmodObj$Obs - mean(nullmodObj$Sim))/sd(nullmodObj$Sim),digits=5), "\n")
  
  #if(!is.null(Output.File)) close(outfile)
}


#' plot a histogram null model
#' @description plot a histogram of a generic null model object
#' @param x the null model to plot
#' @param ... Other variables to be passed on to base plotting
#' @details the valid types for size are "hist" to show a histogram and "size" to show a sample size null model.
#' @export



plot.nullmod <- function(x,...)
{
  nullmodObj <- x
    par(mfrow=c(1,1))
    opar <- par(no.readonly=TRUE)
    par(cex=1, cex.axis = 1.5,
        cex.main=1,cex.lab=1.6)
    par (mar=c(5,6,4,2)+0.1)
    hist(nullmodObj$Sim, breaks=20, col="royalblue3",
         
         xlab="Simulated Metric",ylab="Frequency",main="",
         xlim=range(c(nullmodObj$Sim,nullmodObj$Obs)))
    
    abline(v=nullmodObj$Obs,col="red",lty="solid",lwd=2.5)
    abline(v=quantile(nullmodObj$Sim,c(0.05,0.95)),col="black",lty="dashed",lwd=2.5)
    abline(v=quantile(nullmodObj$Sim,c(0.025,0.975)),col="black",lty="dotted",lwd=2.5)
    mtext(as.character(date()),side=3,adj=1,line=3)
  }


