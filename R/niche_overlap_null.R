#'Niche overlap 
#'@description Create a niche overlap null model
#'@param speciesData a dataframe <put some guidelines in here>
#'@param algo the algorithm to use, must be "ra1", "ra2", "ra3", "ra4"
#'@param metric the metric used to caluclate the null model: choices are "pianka", "czekanowski", "pianka_var", "czekanowski_var", "pianka_skew", "czekanowski_skew"; default is pianka
#'@param nReps the number of replicates to run the null model.
#'@param rowNames Does your dataframe have row names? If yes, they are stripped, otherwise FALSE for data that has no row names
#'@param saveSeed TRUE or FALSE.  If TRUE the current seed is saved so the simulation can be repeated
#'@param algoOpts a list containing all the options for the specific algorithm you want to use.  Must match the algorithm given in the `algo` argument
#'@param metricOpts a list containing all the options for the specific metric you want to use.  Must match the metric given in the `metric` argument
#'@examples \dontrun{
#' ## Load MacAruthur warbler data
#' data(dataMacWarb)
#' 
#' ## Run the null model
#' warbMod <- niche_null_model(dataMacWarb)
#' ## Summary and plot info
#' summary(warbMod)
#' plot(warbMod)
#' plot(warbMod,type="niche")
#'  
#'}
#'
#'@export

niche_null_model <- function(speciesData, algo = "ra3", metric = "pianka", nReps = 1000, rowNames = TRUE,algoOpts = list(),metricOpts = list(),saveSeed=TRUE){
  aChoice <- c("ra1","ra2","ra3","ra4")
  mChoice<- c("pianka", "czekanowski", "pianka_var", "czekanowski_var", "pianka_skew", "czekanowski_skew")
  
  algo <- match.arg(algo,choices = aChoice)
  metric <- match.arg(metric,choices = mChoice)
  
  params <- list(speciesData = speciesData, algo = algo, metric = metric, nReps = nReps, rowNames = rowNames, saveSeed = saveSeed,algoOpts = algoOpts,metricOpts = metricOpts)
  output <- do.call(null_model_engine,params)
  class(output) <- "nichenullmod"
  return(output)
  
}


#' Generic function for calculating null model summary statistics
#' @description Takes as input a list of Null.Model.Out, with Obs, Sim, Elapsed Time, and Time Stamp values
#' @export

summary.nichenullmod <- function(object,...)
{ 
 
  nullmodObj <- object 
  
  cat("Time Stamp: " , nullmodObj$Time.Stamp,   "\n")  
  cat("Reproducible: ",nullmodObj$SaveSeed,  "\n")
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




#' Null Model Plot function
#' @description Generic function for plotting a histogram of simulated values Takes as input a list of Null.Model.Out, with Obs and Sim values
#' @export



plot.nichenullmod <- function(x, type = "hist",...)
{
  nullmodObj <- x
  if(type == "hist"){

  opar <- par(no.readonly=TRUE)
  par(cex=1, cex.axis = 1.5,
      cex.main=1,cex.lab=1.6)
  par (mar=c(5,6,4,2)+0.1,mfrow=c(1,1))
  #------------------------------------------------------
  hist(nullmodObj$Sim, breaks=20, col="royalblue3",
       
       xlab="Simulated Metric",ylab="Frequency",main="",
       xlim=range(c(nullmodObj$Sim,nullmodObj$Obs)))
  abline(v=nullmodObj$Obs,col="red",lty="solid",lwd=2.5)
  abline(v=quantile(nullmodObj$Sim,c(0.05,0.95)),col="black",lty="dashed",lwd=2.5)
  abline(v=quantile(nullmodObj$Sim,c(0.025,0.975)),col="black",lty="dotted",lwd=2.5)
  mtext(as.character(date()),side=3,adj=1,line=3)
  }
  
  if(type == "niche"){
    opar<- par(no.readonly=TRUE)
    par(mfrow=c(2,1))
    Data <- nullmodObj$Data/rowSums(nullmodObj$Data)
    plot(rep(1:ncol(Data),times = nrow(Data)),
         rep(1:nrow(Data),each=ncol(Data)),
         xlab="Resource Category",ylab="Species",cex=10*sqrt(t(Data)/pi),col="red3",lwd=2,
         main="Observed Utilization Matrix",col.main="red3",cex.main=1.5)
  mtext(as.character(nullmodObj$Time.Stamp),side=3,adj=1,line=3)
    
    One.Null.Matrix <- nullmodObj$Randomized.Data
    One.Null.Matrix <- One.Null.Matrix/rowSums(One.Null.Matrix)
    plot(rep(1:ncol(One.Null.Matrix),times = nrow(One.Null.Matrix)),
         rep(1:nrow(One.Null.Matrix),each=ncol(One.Null.Matrix)),
         xlab="Resource Category",ylab="Species",cex=10*sqrt(t(One.Null.Matrix)/pi),col="royalblue3",lwd=2,
         main="Simulated Utilization Matrix",col.main="royalblue3",cex.main=1.5)
    par(opar)
  }
  
}
