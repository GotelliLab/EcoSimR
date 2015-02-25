#'Co-Occurrence Null model 
#'@description Create a Co-Occurrence null model
#'@param speciesData a dataframe in which rows are species, columns are sites,
#' and the entries indicate the absence (0) or presence (1) of a species in a 
#' site. Empty rows and empty columns should not be included in the matrix.
#'@param algo the algorithm to use, must be "sim1", "sim2", "sim3", "sim4", "sim5", "sim6", "sim7", "sim8", "sim9", "sim10"
#'@param metric the metric used to caluclate the null model: choices are "species_combo", "checker", "c_score", "c_score_var", "c_score_skew", "v_ratio"; default is "c_score"
#'@param nReps the number of replicates to run the null model.
#'@param rowNames Does your dataframe have row names? If yes, they are stripped, otherwise FALSE for data that has no row names
#'@param saveSeed TRUE or FALSE.  If TRUE the current seed is saved so the simulation can be repeated
#'@param burn_in The number of burn_in iterations to use with the simFast algorithm
#'@param algoOpts a list containing all the options for the specific algorithm you want to use.  Must match the algorithm given in the `algo` argument
#'@param metricOpts a list containing all the options for the specific metric you want to use.  Must match the metric given in the `metric` argument
#'@examples \dontrun{
#' 
#' ## Run the null model
#' finchMod <- cooc_null_model(dataWiFinches, algo="sim9")
#' ## Summary and plot info
#' summary(finchMod)
#' plot(finchMod,type="burn_in")
#' 
#' ## Example that is repeatable with a saved seed
#' finchMod <- cooc_null_model(dataWiFinches, algo="sim1",saveSeed = TRUE)
#' mean(finchMod$Sim)
#' ## Run the model with the seed saved
#' 
#' finchMod <- cooc_null_model(dataWiFinches, algo="sim1",saveSeed=T)
#' ## Check model output
#' mean(finchMod$Sim)
#' 
#' reproduce_model(finchMod$Sim)
#' 
#' finchMod <- cooc_null_model(dataWiFinches, algo="sim1")
#' ## Check model output is the same as before
#' mean(finchMod$Sim)
#' reproduce_model(finchMod$Sim)
#' 
#' 
#'}
#'
#'@export

cooc_null_model <- function(speciesData, algo = "sim1", metric = "c_score", nReps = 1000, rowNames = TRUE, saveSeed = FALSE, burn_in = 0,algoOpts = list(),metricOpts = list()){
  aChoice <- c(paste("sim",c(1:10),sep=""))
  mChoice <- c("species_combo", "checker", "c_score", "c_score_var", "c_score_skew", "v_ratio")

  algo <- match.arg(algo,choices = aChoice)
  metric <- match.arg(metric,choices = mChoice)
  ## Control behavior of whether or not sim9fast is used.
  if(algo != "sim9"){
  params <- list(speciesData = speciesData, algo = algo, metric = metric, nReps = nReps, rowNames = rowNames, saveSeed
 = saveSeed
  ,algoOpts = algoOpts,metricOpts = metricOpts)
  output <- do.call(null_model_engine,params)
  output$burn.in <- burn_in
  class(output) <- "coocnullmod"
  return(output)
  } else if(algo == "sim9"){
    params <- list(speciesData = speciesData,algo = algo, metric = metric, nReps = nReps, rowNames = rowNames, saveSeed = saveSeed, burn_in = burn_in)
    output <- do.call(sim9,params)
    class(output) <- "coocnullmod"
    return(output)
  }
  
  
}


#' Generic function for calculating null model summary statistics
#' @description Takes as input a list of Null.Model.Out, with Obs, Sim, Elapsed Time, and Time Stamp values
#' @export

summary.coocnullmod <- function(object,...)
{ 
  nullmodObj <- object
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



plot.coocnullmod <- function(x, type = "hist",...)
{
  nullmodObj <- x 
  
  if(type == "cooc"){
  Date.Stamp=date()
  par(mfrow=c(1,2))
  if(nullmodObj$Algorithm!="sim9"){
  Fun.Alg <- get(nullmodObj$Algorithm)
  } else {
    Fun.Alg <- sim9_single
  }
  
  One.Null.Matrix <- Fun.Alg(nullmodObj$Data)
  
  # reverse the matrix rows for plotting consistency
  m <- One.Null.Matrix
  m <- m[rev(1:nrow(m)),]
  
  # setup plotting space
  
  plot(m,xlim=c(0,ncol(m)),ylim=c(0,nrow(m)),type="n",ann=FALSE,axes=FALSE)
  mtext("Sites",side=1,font=2)
  mtext("Species",side=2,font=2)
  mtext("Simulated",side=3,font=2,col="royalblue3")
  # define coordinate vectors
  yrec <- rep(0:(nrow(m)-1),ncol(m))
  xrec <- rep(0:(ncol(m)-1),each=nrow(m))
  
  # Set up color labels
  Plot.cols <- c("white","royalblue3")
  Color.Vector <- Plot.cols[as.integer(m)+1]
  
  # Plot and fill rectangles
  rect(xrec,yrec,xrec+1,yrec+1,col=Color.Vector)
  
  
  mtext(as.character(Date.Stamp),side=3,adj=1,line=3)
  # reverse the matrix rows for plotting consistency
  m <- nullmodObj$Data
  m <- m[rev(1:nrow(m)),]
  
  # setup plotting space
  
  plot(m,xlim=c(0,ncol(m)),ylim=c(0,nrow(m)),type="n",
       ann=FALSE,axes=FALSE)
  mtext("Sites",side=1,font=2)
  mtext("Species",side=2,font=2)
  mtext("Observed",side=3,font=2,col="red3")
  # define coordinate vectors
  yrec <- rep(0:(nrow(m)-1),ncol(m))
  xrec <- rep(0:(ncol(m)-1),each=nrow(m))
  
  # Set up color labels
  Plot.cols <- c("white","red3")
  Color.Vector <- Plot.cols[as.integer(m)+1]
  
  # Plot and fill rectangles
  rect(xrec,yrec,xrec+1,yrec+1,col=Color.Vector)
}

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

if(type=="burn_in"){
  if(is.na(nullmodObj$burn.in)){
    warning("You can only create a burn_in plot for a model run with the 'simFast' algorithm")
   
  }
  par(mfrow=c(1,1))
  v <- nullmodObj$burn.in.metric
  z <- nullmodObj$Obs
  v <- c(z,v)
  plot(x=1:length(v),y=v,xlab="Iteration",ylab="Index",
       las=1,type="l",col="royalblue3")
  abline(h=z,col="red3")
  lines(lowess(1:length(v),v), col="gray",lwd=4) # lowess line (x,y) 
  
}

}
