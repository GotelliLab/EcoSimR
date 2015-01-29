#'Co-Occurrence Null model 
#'@description Create a Co-Occurrence null model
#'@param species_data a dataframe <put some guidelines in here>
#'@param algo the algorithm to use, must be "sim1", "sim2", "sim3", "sim4", "sim5", "sim6", "sim7", "sim8", "sim9", "sim10"
#'@param metric the metric used to caluclate the null model: choices are "Species.Combo", "Checker", "C.Score", "C.Score.var", "C.Score.skew", "V.Ratio"; default is "C.Score"
#'@param n.reps the number of replicates to run the null model.
#'@param rowNames Does your dataframe have row names? If yes, they are stripped, otherwise FALSE for data that has no row names
#'@param random.seed Choose a seed to start your random number.  0 will choose a random seed, otherwise set the seed with any integer.
#'@param burnin The number of burnin iterations to use with the simFast algorithm
#'@examples \dontrun{
#' 
#' ## Run the null model
#' finchMod <- cooc_null_model(wiFinches, algo="sim1")
#' ## Summary and plot info
#' summary(finchMod)
#' plot(finchMod,type="burnin")
#'  
#'}
#'
#'@export

cooc_null_model <- function(species_data, algo = "simFast", metric = "C.Score", n.reps = 1000, rowNames = TRUE, random.seed = 0, burnin = 0,algoOpts = list(),metricOpts = list(), row.names = TRUE){
  mChoice <- c("Species.Combo", "Checker", "C.Score", "C.Score.var", "C.Score.skew", "V.Ratio")
  aChoice <- c(paste("sim",1:10,sep=""),"simFast")
  mFunc <- c("species_combo", "checker", "c_score", "c_score_var", "c_score_skew", "v_ratio")

  algo <- match.arg(algo,choices = aChoice)
  metric <- match.arg(metric,choices = mChoice)
  metric <- mFunc[which(mChoice==metric)]
  ## Control behavior of whether or not sim9fast is used.
  if(algo != "simFast"){
  params <- list(species_data = species_data, algo = algo, metric = metric, n.reps = n.reps, rowNames = rowNames, random.seed
 = random.seed
  ,algoOpts = algoOpts,metricOpts = metricOpts)
  output <- do.call(null_model_engine,params)
  output$burn.in <- burnin
  class(output) <- "coocnullmod"
  return(output)
  } else if(algo == "simFast"){
    params <- list(species_data = species_data,algo = algo, metric = metric, n.reps = n.reps, row.names = row.names, random.seed = random.seed, burnin = burnin)
    output <- do.call(sim9.fast,params)
    class(output) <- "coocnullmod"
    return(output)
  }
  
  
}


#' Generic function for calculating null model summary statistics
#' @description Takes as input a list of Null.Model.Out, with Obs, Sim, Elapsed Time, and Time Stamp values
#' @export

summary.coocnullmod <- function(nullmodObj)
{ 
  
  
  
  
  #if (!is.null(Output.File)) outfile <- file(p$Output.File, "w") else outfile <-""
  
  cat("Time Stamp: " , nullmodObj$Time.Stamp,   "\n") 
  # cat("Data File: ", p$Data.File,  "\n")
  #  cat("Output File: ", p$Output.File,  "\n") 
  cat("Random Number Seed: ",nullmodObj$RandomInteger,  "\n")
  cat("Number of Replications: ",nullmodObj$n.reps,  "\n")
  cat("Elapsed Time: ", nullmodObj$Elapsed.Time, "\n")
  cat("Metric: ", nullmodObj$MetricOut,  "\n")
  cat("Algorithm: ", nullmodObj$AlgorithmOut,  "\n") 
  
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



plot.coocnullmod <- function(nullmodObj, type = "hist")
{
  if(type == "cooc"){
  Date.Stamp=date()
  par(mfrow=c(1,2))
  if(is.na(nullmodObj$burn.in)){
  Fun.Alg <- eval(parse(text = nullmodObj$Algorithm))
  } else {
    Fun.Alg <- sim9.single
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
  par (mar=c(5,6,4,2)+0.1)
  #------------------------------------------------------
  hist(nullmodObj$Sim, breaks=20, col="royalblue3",
       
       xlab="Simulated Metric",ylab="Frequency",main="",
       xlim=range(c(nullmodObj$Sim,nullmodObj$Obs)))
  abline(v=nullmodObj$Obs,col="red",lty="solid",lwd=2.5)
  abline(v=quantile(nullmodObj$Sim,c(0.05,0.95)),col="black",lty="dashed",lwd=2.5)
  abline(v=quantile(nullmodObj$Sim,c(0.025,0.975)),col="black",lty="dotted",lwd=2.5)
  mtext(as.character(date()),side=3,adj=1,line=3)
  
}

if(type=="burnin"){
  if(is.na(nullmodObj$burn.in)){
    warning("You can only create a burnin plot for a model run with the 'simFast' algorithm")
   
  }
  v <- nullmodObj$burn.in.metric
  z <- nullmodObj$Obs
  v <- c(z,v)
  plot(x=1:length(v),y=v,xlab="Iteration",ylab="Index",
       las=1,type="l",col="royalblue3")
  abline(h=z,col="red3")
  lines(lowess(1:length(v),v), col="gray",lwd=4) # lowess line (x,y) 
  
}

}
