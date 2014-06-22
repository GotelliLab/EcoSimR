#'Co-Occurrence Null model 
#'@description Create a Co-Occurrence null model
#'@param species_data a dataframe <put some guidelines in here>
#'@param algo the algorithm to use, must be  "spcombo" "checker" "C.Score" "C.Score.var" "C.Score.skew" "V.Ratio"
#'@param metric the metric used to caluclate the null model: choices are "Pianka", "Czekanowski", "Pianka.var", "Czekanowski.var", "Pianka.skew", "Czekanowski.skew"; default is Pianka
#'@param n.reps the number of replicates to run the null model.
#'@param row.names Does your dataframe have row names? If yes, they are stripped, otherwise FALSE for data that has no row names
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

cooc_null_model <- function(species_data, algo = "simFast", metric = "C.Score", n.reps = 1000, row.names = TRUE, random.seed = 0){
  m.choice <- c("Species.Combo", "Checker", "C.Score", "C.Score.var", "C.Score.skew", "V.Ratio")
  a.choice <- c(paste("sim",1:10,sep=""),"simFast")
  
  algo <- match.arg(algo,choices = a.choice)
  metric <- match.arg(metric,choices = m.choice)
  params <- list(species_data = species_data, algo = algo, metric = metric, n.reps = n.reps, row.names = row.names, random.seed = random.seed)
  output <- do.call(null_model_engine,params)
  class(output) <- "nichenullmod"
  return(output)
  
}


#' Generic function for calculating null model summary statistics
#' @description Takes as input a list of Null.Model.Out, with Obs, Sim, Elapsed Time, and Time Stamp values
#' @export

summary.nichenullmod <- function(nullmodObj)
{ 
  
  
  
  
  #if (!is.null(Output.File)) outfile <- file(p$Output.File, "w") else outfile <-""
  
  cat("Time Stamp: " , nullmodObj$Time.Stamp,   "\n") 
  # cat("Data File: ", p$Data.File,  "\n")
  #  cat("Output File: ", p$Output.File,  "\n") 
  cat("Random Number Seed: ",nullmodObj$RandomInteger,  "\n")
  cat("Number of Replications: ",nullmodObj$N.Reps,  "\n")
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



plot.nichenullmod <- function(nullmodObj)
{
  
  
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
