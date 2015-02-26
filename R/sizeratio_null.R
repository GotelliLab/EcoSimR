#'Size Ratio 
#'@description Create a size Ratio null model
#'@param speciesData a dataframe <put some guidelines in here>
#'@param algo the algorithm to use, must be "size_uniform", "size_uniform_user", "size_source_pool", "size_gamma"
#'@param metric the metric used to caluclate the null model: choices are "min_diff", "min_ratio", "var_diff", "var_ratio"; default is Var.Ratio
#'@param nReps the number of replicates to run the null model.
#'@param rowNames Does your dataframe have row names? If yes, they are stripped, otherwise FALSE for data that has no row names
#'@param saveSeed TRUE or FALSE.  If TRUE the current seed is saved so the simulation can be repeated
#'@param algoOpts a list containing all the options for the specific algorithm you want to use.  Must match the algorithm given in the `algo` argument
#'@param metricOpts a list containing all the options for the specific metric you want to use.  Must match the metric given in the `metric` argument
#'@examples \dontrun{
#' ## Run the null model
#' rodentMod <- size_null_model(dataRodents)
#' ## Summary and plot info
#' summary(rodentMod)
#' plot(rodentMod,type="hist")
#' plot(rodentMod,type="size")
#' 
#' ##  Uniform Size model with user inputs
#' rodentMod2 <- size_null_model(dataRodents,algo="size_uniform_user",
#' algoOpts = list(userLow = 3,userHigh=15))
#' summary(rodentMod2)
#' plot(rodentMod2,type="hist")
#' plot(rodentMod2,type="size")
#' 
#' ### Source pool model
#' 
#' rodentMod_sp <- size_null_model(dataRodents,algo="size_source_pool",
#' algoOpts = list(sourcePool = runif(dim(dataRodents)[1],1,15)))
#' 
#' summary(rodentMod_sp)
#' plot(rodentMod_sp,type="hist")
#' plot(rodentMod_sp,type="size")
#' 
#'}
#'
#'@export

size_null_model <- function(speciesData, algo = "size_uniform", metric = "var_ratio", nReps = 1000, rowNames = TRUE, saveSeed = FALSE, algoOpts = list(), metricOpts = list()){

  mChoice<- c("min_diff", "min_ratio", "var_diff", "var_ratio")
  aChoice <- c("size_uniform", "size_uniform_user", "size_source_pool", "size_gamma")
  
  algo <- match.arg(algo,choices = aChoice)
  metric <- match.arg(metric,choices = mChoice)
  
  params <- list(speciesData = speciesData, algo = algo, metric = metric, nReps = nReps, rowNames = rowNames, saveSeed = saveSeed, algoOpts = algoOpts, metricOpts = metricOpts)
  output <- do.call(null_model_engine,params)
  class(output) <- "sizenullmod"
  return(output)
  
}


#' Generic function for calculating null model summary statistics
#' @description Takes as input a list of Null.Model.Out, with Obs, Sim, Elapsed Time, and Time Stamp values
#' @export

summary.sizenullmod <- function(object,...)
{ 
  nullmodObj <- object
  #if (!is.null(Output.File)) outfile <- file(p$Output.File, "w") else outfile <-""
  
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


plot.sizenullmod <- function(x, type = "hist",...)
{
cat("in here")
nullmodObj <- x
if(type == "hist"){
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

if(type=="size"){


  opar<- par(no.readonly=TRUE)
  par(mfrow=c(2,1))

  One.Null.Vector <- nullmodObj$Randomized.Data
  
  limits <- range(c(nullmodObj$Data,One.Null.Vector))
  x.lims <- c((limits[1]-0.1*limits[2]),(1.1*limits[2]))
  if(x.lims[1]<0)x.lims[1] <- 0
  plot(nullmodObj$Data,xlim=x.lims,ylim=c(0,5),axes=TRUE,ann=TRUE,ylab="",
       xlab="Body Size",cex.lab=1.5,cex.axis=1.5,yaxt="n",bty="n") 
  lines(x=x.lims,y=c(2,2))
  points(x=nullmodObj$Data,y=rep(2,length(nullmodObj$Data)),pch=21,bg="red3",cex=2)
  lines(x=x.lims,y=c(4,4))
  points(x=One.Null.Vector,y=rep(4,length(One.Null.Vector)),pch=21,bg="royalblue3",cex=2)
  text.center <- mean(x.lims)
  text(text.center,1,"Observed",font=2,cex=1.2)
  text(text.center,3,"Simulated",font=2,cex=1.2)
  
  
  
  mtext(as.character(nullmodObj$Time.Stamp),side=3,adj=1,line=3)
  Segs.Obs <-sort(diff(sort(nullmodObj$Data)))
  
  Segs.Sim <- sort(diff(sort(One.Null.Vector)))
  
  Segs.Mat <- matrix(c(Segs.Obs,Segs.Sim),ncol=2)
  
  barplot(Segs.Mat,beside=TRUE,col=rep(c("red3","royalblue3"),each=nrow(Segs.Mat)),
          names.arg=c("Observed","Simulated"),ylab="Body Size Difference",font=2,
          cex.names=1.5,sub="Sorted Size Differences",cex.sub=1.5,font.sub=2)
  
  
  par(opar)

}

}

