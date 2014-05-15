##############################################
## EcoSimR: R code for Null Model Analysis
## 
##
##
##############################################
## EcoSimR Graphics Source File
## Nicholas J. Gotelli & Aaron M. Ellison
##
##
##############################################
## Version 1.00
## 12 May 2013
##############################################

##########################################################################################
## Niche.Overlap.Plot function
## Function for plotting observed and simulated utilization matrix
## Takes as input an input data matrix, a randomization algorithm, and the time stamp
##
##

Niche.Overlap.Plot <- function (Data=matrix(rpois(80,1),nrow=5), Algorithm="RA3", Date.Stamp=date(), Plot.Output="screen")
{
  
opar<- par(no.readonly=TRUE)
if (Plot.Output == "file") par(mfrow=c(2,1))
Data <- Data/rowSums(Data)
plot(rep(1:ncol(Data),times = nrow(Data)),
     rep(1:nrow(Data),each=ncol(Data)),
     xlab="Resource Category",ylab="Species",cex=10*sqrt(t(Data)/pi),col="red",lwd=2,
     main="Observed Utilization Matrix",col.main="red",cex.main=1.5)
if (Plot.Output=="file") mtext(as.character(Date.Stamp),side=3,adj=1,line=3)

Fun.Alg <- get(Algorithm)
One.Null.Matrix <- Fun.Alg(Data)
One.Null.Matrix <- One.Null.Matrix/rowSums(One.Null.Matrix)
plot(rep(1:ncol(One.Null.Matrix),times = nrow(One.Null.Matrix)),
     rep(1:nrow(One.Null.Matrix),each=ncol(One.Null.Matrix)),
     xlab="Resource Category",ylab="Species",cex=10*sqrt(t(One.Null.Matrix)/pi),col="royalblue3",lwd=2,
     main="Simulated Utilization Matrix",col.main="royalblue3",cex.main=1.5)
par(opar)
}