##############################################
## EcoSimR - R Code For Null Model Analysis
## 
##
##
##############################################
## EcoSimR Algorithms Source File
## Nicholas J. Gotelli & Aaron M. Ellison
##
##
##############################################
## Version 1.00
## 12 May 2013
#############################################

#' VectorSample function
#' @description Takes an input binary vector and a weight vector Reassigns 1s randomly in proportion to vector weights
#' @export

vector_sample <- function(speciesData,weights) 

{
x <- mat.or.vec(length(speciesData),1)                   # creates a vector of 0s
x[sample(1:length(speciesData),size=sum(speciesData),prob=weights)] <- 1  # fills with 1s, sampling with weights
return(x)
}


#' Sim1 Function
#' @description Randomizes a binary matrix speciesData by reshuffling all of its elements equiprobably
#' @export

sim1 <- function(speciesData) 

{
matrix(sample(speciesData), ncol=ncol(speciesData))
}



#' Sim2 Function
#' @description Randomizes a binary matrix speciesData by reshuffling elements within each row equiprobably
#' @export


sim2 <- function(speciesData) 

{
t(apply(speciesData,1,sample))
}


#' Sim3 Function
#' @description Randomizes a binary matrix speciesData by reshuffling elements within each column equiprobably
#' @export

sim3 <- function(speciesData) 

{
apply(speciesData,2,sample)
}



#' Sim4 Function
#' @description Randomizes a binary matrix speciesData by reshuffling elements within each row Sampling weights for each column are proportional to column sums Makes a call to the VectorSample
#' @export

sim4 <- function(speciesData) 

{
t(apply(speciesData,1,vector_sample,w=colSums(speciesData)))
}


#' Sim5 Function
#' @description Randomizes a matrix speciesData by reshuffling elements within each columnSampling weights for each row are proportional to row sums Makes a call to the VectorSample.
#' @export


sim5 <- function(speciesData) 

{
apply(speciesData,2,vector_sample,w=rowSums(speciesData))
}


#' Sim6 Function
#' Randomizes a binary matrix speciesData by reshuffling all elements Rows are equiprobable, columns proportional to column sums Makes a call to the VectorSample.
#' @export

sim6 <- function(speciesData) 

{
matrixWeights <- outer(rep(1,nrow(speciesData)),colSums(speciesData))
matrix(vector_sample(speciesData, w= matrixWeights),ncol=ncol(speciesData))
}


#' Sim7 Function
#' @description Randomizes a binary matrix speciesData by reshuffling all elements Columns are equiprobable, rows proportional to row sums Makes a call to the VectorSample
#' @export


sim7 <- function(speciesData) 

{
matrixWeights <- outer(rowSums(speciesData),rep(1,ncol(speciesData)))
matrix(vector_sample(speciesData, w=matrixWeights),ncol=ncol(speciesData))
}



#' Sim8 Function
#' @description Randomizes a binary matrix speciesData by reshuffling all elements Columns are proportional to column sums, rows proportional to row sums Makes a call to the VectorSample
#' @export


sim8 <- function(speciesData) 

{
matrixWeights <- outer(rowSums(speciesData),colSums(speciesData))
matrix(vector_sample(speciesData,w=matrixWeights),ncol=ncol(speciesData))
}



###############################################
######################################################################

#' SIM9 function Takes a binary presence absence matrix returns a new matrix with same number of rows and columns uses swapping function from Denitz
# Not used preserved for posterity
# sim9 <- function(speciesData=matrix(rbinom(100,1,0.5),nrow=10)) 
  
#{
#  speciesData <- speciesData[which(rowSums(speciesData)>0),] # make calculation on submatrix with no missing species
#  
#  Burn.In <- max(c(10*nrow(speciesData),1000)) # set the burn-in
#  # select two random rows and create submatrix
#  for(i in 1:Burn.In)
#  {
#    ran.rows <- sample(nrow(speciesData),2)
#    
#    speciesData.pair <- speciesData[ran.rows,]
#    
#    
#    
#    # find columns if any in pair for which colsum =1; these can be swapped
#    
#    Sum.Is.One <- which(colSums(speciesData.pair)==1)
#    
#    if(length(Sum.Is.One)>1)
#    {
#      # swap them in speciesData.pair.swapped
#      speciesData.pair.swapped <- speciesData.pair
#      speciesData.pair.swapped[,Sum.Is.One] <- speciesData.pair[,sample(Sum.Is.One)]
#      
#      # return the two swapped rows to the original speciesData matrix
#      
#      speciesData[ran.rows,] <- speciesData.pair.swapped
#      
#    }
#  }
#  return(speciesData)
#}


#' Sim10 Function
#' @description Randomizes a binary matrix speciesData by reshuffling all elements Rows & columns proportional to supplied row and column weights Makes a call to the VectorSample
#' @export


sim10 <- function(speciesData,rowWeights,colWeights) 

{
matrixWeights <- outer(rowWeights,colWeights)
matrix(vector_sample(speciesData, w=matrixWeights),ncol=ncol(speciesData))
}




#' RA1 Function
#' @description Randomizes a numeric utilization matrix speciesData by replacing all elements with a random uniform [0,1]
#' @export

ra1 <- function(speciesData=matrix(rpois(80,1),nrow=10)){

matrix(runif(prod(dim(speciesData))),ncol=ncol(speciesData))
}

#' RA2 Function
#' @description Randomizes a numeric utilization matrix speciesData by replacing non-zero elements with a random uniform [0,1]
#' @export

ra2 <- function(speciesData=matrix(rpois(80,1),nrow=10)) 

{
z <- which(speciesData > 0)
RM <- speciesData
RM[z] <- runif(length(z))
return(RM)
}

   

#' RA3 Function
#' @description Randomizes a numeric utilization matrix speciesData by reshuffling the elements within each row
#' @export

ra3 <- function(speciesData=matrix(rpois(80,1),nrow=10)) 

{
RM <- apply(speciesData,1,sample)
RM <- t(as.matrix(RM,nrow=nrow(speciesData),ncol=ncol(speciesData)))
return(RM)
}

   

#' RA4 Function
#' @description Randomizes a numeric utilization matrix speciesData by reshuffling the non-zero elements within each row
#' @export

ra4 <- function(speciesData=matrix(rpois(80,1),nrow=10)) 
  
{
    #.....................................
    ## Nonzero row shuffle function for RA4
    NonZeroRowShuffle <- function(vec=runif(10)) {
      nonzero <- which(vec > 0)
      shuffledvec <- vec
      shuffledvec[nonzero] <- vec[sample(nonzero)]
      return(shuffledvec)
    }
    #....................................
    
    RM <- t(apply(speciesData,1,NonZeroRowShuffle))
    return(RM)
    
  }


#'Uniform size algorithm
#'@description Function to randomize uniformly body sizes within  observed limits (classic Barton-David test)
#'@export
uniform_size <- function(speciesData=runif(20)) {
  
  endpoints <- c(min(speciesData),max(speciesData))  # save max and min boundaries
  sim <- runif(n=(length(speciesData)-2),min=min(speciesData),max=max(speciesData))
  randomVec <- c(endpoints,sim)
  return(randomVec)
}

#' User defined size limits
#' @description Function to randomize uniformly body sizes within user-defined limits. Note that all n species are randomized in this algorithm
#' @export
uniform_size_user <- function(speciesData=runif(n=20),userLow=0.9*min(speciesData),
                              userHigh=1.1*max(speciesData)){
#  if(!is.null(Param.List$Special)){User.low <- Param.List$Special[1]
#                                    User.high <- Param.List$Special[2]}
  randomVec <- runif(n=length(speciesData),min=userLow,max=userHigh)
  
  return(randomVec)
}

#' Source pool of body sizes
#' @description Function to randomize body sizes by drawing from a user-defined source pool. Species are drawn without replacement, and there is a specified probability vector for the source pool species
#' @export
source_pool_draw <- function(speciesData=21:30,sourcePool=
                               runif(n=2*length(speciesData),min=10,max=50),
                             speciesProbs=rep(1,length(sourcePool))) {
  
  #if(!is.null(Param.List$Special)){Source.Pool <- Param.List$Special[[1]]
  #                                 probs <- Param.List$Special[[2]]}
  
  speciesDraw <- sample(sourcePool,size=length(speciesData),replace=FALSE,
                         prob=speciesProbs)
  
  return(speciesDraw)
}

#' Gamma size
#' @description Function to estimate size distribution as a gamma function Gamma function parameters estimated by method of moments from unidentified pdf by Murphy 2007 http://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall07/reading/paramEst.pdf
#' @export
gamma_size <- function (speciesData=runif(20)) {
  a <- mean(speciesData)^2/var(speciesData) # use for shape parameter
  b <- mean(speciesData)/var(speciesData) # use for rate parameter
  gammaDraw <- rgamma(length(speciesData),shape=a,rate=b)
  return(gammaDraw)
}
