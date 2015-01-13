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

vector_sample <- function(v,w) 

{
x <- mat.or.vec(length(v),1)                   # creates a vector of 0s
x[sample(1:length(v),size=sum(v),prob=w)] <- 1  # fills with 1s, sampling with weights
return(x)
}


#' Sim1 Function
#' @description Randomizes a binary matrix m by reshuffling all of its elements equiprobably
#' @export

sim1 <- function(m) 

{
matrix(sample(m), ncol=ncol(m))
}



#' Sim2 Function
#' @description Randomizes a binary matrix m by reshuffling elements within each row equiprobably
#' @export


sim2 <- function(m) 

{
t(apply(m,1,sample))
}


#' Sim3 Function
#' @description Randomizes a binary matrix m by reshuffling elements within each column equiprobably
#' @export

sim3 <- function(m) 

{
apply(m,2,sample)
}



#' Sim4 Function
#' @description Randomizes a binary matrix m by reshuffling elements within each row Sampling weights for each column are proportional to column sums Makes a call to the VectorSample
#' @export

sim4 <- function(m) 

{
t(apply(m,1,vector_sample,w=colSums(m)))
}


#' Sim5 Function
#' @description Randomizes a matrix m by reshuffling elements within each columnSampling weights for each row are proportional to row sums Makes a call to the VectorSample.
#' @export


sim5 <- function(m) 

{
apply(m,2,vector_sample,w=rowSums(m))
}


#' Sim6 Function
#' Randomizes a binary matrix m by reshuffling all elements Rows are equiprobable, columns proportional to column sums Makes a call to the VectorSample.
#' @export

sim6 <- function(m) 

{
Matrix.Weights <- outer(rep(1,nrow(m)),colSums(m))
matrix(vector_sample(m, w=Matrix.Weights),ncol=ncol(m))
}


#' Sim7 Function
#' @description Randomizes a binary matrix m by reshuffling all elements Columns are equiprobable, rows proportional to row sums Makes a call to the VectorSample
#' @export


sim7 <- function(m) 

{
Matrix.Weights <- outer(rowSums(m),rep(1,ncol(m)))
matrix(vector_sample(m, w=Matrix.Weights),ncol=ncol(m))
}



#' Sim8 Function
#' @description Randomizes a binary matrix m by reshuffling all elements Columns are proportional to column sums, rows proportional to row sums Makes a call to the VectorSample
#' @export


sim8 <- function(m) 

{
Matrix.Weights <- outer(rowSums(m),colSums(m))
matrix(vector_sample(m,w=Matrix.Weights),ncol=ncol(m))
}



###############################################
######################################################################

#' SIM9 function Takes a binary presence absence matrix returns a new matrix with same number of rows and columns uses swapping function from Denitz

sim9 <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  Burn.In <- max(c(10*nrow(m),1000)) # set the burn-in
  # select two random rows and create submatrix
  for(i in 1:Burn.In)
  {
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
  }
  return(m)
}


#' Sim10 Function
#' @description Randomizes a binary matrix m by reshuffling all elements Rows & columns proportional to supplied row and column weights Makes a call to the VectorSample
#' @export


sim10 <- function(m,Row.Weights,Col.Weights) 

{
Matrix.Weights <- outer(Row.Weights,Col.Weights)
matrix(vector_sample(m, w=Matrix.Weights),ncol=ncol(m))
}




#' RA1 Function
#' @description Randomizes a numeric utilization matrix m by replacing all elements with a random uniform [0,1]
#' @export

ra1 <- function(m=matrix(rpois(80,1),nrow=10)){

matrix(runif(prod(dim(m))),ncol=ncol(m))
}

#' RA2 Function
#' @description Randomizes a numeric utilization matrix m by replacing non-zero elements with a random uniform [0,1]
#' @export

ra2 <- function(m=matrix(rpois(80,1),nrow=10)) 

{
z <- which(m > 0)
RM <- m
RM[z] <- runif(length(z))
return(RM)
}

   

#' RA3 Function
#' @description Randomizes a numeric utilization matrix m by reshuffling the elements within each row
#' @export

ra3 <- function(m=matrix(rpois(80,1),nrow=10)) 

{
RM <- apply(m,1,sample)
RM <- t(as.matrix(RM,nrow=nrow(m),ncol=ncol(m)))
return(RM)
}

   

#' RA4 Function
#' @description Randomizes a numeric utilization matrix m by reshuffling the non-zero elements within each row
#' @export

ra4 <- function(m=matrix(rpois(80,1),nrow=10)) 
  
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
    
    RM <- t(apply(m,1,NonZeroRowShuffle))
    return(RM)
    
  }


#'Uniform size algorithm
#'@description Function to randomize uniformly body sizes within  observed limits (classic Barton-David test)
#'@export
uniform_size <- function(v=runif(20)) {
  
  Endpoints <- c(min(v),max(v))  # save max and min boundaries
  Sim <- runif(n=(length(v)-2),min=min(v),max=max(v))
  Random.Vec <- c(Endpoints,Sim)
  return(Random.Vec)
}

#' User defined size limits
#' @description Function to randomize uniformly body sizes within user-defined limits. Note that all n species are randomized in this algorithm
#' @export
uniform_size_user <- function(v=runif(n=20),User.low=0.9*min(v),
                              User.high=1.1*max(v)){
  if(!is.null(Param.List$Special)){User.low <- Param.List$Special[1]
                                   User.high <- Param.List$Special[2]}
  Random.Vec <- runif(n=length(v),min=User.low,max=User.high)
  
  return(Random.Vec)
}

#' Source pool of body sizes
#' @description Function to randomize body sizes by drawing from a user-defined source pool. Species are drawn without replacement, and there is a specified probability vector for the source pool species
#' @export
source_pool_draw <- function(v=21:30,Source.Pool=
                               runif(n=2*length(v),min=10,max=50),
                             Species.Probs=rep(1,length(Source.Pool))) {
  if(!is.null(Param.List$Special)){Source.Pool <- Param.List$Special[[1]]
                                   probs <- Param.List$Special[[2]]}
  
  Species.Draw <- sample(Source.Pool,size=length(Data.Matrix),replace=FALSE,
                         prob=Species.Probs)
  
  return(Species.Draw)
}

#' Gamma size
#' @description Function to estimate size distribution as a gamma function Gamma function parameters estimated by method of moments from unidentified pdf by Murphy 2007 http://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall07/reading/paramEst.pdf
#' @export
gamma_size <- function (v=runif(20)) {
  a <- mean(v)^2/var(v) # use for shape parameter
  b <- mean(v)/var(v) # use for rate parameter
  Gamma.Draw <- rgamma(length(v),shape=a,rate=b)
  return(Gamma.Draw)
}
