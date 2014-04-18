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

VectorSample <- function(v,w) 

{
x <- mat.or.vec(length(v),1)                   # creates a vector of 0s
x[sample(1:length(v),size=sum(v),prob=w)] <- 1  # fills with 1s, sampling with weights
return(x)
}


#' Sim1 Function
#' @description Randomizes a binary matrix m by reshuffling all of its elements equiprobably
#' @export


Sim1 <- function(m) 

{
matrix(sample(m), ncol=ncol(m))
}



#' Sim2 Function
#' @description Randomizes a binary matrix m by reshuffling elements within each row equiprobably
#' @export


Sim2 <- function(m) 

{
t(apply(m,1,sample))
}


#' Sim3 Function
#' @description Randomizes a binary matrix m by reshuffling elements within each column equiprobably
#' @export

Sim3 <- function(m) 

{
apply(m,2,sample)
}



#' Sim4 Function
#' @description Randomizes a binary matrix m by reshuffling elements within each row Sampling weights for each column are proportional to column sums Makes a call to the VectorSample
#' @export



Sim4 <- function(m) 

{
t(apply(m,1,VectorSample,w=colSums(m)))
}


#' Sim5 Function
#' @description Randomizes a matrix m by reshuffling elements within each columnSampling weights for each row are proportional to row sums Makes a call to the VectorSample.
#' @export


Sim5 <- function(m) 

{
apply(m,2,VectorSample,w=rowSums(m))
}


#' Sim6 Function
#' Randomizes a binary matrix m by reshuffling all elements Rows are equiprobable, columns proportional to column sums Makes a call to the VectorSample.
#' @export

Sim6 <- function(m) 

{
Matrix.Weights <- outer(rep(1,nrow(m)),colSums(m))
matrix(VectorSample(m, w=Matrix.Weights),ncol=ncol(m))
}


#' Sim7 Function
#' @description Randomizes a binary matrix m by reshuffling all elements Columns are equiprobable, rows proportional to row sums Makes a call to the VectorSample
#' @export


Sim7 <- function(m) 

{
Matrix.Weights <- outer(rowSums(m),rep(1,ncol(m)))
matrix(VectorSample(m, w=Matrix.Weights),ncol=ncol(m))
}



#' Sim8 Function
#' @description Randomizes a binary matrix m by reshuffling all elements Columns are proportional to column sums, rows proportional to row sums Makes a call to the VectorSample
#' @export


Sim8 <- function(m) 

{
Matrix.Weights <- outer(rowSums(m),colSums(m))
matrix(VectorSample(m,w=Matrix.Weights),ncol=ncol(m))
}


#' Sim10 Function
#' @description Randomizes a binary matrix m by reshuffling all elements Rows & columns proportional to supplied row and column weights Makes a call to the VectorSample
#' @export


Sim10 <- function(m,Row.Weights,Col.Weights) 

{
Matrix.Weights <- outer(Row.Weights,Col.Weights)
matrix(VectorSample(m, w=Matrix.Weights),ncol=ncol(m))
}


#' RA1 Function
#' @description Randomizes a numeric utilization matrix m by replacing all elements with a random uniform [0,1]
#' @export

RA1 <- function(m=matrix(rpois(80,1),nrow=10)){

matrix(runif(prod(dim(m2))),ncol=ncol(m2))
}

#' RA2 Function
#' @description Randomizes a numeric utilization matrix m by replacing non-zero elements with a random uniform [0,1]
#' @export

RA2 <- function(m=matrix(rpois(80,1),nrow=10)) 

{
z <- which(m > 0)
RM <- m
RM[z] <- runif(length(z))
return(RM)
}

   

#' RA3 Function
#' @description Randomizes a numeric utilization matrix m by reshuffling the elements within each row
#' @export

RA3 <- function(m=matrix(rpois(80,1),nrow=10)) 

{
RM <- apply(m,1,sample)
RM <- t(as.matrix(RM,nrow=nrow(m),ncol=ncol(m)))
return(RM)
}

   

#' RA4 Function
#' @description Randomizes a numeric utilization matrix m by reshuffling the non-zero elements within each row
#' @export

RA4 <- function(m=matrix(rpois(80,1),nrow=10)) 
  
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
#############################################


##
## 
#############################################