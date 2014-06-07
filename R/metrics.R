##############################################
## EcoSimR: R code for Null Model Analysis
## Metrics for 
## Nicholas J. Gotelli & Aaron M. Ellison
##
##
##############################################
## Version 1.00
## 12 May 2013
#############################################


#' Pianka function
#' @description Takes a niche utilization matrix returns Pianka's niche overlap index
#' @export
#'

Pianka <- function(m=matrix(rpois(80,1),nrow=10)) 

{
m <- m/rowSums(m)
pairwise <- cbind(t(combn(nrow(m),2)),0)	# set up pairwise species list
for (i in 1:nrow(pairwise)) 
	pairwise[i,3] <- sum(m[pairwise[i,1],]*m[pairwise[i,2],])/
	sqrt(sum(m[pairwise[i,1],]^2)*sum(m[pairwise[i,2],]^2))
return(mean(pairwise[,3]))
}


#' Czekanowski function
#' @description Takes a niche utilization matrix returns Czekanowski niche overlap index
#'@export
#'

Czekanowski <- function(m=matrix(rpois(80,1),nrow=10)) 

{
m <- m/rowSums(m)
pairwise <- cbind(t(combn(nrow(m),2)),0)	# set up pairwise species list
for (i in 1:nrow(pairwise)) 
	pairwise[i,3] <- 1 - 0.5*sum(abs((m[pairwise[i,1],] - m[pairwise[i,2],])))

return(mean(pairwise[,3]))
}


#' Pianka variance function
#' @description Takes a niche utilization matrix returns variance of Pianka's niche overlap index
#' @export

Pianka.var <- function(m=matrix(rpois(80,1),nrow=10)) 
  
{
  m <- m/rowSums(m)
  pairwise <- cbind(t(combn(nrow(m),2)),0)  # set up pairwise species list
  for (i in 1:nrow(pairwise)) 
    pairwise[i,3] <- sum(m[pairwise[i,1],]*m[pairwise[i,2],])/
      sqrt(sum(m[pairwise[i,1],]^2)*sum(m[pairwise[i,2],]^2))
  
  return(var(pairwise[,3]))
}

##
##
#' Czekanowski variance function
#' @description Takes a niche utilization matrix returns variance Czekanowski niche overlap index
#' @export

Czekanowski.var <- function(m=matrix(rpois(80,1),nrow=10)) 
  
{
  m <- m/rowSums(m)
  pairwise <- cbind(t(combn(nrow(m),2)),0)	# set up pairwise species list
  for (i in 1:nrow(pairwise)) 
    pairwise[i,3] <- 1 - 0.5*sum(abs((m[pairwise[i,1],] - m[pairwise[i,2],])))
  
  return(var(pairwise[,3]))
}


#' Pianka skewness function
#' @description Takes a niche utilization matrix returns skewness of Pianka's niche overlap index
#' @export

Pianka.skew <- function(m=matrix(rpois(80,1),nrow=10)) 
  
{
  m <- m/rowSums(m)
  pairwise <- cbind(t(combn(nrow(m),2)),0)  # set up pairwise species list
  for (i in 1:nrow(pairwise)) 
    pairwise[i,3] <- sum(m[pairwise[i,1],]*m[pairwise[i,2],])/
    sqrt(sum(m[pairwise[i,1],]^2)*sum(m[pairwise[i,2],]^2))
  
    m3 <- mean((pairwise[,3]-mean(pairwise[,3]))^3)
    Pianka.skew <- m3/(sd(pairwise[,3])^3)
    
  
  return(Pianka.skew)
}


#' Czekanowski skew function
#' @description Takes a niche utilization matrix returns skew of the Czekanowski niche overlap index
#' @export

Czekanowski.skew <- function(m=matrix(rpois(80,1),nrow=10)) 
  
{
  m <- m/rowSums(m)
  pairwise <- cbind(t(combn(nrow(m),2)),0)  # set up pairwise species list
  for (i in 1:nrow(pairwise)) 
    pairwise[i,3] <- 1 - 0.5*sum(abs((m[pairwise[i,1],] - m[pairwise[i,2],])))
  
  m3 <- mean((pairwise[,3]-mean(pairwise[,3]))^3)
  Czekanowski.skew <- m3/(sd(pairwise[,3])^3)
  
  
  return(Czekanowski.skew)
}

##
##


# Co-occurrence algorithms and metrics for testing
# 25 May 2013
# NJG

# ###############################################
#' @description Function to Calculate Number Of Species Combinations
#' @export
#' 
Species.Combo <- function(m=matrix(rbinom(100,1,0.5),nrow=10))
{
  
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  return(ncol(unique(m,MARGIN=2))) # number of columns in submatrix of uniques
  
}
#############################################
#' @description Checkerboard function function Takes a binary presence absence matrix returns number of checkerboard pairs
#' @export


Checker <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  pairwise <- cbind(t(combn(nrow(m),2)),0) # set up pairwise species list
  
  
  shared <- mat.or.vec(1,nrow(pairwise))
  
  for (i in 1:nrow(pairwise)) 
  {
    shared[i] <- sum(m[pairwise[i,1],]==1 & m[pairwise[i,2],]==1)
  }
  
  
  
  return(sum(shared==0)) # return number of pairs with no shared sites
  
}

#' @description C-score function Takes a binary presence absence matrix returns Stone and Roberts C-score
#' @export

C.Score <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  pairwise <- cbind(t(combn(nrow(m),2)),0) # set up pairwise species list
  
  
  C.score <- mat.or.vec(nrow(pairwise),1)
  shared <- mat.or.vec(nrow(pairwise),1)
  
  for (i in 1:nrow(pairwise)) 
  {
    shared[i] <- sum(m[pairwise[i,1],]==1 & m[pairwise[i,2],]==1)
    C.score[i] <- (sum(m[pairwise[i,1],]) - shared[i])*
      (sum(m[pairwise[i,2],]) - shared[i])
    
    
  }
  
  return(mean(C.score)) # return average C-score
  
}


#' @description C-score variance function Takes a binary presence absence matrix returns variance of Stone and Roberts C-score
#' @export


C.Score.var <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  pairwise <- cbind(t(combn(nrow(m),2)),0) # set up pairwise species list
  
  
  C.score <- mat.or.vec(nrow(pairwise),1)
  shared <- mat.or.vec(nrow(pairwise),1)
  
  for (i in 1:nrow(pairwise)) 
  {
    shared[i] <- sum(m[pairwise[i,1],]==1 & m[pairwise[i,2],]==1)
    C.score[i] <- (sum(m[pairwise[i,1],]) - shared[i])*
      (sum(m[pairwise[i,2],]) - shared[i])
    
    
  }
  
  
  
  return(var(C.score))  # return variance of pairwise C-score
  
}

#' @description C-score skew function Takes a binary presence absence matrix returns skewness of Stone and Roberts C-score
#' @export
#' 
C.Score.skew <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  pairwise <- cbind(t(combn(nrow(m),2)),0) # set up pairwise species list
  
  
  C.score <- mat.or.vec(nrow(pairwise),1)
  shared <- mat.or.vec(nrow(pairwise),1)
  
  for (i in 1:nrow(pairwise)) 
  {
    shared[i] <- sum(m[pairwise[i,1],]==1 & m[pairwise[i,2],]==1)
    C.score[i] <- (sum(m[pairwise[i,1],]) - shared[i])*
      (sum(m[pairwise[i,2],]) - shared[i])
    
    
  }
  
  m3 <- mean((C.score-mean(C.score))^3)
  C.score.skew <- m3/(sd(C.score)^3)
  
  
  return(C.score.skew)  # return skewness of pairwise C-score
  
}
######################################################################
## Schluter's V-ratio function
## Takes a binary presence absence matrix
## returns Schluter's V-ratio index
##
##
V.Ratio <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  
  v <- var(colSums(m))/sum(apply(m,1,var))
  return(v)
}

