##############################################
## EcoSimR: R code for Null Model Analysis
##
##
##
##############################################
## EcoSimR Metrics Source File
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

