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


#' Pianka niche overlap
#' @description Takes a resource utilization matrix as input and
#' returns the average pairwise Pianka's niche overlap index.
#' @details Pianka's niche overlap index is averaged 
#' over each unique species pair. The index is symmetric, 
#' with a normalization term in the denominator for the overlap 
#' between species 1 and 2. Values of Pianka's niche overlap index close to 0.0
#' reflect usage of exclusive resource categories, whereas values close to 1.0
#' reflect similar resource utilization spectra.
#' \deqn{O_{jk} = O_{kj} = \frac{\sum_{n}^{i}p_{ij}p_{jk}}{\sqrt{\sum_{n}^{i}p_{ij}^2\sum_{n}^{i}p_{ik}^2}}}{O_jk = O_kj = sum(p_ij*p_jk) / sqrt(sum((p_ij)^2)sum((p_jk)^2))}
#' 
#' @param m A matrix of resource utilization values. 
#' @return Returns the average pairwise niche overlap.
#' @references Pianka, E. 1973. The structure of lizard communities.
#' Annual Review of Ecology and Systematics 4:53-74.
#' 
#' Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural assemblages
#' of desert lizards and tropical fishes. Ecological Monographs 60: 27-55.
#' 
#' @note The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @seealso \code{\link{czekanowski}} niche overlap index.
#' @examples 
#' obsOverlap <- pianka(m=matrix(rpois(40,0.5),nrow=8))
#' @export
#'

pianka <- function(m=matrix(rpois(80,1),nrow=10)) 

{
m <- m/rowSums(m)
pairwise <- cbind(t(combn(nrow(m),2)),0)	# set up pairwise species list
for (i in 1:nrow(pairwise)) 
	pairwise[i,3] <- sum(m[pairwise[i,1],]*m[pairwise[i,2],])/
	sqrt(sum(m[pairwise[i,1],]^2)*sum(m[pairwise[i,2],]^2))
return(mean(pairwise[,3]))
}


#' Czekanowski niche overlap
#' @description Takes a resource utilization matrix as input and
#' returns the average pairwise Czekanowki niche overlap index.
#' @details The Czekanowski niche overlap index is averaged 
#' over each unique species pair. The index measures the area of intersection
#' of the resource utilization histograms of each species pair. 
#' Values of Czekanowski niche overlap index close to 0.0
#' reflect usage of exclusive resource categories, whereas values close to 1.0
#' reflect similar resource utilization spectra.
#' 
#' @param m A matrix of resource utilization values. 
#' @return Returns the average pairwise niche overlap.
#' @references Feinsinger, P., E.E. Spears, and R. Poole. 1981.
#' A simple measure of niche breadth. Ecology 62: 27-32. 
#' 
#' Winemiller, K.O. and E.R. Pianka. 1990. Organization in natural assemblages
#' of desert lizards and tropical fishes. Ecological Monographs 60: 27-55.
#' 
#' @note The resource utilization matrix (rows = species, columns = discrete
#' resource categories) may include zeroes, but no negative numbers or missing
#' values. Relative resource within a species is first calculated, so the rows
#' need not sum to 1.0.
#' @seealso \code{\link{pianka}} niche overlap index.
#' @examples 
#' obsOverlap <- pianka(m=matrix(rpois(40,0.5),nrow=8))
#' @export
#'

czekanowski <- function(m=matrix(rpois(80,1),nrow=10)) 

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

pianka_var <- function(m=matrix(rpois(80,1),nrow=10)) 
  
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

czekanowski_var <- function(m=matrix(rpois(80,1),nrow=10)) 
  
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

pianka_skew <- function(m=matrix(rpois(80,1),nrow=10)) 
  
{
  m <- m/rowSums(m)
  pairwise <- cbind(t(combn(nrow(m),2)),0)  # set up pairwise species list
  for (i in 1:nrow(pairwise)) 
    pairwise[i,3] <- sum(m[pairwise[i,1],]*m[pairwise[i,2],])/
    sqrt(sum(m[pairwise[i,1],]^2)*sum(m[pairwise[i,2],]^2))
  
    m3 <- mean((pairwise[,3]-mean(pairwise[,3]))^3)
    piankaSkew <- m3/(sd(pairwise[,3])^3)
    
  
  return(piankaSkew)
}


#' Czekanowski skew function
#' @description Takes a niche utilization matrix returns skew of the Czekanowski niche overlap index
#' @export
czekanowski_skew <- function(m=matrix(rpois(80,1),nrow=10)) 
  
{
  m <- m/rowSums(m)
  pairwise <- cbind(t(combn(nrow(m),2)),0)  # set up pairwise species list
  for (i in 1:nrow(pairwise)) 
    pairwise[i,3] <- 1 - 0.5*sum(abs((m[pairwise[i,1],] - m[pairwise[i,2],])))
  
  m3 <- mean((pairwise[,3]-mean(pairwise[,3]))^3)
  czekanowskiSkew <- m3/(sd(pairwise[,3])^3)
  
  
  return(czekanowskiSkew)
}

##
##


# Co-occurrence algorithms and metrics for testing
# 25 May 2013
# NJG

#' Species combinations
#' @description Function to Calculate Number Of Species Combinations in a matrix
#' @export
species_combo <- function(m=matrix(rbinom(100,1,0.5),nrow=10))
{
  
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  return(ncol(unique(m,MARGIN=2))) # number of columns in submatrix of uniques
  
}
#' Checkerboard function
#' @description Takes a binary presence absence matrix returns number of checkerboard pairs
#' @export
checker <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
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
#' calculate c-score
#' @description C-score function Takes a binary presence absence matrix returns Stone and Roberts C-score
#' @export

c_score <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  pairwise <- cbind(t(combn(nrow(m),2)),0) # set up pairwise species list
  
  
  cScore <- mat.or.vec(nrow(pairwise),1)
  shared <- mat.or.vec(nrow(pairwise),1)
  
  for (i in 1:nrow(pairwise)) 
  {
    shared[i] <- sum(m[pairwise[i,1],]==1 & m[pairwise[i,2],]==1)
    cScore[i] <- (sum(m[pairwise[i,1],]) - shared[i])*
      (sum(m[pairwise[i,2],]) - shared[i])
    
    
  }
  
  return(mean(cScore)) # return average C-score
  
}

#' Calculate variance of the c-score
#' @description C-score variance function Takes a binary presence absence matrix returns variance of Stone and Roberts C-score
#' @export


c_score_var <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  pairwise <- cbind(t(combn(nrow(m),2)),0) # set up pairwise species list
  
  
  cScore <- mat.or.vec(nrow(pairwise),1)
  shared <- mat.or.vec(nrow(pairwise),1)
  
  for (i in 1:nrow(pairwise)) 
  {
    shared[i] <- sum(m[pairwise[i,1],]==1 & m[pairwise[i,2],]==1)
    cScore[i] <- (sum(m[pairwise[i,1],]) - shared[i])*
      (sum(m[pairwise[i,2],]) - shared[i])
    
    
  }
  
  
  
  return(var(cScore))  # return variance of pairwise C-score
  
}

#' Calculate the skew of the c-score
#' @description C-score skew function Takes a binary presence absence matrix returns skewness of Stone and Roberts C-score
#' @export
#' 
c_score_skew <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
  
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  pairwise <- cbind(t(combn(nrow(m),2)),0) # set up pairwise species list
  
  
  cScore <- mat.or.vec(nrow(pairwise),1)
  shared <- mat.or.vec(nrow(pairwise),1)
  
  for (i in 1:nrow(pairwise)) 
  {
    shared[i] <- sum(m[pairwise[i,1],]==1 & m[pairwise[i,2],]==1)
    cScore[i] <- (sum(m[pairwise[i,1],]) - shared[i])*
      (sum(m[pairwise[i,2],]) - shared[i])
    
    
  }
  
  m3 <- mean((cScore-mean(cScore))^3)
  cScoreSkew <- m3/(sd(cScore)^3)
  
  
  return(cScoreSkew)  # return skewness of pairwise C-score
  
}

#' Schluter's V-ratio function
#' @description Takes a binary presence absence matrix returns Schluter's V-ratio index
#' @export
v_ratio <- function(m=matrix(rbinom(100,1,0.5),nrow=10)) 
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  
  
  v <- var(colSums(m))/sum(apply(m,1,var))
  return(v)
}


#' Min size difference
#' @description Function to calculate the minimum absolute size difference
#' @export
#' 
min_diff <- function(m=runif(20)) {
  m <- sort(m)
  md <- min(diff(m))
  return(md)
}

#' Min size ratio
#' @description Function to calculate the minimum size ratio
#' @export

min_ratio <- function(m=runif(20)) {
  m <- sort(log(m))
  mr <- min(diff(m))
  return(mr)
}

#' Variance in size
#' @description Function to calculate the variance in size differences
#' @export

var_diff <- function(m=runif(20)){
  m <- sort(m)
  vd <- var(diff(m))
  return(vd)
}

#' Variance in ratio
#' @description  Function to calculate the variance in size ratios
#' @export

var_ratio <- function(m=runif(20)){
  m <- sort(log(m))
  vr <- var(diff(m))
  return(vr)
}
