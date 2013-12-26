Dist_test <- function(m=matrix(rpois(nrow(distances)*5,1),ncol=nrow(distances))) 
  
{
  #calculate pairwise Bray dissimilarity
  sp_dist<-as.matrix(vegdist(t(m)))
  
  #return slope of species vs. distance dissimilarity
  dist_slope<-lm(c(distances)~c(sp_dist))$coefficients[2]
  names(dist_slope)<-"dist_slope"
  
  return(dist_slope)
}
