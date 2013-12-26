Sorensen <- function(m=matrix(rpois(80,1),nrow=40)) 
  
{
  # Function calculates Sorensen similarity index
  # 0 is no species shared, 1 is all species shared
  # NOTE - works only for pairs of sites

  if(ncol(m)!=2) {
    print("Error - can only compare pairs of sites. Repeat with a matrix with only two columns.")
  } else {
  
    #shared species
    C <- sum(m[,1]*m[,2]>0)
    
    #species in each plot
    A <- sum(m[,1]>0)
    B <- sum(m[,2]>0)
    
    #calculate index
    sim <- (2*C)/(A+B)
    
    return(sim)
  }
}
