Hill_Shannon <- function(m=matrix(rpois(80,1),nrow=40)) 
  
{
  # Function calculates beta diversity based on
  # second order Hill number (Shannon) decomposition
  # sensu Jost 2007 Ecology 88(10) and Death 2012 Ecology 93(10)
  
  # alpha diversity
  pa <- t(t(m)/colSums(m))
  Ha <- -pa*log(pa)
  Ha[!is.finite(Ha)] <- 0
  Ha <- sum(Ha)/nrow(m)
  
  # gamma diversity
  pg <- rowMeans(m)/sum(m)
  Hg <- -pg*log(pg)
  Hg[!is.finite(Hg)] <- 0
  Hg <- sum(Hg)
  
  # beta diversity (as D, not as H - i.e., is transformed to diversity equivalent)
  eHb <- exp(Hg)/exp(Ha)
  
  return(eHb)
}
