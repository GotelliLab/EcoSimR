empbayes<-function(commat, bytrt, nbrks=21, niter=1000, nburn=1000) {
  #Implements empirical Bayes method from Gotelli and Ulrich 2010 (Oecologia 162:463-77)
  #commat is community matrix with row=site, col=sp  
  #bytrt is list of separate treatments (corresponding to each site) -
  #A separate analysi will be run for each bytrt category
  #nbrks is number of breaks for histogram binning
  
  
  if(sum(dim(matrix(bytrt))>1)>1) {
    print("error - function not defined for bytrt of dimension >1")
    return()
  }
  
  require(EcoSimR)
  c_score_quantiles<-function(m = matrix(rbinom(100, 1, 0.5), nrow = 10),  bksq=seq(0,1,length=20)) {
    #From EcoSimR, Gotelli and Ellison 2013
    m <- m[which(rowSums(m) > 0), ]
    pairwise <- cbind(t(combn(nrow(m), 2)), 0)
    C.score <- mat.or.vec(nrow(pairwise), 1)
    shared <- mat.or.vec(nrow(pairwise), 1)
    
    maxC<-numeric(nrow(pairwise))
    for (i in 1:nrow(pairwise)) {
      shared[i] <- sum(m[pairwise[i, 1], ] == 1 & m[pairwise[i, 
                                                             2], ] == 1)
      C.score[i] <- (sum(m[pairwise[i, 1], ]) - shared[i]) * 
        (sum(m[pairwise[i, 2], ]) - shared[i])
      maxC[i]<-sum(m[pairwise[i, 1], ])*sum(m[pairwise[i, 2], ])
    }
    nsp<-nrow(m)
    C.score<-C.score/maxC
    cvecout<-table(cut(C.score, breaks=bksq))
    
    spnames<-cbind(rownames(m)[pairwise[,1]],rownames(m)[pairwise[,2]])
    
    return(list(cvecout=cvecout, C.score=C.score, spnames=spnames))
  }
  
  matbyList<-split.data.frame(commat, bytrt)
  bksq<-seq(-0.001, 1, length=nbrks) #Sequence for c-score breaks
  
  #Arrays for storing output
  qtout025<-matrix(nrow=length(matbyList), ncol=length(bksq)-1)
  qtout500<-matrix(nrow=length(matbyList), ncol=length(bksq)-1)
  qtout975<-matrix(nrow=length(matbyList), ncol=length(bksq)-1)
  qtout950<-matrix(nrow=length(matbyList), ncol=length(bksq)-1)
  
  ctable<-matrix(nrow=length(matbyList), ncol=length(bksq)-1)
  spparilst<-NULL
  cobslst<-NULL
  csimlst<-NULL
  
  empbayes_ci<-numeric(length(matbyList))
  empbayes_mu<-numeric(length(matbyList))
  c_ci<-numeric(length(matbyList))
  
  csimdi<-NULL
  cdisdi<-NULL
  
  for(i in 1:length(matbyList)) {
    print(paste(names(matbyList)[[i]], i/length(matbyList)))
    mtmp<-t(matbyList[[i]][,colSums(matbyList[[i]])>0])
    nsp<-ncol(mtmp)
    
    #Get C stat
    ctmp<-c_score_quantiles(mtmp, bksq)
    ctable[i,]<-ctmp$cvecout
    cobs<-ctmp$C.score
    
    csim<-matrix(nrow=niter, ncol=length(cobs))
    
    spparilst[[i]]<-ctmp$spnames
    
    #Get simulated null
    ctable_null<-matrix(ncol=length(ctable[i,]), nrow=niter)
    for(j in 1:niter) {
      m2<-sim2(mtmp)
      for(k in 1:nburn) {
        m2<-sim2(m2)
      }
      ctmp<-c_score_quantiles(m2, bksq)
      ctable_null[j,]<-ctmp$cvecout
      csim[j,]<-ctmp$C.score
      
      if(j/100==floor(j/100)) {
        print(j/niter)
      }
    }
    
    #Save outputs
    cobslst[[i]]<-cobs
    csimlst[[i]]<-csim
    
    #Calculate null quantiles for bins
    qt<-apply(ctable_null, 2, function(x) quantile(x, c(0.025, 0.5, 0.975, 0.95)))
    qtout025[i,]<-qt[1,]
    qtout500[i,]<-qt[2,]
    qtout975[i,]<-qt[3,]
    qtout950[i,]<-qt[4,]
    
    plot(c(0, 1), c(0, max(c(qt, ctable[i,]))),
         xlab="C-Stat", ylab="freq", type="n", main=names(matbyList)[i])
    xpos<-bksq[-length(bksq)]+diff(bksq)/2
    
    polygon(c(xpos, rev(xpos)), c(qtout025[i,], rev(qtout975[i,])), border=NA, col="grey")
    lines(xpos, qtout500[i,], lwd=2, col="darkgrey")
    lines(xpos, ctable[i,], lwd=2) 
    
    #Calculate emp bayes
    empbayes_ci[i]<-sum(ctable[i,]>=qt[3,])
    empbayes_mu[i]<-sum(ctable[i,]>=qt[2,])
    
    #Calculate C ci
    cqt<-apply(csim, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
    #too similar
    csimdi[[i]]<-list(which(cobs<=cqt[1,]))
    csimdi[[i]]<-csimdi[[i]][[1]]
    #not similar enough
    cdisdi[[i]]<-list(which(cobs>=cqt[3,]))
    cdisdi[[i]]<-cdisdi[[i]][[1]]
    
    c_ci[i]<-length(c(csimdi[[i]], cdisdi[[i]]))
  }
  
  datout<-list(cobslst=cobslst, csimlst=csimlst, 
               ctable=ctable, qtout025=qtout025, qtout500=qtout500, qtout975=qtout975, qtout950=qtout950,
               empbayes_ci=empbayes_ci, empbayes_mu=empbayes_mu,
               csimdi=csimdi, cdisdi=cdisdi, c_ci=c_ci,
               spparilst=spparilst, datnames=names(matbyList),
               bksq=bksq)
  
  #cobslst and csimlst are observed and simulated c-scores for species at each site
  #catable is expected number of c-scores in each bin. qtout__ is the __ quantile for expected number of entries in each bin
  #embbayes_ci and empbayes_mu are number of species in each bin that exceed the ci and mean criteria for the test
  #csimdi and cdisdi show species pairs that have cooccurrences that are more similar or more dissimlar than expected by chance
  #c_ci is total number of species pairs in each bin that co-occur more or less frequently than expected by chance
  #spparilst is a lookup table for numerical values explaining which species pair each corresponds to
  #datnames is names of the treatment categories
  #bksq is vector of breaks used for the binning in the analysis
  
  return(datout)
}