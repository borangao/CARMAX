#' @keywords internal
#' Utility function for processing z and ld matrices
#'
#' This function takes summary statistics and LD matrices as input and computes 
#' a joint LD matrix used for the CARMA-X fine-mapping model.
#'
#' @param z.list A list of z-score summary statistics for different ancestries.
#' @param ld.list A list of LD matrices for different ancestries.
#' @return A combined LD matrix for all ancestries.
utility<-function(z.list,ld.list){

  p<-nrow(z.list[[1]])
  
  z<-as.matrix(c(z.list[[1]],z.list[[2]]))
  
  S0<-matrix(0,2*p,2*p)
  
  S0[1:p,1:p]<-ld.list[[1]]
  S0[1:p+p,1:p+p]<-ld.list[[2]]
  nu0=2*p+2
  S<-S0+z%*%t(z)/max(z^2)
  S<-S+diag(0.001,nrow(S))
  
  ES<-S/(nu0+1-1-2*p)
  DS<-diag(ES)
  ds_inv<-diag(1/sqrt(DS))
  ESR<-t(ds_inv)%*%ES%*%ds_inv
  R<-ESR
  diag(R)<-1
  s.off<-R[1:p,1:p+p]
  S0[1:p,1:p+p]<-s.off
  S0[1:p+p,1:p]<-t(s.off)
  ld.all<-S0
  return(R)
}
