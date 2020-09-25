#' Unbiased variance of a single ROC curve.
#'
#' Calculate unbiased variance of a single ROC curve.
#' @usage uvar(ROC)
#' @param ROC An \code{ROC} object obtained by ROC.
#' @return A numeric value.
#' @references Lu, Y. and Shao, Y. (2020). Preprint
#' @section Warning:
#' Variance of AUC=1 is always 0 and can be misleading.
#' @examples
#' set.seed(123)
#' ## Generate data:
#' Data=SimuCaseCont(100,100,0.8,0.9,0.9,0.9,"c","t",marg1="Normal",marg2="Normal")
#' ## Construct ROC curves
#' roc1=ROC(Data[,3],Data[,1],case_score_ind = "higher")
#' roc2=ROC(Data[,3],Data[,2],case_score_ind = "higher",ci=TRUE,level=0.95)
#' ## Calculate the unbiased variance of two ROC curves
#' uvar(roc1)
#' uvar(roc2)
#' @export

uvar<-function(ROC) {

  K = ROC$K
  m=nrow(K) # number of case
  n=ncol(K) # number of control
  AUC = mean(K)

  s1=sum(K^2)
  m1=s1/(m*n)
  s2=sum(colSums(K)^2)
  m2=(s2-s1)/(m*n*(m-1))
  s3=sum(rowSums(K)^2)
  m3=(s3-s1)/(m*n*(n-1))
  s4=sum(K)^2
  m4=(s4-s3-s2+s1)/(m*n*(m-1)*(n-1))
  c1=1/(m*n)
  c2=(m-1)/(m*n)
  c3=(n-1)/(m*n)
  c4=(m-1)*(n-1)/(m*n)
  variance=c1*m1+c2*m2+c3*m3-(1-c4)*m4
  variance=max(variance,1e-9)
  if (variance==1e-9 & ROC$auc==1) {
    variance=0
    warning("Variance of a ROC curve with AUC == 1 is always 0 and can be misleading.")
  }
  return(variance)
}
