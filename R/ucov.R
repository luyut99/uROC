#' Unbiased covariance of AUCs under two ROC curves
#'
#' Calculate unbiased covariance of AUCs under two ROC curves.
#' @usage ucov(ROC1, ROC2)
#' @param ROC1,ROC2 Two \code{ROC} objects referring two ROC curves.
#' @return A numeric value.
#' @references Lu, Y. and Shao, Y. (2020). ucompROC: A new powerful test to compare correlated ROC curves.
#' @section Warning:
#' Covariance of two AUCs with any AUC=1 is always 0 and can be misleading.
#' @examples
#' library(uROC)
#'
#' set.seed(123)
#' ## Generate data:
#' Data=SimuCaseCont(100,100,0.8,0.9,0.9,0.9,"c","t",marg1="Normal",marg2="Normal")
#' ## Construct ROC curves
#' roc1=ROC(Data[,3],Data[,1],case_score_ind = "higher")
#' roc2=ROC(Data[,3],Data[,2],case_score_ind = "higher",ci=TRUE,level=0.95)
#' ## Compute the covariance of two ROC curves
#' ucov(roc1,roc2)
#' @export

ucov<-function(ROC1,ROC2) {

  K1 = ROC1$K
  K2 = ROC2$K
  m=nrow(K1)
  n=ncol(K2)
  s1=sum(K1*K2)
  m1=s1/(m*n)
  s2=sum(colSums(K1)*colSums(K2))
  m2=(s2-s1)/(m*n*(m-1))
  s3=sum(rowSums(K1)*rowSums(K2))
  m3=(s3-s1)/(m*n*(n-1))
  s4=sum(K1)*sum(K2)
  m4=(s4-s3-s2+s1)/(m*n*(m-1)*(n-1))
  c1=1/(m*n)
  c2=(m-1)/(m*n)
  c3=(n-1)/(m*n)
  c4=(m-1)*(n-1)/(m*n)
  covariance=c1*m1+c2*m2+c3*m3-(1-c4)*m4
  covariance=max(covariance,1e-9)
  if ((covariance==1e-9 & ROC1$auc==1) | (covariance==1e-9 & ROC2$auc==1)) {
    covariance=0
    warning("Covariance of two ROC curves with any AUC == 1 is always 0 and can be misleading.")
  }
  return(covariance)
}
