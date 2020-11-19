#' Confidence interval of an AUC using unbiased variance estimator
#'
#' Calculate the two-sided confidence interval of the AUC under a single ROC curve using normal approximation
#' with unbiased AUC variance estimator.
#' @usage ci_AUC(ROC,level=NULL)
#' @param ROC An \code{ROC} object obtained by ROC.
#' @param level A numeric scalar in (0, 1). Indicating the level of confidence interval. Default is 0.95.
#' @return A vector of length 2 including the lower bound and upper bound of the confidence interval.
#' @references Lu, Y. and Shao, Y. (2020). ucompROC: A new powerful test to compare correlated ROC curves.
#' @seealso \code{\link{ROC}}, \code{\link{uvar}}
#' @section Warning:
#' Confidence interval for AUC=1 has width 0 and can be misleading.
#' @examples
#' library(uROC)
#' set.seed(123)
#' ## Generate data
#' Data=SimuCaseCont(100,100,0.8,0.9,0.9,0.9,"c","t",marg1="Normal",marg2="Normal")
#' ## Construct ROC curves
#' roc1=ROC(Data[,3],Data[,1],case_score_ind = "higher")
#' roc2=ROC(Data[,3],Data[,2],case_score_ind = "higher")
#' ## 95% CI of roc1
#' ci_AUC(roc1,level=0.95)
#' ## 90% CI of roc2
#' ci_AUC(roc2,level=0.9)
#' @export

ci_AUC<-function(ROC,level=NULL) {
  if (ROC$auc==1) {
    warning("Confidence interval of a ROC curve with AUC == 1 has width 0 and can be misleading.")
  }
  if (missing(level) | is.null(level)) {
    level=0.95
  } else if (level<0 | level>1) stop("Level must be in (0, 1).")
  std=sqrt(suppressWarnings(uvar(ROC)))
  lb=ROC$auc-abs(qnorm((1-level)/2))*std
  ub=ROC$auc+abs(qnorm((1-level)/2))*std
  ci=round(c(lb,ub),3)
  return(ci)
}
