#' Confidence interval of AUC difference using unbiased variance estimator
#'
#' Calculate the confidence interval of the AUC difference using normal approximation with the unbiased variance estimator
#' @usage ci_diff(ROC1,ROC2,type=c("two.sided","less","greater"),level=NULL)
#' @param ROC1,ROC2 Two \code{ROC} objects referring two ROC curves.
#' @param type Character. Type of the confidence interval. Must be one of 'two.sided', 'less', 'greater'. 'less',
#' 'greater' refer to one-sided confidence interval corresponding to the alternatives. The first letter is
#' sufficient.
#' @param level A numeric scalar in (0, 1). Indicating the level of confidence interval. Default is 0.95.
#' @return A vector of length 2 including the lower bound and upper bound of the confidence interval.
#' @references Lu, Y. and Shao, Y. (2020). ucompROC: A new powerful test to compare correlated ROC curves.
#' @seealso \code{\link{ucompROC}}
#' @section Warning:
#' When \eqn{AUC_1=AUC_2}, the variance of the AUC difference may be 0, so the confidence interval should always
#' have width 0 and can be misleading.
#' @examples
#' library(uROC)
#' set.seed(321)
#' ## Generate data
#' Data=SimuCaseCont(30,30,0.9,0.9,0.87,0.90,"c","t",marg1="Normal",marg2="Normal")
#' ## Construct ROC curves
#' uroc1=ROC(Data[,3],Data[,1],case_score_ind = "higher")
#' uroc2=ROC(Data[,3],Data[,2],case_score_ind = "higher")
#' ## Two-sided confidence interval
#' ci_diff(uroc1,uroc2,type = "two.sided")
#'
#' ## One-sided confidence interval
#' ci_diff(uroc1,uroc2,type="less",level=0.95)
#' @export

ci_diff<-function(ROC1,ROC2,type=c("two.sided","less","greater"),level=NULL){
  if (missing(type) | is.null(type)) {
    type='two.sided' # Default is unbiased
  } else {
    type <- match.arg(type,c("two.sided","less","greater"))
  }
  if (missing(level) | is.null(level)) {
    level=0.95
  } else if (level<0 | level>1) stop("Level must be in (0, 1).")
  AUC_diff=ROC1$auc-ROC2$auc
  variance=suppressWarnings(uvar(ROC1)+uvar(ROC2)-2*ucov(ROC1,ROC2))
  if (variance==0) {warning("Variance of AUC difference is 0. Confidence interval can be misleading.")
    std=0
  } else {
    std=sqrt(variance)
  }
  if (type=="two.sided") {
    lb=AUC_diff-abs(qnorm((1-level)/2))*std
    ub=AUC_diff+abs(qnorm((1-level)/2))*std
  } else if (type=="less") {
    lb=-0.5
    ub=AUC_diff+abs(qnorm(level))*std
  } else {
    lb=AUC_diff-abs(qnorm(level))*std
    ub=0.5
  }
  # print(paste0(level*100,"% CI: ",round(lb,3),"-",round(ub,3)," (unbiased)"))
  ci<-round(c(lb,ub),3)
  return(ci)
}
