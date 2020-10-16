#' Hypothesis test for comparing AUCs of two correlated ROC curves using unbiased variance estimator
#'
#' Test whether the AUCs of two correlated ROC curves are different.
#' @usage ucompROC(ROC1,ROC2,alternative=c("two.sided","less","greater"),
#' var_diff=FALSE,ci=FALSE,...)
#' @param ROC1,ROC2 Two \code{ROC} objects referring two ROC curves.
#' @param alternative Character. Indicate the alternative hypothesis. Must be one of 'two.sided', 'less' or 'greater'.
#' Default is 'two.sided'. The first letter is sufficient.
#' @param var_diff Logical. Whether the unbiased variance of AUC difference is calculated and exported. Default is FALSE.
#' @param ci Logical. Whether the confidence interval of AUC difference is calculated and exported. Default is FALSE.
#' If TRUE, the type of the confidence interval will be corresponding to the alternative, and the level could be specified.
#' Default of level is 0.95 if not specified.
#' @param ... Optional arguments passed to \code{\link{ci_diff}}, such as level.
#' @details The test statistic is
#' \deqn{Z=\frac{AUC_1-AUC_2}{var(AUC_1-AUC_2)},} where \eqn{var(AUC_1-AUC_2)=var(AUC_1)+var(AUC_2)
#' -2cov(AUC_1, AUC_2)}, variance and covariance are estimated by the unbiased estimator in this test. Under null
#' hypothesis, the asymptotic distribution of \eqn{Z} is normal distribution \eqn{N(0, 1)}.
#'
#' @return A list of result including the following component:
#' \describe{
#'   \item{null_hypo}{The AUC difference under null hypothesis}
#'   \item{alternative}{The alternative hypothesis}
#'   \item{estimates}{The AUCs of the two ROC curves}
#'   \item{AUC_diff}{The AUC difference}
#'   \item{Z}{Test statistic}
#'   \item{p.value}{P-value}
#'   \item{var_diff}{Unbiased variance estimate of the AUC difference if var_diff=TRUE}
#'   \item{ci_diff}{Confidence interval of the AUC difference calculated by normal approximation with
#'   the unbiased variance estimator if ci_diff=TRUE}
#' }
#'
#' @references Lu, Y. and Shao, Y. (2020). ucompROC: A new powerful test to compare correlated ROC curves.
#' @seealso \code{\link{ROC}}, \code{\link{uvar}}, \code{\link{ucov}}, \code{\link{ci_diff}}
#' @section Warning:
#' If \eqn{AUC_1 = AUC_2}, variance of AUC difference may be 0 and misleading. But this test is still valid.
#' @examples
#' set.seed(321)
#' ## Generate data
#' ## m=n=30, tau1=0.9, tau0=0.9, AUC_1=0.87, AUC_2=0.90, Clayton and t copula for cases and controls
#' ## Normal margins for both biomarkers
#' Data=SimuCaseCont(30,30,0.9,0.9,0.87,0.90,"c","t",marg1="Normal",marg2="Normal")
#'
#' ## New test using unbiased variance estimator, result in p-value<0.05
#' uroc1=ROC(Data[,3],Data[,1],case_score_ind = "higher")
#' uroc2=ROC(Data[,3],Data[,2],case_score_ind = "higher")
#' Utest=ucompROC(uroc1,uroc2,alternative = "two.sided",var_diff = TRUE)
#' Utest
#'
#' ## DeLong test, result in p-value>0.05
#' \dontrun{
#' library(pROC)
#' droc1=suppressMessages(roc(Data[,3],Data[,1],direction = "<"))
#' droc2=suppressMessages(roc(Data[,3],Data[,2],direction = "<"))
#' Dtest=roc.test(droc1,droc2,method="delong",paired = TRUE)
#' Dtest
#'
#' ## The variance estimated by our test
#' Utest$var_diff
#'
#' ## The variance used in DeLong test, which is 41% greather than the unbiased variance
#' var(droc1)+var(droc2)-2*cov(droc1,droc2)}
#' @export

ucompROC<-function(ROC1,ROC2,alternative=c("two.sided","less","greater"),
                   var_diff=FALSE,ci=FALSE,...) { #paired=NULL,
  if (missing(alternative) | is.null(alternative)) {
    alternative='two.sided' # Default is two.sided
  } else {
    alternative <- match.arg(alternative,c("two.sided","less","greater"))
  }

  # if (is.null(ROC1$direction) | is.null(ROC2$direction)) {
  #   stop("Directions of both curves have to be specified.")
  # }
  null_hypo=0
  names(null_hypo)="difference in AUC"

  estimates=c(ROC1$auc,ROC2$auc)
  names(estimates)=c("AUC of ROC 1", "AUC of ROC 2")
  AUC_diff=ROC1$auc-ROC2$auc
  varAUC_diff=suppressWarnings(uvar(ROC1)+uvar(ROC2)-2*ucov(ROC1,ROC2))
  # Calculate AUCs
  varAUC_diff=ifelse((AUC_diff==0 & varAUC_diff==0) | varAUC_diff<=0, 1e-12,varAUC_diff)
  Z=AUC_diff/sqrt(varAUC_diff)
  if (alternative=="two.sided") {
    p.value=2*pnorm(-abs(Z))
  } else if (alternative=="less") {
    p.value=pnorm(Z)
  } else {
    p.value=1-pnorm(Z)
  }
  result<-list()
  result$null_hypo=null_hypo
  result$alternative=alternative
  result$AUC_estimates=estimates
  result$AUC_difference=AUC_diff
  result$Z=Z
  result$p.value=p.value
  if (var_diff==T) {
    if (varAUC_diff==1e-12) {
      varAUC_diff=0
      warning("If AUC1 = AUC2, variance of AUC difference may be 0 and misleading.")
    }
    result$var_diff=varAUC_diff
  }
  if (ci==T) {
    result$ci_diff=ci_diff(ROC1,ROC2,type=alternative,...)
  }
  result$method="Asymptotic Z test for two correlated ROC curves using unbiased variance"
  return(result)
}
