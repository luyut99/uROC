#' Construct a ROC curve
#'
#' It builds a ROC curve given the predictor and binary response and saves it to a list including sensitivity,
#' specificity, auc, etc. Unbiased variance of auc and confidence interval can be selected to output.
#' The result can then be passed to ucompROC for testing two correlated ROC curves.
#' @usage ROC(response,predictor,case_score_ind,uvar=FALSE,ci=FALSE,level=NULL)
#' @param response A factor, numeric or character vector of binary responses encoded with 0(controls) and 1(cases).
#' @param predictor A numeric or ordered vector of the same length as response, containing the predicted value of each observation.
#' @param case_score_ind Character variable. Indicating the relation between case's score and control's score. 'higher': cases are
#' associated with higher predictor's value. 'lower': cases are associated with lower predictor's value.
#' @param uvar Logical. Whether the variance of AUC calculated by the unbiased variance estimator is exported. Default is FALSE.
#' @param ci Logical. Whether the CI of the AUC based on unbiased variance estimator is exported. Default is FALSE.
#' @param level A numeric scalar in (0, 1). Indicating the level of confidence interval. Default is 0.95.
#' @details Sensitivity and specificity are calculated by
#' \deqn{sensitivity=\frac{TP}{TP+FN}, specificity=\frac{TN}{TN+FP}} for each given cut-off. AUC is calculated by
#' Mann-Whitney U-statistic
#' \deqn{AUC=\frac{1}{mn}\sum_{i=1}^m\sum_{j=1}^n \psi(X_i, Y_j),} where \eqn{X_i} is predictor's value in case,
#' \eqn{Y_j} is predictor's value in control. When \code{uvar} and \code{ci} are set to be TRUE, \code{uvar} and
#' \code{ci.AUC} will be called, accordingly.
#' @return A list of class "\code{ROC}" including the following components:
#' \describe{
#'   \item{cases}{Original scores in cases}
#'   \item{controls}{Original scores in controls}
#'   \item{case_score_ind}{The case_score_ind used to construct ROC curve}
#'   \item{sensitivity, specificity}{The corresponding sensitivities, specificities at each cutoff}
#'   \item{auc}{The empirical AUC of this ROC curve}
#'   \item{K}{The kernel matrix used to calculate AUC, which would be used to compute variance and covariance}
#'   \item{auc_uvar}{The unbiased variance of AUC if uvar=TRUE}
#'   \item{auc_ci}{The confidence interval of this AUC derived from unbiased variance if ci=TRUE}
#' }
#' @references Hanley, J. A., & McNeil, B. J. (1982). The meaning and use of the area under a receiver operating
#' characteristic (ROC) curve. \emph{Radiology}, 143(1), 29-36.
#' @references Metz, C. E. (1978). Basic principles of ROC analysis. \emph{Seminars in nuclear medicine}. 8(4), 283-298.
#' @references Lu, Y. Shao, Y. (2020). ucompROC: A new powerful test to compare correlated ROC curves.
#' @seealso \code{\link{uvar}}, \code{\link{ci_AUC}}
#' @examples
#' library(uROC)
#'
#' set.seed(123)
#' ## Generate data: m=n=100, tau1=0.8, tau0=0.9, AUC_1=AUC_2=0.9, Clayton copula in cases,
#' ## t copula in controls, both margins are normal
#' Data=SimuCaseCont(100,100,0.8,0.9,0.9,0.9,"c","t",marg1="Exp",marg2="Exp")
#' ## Construct ROC curve by default
#' roc1=ROC(Data[,3],Data[,1],case_score_ind = "higher",uvar=TRUE)
#' ## Constrct ROC curve with 95% CI
#' roc2=ROC(Data[,3],Data[,2],case_score_ind = "higher",ci=TRUE,level=0.95)
#' roc1
#' roc2
#' @export

ROC<-function(response,predictor,case_score_ind=c("higher","lower"),uvar=FALSE,ci=FALSE,level=NULL) {
  if (is.null(case_score_ind) | missing(case_score_ind)) {
    stop("The relation between case's score and control's score need to be specified.")
  } else {
    case_score_ind<-match.arg(case_score_ind,c("higher","lower"))
  }

  # Create roc object by pROC, calculate auc, sensitivity, specificity, etc
  case=predictor[response==1]
  cont=predictor[response==0]

  if (case_score_ind=="higher") {
    response_sort=response[order(predictor,decreasing = T)]
    kern_func = function(a,b) (sign(a-b)+1)/2 # AUC kernel function
  } else {
    response_sort=response[order(predictor,decreasing = F)]
    kern_func = function(a,b) (sign(b-a)+1)/2 # AUC kernel function
  }

  sen=cumsum(response_sort==1)/length(case) # Sensitivity
  spe=1-cumsum(response_sort==0)/length(cont) # Specificity
  # TP=cumsum(response_sort==1) # True positive
  # TN=length(cont)-cumsum(response_sort==0) # True negative

  K = outer(case,cont,kern_func)
  AUC = mean(K)

  ROC<-list()
  ROC$cases<-case
  ROC$controls<-cont
  ROC$case_score_ind<-case_score_ind
  ROC$sensitivity<-sen
  ROC$specificity<-spe
  ROC$auc<-AUC
  ROC$K<-K
  # Whether calculate AUC variance
  if (uvar) {
    ROC$auc_uvar=uvar(ROC)
  }
  # Whether calculate AUC CI
  if (ci) {
    if (missing(level) | is.null(level)) {
      level=0.95 # stop("Level for CI need to be specified.")
    } else if (level<0 | level>1) stop("Level must be in (0, 1).")
    ROC$auc_ci=ci_AUC(ROC,level=level)
  }
  class(ROC)<-"ROC"
  return(ROC)
}
