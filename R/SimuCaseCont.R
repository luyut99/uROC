#' Random generation of two correlated risk scores under case-control settings
#'
#' Generate two correlated risk score with given AUCs and correlations under case-control setting using copula theory.
#' @usage SimuCaseCont(m,n,tau1,tau0,AUC_1,AUC_2,
#'   cop1=c("g","t","c"),cop0=c("g","t","c"),marg1=c("Normal","Exp"),marg2=c("Normal","Exp"))
#' @param m,n Integers. Specify the number of cases(m), number of controls(n).
#' @param tau1,tau0 Numeric. The Kendall's tau correlation between the two risk scores in cases(tau1)
#' and controls(tau0) separately. Should be specified in (0, 1).
#' @param AUC_1,AUC_2 The AUC of each of the risk scores. Should be specified (0.5, 1)
#' @param cop1,cop0 Character. Specify the copula for cases(cop1) and controls(cop0). Choosing from
#' 'g'(Gaussian copula), 't'(t copula), and 'c'(Clayton copula).
#' @param marg1,marg2 Character. Marginal distributions for two risk scores, respectively. Both choose from "Normal"
#' or "Exp"(exponential).
#' @return A \eqn{(m+n) \times 3} data matrix. The first two columns are the generated values for risk score 1
#' and 2, accordingly, and the last column is the binary response.
#' @references Nagler T., Schepsmeier U., Stoeber J. , Brechmann E C., Graeler B. and Erhardt T. (2019).
#' \emph{VineCopula: Statistical Inference of Vine Copulas.} R package version 2.3.0.
#' https://CRAN.R-project.org/package=VineCopula
#' @examples
#' library(uROC)
#'
#' set.seed(123)
#' ## Generate data: m=n=100, tau1=0.8, tau0=0.9, AUC_1=AUC_2=0.9, Clayton copula in cases,
#' ## t copula in controls, both normal margins for score 1 and 2
#' Data=SimuCaseCont(100,100,0.8,0.9,0.9,0.9,"c","t",marg1="Normal",marg2="Normal")
#' Data
#' ## Check Kendall's tau correlation
#' cor(Data[1:100,1],Data[1:100,2],method="kendall")
#' cor(Data[101:200,1],Data[101:200,2],method="kendall")
#' ## Check AUC
#' roc1=ROC(Data[,3],Data[,1],case_score_ind = "higher")
#' roc2=ROC(Data[,3],Data[,2],case_score_ind = "higher")
#' roc1$auc
#' roc2$auc
#' @importFrom VineCopula BiCopTau2Par BiCopSim
#' @importFrom stats cor qnorm qexp pnorm
#' @export

SimuCaseCont <- function(m,n,tau1,tau0,AUC_1,AUC_2,cop1=c("g","t","c"),cop0=c("g","t","c"),marg1=c("Normal","Exp"),marg2=c("Normal","Exp")) {
  if (tau0<0 | tau0>1 | tau1<0 | tau1>1) stop("tau0, tau1 should be specified in (0, 1).")
  if (AUC_1<0.5 | AUC_1>1 | AUC_2<0.5 | AUC_2>1) stop("AUC_1, AUC_2 should be specified in (0.5, 1).")
  cop0<-match.arg(cop0,c("g","t","c"))
  cop1<-match.arg(cop1,c("g","t","c"))
  cop0n<-ifelse(cop0=="g",1,ifelse(cop0=="t",2,3))
  cop1n<-ifelse(cop1=="g",1,ifelse(cop1=="t",2,3))
  marg1<-match.arg(marg1,c("Normal","Exp"))
  marg2<-match.arg(marg2,c("Normal","Exp"))
  Z0=rep(0,n)
  Z1=rep(1,m)
  Z=c(Z0,Z1)
  ## Given correlation, generate bivariate F(x|z), G(y|z) conditioning on Z##
  if (cop0n==2) {
    rho0=BiCopTau2Par(cop0n,tau0)
    uv0=BiCopSim(n,cop0n,rho0,3)
  } else {
    rho0=BiCopTau2Par(cop0n,tau0)
    uv0=BiCopSim(n,cop0n,rho0)
  }
  if (cop1n==2) {
    rho1=BiCopTau2Par(cop1n,tau1)
    uv1=BiCopSim(m,cop1n,rho1,3)
  } else {
    rho1=BiCopTau2Par(cop1n,tau1)
    uv1=BiCopSim(m,cop1n,rho1)
  }

  ## Given AUC, assign conditional margin and generate scores X, Y based on Z ##
  ## Two margins: Normal and Exponential ##
  if (marg1=="Normal") {
    mu_cont=0
    mu_case=sqrt(2)*qnorm(AUC_1)
    X_cont=qnorm(uv0[,1],mu_cont)
    X_case=qnorm(uv1[,1],mu_case)
  } else {
    b1=(1-AUC_1)/AUC_1
    X_cont=qexp(uv0[,1],1)
    X_case=qexp(uv1[,1],b1)
  }

  if (marg2=="Normal") {
    nu_cont=0
    nu_case=sqrt(2)*qnorm(AUC_2)
    Y_cont=qnorm(uv0[,2],nu_cont)
    Y_case=qnorm(uv1[,2],nu_case)
  } else {
    b=(1-AUC_2)/AUC_2
    Y_cont=qexp(uv0[,2],1)
    Y_case=qexp(uv1[,2],b)
  }
  X=c(X_cont,X_case)
  Y=c(Y_cont,Y_case)
  Data=cbind(X,Y,Z)
  return(Data)
}
