\name{pcd}
\alias{pcd}
\docType{data}
\title{Pancreatic cancer dataset}
\description{
Data of 51 'control' patients with pancreatitis and 90 'case' patients with pancreatic cancer that were studied
#' at the Mayo Clinic with a cancer antigen (CA125) and with a carbohydrate antigen (CA19-9).
}
\usage{pcd}
\format{
  A data frame with 141 observations on the following 3 variables.
  \describe{
    \item{\code{y1}}{CA19-9}
    \item{\code{y2}}{CA125}
    \item{\code{d}}{Status indicator. 'd=1' refers to pancreatic case, 'd=0' refers to control}
  }
}
\source{
\href{https://research.fhcrc.org/diagnostic-biomarkers-center/en/datasets.html}{https://research.fhcrc.org/diagnostic-biomarkers-center/en/datasets.html}
}
\references{
Wieand S, Gail MH, James BR, and James KL. A family of nonparametric statistics for comparing
#' diagnostic markers with paired or unpaired data. Biometrika 76(3):585-92. 1989.(\href{https://www.jstor.org/stable/2336123}{https://www.jstor.org/stable/2336123})
}
\examples{
library(uROC)

data(pcd)
## Calculate ROC and AUC for each biomarker
roc1=ROC(pcd$d,pcd$y1,case_score_ind="higher",uvar=TRUE)
roc2=ROC(pcd$d,pcd$y2,case_score_ind="higher",uvar=TRUE)
## Compare the two biomarkers by hypothesis testing
ucompROC(roc1,roc2,alternative="two.sided",var_diff=TRUE,ci=TRUE)
}
\keyword{datasets}
