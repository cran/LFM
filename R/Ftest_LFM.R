#' @name Ftest_LFM
#' @title Apply the Farmtest method to the Laplace factor model
#' @description This function simulates data from a Lapalce factor model and applies the FarmTest
#' for multiple hypothesis testing. It calculates the false discovery rate (FDR)
#' and power of the test.
#' @param data A matrix or data frame of simulated or observed data from a Laplace factor model.
#' @param p1 The proportion of non-zero hypotheses.
#' @return A list containing the following elements:
#' \item{FDR}{The false discovery rate, which is the proportion of false positives among all discoveries (rejected hypotheses).}
#' \item{Power}{The statistical power of the test, which is the probability of correctly rejecting a false null hypothesis.}
#' \item{PValues}{A vector of p-values associated with each hypothesis test.}
#' \item{RejectedHypotheses}{The total number of hypotheses that were rejected by the FarmTest.}
#' @examples
#' library(LaplacesDemon)
#' library(MASS)
#' n=1000
#' p=10
#' m=5
#' mu=t(matrix(rep(runif(p,0,1000),n),p,n))
#' mu0=as.matrix(runif(m,0))
#' sigma0=diag(runif(m,1))
#' F=matrix(mvrnorm(n,mu0,sigma0),nrow=n)
#' A=matrix(runif(p*m,-1,1),nrow=p)
#' lanor <- rlaplace(n*p,0,1)
#' epsilon=matrix(lanor,nrow=n)
#' D=diag(t(epsilon)%*%epsilon)
#' data=mu+F%*%t(A)+epsilon
#' p1=40
#' results <- Ftest_LFM(data, p1)
#' print(results$FDR)
#' print(results$Power)
#' @export
#' @importFrom FarmTest farm.test
#'
Ftest_LFM <- function(data, p1) {

  # Apply FarmTest
  output <- farm.test(data)

  # Calculate FDR and power
  fdr <- mean(output$reject > p1) / length(output$reject)
  power <- mean(output$reject <= p1) / p1

  # Extract p-values and rejected hypotheses
  pValues <- output$p.value
  RejectedHypotheses <- sum(output$reject)

  # Return the results as a list
  return(list(FDR = fdr,
              Power = power,
              PValues = pValues,
              RejectedHypotheses = RejectedHypotheses))
}
