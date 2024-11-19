#' @name SAPC_LFM
#' @title Apply the SAPC method to the Laplace factor model
#' @description This function calculates several metrics for the SAPC method,
#' including the estimated factor loadings and uniquenesses, and various
#' error metrics comparing the estimated matrices with the true matrices.
#' @param data The data used in the SAPC analysis.
#' @param m The number of common factors.
#' @param A The true factor loadings matrix.
#' @param D The true uniquenesses matrix.
#' @param p The number of variables.
#' @return A list of metrics including:
#' \item{Asa}{Estimated factor loadings matrix obtained from the SAPC analysis.}
#' \item{Dsa}{Estimated uniquenesses vector obtained from the SAPC analysis.}
#' \item{MSESigmaA}{Mean squared error of the estimated factor loadings (Asa) compared to the true loadings (A).}
#' \item{MSESigmaD}{Mean squared error of the estimated uniquenesses (Dsa) compared to the true uniquenesses (D).}
#' \item{LSigmaA}{Loss metric for the estimated factor loadings (Asa), indicating the relative error compared to the true loadings (A).}
#' \item{LSigmaD}{Loss metric for the estimated uniquenesses (Dsa), indicating the relative error compared to the true uniquenesses (D).}
#' @importFrom SOPC SAPC
#' @examples
#' library(SOPC)
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
#' results <- SAPC_LFM(data, m, A, D, p)
#' print(results)
#' @export
#' @importFrom matrixcalc frobenius.norm
#' @importFrom SOPC SAPC
SAPC_LFM <- function(data, m, A, D, p) {
  Asa = SAPC(data, m = m, eta = 0.8)$Asa
  Dsa = SAPC(data, m = m, eta = 0.8)$Dsa
  MSESigmaA = frobenius.norm(Asa - A)^2 / (p^2)
  MSESigmaD = frobenius.norm(Dsa - D)^2 / (p^2)
  LSigmaA = frobenius.norm(Asa - A)^2 / frobenius.norm(A)^2
  LSigmaD = frobenius.norm(Dsa - D)^2 / frobenius.norm(D)^2

  return(list('Asa' = Asa,
              'Dsa' = Dsa,
              'MSESigmaA' = MSESigmaA,
              'MSESigmaD' = MSESigmaD,
              'LSigmaA' = LSigmaA,
              'LSigmaD' = LSigmaD))
}
