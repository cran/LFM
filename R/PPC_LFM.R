#' @name PPC_LFM
#' @title Apply the PPC method to the Laplace factor model
#' @description This function computes Projected Principal Component Analysis (PPC) for the provided input data, estimating factor loadings and uniquenesses. It calculates mean squared errors and loss metrics for the estimated values compared to true values.
#' @param data A matrix of input data.
#' @param m The number of principal components.
#' @param A The true factor loadings matrix.
#' @param D The true uniquenesses matrix.
#' @param p The number of variables.
#' @return A list containing:
#' \item{Ap}{Estimated factor loadings.}
#' \item{Dp}{Estimated uniquenesses.}
#' \item{MSESigmaA}{Mean squared error for factor loadings.}
#' \item{MSESigmaD}{Mean squared error for uniquenesses.}
#' \item{LSigmaA}{Loss metric for factor loadings.}
#' \item{LSigmaD}{Loss metric for uniquenesses.}
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
#' results <- PPC_LFM(data, m, A, D, p)
#' print(results)
#' @export
#' @importFrom SOPC PPC
#' @importFrom matrixcalc frobenius.norm
PPC_LFM <- function(data, m, A, D, p) {
  Ap = PPC(data, m = m, eta = 0.8)$Ap
  Dp = PPC(data, m = m, eta = 0.8)$Dp
  MSESigmaA = frobenius.norm(Ap - A)^2 / (p^2)
  MSESigmaD = frobenius.norm(Dp - D)^2 / (p^2)
  LSigmaA = frobenius.norm(Ap - A)^2 / frobenius.norm(A)^2
  LSigmaD = frobenius.norm(Dp - D)^2 / frobenius.norm(D)^2

  return(list(Ap = Ap,
              Dp = Dp,
              MSESigmaA = MSESigmaA,
              MSESigmaD = MSESigmaD,
              LSigmaA = LSigmaA,
              LSigmaD = LSigmaD))
}

