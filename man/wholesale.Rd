\name{wholesale}
\alias{wholesale}
\docType{data}
\title{Wholesale Customers Data}
\description{
This dataset contains the annual spending amounts of wholesale customers on various product categories, along with their channel and region information.
}
\usage{wholesale}
\format{
A data frame with 440 rows and 8 columns.
\itemize{
\item \code{FRESH}: Annual spending (m.u.) on fresh products.
\item \code{MILK}: Annual spending (m.u.) on milk products.
\item \code{GROCERY}: Annual spending (m.u.) on grocery products.
\item \code{FROZEN}: Annual spending (m.u.) on frozen products.
\item \code{DETERGENTS_PAPER}: Annual spending (m.u.) on detergents and paper products.
\item \code{DELICATESSEN}: Annual spending (m.u.) on delicatessen products.
\item \code{CHANNEL}: Customers' channel - Horeca (Hotel/Restaurant/Café) or Retail channel (Nominal).
\item \code{REGION}: Customers' region - Lisbon, Oporto or Other (Nominal).
}
}
\details{
The dataset includes both continuous spending variables and nominal classification variables (CHANNEL and REGION). Descriptive statistics for the continuous variables are provided, including minimum, maximum, mean, and standard deviation.
}
\note{
The data can be used for clustering analysis to segment customers based on their purchasing behavior.
}
\source{
UCI Machine Learning Repository
}
\examples{
# Load the dataset
data(wholesale)
}
\keyword{datasets}
