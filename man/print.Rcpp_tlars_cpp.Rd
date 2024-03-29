% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.Rcpp_tlars_cpp.R
\name{print.Rcpp_tlars_cpp}
\alias{print.Rcpp_tlars_cpp}
\title{Prints a summary of the results stored in a C++ object of class tlars_cpp.}
\usage{
\method{print}{Rcpp_tlars_cpp}(x, ...)
}
\arguments{
\item{x}{Object of the class tlars_cpp. See \link{tlars_cpp} for details.}

\item{...}{Ignored. Only added to keep structure of generic \link{print} function.}
}
\value{
Prints a summary of the results stored in a C++ object of class tlars_cpp.
}
\description{
Prints a summary of the results stored in a C++ object of class tlars_cpp
(see \link{tlars_cpp} for details), i.e., selected variables, computation time,
and number of included dummies.
}
\examples{
data("Gauss_data")
X <- Gauss_data$X
y <- drop(Gauss_data$y)
p <- ncol(X)
n <- nrow(X)
num_dummies <- p
dummies <- matrix(stats::rnorm(n * p), nrow = n, ncol = num_dummies)
XD <- cbind(X, dummies)
mod_tlars <- tlars_model(X = XD, y = y, num_dummies = num_dummies)
tlars(model = mod_tlars, T_stop = 3, early_stop = TRUE)
print(mod_tlars)
}
\seealso{
\link{tlars_cpp}.
}
