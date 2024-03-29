% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tlars.R
\name{tlars}
\alias{tlars}
\title{Executes the Terminating-LARS (T-LARS) algorithm}
\usage{
tlars(model, T_stop = 1, early_stop = TRUE, info = TRUE)
}
\arguments{
\item{model}{Object of the class tlars_cpp.}

\item{T_stop}{Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.}

\item{early_stop}{Logical. If TRUE, then the forward selection process is stopped after T_stop dummies have been included. Otherwise
the entire solution path is computed.}

\item{info}{If TRUE information about the T-LARS step are printed.}
}
\value{
No return value. Executes the T-LARS algorithm and includes the results in the associated object of class tlars_cpp.
}
\description{
Modifies the generic tlars_cpp model by executing the T-LARS algorithm and including the results in the tlars_cpp model.
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
beta <- mod_tlars$get_beta()
beta
}
