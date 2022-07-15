## ---- include = FALSE---------------------------------------------------------
# Store user's options()
old_options <- options()

library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  fig.retina = 2,
  out.width = "85%",
  dpi = 96
  # pngquant = "--speed=1"
)
options(width = 80)

## ---- eval=FALSE--------------------------------------------------------------
#  library(tlars)
#  help(package = "tlars")
#  ?tlars
#  ?tlars_model
#  ?tlars_cpp
#  ?plot.Rcpp_tlars_cpp
#  ?print.Rcpp_tlars_cpp
#  ?Gauss_data

## ---- eval=FALSE--------------------------------------------------------------
#  citation("tlars")

## -----------------------------------------------------------------------------
library(tlars)

# Setup
n <- 150 # Number of observations
p <- 300 # Number of variables
num_act <- 5 # Number of true active variables
beta <- c(rep(1, times = num_act), rep(0, times = p - num_act)) # Coefficient vector
true_actives <- which(beta > 0) # Indices of true active variables
num_dummies <- p # Number of dummy predictors (or dummies)

# Generate Gaussian data
set.seed(123)
X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
y <- X %*% beta + stats::rnorm(n)

## -----------------------------------------------------------------------------
set.seed(1234)
dummies <- matrix(stats::rnorm(n * num_dummies), nrow = n, ncol = num_dummies)
XD <- cbind(X, dummies)

## -----------------------------------------------------------------------------
mod_tlars <- tlars_model(X = XD, y = y, num_dummies = num_dummies)

## ----terminated_solution_path_T_3, echo=TRUE, fig.align='center', message=TRUE, fig.width = 10, fig.height = 6, out.width = "90%"----
tlars(model = mod_tlars, T_stop = 1, early_stop = TRUE) # Perform one T-LARS step on object "mod_tlars"
print(mod_tlars) # Print information about the results of the performed T-LARS steps
plot(mod_tlars) # Plot the terminated solution path

## ----terminated_solution_path_T_5, echo=TRUE, fig.align='center', message=TRUE, fig.width = 10, fig.height = 6, out.width = "90%"----
# Numerical zero
eps <- .Machine$double.eps

# Perform one additional T-LARS step (going from T_stop = 1 to T_stop = 2) on object "mod_tlars"
tlars(model = mod_tlars, T_stop = 2, early_stop = TRUE)
print(mod_tlars)
plot(mod_tlars)

# Coefficient vector corresponding to original and dummy variables at the terminal T-LARS step 
beta_hat <- mod_tlars$get_beta()

selected_var <- which(abs(beta_hat[seq(p)]) > eps) # Indices of selected original variables
selected_dummies <- p + which(abs(beta_hat[seq(p + 1, ncol(XD))]) > eps) # Indices of selected dummy variables

FDP <- length(setdiff(selected_var, true_actives)) / max(1, length(selected_var)) # False discovery proportion (FDP)
TPP <- length(intersect(selected_var, true_actives)) / max(1, length(true_actives)) # True positive proportion (TPP)

selected_var
selected_dummies
FDP
TPP

## ----terminated_solution_path_T_10, echo=TRUE, fig.align='center', message=TRUE, fig.width = 10, fig.height = 6, out.width = "90%"----
# Perform three additional T-LARS steps (going from T_stop = 2 to T_stop = 5) on object "mod_tlars"
tlars(model = mod_tlars, T_stop = 5, early_stop = TRUE)
print(mod_tlars)
plot(mod_tlars)

# Coefficient vector corresponding to original and dummy variables at the terminal T-LARS step 
beta_hat <- mod_tlars$get_beta()

selected_var <- which(abs(beta_hat[seq(p)]) > eps) # Indices of selected original variables
selected_dummies <- p + which(abs(beta_hat[seq(p + 1, ncol(XD))]) > eps) # Indices of selected dummy variables

FDP <- length(setdiff(selected_var, true_actives)) / max(1, length(selected_var)) # False discovery proportion (FDP)
TPP <- length(intersect(selected_var, true_actives)) / max(1, length(true_actives)) # True positive proportion (TPP)

selected_var
selected_dummies
FDP
TPP

## -----------------------------------------------------------------------------
# Setup
n <- 100 # number of observations
p <- 300 # number of variables

# Parameters
num_act <- 10 # number of true active variables
beta <- rep(0, times = p) # coefficient vector (all zeros first)
beta[sample(seq(p), size = num_act, replace = FALSE)] <- 3 # coefficient vector (active variables with non-zero coefficients)
true_actives <- which(beta > 0) # indices of true active variables
num_dummies <- p # number of dummies
T_vec <- c(1, 2, 5, 10, 20, 50, 100) # stopping points, i.e, number of included dummies before terminating the solution path
MC <- 500 # number of Monte Carlo runs per stopping point

# Initialize results vectors
FDP <- matrix(NA, nrow = MC, ncol = length(T_vec))
TPP <- matrix(NA, nrow = MC, ncol = length(T_vec))

# Numerical zero
eps <- .Machine$double.eps

# Seed
set.seed(12345)

# Run simulations
for (t in seq_along(T_vec)) {
  for (mc in seq(MC)) {
    # Generate Gaussian data
    X <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    y <- X %*% beta + stats::rnorm(n)

    # Generate dummy matrix and append it to X
    dummies <- matrix(stats::rnorm(n * p), nrow = n, ncol = num_dummies)
    XD <- cbind(X, dummies)

    # Create object of class tlars_cpp
    mod_tlars <- tlars_model(X = XD, y = y, num_dummies = num_dummies, type = "lar", info = FALSE)

    # Run T-LARS steps
    tlars(model = mod_tlars, T_stop = t, early_stop = TRUE, info = FALSE)
    beta_hat <- mod_tlars$get_beta()
    selected_var <- which(abs(beta_hat[seq(p)]) > eps)

    # Results
    FDP[mc, t] <- length(setdiff(selected_var, true_actives)) / max(1, length(selected_var))
    TPP[mc, t] <- length(intersect(selected_var, true_actives)) / max(1, length(true_actives))
  }
}

# Compute estimates of FDR and TPR by averaging FDP and TPP over MC Monte Carlo runs
FDR <- colMeans(FDP)
TPR <- colMeans(TPP)

## ----FDR_and_TPR, echo=FALSE, fig.align='center', message=FALSE, fig.width = 10, fig.height = 5, out.width = "90%"----
# Plot results
library(ggplot2)
library(patchwork)

plot_data = data.frame(T_vec = T_vec,
                       FDR = 100 * FDR,
                       TPR = 100 * TPR) # data frame containing data to be plotted (FDR and TPR in %)

# FDR vs. T
FDR_vs_T = 
  ggplot(plot_data, aes(x = T_vec, y = FDR)) +
    labs(x = "T",
         y = "FDR") +
    scale_x_continuous(breaks = T_vec[-2], minor_breaks = c(2), limits = c(T_vec[1], T_vec[length(T_vec)])) +
    scale_y_continuous(breaks = seq(0, 100, by = 10), minor_breaks = c(), limits = c(0, 100)) +
    geom_line(size = 1.5, colour = "#336C68") +
    geom_point(size = 2.5, colour = "#336C68") +
    theme_bw(base_size = 16) +
    theme(panel.background = element_rect(fill = "white", color = "black", size = 1)) +
    coord_fixed(ratio =  0.85 * T_vec[length(T_vec)] / (100 - 0))

# TPR vs. T
TPR_vs_T = 
  ggplot(plot_data, aes(x = T_vec, y = TPR)) +
    labs(x = "T",
         y = "TPR") +
    scale_x_continuous(breaks = T_vec[-2], minor_breaks = c(2), limits = c(T_vec[1], T_vec[length(T_vec)])) +
    scale_y_continuous(breaks = seq(0, 100, by = 10), minor_breaks = c(), limits = c(0, 100)) +
    geom_line(size = 1.5, colour = "#336C68") +
    geom_point(size = 2.5, colour = "#336C68") +
    theme_bw(base_size = 16) +
    theme(panel.background = element_rect(fill = "white", color = "black", size = 1)) +
    coord_fixed(ratio =  0.85 * T_vec[length(T_vec)] / (100 - 0))

FDR_vs_T + TPR_vs_T

# TPR vs. FDR
TPR_vs_FDR = 
  ggplot(plot_data, aes(x = FDR, y = TPR)) +
    labs(x = "FDR",
         y = "TPR") +
    geom_line(size = 1.5, colour = "#336C68") +
    geom_point(size = 2.5, colour = "#336C68") +
    theme_bw(base_size = 16) +
    theme(panel.background = element_rect(fill = "white", color = "black", size = 1))

TPR_vs_FDR

## ---- include = FALSE---------------------------------------------------------
# Reset user's options()
options(old_options)

