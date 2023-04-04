###############################################################################
#
#   Calculate statistical significancy for the Activity Zscores:
#   1) Uniform+Normal mixture model to infer Gaussian parameters for the bg
#   2) Standardize Z_{m,s} with the MLE of parameters
#   3) Calculate Bonferroni-corrected p-values
#   4) Combine standardized per-sample Z_{m,s} into Z_m as in ISMARA model
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 20-01-2020
#   LICENSE: Apache_2.0
#
###############################################################################

# by default: suppress warnings
options(warn = -1)

# load libraries
suppressPackageStartupMessages(library("optparse"))

# list the command-line arguments
option_list <- list(
  make_option(c("--input"),
    action = "store",
    dest = "input",
    type = "character",
    help = "Path to the input Z-scores table."
  ),
  make_option(c("--sorting-strategy"),
    action = "store",
    dest = "sorting_strategy",
    type = "character",
    help = "Sorting strategy for the output table; options: avg OR max."
  ),
  make_option(c("--output"),
    action = "store",
    dest = "output",
    type = "character",
    help = "Path for the updated Z-scores table."
  ),
  make_option(c("--help"),
    action = "store_true",
    dest = "help",
    type = "logical",
    default = FALSE,
    help = "Show this information and exit."
  ),
  make_option(c("--verbose"),
    action = "store_true",
    dest = "verbose",
    type = "logical",
    default = FALSE,
    help = "Run in verbose mode."
  )
)

# parse command-line arguments
opt_parser <- OptionParser(
  usage = "",
  option_list = option_list,
  add_help_option = FALSE,
  description = ""
)
opt <- parse_args(opt_parser)

# if verbose flag was set: print warnings
if (opt$verbose) {
  options(warn = 0)
}

# -----------------------------------------------------------------------------


# Mixture model:
# Probability to get a particular z-score given rho, mu, sigma
P_z <- function(z, mu, sigma, rho, z_max, z_min) {
  # uniform term:
  u <- rho / (z_max - z_min)
  # gaussian term:
  n <- (1 - rho) * exp(-(z - mu)^2 / (2 * sigma^2)) / (sqrt(2 * pi) * sigma)
  return(u + n)
}


# Log-likelihood is the sum over all z-scores
LL <- function(z, mu, sigma, rho, z_max, z_min) {
  return(sum(sapply(z, function(x) {
    log(P_z(x, mu, sigma, rho, z_max, z_min))
  })))
}


# Posterior for motif coming from the 'foreground':
# ratio of the uniform term to the uniform+gaussian
Posterior_for_fg <- function(z, mu, sigma, rho, z_max, z_min) {
  # uniform term:
  u <- rho / (z_max - z_min)
  # gaussian term:
  n <- (1 - rho) * exp(-(z - mu)^2 / (2 * sigma^2)) / (sqrt(2 * pi) * sigma)
  return(u / (u + n))
}


# Posterior for motif coming from the 'background':
# ratio of the gaussian term to the uniform+gaussian
Posterior_for_bg <- function(z, mu, sigma, rho, z_max, z_min) {
  # uniform term:
  u <- rho / (z_max - z_min)
  # gaussian term:
  n <- (1 - rho) * exp(-(z - mu)^2 / (2 * sigma^2)) / (sqrt(2 * pi) * sigma)
  return(n / (u + n))
}


# Expectation-Maximisation to infer the parameters of the mixture model
# C++ do-while is simulated by R repeat-if-break structure
EM <- function(Zscores) {

  # define the requred variables:
  z_max <- max(Zscores)
  z_min <- min(Zscores)
  n_motifs <- length(Zscores)
  Likelihood_eps <- 0.01
  LL_conv_test <- log(1 + Likelihood_eps)

  # initial values for the parameters:
  mu <- mean(Zscores)
  sigma <- sd(Zscores)
  rho <- 0.01

  # calculate the initial log-likelihood
  LL_current <- LL(Zscores, mu, sigma, rho, z_max, z_min)

  repeat {
    LL_old <- LL_current
    mu_old <- mu
    sigma_old <- sigma
    rho_old <- rho

    # update on rho:
    rho <- sum(sapply(Zscores, function(x) {
      Posterior_for_fg(x, mu_old, sigma_old, rho_old, z_max, z_min) / n_motifs
    }))

    # update on sigma
    sigma2 <- sum(sapply(
      Zscores,
      function(x) {
        (1 - Posterior_for_fg(x, mu_old, sigma_old, rho_old, z_max, z_min)) * (x - mu_old)^2
      }
    )) /
      sum(sapply(
        Zscores,
        function(x) {
          (1 - Posterior_for_fg(x, mu_old, sigma_old, rho_old, z_max, z_min))
        }
      ))
    sigma <- sqrt(sigma2)

    # update on mu
    mu <- sum(sapply(
      Zscores,
      function(x) {
        (1 - Posterior_for_fg(x, mu_old, sigma_old, rho_old, z_max, z_min)) * x
      }
    )) /
      sum(sapply(
        Zscores,
        function(x) {
          (1 - Posterior_for_fg(x, mu_old, sigma_old, rho_old, z_max, z_min))
        }
      ))

    # calculate the log-likelihood with the new estimates
    LL_current <- LL(Zscores, mu, sigma, rho, z_max, z_min)

    # iterave until the LL converge i.e. while L[n] / L[n-1] >= 1 + eps
    if (LL_current - LL_old < LL_conv_test) {
      break
    }
  }
  return(c(LL_current, mu, sigma, rho, z_max, z_min))
}


###############################################################################
# MAIN
###############################################################################

input_table <-
  as.matrix(read.table(opt$input,
    h = TRUE,
    as.is = T,
    sep = "\t",
    row.names = 1,
    comment.char = "",
    check.names = FALSE
  ))

# -----------------------------------------------------------------------------

# get the rows for which we did actually fit
input_table_fitted <- input_table[complete.cases(input_table), ]

# get the motifs that were filtered out
filtered_out_motifs <-
  rownames(input_table[!(rownames(input_table) %in% rownames(input_table_fitted)), , drop = FALSE])

# get the rows for which we did not fit
input_table_not_fitted <- input_table[filtered_out_motifs, , drop = FALSE]

# get the number of motifs for which we did fit
n_motifs <- nrow(input_table_fitted)

# -----------------------------------------------------------------------------

standard_Z_columns <- c() # needed to remember the new columns

original_columns <- colnames(input_table_fitted)
for (i in original_columns) {

  # for every Z-score column:
  if (startsWith(i, "zscores")) {
    standard_Z_columns <- c(standard_Z_columns, paste("standard", i, sep = "_"))

    # Optimize the Mixture Model parameters
    scores <- as.numeric(input_table_fitted[, i])
    summary <- EM(scores)
    scores_LL <- summary[1]
    scores_mu <- summary[2]
    scores_sigma <- summary[3]
    scores_rho <- summary[4]
    scores_z_max <- summary[5]
    scores_z_min <- summary[6]

    # standardize the Z-scores according to the MLE of Gaussian parameters
    temp_standard_Zscores <- sapply(scores, function(x) {
      (x - scores_mu) / scores_sigma
    })
    # update the dataframe with the standardized Z-scores
    input_table_fitted <- cbind(input_table_fitted, temp_standard_Zscores)
    colnames(input_table_fitted)[colnames(input_table_fitted) == "temp_standard_Zscores"] <- paste("standard", i, sep = "_")
    # update the dataframe with the p-values from the standard Z-scores
    temp_standard_pvals <- sapply(temp_standard_Zscores, function(x) {
      pnorm(-abs(x))
    })
    input_table_fitted <- cbind(input_table_fitted, temp_standard_pvals)
    colnames(input_table_fitted)[colnames(input_table_fitted) == "temp_standard_pvals"] <- sub("zscores", "pval", i)
    # update the dataframe with Bonferroni-corrected p-values
    temp_pvals_corr <- sapply(temp_standard_pvals, function(x) {
      p.adjust(x, method = "bonferroni", n = n_motifs)
    })
    input_table_fitted <- cbind(input_table_fitted, temp_pvals_corr)
    colnames(input_table_fitted)[colnames(input_table_fitted) == "temp_pvals_corr"] <- sub("zscores", "corr_pval", i)
  }
}

standard_Zscores_subtable <- type.convert(input_table_fitted[, standard_Z_columns])
input_table_fitted <- data.frame(input_table_fitted, check.names = FALSE)

# Calculate the combined per-motif Z-score
input_table_fitted["combined.standard.Zscore"] <- apply(standard_Zscores_subtable, 1, function(x) {
  sqrt(sum(x^2) / length(standard_Z_columns))
})

# Select the max. absolute per-motif Z-score
input_table_fitted["max.abs.standard.Zscore"] <- apply(standard_Zscores_subtable, 1, function(x) {
  max(abs(x))
})

# sort the final table
if (opt$sorting_strategy == "avg") {
  input_table_fitted <-
    input_table_fitted[order(-input_table_fitted["combined.standard.Zscore"]), ]
} else {
  stopifnot(opt$sorting_strategy == "max")
  input_table_fitted <-
    input_table_fitted[order(-input_table_fitted["max.abs.standard.Zscore"]), ]
}

# -----------------------------------------------------------------------------

# if neccessary - append NA rows with not-fitted motifs
input_table_not_fitted <- data.frame(input_table_not_fitted, check.names = FALSE)
if (nrow(input_table_not_fitted) > 0) {
  input_table_not_fitted["combined.standard.Zscore"] <- NA
  input_table_not_fitted["max.abs.standard.Zscore"] <- NA
  for (i in standard_Z_columns) {
    input_table_not_fitted[i] <- NA
  }
  for (i in standard_Z_columns) {
    input_table_not_fitted[sub("standard_zscores", "pval", i)] <- NA
  }
  for (i in standard_Z_columns) {
    input_table_not_fitted[sub("standard_zscores", "corr_pval", i)] <- NA
  }
  output_table <- rbind(input_table_fitted, input_table_not_fitted)
} else {
  output_table <- input_table_fitted
}

# -----------------------------------------------------------------------------

# Save the results to tsv file
write.table(output_table, file = opt$output, sep = "\t", quote = FALSE)
