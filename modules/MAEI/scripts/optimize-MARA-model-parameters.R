###############################################################################
#
#   Optimize MARA model parameters (find maximum likelihood estimates)
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: University_of_Basel
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 20-01-2020
#   LICENSE: GPL_3.0
#   LICENSE: https://www.gnu.org/licenses/gpl-3.0.html
#
###############################################################################

# by default: suppress warnings
options(warn = -1)

# load libraries
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("Rcpp"))

# list the command-line arguments
option_list <- list(
  make_option(c("--cpp_model"),
    action = "store",
    dest = "cpp_model",
    type = "character",
    help = "Path to the cpp file with MARA model functions."
  ),
  make_option(c("--inclusion_table"),
    action = "store",
    dest = "inclusion_table",
    type = "character",
    help = "Filtered and sorted exon inclusion table."
  ),
  make_option(c("--sitecount_matrix"),
    action = "store",
    dest = "sitecount_matrix",
    type = "character",
    help = "Filtered and sorted sitecount matrix."
  ),
  make_option(c("--out1"),
    action = "store",
    dest = "out1",
    type = "character",
    help = "Path for outfile 1 with optimized parameters."
  ),
  make_option(c("--out2"),
    action = "store",
    dest = "out2",
    type = "character",
    help = "Path for outfile 2 with optimized parameters."
  ),
  make_option(c("--min_motif_fraction"),
    action = "store",
    dest = "min_motif_fraction",
    default = 0.00,
    type = "double",
    help = "Minimal fraction of exons which need to bind a specific motif in order for it to be considered."
  ),
  make_option(c("--average_expressions"),
    action = "store",
    dest = "average_expressions",
    default = FALSE,
    type = "logical",
    help = "Average TPM expression across samples."
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

###############################################################################
# MAIN
###############################################################################

# compile cpp functions
# (by default all the files are in the same directory)
sourceCpp(opt$cpp_model)

# load sitecount matrix
N_matrix <-
  as.matrix(read.table(opt$sitecount_matrix,
    h = TRUE,
    as.is = T,
    sep = "\t",
    row.names = 1,
    comment.char = "",
    check.names = FALSE
  ))

# load inclusion table
inclusion_table <-
  as.matrix(read.table(opt$inclusion_table,
    h = TRUE,
    as.is = T,
    sep = "\t",
    row.names = 1,
    comment.char = "",
    check.names = FALSE
  ))

# -----------------------------------------------------------------------------

n_samples <- ncol(inclusion_table) / 2

i_matrix <- inclusion_table[, 1:n_samples, drop = FALSE]
t_matrix <- inclusion_table[, (n_samples + 1):(2 * n_samples), drop = FALSE]

samples <- gsub("_included", "", colnames(i_matrix))

# -----------------------------------------------------------------------------

# the hyperparameter to control the weight of gene expression on the likelihood
# if t_crit == median avg gene expression it means that for the upper half of
# expressed genes the weight of the expression on the likelihood is diminished
#
#
#
if (opt$average_expressions) {
  t_crit <- median(rowMeans(t_matrix))
} else {
  t_crit <- apply(t_matrix, 2, median)
}

# -----------------------------------------------------------------------------

# calculate the 'binding fraction' for each motif
motifs <- colnames(N_matrix)
binding_fraction_values <-
  apply(N_matrix, 2, function(x) {
    (sum(x > 0)) / length(x)
  })
binding_fraction_matrix <-
  matrix(c(motifs, binding_fraction_values), ncol = 2, byrow = FALSE)
colnames(binding_fraction_matrix) <- c("motif", "fraction")
rownames(binding_fraction_matrix) <- binding_fraction_matrix[, "motif"]
binding_fraction_matrix <-
  subset(binding_fraction_matrix, select = c("fraction"), drop = FALSE)

# Filter out motifs that are not present in any exon (always):
cols_to_keep <-
  apply(N_matrix, 2, function(x) {
    sum(x > 0) > 0
  })
N_filtered <- N_matrix[, cols_to_keep, drop = FALSE]

# Filter out motifs that are present below the min fraction (optional):
if (opt$min_motif_fraction != 0.0) {
  cols_to_keep <-
    apply(N_filtered, 2, function(x) {
      (sum(x > 0)) / length(x) >= opt$min_motif_fraction
    })
  N_filtered <- N_filtered[, cols_to_keep, drop = FALSE]
}

# get the motifs that are filtered out
filtered_out_motifs <-
  colnames(N_matrix)[!(colnames(N_matrix) %in% colnames(N_filtered))]

# -----------------------------------------------------------------------------

# Intersect exons of the sitecount matrix and the inclusion table

inclusion_table_exons <- rownames(inclusion_table)
sitecount_matrix_exons <- rownames(N_filtered)
exon_intersection <- intersect(inclusion_table_exons, sitecount_matrix_exons)

intersected_i_matrix <-
  i_matrix[match(exon_intersection, rownames(i_matrix)), , drop = FALSE]
intersected_t_matrix <-
  t_matrix[match(exon_intersection, rownames(t_matrix)), , drop = FALSE]
intersected_N_matrix <-
  N_filtered[match(exon_intersection, rownames(N_filtered)), , drop = FALSE]

# -----------------------------------------------------------------------------

# Expectation-Maximisation optimisation: MARA Model

model_results <- fit_model_parameters(
  intersected_i_matrix,
  intersected_t_matrix,
  intersected_N_matrix,
  t_crit,
  as.integer(opt$average_expressions)
)
A_b_table <- model_results[[1]]
c_table <- model_results[[2]]
c_table <- data.frame(c_table)

# -----------------------------------------------------------------------------

A_b_columns <-
  c(
    paste("A", samples, sep = "_"),
    paste("stdA", samples, sep = "_"),
    paste("b", samples, sep = "_"),
    paste("stdb", samples, sep = "_"), "LL"
  )

# Add the original row/column names from the input matrices:
rownames(A_b_table) <- colnames(intersected_N_matrix)
colnames(A_b_table) <- A_b_columns
rownames(c_table) <- paste("c", rownames(intersected_N_matrix), sep = "_")

# important: update the results table for columns which were not fitted
# if a motif was not fitted activity=0 and we are certain about it => std=0
if (length(filtered_out_motifs) > 0) {
  missing <- data.frame(matrix(ncol = 4 * n_samples + 1, nrow = length(filtered_out_motifs), data = 0))
  colnames(missing) <- A_b_columns
  rownames(missing) <- filtered_out_motifs
  A_b_table <- rbind(A_b_table, missing)
}

# Add a column with the motifs binding fraction information
binding_fraction_matrix <-
  binding_fraction_matrix[row.names(A_b_table), , drop = FALSE]
A_b_table <- cbind(A_b_table, binding_fraction_matrix)

# Save the results to tsv file
write.table(A_b_table,
  file = opt$out1,
  sep = "\t",
  quote = FALSE
)
write.table(c_table,
  file = opt$out2,
  sep = "\t",
  col.names = FALSE,
  quote = FALSE
)
