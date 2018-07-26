#!/usr/bin/Rscript
#===============================================================================
# selex_multinomial.R
#===============================================================================

# A function for building a multinomial regression model on read counts from
# SELEX results




# Libraries ====================================================================

library(nnet)
library(parallel)




# Constants ====================================================================

default_weights <- c(8, 1, 2, 4, 8)
example_counts <- data.frame(
  a=c(25, 50, 75, 95, 95),
  c=c(25, 20, 8, 2, 2),
  g=c(25, 15, 7, 1, 1),
  t=c(25, 15, 7, 1, 1)
)



# Functions ====================================================================

preprocess_counts <- function(counts, ref_allele = NULL) {
  counts <- data.matrix(counts)
  if (!is.null(ref_allele)) {
    alleles <- colnames(counts)
    alt_alleles <- alleles[alleles != ref_allele]
    counts <- counts[,c(ref_allele, alt_alleles)]
  }
  counts
}

selex_multinom <- function(
  counts,
  weights = default_weights,
  ref_allele = NULL
) {
  counts <- preprocess_counts(counts, ref_allele = ref_allele)
  cycle <- c(0, 1, 2, 3, 4)
  fit <- multinom(counts ~ cycle, weights = weights)
  fit[["counts"]] <- counts
  fit[["ref.allele"]] <- colnames(counts)[[1]]
  fit[["input.weights"]] <- weights
  fit[["coefficients"]] <- summary(fit)[["coefficients"]]
  fit[["standard.errors"]] <- summary(fit)[["standard.errors"]]
  fit[["entropy"]] <- -1 * rowSums(
    fit[["fitted.values"]] * log2(fit[["fitted.values"]])
  )
  fit[["info.content"]] <- 2 - fit[["entropy"]]
  fit
}

contrast <- function(fit, ref_allele, alt_allele) {
  if (ref_allele == fit[["ref.allele"]]) {
    fit[["coefficients"]][alt_allele,]
  } else if (alt_allele == fit[["ref.allele"]]) {
    -1 * fit[["coefficients"]][ref_allele,]
  } else {
    fit[["coefficients"]][alt_allele,] - fit[["coefficients"]][ref_allele,]
  }
}

randomize_counts <- function(counts) {
  t(
    sapply(
      rowSums(counts),
      function(n_reads) rmultinom(1, n_reads, colSums(counts))
    )
  )
}

sample_coefficients <- function(
  counts,
  weights = default_weights,
  n = 1,
  cores = detectCores()
) {
  matrix(
    unlist(
      mclapply(
        lapply(1:n, function(x) randomize_counts(counts)),
        function(c) selex_multinom(c, weights = weights)[["coefficients"]][4:6],
        mc.cores = cores
      )
    ),
    nrow = 3
  )
}

estimate_standard_errors <- function(
  counts,
  weights = default_weights,
  n = 100,
  cores = detectCores()
) {
  sink("/dev/null")
  sample <- sample_coefficients(counts, weights = weights, n = n, cores = cores)
  sink()
  sapply(1:3, function(row) sd(sample[row,]))
}

estimate_z_scores <- function(
  fit,
  estimated_se = NULL,
  n = 100,
  cores = detectCores()
) {
  if (is.null(estimated_se)) {
    estimated_se <- estimate_standard_errors(
      fit[["counts"]],
      weights = fit[["input.weights"]],
      n = n,
      cores = cores
    )
  }
  setNames(fit[["coefficients"]][4:6] / estimated_se, fit[["lab"]][2:4])
}

two_tailed_z_test <- function(z) {
  (1 - pnorm(abs(z))) * 2
}

estimate_pvals <- function(
  fit,
  estimated_se = NULL,
  n = 100,
  cores = detectCores()
) {
  p.adjust(
    two_tailed_z_test(estimate_z_scores(fit, estimated_se, n, cores = cores)),
    method = "bonferroni"
  )
}
