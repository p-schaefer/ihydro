#' @noRd
.rem_na_helper <- function(..., rem_one = FALSE) {
  args <- list(...)
  nms <- match.call(expand.dots = FALSE)$`...`
  names(args) <- nms
  out <- tibble::as_tibble(args)
  out <- out[!apply(out, 1, function(row) any(is.na(row))), , drop = FALSE]
  if (rem_one) {
    out <- out[out$cnt != 1, ]
  }
  as.list(out)
}

#' Combine unweighted means across chunks (variance input)
#' @noRd
combine_lumped_mean <- function(sub_mean, sub_n, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_n, cnt)
  sub_mean <- cln$sub_mean
  sub_n <- cln$sub_n
  sum(sub_mean * sub_n) / sum(sub_n)
}

#' Combine unweighted population variances across chunks (denominator N, variance input)
#' @noRd
combine_lumped_var <- function(sub_mean, sub_var, sub_n, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_var, sub_n, cnt, rem_one = TRUE)
  sub_mean <- cln$sub_mean
  sub_var <- cln$sub_var
  sub_n <- cln$sub_n

  mu_pooled <- sum(sub_mean * sub_n) / sum(sub_n)
  sum(sub_n * (sub_var + (sub_mean - mu_pooled)^2)) / sum(sub_n)
}

#' Combine unweighted population SD across chunks (denominator N, variance input)
#' @noRd
combine_lumped_sd <- function(sub_mean, sub_var, sub_n, cnt) {
  suppressWarnings(sqrt(combine_lumped_var(sub_mean, sub_var, sub_n, cnt)))
}

#' Combine unweighted sample variances across chunks (denominator N-1, variance input)
#' @noRd
combine_lumped_sample_var <- function(sub_mean, sub_var, sub_n, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_var, sub_n, cnt, rem_one = TRUE)
  sub_mean <- cln$sub_mean
  sub_var <- cln$sub_var
  sub_n <- cln$sub_n

  mu_pooled <- sum(sub_mean * sub_n) / sum(sub_n)
  # unbiased sample variance for each chunk
  sub_sample_var <- sub_var * sub_n / (sub_n - 1)
  numerator <- sum((sub_n - 1) * sub_sample_var) +
    sum(sub_n * (sub_mean - mu_pooled)^2)
  denominator <- sum(sub_n) - 1
  numerator / denominator
}

#' Combine unweighted sample SD across chunks (denominator N-1, variance input)
#' @noRd
combine_lumped_sample_sd <- function(sub_mean, sub_var, sub_n, cnt) {
  suppressWarnings(sqrt(combine_lumped_sample_var(sub_mean, sub_var, sub_n, cnt)))
}

#' Combine weighted means across chunks (variance input)
#' @noRd
combine_weighted_mean <- function(sub_mean, sub_wt, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_wt, cnt)
  sub_mean <- cln$sub_mean
  sub_wt <- cln$sub_wt
  sum(sub_mean * sub_wt) / sum(sub_wt)
}

#' Combine weighted population variances across chunks (denominator W, variance input)
#' @noRd
combine_weighted_var <- function(sub_mean, sub_var, sub_wt, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_var, sub_wt, cnt, rem_one = TRUE)
  sub_mean <- cln$sub_mean
  sub_var <- cln$sub_var
  sub_wt <- cln$sub_wt

  mu_pooled <- sum(sub_mean * sub_wt) / sum(sub_wt)
  out <- sum(sub_wt * (sub_var + (sub_mean - mu_pooled)^2)) / sum(sub_wt)

  if (is.na(out) || all(out == 0)) {
    out <- NA_real_
  }
  return(out)
}

#' Combine weighted population SD across chunks (denominator W, variance input)
#' @noRd
combine_weighted_sd <- function(sub_mean, sub_var, sub_wt, cnt) {
  suppressWarnings(sqrt(combine_weighted_var(sub_mean, sub_var, sub_wt, cnt)))
}

#' Combine weighted sample variances across chunks (unbiased, denominator W-1, variance input)
#' @noRd
combine_weighted_sample_var <- function(sub_mean, sub_var, sub_wt, cnt) {
  cln <- .rem_na_helper(sub_mean, sub_var, sub_wt, cnt, rem_one = TRUE)
  sub_mean <- cln$sub_mean
  sub_var <- cln$sub_var
  sub_wt <- cln$sub_wt

  mu_pooled <- sum(sub_mean * sub_wt) / sum(sub_wt)
  # unbiased sample variance for each chunk
  sub_sample_var <- sub_var * sub_wt / (sub_wt - 1)
  W <- sum(sub_wt)
  numerator <- sum((sub_wt - 1) * sub_sample_var) +
    sum(sub_wt * (sub_mean - mu_pooled)^2)
  denominator <- W - 1
  out <- numerator / denominator
  if (is.na(out) || all(out == 0)) {
    out <- NA_real_
  }
  return(out)
}

#' Combine weighted sample SD across chunks (unbiased, denominator W-1, variance input)
#' @noRd
combine_weighted_sample_sd <- function(sub_mean, sub_var, sub_wt, cnt) {
  suppressWarnings(sqrt(combine_weighted_sample_var(sub_mean, sub_var, sub_wt, cnt)))
}
