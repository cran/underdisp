## ---------------------------------------------------------------------------
## The unified first-difference contract ("ud_fd"): every first_difference()
## method in the package returns the same object -- a data.frame of class
## c("ud_fd", "data.frame") with columns
##   component | from | to | diff | lower | upper | method
## Single-equation fits return one row (the requested quantity); two-part fits
## return one row per margin LEVEL: the binary-stage probability, the
## count-stage mean, and the marginal mean, each evaluated at `from` and `to`
## (so the rows are quantities, and `stage=` controls WHERE the covariate
## moves: in both equations, or in one holding the other at its reference).
## `method` records the uncertainty source per row: "bootstrap (stored)",
## "bootstrap (refit, B=...)", "delta", or "none" (bounds NA).
## Unknown arguments error loudly (no silent ... swallowing).
## ---------------------------------------------------------------------------

## strict-dots guard: called by every method as .fd_dots(...)
.fd_dots <- function(...) {
  if (...length() > 0L) {
    nm <- names(list(...))
    nm <- if (is.null(nm)) "<positional>" else ifelse(nzchar(nm), nm, "<positional>")
    stop("unused argument(s) for first_difference(): ", paste(nm, collapse = ", "),
         ". See ?first_difference for this model's arguments.", call. = FALSE)
  }
  invisible(NULL)
}

## constructor: assemble + class the unified return
.ud_fd <- function(component, from, to, lower = NA_real_, upper = NA_real_,
                   method = "none") {
  out <- data.frame(component = component, from = from, to = to,
                    diff = to - from,
                    lower = rep_len(lower, length(component)),
                    upper = rep_len(upper, length(component)),
                    method = rep_len(method, length(component)),
                    row.names = NULL, stringsAsFactors = FALSE)
  class(out) <- c("ud_fd", "data.frame")
  out
}

## stage normalization: NULL = all rows; otherwise one canonical component.
## Synonyms accepted across families: participation/zero (binary stage),
## intensity/count (count stage), marginal.
.fd_stage <- function(stage, binary_label) {
  if (is.null(stage)) return(NULL)
  s <- match.arg(stage, c("participation", "zero", "intensity", "count", "marginal"))
  if (s %in% c("participation", "zero")) binary_label
  else if (s %in% c("intensity", "count")) "count-stage"
  else "marginal"
}

#' @export
print.ud_fd <- function(x, digits = 3, ...) {
  y <- as.data.frame(x)
  num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], round, digits)
  if (all(is.na(y$lower))) y$lower <- y$upper <- NULL
  if (all(y$method == "none")) y$method <- NULL
  print.data.frame(y, row.names = FALSE)
  invisible(x)
}

## ---------------------------------------------------------------------------
## The unified rate/odds-ratio contract ("ud_irr"): every irr() method returns
## a data.frame of class c("ud_irr","data.frame") with columns
##   term | equation ("count" or "binary") | ratio ("IRR" or "OR") |
##   estimate | lower | upper | method
## Two-part fits stack both equations; single-equation fits return the count
## rows only. `method` records the interval source per row.
## ---------------------------------------------------------------------------

.ud_irr <- function(term, equation, ratio, estimate,
                    lower = NA_real_, upper = NA_real_, method = "none") {
  out <- data.frame(term = term, equation = rep_len(equation, length(term)),
                    ratio = rep_len(ratio, length(term)),
                    estimate = estimate,
                    lower = rep_len(lower, length(term)),
                    upper = rep_len(upper, length(term)),
                    method = rep_len(method, length(term)),
                    row.names = NULL, stringsAsFactors = FALSE)
  class(out) <- c("ud_irr", "data.frame")
  out
}

## Wald rows on the ratio scale from coefficients + SEs (NULL se -> no interval)
.ud_irr_wald <- function(b, se, level, equation, ratio, method = "wald") {
  if (is.null(se) || !all(is.finite(se))) {
    .ud_irr(names(b), equation, ratio, exp(unname(b)))
  } else {
    z <- stats::qnorm(1 - (1 - level) / 2)
    .ud_irr(names(b), equation, ratio, exp(unname(b)),
            lower = exp(unname(b) - z * se), upper = exp(unname(b) + z * se),
            method = method)
  }
}

#' @export
print.ud_irr <- function(x, digits = 3, ...) {
  y <- as.data.frame(x)
  num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], round, digits)
  if (all(is.na(y$lower))) y$lower <- y$upper <- NULL
  if (all(y$method == "none")) y$method <- NULL
  print.data.frame(y, row.names = FALSE)
  invisible(x)
}

## delta-method interval for a scalar function fdfun(par) with covariance vc;
## returns c(lower, upper) or c(NA, NA)
.fd_delta_ci <- function(fdfun, par, vc, level) {
  if (is.null(vc)) return(c(NA_real_, NA_real_))
  fd <- fdfun(par)
  g  <- tryCatch(numDeriv::grad(fdfun, par), error = function(e) NULL)
  if (is.null(g)) return(c(NA_real_, NA_real_))
  se <- sqrt(max(as.numeric(t(g) %*% vc %*% g), 0))
  z  <- stats::qnorm(1 - (1 - level) / 2)
  c(fd - z * se, fd + z * se)
}
