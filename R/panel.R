## ---------------------------------------------------------------------------
## Model-agnostic panel tools that compose with any estimator in the package
## (cpb / gec / count_reg and their FE variants). Fixed effects (fe / two-way fe)
## and robust/clustered standard errors are built into the estimators; this file
## adds the Mundlak / correlated-random-effects device.
## ---------------------------------------------------------------------------

#' Mundlak (correlated random effects) device
#'
#' Augments a data frame with the unit-level means of the time-varying numeric
#' covariates in `formula`, and returns the augmented formula and data. Fitting
#' any of the package's estimators on the result implements the Mundlak /
#' correlated-random-effects specification: the coefficients on the original
#' covariates recover the within-unit (fixed-effects-consistent) effects, while
#' the coefficients on the unit means capture (and test) the correlation between
#' the covariates and the unit effect. It is a lighter, between-variation-preserving
#' alternative to full fixed effects, and it is model-agnostic --- the same device
#' feeds [cpb()], [gec()], or [count_reg()].
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param unit Column name identifying the panel unit.
#' @param suffix Suffix for the added unit-mean columns (default `"_mean"`).
#' @return A list with `formula` (the augmented formula), `data` (the augmented
#'   data frame), and `added` (the names of the unit-mean columns).
#' @examples
#' set.seed(1)
#' d <- data.frame(y = rpois(200, 3), x = rnorm(200), unit = factor(rep(1:20, each = 10)))
#' m <- mundlak(y ~ x, d, unit = "unit")
#' fit <- count_reg(m$formula, m$data, family = "poisson")
#' @seealso [count_reg()], [cpb()]
#' @export
mundlak <- function(formula, data, unit, suffix = "_mean") {
  if (!unit %in% names(data)) stop("'unit' must be a column in 'data'.")
  uf <- data[[unit]]
  rhs <- all.vars(formula[-2L])
  added <- character(0); aug <- data
  for (v in intersect(rhs, names(data))) {
    if (is.numeric(data[[v]]) && v != unit) {
      nm <- paste0(v, suffix)
      if (nm %in% names(data))
        stop("column '", nm, "' already exists in 'data'; choose a different 'suffix'.")
      aug[[nm]] <- stats::ave(data[[v]], uf, FUN = function(z) mean(z, na.rm = TRUE))
      added <- c(added, nm)
    }
  }
  if (!length(added)) warning("No time-varying numeric covariates found to average.")
  list(formula = stats::reformulate(c(labels(stats::terms(formula)), added),
                                    response = all.vars(formula)[1L]),
       data = aug, added = added)
}
