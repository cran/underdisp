## Optional integration with texreg (screenreg / texreg / htmlreg) for formatted
## regression tables. texreg is in Suggests; its extract() methods are registered
## on load only when texreg is installed, and a registration failure is reported
## without blocking the package from loading. modelsummary users are served by
## the broom tidy()/glance() methods (see broom.R) and need nothing here.
## stargazer dispatches on a hardcoded list of model classes and cannot be
## extended to a new class; we point stargazer users to modelsummary or texreg
## instead.

.extract_underdisp <- function(names, coef, se, ll, df, nobs,
                               extra_names = character(0), extra = numeric(0),
                               extra_dec = logical(0), model.name = character(0)) {
  texreg::createTexreg(
    coef.names  = names, coef = coef, se = se,
    pvalues     = 2 * stats::pnorm(-abs(coef / se)),
    gof.names   = c("AIC", "Log Likelihood", "Num. obs.", extra_names),
    gof         = c(-2 * ll + 2 * df, ll, nobs, extra),
    gof.decimal = c(TRUE, TRUE, FALSE, extra_dec),
    model.name  = model.name)
}

## The dispersion rows of a table: the parameter under its own name, with the side of its Poisson value that is
## underdispersion (the families' parameters point different ways; see .disp_info in broom.R); the negative
## binomial reports 1/theta, which is 0 at the Poisson boundary rather than an unreadable theta.
.tex_disp <- function(model) {
  di <- .disp_info(model)
  if (is.na(di$name)) return(list(name = character(0), value = numeric(0), decimal = logical(0)))
  if (!nzchar(di$rule)) return(list(name = "Overdispersion (1/theta)", value = 1 / di$value, decimal = TRUE))
  nm <- strsplit(di$name, ", ", fixed = TRUE)[[1L]]
  list(name = paste0("Dispersion (", nm, "; ", di$rule, ")"), value = di$value, decimal = rep(TRUE, length(nm)))
}

## A table column from the model's tidy() rows (the same standard errors, and the same blank standard errors for
## separated terms) with the fit statistics and the dispersion rows; the column is named by its model and its count
## equation's coefficient scale (a table's custom.model.names replace it)
.tex <- function(model, tab, extra_names = character(0), extra = numeric(0), extra_dec = logical(0)) {
  ll <- stats::logLik(model); d <- .tex_disp(model)
  .extract_underdisp(tab$term, tab$estimate, tab$std.error, as.numeric(ll), attr(ll, "df"), stats::nobs(model),
                     c(d$name, extra_names), c(d$value, extra), c(d$decimal, extra_dec),
                     model.name = paste0(.model_short(model), ": ", .coef_scale(model)$scale))
}

utils::globalVariables(c(".w", ".d"))   # columns built for the participation glm() calls

.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("texreg", quietly = TRUE)) return(invisible())
  tryCatch({
    ## the model classes are S3; texreg dispatches its S4 generic on them, so
    ## they are declared to the methods package in this namespace
    for (cl in c("cpb", "cpb_fe", "hurdle_cpb", "zi_cpb", "gec", "gec_fe", "hurdle_gec",
                 "zi_gec", "count_reg", "hurdle_count", "zi_count"))
      methods::setOldClass(cl, where = asNamespace(pkgname))    # one call each: no inheritance chain
    ## the methods must attach to texreg's generic, not to a generic of that name
    ## looked up from this namespace (which does not import texreg)
    extract <- methods::getGeneric("extract", where = asNamespace("texreg"))
    methods::setMethod(extract, "cpb", function(model, ...) .tex(model, tidy.cpb(model)))
    methods::setMethod(extract, "cpb_fe", function(model, ...)
      .tex(model, tidy.cpb_fe(model), "Num. units", model$n_units, FALSE))
    methods::setMethod(extract, "hurdle_cpb", function(model, ...) .tex(model, tidy.hurdle_cpb(model)))
    methods::setMethod(extract, "zi_cpb", function(model, ...) .tex(model, tidy.zi_cpb(model)))
    methods::setMethod(extract, "gec", function(model, ...) .tex(model, tidy.gec(model)))
    methods::setMethod(extract, "gec_fe", function(model, ...)
      .tex(model, tidy.gec(model), "Num. units", model$n_units, FALSE))
    methods::setMethod(extract, "hurdle_gec", function(model, ...) .tex(model, tidy.hurdle_gec(model)))
    methods::setMethod(extract, "count_reg", function(model, ...) .tex(model, tidy.count_reg(model)))
    methods::setMethod(extract, "hurdle_count", function(model, ...) .tex(model, tidy.hurdle_count(model)))
    methods::setMethod(extract, "zi_count", function(model, ...) .tex(model, tidy.zi_count(model)))
    methods::setMethod(extract, "zi_gec", function(model, ...) .tex(model, tidy.zi_gec(model)))
  }, error = function(e) assign("texreg_error", conditionMessage(e), envir = .ud_state))
}

## a registration failure is reported when the package is attached (a load stays silent)
.ud_state <- new.env(parent = emptyenv())
.onAttach <- function(libname, pkgname) {
  if (!is.null(.ud_state$texreg_error))
    packageStartupMessage("underdisp: the texreg table methods could not be registered (", .ud_state$texreg_error, ").")
}
