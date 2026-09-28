## Optional integration with texreg (screenreg / texreg / htmlreg) for formatted
## regression tables. texreg is in Suggests; its extract() methods are registered
## on load only when texreg is installed, and wrapped so a registration hiccup can
## never block the package from loading. modelsummary users are served by the
## broom tidy()/glance() methods (see broom.R) and need nothing here. stargazer
## dispatches on a hardcoded list of model classes and cannot be extended to a new
## class; we point stargazer users to modelsummary or texreg instead.

.extract_underdisp <- function(names, coef, se, ll, df, nobs,
                               extra_names = character(0), extra = numeric(0),
                               extra_dec = logical(0)) {
  texreg::createTexreg(
    coef.names  = names, coef = coef, se = se,
    pvalues     = 2 * stats::pnorm(-abs(coef / se)),
    gof.names   = c("AIC", "Log Likelihood", "Num. obs.", extra_names),
    gof         = c(-2 * ll + 2 * df, ll, nobs, extra),
    gof.decimal = c(TRUE, TRUE, FALSE, extra_dec))
}

## The dispersion row of a count-family table: the family's shape parameter on
## its natural scale under its own name (the link is stripped from `shape_name`);
## the negative binomial reports 1/theta, which is 0 at the Poisson boundary
## rather than an unreadable theta.
.count_disp_row <- function(fam, theta) {
  if (!fam$nshape) return(list(name = character(0), value = numeric(0), decimal = logical(0)))
  nm <- sub("^[a-z]+\\((.*)\\)$", "\\1", fam$shape_name)
  if (fam$tag == "negbin") list(name = "Overdispersion (1/theta)", value = 1 / theta, decimal = TRUE)
  else list(name = paste0("Dispersion (", nm, ")"), value = theta, decimal = TRUE)
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
    methods::setMethod(extract, "cpb", function(model, ...) {
      se <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      .extract_underdisp(names(model$coefficients), as.numeric(model$coefficients),
                         as.numeric(se), model$loglik, model$df, model$n)
    })
    methods::setMethod(extract, "cpb_fe", function(model, ...) {
      se <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      .extract_underdisp(names(model$coefficients), as.numeric(model$coefficients),
                         as.numeric(se), model$loglik, model$df, model$n,
                         "Num. units", model$n_units, FALSE)
    })
    methods::setMethod(extract, "hurdle_cpb", function(model, ...) {
      ps  <- summary(model$participation)$coefficients
      ie  <- model$intensity$coefficients
      ise <- if (is.null(model$intensity$se.beta)) rep(NA_real_, length(ie)) else model$intensity$se.beta
      idf <- if (!is.null(model$intensity$df)) model$intensity$df else length(ie) + 1L
      df  <- length(stats::coef(model$participation)) + idf
      ll  <- as.numeric(stats::logLik(model$participation)) + model$intensity$loglik
      .extract_underdisp(c(paste0("participation:", rownames(ps)), paste0("intensity:", names(ie))),
                         c(ps[, 1], as.numeric(ie)), c(ps[, 2], as.numeric(ise)), ll, df, model$n)
    })
    methods::setMethod(extract, "zi_cpb", function(model, ...) {
      cse <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      zse <- if (is.null(model$se.zero)) rep(NA_real_, length(model$zero_coef))    else model$se.zero
      .extract_underdisp(c(paste0("count:", names(model$coefficients)),
                           paste0("zero:",  names(model$zero_coef))),
                         c(as.numeric(model$coefficients), as.numeric(model$zero_coef)),
                         c(as.numeric(cse), as.numeric(zse)), model$loglik, model$df, model$n)
    })
    methods::setMethod(extract, "gec", function(model, ...) {
      se <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      .extract_underdisp(names(model$coefficients), as.numeric(model$coefficients), as.numeric(se),
                         model$loglik, model$df, model$n, "Dispersion", model$delta, TRUE)
    })
    methods::setMethod(extract, "gec_fe", function(model, ...) {
      se <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      .extract_underdisp(names(model$coefficients), as.numeric(model$coefficients), as.numeric(se),
                         model$loglik, model$df, model$n, c("Dispersion", "Num. units"),
                         c(model$delta, model$n_units), c(TRUE, FALSE))
    })
    methods::setMethod(extract, "hurdle_gec", function(model, ...) {
      ps  <- summary(model$participation)$coefficients
      ie  <- model$intensity$coefficients
      ise <- if (is.null(model$intensity$se.beta)) rep(NA_real_, length(ie)) else model$intensity$se.beta
      df  <- length(stats::coef(model$participation)) + model$intensity$df
      ll  <- as.numeric(stats::logLik(model$participation)) + model$intensity$loglik
      .extract_underdisp(c(paste0("participation:", rownames(ps)), paste0("intensity:", names(ie))),
                         c(ps[, 1], as.numeric(ie)), c(ps[, 2], as.numeric(ise)), ll, df, model$n,
                         "Dispersion", model$delta, TRUE)
    })
    methods::setMethod(extract, "count_reg", function(model, ...) {
      se <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      fam <- .count_fam(model$family); disp <- .count_disp_row(fam, model$theta)
      .extract_underdisp(names(model$coefficients), as.numeric(model$coefficients), as.numeric(se),
                         model$loglik, model$df, model$n, disp$name, disp$value, disp$decimal)
    })
    methods::setMethod(extract, "hurdle_count", function(model, ...) {
      ps  <- summary(model$participation)$coefficients
      ie  <- model$intensity$coefficients
      ise <- if (is.null(model$intensity$se.beta)) rep(NA_real_, length(ie)) else model$intensity$se.beta
      disp <- .count_disp_row(.count_fam(model$family), model$theta)
      .extract_underdisp(c(paste0("participation:", rownames(ps)), paste0("intensity:", names(ie))),
                         c(ps[, 1], as.numeric(ie)), c(ps[, 2], as.numeric(ise)), model$loglik, model$df, model$n,
                         disp$name, disp$value, disp$decimal)
    })
    methods::setMethod(extract, "zi_count", function(model, ...) {
      cse <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      zse <- if (is.null(model$se.zero)) rep(NA_real_, length(model$zero.coefficients)) else model$se.zero
      disp <- .count_disp_row(.count_fam(model$family), model$theta)
      .extract_underdisp(c(paste0("count:", names(model$coefficients)),
                           paste0("zero:",  names(model$zero.coefficients))),
                         c(as.numeric(model$coefficients), as.numeric(model$zero.coefficients)),
                         c(as.numeric(cse), as.numeric(zse)), model$loglik, model$df, model$n,
                         disp$name, disp$value, disp$decimal)
    })
    methods::setMethod(extract, "zi_gec", function(model, ...) {
      cse <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      zse <- if (is.null(model$se.zero)) rep(NA_real_, length(model$zero_coef))    else model$se.zero
      .extract_underdisp(c(paste0("count:", names(model$coefficients)),
                           paste0("zero:",  names(model$zero_coef))),
                         c(as.numeric(model$coefficients), as.numeric(model$zero_coef)),
                         c(as.numeric(cse), as.numeric(zse)), model$loglik, model$df, model$n,
                         "Dispersion", model$delta, TRUE)
    })
  }, error = function(e) invisible())
}
