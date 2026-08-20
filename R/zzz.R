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

.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("texreg", quietly = TRUE)) return(invisible())
  tryCatch({
    methods::setMethod("extract", "cpb", function(model, ...) {
      se <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      .extract_underdisp(names(model$coefficients), as.numeric(model$coefficients),
                         as.numeric(se), model$loglik, model$df, model$n)
    })
    methods::setMethod("extract", "cpb_fe", function(model, ...) {
      se <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      .extract_underdisp(names(model$coefficients), as.numeric(model$coefficients),
                         as.numeric(se), model$loglik, model$df, model$n,
                         "Num. units", model$n_units, FALSE)
    })
    methods::setMethod("extract", "hurdle_cpb", function(model, ...) {
      ps  <- summary(model$participation)$coefficients
      ie  <- model$intensity$coefficients
      ise <- if (is.null(model$intensity$se.beta)) rep(NA_real_, length(ie)) else model$intensity$se.beta
      idf <- if (!is.null(model$intensity$df)) model$intensity$df else length(ie) + 1L
      df  <- length(stats::coef(model$participation)) + idf
      ll  <- as.numeric(stats::logLik(model$participation)) + model$intensity$loglik
      .extract_underdisp(c(paste0("part: ", rownames(ps)), paste0("int: ", names(ie))),
                         c(ps[, 1], as.numeric(ie)), c(ps[, 2], as.numeric(ise)), ll, df, model$n)
    })
    methods::setMethod("extract", "zi_cpb", function(model, ...) {
      cse <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      zse <- if (is.null(model$se.zero)) rep(NA_real_, length(model$zero_coef))    else model$se.zero
      .extract_underdisp(c(paste0("count: ", names(model$coefficients)),
                           paste0("zero: ",  names(model$zero_coef))),
                         c(as.numeric(model$coefficients), as.numeric(model$zero_coef)),
                         c(as.numeric(cse), as.numeric(zse)), model$loglik, model$df, model$n)
    })
    methods::setMethod("extract", "gec", function(model, ...) {
      se <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      .extract_underdisp(names(model$coefficients), as.numeric(model$coefficients), as.numeric(se),
                         model$loglik, model$df, model$n, "Dispersion", model$delta, TRUE)
    })
    methods::setMethod("extract", "gec_fe", function(model, ...) {
      se <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      .extract_underdisp(names(model$coefficients), as.numeric(model$coefficients), as.numeric(se),
                         model$loglik, model$df, model$n, c("Dispersion", "Num. units"),
                         c(model$delta, model$n_units), c(TRUE, FALSE))
    })
    methods::setMethod("extract", "hurdle_gec", function(model, ...) {
      ps  <- summary(model$participation)$coefficients
      ie  <- model$intensity$coefficients
      ise <- if (is.null(model$intensity$se.beta)) rep(NA_real_, length(ie)) else model$intensity$se.beta
      df  <- length(stats::coef(model$participation)) + model$intensity$df
      ll  <- as.numeric(stats::logLik(model$participation)) + model$intensity$loglik
      .extract_underdisp(c(paste0("part: ", rownames(ps)), paste0("int: ", names(ie))),
                         c(ps[, 1], as.numeric(ie)), c(ps[, 2], as.numeric(ise)), ll, df, model$n,
                         "Dispersion", model$delta, TRUE)
    })
    methods::setMethod("extract", "zi_gec", function(model, ...) {
      cse <- if (is.null(model$se.beta)) rep(NA_real_, length(model$coefficients)) else model$se.beta
      zse <- if (is.null(model$se.zero)) rep(NA_real_, length(model$zero_coef))    else model$se.zero
      .extract_underdisp(c(paste0("count: ", names(model$coefficients)),
                           paste0("zero: ",  names(model$zero_coef))),
                         c(as.numeric(model$coefficients), as.numeric(model$zero_coef)),
                         c(as.numeric(cse), as.numeric(zse)), model$loglik, model$df, model$n,
                         "Dispersion", model$delta, TRUE)
    })
  }, error = function(e) invisible())
}
