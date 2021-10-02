# This part of function estimate partial effect and outputs the following:
# (1) Partial Effects
# (2) Partial Effects on the specified subcat group
# (3) Sampling weights for the subcat group

# One caution for compare: the vector should be compare = c("a", "b") where "a"
# and "b" denote the two levels that the users are interested in. Must be
# specified as character otherwise the code breaks down.
# If variable type is categorical, varT should be char factor. For example:
# varT = "occ2" subcat allows users to specify partial effect they are
# interested in. The default is treated group, but can be more general like
# treated group conditional on education, income bracket, etc.
#
#' @importFrom SparseM as.matrix.csr
#' @importFrom quantreg sfn.control rq
#' @importFrom stats lm binomial glm predict
peestimate <- function(fm, data, samp_weight = NULL, var_type, var, compare,
                       method, subcat = NULL, taus = NULL) {
  # --------------------------- Stopping conditions ---------------------------
  if (method == "QR" && is.null(taus)) {
    stop("Please specify quantile indexes for quantile regression.")
  }
  if (var_type != "binary" && var_type != "categorical" && var_type !=
      "continuous") {
    stop("Please specify a correct variable type")
  }
  if (is.null(var) || !is.character(var) || length(var) != 1) {
    stop("Please specify the variable T in interest as a character.")
  }
  if (var_type == "categorical" && !is.factor(data[, var])) {
    stop("Please first declare the categorical variable T as a factor.")
  }
  if (var_type == "categorical" && (is.null(compare) || length(compare) != 2 ||
                                    !is.character(compare))) {
    stop("Please specify which two category levels (char) to compare with.")
  }
  if (!is.null(subcat) && !is.logical(subcat)) {
    stop("Please provide a correct specification of sub category group.")
  }
  # ------Taking Care of Sampling Weights -----------
  # This is to ensure samp_weight can be passed to regression
  environment(fm) <- environment()
  n <- dim(data)[1]
  if (is.null(samp_weight)) samp_weight <- rep(1, n)
  samp_weight_t <- samp_weight[subcat]
  # ------1. Fit a model under a chosen regression method--------------
  if (method == "logit") {
  #  model_fit <- glm(fm, data = data, weights = samp_weight,
  #                   family = binomial(link = "logit"))
    model_fit <- glm(fm, data = data, weights = samp_weight,
                     family = binomial(link = "logit"))
  }
  if (method == "probit") {
    model_fit <- glm(fm, data = data, weights = samp_weight,
                     family = binomial(link = "probit"))
  }
  if (method == "ols") model_fit <- lm(fm, data = data, weights = samp_weight)
  if (method == "QR" && !is.null(taus)) {
    model_fit <- rq(fm, data = data, weights = samp_weight, tau = taus,
                    method = "sfn", control = list(tmpmax = 1e6))
  }
  # ------2. Estimated Partial Effects under a chosen variable T type--------
  d0 <- d1 <- data
  # Three types of variably in interest
  if (var_type == "binary") {
    d0[, var] <- rep(0, n)
    d1[, var] <- rep(1, n)
  } else if (var_type == "categorical") {
    d0[, var] <- rep(compare[1], n)
    d1[, var] <- rep(compare[2], n)
  } else if (var_type == "continuous") {
    ######################################
    # Acknowledgement: This part is inspired by Thomas Leeper's code of
    # `margins` package. His code is available on Github.
    ######################################
    # set a small number
    eps <- 1e-7
    # set value of h based on eps to deal with machine precision
    setstep <- function(x) x + (max(abs(x), 1, na.rm = TRUE) * sqrt(eps)) - x
    d0[, var] <- d0[, var] - setstep(d0[, var])
    d1[, var] <- d1[, var] + setstep(d1[, var])
  } # For binary and cat variable in interest, PE is just discrete difference
  if (var_type == "binary" || var_type == "categorical") {
    if (method == "logit" || method == "probit") {
      pe_est <- predict(model_fit, newdata = d1, type = "response") -
        predict(model_fit, newdata = d0, type = "response")
    } else if (method == "QR" || method == "ols") {
      pe_est <- predict(model_fit, newdata = d1) -
        predict(model_fit, newdata = d0)
    }
  } else if (var_type == "continuous") {
    # Use central difference numerical differentation
    if (method == "logit" || method == "probit") {
      pred <- (predict(model_fit, newdata = d1, type = "response") -
              predict(model_fit, newdata = d0, type = "response")) / (d1 - d0)
      pe_est <- pred[, var]
    } else if (method == "QR" || method == "ols") {
      denominator <- d1 - d0
      denominator <- denominator[[var]]
      pe_est <- (predict(model_fit, newdata = d1) -
                   predict(model_fit, newdata = d0)) / denominator
    }
  }
  # ----- 3. Return Output as a List
  # First set PE_est as matrix for methods other than QR
  # Cannot set QR PE as matrix here: otherwise the which function later can
  # only subset the first index
  if (method != "QR") pe_est <- as.matrix(pe_est)
  # subcat is not null, then output PE_est, PEsub.est and sub sample weight
  if (!is.null(subcat)) {
    if (method != "QR") {
      pet_est <- as.matrix(pe_est[subcat, ])
    } else {
      # This is to take care of the quantile regression (For taus >= 2, output
      # is a matrix, need to conver to a giant col vector for future steps)
      pet_est <- matrix(pe_est[which(subcat), ], ncol = 1, byrow = FALSE)
      # samp_weight_t is repetitive by design
      samp_weight_t <- matrix(samp_weight_t, ncol = 1, nrow = nrow(pet_est),
                              byrow = FALSE)
      pe_est <- matrix(pe_est, ncol = 1, byrow = FALSE)
    }
    out <- list(pe_est = pe_est, pesub_est = pet_est,
                samp_weight_sub = samp_weight_t)
    return(out)
  } else {
    # otherwise just output PE_est
    if (method == "QR") pe_est <- matrix(pe_est, ncol = 1, byrow = FALSE)
    out <- list(pe_est = pe_est)
    return(out)
  }
}
