#' Inference on Most and Least Affected Groups
#'
#' \code{subpop} conducts set inference on the groups of most and least
#' affected. When \code{subgroup = NULL}, output is for whole sample. Otherwise
#' the results are subgroup. The output of \code{subpop} is a list
#' containing six components: \code{cs_most}, \code{cs_least}, \code{u},
#' \code{subgroup}, \code{most} and \code{least}. As the names
#' indicate, \code{cs_most} and \code{cs_least} denote the confidence sets for
#' the most and least affected units. \code{u} stores the u-th most and least
#' affected index. \code{subgroup} stores the indicators for subpopulations.
#' \code{most} and \code{least} store the data of the most and
#' least affected groups. The confidence sets can be visualized using the
#' \code{\link{plot.subpop}} command while the two groups can be tabulated via
#' the \code{\link{summary.subpop}} command.
#'
#' @param fm          Regression formula
#' @param data        The data in use
#' @param method      Models to be used for estimating partial effects. Four
#'                    options: \code{"logit"} (binary response),
#'                    \code{"probit"} (binary response), \code{"ols"}
#'                    (interactive linear with additive errors), \code{"QR"}
#'                    (linear model with non-additive errors). Default is
#'                    \code{"ols"}.
#' @param var_type    The type of parameter in interest. Three options:
#'                    \code{"binary"}, \code{"categorical"},
#'                    \code{"continuous"}. Default is \code{"binary"}.
#' @param var         Variable T in interset. Should be a character.
#' @param compare     If parameter in interest is categorical, then user needs
#'                    to specify which two category to compare with. Should be
#'                    a 1 by 2 character vector. For example, if the two levels
#'                    to compare with is 1 and 3, then \code{c=("1", "3")},
#'                    which will calculate partial effect from 1 to 3. To use
#'                    this option, users first need to specify \code{var} as a
#'                    factor variable.
#' @param subgroup    Subgroup in interest. Default is \code{NULL}.
#'                    Specifcation should be a logical variable. For example,
#'                    suppose data contains indicator variable for women
#'                    (female if 1, male if 0). If users are interested in
#'                    women SPE, then users should specify
#'                    \code{subgroup = data[, "female"] == 1}.
#' @param samp_weight Sampling weight of data. Input should be a n by 1 vector,
#'                    where n denotes sample size. Default is \code{NULL}.
#' @param taus        Indexes for quantile regression.
#'                    Default is \code{c(5:95)/100}.
#' @param u           Percentile of most and least affected. Default is set to
#'                    be 0.1.
#' @param alpha       Size for confidence interval. Shoule be between 0 and 1.
#'                    Default is 0.1
#' @param b           Number of bootstrap draws. Default is set to be 500.
#' @param parallel    Whether the user wants to use parallel computation.
#'                    The default is \code{FALSE} and only 1 CPU will be used.
#'                    The other option is \code{TRUE}, and user can specify
#'                    the number of CPUs in the \code{ncores} option.
#' @param ncores      Number of cores for computation. Default is set to be
#'                    \code{detectCores()}, which is a function from package
#'                    \code{parallel} that detects the number of CPUs on the
#'                    current host. For large dataset, parallel computing is
#'                    highly recommended since bootstrap is time-consuming.
#' @param seed        Pseudo-number generation for reproduction. Default is 1.
#' @param boot_type   Type of bootstrap. Default is \code{"nonpar"}, and the
#'                    package implements nonparametric bootstrap. The
#'                    alternative is \code{"weighted"}, and the package
#'                    implements weighted bootstrap.
#'
#' @examples
#' data("mortgage")
#' ### Regression Specification
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec +
#'    ltv_med + ltv_high + denpmi + selfemp + single + hischl
#' ### Issue the subpop command
#' set_b <- subpop(fm, data = mortgage, method = "logit", var = "black",
#' u = 0.1, alpha = 0.1, b = 50)
#'
#' @importFrom Hmisc wtd.quantile
#' @importFrom boot boot
#' @importFrom stats quantile rexp qnorm
#' @importFrom parallel detectCores
#' @importFrom pbapply setpb startpb closepb
#' @export
subpop <- function(fm, data, method = c("ols", "logit", "probit", "QR"),
                   var_type = c("binary", "continuous", "categorical"),
                   var, compare, subgroup = NULL, samp_weight = NULL,
                   taus = c(5:95)/100, u = 0.1, alpha = 0.1, b = 500, seed = 1,
                   parallel = FALSE, ncores = detectCores(),
                   boot_type = c("nonpar", "weighted")) {
  #---------------------- Stopping Conditions ---------------------------
  if (alpha >= 1 || alpha <= 0) stop("Please specify a correct size for
                                     hypothesis testing between 0 and 1.")
  if (u >= 1 || u <= 0) stop("Please provide a group classification
                             quantile between 0 and 1.")
  # ------ Replace Null samp_weight specification
  if (is.null(samp_weight)) samp_weight <- rep(1, dim(data)[1])
  # ------ Matching Arguments
  method <- match.arg(method)
  var_type <- match.arg(var_type)
  boot_type <- match.arg(boot_type)
  # ------ 1. Call to estimate PE and PEsub
  output <- suppressWarnings(peestimate(fm, data, samp_weight, var_type, var, compare, method,
                       subgroup, taus))
  pe_est <- output$pe_est
  # ----- 2. u-most and least Affected Groups
  if (method != "QR") {
    if (is.null(subgroup)) {
      # Threshold Values for u-most/least Affected for WHOLE sample
      effect_high <- wtd.quantile(pe_est, samp_weight, 1 - u)
      effect_low <- wtd.quantile(pe_est, samp_weight, u)
      # Get weight for the two groups. If weight is NULL, it does nothing
      # Two affected groups
      high_affected <- data[pe_est >= effect_high, ]
      low_affected <- data[pe_est <= effect_low, ]
    } else {
      # SUB sample
      pesub_est <- output$pesub_est
      pesub_w <- output$samp_weight_sub
      # Threshold Values for u-most/least Affected
      effect_sub_high <- wtd.quantile(pesub_est, pesub_w, 1 - u)
      effect_sub_low <- wtd.quantile(pesub_est, pesub_w, u)
      subdata <- data[subgroup, ]
      high_affected <- subdata[pesub_est >= effect_sub_high, ]
      low_affected <- subdata[pesub_est <= effect_sub_low, ]
    }
  } # QR requires special treatment due to stacking of quantile indices
  if (method == "QR") {
    if (is.null(subgroup)) {
      # This is where QR is different.
      effect_high <- wtd.quantile(pe_est, matrix(samp_weight, ncol = 1,
                                                 nrow = nrow(pe_est),
                                                 byrow = FALSE), 1 - u)
      effect_low <- wtd.quantile(pe_est, matrix(samp_weight, ncol = 1,
                                                nrow = nrow(pe_est),
                                                byrow = FALSE), u)
      # Kronecker function doesn't work for all data types, so I wrote this
      # alternative
      mesh <- data[rep(1:nrow(data), times = length(taus)), ]
      # Two Affected Groups
      high_affected <- mesh[pe_est >= effect_high, ]
      low_affected <- mesh[pe_est <= effect_low, ]
    } else {# Subsample: use PEsub_est
      pesub_est <- output$pesub_est
      pesub_w <- output$samp_weight_sub
      # Threshold Values for u-most/least Affected
      effect_sub_high <- wtd.quantile(pesub_est, pesub_w, 1 - u)
      effect_sub_low <- wtd.quantile(pesub_est, pesub_w, u)
      subdata <- data[subgroup, ]
      mesh <- subdata[rep(1:nrow(subdata), times = length(taus)), ]
      # Two Affected Groups
      high_affected <- mesh[pesub_est >= effect_sub_high, ]
      low_affected <- mesh[pesub_est <= effect_sub_low, ]
    }
  }
  # ---------------------------- 3. Bootstrap Samples -------------------------
  # set a bootstrap counting variable for the purpose of showing a progress bar
  rep_count <- 1
  # Resampling function for weighted bootstrap
  data_rg <- function(data, mle){
    n <- dim(data)[1]
    # Exponential weights
    multipliers  <- rexp(n)
    # Sampling weight of data
    weight <- samp_weight * multipliers / sum(multipliers) * 20000
    data$.w <- weight
    return(data)
  }
  # Resampling function for nonparametric bootstrap
  data_non <- function(data, mle){
    n <- dim(data)[1]
    multipliers <- as.vector(table(factor(sample(n,n,replace = T),
                                          levels = c(1:n))))
    # Sampling weight of data.bs
    weight <- (multipliers/sum(multipliers)) * samp_weight * 20000
    data$.w <- weight
    return(data)
  }
  # Function that computes bootstrap statistics in each draw
  boot_stat_weight <- function(data){
    # set up a progress bar to document the bootstrap progress
    setpb(pb, rep_count)
    rep_count <<- rep_count + 1
    out_bs <- suppressWarnings(peestimate(fm, data, samp_weight = data$.w, var_type, var,
                         compare, method, subgroup, taus))
    pe_est_bs <- out_bs$pe_est
    if (method == "QR") {
      effect_high_bs <- wtd.quantile(pe_est_bs, matrix(data$.w, ncol = 1,
                                                       nrow = nrow(pe_est),
                                                       byrow = FALSE), 1 - u)
      effect_low_bs <- wtd.quantile(pe_est_bs, matrix(data$.w, ncol = 1,
                                                      nrow = nrow(pe_est),
                                                      byrow = FALSE), u)
    } else {
      effect_high_bs <- wtd.quantile(pe_est_bs, data$.w, 1 - u)
      effect_low_bs <- wtd.quantile(pe_est_bs, data$.w, u)
    }
    if (!is.null(subgroup)) {
      pesub_est_bs <- out_bs$pesub_est
      pesub_w_bs <- out_bs$samp_weight_sub
      effect_sub_high_bs <- wtd.quantile(pesub_est_bs, pesub_w_bs, 1 - u)
      effect_sub_low_bs <- wtd.quantile(pesub_est_bs, pesub_w_bs, u)
      return(c(effect_sub_high_bs, effect_sub_low_bs, pesub_est_bs))
    } else {
      return(c(effect_high_bs, effect_low_bs, pe_est_bs))
    }
  }
  # ----- 4. Conduct Bootstrap
  set.seed(seed)
  if (parallel == FALSE) ncores <- 1
  if (boot_type == "nonpar") {
    data$.w <- samp_weight
    cat(paste("Using", ncores, "CPUs now.\n"))
    pb <- startpb(min = 0, max = b)
    result_boot <- boot(data = data, statistic = boot_stat_weight,
                        sim = "parametric", ran.gen = data_non, mle = 0,
                        parallel = "multicore", ncpus = ncores, R = b)
    closepb(pb)
    data$.w <- NULL
  } else if (boot_type == "weighted") {
    data$.w <- samp_weight
    cat(paste("Using", ncores, "CPUs now.\n"))
    pb <- startpb(min = 0, max = b)
    result_boot <- boot(data = data, statistic = boot_stat_weight,
                        sim = "parametric", ran.gen = data_rg, mle = 0,
                        parallel = "multicore", ncpus = ncores, R = b)
    closepb(pb)
    data$.w <- NULL
  }

  # -----5. Analysis for the subpopulation -=
  if (is.null(subgroup)) {
    # (a) PE
    # MOST Affected Group confidence set
    is_he <- submost(pe_est, result_boot$t[, 3:(length(pe_est) + 2)],
                     result_boot$t[, 1], effect_high, alpha, b)
    # LEAST Affected Group confidence set
    is_le <- subleast(pe_est, result_boot$t[, 3:(length(pe_est) + 2)],
                      result_boot$t[, 2], effect_low, alpha, b)
    output <- list(cs_most = is_he, cs_least = is_le, u = u,
                   subgroup = subgroup, most = high_affected,
                   least = low_affected)
  } else{
    # (b) PEsub
    is_he_sub <- submost(pesub_est, result_boot$t[, 3:(length(pesub_est) + 2)],
                         result_boot$t[, 1], effect_sub_high, alpha, b)
    is_le_sub <- subleast(pesub_est, result_boot$t[, 3:(length(pesub_est) + 2)],
                          result_boot$t[, 2], effect_sub_low, alpha, b)
    output <- list(cs_most = is_he_sub, cs_least = is_le_sub, u = u,
                   subgroup = subgroup, most = high_affected,
                   least = low_affected)
  }
  # claim output as a class
  output <- structure(output, class = "subpop")
  return(output)
}

# -----Two Auxiliary Functions
# Implementing algorithm: Output is confidence set indicator
# Most affected group
submost <- function(est_pe, bs_pe, bs_u, est_u, alpha, b){
  # Find the min (to implement the sup condition as in the paper)
  is_0 <- rank(abs(est_pe - est_u)) == min(rank(abs(est_pe - est_u)))
  draws <- bs_pe - matrix(bs_u, nrow = b, ncol = length(est_pe)) -
    matrix(est_pe - est_u, nrow = b, ncol = length(est_pe), byrow = TRUE)
  bse <- (apply(draws, 2, quantile, .75, na.rm = TRUE) -
            apply(draws, 2, quantile, .25, na.rm = TRUE))/(qnorm(0.75) - qnorm(.25))
  #bm <- apply(draws, 2, mean) # bias estimator
  zs <- apply(-draws[, is_0] / matrix(bse[is_0], nrow = b,
                                      ncol = length(bse[is_0]), byrow = TRUE),
              1, max, na.rm = TRUE)
  crt <- quantile(zs, 1 - alpha)  #critical value
  is_he <- -(est_pe - est_u) / bse <= crt
  return(is_he)
}
# Least affected group
subleast <- function(est_pe, bs_pe, bs_u, est_u, alpha, b) {
  # Find the min (to implement the sup condition as in the paper)
  is_0 <- rank(abs(est_pe - est_u)) == min(rank(abs(est_pe - est_u)))
  draws <- bs_pe - matrix(bs_u, nrow = b, ncol = length(est_pe)) -
    matrix(est_pe - est_u, nrow = b, ncol = length(est_pe), byrow = TRUE)
  bse <- (apply(draws, 2, quantile, .75, na.rm = TRUE) -
            apply(draws, 2, quantile, .25, na.rm = TRUE))/(qnorm(0.75) - qnorm(.25))
  #bm <- apply(draws, 2, mean) # bias estimator
  zs <- apply(draws[, is_0] / matrix(bse[is_0], nrow = b,
                                     ncol = length(bse[is_0]), byrow = TRUE),
              1, max, na.rm = TRUE)
  crt <- quantile(zs, 1 - alpha)  #critical value
  is_le <- (est_pe - est_u) / bse <= crt
  return(is_le)
}

# ----- Summary of subpop output
#' Return the output of \code{\link{subpop}} function.
#'
#' The \code{\link{subpop}} function stores the most and least affected groups.
#' This command allows users to see these two groups and their corresponding
#' characteristics. The command also allows users to check the summary
#' statistics of variables in interest, which can be useful for plotting the
#' projections plot via the \code{\link{plot.subpop}} method.
#'
#' @param object    Output of \code{\link{subpop}} command.
#' @param vars      The variables that users want to see the summary
#'                  statistics. The default is \code{NULL} and the command
#'                  shows all variables. The summary statistcs include min
#'                  1st quartile, median, mean, 3rd quartile and the max.
#' @param ...       additional arguments affecting the summary produced.
#'
#' @examples
#' data("mortgage")
#' ### Regression Specification
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec +
#'    ltv_med + ltv_high + denpmi + selfemp + single + hischl
#' ### Issue the subpop command
#' set_b <- subpop(fm, data = mortgage, method = "logit", var = "black",
#' u = 0.1, alpha = 0.1, b = 50)
#' ### Produce summary of two variables
#' groups <- summary(set_b, vars = c("p_irat", "hse_inc"))
#' @export
summary.subpop <- function(object, vars = NULL, ...) {
  # Retrieve the two affected groups from the subpop object
  high_affected <- object$most
  low_affected <- object$least
  # report all variables if users don't specify a particular interest
  if (is.null(vars)) vars <- colnames(high_affected)
  # an axiliary function that reports summary statistics
  sumstats <- function(x) {
    min <- min(x, na.rm = TRUE)
    max <- max(x, na.rm = TRUE)
    median <- median(x, na.rm = TRUE)
    mean <- mean(x, na.rm = TRUE)
    first_quartile <- quantile(x, 0.25, na.rm = TRUE)
    third_quartile <- quantile(x, 0.75, na.rm = TRUE)
    output <- matrix(rbind(min, first_quartile, median, mean, third_quartile, max))
    #rownames(output) <- c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max")
    return(output)
  }
  info_m <- matrix(as.numeric(unlist(high_affected[,vars])), ncol = length(vars))
  info_l <- matrix(as.numeric(unlist(low_affected[,vars])), ncol = length(vars))
  info_m <- apply(info_m, 2, sumstats)
  info_l <- apply(info_l, 2, sumstats)
  rownames(info_m) <- rownames(info_l) <- c("Min", "1st Quartile", "Median",
                                            "Mean", "3rd Quartile", "Max")
  colnames(info_m) <- colnames(info_l) <- vars
  output <- list(most_affected = high_affected, least_affected = low_affected,
                 stats_most = info_m, stats_least = info_l)
  return(output)
}

# ----- Plotting (2-dimensional projection plots of two specified variables)
#'
#' Plot 2-dimensional projections of variables in interest.
#'
#' Takes output from \code{\link{subpop}} command as inputs and plots
#' 2-dimensional projection plots of two specified variables. If a
#' variable in interest is of type factor, then the user must put it on
#' the y-axis. If the variable on the y-coordinate is a factor, range of
#' y-axis is set to be the factor level. Otherwise, users can use
#' \code{\link{summary.subpop}} to know the ranges of variables in the
#' two groups.
#'
#' @param x           Output of \code{\link{subpop}} command.
#' @param varx        The name of the variable to be plotted on the x-axis.
#' @param vary        The name of the variable name to be plotted on the
#'                    y-axis.
#' @param xlim        The range of x-axis. Default is \code{NULL}.
#' @param ylim        The range of y-axis. Default is \code{NULL}. If the
#'                    variable on the y-coordinate is a factor, the default
#'                    will set it to be the factor level, and users don't
#'                    need to specify \code{ylim}.
#' @param main        Main title of the plot. Default is \code{NULL}.
#' @param sub         Sub title of the plot. Default is NULL.
#' @param xlab        x-axis label. Default is \code{NULL}.
#' @param ylab        y-axis label. Default is \code{NULL}.
#' @param overlap     Whether user wants to allow observations included in both
#'                    confidence sets. Default is \code{FALSE}, and the plot
#'                    drops the overlapped observations.
#' @param ...         Graphics parameters to be passed to the plotting
#'                    routines.
#'
#' @examples
#' data("mortgage")
#' ### Regression Specification
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec +
#'    ltv_med + ltv_high + denpmi + selfemp + single + hischl
#' ### Issue the subpop command
#' set_b <- subpop(fm, data = mortgage, method = "logit", var = "black",
#' u = 0.1, alpha = 0.1, b = 50)
#' ### Plotting
#' plot(set_b, varx = mortgage$p_irat, vary = mortgage$hse_inc,
#'      xlim = c(0, 1.5), ylim = c(0, 1.5), xlab = "Debt/Income",
#'      ylab = "Housing expenses/Income", overlap = TRUE)
#' @importFrom graphics plot polygon lines legend
#' @export
plot.subpop <- function(x, varx, vary, xlim = NULL, ylim = NULL, main = NULL,
                        sub = NULL, xlab = NULL, ylab = NULL, overlap = FALSE,
                        ...) {
  if (is.factor(varx)) stop("Variable on the x-axis cannot be a factor.")
  is_he <- x$cs_most
  is_le <- x$cs_least
  subgroup <- x$subgroup
  u <- x$u
  if (!is.factor(vary)) {
    plot(varx, vary, type = "n", xlim = xlim, ylim = ylim,
         log = "", main, sub, xlab, ylab, col = 4, pch = 20, lwd = 2)
  } else {
    plot(varx, vary, type = "n", xlim = xlim, ylim = NULL,
         log = "", main, sub, xlab, ylab, col = 4, pch = 20, lwd = 2,
         yaxt = 'n')
    axis(side = 2, at = c(1:nlevels(vary)),
         labels = substr(levels(vary), 1, 3))
  }
  if (is.null(subgroup)) {
    # Full sample
    if (overlap == TRUE) {
      points(varx[is_he], vary[is_he], col = 4, pch = 20, lwd = 2)
      points(varx[is_le], vary[is_le], col = 'lightblue1', pch = 1, lwd = 2)
    } else {
      points(varx[is_he & !is_le], vary[is_he & !is_le], col = 4,
             pch = 20, lwd = 2)
      points(varx[is_le & !is_he], vary[is_le & !is_he], col = 'lightblue1',
             pch = 1, lwd = 2)
    }
  } else {
    # Sub sample
    subx <- varx[subgroup]
    suby <- vary[subgroup]
    if (overlap == TRUE) {
      points(subx[is_he], suby[is_he], col = 4, pch = 20, lwd = 2)
      points(subx[is_le], suby[is_le], col = 'lightblue1', pch = 1, lwd = 2)
    } else {
      points(subx[is_he & !is_le], suby[is_he & !is_le], col = 4,
             pch = 20, lwd = 2)
      points(subx[is_le & !is_he], suby[is_le & !is_he], col = 'lightblue1',
             pch = 1, lwd = 2)
    }
  }
  legend('topleft', c(paste0(u*100, "% Most"), paste0(u*100, "% Least")),
         col = c(4, 'lightblue1'), pch = c(20, 1), horiz = F, bty = 'n')
}
