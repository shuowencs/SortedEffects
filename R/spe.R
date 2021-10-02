#' Empirical Sorted Partial Effects (SPE) and Inference
#'
#' \code{spe} conducts SPE estimation and inference at user-specifed quantile
#' index. The bootstrap procedures follows algorithm 2.1 as in Chernozhukov,
#' Fernandez-Val and Luo (2018). All estimates are bias-corrected and all
#' confidence bands are monotonized. For graphical results, please use
#' \code{\link{plot.spe}}.
#'
#' @return The output is a list with 4 components: (1) \code{spe} stores spe
#' estimates, the upper and lower confidence bounds, and standard errors;
#' (2) \code{ape} stores ape estimates, the upper and lower confidence bounds,
#' and the standard error; (3) \code{us} stores percentile index as in \
#' code{spe} command; (4) \code{alpha} stores significance level as in
#' \code{spe} command.
#'
#' @param fm          Regression formula.
#' @param data        Data in use.
#' @param method      Models to be used for estimating partial effects. Four
#'                    options: \code{"logit"} (binary response),
#'                    \code{"probit"} (binary response), \code{"ols"}
#'                    (interactive linear with additive errors), \code{"QR"}
#'                    (linear model with non-additive errors). Default is
#'                    \code{"ols"}.
#' @param var_type    The type of parameter in interest. Three options:
#'                    \code{"binary"}, \code{"categorical"},
#'                    \code{"continuous"}. Default is \code{"binary"}.
#' @param var         Variable T in interset. Should be a character type.
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
#'                    women SPE, then users should specify \code{subgroup =
#'                    data[, "female"] == 1}.
#' @param samp_weight Sampling weight of data. Input should be a n by 1 vector,
#'                    where n denotes sample size. Default is \code{NULL}.
#' @param us          Percentile of interest for SPE. Should be a vector of
#'                    values between 0 and 1. Default is \code{c(1:9)/10}.
#' @param alpha       Size for confidence interval. Shoule be between 0 and 1.
#'                    Default is 0.1
#' @param taus        Indexes for quantile regression. Default is
#'                    \code{c(5:95)/100}.
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
#' @param bc          Whether want the estimate to be bias-corrected. Default
#'                    is \code{TRUE}. If \code{FALSE} uncorrected estimate and
#'                    corresponding confidence bands will be reported.
#' @param boot_type   Type of bootstrap. Default is \code{"nonpar"}, and the
#'                    package implements nonparametric bootstrap. The other
#'                    alternative is \code{"weighted"}, and the package
#'                    implements weighted bootstrap.
#' @examples
#' data("mortgage")
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec + ltv_med +
#' ltv_high + denpmi + selfemp + single + hischl
#'
#' test <- spe(fm = fm, data = mortgage, var = "black", method = "logit",
#' us = c(2:98)/100, b = 50)
#'
#' @importFrom Hmisc wtd.quantile
#' @importFrom boot boot
#' @importFrom stats quantile rexp qnorm
#' @importFrom parallel detectCores
#' @importFrom pbapply setpb startpb closepb
#' @export
spe <- function(fm, data, method = c("ols", "logit", "probit", "QR"),
                var_type = c("binary", "continuous", "categorical"), var,
                compare, subgroup = NULL, samp_weight = NULL, us = c(1:9)/10,
                alpha = 0.1, taus = c(5:95)/100, b = 500, parallel = FALSE,
                ncores = detectCores(), seed = 1, bc = TRUE,
                boot_type = c("nonpar", "weighted")) {
  # ----- Stopping Condition
  if (alpha >= 1 || alpha <= 0) stop("Please specify a correct size for
                                     hypothesis testing between 0 and 1.")
  # ----- Replace Null samp_weight Specification
  if (is.null(samp_weight)) samp_weight <- rep(1, dim(data)[1])
  samp_weight <- samp_weight/mean(samp_weight) # renormalize
  # ----- Matching Arguments
  method <- match.arg(method)
  var_type <- match.arg(var_type)
  boot_type <- match.arg(boot_type)
  # ----- 1. Call to estimate PE
  output <- suppressWarnings(peestimate(fm, data, samp_weight, var_type, var, compare, method,
                       subgroup, taus))
  pe_est <- output$pe_est
  # ----- 2. Now get estimated SPE (full and subgroup samples)
  if (method != "QR") {
    spe_est <- wtd.quantile(pe_est, samp_weight, us)
  } else {
    spe_est <- wtd.quantile(pe_est, matrix(samp_weight, ncol = 1,
                                           nrow = nrow(pe_est),
                                           byrow = FALSE), us)
  }
  if (!is.null(subgroup)) {
    pesub_est <- output$pesub_est
    pesub_w <- output$samp_weight_sub
    spesub_est <- wtd.quantile(pesub_est, pesub_w, us)
  }
  # -----3. Bootstrap Samples
  #  statistic in one bootstrap for SPE
  #  Accommodate two types of bootstraps:
  #  (1) Nonparametric: draw multinomial weights
  #  (2) Weighted:      draw exponential weights

  # set a bootstrap counting variable for the purpose of showing a progress bar
  rep_count <- 1

  # The resampling function for weighted bootstrap
  data_rg <- function(data, mle) {
    n <- dim(data)[1]
    # Exponential weights
    multipliers  <- rexp(n)
    # Sampling weight of data.bs
    weight <- (multipliers/sum(multipliers)) * samp_weight
    data$.w <- weight/mean(weight)
    return(data)
  }
  # Resampling function for nonparametric bootstrap
  # Note: nonpar is a special type of weighted bootstrap with multinomial weight
  data_non <- function(data, mle) {
    n <- dim(data)[1]
    multipliers <- as.vector(table(factor(sample(n, n, replace = T),
                                          levels = c(1:n))))
    # Sampling weight of data.bs
    weight <- (multipliers/sum(multipliers)) * samp_weight
    # Divide by the mean of weight to stablize bootstrap estimate
    weight <- weight/mean(weight)
    data$.w <- weight
    return(data)
  }
  # Function that computes bootstrap statistics in each draw
  stat_boot_weight <- function(data) {
    setpb(pb, rep_count)
    rep_count <<- rep_count + 1
    output_bs <- suppressWarnings(peestimate(fm, data, samp_weight = data$.w, var_type, var,
                            compare, method, subgroup, taus))
    est_pe_bs <- output_bs$pe_est
    est_ape_bs <- mean(est_pe_bs) # scalar, length = 1
    if (method != "QR") {
      est_spe_bs <- wtd.quantile(est_pe_bs, samp_weight, us)
    } else {
      est_spe_bs <- wtd.quantile(est_pe_bs, matrix(samp_weight, ncol = 1,
                                                   nrow = nrow(pe_est),
                                                   byrow = FALSE), us)
    } # length = length(us)
    if (!is.null(subgroup)) {
      est_pesub_bs <- output_bs$pesub_est
      pesub_w_bs <- output_bs$samp_weight_sub
      est_apesub_bs <- mean(est_pesub_bs) # length = 1
      est_spesub_bs <- wtd.quantile(est_pesub_bs, pesub_w_bs, us) # length(us)
      return(c(est_apesub_bs, est_spesub_bs)) # Each has length: 1+length(us)
    } else {
      return(c(est_ape_bs, est_spe_bs))
    }
  }

  # Use boot command
  set.seed(seed)
  if (parallel == FALSE) ncores <- 1
  if (boot_type == "nonpar") {
    data$.w <- samp_weight
    # print a message showing how many cores are used
    cat(paste("Using", ncores, "CPUs now.\n"))
    # set up a progress bar
    pb <- startpb(min = 0, max = b)
    result_boot <- boot(data = data, statistic = stat_boot_weight,
                        sim = "parametric", ran.gen = data_non, mle = 0,
                        parallel = "multicore", ncpus = ncores, R = b)
    data$.w <- NULL
    closepb(pb)
  } else if (boot_type == "weighted") {
    data$.w <- samp_weight
    cat(paste("Using", ncores, "CPUs now.\n"))
    pb <- startpb(min = 0, max = b)
    result_boot <- boot(data = data, statistic = stat_boot_weight,
                        sim = "parametric", ran.gen = data_rg, mle = 0,
                        parallel = "multicore", ncpus = ncores, R = b)
    data$.w <- NULL
    closepb(pb)
  }
  # ------ 4. Inference
  ##############################
  ### Full Sample
  ##############################
  # SPE inference
  draws_spe_bs <- result_boot$t[, 2:(length(us) + 1)]
  # bundle function implements algorithm 2.1, remark 2.2 and 2.3
  inf_spe <- bundle(draws_spe_bs, spe_est, alpha)
  # APE inference
  est_ape <- mean(pe_est)
  draws_ape_bs <- result_boot$t[, 1]
  inf_ape <- ape(draws_ape_bs, est_ape, alpha)
  ##############################
  ### Subgroup Sample (if subgroup is not NULL)
  ##############################
  if (!is.null(subgroup)) {
    # SPE
    draws_spesub_bs <- result_boot$t[, 2:(length(us) + 1)]
    inf_spesub <- bundle(draws_spesub_bs, spesub_est, alpha)
    # APE
    est_apesub <- mean(pesub_est)
    draws_apesub_bs <- result_boot$t[, 1]
    inf_apesub <- ape(draws_apesub_bs, est_apesub, alpha)
  }
  # ----- 5. Return Results
  # Depends on whether subgroup is NULL and whether bias correction is wanted
  # The resulting outputs are APE(sub) & SPE(sub) with correpsonding confidence
  # intervals. All are stored in an "inf" bundle list
  if (!is.null(subgroup)) {
    if (bc == TRUE) {
      output <- list(spe = inf_spesub[c(1:3, 7)], ape = inf_apesub[c(1:3, 7)],
                     us = us, alpha = alpha)
    } else {
      output <- list(spe = inf_spesub[4:7], ape = inf_apesub[4:7], us = us,
                     alpha = alpha)
    }
  } else {
    if (bc == TRUE) {
      output <- list(spe = inf_spe[c(1:3, 7)], ape = inf_ape[c(1:3, 7)],
                     us = us, alpha = alpha)
    } else {
      output <- list(spe = inf_spe[4:7], ape = inf_ape[4:7], us = us,
                     alpha = alpha)
    }
  }
  # claim output as an object of class "spe" so as to use S3
  output <- structure(output, class = "spe")
  return(output)
}

# ------- Auxiliary Functions ----------------
# Implementing algorithm 2.1 and get (bias corrected) estimate and confidence
# bands for sorted effects
bundle <- function(bs, est, alpha) {
  z_bs <- bs - matrix(est, nrow = nrow(bs), ncol = ncol(bs), byrow = TRUE)
  sigma <- (apply(z_bs, 2, quantile, .75, na.rm = TRUE) -
              apply(z_bs, 2, quantile, .25, na.rm = TRUE)) / (qnorm(0.75) -
                                                              qnorm(.25))
  t_hat <- apply(abs(z_bs / matrix(sigma, nrow = nrow(bs), ncol = ncol(bs),
                                   byrow = TRUE)), 1, max)
  crt <- quantile(t_hat, 1 - alpha)
  # bias-correction
  est_bc <- sort(2*est - apply(bs, 2, mean))
  ubound_est_bc <- sort(est_bc + crt * sigma)
  lbound_est_bc <- sort(est_bc - crt * sigma)
  # uncorrected
  ubound_est <- sort(est + crt * sigma)
  lbound_est <- sort(est - crt * sigma)
  out <- list(est_bc = est_bc, ubound_est_bc = ubound_est_bc,
              lbound_est_bc = lbound_est_bc, est = est, ubound_est = ubound_est,
              lbound_est = lbound_est, sigma_spe = sigma)
  return(out)
}

# Get (biased corrected) average partial effects and confidence intervals
ape <- function(bs, est, alpha) {
  sigma <- (quantile(bs, .75, na.rm = TRUE) -
              quantile(bs, .25, na.rm = TRUE)) / (qnorm(0.75) - qnorm(.25))
  mzs <- abs(bs - est) / sigma
  crt2 <- quantile(mzs, 1 - alpha)
  # bias-correction
  est_bc <- 2 * est - mean(bs)
  ubound_ape_bc <- est_bc + crt2 * sigma
  lbound_ape_bc <- est_bc - crt2 * sigma
  # uncorrected
  ubound_ape <- est + crt2 * sigma
  lbound_ape <- est - crt2 * sigma
  out <- list(est_bc = est_bc, ubound_est_bc = ubound_ape_bc,
              lbound_est_bc = lbound_ape_bc, est = est, ubound_ape = ubound_ape,
              lbound_ape = lbound_ape, sigma_ape = sigma)
  return(out)
}

# ------- Plotting Function ----------------
# Since plot() is a generic, we write a new plotting method for the SPE class.
# --Plotting (Output: average effect, sorted effect, correspondent confidence bands) --
#' Plot output of \code{\link{spe}} command. The x-axis limits are set to the
#' specified range of percentile index.
#' @param x           Output of \code{\link{spe}} command.
#' @param ylim        y-axis limits. Default is NULL.
#' @param main        Main title of the plot. Defualt is NULL.
#' @param sub         Sub title of the plot. Default is NULL.
#' @param xlab        x-axis label. Default is "Percentile Index".
#' @param ylab        y-axis label. Default is "Sorted Effects".
#' @param ...         graphics parameters to be passed to the plotting
#'                    routines.
#' @examples
#' data("mortgage")
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec + ltv_med +
#' ltv_high + denpmi + selfemp + single + hischl
#' test <- spe(fm = fm, data = mortgage, var = "black", method = "logit",
#' us = c(2:98)/100, b = 50)
#'
#' plot(x = test, main="APE and SPE of Being Black on the prob of
#' Mortgage Denial", sub="Logit Model", ylab="Change in Probability")
#'
#' @importFrom graphics plot polygon lines abline points legend
#' @export
plot.spe <- function(x, ylim = NULL, main = NULL, sub = NULL,
                     xlab = "Percentile Index", ylab = "Sorted Effects",
                     ...) {
  # SPE and confidence bands
  SE <- x$spe
  AE <- x$ape
  us <- x$us
  alpha <- x$alpha
  xlim <- range(us)
  plot(us, SE[[1]], type = "l", xlim, ylim, log = "", main, sub, xlab, ylab,
       col = 4, lwd = 2)
  polygon(c(us, rev(us)),c(SE[[2]], rev(SE[[3]])), density = 60, border = F,
          col = 'light blue', lty = 1, lwd = 1)
  lines(us, SE[[1]], lwd = 2, col = 4 )
  # APE and CI
  abline(h = AE[[1]], col = 1)
  abline(h = AE[[2]], col = 1, lty = 2)
  abline(h = AE[[3]], col = 1, lty = 2)
  points(c(min(us),.2,.4,.6,.8,max(us)), rep(AE[[1]], 6), col = 1, pch = 15)
  legend(x = "topleft", col = c(4, 1, "light blue", 1), lwd = c(1, 1, 5, 1),
         lty = c(1, 1, 1, 2), pch = c(NA, 15, NA, NA), pt.cex = c(2, 1),
         bty = 'n', legend = c("SPE","APE", paste0((1 - alpha)*100,"% CB(SPE)"),
                    paste0((1 - alpha)*100,"% CB(APE)")))
}

# ------- Summary Function ----------------
#' Tabulate the output of \code{\link{spe}} function.
#'
#' The option \code{result} allows user to tabulate either sorted estimates or
#' average estimates. For sorted estimates, the table shows user-specified
#' quantile indices, sorted estimates, standard errors, point-wise confidence
#' intervals, and uniform confidence intervals. For average estimates, the
#' table shows average estiamtes, standard errors, and confidence intervals.
#'
#' @param object   The output of \code{\link{spe}} function.
#' @param result   Whether the user wants to see the sorted or the average
#'                 estimates. Default is \code{sorted}, which shows the
#'                 sorted estimates.
#' @param ...      additional arguments affecting the summary produced.
#' @examples
#' data("mortgage")
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec + ltv_med +
#' ltv_high + denpmi + selfemp + single + hischl
#' test <- spe(fm = fm, data = mortgage, var = "black", method = "logit",
#' us = c(2:98)/100, b = 50)
#' summary(test)
#' @export
summary.spe <- function(object, result = c("sorted", "average"), ...) {
  spe <- object$spe
  ape <- object$ape
  us <- object$us
  alpha <- object$alpha
  result <- match.arg(result)
  if (result == "sorted") {
    table <- matrix(0, nrow = length(us), ncol = 6)
    table[, 1] <- unlist(spe[1])
    table[, 2] <- unlist(spe[4])
    table[, 3] <- table[, 1] - qnorm(1 - alpha/2)*table[, 2]
    table[, 4] <- table[, 1] + qnorm(1 - alpha/2)*table[, 2]
    table[, 5] <- unlist(spe[3])
    table[, 6] <- unlist(spe[2])
    rownames(table) <- us
    colnames(table) <- c("Est", "SE", paste0((1 - alpha)*100,"% PLB"),
                         paste0((1 - alpha)*100,"% PUB"),
                         paste0((1 - alpha)*100,"% ULB"),
                         paste0((1 - alpha)*100,"% UUB"))
  } else {
    table <- matrix(0, nrow = 1, ncol = 4)
    table[, 1] <- unlist(ape[1])
    table[, 2] <- unlist(ape[4])
    table[, 3] <- unlist(ape[3])
    table[, 4] <- unlist(ape[2])
    rownames(table) <- "APE"
    colnames(table) <- c("Est", "SE", paste0((1 - alpha)*100,"% LB"),
                         paste0((1 - alpha)*100,"% UB"))
  }
  return(table)
}
