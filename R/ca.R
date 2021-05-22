#' Empirical Classification Analysis (CA) and Inference
#'
#' \code{ca} conducts CA estimation and inference on user-specified objects of
#' interest: first (weighted) moment or (weighted) distribution. Users can use
#' \code{t} to specify variables in interest. When object of interest is
#' moment, use \code{cl} to specify whether want to see averages or difference
#' of the two groups.
#'
#' All estimates are bias-corrected and all confidence bands are monotonized.
#' The bootstrap procedures follow algorithm 2.2 as in Chernozhukov,
#' Fernandez-Val and Luo (2018).
#'
#' @return
#' If \code{subgroup = NULL}, all outputs are whole sample. Otherwise output
#' are subgroup results. When \code{interest = "moment"}, the output is a list
#' showing
#' \itemize{
#'   \item \code{est} Estimates of variables in interest.
#'   \item \code{bse} Bootstrap standard errors.
#'   \item \code{joint_p} P-values that are adjusted for multiplicity to
#'   account for joint testing for all variables.
#'   \item \code{pointwise_p} P-values that doesn't adjust for join testing
#' }
#' If users have further specified \code{cat} (e.g., \code{!is.null(cat)}), the
#' fourth component will be replaced with \code{p_cat}: P-values that are a
#' djusted for multiplicity to account for joint testing for all variables
#' within a category. Users can use \code{\link{summary.ca}} to tabulate the
#' results.
#'
#' When \code{interest = "dist"}, the output is a list of two components:
#' \itemize{
#'   \item \code{infresults} A list that stores estimates, upper and lower
#'   confidence bounds for all variables in interest for least and most
#'   affected groups.
#'   \item \code{sortvar} A list that stores sorted and unique variables in
#'   interest.
#' }
#' We recommend using \code{\link{plot.ca}} command for result visualization.
#'
#' @param fm          Regression formula
#' @param data        The data in use: full sample or subpopulation in interset
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
#'                    this option, users first need to specify \code{var} as
#'                    a factor variable.
#' @param subgroup    Subgroup in interest. Default is \code{NULL}.
#'                    Specifcation should be a logical variable. For example,
#'                    suppose data contain indicator variable for women (female
#'                    if 1, male if 0). If users are interested in women SPE,
#'                    then users should specify
#'                    \code{subgroup = data[, "female"] == 1}.
#' @param samp_weight Sampling weight of data. Input should be a n by 1 vector,
#'                    where n denotes sample size. Default is \code{NULL}.
#' @param taus        Indexes for quantile regression. Default is
#'                    \code{c(5:95)/100}.
#' @param u           Percentile of most and least affected. Default is set to
#'                    be 0.1.
#' @param interest    Generic objects in the least and most affected
#'                    subpopulations. Two options:
#'                    (1) \code{"moment"}: weighted mean of Z in the
#'                    u-least/most affected subpopulation.
#'                    (2) \code{"dist"}: distribution of Z in the u-least/most
#'                    affected subpopulation.
#'                    Default is \code{interest = "moment"}.
#' @param cl          If \code{moment = "interest"}, \code{cl} allows the user
#'                    to get the variables of interest (specified in \code{t}
#'                    option) of the most and least affected groups. The
#'                    default is \code{"both"}, which shows the variables of
#'                    the two groups; the alternative is \code{"diff"}, which
#'                    shows the difference of the two groups. The user can
#'                    use the \code{\link{summary.ca}} to tabulate the
#'                    results, which also contain the standard errors and p-
#'                    values. If \code{interest = "dist"}, this option doesn't
#'                    have any bearing and user can leave it to be the default
#'                    value.
#' @param t           An index for ca object. Should be a 1 by ncol(data)
#'                    indicator vector. Users can either (1) specify names of
#'                    variables of interest directly, or (2) use 1 to indicate
#'                    the variable of interest. For example, total number of
#'                    variables is 5 and interested in the 1st and 3rd vars,
#'                    then specify \code{t = c(1, 0, 1, 0, 0)}.
#' @param cat         P-values in classification analysis are adjusted for
#'                    multiplicity to account for joint testing of zero
#'                    coefficients on for all variables within a category.
#'                    Suppose we have selected specified 3 variables in
#'                    interest: \code{t = c("a", "b", "c")}. Without loss of
#'                    generality, assume \code{"a"} is not a factor, while
#'                    \code{"b"} and \code{"c"} are two factors. Then users
#'                    need to specify as \code{cat = c("b", "c")}. Default is
#'                    \code{NULL}.
#' @param alpha       Size for confidence interval. Shoule be between 0 and 1.
#'                    Default is 0.1
#' @param b           Number of bootstrap draws. Default is 500.
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
#' @param range_cb    When \code{interest = "dist"}, we sort and unique
#'                    variables in interest to estimate weighted CDF. For large
#'                    dataset there can be memory problem storing very many of
#'                    observations, and thus users can provide a Sort value and
#'                    the package will sort and unique based on the weighted
#'                    quantile of Sort. If users don't want this feature, set
#'                    \code{range_cb = NULL}. Default is \code{c(1:99)/100}.
#' @param boot_type   Type of bootstrap. Default is \code{"nonpar"}, and the
#'                    package implements nonparametric bootstrap. The
#'                    alternative is \code{"weighted"}, and the package
#'                    implements weighted bootstrap.
#'
#' @examples
#' data("mortgage")
#' ### Regression Specification
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec +
#' ltv_med + ltv_high + denpmi + selfemp + single + hischl
#' ### Specify characteristics of interest
#' t <- c("deny", "p_irat", "black", "hse_inc", "ccred", "mcred", "pubrec",
#' "denpmi", "selfemp", "single", "hischl", "ltv_med", "ltv_high")
#' ### issue ca command
#' CA <- ca(fm = fm, data = mortgage, var = "black", method = "logit",
#' cl = "diff", t = t, b = 50, bc = TRUE)
#' @importFrom boot boot
#' @importFrom rlist list.cbind
#' @importFrom Hmisc wtd.quantile wtd.mean
#' @importFrom stats quantile rexp qnorm pnorm
#' @importFrom dummies dummy.data.frame
#' @importFrom parallel detectCores
#' @importFrom pbapply setpb startpb closepb
#' @export
ca <- function(fm, data, method = c("ols", "logit", "probit", "QR"),
               var_type = c("binary", "continuous", "categorical"), var,
               compare, subgroup = NULL, samp_weight = NULL, taus = c(5:95)/100,
               u = 0.1, interest = c("moment", "dist"),
               t = c(1, 1, rep(0, dim(data)[2] - 2)), cl = c("both", "diff"),
               cat = NULL, alpha = 0.1, b = 500, parallel = FALSE,
               ncores = detectCores(), seed = 1, bc = TRUE,
               range_cb = c(1:99)/100, boot_type = c("nonpar", "weighted"))
  {
  # ----- Stopping Conditions
  if (alpha >= 1 || alpha <= 0) stop("Please specify a correct size for
                                     hypothesis testing between 0 and 1.")
  if (u >= 1 || u <= 0) stop("Please provide a group classification quantile
                             between 0 and 1.")
  if (typeof(t) == "double" && length(t) != dim(data)[2]) {
    stop("Number of variable indicators don't match.")
  }
  # ----- Replace Null Weight Specification
  if (is.null(samp_weight)) samp_weight <- rep(1, dim(data)[1])
  # ----- Matching Arguments
  method <- match.arg(method)
  var_type <- match.arg(var_type)
  boot_type <- match.arg(boot_type)
  interest <- match.arg(interest)
  temp <- match.arg(cl)
  if (temp == "both") cl <- matrix(c(1, 0, 0, 1), nrow = 2)
  if (temp == "diff") cl <- matrix(c(1, -1), nrow = 2)
  # ---- Index Configuration
  # t is variable names, convert to index
  if (typeof(t) == "character") {
    posit <- match(t, colnames(data))
    ind <- rep(0, dim(data)[2])
    ind[posit] <- 1
    t <- ind
  }
  # Now that t is an index, use it to get prod and vars
  n2 <- dim(data)[2]
  col.num <- c(1:n2)
  # prod is index vector what variables are selected
  prod <- col.num * t
  # The new data is data[, prod]
  # Number of characateristics in interest
  m <- dim(data[, prod])[2]
  if (is.null(m)) m <- 1
  # Extract Names of variables
  vars <- colnames(data)[which(prod != 0)]
  # ----- Auxiliary function
  # The following function creates a dummy matrix for inputted factors
  subsetting <- function(data, vars) {
    # check if the inputted variable is a factor
    if (length(vars) == 1) {
      temp <- matrix(is.factor(data[, vars]))
    } else {
      temp <- matrix(sapply(data[, vars], is.factor), nrow = 1)
    }
    colnames(temp) <- vars
    if (length(which(temp)) == 0) {
      # this means there is no factor inputs
      # then it just subset characteristics of interest
      dat <- data[, vars]
    } else {
      # there is at least one factor input
      # create dummy matrix for inputted factors
      dummydata <- dummy.data.frame(data, names = colnames(temp)[which(temp)],
                                    sep = "_", all = FALSE)
      # create
      otherdata <- data[, colnames(temp)[-which(temp)]]
      dat <- cbind(otherdata, dummydata)
    }
    return(dat)
  }
  # ----- 1. Call to estimate PE
  output <- suppressWarnings(peestimate(fm, data, samp_weight, var_type, var,
                                        compare, method, subgroup, taus))
  pe_est <- output$pe_est
  # if(is.null(subgroup)) is equivalent to if(groupCA=="full")
  # if(!is.null(subgroup)) is equivalent to if(groupCA=="subgroup")
  # ----- 2. u-most and least Affected Groups
  if (method != "QR") {
    if (is.null(subgroup)) {
      # Threshold Values for u-most/least Affected for WHOLE sample
      effect_high <- wtd.quantile(pe_est, samp_weight, 1 - u)
      effect_low <- wtd.quantile(pe_est, samp_weight, u)
      if (interest == "moment") {
        # subset data for interest == moment
        mat <- subsetting(data, vars)
        mat$.w <- samp_weight
        high_affected <- mat[pe_est >= effect_high, ]
        low_affected <- mat[pe_est <= effect_low, ]
      } else if (interest == "dist") {
        # Get weight for the two groups. If weight is NULL, it does nothing
        data$.w <- samp_weight
        # Two affected groups for dist purpose
        high_affected <- data[pe_est >= effect_high, ]
        low_affected <- data[pe_est <= effect_low, ]
      }
    } else {
      # SUB sample
      pesub_est <- output$pesub_est
      pesub_w <- output$samp_weight_sub
      # Threshold Values for u-most/least Affected
      effect_high <- wtd.quantile(pesub_est, pesub_w, 1 - u)
      effect_low <- wtd.quantile(pesub_est, pesub_w, u)
      data$.w <- samp_weight
      subdata <- data[subgroup, ]
      if (interest == "moment") {
        vars2 <- c(vars, ".w")
        mat_sub <- subsetting(subdata, vars2)
        high_affected <- mat_sub[pesub_est >= effect_high, ]
        low_affected <- mat_sub[pesub_est <= effect_low, ]
      } else if (interest == "dist") {
        high_affected <- subdata[pesub_est >= effect_high, ]
        low_affected <- subdata[pesub_est <= effect_low, ]
      }
    }
  }
  # QR requires special treatment due to stacking of quantile indices
  if (method == "QR") {
    if (is.null(subgroup)) {
      # This is where QR is different.
      effect_high <- wtd.quantile(pe_est, matrix(samp_weight, ncol = 1,
                                                 nrow = nrow(pe_est),
                                                 byrow = FALSE), 1 - u)
      effect_low <- wtd.quantile(pe_est, matrix(samp_weight, ncol = 1,
                                                nrow = nrow(pe_est),
                                                byrow = FALSE), u)
      if (interest == "moment") {
        mat <- subsetting(data, vars)
        mat$.w <- samp_weight
        # Kronecker function doesn't work for all data types, so I wrote this
        # alternative
        mesh <- mat[rep(1:nrow(mat), times = length(taus)), ]
        high_affected <- mesh[pe_est >= effect_high, ]
        low_affected <- mesh[pe_est <= effect_low, ]
      } else if (interest == "dist") {
        # Obtain sampling weight of each group
        data$.w <- samp_weight
        mesh <- data[rep(1:nrow(data), times = length(taus)), ]
        # Two Affected Groups
        high_affected <- mesh[pe_est >= effect_high, ]
        low_affected <- mesh[pe_est <= effect_low, ]
      }
    } else {# Subsample: use PEsub_est
      pesub_est <- output$pesub_est
      pesub_w <- output$samp_weight_sub
      # Threshold Values for u-most/least Affected
      effect_high <- wtd.quantile(pesub_est, pesub_w, 1 - u)
      effect_low <- wtd.quantile(pesub_est, pesub_w, u)
      data$.w <- samp_weight
      subdata <- data[subgroup, ]
      if (interest == "moment") {
        vars2 <- c(vars, ".w")
        mat_sub <- subsetting(subdata, vars2)
        mesh <- mat_sub[rep(1:nrow(mat_sub), times = length(taus)), ]
        high_affected <- mat_sub[pesub_est >= effect_high, ]
        low_affected <- mat_sub[pesub_est <= effect_low, ]
      } else if (interest == "dist") {
        mesh <- subdata[rep(1:nrow(subdata), times = length(taus)), ]
        # Two Affected Groups
        high_affected <- mesh[pesub_est >= effect_high, ]
        low_affected <- mesh[pesub_est <= effect_low, ]
      }
    }
  }
  # Obtain the sampling weight for each group (If samp_weight = NULL, then
  # these two vars are NULL)
  weight_high <- high_affected$.w
  weight_low <- low_affected$.w
  # Remove the .w
  high_affected$.w <- NULL
  low_affected$.w <- NULL
  data$.w <- NULL
  ##### Classicification Analysis Starts Here #####
  # ----- 3. Specify object of interest
  if (interest == "moment") {
    # weighted mean of var in interest for the u-least/most affected groups
    interest_high <- apply(as.matrix(high_affected), 2, wtd.mean,
                           weights = weight_high, na.rm = TRUE)
    interest_low <- apply(as.matrix(low_affected), 2, wtd.mean,
                          weights = weight_low, na.rm = TRUE)
    varname <- colnames(as.matrix(high_affected))
    # Average of var in interest for the two groups (m * 2), where m is number
    # of variables in interest
    est_interest <- cbind(interest_high, interest_low)
    #est_interest <- cbind(interest_low, interest_high)
    # Hypothesis in interest (m * L matrix)
    # Component (a, b) meaning: b-th hypothesis statistics for a-th variable
    # in interest
    H <- est_interest %*% cl
    # Reshape: 1 * mL
    H <- matrix(H, nrow = 1)
    # Assign names
    m <- length(varname)
    colnames(H) <- rep(varname, length(H)/m)
  } else if (interest == "dist") {
    # weighted CDF of var in interest for the u-least/most affected groups
    # Sort and unique eliminate repetitive obs. Note the list "temp" is ragged.
    # temp0 will be used at the end of the function
    # Depending on whether user provides Sort
    if (is.null(range_cb)) {
      temp0 <- temp <- lapply(as.data.frame(data[, prod]), function(x)
        sort(unique(x)))
    } else {
      temp0 <- temp <- lapply(as.data.frame(data[, prod]), function(x)
        sort(unique(wtd.quantile(x, samp_weight, range_cb))))
    }
    names(temp0) <- names(temp) <- vars
    # Double each list and then split each list into a k by 2 matrix, where k
    # denotes number of elements for each var
    temp <- lapply(temp, rep, 2)
    temp <- lapply(temp, matrix, ncol = 2)
    # I hate for loops in general, but this loop is short: it assigns names to
    # each matrix within the list. In general, I find it easier to use for loop
    # if the purpose is to modify part of dataframe
    for (x in 1:length(temp)) colnames(temp[[x]]) <- rep(names(temp)[[x]], 2)
    # Get the trim tails (for inference and plotting)
    trim <- function(x){
      varname <- colnames(x)[1]
      notails <- (x[, 1] >= wtd.quantile(unlist(data[varname]), samp_weight, 0.02) &
                    x[, 1] <= wtd.quantile(unlist(data[varname]), samp_weight, 0.98))
      notails <- rep(notails, 2)
    }
    trimtails <- lapply(temp, trim)
    # The following is a function to calculate weighted cdf
    fun <- function(a){ # a is a name component in the list
      if (is.null(samp_weight)) {
        # Take care of no samp weight scenario
        weight_low <- rep(1, length(unlist(low_affected[, colnames(a)[1]])))
        weight_high <- rep(1, length(unlist(high_affected[, colnames(a)[1]])))
      }
      cdf_low <- lapply(a[, 1], wcdf, x = unlist(low_affected[, colnames(a)[1]]),
                        w = weight_low) # low group cdf
      cdf_high <- lapply(a[, 2], wcdf, x = unlist(high_affected[, colnames(a)[[1]]]),
                         w = weight_high) # high group cdf
      output <- cbind(cdf_low, cdf_high)
    }
    # Call 'fun' to calculate wcdfs for both groups for all vars in interest
    # (least group 1st col, most group 2nd col)
    cdf_bundle <- lapply(temp, fun)
    index <- lengths(cdf_bundle)
    # Reshape
    H <- lapply(cdf_bundle, matrix, nrow = 1)
    H <- list.cbind(H) # estimates of weighted cdf
    trimtails <- lapply(trimtails, matrix, nrow = 1)
    trimtails <- list.cbind(trimtails) # long vector of trimtail logics
  }
  #  ----- 4. The Bootstrap Statistics Function
  # set a bootstrap counting variable for the purpose of showing a progress bar
  rep_count <- 1
  # DGP for weighted bootstrap
  data_rg <- function(data, mle){
    n <- dim(data)[1]
    # Exponential weights
    multipliers  <- rexp(n)
    # Sampling weight of data.bs
    # Multiply 20000 to enlarge the weight (otherwise the effect_high and
    # effect_low in stat.boot.CA.weight
    # function would be the same. Same strategy applies to u.Subpop. For SPE it
    # doesn't matter since we plot the entire quantile index.)
    weight <- (multipliers/sum(multipliers)) * samp_weight * 20000
    data$.w <- weight
    return(data)
  }
  # DGP for nonparametric bootstrap
  data_non <- function(data, mle){
    n <- dim(data)[1]
    multipliers <- as.vector(table(factor(sample(n, n, replace = T),
                                          levels = c(1:n))))
    # Sampling weight of data.bs
    weight <- (multipliers/sum(multipliers)) * samp_weight * 20000
    data$.w <- weight
    return(data)
  }
  # Function that computes bootstrap statistics in each draw
  stat_boot_CA_weight <- function(data){
    # set up a progress bar to document the bootstrap progress
    setpb(pb, rep_count)
    rep_count <<- rep_count + 1
    output_bs <- suppressWarnings(peestimate(fm, data = data,
                                             samp_weight = data$.w, var_type,
                                             var, compare, method, subgroup,
                                             taus))
    pe_est <- output_bs$pe_est
    if (method != "QR") {
      if (is.null(subgroup)) {
        # Threshold Values for u-most/least Affected
        effect_high <- wtd.quantile(pe_est, data$.w, 1 - u)
        effect_low <- wtd.quantile(pe_est, data$.w, u)
        if (interest == "moment") {
          mat <- subsetting(data, vars)
          mat$.w <- data$.w
          high_affected <- mat[pe_est >= effect_high, ]
          low_affected <- mat[pe_est <= effect_low, ]
        } else if (interest == "dist") {
          # Two affected groups
          high_affected <- data[pe_est >= effect_high, ]
          low_affected <- data[pe_est <= effect_low, ]
        }
      } else {
        pesub_est <- output_bs$pesub_est
        pesub_w <- output_bs$samp_weight_sub
        # Threshold Values for u-most/least Affected
        effect_high <- wtd.quantile(pesub_est, pesub_w, 1 - u)
        effect_low <- wtd.quantile(pesub_est, pesub_w, u)
        subdata <- data[subgroup, ]
        if (interest == "moment") {
          vars2 <- c(vars, ".w")
          mat_sub <- subsetting(subdata, vars2)
          high_affected <- mat_sub[pesub_est >= effect_high, ]
          low_affected <- mat_sub[pesub_est <= effect_low, ]
        } else if (interest == "dist") {
          high_affected <- subdata[pesub_est >= effect_high, ]
          low_affected <- subdata[pesub_est <= effect_low, ]
        }
      }
    } # QR requires special treatment due to stacking of quantile indices
    if (method == "QR") {
      # Full sample: use PE_est
      if (is.null(subgroup)) {
        # Threshold Values for u-most/least Affected
        effect_high <- wtd.quantile(pe_est, matrix(data$.w, ncol = 1,
                                                   nrow = nrow(pe_est),
                                                   byrow = FALSE), 1 - u)
        effect_low <- wtd.quantile(pe_est, matrix(data$.w, ncol = 1,
                                                  nrow = nrow(pe_est),
                                                  byrow = FALSE), u)
        if (interest == "moment") {
          mat <- subsetting(data, vars)
          mat$.w <- data$.w
          mesh <- mat[rep(1:nrow(mat), times = length(taus)), ]
          high_affected <- mesh[pe_est >= effect_high, ]
          low_affected <- mesh[pe_est <= effect_low, ]
        } else if (interest == "dist") {
          mesh <- data[rep(1:nrow(data), times = length(taus)), ]
          # Two Affected Groups
          high_affected <- mesh[pe_est >= effect_high, ]
          low_affected <- mesh[pe_est <= effect_low, ]
        }
      } else {
        pesub_est <- output_bs$pesub_est
        pesub_w <- output_bs$samp_weight_sub
        # Subsample: use pesub_est
        # Threshold Values for u-most/least Affected
        effect_high <- wtd.quantile(pesub_est, pesub_w, 1 - u)
        effect_low <- wtd.quantile(pesub_est, pesub_w, u)
        subdata <- data[subgroup, ]
        if (interest == "moment") {
          vars2 <- c(vars, ".w")
          mat_sub <- subsetting(data, vars2)
          mesh <- mat_sub[rep(1:nrow(mat_sub), times = length(taus)), ]
          high_affected <- mesh[pesub_est >= effect_high, ]
          low_affected <- mesh[pesub_est <= effect_low, ]
        } else if (interest == "dist") {
          mesh <- subdata[rep(1:nrow(subdata), times = length(taus)), ]
          # Two Affected Groups
          high_affected <- mesh[pesub_est >= effect_high, ]
          low_affected <- mesh[pesub_est <= effect_low, ]
        }
      }
    }
    # Obtain the sampling weight for each group
    weight_high <- high_affected$.w
    weight_low <- low_affected$.w
    # Remove the .w
    high_affected$.w <- NULL
    low_affected$.w <- NULL
    data$.w <- NULL
    if (interest == "moment") {
      # weighted mean of var in interest for the u-least/most affected groups
      interest_high <- apply(as.matrix(high_affected), 2, wtd.mean,
                             weights = weight_high, na.rm = TRUE)
      interest_low <- apply(as.matrix(low_affected), 2, wtd.mean,
                            weights = weight_low, na.rm = TRUE)
      varname <- colnames(as.matrix(high_affected))
      # Average of var in interest for the two groups (m * 2), where m is
      # number of variables in interest
      est_interest <- cbind(interest_high, interest_low)
      # est_interest <- cbind(interest_low, interest_high)
      # Hypothesis in interest (m * L matrix)
      # Component (a, b) meaning: b-th hypothesis statistics for a-th variable
      # in interest
      H <- est_interest %*% cl
      # Reshape: 1 * mL
      H <- matrix(H, nrow = 1)
      # Assign names
      m <- length(varname)
      colnames(H) <- rep(varname, length(H)/m)
    }
    if (interest == "dist") {
      fun2 <- function(a, weight_high = weight_high, weight_low = weight_low){
        # a is a name component in the list
        cdf_low <- lapply(a[, 1], wcdf, x = unlist(low_affected[, colnames(a)[1]]), w = weight_low) # low group cdf
        cdf_high <- lapply(a[, 2], wcdf, x = unlist(high_affected[, colnames(a)[1]]), w = weight_high) # high group cdf
        output <- cbind(cdf_low, cdf_high)
      }
      # Call 'fun' to calculate wcdfs for both groups for all vars in interest
      # (least group 1st col, most group 2nd col)
      cdf_bundle <- lapply(temp, fun2)
      index <- lengths(cdf_bundle)
      # Reshape
      H <- lapply(cdf_bundle, matrix, nrow = 1)
      H <- unlist(list.cbind(H)) # estimates of weighted cdf
    }
    out <- H
  }
  # -----  5. Call boot for bootstrap
  set.seed(seed)
  if (parallel == FALSE) ncores <- 1
  if (boot_type == "nonpar") {
    # print a message showing how many cores are used
    data$.w <- samp_weight
    cat(paste("Using", ncores, "CPUs now.\n"))
    # set up a progress bar
    pb <- startpb(min = 0, max = b)
    result_boot <- boot(data = data, statistic = stat_boot_CA_weight,
                        sim = "parametric", ran.gen = data_non, mle = 0,
                        parallel = "multicore", ncpus = ncores, R = b)
    closepb(pb)
    data$.w <- NULL
  } else if (boot_type == "weighted") {
    data$.w <- samp_weight
    cat(paste("Using", ncores, "CPUs now.\n"))
    pb <- startpb(min = 0, max = b)
    result_boot <- boot(data = data, statistic = stat_boot_CA_weight,
                        sim = "parametric", ran.gen = data_rg, mle = 0,
                        parallel = "multicore", ncpus = ncores, R = b)
    closepb(pb)
    data$.w <- NULL
  }
  # -----  6. Inference
  if (interest == "moment") {
    # draws.H is B * mL
    draws_H <- result_boot$t
    colnames(draws_H) <- rep(varname, dim(cl)[2])
    # Zc is B * mL
    Zc <- draws_H - matrix(H, nrow = b, ncol = length(H), byrow = TRUE)
    # sig is 1 * mL (bootstrap standard error) (Pointwise SE)
    sig <- (apply(Zc, 2, quantile, 0.75, na.rm = TRUE) -
              apply(Zc, 2, quantile, 0.25, na.rm = TRUE)) / (qnorm(0.75) - qnorm(0.25))
    # t.tilde is B * mL
    t_tilde <- apply(abs(Zc)/matrix(sig, nrow = b, ncol = length(sig),
                                    byrow = TRUE), 1, max, na.rm = TRUE)
    # Bias correction
    H_bc <- 2 * H - apply(draws_H, 2, mean, na.rm = TRUE) # 1 * mL
    # We calculate the pointwise p-value (non and bc corrected)
    pw_p <- 1 - pnorm(abs(H/sig))
    pw_p_bc <- 1 - pnorm(abs(H_bc/sig))
    # Joint p-values
    stat <- abs(H_bc)/sig
    # Proportions of the B draws of t.tilde that are greater than s
    # 1 * mL. First m columns represent joint p-vals for the m vars for the
    # first hypothesis.
    j_pvals <- sapply(stat, epvals, zs = t_tilde)
    # Not bias-corrected joint p-values
    stat_un <- abs(H)/sig
    j_pvals_un <- sapply(stat_un, epvals, zs = t_tilde)
    # Now work on p-vals accounting for testing all vars within one category
    if (!is.null(cat)) {
      # first tease out characteristics that are not factors
      noncat <- vars[-match(cat, vars)]
      # for noncat subset pointwise p-values from previous results
      if (dim(cl)[2] == 1) {
        noncat_p <- pw_p[, noncat]
        noncat_p_bc <- pw_p_bc[, noncat]
      } else if (dim(cl)[2] == 2) {
        pw_p_most <- matrix(pw_p[, 1:length(varname)], nrow = 1)
        pw_p_least <- matrix(pw_p[, -(1:length(varname))], nrow = 1)
        colnames(pw_p_most) <- colnames(pw_p_least) <- varname
        noncat_p_most <- pw_p_most[, noncat]
        noncat_p_least <- pw_p_least[, noncat]
        pw_p_most_bc <- matrix(pw_p_bc[, 1:length(varname)], nrow = 1)
        pw_p_least_bc <- matrix(pw_p_bc[, -(1:length(varname))], nrow = 1)
        colnames(pw_p_most_bc) <- colnames(pw_p_least_bc) <- varname
        noncat_p_most_bc <- pw_p_most_bc[, noncat]
        noncat_p_least_bc <- pw_p_least_bc[, noncat]
      }
      # now we calculate joint p-values within each factor category
      # for each factor, get the dummy names, subset from the bootstrap matrix
      # and conduct inference
      newfunc <- function(x, bootstrap_mat, H, b, cl, varname) {
        # get the dummy name for each factor
        name <- colnames(subsetting(data, x))
        if (dim(cl)[2] == 1) {
          # get bootstrap matrix and estimates within each category
          boot_mat <- bootstrap_mat[, name]
          H_mat <- H[, name]
        } else if (dim(cl)[2] == 2) {
          # subset only gets the first group, need to split the matrix
          # and get the second group, too
          most <- bootstrap_mat[, 1:length(varname)]
          least <- bootstrap_mat[, -(1:length(varname))]
          most_H <- matrix(H[1:length(varname)], nrow = 1)
          least_H <- matrix(H[-(1:length(varname))], nrow = 1)
          colnames(most_H) <- colnames(least_H) <- varname
          mat_most <- most[, name]
          mat_least <- least[, name]
          H_most <- most_H[, name]
          H_least <- least_H[, name]
        }
        # conduct inference
        inference <- function(h, bt, rownum) {
          Zc <- bt - matrix(h, nrow = rownum, ncol = length(h), byrow = TRUE)
          sig <- (apply(Zc, 2, quantile, 0.75, na.rm = TRUE) -
                    apply(Zc, 2, quantile, 0.25, na.rm = TRUE)) / (qnorm(0.75) - qnorm(0.25))
          t_tilde <- apply(abs(Zc)/matrix(sig, nrow = b, ncol = length(sig),
                                          byrow = TRUE), 1, max, na.rm = TRUE)
          H_bc <- 2 * h - apply(bt, 2, mean, na.rm = TRUE)
          stat <- abs(H_bc)/sig
          # bias corrected categorical pvalues
          cat_pvals <- sapply(stat, epvals, zs = t_tilde)
          stat_un <- abs(h)/sig
          # non bias corrected categorical pvalues
          cat_pvals_un <- sapply(stat_un, epvals, zs = t_tilde)
          out <- cbind(cat_pvals, cat_pvals_un)
        }
        if (dim(cl)[2] == 1) {
          pvalues <- inference(h = H_mat, bt = boot_mat, rownum = b)
        } else if (dim(cl)[2] == 2) {
          pvalues_most <- inference(h = H_most, bt = mat_most, rownum = b)
          pvalues_least <- inference(h = H_least, bt = mat_least, rownum = b)
          pvalues <- cbind(pvalues_most, pvalues_least)
        }
        return(pvalues)
      }
      result <- sapply(cat, newfunc, bootstrap_mat = draws_H, H = H, b = b,
                       cl = cl, varname = varname)
      # Now we tabulate the results
      # placeholder
      if (dim(cl)[2] == 1) {
        mat <- matrix(0, nrow = 1, ncol = 2)
        for (i in 1:length(result)) mat <- rbind(mat, unlist(result[[i]]))
        # the first column is cat pvalues; second column is not bias corrected
        mat <- mat[-1, ]
        # combine with previous noncat pvalues
        mat <- rbind(cbind(noncat_p_bc, noncat_p), mat)
        # mat might have different orders from est, bse, so need to reorder and
        # match.
        # the order of variables in mat should line up with that in varname
        rightorder <- match(varname, rownames(mat))
        mat <- mat[rightorder, ]
        p_cat_bc <- mat[, 1]
        p_cat_un <- mat[, 2]
      } else if (dim(cl)[2] == 2) {
        mat <- matrix(0, nrow = 1, ncol = 4)
        for (i in 1:length(result)) mat <- rbind(mat, unlist(result[[i]]))
        # first two columns for most affected, last two for the least
        mat <- mat[-1, ]
        # combine with previous noncat pvalues
        noncategory <- cbind(noncat_p_most_bc, noncat_p_least_bc,
                             noncat_p_most, noncat_p_least)
        mat <- rbind(noncategory, mat)
        # reorder
        rightorder <- match(varname, rownames(mat))
        mat <- mat[rightorder, ]
        p_cat_bc <- mat[, 1:2]
        p_cat_un <- mat[, 3:4]
      }
      if (bc == TRUE) {
        output <- list(est = H_bc, bse = sig, joint_p = j_pvals,
                       p_cat = p_cat_bc, vars = varname, type = temp)
      } else {
        output <- list(est = H, bse = sig, joint_p = j_pvals_un,
                       p_cat = p_cat_un, vars = varname, type = temp)
      }
    } else {# if cat is null
      if (bc == TRUE) {
        output <- list(est = H_bc, bse = sig, joint_p = j_pvals,
                       pointwise_p = pw_p_bc, vars = varname, type = temp)
      } else {
        output <- list(est = H, bse = sig, joint_p = j_pvals_un,
                       pointwise_p = pw_p, vars = varname, type = temp)
      }
    }
    # claim the object to be class "ca"
    output <- structure(output, class = "ca")
    return(output)
  } else if (interest == "dist") {
    draws_H <- result_boot$t # nrow is B
    # Zc <- mapply("-", draws.H, lapply(H, rep, b))
    Zc <- draws_H - matrix(unlist(H), nrow = b, ncol = length(H), byrow = TRUE)
    bse <- (apply(Zc, 2, quantile, .75, na.rm = TRUE) -
              apply(Zc, 2, quantile, .25, na.rm = TRUE)) / (qnorm(0.75) - qnorm(.25))
    # Bias corrected estimation
    bm <- apply(draws_H, 2, mean, na.rm = TRUE)
    H_bc <- 2*unlist(H) - bm
    # Extract the var in interest by Assigning names to H
    # The following chunk of code is to get A: column of variable names
    indexname <- names(index)
    A <- NULL
    A[indexname] <- list(NULL)
    for (i in 1:length(A)) A[i] <- names(A[i])
    # A small function to assign names
    assign <- function(x, y = index){
      num <- which(names(y) == x)
      x <- rep(names(y)[num], y[num])
    }
    A <- lapply(A, assign)
    A <- lapply(A, matrix, nrow = 1)
    A <- list.cbind(A)
    colnames(Zc) <- A
    # The following function outputs a bundle of bias-corrected estimate, u
    # pper and lower bound for both groups for a variable in interest
    CDFinf <- function(x){
      Zc_sub <- x[1:b, ]
      bse_sub <- bse_sub1 <- x[b + 1, ]
      bse_sub1[bse_sub1 == 0] <- NA
      est_sub <- x[b + 2, ]
      trim_sub <- as.logical(x[b + 3, ])
      tempo <- abs(Zc_sub)/matrix(bse_sub1, nrow = b, ncol = length(bse_sub1),
                                  byrow = TRUE)
      t_sub <- apply(tempo[, trim_sub], 1, max, na.rm = TRUE) # B by 1
      crt_sub <- quantile(t_sub, 1 - alpha) # scalar
      # For Least-Affected Group (bias-corrected)
      least_upper <- est_sub[1:(length(est_sub)/2)] + crt_sub *
        bse_sub[1:(length(bse_sub)/2)]
      least_lower <- est_sub[1:(length(est_sub)/2)] - crt_sub *
        bse_sub[1:(length(bse_sub)/2)]
      # For Most-Affected Group (bias-corrected)
      most_upper <- est_sub[-(1:(length(est_sub)/2))] + crt_sub *
        bse_sub[-(1:(length(bse_sub)/2))]
      most_lower <- est_sub[-(1:(length(est_sub)/2))] - crt_sub *
        bse_sub[-(1:(length(bse_sub)/2))]
      bundle <- list(least = est_sub[1:(length(est_sub)/2)],
                     least_up = least_upper, least_low = least_lower,
                     most = est_sub[-(1:(length(est_sub)/2))],
                     most_up = most_upper, most_low = most_lower)
      # Impose shape restriction
      bundle <- lapply(bundle, function(y) sort(y * (y > 0 & y < 1) + (y >= 1)))
    }
    if (bc == TRUE) {
      # Returns a list storing statistics in interest
      # Zc (row 1:B), bse (row B+1), H.bc (row B+2), trimtails (row B+3)
      info <- rbind(Zc, bse, H_bc, trimtails)
      info_sub <- lapply(split.data.frame(t(info), colnames(info)), t)
      # Now call function "CDFinf" to get the bundled result as a list.
      infresults <- lapply(info_sub, CDFinf)
      # Return output
      output <- list(infresults = infresults, sortvar = temp0, alpha = alpha)
    } else if (bc == FALSE) {
      # For non bias-correction case
      info_n <- rbind(Zc, bse, unlist(H), trimtails) # For non bias-corrected
      info_sub_n <- lapply(split.data.frame(t(info_n), colnames(info_n)), t)
      infresults_n <- lapply(info_sub_n, CDFinf)
      output <- list(infresults = infresults_n, sortvar = temp0, alpha = alpha)
    }
    # claim the output as an object of class "cadist" so as to use S3
    output <- structure(output, class = "ca")
    return(output)
  }
}

# ------ Two Auxiliary Functions
# 1. Function to obtain empirical p-values
epvals <- function(z, zs){
  return(mean(zs > z))
}
# 2. Weighted distribution function estimation
#' @importFrom stats weighted.mean
wcdf <- function(y, x, w){
  Ff <- weighted.mean((x <= y), w)
  return(Ff)
}
# ----- Summary (Moments of specified variables in interest for least/most
# affected groups) ------
#' Return the output of \code{\link{ca}} function.
#'
#' @param object   Output of \code{\link{ca}} command with
#'                 \code{interest = "moments"}.
#' @param ...      additional arguments affecting the summary produced.
#' @examples
#' data("mortgage")
#' ### Regression Specification
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec +
#' ltv_med + ltv_high + denpmi + selfemp + single + hischl
#' ### Specify characteristics of interest
#' t <- c("deny", "p_irat", "black", "hse_inc", "ccred", "mcred", "pubrec",
#' "denpmi", "selfemp", "single", "hischl", "ltv_med", "ltv_high")
#' ### Issue ca command
#' CA <- ca(fm = fm, data = mortgage, var = "black", method = "logit",
#' cl = "both", t = t, b = 50, bc = TRUE)
#' ### Report summary table
#' summary(CA)
#' @export
summary.ca <- function(object, ...) {
  type <- object$type
  n_col <- length(object$est)/length(object$vars)
  est <- matrix(object$est, ncol = n_col)
  se <- matrix(object$bse, ncol = n_col)
  joint_p <- matrix(object$joint_p, ncol = n_col)
  if (type == "both") {
    table <- matrix(0, ncol = 4, nrow = length(object$vars))
    # Least Affected Bias-corrected estimate
    table[, 1] <- est[, 1]
    # Corresponding SE
    table[, 2] <- se[, 1]
    # Most affected
    table[, 3] <- est[, 2]
    # Corresponding SE
    table[, 4] <- se[, 2]
    # assign names to each row
    rownames(table) <- colnames(object$est)[1:length(object$vars)]
    colnames(table) <- c("Most", "SE", "Least", "SE")
  } else {
    table <- matrix(0, ncol = 3, nrow = length(object$vars))
    table[, 1] <- est
    table[, 2] <- se
    table[, 3] <- joint_p
    rownames(table) <- colnames(object$est) # assign names to each row
    colnames(table) <- c("Estimate", "SE", "JP-vals")
  }
  # add p_cat or pointwise_p
  if (!is.null(object$p_cat)) {
    new_p <- matrix(object$p_cat, ncol = n_col)
  } else if (!is.null(object$pointwise_p)) {
    new_p <- matrix(object$pointwise_p, ncol = n_col)
  }
  if (type == "diff") {
    temp <- matrix(new_p, ncol = 1, nrow = length(object$vars))
    if (!is.null(object$p_cat)) colnames(temp) <- "Cat P-vals"
    if (is.null(object$p_cat)) colnames(temp) <- "PW P-vals"
    table <- cbind(table, temp)
  }
  return(table)
}
# ----------------- Plotting (CDFs of a continuous variable in interest for
# least/most affected groups) ------------
#' Distribution plotting
#'
#' Plots distributions and joint uniform confidence bands of variables in
#' interest from \code{\link{ca}} command.
#'
#'
#' @param x           Output of \code{\link{ca}} command with
#'                    \code{interest = "dist"}.
#' @param var         Name of variable for plotting
#' @param main        Main title of the plot. Defualt is NULL.
#' @param sub         Sub title of the plot. Default is NULL.
#' @param xlab        x-axis label. Default is NULL.
#' @param ylab        y-axis label. Default is NULL.
#' @param ...         graphics parameters to be passed to the plotting
#'                    routines.
#'
#' @examples
#' data("mortgage")
#' ### Regression Specification
#' fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec +
#' ltv_med + ltv_high + denpmi + selfemp + single + hischl
#' ### Specify characteristics of interest for plotting
#' t2 <- "p_irat"
#' ### issue ca command
#' CAdist <- ca(fm = fm, data = mortgage, var = "black", method = "logit",
#' t = "p_irat", b = 50, interest = "dist")
#' ### plotting
#' plot(CAdist, var = "p_irat", ylab = "Prob",
#' xlab = "Monthly Debt-to-Income Ratio", sub = "logit model")
#'
#' @importFrom graphics axis legend lines plot polygon
#' @export
plot.ca <- function(x, var, main = NULL, sub = NULL, xlab = NULL,
                    ylab = NULL, ...) {
  alpha = x$alpha
  if (is.na(match(var, names(x$sortvar)))) {
    stop("The variable must be consistent with interest specification in the ca
         command.")
  }
  xvar <- unlist(x$sortvar[var])
  yvars <- unlist(x$infresults[var])
  q <- length(xvar)
  plot(xvar, yvars[(3*q + 1):(4*q)], type = "n", xlim = range(xvar),
       ylim = c(0, 1), log = "", main, sub, xlab, ylab, col = 4, lwd = 2)
  # CB of least affected
  polygon(c(xvar, rev(xvar)),
          c(yvars[(q + 1):(2*q)], rev(yvars[(2*q + 1):(3*q)])), density = 60,
          border = F, col = 'darkcyan', lty = 1, lwd = 1)
  # CB of most affected
  polygon(c(xvar, rev(xvar)),
          c(yvars[(4*q + 1):(5*q)], rev(yvars[(5*q + 1):(6*q)])),
          density = 60, border = F, col = 'tomato', lty = 1, lwd = 1)
  # Est of least and most affected
  lines(xvar, yvars[1:q], lwd = 2, col = 4)
  lines(xvar, yvars[(3*q + 1):(4*q)], lwd = 2, col = 2)
  legend(x = "topleft", col = c(4, 2, "darkcyan","tomato"),
         lwd = c(1, 1, 5, 5), lty = c(1, 1, 1, 1), bty = 'n',
         legend = c("Least Affected","Most Affected",
                    paste0((1 - alpha)*100, "% CB"),
                    paste0((1 - alpha)*100, "% CB")))
}
