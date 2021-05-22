
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.org/shuowencs/SortedEffects.svg?branch=master)](https://travis-ci.org/shuowencs/SortedEffects)

# SortedEffects

The R package **SortedEffects** implements econometric methods in
Chernozhukov, Fernandez-Val and Luo (2018) and reports *sorted partial
effects*, a collection of estimated partial effects sorted in increasing
order and indexed by percentiles. In contrast to the classical average
partial effect, sorted effect curves completely represent range of
heterogeneous partial effects, which are common in nonlinear predictive
models. The package further allows users to conduct classification
analysis with sorted effects and visualize the results. You can install
the current version of with:

``` r
install.packages("devtools")
devtools::install_github("yuqimemeda/SortedEffects")
```

This repository hosts the (commented) source code and data. For more
theoretical details, please refer to the [original
paper](https://www.econometricsociety.org/publications/econometrica/2018/11/01/sorted-effects-method-discovering-heterogeneous-effects-beyond)
.

This document is structured to answer three curiosity-sparked questions:

  - What do you mean by (heterogeneous) partial effects and
    classification analysis?
  - What function commands does **SortedEffects** offer?
  - How do I use **SortedEffects** for my research?

I will use racial-based discrimination in mortgage lending as an
application to illustrate the basic functionality of **SortedEffects**.
For more detailed explanations of function commands, refer to the
companion vignette.

## Econometrics Background

Many economic inquiries reduce to the following question: holding all
else (control variables) fixed, how does change in a key variable affect
the outcome? Depending on the assumptions we impose, such effects have
different names like treatment or marginal effects. To ease conversation
we call them **partial effects**. To answer the question, economists
often work with nonlinear models. In particular, consider the following
two kinds of nonlinearity: nonlinearity in key variables of interest or
in parameters. The first kind includes OLS and quantile regressions in
which the variable of interest is interacted with control variables. The
second kind includes generalized linear models like Logit and Probit. It
turns out that the partial effects derived under nonlinear settings vary
among obervational units. We call these different effects as
**heterogenous partial effects**.

Before jumping to a concrete example, let’s fix some notations. Let
\(Y\) denote outcome variable of interest and \(X=(T, W)\) denote
regressors. In particular, \(T\) is the variable in interest while \(W\)
means control variables. We link between \(Y\) and \(X\) with a
predictive function \(g(X)\). With more assumptions \(g\) can represent
causal relationship, but it’s not necessary for our purposes.

The economic question is: holding \(W\) fixed, how do changes in \(T\)
affect \(Y\)? The mathematical expression for partial effects is
\[\Delta(x) = g(t_{1}, w)-g(t_{0}, w)\]

for discrete \(T\) and \[\Delta(x) = \partial_{t}g(t, w)\] for \(T\)
being continuous. It turns out that if we work with nonlinear models,
\(\Delta(x)\) depends on values of \(w\).

Here is a concrete probit example:
\(Y=\Phi(T \beta_{1}+W \beta')+\epsilon\), where \(\Phi\) denotes CDF of
standard normal distribution. If \(T\) is continuous, then the partial
effect is \[\Delta(x)=\phi(T\beta_{1}+W\beta')\beta_{1}.\] It is obvious
that \(\Delta(x)\) depends on values of \(W\). Since sample observations
will have different values of \(W\), \(\Delta(x)\) features
heterogeneity. Because \(\Delta(x)\) is unobservable in data, we use
econometric techniques to get an estimate. However, no estimation is
accurate. Therefore we need to quantify the uncertainty of estimation,
and in econometrics this is called **inference**.

It turns out hard to do inference on the estimate of partial effects, so
in applied economic work researchers often report the average partial
effect as an alternative (averaged over the empirical distribution of
\(x\)). But by definition such reports neglect heterogeneity.

Chernozhukov, Fernandez-Val and Luo (forthcoming, *Econometrica*)
provide theoretical underpinnings to conduct estimation and inference on
the partial effect values at any quantile in interest. For example, if
we want to know the partial effects from 2% to 98%, then the method will
report estimate of all these effects and quantification of estimation
uncertainty (confidence bands). Because the estimated partial effects
are sorted in an increasing order, we call them **sorted partial effect
(SPE)**.

With SPE, we can classify observational units into groups based on how
greatly they are affected. For example, researchers might be curious
whether mostly affected observational units have different
characteristics compared to peers who are least affected.

## Example 1: Racial-Based Discrimination in Mortgage Lending

We use data on mortgage applications in Boston from 1990. The Federal
Reserve Bank of Boston collected these data in relation to the Home
Mortgage Disclosure Act (HMDA), which was passed to monitor minority
access to the mortgage market. Providing better access to credit markets
can arguably help the disadvantaged groups escape poverty traps. To
access the data, type in the following command

``` r
data("mortgage")
```

The sample includes 2380 observations corresponding to 2041 white
applicants and 339 black applicants.

We estimate a binary response model where the outcome variable \(Y\) is
an indicator for mortgage denial, the key covariate \(T\) is an
indicator for the applicant being black, and the controls \(W\) contain
financial and other characteristics of the applicant that banks take
into account in the mortgage decisions. The regression formula is
specified as

``` r
fm <- deny ~ black + p_irat + hse_inc + ccred + mcred + pubrec + ltv_med + ltv_high + denpmi + selfemp + 
      single + hischl
```

We use `SPE` command to estimate SPE and APE.

``` r
test <- spe(fm = fm, data = mortgage, var = "black", method = "logit",
            us = c(2:98)/100, b = 500, bc = TRUE)
```

In this line of code, `var.T = "black"` tells the package that \(T\) is
`black`, while `us` specifies the SPE quantile index to range from 2% to
98%. `B` refers to the number of bootstrap for inference. For accuracy
we recommend setting `B` to be 500. The result `ex` is a list containing
SPE and APE with corresponding confidence bands, and for visualization
we can use the `plot.SPE` command as follows:

``` r
plot(x = test, ylim = c(0, 0.25), ylab = "Change in Probability",
     main = "APE and SPE of Being Black on the Prob of Mortgage Denial",
     sub = "Logit Model")
```

We see that the SPE ranges from 2% to 14% while APE is 6%. It’s clear
that reporting APE masks rich heterogeneity of the partial effects.

Now that we know being black affects mortgage denial probability
differently among applicants, we can ask a follow-up question: what are
the average characteristics of most and least affected applicants? Such
questions are referred to as **classification analysis**: first classify
observational units into most/least affected groups and then study
objects of interest like mean and distributions.

To conduct estimation and inference on classification analysis, we use
the `u.CA` command. We first choose the characteristics in interest:
Deny, Black, Debt-to-income, Expenses-to-income, Bad consumer credit,
Bad mortgage credit, Credit problems, Denied mortgage insurance, Medium
loan-to-house, High loan-to-house, Self employed, Single and High school
grad. The corresponding specification is

``` r
t <- c("deny", "p_irat", "black", "hse_inc", "ccred", "mcred", "pubrec",
       "denpmi", "selfemp", "single", "hischl", "ltv_med", "ltv_high")
```

Then we issue to `u.CA` command:

``` r
CA <- ca(fm = fm, data = mortgage, var = "black", method = "logit",
         cl = "both", t = t, b = 500, bc = TRUE)
```

The specification `cl=matrix(c(1,0,0,1), nrow=2)` means we tabulate the
mean of the characteristics for two groups. The result `CA` is a list
containing the estimate, standard errors and joint p-values. We can
tabulate the result as follows:

``` r
summary(CA)
```

We can see that the characteristics are different for these two groups.
