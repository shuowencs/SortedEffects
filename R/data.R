#' Mortgage Denial
#'
#' @format Contains the data on mortgage application in Boston from 1990,
#' (Munnell et al., 1996.)
#' We obtain the data from the companion website of Stock and Watson (2011).
#' The file contains the following variables:
#' \describe{
#' \item{deny}{indicator for mortgage application denied}
#' \item{p_irat}{monthly debt to income ratio}
#' \item{black}{indicator for black applicant}
#' \item{hse_inc}{monthly housing expenses to income ratio}
#' \item{loan_val}{loan to assessed property value ratio}
#' \item{ccred}{consumer credit score with 6 categories. 1 if no "slow"
#' payments or delinquencies, 2 if one or two "slow" payments or delinquencies,
#' 3 if more than two "slow" payments or delinguencies, 4 if insufficient
#' credit history for determination, 5 if delinquent credit history with
#' payment 60 days overdue, and 6 if delinquent credit history with payments
#' 90 days overdue.}
#' \item{mcred}{mortgage credit score with 4 categories. 1 if no late mortgage
#' payments, 2 if no mortgage payment history, 3 if one or two late mortage
#' payments, and 4 if more than two late mortgages payments}
#' \item{pubrec}{indicator for any public record of credit problems: bankruptcy
#' , charge-offs, collection actions}
#' \item{denpmi}{indicator for applicant applied for mortgage insurance and was
#' denied}
#' \item{selfemp}{indicator for self-employed applicant}
#' \item{single}{indicator for single applicant}
#' \item{hischl}{indicator for high school graduated applicant}
#' \item{probunmp}{1989 Massachusetts unemployment rate in the applicant's
#' history}
#' \item{condo}{indicator for unit is a condominium}
#' \item{ltv_med}{indicator for medium loan to property value ratio [.80, .95]}
#' \item{ltv_high}{indicator for high loan to property value ratio >.95}
#'}
#'@source Munnell, Alicia, Geoffrey Tootell, Lynn Browne, and James McEneaney,
#'"Mortgage Lending in Boston: Interpreting HMDA Data", The American Economic
#'Review, 1996.
"mortgage"

#' Wage Data
#'
#' @format Consists of white, non-hispanic individuals aging from 25 to 64 and
#' working more than 35 hours per week during at least 50 weeks of the year.
#' Excludes self-employed, individuals living in group quarters; individuals in
#' the military, agricultural or private household sectors; individuals with
#' inconsistent reports on earnings and employment status; individuals with
#' allocated or missing information in any of the variables used in the
#' analysis; individuals with hourly wage rate below $3. Contains 32,523
#' workers including 18,137 men and 14,382 women. The file contains the
#' following variables:
#' \describe{
#' \item{lnw}{log of hourly wages}
#' \item{weight}{CPS sampling weight}
#' \item{female}{gender indicator: 1 if female}
#' \item{exp1}{max(age-years of educ-7, 0)}
#' \item{exp2}{exp1^2/100}
#' \item{exp3}{exp1^3/100}
#' \item{exp4}{exp1^4/100}
#' \item{occ}{Aggregated occupation with 5 categories: managers, service,
#' sales, construction and production.}
#' \item{ind}{Aggregated industry with 12 categories: minery, construction,
#' manufacture, retail, transport, information, finance, professional,
#' education, leisure, services, public.}
#' \item{educ}{Education attainment with 5 categories: lhs (less than high
#' school graduate, years of educ < 12), hsg (high school graduate: years of
#' educ = 12), sc (some college: 13<=years of educ<=15), cg (college:
#' 16<=years of educ<=17), ad (advanced degree: years of educ>=18).}
#' \item{ms}{Marital Status with 5 categories: married, widowed, separated,
#' divorced, and nevermarried.}
#' \item{region}{Regions with 4 categories: mw (midwest), so (south),
#' we (west), ne (northeast).}
#' }
#' @source U.S. March Supplement of the Current Population Survey (CPS) in 2015.
#'
"wage2015"
