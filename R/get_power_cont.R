####Parent function (one or more sample sizes)
#source("Clean_child_functions.R")


#
#' @title Power Calculations for Continuous EWAS
#' @description
#' Calculates power for EWAS with a continuous outcome for multiple
#' sample sizes and/or correlations. Data sets are only simulated for the
#' non-null tests; p-values are generated directly for the null tests.
#' Rather than specifying the number of data sets to calculate power, you
#' specify the precision level (`MOE`)
#' and a maximum number of data sets (`Nmax`). After 20 data sets, the
#' function terminates when the desired precision level is reached or if the
#' number of tested data sets reaches `Nmax`.
#'
#' @param dm Number of non-null tests.
#' @param Total The total number of tests (null and non-null).
#' @param n Sample size(s) for which power is calculated (accepts a vector).
#' @param fdr_fwer Either the false discovery rate or the family-wise type I
#' error rate, depending on `use_fdr`.
#' @param rho_mu Mean correlation(s) of methylation and phenotype for non-null
#' tests (accepts a vector). If multiple `n` and `rho_mu` are specified, power
#' is calculated for all unique settings under the Cartesian product of these
#' vectors.
#' @param rho_sd Standard deviation of methylation and phenotype for non-null
#' tests. If 0 all correlations are fixed at `rho_mu`.
#' @param Tissue Tissue type of Empirical EWAS to be used for data generation
#' (see details for valid options).
#' @param Nmax The maximum number of datasets used to calculate power
#' @param MOE The target margin of error of a 95% confidence interval
#' for average power. This determines the stopping point of the algorithm.
#' @param test The type of statistical test to be used. `"pearson"`,
#' `"kendall"` and `"spearman"` are valid options.
#' @param use_fdr If `TRUE`, uses `fdr_fwer` as the false discovery rate.
#' If `FALSE`, uses the family-wise type I error rate.
#' @param Suppress_updates If `TRUE`, blocks messages reporting the completion
#' of each unique setting.
#' @param emp_data Reference data set in matrix or data frame format (Beta
#' values with CpG sites as rows, samples as columns). Ignored unless
#' `Tissue="Custom"`.
#' @param phenotype_data A sample of phenotype data to be used in the power
#' calculation(s). Accepts a vector with length > 100, although we recommend
#' that users make this at least as large as their maximum sample size.
#' If left blank, a normal distribution is used to generate correlations.
#' @details
#' Valid options for the `Tissue` argument are `"Saliva"`, `"Lymphoma"`,
#' `"Placenta"`, `"Liver"`, `"Colon"`, `"Blood adult"`, `"Blood 5 year olds"`,
#' `"Blood newborns"`, `"Cord-blood (whole blood)"`, `"Cord-blood (PBMC)"`,
#' `"Adult (PBMC)"`, and `"Sperm"`. All data sets are publicly available on the
#' gene expression omnibus (see
#' \href{https://github.com/jbarth216/EpipwR.data}{EpipwR.data} package for
#' more details) and were identified by Graw, et. al. (2019). Please note that,
#' due to some extreme values in this data, the Lymphoma option will
#' occasionally throw a warning related to data generation. At this time, we
#' recommend using one of the other tissue options. Users who would like to use
#' their own reference data set should set `Tissue="Custom"` and provide the
#' data in matrix or data frame format with `emp_data`. Similarly, users can
#' also now specify their own phenotype data, either from a real data set or by
#' generating samples from a known distribution (i.e., `rt(1000,2)`). Users who
#' would like to take advantage of either of these settings are responsible for
#' the quality and formatting of the data provided.
#'
#' Although this function only covers 3 types of statistical tests, its worth
#' noting that tests run using software packages such as limma will yield the
#' same results as a pearson correlation test in the absence of covariates
#' or any dependence across CpG sites (as is the assumption here).
#' Any users wanting to mimic an analysis done in limma should use
#' `test="pearson"`.
#'
#'
#' @references Graw, S., Henn, R., Thompson, J. A., and Koestler, D. C. (2019).
#' pwrEWAS: A user-friendly tool for comprehensive power estimation for
#' epigenome wide association studies (EWAS). \emph{BMC Bioinformatics},
#' 20(1):218.
#'
#' @return A data frame with rows equal to the number of `n` and `rho_mu`
#' combinations
#' @examples
#' # This examples calculates power for 100 non-null tests out of 10,000 total
#' # with an FDR of 5%.
#' # Sample sizes of 70,80,90,100 and fixed correlations at 0.3,0.4,0.5 are
#' # used.
#' get_power_cont(100,10000,c(70,80,90,100),.05,c(.3,.4,.5))
#' @export
get_power_cont <- function(dm, Total, n, fdr_fwer, rho_mu, rho_sd=0,
                           Tissue="Saliva", Nmax=1000, MOE=.03, test="pearson",
                           use_fdr=TRUE, Suppress_updates=FALSE, emp_data=NULL,
                           phenotype_data=NULL){
 ##Check that all inputs are valid
 gpcont_check(dm,Total,n,fdr_fwer,rho_mu,rho_sd,Tissue,Nmax,MOE,test,use_fdr,
              Suppress_updates, emp_data, phenotype_data)

  ###Start Power Analysis
  n <- sort(n)
  rho <- sort(rho_mu)
  runs <- length(n)*length(rho)
  out <- expand.grid(n, rho)
  colnames(out) <- c("n","rho")
  out$avg_power <- 0
  out$sd_power <- 0
  out$N <- 0

  ##Suppress EH message
  if(Tissue == "Custom"){ab_sets <- abfinder(emp_data)} else{
    suppressMessages(ab_sets <-get_EpipwR_data(Tissue))
  }

  for(i in seq_len(runs)){
    Results <- get_power_findN(dm, Total, out$n[i], fdr_fwer, out$rho[i],
                               rho_sd, ab_sets, Nmax, MOE, Nmin=20, test,
                               use_fdr, phenotype_data)
    out$avg_power[i] <- mean(Results$power)
    out$sd_power[i] <- sd(Results$power)
    out$N[i] <- nrow(Results)
    if(Suppress_updates==FALSE){
      message(c(i," out of ",runs," settings completed"))
    }
  }
  out$se_power <- out$sd_power/sqrt(out$N)
  colnames(out)[c(1,2)] <- c("sample_size","rho_mu")
  out
}

