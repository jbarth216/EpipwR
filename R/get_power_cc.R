####Parent function (one or more sample sizes)



#
#' @title Power Calculations for Case/Control EWAS
#' @description
#' Calculates power for EWAS with a binary outcome for multiple
#' sample sizes and/or effect sizes. Data sets are only simulated for the
#' non-null tests; p-values are generated directly for the null tests.
#' Rather than specifying the number of data sets to calculate power, you specify the precision level (`MOE`)
#' and a maximum number of data sets (`Nmax`). After 20 data sets, the
#' function terminates when the desired precision level is reached or if the number
#' of tested data sets reaches `Nmax`.
#'
#' @param dm Number of non-null tests.
#' @param Total The total number of tests (null and non-null).
#' @param n Sample size(s) for which power is calculated (accepts a vector). This is total sample size.
#' @param fdr_fwer Either the false discovery rate or the family-wise type I error rate, depending on `use_fdr`.
#' @param delta_mu Average effect size (difference in mean methylation between groups) for non-null tests (accepts a vector).
#' If multiple `n` and `delta_mu` are specified, power is calculated for all unique settings under the Cartesian product of these vectors.
#' @param delta_sd Standard deviation of the effect size for non-null tests. If 0, all effect sizes are fixed at `delta_mu`.
#' @param n1_prop Indicates the proportion of the total sample size (`n`) in group 1 (rounded to the nearest integer).
#' @param Tissue Tissue type of Empirical EWAS to be used for data generation (see details for valid options).
#' @param Nmax The maximum number of data sets used to calculate power.
#' @param MOE The target margin of error of a 95% confidence interval for average power. This determines the stopping point of the algorithm.
#' @param test The type of t-test to be used. `"pooled"` indicates a pooled variance t-test while `"WS"` indicates Welch's t-test.
#' @param use_fdr If `TRUE`, uses `fdr_fwer` as the false discovery rate. If `FALSE`, uses the family-wise type I error rate.
#' @param Suppress_updates If `TRUE`, blocks messages reporting the completion of each unique setting.
#' @details
#' Valid options for the `Tissue` argument are `"Saliva"`, `"Lymphoma"`, `"Placenta"`, `"Liver"`,
#' `"Colon"`, `"Blood adult"`, `"Blood 5 year olds"`, `"Blood newborns"`, `"Cord-blood (whole blood)"`,
#' `"Cord-blood (PBMC)"`, `"Adult (PBMC)"`, and `"Sperm"`. All data sets are publicly available on the gene
#' expression omnibus (see \href{https://github.com/jbarth216/EpipwR.data}{EpipwR.data} package for more details) and were identified by Graw, et. al. (2019).
#' Please note that, due to some extreme values in this data set, the Lymphoma option will occasionally throw
#' a warning related to data generation. At this time, we recommend using one of the other tissue options.
#'
#' Unlike pwrEWAS (Graw et al., 2019), EpipwR enforces equality of precision (sum of the parameters) on the distributions
#' of each group rather than equality of variance.
#' @references Graw, S., Henn, R., Thompson, J. A., and Koestler, D. C. (2019). pwrEWAS:
#' A user-friendly tool for comprehensive power estimation for epigenome wide
#' association studies (EWAS). \emph{BMC Bioinformatics}, 20(1):218.
#'
#' @return A dataframe with rows equal to the number of `n` and `rho_mu` combinations
#' @examples
#' # This examples calculates power for 100 non-null tests out of 10,000 total with an FDR of 5%.
#' # Sample sizes of 40, 50, 60, 70 and fixed effect sizes of 0.01, 0.02, 0.05 are used.
#' # For improved accuracy, MOE=.01 to ensure that the 95% confidence interval for average power
#' # has a margin of error no larger than .01 (unless the maximum of 1,000 data sets is exceeded).
#' get_power_cc(100,10000,c(40,50,60,70),.05,c(.01,.02,.05), MOE=.01)
#' @export
get_power_cc <- function(dm, Total, n, fdr_fwer, delta_mu, delta_sd=0, n1_prop=0.5, Tissue="Saliva", Nmax=1000, MOE=.03, test="pooled", use_fdr=TRUE,
                      Suppress_updates=FALSE){
  n1 <- round(n*n1_prop)
  n2 <- n - n1
  tiss <- c("Saliva",
            "Lymphoma",
            "Placenta",
            "Liver",
            "Colon",
            "Blood adult",
            "Blood 5 year olds",
            "Blood newborns",
            "Cord-blood (whole blood)",
            "Cord-blood (PBMC)",
            "Adult (PBMC)",
            "Sperm")

  if(!is.numeric(dm) || length(dm)>1 || dm < 1 || floor(dm)!= dm){stop(
    "The number of significant associations (dm) must be a single positive integer"
    )}
  if(!is.numeric(Total) || length(Total)>1 || Total < 1 || floor(Total)!= Total){stop(
    "The total number of CpG sites (Total) must be a single positive integer"
  )}
  #if(!is.numeric(n1_prop) || length(n1_prop)>1 || n1_prop >= 1 || n1_prop <= 0){stop(
  #  "The proportion of group 1 (n1_prop) must be a single value between 0 and 1"
  #)}

  for(i in 1:length(n)){
  if(!is.numeric(n[i]) || n[i] < 4 || floor(n[i])!= n[i]){stop(
    "All sample sizes (n) must be positive integers greater than 3"
  )} else if(n1[i] < 2 || n2[i] < 2){stop(
    "All sample sizes (n) must be large enough so that each group has at least 2 subjects"
  )}
  }

  if(dm > Total){stop(
    "dm must be less than or equal to N"
    )}
  if(!is.numeric(fdr_fwer) || length(fdr_fwer)>1 || fdr_fwer <= 0 || fdr_fwer >= 1){stop(
    "fdr_fwer must be a single number in (0,1)"
  )}
  for (i in 1:length(delta_mu)){
  if(!is.numeric(delta_mu[i]) || abs(delta_mu[i]) > 1 || abs(delta_mu[i]) <= 0){stop(
    "All correlations (rho_mu) must have absolute value in (0,1]"
  )}
    if(delta_mu[i] < .001 && delta_mu[i]>-.01){warning(
      "|delta_mu| < .001 will be extremely difficult to detect"
    )}
    }
  if(!is.numeric(delta_sd) || length(delta_sd)>1 || delta_sd<0){stop(
    "delta_sd must be a single non-negative number"
    )}
  if(!(Tissue %in% tiss)){stop(
    "Unknown Tissue selected. Valid options are Saliva, Lymphoma, Placenta,
     Liver, Colon, Blood adult, Blood 5 year olds, Blood newborns,
     Cord-blood (whole blood), Cord-blood (PBMC), Adult (PBMC), Sperm"
    )}
  if(!is.numeric(Nmax) || length(Nmax)>1 || Nmax<20 || floor(Nmax)!=Nmax){stop(
    "Nmax must be a positive integer greater than or equal to 20")}
  if(!is.numeric(MOE) || length(MOE) > 1 || MOE<0 || MOE>1){stop(
    "MOE must be a single number in [0,1]"
    )}
  if(!(test %in% c("pooled","WS"))){stop(
    "Invalid test selected. Valid options are 'pooled' or 'WS' (Welch-Satterthwaite)"
    )}

  if(!identical(use_fdr,TRUE) && !identical(use_fdr,FALSE)){stop(
    "use_fdr must be boolean (either TRUE or FALSE)"
  )}



  ###Start Power Analysis
  n <- sort(n)
  delta <- sort(delta_mu)
  runs <- length(n)*length(delta)
  out <- expand.grid(n, delta)
  colnames(out) <- c("n","delta")
  out$avg_power <- 0
  out$sd_power <- 0
  out$N <- 0

  ##Suppress EH message
  suppressMessages(ab_sets <-get_EpipwR_data(Tissue))

    for(i in 1:runs){
      Results <- get_power_cc_findN(dm, Total, out$n[i], fdr_fwer, out$delta[i], delta_sd, ab_sets, test, n1_prop, Nmax, MOE, Nmin=20, use_fdr)
      out$avg_power[i] <- mean(Results$power)
      out$sd_power[i] <- sd(Results$power)
      out$N[i] <- nrow(Results)
      if(Suppress_updates==F){
        print(paste0(i," out of ",runs," settings completed"))
      }
    }

  out$se_power <- out$sd_power/sqrt(out$N)
  colnames(out)[1:2] <- c("sample_size","delta_mu")
  out
}

