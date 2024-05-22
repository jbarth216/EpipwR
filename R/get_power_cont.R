####Parent function (one or more sample sizes)
#source("Clean_child_functions.R")


#
#' @title Power Calculations for Continuous EWAS
#' @description
#' Calculates power for EWAS with a continuous outcome for multiple
#' sample sizes and/or correlations. Data sets are only simulated for the
#' non-null tests; p-values are generated directly for the null tests.
#' Rather than specifying the number of data sets to calculate power, you specify the precision level (`MOE`)
#' and a maximum number of data sets (`Nmax`). After 20 data sets, the
#' function terminates when the desired precision level is reached or if the number
#' of tested data sets reaches `Nmax`.
#'
#' @param dm Amount of non-null tests.
#' @param Total The total number of tests (null and non-null).
#' @param n Sample size(s) for which power is calculated (accepts a vector).
#' @param fdr_fwer Either the false discovery rate or the family-wise type I error rate, depending on `use_fdr`.
#' @param rho_mu Mean correlation(s) of methylation and phenotype for non-null tests (accepts a vector). If multiple `n` and `rho_mu` are specified, power is calculated for all unique settings under the Cartesian product of these vectors.
#' @param rho_sd Standard deviation of methylation and phenotype for non-null tests. If 0 all correlations are fixed at `rho_mu`.
#' @param Tissue Tissue type of Empirical EWAS to be used for data generation (see details for valid options).
#' @param Nmax The maximum number of datasets used to calculate power
#' @param MOE The target margin of error of a 95% confidence interval for average power. This determines the stopping point of the algorithm.
#' @param test The type of statistical test to be used. `"pearson"`, `"kendall"` and `"spearman"` are valid options.
#' @param use_fdr If `TRUE`, uses `fdr_fwer` as the false discovery rate. If `FALSE`, uses the family-wise type I error rate.
#' @param Suppress_updates If `TRUE`, blocks messages reporting the completion of each unique setting.
#' @details
#' Valid options for the `Tissue` argument are `"Saliva"`, `"Lymphoma"`, `"Placenta"`, `"Liver"`,
#' `"Colon"`, `"Blood adult"`, `"Blood 5 year olds"`, `"Blood newborns"`, `"Cord-blood (whole blood)"`,
#' `"Cord-blood (PBMC)"`, `"Adult (PBMC)"`, and `"Sperm"`. All data sets are publicly available on the gene
#' expression omnibus (see [Epipwr.data] package for more details) and were identified by Graw, et. al. (2019).
#' @references Graw, S., Henn, R., Thompson, J. A., and Koestler, D. C. (2019). pwrEWAS:
#' A user-friendly tool for comprehensive power estimation for epigenome wide
#' association studies (EWAS). \emph{BMC Bioinformatics}, 20(1):218.
#'
#' @return A dataframe with rows equal to the number of `n` and `rho_mu` combinations
#' @examples
#' # This examples calculates power for 100 non-null tests out of 10,000 total with an FDR of 5%.
#' # Sample sizes of 70,80,90,100 and fixed correlations at 0.3,0.4,0.5 are used.
#' get_power_cont(100,10000,c(70,80,90,100),.05,c(.3,.4,.5))
#' @export
get_power_cont <- function(dm, Total, n, fdr_fwer, rho_mu, rho_sd=0, Tissue="Saliva", Nmax=1000, MOE=.03, test="pearson", use_fdr=TRUE,
                      Suppress_updates=FALSE){
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
    "The total numnber of CpG sites (Total) must be a single positive integer"
  )}
  for(i in 1:length(n)){
  if(!is.numeric(n[i]) || n[i] < 2 || floor(n[i])!= n[i]){stop(
    "All sample sizes (n) must be positive integers greater than 1"
  )}
  }
  if(dm > Total){stop(
    "dm must be less than or equal to N"
    )}
  if(!is.numeric(fdr_fwer) || length(fdr_fwer)>1 || fdr_fwer <= 0 || fdr_fwer >= 1){stop(
    "fdr_fwer must be a single number in (0,1)"
  )}
  for (i in 1:length(rho_mu)){
  if(!is.numeric(rho_mu[i]) || abs(rho_mu[i]) > 1 || abs(rho_mu[i]) <= 0){stop(
    "All correlations (rho_mu) must have absolute value in (0,1]"
  )}
    if(rho_mu[i] < .05 && rho_mu[i]>-.05){warning(
      "|rho_mu| < .05 will be extremely difficult to detect"
    )}
    }
  if(!is.numeric(rho_sd) || length(rho_sd)>1 || rho_sd<0){stop(
    "rho_sd must be a single non-negative number"
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
  if(!(test %in% c("pearson","kendall","spearman"))){stop(
    "Invalid test selected. Valid options are 'pearson', 'spearman' or 'kendall'"
    )}

  if(!identical(use_fdr,TRUE) && !identical(use_fdr,FALSE)){stop(
    "use_fdr must be boolean (either TRUE or FALSE)"
  )}



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
  suppressMessages(ab_sets <-get_EpipwR_data(Tissue))

  if(test=="pearson" & use_fdr){
    for(i in 1:runs){
      Results <- get_power_findN(dm, Total, out$n[i], fdr_fwer, out$rho[i], rho_sd, ab_sets, Nmax, MOE)
      out$avg_power[i] <- mean(Results$power)
      out$sd_power[i] <- sd(Results$power)
      out$N[i] <- nrow(Results)
      if(Suppress_updates==F){
        print(paste0(i," out of ",runs," settings completed"))
      }
    } } else if(test=="pearson" & !(use_fdr)) {
    for(i in 1:runs){
      Results <- get_power_findN_bf(dm, Total, out$n[i], fdr_fwer, out$rho[i], rho_sd, ab_sets, Nmax, MOE)
      out$avg_power[i] <- mean(Results$power)
      out$sd_power[i] <- sd(Results$power)
      out$N[i] <- nrow(Results)
      if(Suppress_updates==F){
        print(paste0(i," out of ",runs," settings completed"))
      }
    }
  } else if(test!="pearson" & use_fdr){
    for(i in 1:runs){
      Results <- get_power_np(dm, Total, out$n[i], fdr_fwer, out$rho[i], rho_sd, ab_sets, Nmax, MOE, Nmin=20, test)
      out$avg_power[i] <- mean(Results$power)
      out$sd_power[i] <- sd(Results$power)
      out$N[i] <- nrow(Results)
      if(Suppress_updates==F){
        print(paste0(i," out of ",runs," settings completed"))
      }
    }} else{
    for(i in 1:runs){
      Results <- get_power_np_bf(dm, Total, out$n[i], fdr_fwer, out$rho[i], rho_sd, ab_sets, Nmax, MOE, Nmin=20, test)
      out$avg_power[i] <- mean(Results$power)
      out$sd_power[i] <- sd(Results$power)
      out$N[i] <- nrow(Results)
      if(Suppress_updates==F){
        print(paste0(i," out of ",runs," settings completed"))
      }
    }
  }

  out$se_power <- out$sd_power/sqrt(out$N)
  colnames(out)[1:2] <- c("sample_size","rho_mu")
  out
}

