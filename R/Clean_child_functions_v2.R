#' @importFrom ExperimentHub ExperimentHub
#' @import EpipwR.data
#' @importFrom stats cor cor.test pbeta pnorm pt qbeta qnorm rbeta rnorm
#' runif sd var


generate_null_pvals <- function(Total,lim){
  output <- numeric(lim)
  draws<-runif(lim)
  output[1] <- rbeta(1,1,Total)
  for (i in 2:lim){
    output[i] <- 1-((1-output[i-1])^(Total-i+1)-
                      draws[i]*(1-output[i-1])^(Total-i+1))^(1/(Total-i+1))
  }
  output
}

rtnorm <- function(n,mu,sig,lb,ub){
  CDF_lb <- pnorm(lb,mu,sig)
  CDF_ub <- pnorm(ub,mu,sig)
  u<-runif(n,CDF_lb,CDF_ub)
  qnorm(u,mu,sig)
}

get_cor_pval <- function(M, Y, np_method){
  suppressWarnings(cor.test(M, Y, method=np_method)$p.value)
}

get_EpipwR_data <- function(Tissue){
  hub <- ExperimentHub::ExperimentHub()
  dat <- switch(Tissue,
                     "Saliva" = hub[["EH9499"]],
                     "Lymphoma" = hub[["EH9500"]],
                     "Placenta" = hub[["EH9501"]],
                     "Liver" = hub[["EH9502"]],
                     "Colon" = hub[["EH9503"]],
                     "Blood adult" = hub[["EH9504"]],
                     "Blood 5 year olds" = hub[["EH9505"]],
                     "Blood newborns" = hub[["EH9506"]],
                     "Cord-blood (whole blood)" = hub[["EH9507"]],
                     "Cord-blood (PBMC)" = hub[["EH9508"]],
                     "Adult (PBMC)" = hub[["EH9509"]],
                     "Sperm" = hub[["EH9510"]],
                     stop("Tissue type not found")
  )
  if(Tissue == "Blood newborns"){dat <- dat[dat[,2]>0 & dat[,1]>0,]}
  return(dat)
}

how_many_null_pvals <- function(dm, Total, sig){
  l <- 1
  while(l < dm){
    if(pbeta((dm+l)*sig/Total,l,Total-dm-l+1)<(1*10^-(10))){
      break
    } else{l <- l + 1}
  }
  ###The minimum number of potential false positives to consider is 5
  ###(unless there are less than 5 non-associated tests)
  l <- min(max(l,5),Total-dm)
}

generate_non_null_pvals <- function(n, dm, par, ic, test){
  ###Start by generating correlated normal samples
  Y <- rnorm(n)
  X <- rnorm(n*dm)
  X1 <- X*sqrt(1-ic^2) + Y*ic
  ##convert X1 to the selected beta distribution using percentiles
  ds <- qbeta(pnorm(X1),rep(par[,1],rep(n,dm)),rep(par[,2],rep(n,dm)))
  df <- data.frame(M = log(ds/(1-ds)),cpg = rep(seq_len(dm),rep(n,dm)))
  if(test=="pearson"){
    ##calculate correlations
    r<-as.vector(tapply(df$M,df$cpg,cor,y=Y))
    ##calculate p-values
    2*(1-pt(abs(r/sqrt((1-r^2)/(n-2))),n-2))
  } else {
    as.vector(tapply(df$M,INDEX=df$cpg,get_cor_pval,Y=Y,np_method=test))
  }
}

find_BH_cutoff <- function(p,sig,dm,Total,l){
  ##Generate smallest "l" p-values from true null tests
  if(l==Total-dm){nsp <- runif(Total-dm)} else{
    nsp <- generate_null_pvals(Total-dm,l) }

  ##Group all p-values together
  all <- c(p,nsp)
  all <- all[order(all)]

  ##Determine the BH cutoff
  cutoffs <- sig*(seq(length(all))/Total)
  max(c(sig/Total,cutoffs[all < cutoffs]),na.rm=TRUE)
}


generate_non_null_pvals_cc <- function(dm, n1, n2, par, par2, test){
  ###Generate data
  n <- n1 + n2
  M1 <- logit(rbeta(dm*n1,par[,1],par[,2]))
  M2 <- logit(rbeta(dm*n2,par2[,1],par2[,2]))
  cpg1 <- rep(seq_len(dm),n1)
  cpg2 <- rep(seq_len(dm),n2)


  ##calculate test statistic (pooled variance and WS)
  s1 <- tapply(M1,cpg1,var)
  s2 <- tapply(M2,cpg2,var)

  if(test=="pooled"){
    TS<-(tapply(M1,cpg1,mean) - tapply(M2,cpg2,mean))/sqrt(
          (1/n1 + 1/n2)*(s1*(n1-1) + s2*(n2-1))/(n-2))
    2*pt(-1*abs(TS), n-2)

  } else if (test=="WS"){
    TS<-(tapply(M1,cpg1,mean) - tapply(M2,cpg2,mean))/sqrt(s1/n1 + s2/n2)
    df <- ((s1/n1 + s2/n2)^2)/((s1^2)/((n1-1)*n1^2) + (s2^2)/((n2-1)*n2^2))
    2*pt(-1*abs(TS),df)
  }
}

logit <- function(x){
  log(x/(1-x))
}

get_ab <- Vectorize(function(mu, var){
  a <- mu^2*(1-mu)/var - mu
  b <- a*(1-mu)/mu
  c(a,b)
}, vectorize.args = c("mu","var"))


get_power_findN <- function(dm, Total, n, sig, rho_mu, rho_sd=0, ab_sets,
                            Nmax=1000, MOE=.01, Nmin=20, test="pearson",
                            use_fdr=TRUE){

  #Start by generating the maximum number of parameters needed
  N <- Nmax
  ab<-sample(seq_len(dim(ab_sets)[1]),dm*N,TRUE)
  par1 <- ab_sets[ab,]

  ##Find a very conservative number of false-positive p-values for
  ##FDR correction
  l <- how_many_null_pvals(dm, Total, sig)

  ##Set Upper/Lower bounds if rho_sd is not 0
  if(rho_sd==0){ic<-rho_mu} else if(rho_mu>0){lb<-.03
    ub <- 1} else{lb <- -1
    ub <- -.03}
  Collect <- data.frame(power=numeric(N),cutoff=0)

  ##Begin loop: for each dataset, generate data and calculate power
  for(i in seq_len(N)){
    ##If correlation isn't fixed, generate rho ("ic")
    if(rho_sd>0){
      ic <- rtnorm(n,rho_mu,rho_sd,lb,ub)
    }
    ###Generate data and p-values for non-null tests
    rng <- ((i-1)*dm+1):(i*dm)
    par <- matrix(par1[rng,],length(rng),2)
    p <- generate_non_null_pvals(n, dm, par, ic, test)

    if(use_fdr==TRUE){
      cutoff <- find_BH_cutoff(p, sig, dm, Total, l) } else{
      cutoff <- sig/Total
    }

    ##Calculate power and store the cutoff
    Collect[i,1] <- length(p[p<cutoff])/dm
    Collect[i,2] <- cutoff

    ##breaks the loop if N is above 20 and if the 95% CI for power has a
    ##small enough MOE
    if(i>=Nmin && sd(Collect$power[seq_len(i)])<(sqrt(i)*MOE/1.96)){break}
  }
  Collect[Collect$cutoff>0,]
}



get_power_cc_findN <- function(dm, Total, n, sig, delta_mu, delta_sd=0,
                               ab_sets, test="pooled", n1_prop=0.5, Nmax=1000,
                               MOE=.01, Nmin=20, use_fdr=TRUE){
  #Start by generating the maximum number of parameters needed
  N <- Nmax
  ab<-sample(seq_len(dim(ab_sets)[1]),dm*N,TRUE)
  par1 <- ab_sets[ab,]
  ##Find a very conservative number of false-positive p-values required
  l <- how_many_null_pvals(dm, Total, sig)
  Collect <- data.frame(power=numeric(N),cutoff=0)
  ##Begin loop: for each dataset, generate data and calculate power
  for(i in seq_len(N)){
    ####Establish sample sizes
    n1 <- round(n*n1_prop)
    n2 <- n-n1
    ###Generate data and p-values for non-null tests
    rng <- ((i-1)*dm+1):(i*dm)
    par <- matrix(par1[rng,],length(rng),2)
    mu <- par[,1]/(par[,1] + par[,2])
    var <- mu*(1-mu)/(par[,1] + par[,2] + 1)
    ##FOR NOW: FINISH ONLY FOR DELTA_SD > 0
    if(delta_sd>0){
      lb <- ifelse(delta_mu>0,0,-0.5)
      ub <- ifelse(delta_mu>0,0.5,0)
      del <- rtnorm(dm,delta_mu*ifelse(delta_mu > lb & delta_mu<ub,1,-1),
                    delta_sd,lb,ub)
    } else{ del <- delta_mu}

    mu_new <- ifelse(mu+del<1 & mu+del>0,mu+del,mu-del)
    par2 <- matrix(0,length(rng),2)
    par2[,1] <- mu_new*apply(par,1,sum)
    par2[,2] <- apply(par,1,sum) - par2[,1]
    p <- generate_non_null_pvals_cc(dm, n1, n2, par, par2, test)
    if(use_fdr==TRUE){
        cutoff <- find_BH_cutoff(p, sig, dm, Total, l) } else{
        cutoff <- sig/Total
      }
    ##Calculate power and store the cutoff
    Collect[i,1] <- length(p[p<cutoff])/dm
    Collect[i,2] <- cutoff
    ##breaks the loop if N is above 20 and if the 95% CI for power has
    ##a small enough MOE
    if(i>=Nmin && sd(Collect$power[seq_len(i)])<(sqrt(i)*MOE/1.96)){break}
  }
  Collect[Collect$cutoff>0,]
}

gpcont_check <- function(dm, Total, n, fdr_fwer, rho_mu, rho_sd=0,
                           Tissue="Saliva", Nmax=1000, MOE=.03, test="pearson",
                           use_fdr=TRUE, Suppress_updates=FALSE){
  tiss <- c("Saliva", "Lymphoma", "Placenta", "Liver", "Colon", "Blood adult",
            "Blood 5 year olds",  "Blood newborns", "Cord-blood (whole blood)",
            "Cord-blood (PBMC)", "Adult (PBMC)", "Sperm")
  if(!is.numeric(dm) || length(dm)>1 || dm < 1 || floor(dm)!= dm){stop(
    "The number of significant associations (dm) must be a single positive
    integer")}
  if(!is.numeric(Total) || length(Total)>1 || Total < 1 ||
     floor(Total)!= Total){stop(
       "The total numnber of CpG sites (Total) must be a single positive
       integer")
    }
  for(i in seq_len(length(n))){
    if(!is.numeric(n[i]) || n[i] < 2 || floor(n[i])!= n[i]){stop(
      "All sample sizes (n) must be positive integers greater than 1")}
  }
  if(dm > Total){stop("dm must be less than or equal to N")}
  if(!is.numeric(fdr_fwer) || length(fdr_fwer)>1 || fdr_fwer <= 0 ||
     fdr_fwer >= 1){stop("fdr_fwer must be a single number in (0,1)")}
  for (i in seq_len(length(rho_mu))){
    if(!is.numeric(rho_mu[i]) || abs(rho_mu[i]) > 1 || abs(rho_mu[i]) <= 0){
      stop("All correlations (rho_mu) must have absolute value in (0,1]")}
    if(rho_mu[i] < .05 && rho_mu[i]>-.05){warning(
      "|rho_mu| < .05 will be extremely difficult to detect")}
  }
  if(!is.numeric(rho_sd) || length(rho_sd)>1 || rho_sd<0){stop(
    "rho_sd must be a single non-negative number"
  )}
  if(!(Tissue %in% tiss)){stop(
    "Unknown Tissue selected. Valid options are Saliva, Lymphoma, Placenta,
     Liver, Colon, Blood adult, Blood 5 year olds, Blood newborns,
     Cord-blood (whole blood), Cord-blood (PBMC), Adult (PBMC), Sperm"
  )}
  if(!is.numeric(Nmax) || length(Nmax)>1 || Nmax<20 || floor(Nmax)!=Nmax){
    stop("Nmax must be a positive integer greater than or equal to 20")
    }
  if(!is.numeric(MOE) || length(MOE) > 1 || MOE<0 || MOE>1){stop(
    "MOE must be a single number in [0,1]"
  )}
  if(!(test %in% c("pearson","kendall","spearman"))){stop(
    "Invalid test selected. Valid options are 'pearson', 'spearman' or
    'kendall'"
  )}
  if(!identical(use_fdr,TRUE) && !identical(use_fdr,FALSE)){stop(
    "use_fdr must be boolean (either TRUE or FALSE)"
  )}
}

gpcc_check1 <- function(dm, Total, n, fdr_fwer, delta_mu, delta_sd=0,
                       n1_prop=0.5, Tissue="Saliva", Nmax=1000, MOE=.03,
                       test="pooled", use_fdr=TRUE, Suppress_updates=FALSE){

  if(!is.numeric(dm) || length(dm)>1 || dm < 1 || floor(dm)!= dm){stop(
    "The number of significant associations (dm) must be a single positive
    integer"
  )}
  if(!is.numeric(Total) || length(Total)>1 || Total < 1 ||
     floor(Total)!= Total){stop(
    "The total number of CpG sites (Total) must be a single positive integer"
  )}
  if(dm > Total){stop("dm must be less than or equal to N")}
  if(!is.numeric(fdr_fwer) || length(fdr_fwer)>1 || fdr_fwer <= 0 ||
     fdr_fwer >= 1){stop("fdr_fwer must be a single number in (0,1)")}
  if(!is.numeric(delta_sd) || length(delta_sd)>1 || delta_sd<0){stop(
    "delta_sd must be a single non-negative number"
  )}
  if(!is.numeric(Nmax) || length(Nmax)>1 || Nmax<20 || floor(Nmax)!=Nmax){stop(
    "Nmax must be a positive integer greater than or equal to 20")}
  if(!is.numeric(MOE) || length(MOE) > 1 || MOE<0 || MOE>1){stop(
    "MOE must be a single number in [0,1]"
  )}
  if(!(test %in% c("pooled","WS"))){stop(
    "Invalid test selected. Valid options are 'pooled' or 'WS'
    (Welch-Satterthwaite)"
  )}
  if(!identical(use_fdr,TRUE) && !identical(use_fdr,FALSE)){stop(
    "use_fdr must be boolean (either TRUE or FALSE)"
  )} }

gpcc_check2 <- function(dm, Total, n, fdr_fwer, delta_mu, delta_sd=0,
                        n1_prop=0.5, Tissue="Saliva", Nmax=1000, MOE=.03,
                        test="pooled", use_fdr=TRUE, Suppress_updates=FALSE){
  n1 <- round(n*n1_prop)
  n2 <- n - n1
  tiss <- c("Saliva", "Lymphoma", "Placenta", "Liver", "Colon", "Blood adult",
            "Blood 5 year olds",  "Blood newborns", "Cord-blood (whole blood)",
            "Cord-blood (PBMC)", "Adult (PBMC)", "Sperm")

  if(!is.numeric(n1_prop) || length(n1_prop)>1 || n1_prop >= 1 ||
     n1_prop <= 0){stop(
       "The proportion of group 1 (n1_prop) must be a single value
       between 0 and 1"
     )}

  for(i in seq_len(length(n))){
    if(!is.numeric(n[i]) || n[i] < 4 || floor(n[i])!= n[i]){stop(
      "All sample sizes (n) must be positive integers greater than 3"
    )} else if(n1[i] < 2 || n2[i] < 2){stop(
      "All sample sizes (n) must be large enough so that each group has at
      least 2 subjects"
    )}
  }

  for (i in seq_len(length(delta_mu))){
    if(!is.numeric(delta_mu[i]) || abs(delta_mu[i]) > 1 ||
       abs(delta_mu[i]) <= 0){stop(
         "All correlations (rho_mu) must have absolute value in (0,1]"
       )}
    if(delta_mu[i] < .001 && delta_mu[i]>-.01){warning(
      "|delta_mu| < .001 will be extremely difficult to detect"
    )}
  }

  if(!(Tissue %in% tiss)){stop(
    "Unknown Tissue selected. Valid options are Saliva, Lymphoma, Placenta,
     Liver, Colon, Blood adult, Blood 5 year olds, Blood newborns,
     Cord-blood (whole blood), Cord-blood (PBMC), Adult (PBMC), Sperm"
  )} }
