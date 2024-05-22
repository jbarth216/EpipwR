###Clean_child_functions
#' @importFrom ExperimentHub ExperimentHub
#' @import EpipwR.data


beta_sampling_cond <- function(Total,lim){
  output <- numeric(lim)
  draws<-runif(lim)
  output[1] <- rbeta(1,1,Total)
  for (i in 2:lim){
    output[i] <- 1-((1-output[i-1])^(Total-i+1)-draws[i]*(1-output[i-1])^(Total-i+1))^(1/(Total-i+1))
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
  if(Tissue == "Blood newborns"){dat <- dat[dat[,2]>0,]}
  return(dat)
}


get_power_findN <- function(dm, Total, n, sig, rho_mu, rho_sd=0, ab_sets, Nmax=1000, MOE=.01, Nmin=20){

  #Start by generating the maximum number of parameters needed
  N <- Nmax
  ab<-sample(1:dim(ab_sets)[1],dm*N,T)
  par1 <- ab_sets[ab,]

  ##Find a very conservative number of false-positive p-values required
  l <- 1
  while(l < dm){
    if(pbeta((dm+l)*sig/Total,l,Total-dm-l+1)<(1*10^-(10))){
      break
    } else{l <- l + 1}
  }
  ###The minimum number of potential false positives to consider is 5
  l <- max(l,5)

  ##Set Upper/Lower bounds if rho_sd is not 0
  if(rho_sd==0){ic<-rho_mu} else if(rho_mu>0){lb<-.03
    ub <- 1} else{lb <- -1
    ub <- -.03}
  Collect <- data.frame(power=numeric(N),cutoff=0)

  ##Begin loop: for each dataset, generate data and calculate power
  for(i in 1:N){
    ##If correlation isn't fixed, generate rho ("ic")
    if(rho_sd>0){
      ic <- rtnorm(n,rho_mu,rho_sd,lb,ub)
    }

    ###Generate data and p-values for non-null tests
      rng <- ((i-1)*dm+1):(i*dm)
      par<-matrix(par1[rng,],length(rng),2)
      ###Start by generating correlated normal samples
      Y <- rnorm(n)
      X <- rnorm(n*dm)
      X1 <- X*sqrt(1-ic^2) + Y*ic
      ##convert X1 to the selected beta distribution using percentiles
      ds <- qbeta(pnorm(X1),rep(par[,1],rep(n,dm)),rep(par[,2],rep(n,dm)))
      df <- data.frame(M = log(ds/(1-ds)),cpg = rep(1:dm,rep(n,dm)))
      ##calculate correlations
      r<-as.vector(tapply(df$M,df$cpg,cor,y=Y))
      ##calculate p-values
      p <- 2*(1-pt(abs(r/sqrt((1-r^2)/(n-2))),n-2))

    ##Generate smallest "l" p-values from true null tests
    nsp <- beta_sampling_cond(Total-dm,l)

    ##Group all p-values together
    all <- c(p,nsp)
    all <- all[order(all)]

    ##Determine the BH cutoff
    cutoffs <- sig*(seq(length(all))/Total)
    use_this <- max(c(sig/Total,cutoffs[all < cutoffs]),na.rm=T)

    ##Calculate power and store the cutoff
    Collect[i,1] <- length(p[p<use_this])/dm
    Collect[i,2] <- use_this

    ##breaks the loop if N is above 20 and if the 95% CI for power has a small enough MOE
    if(i>=Nmin && sd(Collect$power[1:i])<(sqrt(i)*MOE/1.96)){break}
  }
  Collect[Collect$cutoff>0,]
}

get_power_np <- function(dm, Total, n, sig, rho_mu, rho_sd=0, ab_sets, Nmax=1000, MOE=.01, Nmin=20, np_method="spearman"){
  N <- Nmax
  ab<-sample(1:dim(ab_sets)[1],dm*N,T)
  par1 <- ab_sets[ab,]

  ##Find max number of false-positive pvals needed
  l <- 1
  while(l < dm){
    if(pbeta((dm+l)*sig/Total,l,Total-dm-l+1)<(1*10^-(10))){
      break
    } else{l <- l + 1}
  }
  l <- max(l,5)
  if(rho_sd==0){ic<-rho_mu} else if(rho_mu>0){lb<-.03
  ub <- 1} else{lb <- -1
  ub <- -.03}
  Collect <- data.frame(power=numeric(N),cutoff=0)
  for(i in 1:N){
    if(rho_sd>0){
      ic <- rtnorm(n,rho_mu,rho_sd,lb,ub)
    }
    rng <- ((i-1)*dm+1):(i*dm)
    par<-matrix(par1[rng,],length(rng),2)
    Y <- rnorm(n)
    X <- rnorm(n*dm)
    X1 <- X*sqrt(1-ic^2) + Y*ic
    ds <- qbeta(pnorm(X1),rep(par[,1],rep(n,dm)),rep(par[,2],rep(n,dm)))
    df <- data.frame(M = log(ds/(1-ds)),cpg = rep(1:dm,rep(n,dm)))

    p<-as.vector(tapply(df$M,INDEX=df$cpg,get_cor_pval,Y=Y,np_method=np_method))
    nsp <- beta_sampling_cond(Total-dm,l)
    all <- c(p,nsp)
    all <- all[order(all)]
    cutoffs <- sig*(seq(length(all))/Total)
    use_this <- max(c(sig/Total,cutoffs[all < cutoffs]),na.rm=T)
    Collect[i,1] <- length(p[p<use_this])/dm
    Collect[i,2] <- use_this
    if(i>=Nmin && sd(Collect$power[1:i])<(sqrt(i)*MOE/1.96)){break}
  }
  Collect[Collect$cutoff>0,]
}


get_power_findN_bf <- function(dm, Total, n, sig, rho_mu, rho_sd=0, ab_sets, Nmax=1000, MOE=.01, Nmin=20){
  N <- Nmax
  ab<-sample(1:dim(ab_sets)[1],dm*N,T)
  par1 <- ab_sets[ab,]

  if(rho_sd==0){ic<-rho_mu} else if(rho_mu>0){lb<-.03
  ub <- 1} else{lb <- -1
  ub <- -.03}
  Collect <- data.frame(power=numeric(N),cutoff=0)
  for(i in 1:N){
    if(rho_sd>0){
      ic <- rtnorm(n,rho_mu,rho_sd,lb,ub)
    }
    rng <- ((i-1)*dm+1):(i*dm)
    par<-matrix(par1[rng,],length(rng),2)##Add this to other functions so dm=1 is ok
    Y <- rnorm(n)
    X <- rnorm(n*dm)
    X1 <- X*sqrt(1-ic^2) + Y*ic
    ds <- qbeta(pnorm(X1),rep(par[,1],rep(n,dm)),rep(par[,2],rep(n,dm)))
    df <- data.frame(M = log(ds/(1-ds)),cpg = rep(1:dm,rep(n,dm)))

    r<-as.vector(tapply(df$M,df$cpg,cor,y=Y))
    p <- 2*(1-pt(abs(r/sqrt((1-r^2)/(n-2))),n-2))
    use_this <- sig/Total
    Collect[i,1] <- length(p[p<use_this])/dm
    Collect[i,2] <- use_this
    if(i>=Nmin && sd(Collect$power[1:i])<(sqrt(i)*MOE/1.96)){break}
  }
  Collect[Collect$cutoff>0,]
}

get_power_np_bf <- function(dm, Total, n, sig, rho_mu, rho_sd=0, ab_sets, Nmax=1000, MOE=.01, Nmin=20, np_method="spearman"){
  N <- Nmax
  ab<-sample(1:dim(ab_sets)[1],dm*N,T)
  par1 <- ab_sets[ab,]

  if(rho_sd==0){ic<-rho_mu} else if(rho_mu>0){lb<-.03
  ub <- 1} else{lb <- -1
  ub <- -.03}
  Collect <- data.frame(power=numeric(N),cutoff=0)
  for(i in 1:N){
    if(rho_sd>0){
      ic <- rtnorm(n,rho_mu,rho_sd,lb,ub)
    }
    rng <- ((i-1)*dm+1):(i*dm)
    par<-matrix(par1[rng,],length(rng),2)##Add this to other functions so dm=1 is ok
    Y <- rnorm(n)
    X <- rnorm(n*dm)
    X1 <- X*sqrt(1-ic^2) + Y*ic
    ds <- qbeta(pnorm(X1),rep(par[,1],rep(n,dm)),rep(par[,2],rep(n,dm)))
    df <- data.frame(M = log(ds/(1-ds)),cpg = rep(1:dm,rep(n,dm)))

    p<-as.vector(tapply(df$M,INDEX=df$cpg,get_cor_pval,Y=Y,np_method=np_method))
    use_this <- sig/Total
    Collect[i,1] <- length(p[p<use_this])/dm
    Collect[i,2] <- use_this
    if(i>=Nmin && sd(Collect$power[1:i])<(sqrt(i)*MOE/1.96)){break}
  }
  Collect[Collect$cutoff>0,]
}



