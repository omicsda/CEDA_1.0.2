
#' Median normalization of sgRNA counts
#'
#' This function adjusts sgRNA counts by the median ratio method.
#' The normalized sgRNA read counts are calculated as the raw read counts 
#' devided by a size factor. The size factor is calcuated as the median of 
#' all size factors caculated from negative control sgRNAs (eg., sgRNAs 
#' corresponding to non-targeting or non-essential genes).
#'
#' @param data A numeric matrix containing raw read counts of sgRNAs
#'   with rows corresponding to sgRNAs and columns correspondings to samples.
#' @param control A numeric matrix containing raw read counts of negative 
#'   control sgRNAs with rows corresponding to sgRNAs and columns 
#'   corresponding to samples. Sample ordering is the same as in data.
#' @return A list with two elements: 1) size factors of all samples; 
#'   2) normalized counts of sgRNAs. 
#' @examples
#' count <- matrix(rnbinom(5000 * 6, mu=500, size=3), ncol = 6)
#' colnames(count) = paste0("sample", 1:6)
#' rownames(count) = paste0("sgRNA", 1:5000)
#' control <- count[1:100,]
#' normalizedcount <- medianNormalization(count, control)
#'
#' @importFrom stats median
#' 
#' @export
medianNormalization <- function(data, control) {
  gm <- exp(rowMeans(log(control+1))) 
  f <- apply(control, 2, function(u) median((u/gm)[gm > 0]))
  norm <- sweep(data,2,f,FUN="/")
  return(list(f=f, count=norm))
}



#' Modeling CRISPR screen data by R package limma
#'
#' The lmFit function in R package limma is employed for group comparisons.
#'  
#' @param data A numeric matrix containing log2 expression levels of sgRNAs
#'   with rows corresponding to sgRNAs and columns corresponding to samples.
#' @param design A design matrix with rows corresponding to samples and
#'   columns corresponding to coefficients to be estimated.
#' @param contrast.matrix A matrix with columns corresponding to contrasts.
#' @return A data frame with rows corresponding to sgRNAs and columns
#'   corresponding to limma results
#' @examples
#' y <- matrix(rnorm(1000*6),1000,6)
#' condition <- gl(2,3,labels=c("Treatment","Baseline"))
#' design <- model.matrix(~ 0 + condition)
#' contrast.matrix <- makeContrasts("conditionTreatment-conditionBaseline",levels=design)
#' limma.fit <- runLimma(y,design,contrast.matrix)
#'
#' @importFrom limma contrasts.fit makeContrasts lmFit eBayes
#'
#' @export
runLimma <- function(data, design, contrast.matrix) {
  lmfit <- limma::lmFit(data,design)
  lmfit.eBayes <- limma::eBayes(contrasts.fit(lmfit, contrast.matrix))
  results <- data.frame(lmfit.eBayes$coef,
                        lmfit.eBayes$stdev.unscaled*lmfit.eBayes$sigma,
                        lmfit.eBayes$p.value)
  names(results) <- c("lfc","se","p")
  return(results)
}



#' Modeling CRISPR data with a permutation test between conditions 
#' by R package limma
#'
#' The lmFit function in R package limma is employed for group comparisons
#' under permutations.
#' 
#' @param data A numeric matrix containing log2 expression level of sgRNAs
#'   with rows corresponding to sgRNAs and columns to samples.
#' @param design A design matrix with rows corresponding to samples and
#'   columns to coefficients to be estimated.
#' @param contrast.matrix A matrix with columns corresponding to contrasts.
#' @param nperm Number of permutations
#' @return A numeric matrix containing log2 fold changes with permutations
#' @examples
#' y <- matrix(rnorm(1000*6),1000,6)
#' condition <- gl(2,3,labels=c("Control","Baseline"))
#' design <- model.matrix(~ 0 + condition)
#' contrast.matrix <- makeContrasts("conditionControl-conditionBaseline",levels=design)
#' fit <- permuteLimma(y,design,contrast.matrix,20)
#'
#' @export
permuteLimma <- function(data, design, contrast.matrix, nperm) {
  n.rna <- dim(data)[1]
  beta.null <- matrix(0,n.rna,nperm)
  ns.grp <- dim(design)[1]/2
  for (j in 1:nperm)
  {
    n.floor <- floor(ns.grp/2)
    n.ceiling <- ceiling(ns.grp/2)
    col.grp1 <- c(sample(1:ns.grp,n.floor),sample((ns.grp+1):(2*ns.grp),n.ceiling))
    col.grp2 <- setdiff(1:(2*ns.grp),col.grp1)
    col.new <- c(col.grp1,col.grp2)
    limma.fit.null <- runLimma(data[,col.new],design,contrast.matrix)
    beta.null[,j] <- limma.fit.null$lfc
  }
  return(beta.null)
}



#' Fitting multi-component normal mixture models by R package mixtools
#'
#' The function normalmixEM in R package mixtools is employed for 
#' fitting multi-component normal mixture models.
#'
#' @param x A numeric vector
#' @param k0 Number of components in the normal mixture model
#' @param mean_constr A constrain on means of components
#' @param sd_constr A constrain on standard deviations of components
#' @param npara Number of parameters
#' @param d0 Number of times for fitting mixture model using different 
#'   starting values
#' @return Normal mixture model fit and BIC value of the log-likelihood
#'
#' @importFrom mixtools normalMixEM
EMFit <- function(x, k0, mean_constr, sd_constr, npara, d0) {
  for (i in 1:d0)
  {
    EM.fit.temp <- mixtools::normalmixEM(x,k=k0,mean.constr=mean_constr,sd.constr=sd_constr)
    if (i==1)
    {
      EM.fit <- EM.fit.temp
    } else
    {
      if (EM.fit.temp$loglik > EM.fit$loglik)
      {
        EM.fit <- EM.fit.temp
      }
    }
  }
  BIC <- -2*EM.fit$loglik + npara*log(length(x))
  return(list(EM.fit=EM.fit,BIC=BIC))
}



#' Performing empirical Bayes modeling on limma results
#' 
#' This function perform an empirical Bayes modeling on log fold ratios
#' and return the posterior log fold ratios.
#' 
#' @param data A numeric matrix containing limma results and log2 gene 
#'   expression levels that has a column nameed 'lfc' and a column 
#'   named 'exp.level.log2'
#' @param theta0 Standard deviation of log2 fold changes under permutations
#' @return A numeric matrix containing limma results, RNA expression levels,
#'   posterior log2 fold ratio, log p-values, and estimates of mixture model
#' @examples
#' \donttest{
#' nmm.fit <- normalMM(data,theta0)
#' }
#'
#' @importFrom stats dnorm pnorm
#'
#' @export
normalMM <- function(data,theta0) {
  eta <- 0.5
  d <- 10
  n.b <- 5
  mu.mat <- matrix(0,n.b,3)
  sigma.mat <- matrix(0,n.b,3)
  lambda.mat <- matrix(0,n.b,3)
  beta.cutoff.mat <- matrix(0,n.b,2)
  xs <- min(data$exp.level.log2)
  xe <- max(data$exp.level.log2)+0.1
  bl <- (xe-xs)/(n.b+1)
  binter <- c(xs,xs+bl,xs+2*bl,xs+3*bl,xs+4*bl,xe)
  for(b in 1:n.b)
  {
    data$exp.level.log2.b[(data$exp.level.log2 >= binter[b]) & (data$exp.level.log2 < binter[b+1])] = b
  }
  for (b in 1:n.b)
  {
    x <- data$lfc[data$exp.level.log2.b==b]
    EMfit.3mm <- EMFit(x,k0=3,mean_constr=c(NA,0,NA),sd_constr=c(NA,theta0,NA),npara=4,d0=d)
    if ( (max(EMfit.3mm$EM.fit$mu) < eta) & (min(EMfit.3mm$EM.fit$mu) < -eta) )
    {
      EMfit.3mm <- EMFit(x,k0=3,mean_constr=c(NA,0,eta),sd_constr=c(NA,theta0,NA),npara=3,d0=d)
      if (min(EMfit.3mm$EM.fit$mu) > -eta)
      {
        EMfit.3mm <- EMFit(x,k0=3,mean_constr=c(-eta,0,eta),sd_constr=c(NA,theta0,NA),npara=2,d0=d)
      }
    }
    if ( (max(EMfit.3mm$EM.fit$mu) > eta) & (min(EMfit.3mm$EM.fit$mu) > -eta) )
    {
      EMfit.3mm <- EMFit(x,k0=3,mean_constr=c(-eta,0,NA),sd_constr=c(NA,theta0,NA),npara=3,d0=d)
      if (max(EMfit.3mm$EM.fit$mu) < eta)
      {
        EMfit.3mm <- EMFit(x,k0=3,mean_constr=c(-eta,0,eta),sd_constr=c(NA,theta0,NA),npara=2,d0=d)
      }
    }
    if ( (max(EMfit.3mm$EM.fit$mu) < eta) & (min(EMfit.3mm$EM.fit$mu) > -eta) )
    {
      EMfit.3mm <- EMFit(x,k0=3,mean_constr=c(-eta,0,eta),sd_constr=c(NA,theta0,NA),npara=2,d0=d)
    }
    EMfit.2mm <- EMFit(x,k0=2,mean_constr=c(NA,0),sd_constr=c(NA,theta0),npara=2,d0=d)
    if ( (EMfit.2mm$EM.fit$mu[1] > 0) &  (EMfit.2mm$EM.fit$mu[1] < eta) )
    {
      EMfit.2mm <- EMFit(x,k0=2,mean_constr=c(eta,0),sd_constr=c(NA,theta0),npara=1,d0=d)
    }
    if ( (EMfit.2mm$EM.fit$mu[1] < 0) &  (EMfit.2mm$EM.fit$mu[1] > -eta) )
    {
      EMfit.2mm <- EMFit(x,k0=2,mean_constr=c(-eta,0),sd_constr=c(NA,theta0),npara=1,d0=d)
    }
    BIC.1mm <- -2*sum(dnorm(x,0,theta0,log=TRUE)) + log(length(x))
    BIC.all <- data.frame(EMfit.3mm$BIC, EMfit.2mm$BIC, BIC.1mm)
    if (which(rank(BIC.all)==1)==1)
    {
      xorder <- sort(EMfit.3mm$EM.fit$x,index=T)$ix
      x.ordered <- x[xorder]
      posterior.xorder <- EMfit.3mm$EM.fit$posterior[xorder,]
      n1 <- which(posterior.xorder[,1]<posterior.xorder[,2])[1]-1
      if (length(n1)==0)
      {
        x.comp1.cutoff <- min(x.ordered)-1
      } else
      {
        x.comp1.cutoff <- x.ordered[n1]
      }
      x.ordered.pos <- x.ordered[x.ordered>0]
      n3 <- which(posterior.xorder[x.ordered>0,3]>posterior.xorder[x.ordered>0,2])[1]
      if (length(n3)==0)
      {
        x.comp3.cutoff <- max(x.ordered)+1
      } else
      {
        x.comp3.cutoff <- x.ordered.pos[n3]
      }
      x.cutoff <- c(x.comp1.cutoff,x.comp3.cutoff)
      null.posterior <- EMfit.3mm$EM.fit$posterior[,2]
      beta.cutoff.mat[b,] <- x.cutoff
      mu.mat[b,] <- EMfit.3mm$EM.fit$mu
      sigma.mat[b,] <- EMfit.3mm$EM.fit$sigma
      lambda.mat[b,] <- EMfit.3mm$EM.fit$lambda
    }
    if (which(rank(BIC.all)==1)==2)
    {
      xorder <- sort(EMfit.2mm$EM.fit$x,index=T)$ix
      x.ordered <- x[xorder]
      posterior.xorder <- EMfit.2mm$EM.fit$posterior[xorder,]
      if(EMfit.2mm$EM.fit$mu[1]<0)
      {
        n1 <- which(posterior.xorder[,1]<posterior.xorder[,2])[1]-1
        if (length(n1)==0)
        {
          x.comp1.cutoff <- min(x.ordered)-1
        } else
        {
          x.comp1.cutoff <- x.ordered[n1]
        }
      } else
      {
        n1 <- which(posterior.xorder[,1]>posterior.xorder[,2])[1]
        if (length(n1)==0)
        {
          x.comp1.cutoff <- max(x.ordered)+1
        } else
        {
          x.comp1.cutoff <- x.ordered[n1]
        }
      }
      x.cutoff <- c(x.comp1.cutoff)
      null.posterior <- EMfit.3mm$EM.fit$posterior[,2]

      beta.cutoff.mat[b,1] <- x.cutoff
      mu.mat[b,1:2] <- EMfit.2mm$EM.fit$mu
      sigma.mat[b,1:2] <- EMfit.2mm$EM.fit$sigma
      lambda.mat[b,1:2] <- EMfit.2mm$EM.fit$lambda
    }
    if (b==1) {null.posterior.b1 <- null.posterior}
    if (b==2) {null.posterior.b2 <- null.posterior}
    if (b==3) {null.posterior.b3 <- null.posterior}
    if (b==4) {null.posterior.b4 <- null.posterior}
    if (b==5) {null.posterior.b5 <- null.posterior}
  }
  data.b1 <- data[data$exp.level.log2.b==1,]
  data.b1$null.posterior <- null.posterior.b1
  data.b2 <- data[data$exp.level.log2.b==2,]
  data.b2$null.posterior <- null.posterior.b2
  data.b3 <- data[data$exp.level.log2.b==3,]
  data.b3$null.posterior <- null.posterior.b3
  data.b4 <- data[data$exp.level.log2.b==4,]
  data.b4$null.posterior <- null.posterior.b4
  data.b5 <- data[data$exp.level.log2.b==5,]
  data.b5$null.posterior <- null.posterior.b5
  b <- 1
  if (beta.cutoff.mat[b,1]==0 & beta.cutoff.mat[b,2]==0)
  {
    data.b1$lfc_posterior <- (data.b1$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b1$se)^2)/((data.b1$se)^2 + (sigma.mat[b,2])^2)
  }
  if (beta.cutoff.mat[b,1]<0 & beta.cutoff.mat[b,2]>0)
  {
    data.b1.comp2 <- data.b1[data.b1$lfc <= beta.cutoff.mat[b,1],]
    data.b1.comp0 <- data.b1[data.b1$lfc >  beta.cutoff.mat[b,1] & data.b1$lfc < beta.cutoff.mat[b,2],]
    data.b1.comp1 <- data.b1[data.b1$lfc >= beta.cutoff.mat[b,2],]
    data.b1.comp2$lfc_posterior <- (data.b1.comp2$lfc * (sigma.mat[b,1])^2 + mu.mat[b,1] * (data.b1.comp2$se)^2)/((data.b1.comp2$se)^2 + (sigma.mat[b,1])^2)
    data.b1.comp0$lfc_posterior <- (data.b1.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b1.comp0$se)^2)/((data.b1.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b1.comp1$lfc_posterior <- (data.b1.comp1$lfc * (sigma.mat[b,3])^2 + mu.mat[b,3] * (data.b1.comp1$se)^2)/((data.b1.comp1$se)^2 + (sigma.mat[b,3])^2)
    data.b1 <- rbind(data.b1.comp0,data.b1.comp1,data.b1.comp2)
  }
  if (beta.cutoff.mat[b,1]<0 & beta.cutoff.mat[b,2]==0)
  {
    data.b1.comp2 <- data.b1[data.b1$lfc <= beta.cutoff.mat[b,1],]
    data.b1.comp0 <- data.b1[data.b1$lfc >  beta.cutoff.mat[b,1],]
    data.b1.comp2$lfc_posterior <- (data.b1.comp2$lfc * (sigma.mat[b,1])^2 + mu.mat[b,1] * (data.b1.comp2$se)^2)/((data.b1.comp2$se)^2 + (sigma.mat[b,1])^2)
    data.b1.comp0$lfc_posterior <- (data.b1.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b1.comp0$se)^2)/((data.b1.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b1 <- rbind(data.b1.comp0,data.b1.comp2)
  }
  if (beta.cutoff.mat[b,1]>0 & beta.cutoff.mat[b,2]==0)
  {
    data.b1.comp0 <- data.b1[data.b1$lfc < beta.cutoff.mat[b,1],]
    data.b1.comp1 <- data.b1[data.b1$lfc >= beta.cutoff.mat[b,1],]
    data.b1.comp0$lfc_posterior <- (data.b1.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b1.comp0$se)^2)/((data.b1.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b1.comp1$lfc_posterior <- (data.b1.comp1$lfc * (sigma.mat[b,3])^2 + mu.mat[b,3] * (data.b1.comp1$se)^2)/((data.b1.comp1$se)^2 + (sigma.mat[b,3])^2)
    data.b1 <- rbind(data.b1.comp0,data.b1.comp1)
  }
  b <- 2
  if (beta.cutoff.mat[b,1]==0 & beta.cutoff.mat[b,2]==0)
  {
    data.b2$lfc_posterior <- (data.b2$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b2$se)^2)/((data.b2$se)^2 + (sigma.mat[b,2])^2)
  }
  if (beta.cutoff.mat[b,1]<0 & beta.cutoff.mat[b,2]>0)
  {
    data.b2.comp2 <- data.b2[data.b2$lfc <= beta.cutoff.mat[b,1],]
    data.b2.comp0 <- data.b2[data.b2$lfc >  beta.cutoff.mat[b,1] & data.b2$lfc < beta.cutoff.mat[b,2],]
    data.b2.comp1 <- data.b2[data.b2$lfc >= beta.cutoff.mat[b,2],]
    data.b2.comp2$lfc_posterior <- (data.b2.comp2$lfc * (sigma.mat[b,1])^2 + mu.mat[b,1] * (data.b2.comp2$se)^2)/((data.b2.comp2$se)^2 + (sigma.mat[b,1])^2)
    data.b2.comp0$lfc_posterior <- (data.b2.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b2.comp0$se)^2)/((data.b2.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b2.comp1$lfc_posterior <- (data.b2.comp1$lfc * (sigma.mat[b,3])^2 + mu.mat[b,3] * (data.b2.comp1$se)^2)/((data.b2.comp1$se)^2 + (sigma.mat[b,3])^2)
    data.b2 <- rbind(data.b2.comp0,data.b2.comp1,data.b2.comp2)
  }
  if (beta.cutoff.mat[b,1]<0 & beta.cutoff.mat[b,2]==0)
  {
    data.b2.comp2 <- data.b2[data.b2$lfc <= beta.cutoff.mat[b,1],]
    data.b2.comp0 <- data.b2[data.b2$lfc >  beta.cutoff.mat[b,1],]
    data.b2.comp2$lfc_posterior <- (data.b2.comp2$lfc * (sigma.mat[b,1])^2 + mu.mat[b,1] * (data.b2.comp2$se)^2)/((data.b2.comp2$se)^2 + (sigma.mat[b,1])^2)
    data.b2.comp0$lfc_posterior <- (data.b2.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b2.comp0$se)^2)/((data.b2.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b2 <- rbind(data.b2.comp0,data.b2.comp2)
  }
  if (beta.cutoff.mat[b,1]>0 & beta.cutoff.mat[b,2]==0)
  {
    data.b2.comp0 <- data.b2[data.b2$lfc < beta.cutoff.mat[b,1],]
    data.b2.comp1 <- data.b2[data.b2$lfc >= beta.cutoff.mat[b,1],]
    data.b2.comp0$lfc_posterior <- (data.b2.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b2.comp0$se)^2)/((data.b2.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b2.comp1$lfc_posterior <- (data.b2.comp1$lfc * (sigma.mat[b,3])^2 + mu.mat[b,3] * (data.b2.comp1$se)^2)/((data.b2.comp1$se)^2 + (sigma.mat[b,3])^2)
    data.b2 <- rbind(data.b2.comp0,data.b2.comp1)
  }
  b <- 3
  if (beta.cutoff.mat[b,1]==0 & beta.cutoff.mat[b,2]==0)
  {
    data.b3$lfc_posterior <- (data.b3$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b3$se)^2)/((data.b3$se)^2 + (sigma.mat[b,2])^2)
  }
  if (beta.cutoff.mat[b,1]<0 & beta.cutoff.mat[b,2]>0)
  {
    data.b3.comp2 <- data.b3[data.b3$lfc <= beta.cutoff.mat[b,1],]
    data.b3.comp0 <- data.b3[data.b3$lfc >  beta.cutoff.mat[b,1] & data.b3$lfc < beta.cutoff.mat[b,2],]
    data.b3.comp1 <- data.b3[data.b3$lfc >= beta.cutoff.mat[b,2],]
    data.b3.comp2$lfc_posterior <- (data.b3.comp2$lfc * (sigma.mat[b,1])^2 + mu.mat[b,1] * (data.b3.comp2$se)^2)/((data.b3.comp2$se)^2 + (sigma.mat[b,1])^2)
    data.b3.comp0$lfc_posterior <- (data.b3.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b3.comp0$se)^2)/((data.b3.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b3.comp1$lfc_posterior <- (data.b3.comp1$lfc * (sigma.mat[b,3])^2 + mu.mat[b,3] * (data.b3.comp1$se)^2)/((data.b3.comp1$se)^2 + (sigma.mat[b,3])^2)
    data.b3 <- rbind(data.b3.comp0,data.b3.comp1,data.b3.comp2)
  }
  if (beta.cutoff.mat[b,1]<0 & beta.cutoff.mat[b,2]==0)
  {
    data.b3.comp2 <- data.b3[data.b3$lfc <= beta.cutoff.mat[b,1],]
    data.b3.comp0 <- data.b3[data.b3$lfc >  beta.cutoff.mat[b,1],]
    data.b3.comp2$lfc_posterior <- (data.b3.comp2$lfc * (sigma.mat[b,1])^2 + mu.mat[b,1] * (data.b3.comp2$se)^2)/((data.b3.comp2$se)^2 + (sigma.mat[b,1])^2)
    data.b3.comp0$lfc_posterior <- (data.b3.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b3.comp0$se)^2)/((data.b3.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b3 <- rbind(data.b3.comp0,data.b3.comp2)
  }
  if (beta.cutoff.mat[b,1]>0 & beta.cutoff.mat[b,2]==0)
  {
    data.b3.comp0 <- data.b3[data.b3$lfc < beta.cutoff.mat[b,1],]
    data.b3.comp1 <- data.b3[data.b3$lfc >= beta.cutoff.mat[b,1],]
    data.b3.comp0$lfc_posterior <- (data.b3.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b3.comp0$se)^2)/((data.b3.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b3.comp1$lfc_posterior <- (data.b3.comp1$lfc * (sigma.mat[b,3])^2 + mu.mat[b,3] * (data.b3.comp1$se)^2)/((data.b3.comp1$se)^2 + (sigma.mat[b,3])^2)
    data.b3 <- rbind(data.b3.comp0,data.b3.comp1)
  }
  b <- 4
  if (beta.cutoff.mat[b,1]==0 & beta.cutoff.mat[b,2]==0)
  {
    data.b4$lfc_posterior <- (data.b4$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b4$se)^2)/((data.b4$se)^2 + (sigma.mat[b,2])^2)
  }
  if (beta.cutoff.mat[b,1]<0 & beta.cutoff.mat[b,2]>0)
  {
    data.b4.comp2 <- data.b4[data.b4$lfc <= beta.cutoff.mat[b,1],]
    data.b4.comp0 <- data.b4[data.b4$lfc >  beta.cutoff.mat[b,1] & data.b4$lfc < beta.cutoff.mat[b,2],]
    data.b4.comp1 <- data.b4[data.b4$lfc >= beta.cutoff.mat[b,2],]
    data.b4.comp2$lfc_posterior <- (data.b4.comp2$lfc * (sigma.mat[b,1])^2 + mu.mat[b,1] * (data.b4.comp2$se)^2)/((data.b4.comp2$se)^2 + (sigma.mat[b,1])^2)
    data.b4.comp0$lfc_posterior <- (data.b4.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b4.comp0$se)^2)/((data.b4.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b4.comp1$lfc_posterior <- (data.b4.comp1$lfc * (sigma.mat[b,3])^2 + mu.mat[b,3] * (data.b4.comp1$se)^2)/((data.b4.comp1$se)^2 + (sigma.mat[b,3])^2)
    data.b4 <- rbind(data.b4.comp0,data.b4.comp1,data.b4.comp2)
  }
  if (beta.cutoff.mat[b,1]<0 & beta.cutoff.mat[b,2]==0)
  {
    data.b4.comp2 <- data.b4[data.b4$lfc <= beta.cutoff.mat[b,1],]
    data.b4.comp0 <- data.b4[data.b4$lfc >  beta.cutoff.mat[b,1],]
    data.b4.comp2$lfc_posterior <- (data.b4.comp2$lfc * (sigma.mat[b,1])^2 + mu.mat[b,1] * (data.b4.comp2$se)^2)/((data.b4.comp2$se)^2 + (sigma.mat[b,1])^2)
    data.b4.comp0$lfc_posterior <- (data.b4.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b4.comp0$se)^2)/((data.b4.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b4 <- rbind(data.b4.comp0,data.b4.comp2)
  }
  if (beta.cutoff.mat[b,1]>0 & beta.cutoff.mat[b,2]==0)
  {
    data.b4.comp0 <- data.b4[data.b4$lfc < beta.cutoff.mat[b,1],]
    data.b4.comp1 <- data.b4[data.b4$lfc >= beta.cutoff.mat[b,1],]
    data.b4.comp0$lfc_posterior <- (data.b4.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b4.comp0$se)^2)/((data.b4.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b4.comp1$lfc_posterior <- (data.b4.comp1$lfc * (sigma.mat[b,3])^2 + mu.mat[b,3] * (data.b4.comp1$se)^2)/((data.b4.comp1$se)^2 + (sigma.mat[b,3])^2)
    data.b4 <- rbind(data.b4.comp0,data.b4.comp1)
  }
  b <- 5
  if (beta.cutoff.mat[b,1]==0 & beta.cutoff.mat[b,2]==0)
  {
    data.b5$lfc_posterior <- (data.b5$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b5$se)^2)/((data.b5$se)^2 + (sigma.mat[b,2])^2)
  }
  if (beta.cutoff.mat[b,1]<0 & beta.cutoff.mat[b,2]>0)
  {
    data.b5.comp2 <- data.b5[data.b5$lfc <= beta.cutoff.mat[b,1],]
    data.b5.comp0 <- data.b5[data.b5$lfc >  beta.cutoff.mat[b,1] & data.b5$lfc < beta.cutoff.mat[b,2],]
    data.b5.comp1 <- data.b5[data.b5$lfc >= beta.cutoff.mat[b,2],]
    data.b5.comp2$lfc_posterior <- (data.b5.comp2$lfc * (sigma.mat[b,1])^2 + mu.mat[b,1] * (data.b5.comp2$se)^2)/((data.b5.comp2$se)^2 + (sigma.mat[b,1])^2)
    data.b5.comp0$lfc_posterior <- (data.b5.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b5.comp0$se)^2)/((data.b5.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b5.comp1$lfc_posterior <- (data.b5.comp1$lfc * (sigma.mat[b,3])^2 + mu.mat[b,3] * (data.b5.comp1$se)^2)/((data.b5.comp1$se)^2 + (sigma.mat[b,3])^2)
    data.b5 <- rbind(data.b5.comp0,data.b5.comp1,data.b5.comp2)
  }
  if (beta.cutoff.mat[b,1]<0 & beta.cutoff.mat[b,2]==0)
  {
    data.b5.comp2 <- data.b5[data.b5$lfc <= beta.cutoff.mat[b,1],]
    data.b5.comp0 <- data.b5[data.b5$lfc >  beta.cutoff.mat[b,1],]
    data.b5.comp2$lfc_posterior <- (data.b5.comp2$lfc * (sigma.mat[b,1])^2 + mu.mat[b,1] * (data.b5.comp2$se)^2)/((data.b5.comp2$se)^2 + (sigma.mat[b,1])^2)
    data.b5.comp0$lfc_posterior <- (data.b5.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b5.comp0$se)^2)/((data.b5.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b5 <- rbind(data.b5.comp0,data.b5.comp2)
  }
  if (beta.cutoff.mat[b,1]>0 & beta.cutoff.mat[b,2]==0)
  {
    data.b5.comp0 <- data.b5[data.b5$lfc < beta.cutoff.mat[b,1],]
    data.b5.comp1 <- data.b5[data.b5$lfc >= beta.cutoff.mat[b,1],]
    data.b5.comp0$lfc_posterior <- (data.b5.comp0$lfc * (sigma.mat[b,2])^2 + mu.mat[b,2] * (data.b5.comp0$se)^2)/((data.b5.comp0$se)^2 + (sigma.mat[b,2])^2)
    data.b5.comp1$lfc_posterior <- (data.b5.comp1$lfc * (sigma.mat[b,3])^2 + mu.mat[b,3] * (data.b5.comp1$se)^2)/((data.b5.comp1$se)^2 + (sigma.mat[b,3])^2)
    data.b5 <- rbind(data.b5.comp0,data.b5.comp1)
  }
  data <- rbind(data.b1, data.b2, data.b3, data.b4, data.b5)
  data$log_p <- log(2)+pnorm(abs(data$lfc_posterior),mean=0,sd=theta0,lower.tail=FALSE,log.p=TRUE)
  data$log_p_noshrink <- log(2)+pnorm(abs(data$lfc),mean=0,sd=theta0,lower.tail=FALSE,log.p=TRUE)
  return(list(data=data,beta.cutoff.mat=beta.cutoff.mat,mu.mat=mu.mat,sigma.mat=sigma.mat,
              lambda.mat=lambda.mat))
}



#' Scatter plot of log2 fold ratios against gene expression levels
#'
#' This function generates a scatter plot of log2 fold ratios of sgRNAs
#' against the corresponding gene expression levels.
#'
#' @param data A numeric matrix from the output of normalMM function
#' @param fdr A level of false discovery rate
#' @param ... Other graphical parameters
#'
#' @return No return value
#' @importFrom ggplot2 aes_string aes geom_point geom_vline theme theme_bw element_blank xlab ylab
#'
#' @export
scatterPlot <- function(data, fdr, ...) {
  b <- 6
  xs <- min(data$exp.level.log2)
  xe <- max(data$exp.level.log2)+0.1
  bl <- (xe-xs)/b
  binter <- c(xs,xs+bl,xs+2*bl,xs+3*bl,xs+4*bl,xe)
  data.fdr <- data[data$null.posterior<=fdr,]
  exp.level.log2 <- lfc <- NULL
  ggplot2::ggplot(data, aes(x = exp.level.log2, y = lfc)) +
    geom_point(size=1,alpha=0.2) +
    xlab("log2(gene expression)") + ylab("log2(FC)") +
    geom_point(data=data.fdr,aes(x=exp.level.log2,y=lfc),size=1,alpha=0.2,color="red3") +
    geom_vline(aes(xintercept=binter[2]),color="cyan") +
    geom_vline(aes(xintercept=binter[3]),color="cyan") +
    geom_vline(aes(xintercept=binter[4]),color="cyan") +
    geom_vline(aes(xintercept=binter[5]),color="cyan") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}



#' Calculating a significance score of a gene based on 
#' the corresponding sgRNAs' p-values of the gene.
#'
#' Code was adapted from R package gscreend.
#'
#' @param pvec A numeric vector of p-values.
#' @return A min value of the kth smallest value based on the beta 
#'   distribution B(k, n-k+1), where the n is the number of probabiliteis 
#'   in the vector. This min value is the significance score of the gene.
alphaBeta <- function(pvec) {
  pvec <- sort(pvec)
  n <- length(pvec)
  min(stats::pbeta(pvec, seq_len(n), n - seq_len(n) + 1))
}



#' Generating the null distribution of the significance score 
#' of a gene.
#'
#' Code was adapted from R package gscreend.
#'
#' @param n An integer representing sgRNA number of a gene.
#' @param p A numeric vector which contains the percentiles of the 
#'   p-values that meet the cut-off (alpha).
#' @param nperm Number of permutation runs.
#' @return A numric vector which contains all the significance scores 
#'   (rho) of genes generated by a permutation test where the sgRNAs are 
#'   randomly assigned to genes.
makeRhoNull <- function(n, p, nperm) {
  rhonull <- lapply(seq_len(nperm), function(x) {
    p_test <- sort.int(sample(p, n, replace = FALSE))
    p_test <- sort(p_test)
    min(stats::pbeta(p_test, seq_len(n), n - seq_len(n) + 1))
  })
  unlist(rhonull)
}



#' Calculating gene level p-values using modified robust rank aggregation
#' (alpha-RRA method) on sgRNAs' p-values
#'
#' Code was adapted from R package gscreend. The alpha-RRA method is 
#' adapted from MAGeCK.
#'
#' @param pvec A numeric vector containing p-values of sgRNAs.
#' @param genes A character string containing gene names corresponding 
#'   to sgRNAs.
#' @param alpha A numeric number denoting the alpha cutoff (i.e. 0.05).
#' @return A list with four elements: 1) a list of genes with their p-values; 
#'   2) a numeric matrix of rho null, each column corresponding to a different 
#'   number of sgRNAs per gene; 3)a numeric vector of rho; 4) a numeric vector 
#'   of number of sgRNAs per gene.
#' @export
calculateGenePval <- function(pvec, genes, alpha) {
  cut.pvec <- pvec <= alpha
  score_vals <- rank(pvec)/length(pvec)
  score_vals[!cut.pvec] <- 1
  rho <- unsplit(vapply(split(score_vals, genes),
                        FUN = alphaBeta,
                        FUN.VALUE = numeric(1)),
                 genes)
  guides_per_gene <- sort(unique(table(genes)))
  permutations = 10 * length(unique(genes))
  rho_nullh <- vapply(guides_per_gene,
                      FUN = makeRhoNull,
                      p = score_vals,
                      nperm = permutations,
                      FUN.VALUE = numeric(permutations))
  pvalue_gene <- lapply(split(rho, genes), function(x) {
    n_sgrnas = length(x)
    mean(rho_nullh[, guides_per_gene == n_sgrnas] <= x[[1]])
  })
  result <- list(pvalue=pvalue_gene, rho_null=rho_nullh, rho=rho, 
                 guides_per_gene=guides_per_gene)
  return(result)
}



#' Calculating gene-level log fold ratios
#' 
#' Log fold ratios of all sgRNAs of a gene are averaged to obtain the 
#' gene level log fold ratio.
#' 
#' @param lfcs A numeric vector containing log fold change of sgRNAs.
#' @param genes A character string containing gene names corresponding to sgRNAs.
#' @return A numeric vector containing log fold ratio of genes.
#'
#' @export
calculateGeneLFC <- function(lfcs, genes) {
  vapply(split(lfcs, genes), FUN = mean, FUN.VALUE = numeric(1))
}





