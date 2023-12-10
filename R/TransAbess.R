#' @importFrom abess abess
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @importFrom Matrix Matrix
#' @importFrom Matrix bdiag
#' @importFrom stats coef
#' @importFrom stats rnorm
#' @import Rcpp
#' @import RcppArmadillo
#' @useDynLib SA23204160
NULL

#' @import boot
#' @import bootstrap
#' @import coda
#' @import DAAG
#' @import dplyr
#' @import foreach
#' @import twosamples
NULL

OracleTransAbess <- function(target, source = NULL, pred.target = NULL, transfer.source.id = NULL) {
  transfer.source.id.ori <- transfer.source.id
  k <- NULL
  p <- ncol(target$x)
  
  if (!is.null(source) && (is.character(transfer.source.id) && transfer.source.id == "all")) { # transfer all source data
    transfer.source.id <- 1:length(source)
  } else if (0 %in% transfer.source.id || is.null(source)) { # don't transfer any source
    transfer.source.id <- 0
  } else { # transfer oracle source data
    transfer.source.id <- transfer.source.id
  }
  
  K_B <- length(transfer.source.id)
  
    all.x <- as.matrix(foreach(k = unique(c(0, transfer.source.id)), .combine = "bdiag") %do% {
      if (k != 0) {
        source[[k]]$x
      } else {
        target$x
      }
    })
    
    all.x <- Matrix(all.x, sparse = TRUE)

    C <- GenerateMatrixC(K_B,p)
    
    all.x.new <- all.x %*% C
    
    all.y <- foreach(k = unique(c(0, transfer.source.id)), .combine = "c") %do% {
      if (k != 0) {
        source[[k]]$y
      } else {
        target$y
      }
    }

    # prevent columns valued 0
    # index <- which(apply(all.x.new,2,function(x) all(x==0)))
    v <- 1:ncol(target$x)
    # if (!is.null(index)){
    #   X <- all.x.new[,-index]
    #   v <- v[-index]
    #   transfer.fit <- abess(x = X, y = all.y, tune.type = 'cv')
    #   wa <- coef(transfer.fit, support.size = transfer.fit[["best.size"]])
    #   beta.hat <- as.vector(wa[1:(length(v)+1)])
    # }
      transfer.fit <- abess(x = all.x.new, y = all.y, tune.type = 'cv')
      wa <- coef(transfer.fit, support.size = transfer.fit[["best.size"]])
      beta.hat <- as.vector(wa[1:(length(v)+1)])
  
  if (!all(is.null(colnames(target$x)))) {
    names(beta.hat) <- c("intercept", colnames(target$x)[v])
  } else {
    names(beta.hat) <- c("intercept", paste0("V", v))
  }
  
  if (0 %in% transfer.source.id) {
    transfer.source.id <- NULL
  }
  
  if (!is.null(pred.target)) {
    loss <- sum((cbind(rep(1,length(pred.target$y)), pred.target$x[,v]) %*% beta.hat - pred.target$y)^2)/2/length(pred.target$y)
  } else {loss <- NA}
  
  obj <- list(beta = beta.hat, transfer.source.id = transfer.source.id, loss = loss)
  class(obj) <- "OracleTransAbess"
  return(obj)
}



#' Fit a transfer learning linear model based on \code{abess} algorithm.
#'
#' @param target target data. Should be a list with elements x and y, where x indicates a predictor matrix with each row/column as a(n) observation/variable, and y indicates the response vector.
#' @param source source data. Should be a list with some sublists, where each of the sublist is a source data set, having elements x and y with the same meaning as in target data.
#' @param pred.target a new target data. Should be a list with elements x and y which we want to predict.
#' @param transfer.source.id transferable source indices. Can be either a subset of \code{{1, ..., length(source)}}, "all" or "auto". Default = \code{"auto"}.
#' \itemize{
#' \item a subset of \code{{1, ..., length(source)}}: only transfer sources with the specific indices.
#' \item "all": transfer all sources.
#' \item "auto": run transferable source detection algorithm to automatically detect which sources to transfer.
#' }
#'
#' @return a list named \code{TransAbess}, containing the parameter estimation \code{beta.hat} and the transferable source indices estimation \code{transfer.source.id} and the prediction error.
#' @export
#' 
#' @examples
#' \dontrun{
#' data <- models(type = 'all', sig.strength = 1, p = 500, s = 10, 
#' h = 5, h.not.in = 100, K = 10, Ka = 5)
#' transfer.source.id <- 1:5
#' fit.OracleTransABESS <- TransAbess(target = data$target, source = data$source, 
#' transfer.source.id = transfer.source.id)
#' }
TransAbess <- function(target, source = NULL, pred.target = NULL, transfer.source.id = 'auto') {
  K <- length(source)
  
  if (!is.null(source) && (is.character(transfer.source.id) && transfer.source.id == "all")) { # transfer all source data
    obj <- OracleTransAbess(target = target, source = source, pred.target = pred.target, transfer.source.id = 1:K)
  } else if (0 %in% transfer.source.id || is.null(source)) { # don't transfer any source
    obj <- OracleTransAbess(target = target, source = source, pred.target = pred.target, transfer.source.id = 0)
  } else if (!is.null(source) && (is.character(transfer.source.id) && transfer.source.id == "auto")) {
    
    I <- sample(1:length(target$y), length(target$y) %/% 2, replace = FALSE)
    target1 = target2 = list()
    target1$x <- target$x[I,]
    target1$y <- target$y[I]
    target2$x <- target$x[-I,]
    target2$y <- target$y[-I]
    
    R.hat <- foreach(k = 1:K, .combine = "c") %do% {
      Delta.hat <- t(source[[k]]$x) %*% source[[k]]$y / length(source[[k]]$y) - t(target1$x) %*% target1$y / length(target1$y)
      sum(Delta.hat^2)
    }
    
    all.beta.hat <- foreach(k = 0:K) %do% {
      if (k == 0){
        G.hat <- 0
      } else {
        G.hat <- which(rank(R.hat) <= rank(R.hat)[k])
      }
      fit <- OracleTransAbess(target = target1, source = source, pred.target = target2, transfer.source.id = G.hat)
      list(beta.hat = fit$beta, transfer.source.id = G.hat, loss = fit$loss)
    }
    source.loss <- foreach(k = 0:K, .combine = 'c') %do% {
      all.beta.hat[[k+1]]$loss
    }
    k <- which.min(source.loss)
    transfer.source.id <- all.beta.hat[[k]]$transfer.source.id
    
    for (i in 0:length(source)) {
      cat(paste0("Source ", paste(all.beta.hat[[i+1]]$transfer.source.id, collapse = ', '), ": ", format(round(source.loss[i+1], 6), nsmall = 6) , "\n"))
    }
    cat("\n")
    cat(paste("Source data set(s)", paste(transfer.source.id, collapse = ", "), "are transferable!\n"))
    
    obj <- OracleTransAbess(target = target, source = source, pred.target = pred.target, transfer.source.id = transfer.source.id)
    
  } else { # transfer oracle source data
    obj <- OracleTransAbess(target = target, source = source, pred.target = pred.target, transfer.source.id = transfer.source.id)
  }

  class(obj) <- "TransAbess"
  return(obj)
}

#' Generate different target data and source data.
#'
#' @param type different data to be generated. Can be either "all" or "source" or "target". 
#' \itemize{
#' \item "all": generate a list with a target data set of size \code{n.target} and \code{K} source data set of size \code{n.source}.
#' \item "source": generate a list with K source data set of size \code{n.source}.
#' \item "target": generate a list with a target data set of size \code{n.target}.
#' }
#' @param beta regression coefficient. Default = \code{c(3, 1.5, 1, 0.5, 2, 2.5, 4, 3.5)}.
#' @param p dimension of \code{beta}. Default is 8.
#' @param s sparse level of \code{beta}.
#' @param sig.strength signal strength in \code{beta}, i.e. the size of nonzero components.
#' @param cov.strength strength of elements in the covariance matrix.
#' @param eps.strength variance of the residuals.
#' @param h transferring level, i.e. the difference of the regression parameter between target data and transferring source data.
#' @param h.not.in to set the regression parameter of source data which can not be transferred.
#' @param K size of source data.
#' @param Ka size of transferring source data.
#' @param n.target sample size of target data.
#' @param n.source sample size of source data.
#'
#' @return a list containing target and source data.
#' @export
#'
#' @examples
#' \dontrun{
#' data <- models(type = 'all', K = 10, Ka = 5)
#' data <- models(type = 'all', sig.strength = 1, p = 500, s = 10, 
#' h = 5, h.not.in = 100, K = 10, Ka = 5)
#' data <- models(type = 'all', sig.strength = 1, beta = c(3, 1.5, 0, 0, 2, 0, 0, 0), K = 10, Ka = 5)
#' }
models <- function(type = c("all", "source", "target"), beta = c(3, 1.5, 1, 0.5, 2, 2.5, 4, 3.5), p = 8, s = NULL, sig.strength = 1, cov.strength = 0.5, eps.strength = 0.3, 
                   h = 2, h.not.in = 6, K = 10, Ka = NULL, n.target = 200, n.source = rep(100, K)) {
  
  target <- list(x = NULL, y = NULL, beta = NULL)
  source <- NULL
  Sigma <- outer(1:p, 1:p, function(x,y){
    cov.strength^(abs(x-y))
  })
  if (!is.null(s)) {
    beta <- c(rep(sig.strength, s), rep(0, p-s))
  }
  if (type == 'all' || type == 'target'){
    R <- chol(Sigma)
    target$x <- matrix(rnorm(n.target*p), nrow = n.target) %*% R
    target$y <- as.numeric(target$x %*% beta + rnorm(n.target))
    target$beta <- beta
  }
  if (type == 'all' || type =='source'){
    eps <- rnorm(p, sd = eps.strength)
    Sigma <- Sigma + eps %*% t(eps)
    R <- chol(Sigma)
    source <- sapply(1:K, function(k){
      if (k <= Ka){
        wk <- beta
        for (i in sample(1:p, size = h, replace = FALSE)) {
          if(i %in% sample(which(beta != 0))) wk[i] <- wk[i] + 0.5*wk[i]*sample(c(-1,1),1)
          else wk[i] <- 0.5*wk[i]*sample(c(-1,1),1)
        }
      } else {
        wk <- beta
        for (i in sample(1:p, size = h.not.in, replace = FALSE)) {
          if(i %in% sample(which(beta != 0))) wk[i] <- wk[i] + 0.5*wk[i]*sample(c(-1,1),1)
          else wk[i] <- 0.5*wk[i]*sample(c(-1,1),1)
        }
      }
      x <- matrix(rnorm(n.source[k]*p), nrow = n.source[k]) %*% R
      y <- as.numeric(x %*% wk + rnorm(n.source[k]))
      list(x = x, y = y, wk = wk)
    }, simplify = FALSE)
  }
  
  if (type == "all") {
    return(list(target = target, source = source))
  } else if (type == "target") {
    return(list(target = target))
  } else {
    return(list(source = source))
  }
  
}

