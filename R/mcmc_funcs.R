#' A running mean function.
#'
#' This function calculates a running mean for columns of a data frame.
#' @param df A data frame of samples. Rows correspond with iterations, and columns with parameters.
#' @param cols column indices for the columns you want to look at.
#' @param burn Number of rows you want to discard as burn in. Defaults to 0.
#' @keywords mean means running
#' @export
#' @examples
#' fake_samples <- data.frame(rnorm(100))
#' plot(runningMeans(fake_samples, 1), type ="l")
runningMeans <- function(df, cols, burn=0){

  start_index <- burn+1
  apply(df[start_index:nrow(df),cols, drop=FALSE],
        2,
        function(vec) cumsum(vec)/(1:length(vec)))
}


#' A acceptance rate function.
#'
#' This function calculates the acceptance rate for samples stored in a data frame.
#' @param df A data frame of samples. Rows correspond with iterations, and columns with parameters.
#' @param burn Number of rows you want to discard as burn in. Defaults to 0.
#' @keywords accept acceptance rate
#' @export
#' @examples
#' fake_samples <- data.frame(rnorm(100))
#' acceptRate(fake_samples)
acceptRate <- function(df, burn=0){
  start <- burn + 1
  mean(abs(diff(df[start:nrow(df),1])) > .000000001)
}


#' A correlogram function.
#'
#' This function plots overlapping correlograms for samples stored in a data frame.
#' @param df A data frame of samples. Rows correspond with iterations, and columns with parameters.
#' @param cols column indices of the parameters you're interested in plotting.
#' @param burn Number of rows you want to discard as burn in. Defaults to 0.
#' @param ... Extra arguments to be passed into stats::acf() (e.g. lag.max )
#' @keywords correlogram acf autocorrelation
#' @export
#' @examples
#' fake_samples <- data.frame(rnorm(100))
#' correlogram(fake_samples, 1, lag.max = 100)
correlogram <- function(df, cols, burn = 0, ...){

  # calculate acfs
  start <- burn + 1
  df <- df[start:nrow(df),cols, drop = FALSE]
  the_acfs <- acf(df, plot=F, ...)$acf

  # plot acfs
  plot(the_acfs[,1,1], col =1, type = "l", ylim = c(-1,1))
  abline(h=0)
  if(dim(the_acfs)[3] > 1){
    for(i in 2:length(cols)){
      lines(the_acfs[,i,i], col = i)
    }
  }
}


#' A log-sum-exp function.
#'
#' This function calculates the function log(sum(exp(vec))) in a numerically stable way.
#' @param x a numeric vector.
#' @keywords log-sum-exp
#' @export
#' @examples
#' log.sum.exp(rep(1e5, 3))
#' log(sum(exp(rep(1e5,3))))
log.sum.exp<- function(x) {
  if ( max(abs(x)) > max(x) ) # too negative
    offset <- min(x)
  else  # too big
    offset <- max(x)
  log(sum(exp(x - offset))) + offset
}


#' A Gelman-Rubin convergence diagnostic function.
#'
#' This function calculates the Gelman-Rubin convergence diagnostic (Rhat) on some samples.
#' It assumes each chain is in a separate csv file, and each file contains the samples for
#' all parameters. For each chain, this function splits each column into two, which means
#' "m" comes out to be twice the number of csv files. Finally, this function assumes each
#' csv has a one-line header.
#' @param chain_files vector of string paths for each file of samples
#' @param burn the number of rows you want to discard as burn in. Defaults to 0.
#' @param ... arguments to be passed in to lapply
#' @keywords rhat Gelman-Rubin
#' @export
#' @examples
#' files <- c('/fake/path/samps1.csv', '/fake/path/samps2.csv')
#' rhats(files, burn = 300)
rhats <- function(chain_files, burn = 0, ...){
  # TODO: option for discarding burnin (or maybe put this in the ...)
  # TODO: check all files are hte same shape
  # TODO: print the first couple of rows to make sure they look ok
  if(burn > 0){
    dfs <- lapply(chain_files,
                  function(path) read.csv(path, header=T)[-1:-burn,],
                  ...)
  }else{
    dfs <- lapply(chain_files,
                  function(path) read.csv(path, header=T),
                  ...)
  }
  m <- 2*length(chain_files)
  n <- nrow(dfs[[1]])
  n <- (n %/% 2)*2 # in case n is odd
  dfs_halved <- lapply(dfs, function(df) df[((n/2+1):n),])
  dfs_halved <- c(dfs_halved, lapply(dfs, function(df) df[(1:(n/2)),]))
  n <- n/2
  rm(dfs)
  gc(dfs)

  # calculate all rhat statistics for each untransformed univariate parameter
  num_params <- ncol(dfs_halved[[1]])
  rhats <- vector(mode="numeric", length = num_params)
  for(param_num in 1:num_params){
    univariate_chains <- sapply(dfs_halved, '[', param_num) # pick out the same col from each df
    psi_bar_dot_js <- sapply(univariate_chains, mean)
    psi_bar_dot_dot <- mean(psi_bar_dot_js)
    ss_j <- sapply(univariate_chains, var)
    w <- mean(ss_j)
    b <- var(psi_bar_dot_js)*n
    varhat_plus <- w*(n-1)/n + b/n
    rhats[param_num] <- sqrt(varhat_plus/w)
  }
  return(rhats)
}


#' A function that removes rows before displaying pairwise scatterplots.
#'
#' If pairs is choking on too much data, this function will remove rows,
#' and then call it. It will also plot two different colors; one color will
#' be given to the first half of all the samples, and the other color will
#' be given to the second half of all the samples.
#' @param df a data frame of samples
#' @param every take every "every"rd row
#' @keywords pairs scatterplot pairwise
#' @export
#' @examples
#' smart_pairs(fake_df, 100)
smart_pairs <- function(df, every){
  n <- nrow(df)
  rows <- seq(1, n, by = every)
  pairs(df[rows,],
        col = rep(seq(1,2),
                  each = length(rows)/2))
}
