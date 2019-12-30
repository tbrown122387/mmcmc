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
  only_one_column <- is.null(nrow(dfs[[1]]))
  m <- 2*length(chain_files)
  if(only_one_column){
    n <- length(dfs[[1]])
  }else {
    n <- nrow(dfs[[1]])
  }
  n <- (n %/% 2)*2 # in case n is odd
  if(only_one_column){
    dfs_halved <- lapply(dfs, function(df) df[((n/2+1):n)])
    dfs_halved <- c(dfs_halved, lapply(dfs, function(df) df[(1:(n/2))]))
  }else {
    dfs_halved <- lapply(dfs, function(df) df[((n/2+1):n),])
    dfs_halved <- c(dfs_halved, lapply(dfs, function(df) df[(1:(n/2)),]))
  }
  n <- n/2
  if(only_one_column){
    num_params <- 1
  }else {
    num_params <- ncol(dfs_halved[[1]])
  }
  rm(dfs)

  # calculate all rhat statistics for each untransformed univariate parameter
  rhats <- vector(mode="numeric", length = num_params)
  for(param_num in 1:num_params){
    if(only_one_column)
      univariate_chains <- dfs_halved
    else
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


#' A function that helps you elicits inverse gamma prior parameters.
#'
#' This function takes some analyst specifications and returns the
#' "best" parameter tuple that satisfies your criteria. The analyst
#' specifies an ideal average for the value, and an interval (along
#' with a given confidence percentage). If there are tuples that
#' possess, both, coverage probabilities below and above the given
#' percentage, this function will return the tuple that yields the
#' given mean, and that is closest to the desired percentage. If there
#' is no "crossing" then NULL is returned. Optionally, this function
#' will plot the coverage probabilities.
#' @param ideal_mean desired mean for your random variable
#' @param percent_ci what percent confidence do you want (e.g. .95)
#' @param lower_ci what is the lower bound on your interval?
#' @param upper_ci what is the upper bound on your interval?
#' @param beta_grid_left what is the lowest beta value for your grid?
#' @param beta_grid_right what is the highest beta value for your grid?
#' @param beta_increment how spaced apart do you want each beta grid element?
#' @param plot whether a plot should be given as a side effect
#' @keywords prior elicit elicitation inverseGamma
#' @export
#' @examples
#' smart_pairs(fake_df, 100)
elicitInvGamma <- function(ideal_mean, percent_ci, lower_ci,
                           upper_ci, beta_grid_left, beta_grid_right,
                           beta_increment, plot = T){
  beta <- seq(beta_grid_left, beta_grid_right, beta_increment)
  #alpha <- ideal_mean*beta/(1-ideal_mean)
  alpha <- (beta + ideal_mean)/ideal_mean
  prob_between <- invgamma::pinvgamma(upper_ci, shape = alpha, rate = beta) -
    invgamma::pinvgamma(lower_ci, shape = alpha, rate = beta)

  # if plot
  if(plot){
    plot(beta, prob_between, type ="l")
    abline(h=.95, col = "red")
  }

  le_df <- data.frame(alpha, beta, prob_between)
  is_cross <- any(le_df[,3] < percent_ci) && any(le_df[,3] > percent_ci)

  if(is_cross){
    le_df$abs_diff <- abs(le_df[,3] - percent_ci)
    final_results <- le_df[le_df$abs_diff == min(le_df$abs_diff),]
    rownames(final_results) <- c()
    return(final_results)
  }else{
    return(NULL)
  }
}


#' A function that helps you elicits beta prior parameters.
#'
#' This function takes some analyst specifications and returns the
#' "best" parameter tuple that satisfies your criteria. The analyst
#' specifies an ideal average for the value, and an interval (along
#' with a given confidence percentage). If there are tuples that
#' possess, both, coverage probabilities below and above the given
#' percentage, this function will return the tuple that a.) yields the
#' given mean, and b.) that is closest to the desired percentage. If there
#' is no "crossing" then NULL is returned. Optionally, this function
#' will plot the coverage probabilities.
#' @param ideal_mean desired mean for your random variable
#' @param percent_ci what percent confidence do you want (e.g. .95)
#' @param lower_ci what is the lower bound on your interval?
#' @param upper_ci what is the upper bound on your interval?
#' @param beta_grid_left what is the lowest beta value for your grid?
#' @param beta_grid_right what is the highest beta value for your grid?
#' @param beta_increment how spaced apart do you want each beta grid element?
#' @param plot whether a plot should be given as a side effect
#' @keywords prior elicit elicitation inverseGamma
#' @export
#' @examples
#' smart_pairs(fake_df, 100)
elicitBeta <- function(ideal_mean, percent_ci, lower_ci,
                           upper_ci, beta_grid_left, beta_grid_right,
                           beta_increment, plot = T){
  beta <- seq(beta_grid_left, beta_grid_right, beta_increment)
  alpha <- ideal_mean*beta/(1-ideal_mean)
  prob_between <- pbeta(upper_ci, shape1 = alpha, shape2 = beta) -
    pbeta(lower_ci, shape1 = alpha, shape2 = beta)

  # if plot
  if(plot){
    plot(beta, prob_between, type ="l")
    abline(h=.95, col = "red")
  }

  le_df <- data.frame(alpha, beta, prob_between)
  is_cross <- any(le_df[,3] < percent_ci) && any(le_df[,3] > percent_ci)

  if(is_cross){
    le_df$abs_diff <- abs(le_df[,3] - percent_ci)
    return(le_df[le_df$abs_diff == min(le_df$abs_diff),])
  }else{
    return(NULL)
  }
}


#' A function that evaluates the log-likelihood for a logistic regression model.
#'
#' This function evaluates the log-likelihood of a logistic regression model
#' in a numerically stable way. It uses the log-sum-exp trick.
#' @param y observed dependent variables, assumed to be coded as {0,1}.
#' @param linkedMeans a numeric vector. The model matrix times the beta vector.
#' @keywords prior elicit elicitation inverseGamma
#' @export
#' @examples
#' xiB <- rnorm(length(y), mean = 100)
#' logLogisticRegCndtlLike(y, xiB)
#' sum(dbinom(y, 1, inv.logit(xiB), TRUE)) #-Inf
logLogisticRegCndtlLike <- function(y, linkedMeans){
  stopifnot(is.vector(linkedMeans))
  stopifnot(is.vector(y))
  sum(
    mapply(function(y,xiTransposeBeta){ -matrixStats::logSumExp(c(0, (1-2*y)*xiTransposeBeta)) },
           y,
           linkedMeans))
}


#' A function that produces a 3D plot.
#'
#' This function produces a 3D plot. x,y are the independent variables,
#' and the z-axis is the dependent.
#' @param lowerFirst lower bound of x axis
#' @param upperFirst upper bound of x axis
#' @param lowerSecond lower bound of y axis
#' @param upperSecond upper bound of y axis
#' @param numGridPointsOnEachAxis how many grid points do you want on each axis
#' @param f the function that takes two scalar arguments (x and y) and produces one scalar argument (z)
#' @param contour do you want a contour plot? (True or False)
#' @param ... extra arguments to be passed to graphics::contour() or graphics::persp() (depending on what contour arg was set to)
#' @keywords plotting 3D 3-D 3d
#' @export
#' @examples
#' plotSurface(-50, 50, 0.0001, 50, 20, eval_log_unnormalized_posterior, F, theta=-120, zlab = "log unnorm dens", xlab = "mu", ylab = "ss")
plotSurface <- function(lowerFirst, upperFirst, lowerSecond, upperSecond,
                        numGridPointsOnEachAxis, f, contour = F, ...)
{
  A <- seq(lowerFirst, upperFirst, length.out = numGridPointsOnEachAxis)
  B <- seq(lowerSecond, upperSecond, length.out = numGridPointsOnEachAxis)
  args <- expand.grid(A,B)
  z <- mapply(f, args[,1], args[,2])
  dim(z) <- c(length(A), length(B))
  if(contour){
    graphics::contour(A, B, z, ...)
  }else{
    graphics::persp(x=A, y=B, z=z, ...)
  }
}


#' A function that produces raw assignment code for a given object.
#'
#' A function that produces raw assignment code that you can copy/paste
#' into a script. The idea is to make it more reproducible.
#' @param d the object you have sitting in memory
#' @param dest the file path destination that the code will get written to
#' @param objectName the desired name of the object produced by the code.
#' Defaults to NULL, which means the name will get generated from how you
#' have the variable stored in your interactive session.
#' @keywords code generation
#' @export
genAssigntmentCode <- function(d, dest, objectName){

  if(is.matrix(d)){
    myString <- paste0(objectName, " <- matrix(c(")
    for(i in 1:(nrow(d)-1)){
      myString <- paste0(myString, paste(d[i,], collapse=", "), ",\n")
    }
    myString <- paste0(myString, paste(d[nrow(d),], collapse=", "))
    myString <- paste0(myString, ")")
    myString <- paste0(myString, ", nrow=", nrow(d),")")

  }else if(is.numeric(d)){
    myString <- paste0(objectName, " <- c(")
    myString <- paste0(myString, paste(d, collapse=", "))
    myString <- paste0(myString, ")")

  }else{
    stop("no implementation for this type of object")
  }

  # write everything out
  fileConn<-file(dest)
  writeLines(myString, fileConn)
  close(fileConn)

}


