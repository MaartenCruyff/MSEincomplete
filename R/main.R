#' Fit log linear models for incomplete data
#'
#' @description `loglinear` performs a log linear analyses with optionally 
#' missing values imputation, and/or latent class analyses and/or multiple 
#' systems estimation.
#'
#' @param formula formula of the model to be fitted. See 'Details'.
#' @param data data frame containing the manifest variables in the model. See 'Details'.
#' @param lambda parameter for the ridge penalty. See 'Details'.
#' @param latlevs vector with the number of levels of the latent 
#' variable(s) in the formula, if any. Defaults to 2 levels.
#' @param stderr computation method for the standard errors. See 'Details'. 
#' @param method imputation method for missing values. Default option 1 uses 
#' matrix algebra, option 2 a loop.
#' @param seed optional seed for the starting values.
#' @param printPars logical for printing the log-linear parameter estimates to 
#' the screen.
#' @param control list with control parameters:
#' * `maxit` the maximum number of iterations of the EM algorithm. Defaults to 2000.
#' * `tol` convergence criterion for the incomplete data log likelihood. Defaults to 1e-5.
#' * `mtol` convergence criterion for the complete data log likelihood in the M-step. 
#' Defaults to 1e-5
#' * `vlat` parameter determining the variance of the random starting values for the latent
#' class probabilities. The value 0 yields approximately equal probabilities, other 
#' values increase the variance.
#' @return A list with the following components:
#' * `formula` the formula of the fitted model.
#' * `lists` names of the population registers, if present.
#' * `lat` names of the latent variables, if present.
#' * `dobs` contingency table of the observed data, including sampling zeros.
#' * `dcmp` contingency table of the complete data estimates.
#' * `pars` parameter estimates, including SE's and z- and p-values.
#' * `probs` tables with latent class probabilities, if applicable.
#' * `hist` iteration history.
#'
#' @details
#' The role of the variables in the formula is determined by their names.  
#' Variables starting with an upper case letter are interpreted as lists, except 
#' for variables named `X`, `Y` and `Z`, which are interpreted  as latent variables. If 
#' there is more than one list in the formula, MSE is performed. The imputation of 
#' missing values is done automatically.
#' 
#' The data can be provided in the form of a contingency table or as individual records.
#' In the former case sampling zeros are automatically included in the table.
#' 
#' The parameter `lambda` can be set if the algorithm breaks down due to a 
#' non-invertible hessian. A value as small as 1e-5 is usually sufficient to 
#' overcome this problem.    
#' 
#' For models with latent classes and/or imputation of missing values the default 
#' option "incomplete" of the argument `stderr` correctly uses the hessian of the incomplete data 
#' log likelihood. For complex models, however, this method is painstakingly slow, 
#' In that case it is advisable to use the option "complete", which underestimates 
#' the standard errors but is much faster.
#' 
#' @examples
#' # An MSE model with latent classes and missing value imputation for a, b, and c.
#' loglinear(formula = Freq ~ X * (a + b + c) + A + B + C, data = lcaimp)
#' 
#'
#' @importFrom stats complete.cases model.matrix na.omit optim pnorm xtabs
#' glm poisson runif update terms coef formula dpois optim rnorm rmultinom
#' @importFrom dplyr mutate left_join %>% across select filter everything
#' filter rename all_of any_of case_when
#' @importFrom Matrix Matrix Diagonal crossprod tcrossprod t
#' @importFrom pracma hessian 
#' @importFrom rlang :=
#' @export

loglinear <- function(formula, data, lambda = 0, latlevs = NULL, method = 1,
                      stderr  = c("incomplete", "complete"), seed = NULL,
                      printPars = TRUE, control   = list())
{

  call  <- match.call()

  if (is.null(seed))
      seed <- sample(1e+5, 1)

  maxit    <- ifelse(exists("maxit", control), control$maxit, 2000)
  tol      <- ifelse(exists("tol", control), control$tol, 1e-5)
  mtol     <- ifelse(exists("mtol", control), control$mtol, 1e-5)
  vlat     <- ifelse(exists("vlat", control), control$vlat, 2)
  
  # update . in rhs of formula ~ . to include all terms
  terms    <- terms(x = formula, data = data)
  formula  <- formula(terms(x = formula, data = data))
  
  # check if data is contingency table 
  ctable   <- attr(terms, "response") == 1
  
  #if so, change name of response to "Freq", update formula and get other variables
  if (ctable) {
    response <- all.vars(formula)[1]
    data     <- rename(data, Freq = !!{response})
    formula  <- update(formula, Freq ~ .)
    vars     <- all.vars(formula)[-1]
  } else {
    vars     <- all.vars(formula)
  }
  
  # get vectors with the lists, latent and manifest variables 
  lists    <- c()
  for (i in 1:26) {
    for (j in 1:length(vars)){
      if (startsWith(vars[j], LETTERS[i])) lists <- c(lists, vars[j])
    }
  }
  nlists   <- length(lists)
  
  lat  <- vars[vars %in% LETTERS[24:26]]
  nlat     <- length(lat)
  if (is.null(latlevs) & nlat > 0)  # set the latent variable levels to the default value
    latlevs <- rep(2, nlat)
  
  lists    <- setdiff(lists, lat)
  covs     <- setdiff(vars, c(lists, lat))
  ncovs    <- length(covs)
  
  # select variables in formula
  if (ctable){
    if (nlists > 0 & ncovs > 0){
      data     <- data[, c(lists, covs, "Freq")] %>% mutate(across(!!{{lists}}, as.factor))
    } else if (nlists > 0 & ncovs == 0){
      data     <- data[, c(lists, "Freq")] %>% mutate(across(!!{{lists}}, as.factor))
    } else if (nlists == 0 & ncovs > 0) {
      data     <- data[, c(covs, "Freq")]
    }
  } else {
    if (nlists > 0 & ncovs > 0){
      data     <- data[, c(lists, covs)]
    } else if (nlists > 0 & ncovs == 0){
      data     <- data[, c(lists)]
    } else if (nlists == 0 & ncovs > 0) {
      data     <- data[, c(covs)]
    }
  }
  
  # determine which covariates are numeric before making the complete contingency table
  numvars <- colnames(data)[sapply(data, class) %in% c("integer", "numeric")] 
  
  # make the complete contingency table
  if (ctable){
    d0 <- as.data.frame(xtabs(Freq ~ ., data = data, addNA = TRUE)) %>%
      filter(!(!complete.cases(.) & Freq == 0))
  } else {
    d0 <- as.data.frame(xtabs( ~ .,     data = data, addNA = TRUE)) %>%
      filter(!(!complete.cases(.) & Freq == 0))
    formula <- update(formula, Freq ~ .)
  }
  
  # convert the numeric covariates back to numeric
  d0   <- mutate(d0, across(!!{numvars}, as.numeric))
  
  # reorder the data complete followed by incomplete
  if (nrow(d0) > sum(complete.cases(d0))) {
    d0 <- rbind(na.omit(d0), d0[!complete.cases(d0), ])
  }
  
  # determine which model is to be fitted
  mse  <- length(lists) > 1
  lca  <- nlat > 0
  imp  <- sum(!complete.cases(d0)) > 0


  ##########
  # models #
  ##########

  if (!imp & !lca)
    ret <- fcomp(formula = formula, d0 = d0, lambda = lambda, maxit = maxit, 
                 tol = tol, lists = lists, mse = mse, imp = imp, lca = lca,
                 printPars = printPars, call = call)
  
  if (imp & !lca & method == 1)
    ret <- imp1(formula = formula, d0 = d0, lambda = lambda,
                stderr = stderr, maxit = maxit, tol = tol, mtol = mtol,
                seed = seed, lists = lists,
                mse = mse, lca = lca, imp = imp,
                printPars = printPars, call = call)

  if (imp & !lca & method == 2)
    ret <- imp2(formula = formula, d0 = d0, lambda = lambda,
                stderr = stderr, maxit = maxit, tol = tol, mtol = mtol,
                seed = seed, lists = lists,
                mse = mse, lca = lca, imp = imp,
                printPars = printPars, call = call)
  
  if (lca & !imp)
    ret <- lca2(formula = formula, d0 = d0, latlevs = latlevs, 
                lambda = lambda, stderr = stderr,  maxit = maxit, tol = tol,
                mtol = mtol, seed = seed, nlat = nlat, lat = lat, lists = lists,
                numvars = numvars, covs = covs, vlat = vlat,
                mse = mse, lca = lca, imp = imp,
                printPars = printPars, call = call)

  if (lca & imp & method == 1)
    ret <- imp1lca(formula = formula, d0 = d0, latlevs = latlevs, 
                   lambda = lambda, stderr = stderr,  maxit = maxit, tol = tol,
                   mtol = mtol, seed = seed, nlat = nlat, lat = lat, lists = lists,
                   numvars = numvars, covs = covs, vlat = vlat,
                   mse = mse, lca = lca, imp = imp,
                   printPars = printPars, call = call)
  
  if (lca & imp & method == 2)
    ret <- imp2lca(formula = formula, d0 = d0, latlevs = latlevs, 
                   lambda = lambda, stderr = stderr,  maxit = maxit, tol = tol,
                   mtol = mtol, seed = seed, nlat = nlat, lat = lat, lists = lists,
                   numvars = numvars, covs = covs, vlat = vlat, 
                   mse = mse, lca = lca, imp = imp,
                   printPars = printPars, call = call)
  
  invisible(ret)

}



#' Bootstrap for MSE models
#'
#' @description Non-parametric bootstrap of the parameter estimates and the fitted 
#' frequencies. 95% confidence intervals are computed with the percentile method 
#' and bias-corrected accelerated bootstrap (BCa) method.
#'
#' @param object object created with \code{\link{loglinear}}.
#' @param B number of bootstrap samples.
#' @param lambda parameter for the ridge penalty \eqn{\lambda\sum\beta^2}.
#' @param cores the number of cores used for parallel computing.
#' @param seed optional seed for reproducibility.
#' @param control a list with the parameters:
#' * `maxit` maximum number of iteration of the EM algorithm.
#' * `tol` stop criterion for maximization of the incomplete data log likelihood.
#' * `mtol` stop criterion for maximization of the complete data log-likelihood.
#' @return a list with the components:
#' * `par` the parameter estimates of the bootstrap sample
#' * `fitted` the fitted population frequencies
#' @importFrom stats quantile sd qnorm
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores clusterSetRNGStream
#' @importFrom foreach foreach %dopar% 
#' @export

bootmse <- function(object, B = 500, lambda = 0, cores = 2, seed = NULL, control = list())
{
  if (length(object$lists) %in% 0:1) stop("Bootstrap only implemented for MSE models")
  
  
  
  cl      <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (!is.null(seed))
    clusterSetRNGStream(cl = cl, iseed = seed)
  
  
  if (!exists("lat", object) & !any(is.na(object$dobs))) {
    
    bootres <- foreach(i = 1:B, .combine = cbind, .inorder = F,
                       .errorhandling = "remove", .packages = "MSEincomplete") %dopar% {
                         bootcomp(object = object, lambda = lambda, control = control)
                       }
    
    bootfits <- bootres[which(rownames(bootres) == 1):nrow(bootres), ]
    bootpars <- bootres[1:(which(rownames(bootres) == 1) - 1),  ]
    ci95     <- cicomp(object = object, lambda = lambda, pseb = colSums(bootfits), control = control)
    
  } else if (exists("lat", object) & !any(is.na(object$dobs))) {
    
    bootres <- foreach(i = 1:B, .combine = cbind, .inorder = F, 
                       .errorhandling = "remove", .packages = "MSEincomplete") %dopar% {
                         bootlca(object = object, lambda = lambda, control = control)
                       }
    bootfits <- bootres[which(rownames(bootres) == 1):nrow(bootres), ]
    bootpars <- bootres[1:(which(rownames(bootres) == 1) - 1),  ]
    ci95     <- cilca(object = object, lambda = lambda, pseb = colSums(bootfits), control = control)
    
    
  } else if (!exists("lat", object) & any(is.na(object$dobs))) {
    
    bootres <- foreach(i = 1:B, .combine = cbind, .inorder = F, 
                       .errorhandling = "remove", .packages = "MSEincomplete") %dopar% {
                         bootimp(object = object, lambda = lambda, control = control)
                       } 
    
    bootfits <- bootres[which(rownames(bootres) == 1):nrow(bootres), ]
    bootpars <- bootres[1:(which(rownames(bootres) == 1) - 1),  ]
    ci95     <- ciimp(object = object, lambda = lambda, pseb = colSums(bootfits), control = control)
    
    
  } else if (exists("lat", object) & any(is.na(object$dobs))) {
    
    bootres <- foreach(i = 1:B, .combine = cbind, .inorder = F, 
                       .errorhandling = "remove", .packages = "MSEincomplete") %dopar% {
                         bootimplca(object = object, lambda = lambda, control = control)
                       }
    bootfits <- bootres[which(rownames(bootres) == 1):nrow(bootres), ]
    bootpars <- bootres[1:(which(rownames(bootres) == 1) - 1),  ]
    ci95     <- ciimplca(object = object, lambda = lambda, pseb = colSums(bootfits), control = control)
    
  }
  
  stopCluster(cl)
  
  invisible(list(ci95 = ci95, bootfits = bootfits, bootpars = bootpars))
  
}


#' Data simulation
#'
#' @description Generates a contingency table for analysis with the function
#' \code{loglinear}.
#'
#' @param lists list with the levels 0,1 for the registers. If a register 
#' contains missing, add \code{NA}. \code{NULL} if the data has no registers.
#' @param covs list with levels of the covariates, If a covariate contains 
#' missings add \code{NA}.\code{NULL} if the data has no covariates.
#' @param N population/sample size, excluding the observations with missing values.
#' @param mean,sd parameters of the normal distribution from which the loglinear
#' parameters are drawn. See 'Details'
#' @param seed optional seed for reproducibility.
#'
#' @return A contingency table
#'
#' @details \code{mean} and \code{sd} apply to the lists
#'
#' @importFrom stats rmultinom as.formula model.matrix rnorm glm poisson
#' @importFrom dplyr filter mutate
#' @importFrom rlang :=
#' @export

simdat  <- function(lists = list(0:1, 0:1, 0:1), 
                    covs = list(c(1, 2, NA), c(1, 2, NA), c(1:2, NA)), 
                    N = 1000,
                    mean = -2, sd = .5, seed = NULL)
{

  if (!is.null(seed)) set.seed(seed)

  levs <- list()

  if (!is.null(lists)){
    nlists <- length(lists)
  for (i in 1:nlists)
    levs[[LETTERS[i]]] <- as.factor(lists[[i]])
  } else {
    nlists <- 0
  }

  if (!is.null(covs)) {
    for (i in (nlists + 1):(nlists + length(covs)))
      levs[[letters[i - nlists]]] <- as.factor(covs[[i - nlists]])
  } else {
    ncovs <- 0
  }

  dpop        <- expand.grid(lapply(1:length(levs), function(i)levs[[i]]), KEEP.OUT.ATTRS = F)
  names(dpop) <- names(levs)

  f   <- as.formula(paste0(" ~ .^", ncol(dpop)))
  D   <- model.matrix(f, dpop)

  if (!is.null(lists)) {
    b <- c(rnorm(nlists, mean = mean, sd = sd),
           rnorm(ncol(D) - nlists - 1, sd = sd / 2))
  } else {
    b <- rnorm(ncol(D) - 1)
  }

  b     <- c(-log(sum(exp(D[, -1 ] %*% b)) + log(N)), b)

  if (nrow(dpop) == sum(complete.cases(dpop))){
    dpop$Freq  <- c(rmultinom(1, N, c(exp(D %*% b))))
  } else {
    dy    <- na.omit(dpop) 
    dz    <- dpop[!complete.cases(dpop), ]
    Freq  <- c(rmultinom(1, N, c(exp(D %*% b))),
               sample(10:0, size = nrow(dz), T, prob = softmax(seq(1, 5, length.out = 10))))
    dpop  <- rbind(dy, dz) %>% mutate(Freq = Freq) %>% 
      filter(!(!complete.cases(.) & Freq == 0))
  }

  dpop
}

