
softmax <- function(x) exp(c(0, x)) / sum(exp(c(0, x)))



##############################################
# complete-data log-likelihoods              #
##############################################

# for imputation with method 1
LLimp1 <- function(beta, D, x, rimp, lambda)
{
  b     <- c(-log(sum(exp(D[, -1] %*% beta))), beta)
  mu    <- as.vector(exp(D %*% b))
  num   <- t(mu * rimp) * x
  x_y   <- crossprod(mu, rimp)

  as.vector(log(x_y) %*% x - lambda * sum(beta^2))
}

# for imputation with method 2
LLimp2 <- function(beta, D, z, y, rimp, lambda)
{
  b     <- c(-log(sum(exp(D[, -1] %*% beta))), beta)
  mu    <- as.vector(exp(D %*% b))
  zhat  <- numeric(length(y))
  for (i in 1:length(z))
    zhat[rimp[[i]]] <-  zhat[rimp[[i]]] + z[i] * mu[rimp[[i]]] / sum(mu[rimp[[i]]])
  
  cpz <- numeric(length(rimp))
  for (i in 1:length(rimp))
    cpz[i] <- sum(mu[rimp[[i]]])
  
  as.vector(sum(y * log(mu)) + sum(z * log(cpz)) - lambda * sum(beta^2))
}

# for lca

LLlca <- function(beta, D, y, nclasses, robs, lambda)
{
  b    <- c(-log(sum(exp(D[, -1] %*% beta))), beta)
  mu   <- as.vector(exp(D %*% b))
  ymat <- matrix(mu, length(y), nclasses)
  yhat <- rowSums(ymat)
  
  sum(y %*% log(yhat)) - lambda * sum(beta^2)
}

# for lca and imp for method 1

LLimp1lca <- function(beta, D, x, n, rimp, lambda)
{
  b     <- c(-log(sum(exp(D[, -1] %*% beta))), beta)
  mu    <- as.vector(exp(D %*% b)) 
  num   <- t(mu * rimp) * x 
  x_y   <- crossprod(mu, rimp)
  
  as.vector(log(x_y) %*% x) - lambda * sum(beta^2)
}


# for lca and imp for method 2

LLimp2lca <- function(beta, D, y, z, nclasses, rimp, lambda)
{
  b      <- c(-log(sum(exp(D[, -1] %*% beta))), beta) 
  pcomp  <- exp(D %*% b)
  mprobs <- matrix(pcomp, length(y), nclasses)
  rsum   <- rowSums(mprobs)
  prmiss <- NULL
  for (i in 1:length(rimp))
    prmiss <- c(prmiss, sum(mprobs[rimp[[i]]]))
  
  sum(c(y, z) * log(c(rsum, prmiss))) - lambda * sum(beta^2)
}

###########################################################
# Compute (conditional) probabilities of latent variables #
###########################################################

get_probs <- function(dx, covs, lists, lat)
{
  nvars <- ncol(dx) - 1
  dx    <- mutate(dx, across(any_of(!!{c(lists, covs)}), as.factor),
                  fitted = fitted / sum(fitted)) 
    
  probs <- list()
  for(j in 1:length(lat)){
    probs[[paste0("P(", colnames(dx[j]), ")")]] <- round(xtabs(paste("fitted ~", colnames(dx[j])), data = dx), 4)
    q     <- NULL
    names <- NULL
    for (k in (length(lat) + 1):(ncol(dx) - 1)){
      form  <- paste("fitted ~ ", colnames(dx)[k], "+", colnames(dx)[j])
      q     <- round(rbind(q, xtabs(form, dx)/matrix(colSums(xtabs(form, dx)),
                                                     nlevels(dx[, k]),
                                                     nlevels(dx[, j]),
                                                     byrow = T)), 4)
      names <- c(names, rep(colnames(dx)[k], nlevels(dx[, k])))
    }
    rownames(q) <- paste(names, "=",  rownames(q))
    colnames(q) <- paste(colnames(dx)[j], "=", 1:nlevels(dx[, j]))
    probs[[paste0("P(.|", colnames(dx)[j],")")]] <- q
  }
  order      <- order(probs[[1]], decreasing = T)
  probs[[1]] <- probs[[1]][order]
  probs[[2]] <- probs[[2]][, order]
  
  qq <- rr <- NULL
  
  for (i in 1:length(covs)) {
    qq <- c(qq, which(startsWith(rownames(probs[[2]]), covs[i])))
  }
  
  for (i in 1:length(lists)) {
    rr <- c(rr, which(startsWith(rownames(probs[[2]]), lists[i])))
  }
  probs[[2]] <- probs[[2]][c(qq, rr), ]
  
  probs
}

#################################################
# obtain starting values for b                  #
#################################################

Newt <- function(D, x, n)
{
  b  <- c(-log(nrow(D)) + log(n), rep(0, ncol(D) - 1))
  for (i in 1:10) {
    mu   <- as.vector(exp(D %*% b))
    db   <- crossprod(D, x) - crossprod(D, mu) - 2 * c(0, b[-1])
    I    <- crossprod(D, Diagonal(x = mu)) %*% D + 2 * Diagonal(x = c(0, rep(2, length(b) - 1)))
    H    <- chol2inv(chol(I))
    b    <- b + .1 * H %*% db
    b[1] <- -log(sum(exp(D[, -1] %*% b[-1]))) + log(n)
  }
  b
}


#################################################
# Newton for incomplete data                    #
#################################################

NewtRaph <- function(beta, D, comp, n, pen, dpen, mtol)
{
  conv <- F
  b    <- c(-log(sum(exp(D[, -1] %*% beta))) + log(n), beta)
  while (!conv) {
    mu   <- as.vector(exp(D %*% b))
    db   <- crossprod(D, comp) - crossprod(D, mu) - pen * b
    I    <- crossprod(D, Diagonal(x = mu)) %*% D + dpen
    H    <- chol2inv(chol(I))
    bt   <- b + H %*% db
    conv <- sum(abs(b - bt)) < mtol
    b    <- bt
  }
  bt
}


##########################
# print output to screen #
##########################

output <- function(call, mse, imp, lca, lambda, iter, conv, n, LL, pse,
                   seed, npar, df, dev, AIC, BIC, dpop, pars, probs, printPars)
{
  print(call)
  cat("\n\n")
  cat("MODEL \n")
  cat("=======================================================")
  cat("\n")
  cat("Multiple systems      =", mse)
  cat("\n")
  cat("Imputation missings   =", imp)
  cat("\n")
  cat("Latent classes        =", lca)
  cat("\n")
  cat("Lambda                =", lambda)
  cat("\n\n")
  cat("STATISTICS \n")
  cat("=======================================================")
  cat("\n")
  cat("Seed                   =", seed)
  cat("\n")
  cat("Number of iterations   =", iter)
  cat("\n")
  cat("Convergence            =", conv)
  cat("\n")
  cat("Sample size            =", n)
  cat("\n")
  cat("Log-likelihood         =", LL)
  cat("\n")
  cat("Number of parameters   =", npar)
  cat("\n")
  cat("Degrees of freedom     =", df)
  cat("\n")
  cat("Deviance               =", ifelse(dev < 1e-6, 0, dev))
  cat("\n")
  cat("AIC                    =", AIC)
  cat("\n")
  cat("BIC                    =", BIC)
  cat("\n\n")
  cat("PARAMETER ESTIMATES \n")
  cat("=======================================================")
  cat("\n\n")
  if (mse){
    cat("population size estimate:\n")
    print(pse)
    cat("\n\n")
  }
  if (printPars){
  cat("Log-linear parameters")
  cat("\n\n")
  print(pars %>% round(4))
  cat("\n\n")
  }
  if (!is.null(probs))  {
    cat("Latent class probabilities")
    cat("\n\n")
    print(probs)
  }
}

###############################################
# Include latent variable(s) in complete data #
###############################################

make_dx <- function(dy, nclasses, lat, nlat, latlevs, numvars)
{
  dx <- NULL
  
  for (j in 1:nclasses)
    dx    <- rbind(dx, dy)
  
  X     <- matrix(0, nrow = nrow(dx), ncol = nlat, dimnames = list(NULL, sort(lat)))
  
  for (i in 1:nlat)
    X[, i]    <- rep(rep(1:latlevs[i], each = nrow(X) / prod(latlevs[1:i])), times = i)
  
  cbind(mutate(as.data.frame(X), across(everything(), .fns = as.factor)), dx) %>%
    mutate(across(!!{numvars}, as.numeric))
}

###############################################
# Include latent variable(s) in complete data #
###############################################

make_dxpop <- function(dpop, nclasses, lat, nlat, latlevs, numvars)
{
  dlat <- NULL
  
  for (j in 1:nclasses)
    dlat    <- rbind(dlat, dpop)
  
  X     <- matrix(0, nrow = nrow(dlat), ncol = nlat, dimnames = list(NULL, sort(lat)))
  
  for (i in 1:nlat)
    X[, i]    <- rep(rep(1:latlevs[i], each = nrow(X) / prod(latlevs[1:i])), times = i)
  
  s <- cbind(mutate(as.data.frame(X), across(everything(), .fns = as.factor)), dlat) %>%
    mutate(across(!!{numvars}, as.numeric),
           Freq = 0) 
}


# deviance of cases with missing values in z and 
# number of groups with identical missing profiles

devz <- function(dmz, z, cpz)
{
  grps <- vector("list", nrow(dmz))
  
  for (i in 1:nrow(dmz))
  {
    obsval    <- dmz[i, !is.na(dmz[i, ]), drop = F]
    grps[[i]] <- colnames(obsval)
  }
  grp <- unique(grps)
  
  gr <- vector("list", length(grp))
  
  for (i in 1:length(grp)){
    for (j in 1:length(grps)){
      if (identical(grp[[i]], grps[[j]])){
        gr[[i]] <- c(gr[[i]], j)
      }
    }
  }
  devz <- NULL
  for (i in 1:length(gr))
  {
    O    <- z[gr[[i]]]
    E    <- sum(z[gr[[i]]]) * cpz[gr[[i]]]
    devz <- c(devz, 2 * sum(O * log(O / E), na.rm = T))
  }
  list(sum(devz), length(grp))
}

