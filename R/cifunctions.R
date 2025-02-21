########################
# acceleration for BCa #
########################

# for MSE data without latent classes and missing values

cicomp <- function(object, lambda, pseb, control)
{
  
  maxit    <- ifelse(exists("maxit", control), control$maxit, 2000)
  tol      <- ifelse(exists("tol", control), control$tol, 1e-5)
  mtol     <- ifelse(exists("mtol", control), control$mtol, 1e-5)
  
  pse0  <- sum(object$dcmp$fitted)
  y     <- object$dobs$Freq
  n     <- sum(y)
  D     <- model.matrix(object$formula, object$dobs)
  
  Dpop  <- model.matrix(object$formula, object$dcmp %>% select(-fitted))
  
  pen   <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen  <- diag(pen)
  b     <- object$pars[, 1]
  
  psej  <- NULL
  
  for (i in 1:length(y))
  {
    if (y[i] > 0) 
    {
      yj    <- y
      yj[i] <- y[i] - 1
      
      b <- NewtRaph(object$pars[-1, 1], D = D, comp = yj, n = n - 1, pen = pen, dpen = dpen, mtol = mtol)
      
      psej <- c(psej, rep(sum(exp(Dpop %*% b)), y[i]))
    }
  }
  
  I    <- (n - 1) * (mean(psej) - psej)
  a    <- sum(I^3) * sum(I^2)^-1.5 / 6
  
  z0   <- qnorm(mean(pseb < pse0))
  u    <- c(.025, .975)
  zu   <- qnorm(u)
  pbca <- pnorm(z0 + (z0 + zu) / (1 - a * (z0 + zu)))
  
  
  pse   <- sum(object$dcmp$fitted)
  ciPM  <- unname(c(quantile(pseb, u)))
  ciBCa <- unname(c(quantile(pseb, pbca)))
  ci    <- matrix(c(pse, ciPM, pse, ciBCa), 2, 3, byrow = T, 
                 dimnames = list(c("Percentile", "BCa"),
                                 c('pse', '2.5%', '97.5%')))
  
  cat("95% Bootstrap Confidence Intervals \n")
  cat("===============================================\n\n")
  print(ci)

  invisible(ci)
}

# for MSE data without latent classes and with missing values

ciimp <- function(object, lambda, pseb, control)
{
  
  maxit    <- ifelse(exists("maxit", control), control$maxit, 2000)
  tol      <- ifelse(exists("tol", control), control$tol, 1e-5)
  mtol     <- ifelse(exists("mtol", control), control$mtol, 1e-5)
  
  pse0 <- sum(object$dcmp$fitted)
  Dpop <- model.matrix(object$formula, object$dcmp %>% select(-fitted))
  
  d0  <- object$dobs
  n   <- sum(object$dobs$Freq)
  dy  <- na.omit(d0)
  dz  <- d0[!complete.cases(d0), ]
  y   <- dy$Freq
  z   <- dz$Freq
  x   <- c(y, z)
  D   <- Matrix(model.matrix(object$formula, dy), sparse = T)
  
  dmy <- dy %>% select(-Freq) %>% mutate(across(!!{object$lists}, ~ as.numeric(.x) - 1)) %>% data.matrix()
  dmz <- dz %>% select(-Freq) %>% mutate(across(!!{object$lists}, ~ as.numeric(.x) - 1)) %>% data.matrix()
  
  tmp  <- matrix(1, nrow(dmy), nrow(dmz))
  for (k in 1:nrow(dmz))
    for (j in 1:ncol(dmy))
      if (!is.na(dmz[k, j]))
        tmp[, k] <- tmp[, k] * (dmy[, j] == dmz[k, j])
  
  rimp  <- cbind(Matrix::Diagonal(nrow(tmp)), tmp)

  pen   <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen  <- diag(pen)
  b     <- object$pars[, 1]
  
  iter  <- 0
  psej  <- NULL
  
  for (i in 1:length(x))
  {
    if (x[i] > 0) 
    {
      xj    <- x
      xj[i] <- xj[i] - 1
      
      conv <- F
      
      while (!conv)
      {
        mu    <- as.vector(exp(D %*% b))
        num   <- t(mu / (n - 1) * rimp) * xj 
        x_y   <- crossprod(mu / (n - 1), rimp)
        LLnew <- as.vector(log(x_y) %*% xj)
        ycomp <- as.vector(x_y^-1 %*% num)
        b     <- NewtRaph(b[-1], D = D, comp = ycomp, n = n - 1, pen = pen, dpen = dpen, mtol = mtol)
        
        if (iter > 0) {
          delta <- LLnew - LL
          conv  <- LLnew - LL < tol
        }
        LL <- LLnew
        
        iter <- iter + 1
        if (iter == maxit) break
      }
      
      psej <- c(psej, rep(sum(exp(Dpop %*% b)), x[i]))
    }
  }
  
  I    <- (n - 1) * (mean(psej) - psej)
  a    <- sum(I^3) * sum(I^2)^-1.5 / 6
  
  z0   <- qnorm(mean(pseb < pse0))
  u    <- c(.025, .975)
  zu   <- qnorm(u)
  pbca <- pnorm(z0 + (z0 + zu) / (1 - a * (z0 + zu)))
  
  
  pse   <- sum(object$dcmp$fitted)
  ciPM  <- unname(c(quantile(pseb, u)))
  ciBCa <- unname(c(quantile(pseb, pbca)))
  ci    <- matrix(c(pse, ciPM, pse, ciBCa), 2, 3, byrow = T, 
                  dimnames = list(c("Percentile", "BCa"),
                                  c('pse', '2.5%', '97.5%')))
  
  cat("95% Bootstrap Confidence Intervals \n")
  cat("===============================================\n\n")
  print(ci)
  
  invisible(ci)
}


# for MSE data with latent classes and without missing values

cilca <- function(object, lambda, pseb, control)
{
  
  maxit    <- ifelse(exists("maxit", control), control$maxit, 2000)
  tol      <- ifelse(exists("tol", control), control$tol, 1e-5)
  mtol     <- ifelse(exists("mtol", control), control$mtol, 1e-5)
  
  pse0 <- sum(object$dcmp$fitted)
  Dpop <- model.matrix(object$formula, object$dcmp %>% rename(Freq = fitted))
  
  dy  <- object$dobs
  y   <- dy$Freq
  n   <- sum(y)
  
  nclasses <- prod(object$latlevs)
  
  dx    <- make_dx(dy = dy, nclasses = nclasses, lat = object$lat, 
                   nlat = length(object$lat), latlevs = object$latlevs, 
                   numvars = object$numvars)
  
  D     <- Matrix(model.matrix(object$formula, dx), sparse = T)  
  
  xcomp <- numeric(length(y) * nclasses)
  robs  <- matrix(1:length(xcomp), length(y))
  
  pen   <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen  <- diag(pen)
  b     <- object$pars[, 1]
  
  iter  <- 0
  psej  <- NULL
  
  for (i in 1:length(y))
  {
    if (y[i] > 0) 
    {
      yj    <- y
      yj[i] <- yj[i] - 1
      
      conv <- F
      iter <- 0
      
      while (!conv)
      {
        mu    <- as.vector(exp(D %*% b))
        ymat  <- matrix(mu, length(y), nclasses)
        yhat  <- rowSums(ymat)
        LLnew <- sum(yj %*% log(yhat / (n - 1))) - lambda * sum(b[-1]^2)
        
        for (i in 1:length(y))
          xcomp[robs[i, ]]  <- ymat[i, ] * yj[i] / yhat[i]
        
        b    <- NewtRaph(b[-1], D = D, comp = xcomp, n = n - 1, pen = pen, dpen = dpen, mtol = mtol)
        
        if (iter > 0)
          conv <- LLnew - LL < tol
        
        LL <- LLnew
        
        iter <- iter + 1
        if (iter == maxit) break
      }
      
      psej <- c(psej, rep(sum(exp(Dpop %*% b)), y[i]))
    }
  }
  
  I    <- (n - 1) * (mean(psej) - psej)
  a    <- sum(I^3) * sum(I^2)^-1.5 / 6
  
  z0   <- qnorm(mean(pseb < pse0))
  u    <- c(.025, .975)
  zu   <- qnorm(u)
  pbca <- pnorm(z0 + (z0 + zu) / (1 - a * (z0 + zu)))
  
  
  pse   <- sum(object$dcmp$fitted)
  ciPM  <- unname(c(quantile(pseb, u)))
  ciBCa <- unname(c(quantile(pseb, pbca)))
  ci    <- matrix(c(pse, ciPM, pse, ciBCa), 2, 3, byrow = T, 
                  dimnames = list(c("Percentile", "BCa"),
                                  c('pse', '2.5%', '97.5%')))
  
  cat("95% Bootstrap Confidence Intervals \n")
  cat("===============================================\n\n")
  print(ci)
  
  invisible(ci)
}

# for MSE data with latent classes and missing values

ciimplca <- function(object, lambda, pseb, control)
{
  
  maxit    <- ifelse(exists("maxit", control), control$maxit, 2000)
  tol      <- ifelse(exists("tol", control), control$tol, 1e-5)
  mtol     <- ifelse(exists("mtol", control), control$mtol, 1e-5)
  
  nclasses <- prod(object$latlevs)
  
  pse0 <- sum(object$dcmp$fitted)
  Dpop <- model.matrix(object$formula, object$dcmp %>% rename(Freq = fitted))
  
  d0  <- object$dobs
  n  <- sum(d0$Freq)
  dy <- na.omit(d0)
  dz <- d0[!complete.cases(d0), ]
  y  <- dy$Freq
  z  <- dz$Freq
  x  <- c(y, z)
  
  rownames(dy) <- 1:nrow(dy)
  rownames(dz) <- 1:nrow(dz)
  
  dx    <- make_dx(dy = dy, nclasses = nclasses, lat = object$lat, 
                   nlat = length(object$lat), latlevs = object$latlevs, 
                   numvars = object$numvars)
  
  D     <- Matrix(model.matrix(object$formula, dx), sparse = T)  
  
  dmy <- dy %>% select(-Freq) %>% mutate(across(!!{object$lists}, ~ as.numeric(.x) - 1)) %>% data.matrix()
  dmz <- dz %>% select(-Freq) %>% mutate(across(!!{object$lists}, ~ as.numeric(.x) - 1)) %>% data.matrix()
  
  tmp <- matrix(1, nrow(dmy), nrow(dmz))
  
  for (k in 1:nrow(dmz))
    for (j in 1:ncol(dmy))
      if (!is.na(dmz[k, j]))
        tmp[, k] <- tmp[, k] * (dmy[, j] == dmz[k, j])
  
  tmp <- rimp <- cbind(Diagonal(nrow(tmp)), tmp)
  
  for (j in 2:nclasses)
    rimp <- rbind(rimp, tmp)
  
  
  pen   <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen  <- diag(pen)
  b     <- object$pars[, 1]
  
  iter  <- 0
  psej  <- NULL
  
  for (i in 1:length(y))
  {
    if (x[i] > 0) 
    {
      xj    <- x
      xj[i] <- xj[i] - 1
      
      conv <- F
      iter <- 0
      
      while (!conv)
      {
        mu    <- as.vector(exp(D %*% b)) 
        num   <- t(mu / (n - 1) * rimp) * xj
        x_y   <- crossprod(mu / (n - 1), rimp)
        LLnew <- as.vector(log(x_y) %*% xj)
        xcomp <- as.vector(x_y^-1 %*% num)
        
        
        if (iter > 0) {
          delta <- LLnew - LL
          conv  <- LLnew - LL < tol
        }
        
        b  <- NewtRaph(b[-1], D = D, comp = xcomp, n = n - 1, pen = pen, dpen = dpen, mtol = mtol)
        
        LL <- LLnew
        
        iter <- iter + 1
        if (iter == maxit) break
      }
      
      psej <- c(psej, rep(sum(exp(Dpop %*% b)), x[i]))
    }
  }
  
  I    <- (n - 1) * (mean(psej) - psej)
  a    <- sum(I^3) * sum(I^2)^-1.5 / 6
  
  z0   <- qnorm(mean(pseb < pse0))
  u    <- c(.025, .975)
  zu   <- qnorm(u)
  pbca <- pnorm(z0 + (z0 + zu) / (1 - a * (z0 + zu)))
  
  
  pse   <- sum(object$dcmp$fitted)
  ciPM  <- unname(c(quantile(pseb, u)))
  ciBCa <- unname(c(quantile(pseb, pbca)))
  ci    <- matrix(c(pse, ciPM, pse, ciBCa), 2, 3, byrow = T, 
                  dimnames = list(c("Percentile", "BCa"),
                                  c('pse', '2.5%', '97.5%')))
  
  cat("95% Bootstrap Confidence Intervals \n")
  cat("===============================================\n\n")
  print(ci)
  
  invisible(ci)
}
