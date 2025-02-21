
#######################
# bootstrap functions #
#######################

# no missings and no latent variables

bootcomp  <- function(object, lambda, control)
{
  maxit  <- ifelse(exists("maxit", control), control$maxit, 2000)
  tol    <- ifelse(exists("tol", control), control$tol, 1e-5)
  mtol   <- ifelse(exists("mtol", control), control$mtol, 1e-5)
  
  d0  <- object$dobs
  D   <- Matrix(model.matrix(object$formula, d0), sparse = T)
  y   <- c(rmultinom(1, sum(d0$Freq), d0$Freq))
  n   <- sum(y)
  
  
  pen  <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen <- Diagonal(x = pen)
  
  b    <- object$pars[, 1]
  b[1] <- c(-log(sum(exp(D[, -1] %*% b[-1]))) + log(n))
  Dy   <- crossprod(D, y)
  
  conv <- F
  iter <- 1
  
  while (!conv) {
    mu    <- as.vector(exp(D %*% b))
    db    <- Dy - crossprod(D, mu) - pen * b
    I     <- crossprod(D, Diagonal(x = mu)) %*% D + dpen
    H     <- chol2inv(chol(I))
    bt    <- b + H %*% db
    conv  <- sum(abs(b - bt)) < tol
    b     <- bt
    iter  <- iter + 1
  }
  
  b      <- as.vector(b)
  fitted <- as.vector(exp(model.matrix(object$formula, object$dcmp %>% select(-fitted)) %*% b))
  
  ret        <- c(b, fitted)
  names(ret) <- c(colnames(D), 1:length(fitted))
  
  ret
  
}

# missings but no latent variables

bootimp  <- function(object, lambda, control)
{
  maxit  <- ifelse(exists("maxit", control), control$maxit, 2000)
  tol    <- ifelse(exists("tol", control), control$tol, 1e-5)
  mtol   <- ifelse(exists("mtol", control), control$mtol, 1e-5)
  
  f0    <- object$dobs$Freq
  xboot <- c(rmultinom(1, sum(f0), f0))
  dboot <- mutate(object$dobs, Freq = xboot)
  dy    <- na.omit(dboot)
  dz    <- dboot[!complete.cases(dboot), ] %>% filter(Freq > 0)
  y     <- dy$Freq
  z     <- dz$Freq
  x     <- c(y, z)
  n     <- sum(x)
  
  dmy <- dy %>% select(-Freq) %>% mutate(across(!!{object$lists}, ~ as.numeric(.x) - 1)) %>% data.matrix()
  dmz <- dz %>% select(-Freq) %>% mutate(across(!!{object$lists}, ~ as.numeric(.x) - 1)) %>% data.matrix()
  
  tmp  <- matrix(1, nrow(dmy), nrow(dmz))
  for (k in 1:nrow(dmz))
    for (j in 1:ncol(dmy))
      if (!is.na(dmz[k, j]))
        tmp[, k] <- tmp[, k] * (dmy[, j] == dmz[k, j])
  
  rimp <- cbind(Matrix::Diagonal(nrow(tmp)), tmp)
  
  D    <- Matrix(model.matrix(object$formula, dy), sparse = T)
  
  b    <- object$pars[, 1]
  b[1] <- sum(-log(sum(exp(D[, -1] %*% b[-1])))) + log(n)
  
  pen  <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen <- Diagonal(x = pen)
  
  iter  <- 0
  conv  <- FALSE
  
  while (!conv)
  {
    mu    <- as.vector(exp(D %*% b))
    num   <- t(mu / n * rimp) * x 
    x_y   <- crossprod(mu / n, rimp)
    LLnew <- as.vector(log(x_y) %*% x)
    ycomp <- as.vector(x_y^-1 %*% num)
    bt    <- NewtRaph(b[-1], D = D, comp = ycomp, n = n, pen = pen, dpen = dpen, mtol = mtol)
    conv  <- sum(abs(b - bt)) < tol
    b     <- bt
    
    if (iter == maxit) break
    iter <- iter + 1
  }
  
  b      <- as.vector(b)
  fitted <- as.vector(exp(model.matrix(object$formula, 
                                     object$dcmp %>% 
                                       select(-fitted)) %*% b))
  
  ret        <- c(b, fitted)
  names(ret) <- c(colnames(D), 1:length(fitted))
  
  ret
  
}



# missings no latent variables

bootlca  <- function(object, lambda, control)
{
  maxit  <- ifelse(exists("maxit", control), control$maxit, 2000)
  tol    <- ifelse(exists("tol", control), control$tol, 1e-5)
  mtol   <- ifelse(exists("mtol", control), control$mtol, 1e-5)
  
  out <- parse(text = paste(object$lists, "== 0", collapse = " & "))
  dy  <- dplyr::filter(object$dobs, !eval(out))
  y   <- c(rmultinom(1, sum(dy$Freq), dy$Freq))
  dy  <- mutate(dy, Freq = y)
  n   <- sum(y)
  
  nclasses <- prod(object$latlevs)
  
  dx    <- make_dx(dy = dy, nclasses = nclasses, lat = object$lat, 
                   nlat = length(object$lat), latlevs = object$latlevs, 
                   numvars = object$numvars)
  
  D     <- Matrix(model.matrix(object$formula, dx), sparse = T)
  
  xcomp <- numeric(length(y) * nclasses)
  robs  <- matrix(1:length(xcomp), length(y))
  
  b     <- object$pars[, 1]
  
  pen   <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen  <- Diagonal(x = pen)
  
  iter  <- 0
  conv  <- FALSE

  while (!conv)
  {
    mu    <- as.vector(exp(D %*% b))
    ymat  <- matrix(mu, length(y), nclasses)
    yhat  <- rowSums(ymat)
    LLnew <- sum(y %*% log(yhat / n)) - lambda * sum(b[-1]^2)
    
    for (i in 1:length(y))
      xcomp[robs[i, ]]  <- ymat[i, ] * y[i] / yhat[i]
    
    b    <- NewtRaph(b[-1], D = D, comp = xcomp, n = n, pen = pen, dpen = dpen, mtol = mtol)
 
    if (iter > 0)
      conv <- LLnew - LL < tol
    
    LL <- LLnew
    
    iter <- iter + 1
    if (iter == maxit) break
  }
  
  b      <- as.vector(b)
  dpop   <- object$dcmp %>% rename(Freq = fitted)
  fitted <- as.vector(exp(model.matrix(object$formula, dpop) %*% b))
  
  ret        <- c(b, fitted)
  names(ret) <- c(colnames(D), 1:length(fitted))

  ret
  
}

# missings and latent variables

bootimplca  <- function(object, lambda, control)
{
  maxit  <- ifelse(exists("maxit", control), control$maxit, 2000)
  tol    <- ifelse(exists("tol", control), control$tol, 1e-5)
  mtol   <- ifelse(exists("mtol", control), control$mtol, 1e-5)
  
  d0    <- object$dobs
  dboot <- mutate(d0, Freq = c(rmultinom(1, sum(d0$Freq), d0$Freq)))
  dy    <- na.omit(dboot)
  dz    <- dboot[!complete.cases(dboot), ] %>% filter(Freq > 0)
  y     <- dy$Freq
  z     <- dz$Freq
  x     <- c(y, z)
  n     <- sum(x)
  
  rownames(dy) <- 1:nrow(dy)
  rownames(dz) <- 1:nrow(dz)
  
  nclasses <- prod(object$latlevs)
  
  dx <- make_dx(dy = dy, nclasses = nclasses, lat = object$lat, 
                nlat = length(object$lat), latlevs = object$latlevs, 
                numvars = object$numvars)
  
  D  <- Matrix(model.matrix(object$formula, dx), sparse = T)
  
  
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
  
  b    <- object$pars[, 1]
  b[1] <- sum(-log(sum(exp(D[, -1] %*% b[-1])))) + log(n)
  
  pen  <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen <- Diagonal(x = pen)
  
  iter  <- 0
  conv  <- FALSE
  
  while (!conv)
  {
    mu    <- as.vector(exp(D %*% b)) 
    num   <- t(mu / n * rimp) * x
    x_y   <- crossprod(mu / n, rimp)
    LLnew <- as.vector(log(x_y) %*% x)
    xcomp <- as.vector(x_y^-1 %*% num)
    
    
    if (iter > 0) 
      conv  <- LLnew - LL < tol
    
    b  <- NewtRaph(b[-1], D = D, comp = xcomp, n = n, pen = pen, dpen = dpen, mtol = mtol)
    
    LL <- LLnew
    
    iter <- iter + 1
    if (iter == maxit) break
  }
  
  b <- as.vector(b)
  

  fitted     <- as.vector(exp(model.matrix(object$formula, object$dcmp %>% rename(Freq = fitted)) %*% b))
  ret        <- c(b, fitted)
  names(ret) <- c(colnames(D), 1:length(fitted))
  
  ret
  
}
