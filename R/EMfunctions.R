# MSE, no latent variables and no missing values

fcomp <- function(formula, d0, lambda, lists, tol,
                  maxit, mse, imp, lca, printPars, call)
{

  if (mse)
  {
    dpop <- d0
    Dpop <- model.matrix(formula, dpop)
    out  <- parse(text = paste(lists, "== 0", collapse = " & "))
    dpop <- mutate(dpop, Freq = case_when(eval(out) ~ 0,
                                          .default = Freq))
    d0   <- filter(d0, !eval(out))
  } 
   
  D   <- Matrix(model.matrix(formula, d0), sparse = T)
  y   <- d0$Freq
  n   <- sum(y)


  pen  <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen <- Diagonal(x = pen)
  
  b    <- Newt(D = D, x = y, n = n)
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
    if (iter == maxit) break
  }
  
  LL   <- sum(y * log(mu / n)) - sum(pen * b^2)
  
  b <- as.vector(b)

  pars <-  data.frame(b = b, se = sqrt(diag(H)), row.names = colnames(D) ) %>%
    mutate(zval = b / se, pval = 2 * pnorm(-abs(zval)))

  npar <- ncol(D)
  df   <- nrow(d0) - npar
  dev  <- 2 * sum(y * log(y / mu), na.rm = T)
  AIC  <- 2 * (npar - LL)
  BIC  <- npar * log(n) - 2 * LL
  
  if (mse)
  {
    fitted  <- exp(Dpop %*% b)
    dpop    <- mutate(dpop, fitted = fitted) 
    Nhat    <- sum(fitted)
    pse     <- data.frame(n = n, Nhat = Nhat, n0 = Nhat - n, row.names = "")
  } else {
    dpop <- pse <- NULL
  }
  
  dcmp <- if (mse) dpop else mutate(d0, fitted = mu)
    
  out <- output(call = call, mse = mse, imp = imp, lca = lca, lambda = lambda, pse = pse,
                iter = iter, conv = conv, n = n, LL = LL, npar = npar, printPars = printPars, 
                df = df, dev = dev, AIC = AIC, BIC = BIC, pars = pars,
                seed = NA, probs = NULL)
  

  invisible(list(formula = formula,
                 lists   = lists,
                 dobs    = d0,
                 dcmp    = dcmp,
                 pars    = pars))
}

# Imputation of missings: method = 1

imp1 <- function(formula, d0, lambda, stderr, maxit, tol, mtol,
                 seed, lists, mse, lca, imp, printPars, call)
{
  set.seed(seed)
  
  if (mse){
    dpop <- na.omit(d0) 
    Dpop <- model.matrix(formula, dpop)
    out  <- parse(text = paste(lists, "== 0", collapse = " & "))
    dpop <- mutate(dpop, Freq = case_when(eval(out) ~ 0,
                                          .default = Freq))
    d0   <- filter(d0, !eval(out))
  } 
  
  n  <- sum(d0$Freq)
  dy <- na.omit(d0)
  dz <- d0[!complete.cases(d0), ]
  y <- dy$Freq
  z <- dz$Freq
  x <- c(y, z)

  dmy <- dy %>% select(-Freq) %>% mutate(across(!!{lists}, ~ as.numeric(.x) - 1)) %>% data.matrix()
  dmz <- dz %>% select(-Freq) %>% mutate(across(!!{lists}, ~ as.numeric(.x) - 1)) %>% data.matrix()

  tmp  <- matrix(1, nrow(dmy), nrow(dmz))
  for (k in 1:nrow(dmz))
    for (j in 1:ncol(dmy))
      if (!is.na(dmz[k, j]))
        tmp[, k] <- tmp[, k] * (dmy[, j] == dmz[k, j])

  rimp  <- cbind(Matrix::Diagonal(nrow(tmp)), tmp)
  D     <- Matrix(model.matrix(formula, dy), sparse = T)

  ycomp <- y + sum(z) * softmax(runif(length(y) - 1, -.1, .1))
  b     <- Newt(D = D, x = ycomp, n = n)
  
  pen  <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen <- Diagonal(x = pen)

  hist  <- NULL
  iter  <- 0
  conv  <- FALSE

  while (!conv)
  {
    mu    <- as.vector(exp(D %*% b))
    num   <- t(mu / n * rimp) * x 
    x_y   <- crossprod(mu / n, rimp)
    LLnew <- as.vector(log(x_y) %*% x)
    ycomp <- as.vector(x_y^-1 %*% num)
    b     <- NewtRaph(b[-1], D = D, comp = ycomp, n = n, pen = pen, dpen = dpen, mtol = mtol)

    if (iter == 0) {
      hist  <- rbind(hist, cbind(iter = iter, LL = LLnew, stopcrit = - LLnew))
    } else {
      delta <- LLnew - LL
      hist  <- rbind(hist, cbind(iter = iter, LL = LLnew, delta))
      conv  <- LLnew - LL < tol
    }
    LL <- LLnew

    iter <- iter + 1
    if (iter == maxit) break
  }

  b  <- as.vector(b)

  if (stderr[1] == "incomplete") {
    H  <- pracma::hessian(LLimp1, x0 = b[-1], D = D, x = x, rimp = rimp, lambda = lambda)
    se <- c(NA, sqrt(diag(solve(-H))))
  } else {
    se <- sqrt(diag(chol2inv(chol(crossprod(D, Diagonal(x = mu)) %*% D))))
  }
  
  pars  <- data.frame(b = b, row.names = colnames(D)) %>%
    mutate(se    = se,
           zval  = b / se,
           pval  = 2 * pnorm(-abs(b / se)))
  


  yfit  <- mu * sum(y) / n
  dy    <- mutate(dy, yfitted = yfit, yzcomp = ycomp, yzfitted = mu) %>% rename(y = Freq)
  dp    <- devz(dmz = dmz, z = z, cpz = x_y[-(1:length(y))])
  dev   <- 2 * sum(y * log(y / yfit), na.rm = T) + dp[[1]]
  npar  <- ncol(D)  + dp[[2]]
  df    <- nrow(d0) - npar
  AIC   <- 2 * (npar - 1 - LL)
  BIC   <- (npar - 1) * log(n) - 2 * LL

  if (mse)
  {
    fitted  <- exp(Dpop %*% b)
    dpop    <- mutate(dpop, fitted = fitted) 
    Nhat    <- sum(fitted)
    pse     <- data.frame(n = n, Nhat = Nhat, n0 = Nhat - n, row.names = "")
  } else {
    pse <- dpop <- NULL
  }
  

  out <- output(call = call, mse = mse, imp = imp, lca = lca, lambda = lambda, pse = pse,
                iter = iter, conv = conv, n = n, LL = LL, npar = npar, df = df,
                dev = dev, AIC = AIC, BIC = BIC, dpop = dpop, pars = pars, printPars = printPars,
                seed = seed, probs = NULL)

  invisible(list(formula = formula,
                 lists   = lists,
                 dobs    = d0,
                 dcmp    = if (mse) dpop else dy,
                 pars    = pars,
                 hist    = data.frame(hist)))

}

# Imputation of missings: method = 2

imp2 <- function(formula, d0, lambda, stderr, maxit, tol, mtol,
                 seed, lists, mse, lca, imp, printPars, call)
{
  set.seed(seed)
  
  if (mse){
    dpop <- na.omit(d0) 
    Dpop <- model.matrix(formula, dpop)
    out  <- parse(text = paste(lists, "== 0", collapse = " & "))
    dpop <- mutate(dpop, Freq = case_when(eval(out) ~ 0,
                                  .default = Freq))
    d0   <- filter(d0, !eval(out))
  } 
  
  n  <- sum(d0$Freq)
  dy <- na.omit(d0)
  dz <- d0[!complete.cases(d0), ]
  dx <- rbind(dy, dz)
  y <- dy$Freq
  z <- dz$Freq
  x <- c(y, z)

  rownames(dy) <- 1:nrow(dy)
  rownames(dz) <- 1:nrow(dz)
  
  dmy <- dy %>% select(-Freq) %>% 
    mutate(across(everything(), ~ as.numeric(.x)),
           across(!!{lists}, ~.x - 1)) 
  dmz <- dz %>% select(-Freq) %>% 
    mutate(across(everything(), ~ as.numeric(.x)),
           across(!!{lists}, ~.x - 1)) 
  
  rimp <- vector("list", nrow(dmz)) 
  grps <- vector("list", nrow(dmz))
  
  for (i in 1:nrow(dmz))
  {
    obsval    <- dmz[i, !is.na(dmz[i, ]), drop = F]
    grps[[i]] <- colnames(obsval)
    ex        <- paste(colnames(obsval), unlist(obsval), sep = " == ", collapse = " & ")
    rimp[[i]] <- as.numeric(rownames(subset(dy, subset = eval(parse(text = ex)))))
  }
  
  D     <- Matrix(model.matrix(formula, data = dy), sparse = T)
  
  ycomp <- y + sum(z) * softmax(runif(length(y) - 1, -.1, .1))
  b     <- Newt(D = D, x = ycomp, n = n)
  
  cpz  <- numeric(length(rimp))
  zhat <- numeric(nrow(dy))
  
  pen  <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen <- Diagonal(x = pen)
  
  hist  <- NULL
  iter  <- 0
  conv  <- FALSE
  delta <- 1

  while (!conv)
  {
    mu    <- as.vector(exp(D %*% b))
    zhat <- 0 * zhat
    for (i in 1:length(z))
      zhat[rimp[[i]]] <-  zhat[rimp[[i]]] + z[i] * mu[rimp[[i]]] / sum(mu[rimp[[i]]])
    
    ycomp <- y + zhat
    
    for (i in 1:length(rimp))
      cpz[i] <- sum(mu[rimp[[i]]]) / n
    
 
    LLnew <- sum(x * log(c(mu / n, cpz)))
    
    
    if (iter == 0) {
      hist  <- rbind(hist, cbind(iter = iter, LL = LLnew, stopcrit = -LLnew))
    } else {
      delta <- LLnew - LL
      hist  <- rbind(hist, cbind(iter, LLnew, delta))
      conv  <- delta < tol
    }

   b  <- NewtRaph(b[-1], D = D, comp = ycomp, n = n, pen = pen, dpen = dpen, mtol = mtol)

   LL <- LLnew

    if (iter == maxit) break
    iter  <- iter + 1
  }
  
  b   <- as.vector(b)

  if (stderr[1] == "incomplete") {
    H  <- pracma::hessian(LLimp2, x0 = b[-1], D = D, z = z, y = y, rimp = rimp, lambda = lambda)
    se <- c(NA, sqrt(diag(solve(-H))))
  } else {
    se <- sqrt(diag(chol2inv(chol(crossprod(D, Diagonal(x = mu * n)) %*% D))))
  }
  pars  <- data.frame(b = b, row.names = colnames(D)) %>%
    mutate(se    = se,
           zval  = b / se,
           pval  = 2 * pnorm(-abs(b / se)))
  
  yhat  <- mu * sum(y) / n
  dy    <- mutate(dy, yfitted = yhat, yzcomp = ycomp, yzfitted = mu) %>% rename(y = Freq)
  dp    <- devz(dmz = dmz, z = z, cpz = cpz)
  dev   <- 2 * sum(y * log(y / yhat), na.rm = T) + dp[[1]]
  npar  <- ncol(D) + dp[[2]]
  df    <- nrow(d0) - npar 
  AIC   <- 2 * (npar - 1 - LL)
  BIC   <- (npar - 1) * log(n) - 2 * LL
  
  if (mse)
  {
    fitted  <- exp(Dpop %*% b)
    dpop    <- mutate(dpop, fitted = fitted) 
    Nhat    <- sum(fitted)
    pse     <- data.frame(n = n, Nhat = Nhat, n0 = Nhat - n, row.names = "")
  } else {
    dpop <- pse <- NULL
  }
  
  out <- output(call = call, mse = mse, imp = imp, lca = lca, lambda = lambda, pse = pse,
                iter = iter, conv = conv, n = n, LL = LL, npar = npar, df = df,
                dev = dev, AIC = AIC, BIC = BIC, dpop = dpop, pars = pars, printPars = printPars,
                seed = seed, probs = NULL)
  
  invisible(list(formula = formula,
                 lists   = lists,
                 dobs    = d0,
                 dcmp    = if (mse) dpop else dy,
                 pars    = pars,
                 hist    = hist))
  
}

# Latent class analysis without missings

lca2 <- function(formula, d0, latlevs, lambda, stderr, maxit, tol, mtol,
                 seed, nlat, lat, lists, numvars, covs, vlat, mse, lca, imp,
                 printPars, call)
{
  
  set.seed(seed)
  
  nclasses <- prod(latlevs)
  prlat    <- softmax(-runif(1, vlat - .5, vlat + .5) * (2:nclasses))
  
  if (mse)
  {
    out   <- parse(text = paste(lists, "== 0", collapse = " & "))
    dpop <- mutate(d0, Freq = case_when(eval(out) ~ 0,
                                        .default = Freq))
    dy    <- filter(d0, !eval(out))
    y     <- dy$Freq
    n     <- sum(y)
  } else {
    dpop <- NULL
    dy   <- d0
    y    <- dy$Freq
    n    <- sum(y)
  }
  
  dx    <- make_dx(dy = dy, nclasses = nclasses, lat = lat, nlat = nlat,
                  latlevs = latlevs, numvars = numvars)
  D     <- Matrix(model.matrix(formula, data = dx), sparse = T)
  xcomp <- rep(y, nclasses) * rep(prlat,  each = nrow(dy))
  b     <- Newt(D = D, x = xcomp, n = n)
  robs  <- matrix(1:length(xcomp), length(y))
  
  pen   <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen  <- Diagonal(x = pen)
  
  hist  <- NULL
  iter  <- 0
  conv  <- FALSE
  delta <- 1

  while (!conv)
  {
    mu    <- as.vector(exp(D %*% b))
    ymat  <- matrix(mu, length(y), nclasses)
    yhat  <- rowSums(ymat)
    LLnew <- sum(y %*% log(yhat / n)) - lambda * sum(b[-1]^2)
    
    for (i in 1:length(y))
      xcomp[robs[i, ]]  <- ymat[i, ] * y[i] / yhat[i]

    b <- NewtRaph(b[-1], D = D, comp = xcomp, n = n, pen = pen, dpen = dpen, mtol = mtol)
    
    if (iter == 0) {
      hist  <- rbind(hist, cbind(iter = iter, LL = LLnew, stopcrit = - LLnew))
    } else {
      delta <- LLnew - LL
      hist  <- rbind(hist, cbind(iter = iter, LL = LLnew, delta))
      conv  <- delta < tol
    }
    LL <- LLnew
    
    iter <- iter + 1
    if (iter == maxit) break
  }
  
  b <- as.vector(b)

  if (stderr[1] == "incomplete") {
    H  <- pracma::hessian(LLlca, b[-1], D = D, y = y, nclasses = nclasses, robs = robs, lambda = lambda)
    se <- c(NA, sqrt(diag(solve(-H))))
  } else {
    se <- sqrt(diag(chol2inv(chol(crossprod(D, Diagonal(x = mu)) %*% D))))
  }
  
  pars  <- data.frame(b = b, row.names = colnames(D)) %>%
    mutate(se    = se,
           zval  = b / se,
           pval  = 2 * pnorm(-abs(b / se)))
  
  
  dx <- mutate(dx, fitted = mu) %>% select(-Freq) 
  dy <- mutate(dy, yfitted = yhat) 
  
  
  dev  <- 2 * sum(xcomp * log(xcomp / mu ), na.rm = T)
  npar <- ncol(D)
  df   <- nrow(dy) - npar
  AIC  <- 2 * (npar - LL)
  BIC  <- npar * log(n) - 2 * LL
  
  if (mse) {
    dxpop  <- make_dxpop(dpop = dpop, nclasses = nclasses, lat = lat, nlat = nlat,
                         latlevs = latlevs, numvars = numvars)
    fitted <- c(exp(model.matrix(formula, dxpop) %*% b))
    dxpop  <- mutate(dxpop, fitted = fitted) %>% select(-Freq)
    Nhat   <- sum(fitted)
    pse    <- data.frame(n = n, Nhat = Nhat, n0 = Nhat - n, row.names = "")
    
  }  else {
    
    pse <- dpop <- NULL
    
  }
  
  dcmp <- if (mse) dxpop else dx
  
  probs <- get_probs(dx = dcmp, covs = covs, lists = lists, lat = lat)
  
  out <- output(call = call, mse = mse, lca = lca, imp = imp, lambda = lambda,
                iter = iter, conv = conv, n = n, LL = LL, npar = npar, df = df,
                dev = dev, AIC = AIC, BIC = BIC, dpop = dpop, pars = pars, 
                printPars = printPars, pse = pse, seed = seed, probs = probs)
  
  invisible(list(formula = formula,
                 numvars = numvars,
                 lists   = lists,
                 lat     = lat,
                 latlevs = latlevs,
                 dobs    = dy,
                 dcmp    = if (mse) dxpop else dx,
                 pars    = pars,
                 probs   = probs,
                 hist    = hist))
}

# Latent class analysis with imputation of missings: method = 1

imp1lca <- function(formula, d0, latlevs, lambda, stderr, maxit, tol, mtol,
                    seed, nlat, lat, lists, numvars, covs, vlat, mse, lca, imp,
                    printPars, call)
{
  set.seed(seed)
  
  nclasses <- prod(latlevs)
  prlat    <- softmax(-runif(1, vlat - .5, vlat + .5) * (2:nclasses))
  
  if (mse){
    dpop <- na.omit(d0) 
    out  <- parse(text = paste(lists, "== 0", collapse = " & "))
    dpop <- mutate(dpop, Freq = case_when(eval(out) ~ 0,
                                          .default = Freq))
    d0   <- filter(d0, !eval(out))
  } 

  n  <- sum(d0$Freq)
  dy <- na.omit(d0)
  dz <- d0[!complete.cases(d0), ]
  y  <- dy$Freq
  z  <- dz$Freq
  x  <- c(y, z)
  
  
  rownames(dy) <- 1:nrow(dy)
  rownames(dz) <- 1:nrow(dz)
  
  dx <- make_dx(dy = dy, nclasses = nclasses, lat = lat, nlat = nlat,
                latlevs = latlevs, numvars = numvars)
  D  <- Matrix(model.matrix(formula, dx), sparse = T)
  

  dmy <- dy %>% select(-Freq) %>% mutate(across(!!{lists}, ~ as.numeric(.x) - 1)) %>% data.matrix()
  dmz <- dz %>% select(-Freq) %>% mutate(across(!!{lists}, ~ as.numeric(.x) - 1)) %>% data.matrix()
  
  tmp <- matrix(1, nrow(dmy), nrow(dmz))
  
  for (k in 1:nrow(dmz))
    for (j in 1:ncol(dmy))
      if (!is.na(dmz[k, j]))
        tmp[, k] <- tmp[, k] * (dmy[, j] == dmz[k, j])
  
  tmp <- rimp <- cbind(Diagonal(nrow(tmp)), tmp)
  
  for (j in 2:nclasses)
    rimp <- rbind(rimp, tmp)

  yhat  <- rep(y, nclasses) * rep(prlat,  each = nrow(dy))  
#  zhat  <- sum(z) * softmax(runif(nclasses * length(y) - 1, -.1, 1))
  zhat  <- sum(z) * rep(prlat, each = length(y)) / (nclasses * length(y))
  xcomp <- yhat + zhat
  
  b     <- Newt(D = D, x = xcomp, n = n)

  pen  <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen <- Diagonal(x = pen)
  
  hist  <- NULL
  iter  <- 0
  conv  <- FALSE
  
  while (!conv)
  {
    mu    <- as.vector(exp(D %*% b)) 
    num   <- t(mu / n * rimp) * x
    x_y   <- crossprod(mu / n, rimp)
    LLnew <- as.vector(log(x_y) %*% x)
    xcomp <- as.vector(x_y^-1 %*% num)


    if (iter == 0) {
      hist  <- rbind(hist, cbind(iter = iter, LL = LLnew, stopcrit = - LLnew))
    } else {
      delta <- LLnew - LL
      hist  <- rbind(hist, cbind(iter = iter, LL = LLnew, delta))
      conv  <- LLnew - LL < tol
    }
    
    b  <- NewtRaph(b[-1], D = D, comp = xcomp, n = n, pen = pen, dpen = dpen, mtol = mtol)

    LL <- LLnew
    
    iter <- iter + 1
    if (iter == maxit) break
  }
  
  b <- as.vector(b)
  
  if (stderr[1] == "incomplete") {
    H  <- pracma::hessian(LLimp1lca, x0 = b[-1], D = D, x = x, n = n, rimp = rimp, lambda = lambda)
    se <- c(NA, sqrt(diag(solve(-H))))
  } else {
    se <- sqrt(diag(chol2inv(chol(crossprod(D, Diagonal(x = mu * n)) %*% D))))
  }
  
  pars  <- data.frame(b = as.vector(b), row.names = colnames(D)) %>%
    mutate(se    = se,
           zval  = b / se,
           pval  = 2 * pnorm(-abs(b / se)))
  
  
  yfit  <- rowSums(matrix(mu, ncol = nclasses) * sum(y) / n) 
  yzcmp <- rowSums(matrix(xcomp, ncol = nclasses))
  yzfit <- rowSums(matrix(mu, ncol = nclasses))
  dy    <- mutate(dy, yfitted = yfit, yzcomp = yzcmp, yzfitted = yzfit) %>% rename(y = Freq)
  dp    <- devz(dmz = dz %>% select(-Freq), z = z, cpz = x_y[-(1:length(y))])
  dev   <- 2 * sum(y * log(y / yfit), na.rm = T) + dp[[1]]
  npar  <- ncol(D) + dp[[2]]
  df    <- nrow(D) + nrow(dz) - npar
  AIC   <- 2 * (npar - 1 - LL)
  BIC   <- (npar - 1) * log(n) - 2 * LL
  
  if (mse)
  {
    dxpop   <- make_dxpop(dpop = dpop, nclasses = nclasses, lat = lat, nlat = nlat, latlevs = latlevs, numvars = numvars)
    Dpop    <- model.matrix(formula, dxpop)
    fitted  <- as.vector(exp(Dpop %*% b))
    dxpop   <- mutate(dxpop, fitted = fitted) %>% select(-Freq)
    Nhat    <- sum(fitted)
    pse     <- data.frame(n = n, Nhat = Nhat, n0 = Nhat - n, row.names = "")
  } else {
    pse <- dpop <- NULL
  }
  
  dcmp  <- if (mse) dxpop else mutate(dx, fitted = mu) %>% select(-Freq)
  probs <- get_probs(dx = dcmp, covs = covs, lists = lists, lat = lat)
  

  out <- output(call = call, mse = mse, imp = imp, lca = lca, lambda = lambda, pse = pse,
                iter = iter, conv = conv, n = n, LL = LL, npar = npar, df = df,
                dev = dev, AIC = AIC, BIC = BIC, dpop = dpop, pars = pars,
                printPars = printPars, seed = seed, probs = probs)
  
  invisible(list(formula = formula,
                 numvars = numvars,
                 lists   = lists,
                 lat     = lat,
                 latlevs = latlevs,
                 dobs    = d0,
                 dcmp    = if (mse) dxpop else dx,
                 pars    = pars,
                 hist    = data.frame(hist)))
  
}

# Latent class analysis with imputation of missings: method = 2

imp2lca <- function(formula, d0, latlevs, lambda, stderr, maxit, tol, mtol,
                 seed, nlat, lat, lists, numvars, covs, vlat, mse, lca, imp,
                 printPars, call)
{
  set.seed(seed)
  
  nclasses <- prod(latlevs)
  prlat    <- softmax(-runif(1, vlat - .5, vlat + .5) * (2:nclasses) / nclasses)
  
  
  if (mse){
    dpop <- na.omit(d0) 
    out  <- parse(text = paste(lists, "== 0", collapse = " & "))
    dpop <- mutate(dpop, Freq = case_when(eval(out) ~ 0,
                                          .default = Freq))
    d0   <- filter(d0, !eval(out))
  } 
  
  dy <- na.omit(d0)
  dz <- d0[!complete.cases(d0), ]
  dx <- make_dx(dy = dy, nclasses = nclasses, lat = lat, nlat = nlat,
                latlevs = latlevs, numvars = numvars)
  
  y  <- dy$Freq
  z  <- dz$Freq
  x  <- c(y, z)
  n  <- sum(d0$Freq)
  
  rownames(dy) <- 1:nrow(dy)
  rownames(dz) <- 1:nrow(dz)
  
  dmy <- dy %>% select(-Freq) %>% 
    mutate(across(everything(), ~ as.numeric(.x)),
           across(!!{lists}, ~.x - 1)) 
  dmz <- dz %>% select(-Freq) %>% 
    mutate(across(everything(), ~ as.numeric(.x)),
           across(!!{lists}, ~.x - 1)) 
  
  rimp <- vector("list", nrow(dmz)) 
  grps <- vector("list", nrow(dmz))
  
  for (i in 1:nrow(dmz))
  {
    obsval    <- dmz[i, !is.na(dmz[i, ]), drop = F]
    grps[[i]] <- colnames(obsval)
    ex        <- paste(colnames(obsval), unlist(obsval), sep = " == ", collapse = " & ")
    tmp       <- as.numeric(rownames(subset(dmy, subset = eval(parse(text = ex)))))
    rows      <- NULL
    for (j in 1:nclasses - 1)
      rows <- c(rows, tmp + j * nrow(dy))
    
    rimp[[i]] <- rows
  }

  D    <- Matrix(model.matrix(formula, dx), sparse = T)
  yhat <- rep(y, nclasses) * rep(prlat,  each = nrow(dy))
#    zhat  <- sum(z) * softmax(runif(nclasses * length(y) - 1, -.1, 1))
  zhat  <- sum(z) * rep(prlat, each = length(y)) / (nclasses * length(y))
  xhat <- yhat + zhat
  b    <- Newt(D = D, x = xhat, n = n)

  robs <- matrix(1:length(yhat), length(y))
  pen  <- lambda * c(0, rep(2, ncol(D) - 1))
  dpen <- Diagonal(x = lambda * c(0, rep(2, ncol(D) - 1)))
  
  hist  <- NULL
  iter  <- 0
  conv  <- FALSE
  delta <- 1
  
  while (!conv)
  {
    mu    <- as.vector(exp(D %*% b))
    pcomp <- mu / sum(mu)

    zhat  <- 0 * zhat
    for (i in 1:length(z))
      zhat[rimp[[i]]] <-  zhat[rimp[[i]]] + pcomp[rimp[[i]]] * z[i] / sum(pcomp[rimp[[i]]])
    
    mprobs <- matrix(pcomp, length(y), nclasses)
    rsum   <- rowSums(mprobs)
    
    for (i in 1:length(y))
      yhat[robs[i, ]]  <- c(mprobs[i, ]) * y[i] / rsum[i]
    
    xcomp <- yhat + zhat
    
    prmiss <- NULL
    for (i in 1:length(rimp))
      prmiss <- c(prmiss, sum(mprobs[rimp[[i]]]))
    
    LLnew <- sum(c(y, z) * log(c(rsum, prmiss))) - lambda * sum(b[-1]^2)

    b  <- NewtRaph(b[-1], D = D, comp = xcomp, n = n, pen = pen, dpen = dpen, mtol = mtol)
    
    if (iter == 0) {
      hist  <- c(iter = iter, LL = LLnew, delta = -LLnew)
    } else {
      delta <- LLnew - LL
      hist  <- rbind(hist, cbind(iter, LLnew, delta))
    }
    
    conv <- delta < tol
    
    LL <- LLnew
    
    if (iter == maxit) break
    iter  <- iter + 1
  }
  
  b <- as.vector(b)

  if (stderr[1] == "incomplete") {
    H  <- pracma::hessian(LLimp2lca, x0 = b[-1], D = D, y = y, z = z, nclasses = nclasses, rimp = rimp, lambda = lambda)
    se <- c(NA, sqrt(diag(solve(-H))))
  } else {
    se <- sqrt(diag(chol2inv(chol(crossprod(D, Diagonal(x = mu)) %*% D))))
  }
  
  pars  <- data.frame(b = b, row.names = colnames(D)) %>%
    mutate(se    = se,
           zval  = b / se,
           pval  = 2 * pnorm(-abs(b / se)))
  
  cpz <- NULL
  for (i in 1:length(rimp))
    cpz <- c(cpz, sum(pcomp[rimp[[i]]])) 
  
  

  yfit  <- rowSums(matrix(mu, ncol = nclasses) * sum(y) / n) 
  yzcmp <- rowSums(matrix(xcomp, ncol = nclasses)) 
  yzfit <- rowSums(matrix(xcomp, ncol = nclasses)) 
  dy    <- mutate(dy, yfitted = yfit, yzcomp = yzcmp, yzfitted = yzfit) %>% rename(y = Freq)
  dp    <- devz(dmz = dmz, z = z, cpz = cpz)
  dev   <- 2 * sum(y * log(y / yfit), na.rm = T) + dp[[1]]
  npar  <- ncol(D) + dp[[2]]
  df    <- nrow(D) + nrow(dz) - npar 
  AIC   <- 2 * (npar - 1 - LL)
  BIC   <- (npar - 1) * log(n) - 2 * LL
  
  if (mse)
  {
    dxpop   <- make_dxpop(dpop = dpop, nclasses = nclasses, lat = lat, nlat = nlat, latlevs = latlevs, numvars = numvars)
    Dpop    <- model.matrix(formula, dxpop)
    fitted  <- as.vector(exp(Dpop %*% b))
    dxpop   <- mutate(dxpop, fitted = fitted) %>% select(-Freq)
    Nhat    <- sum(fitted)
    pse     <- data.frame(n = n, Nhat = Nhat, n0 = Nhat - n, row.names = "")
  } else {
    pse <- dpop <- NULL
  }
  
  dcmp  <- if (mse) dxpop else mutate(dx, fitted = mu) %>% select(-Freq)
  probs <- get_probs(dx = dcmp, covs = covs, lists = lists, lat = lat)
  

  out <- output(call = call, mse = mse, imp = imp, lca = lca, lambda = lambda, pse = pse,
                iter = iter, conv = conv, n = n, LL = LL, npar = npar, df = df,
                dev = dev, AIC = AIC, BIC = BIC, dpop = dpop, pars = pars,
                printPars = printPars, seed = seed, probs = probs)
  
  invisible(list(formula = formula,
                 numvars = numvars,
                 lists   = lists,
                 lat     = lat,
                 latlevs = latlevs,
                 dobs    = d0,
                 dcmp    = if (mse) dxpop else dx,
                 pars    = pars,
                 hist    = hist))
  
}

