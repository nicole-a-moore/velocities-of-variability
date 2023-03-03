## Functions needed from fractal package 
## Author: William Constantine and Donald Percival (Applied Physics Laboratory, University of Washington)

FDSimulate <- function (delta, innovations.var = 1, method = "ce", seed = 0) 
{
  "FDCirculantEmbedding" <- function(delta, innovations.var = 1, 
                                     n.sample = 16, seed = 0) {
    checkScalarType(delta, "numeric")
    checkScalarType(innovations.var, "numeric")
    checkScalarType(n.sample, "integer")
    checkScalarType(seed, "integer")
    seed <- as.integer(abs(seed))
    if (n.sample < 2) 
      stop("n.sample must be an integer greater than unity")
    ncumsum <- max((2 * delta - 1)%/%2 + 1, 0)
    delta <- delta - ncumsum
    fdp.acvs <- lmACF(lmModel("fdp", delta = delta, innov = innovations.var), 
                      lag.max = n.sample, type = "covariance")@data
    fdp.acvs <- c(fdp.acvs, rev(fdp.acvs[2:n.sample]))
    S <- Re(fft(fdp.acvs, inverse = FALSE))
    z <- as.vector(.Call("RS_fractal_bootstrap_circulant_embedding", 
                         S, seed, COPY = rep(FALSE, 2), CLASSES = c("matrix", 
                                                                    "integer"), PACKAGE = "ifultools"))[seq(n.sample)]
    if (!ncumsum) 
      return(z)
    for (i in seq(ncumsum)) z <- cumsum(z)
    z
  }
  checkVectorType(delta, "numeric")
  checkScalarType(method, "character")
  if (!is.numeric(innovations.var)) 
    stop("innovations.var must be a numeric object")
  checkScalarType(seed, "integer")
  seed <- as.integer(abs(seed))
  n.innov <- length(innovations.var)
  n.delta <- length(delta)
  if (n.delta > n.innov) 
    innovations.var <- c(innovations.var, rep(innovations.var[n.innov], 
                                              n.delta - n.innov))
  method <- match.arg(method, c("ce", "cholesky"))
  if (method == "cholesky") {
    z <- as.vector(.Call("RS_fractal_fdp_simulate", delta, 
                         innovations.var, COPY = rep(FALSE, 2), CLASSES = c("numeric", 
                                                                            "numeric"), PACKAGE = "ifultools"))
  }
  else {
    delta.levels <- sort(unique(delta))
    idelta <- match(delta, delta.levels)
    z <- apply(matrix(delta.levels), MARGIN = 1, function(x, 
                                                          n.delta, seed, FDCirculantEmbedding) FDCirculantEmbedding(delta = x, 
                                                                                                                    innovations.var = 1, n.sample = n.delta), n.delta = n.delta, 
               seed = seed, FDCirculantEmbedding = FDCirculantEmbedding)
    z <- z[cbind(seq(n.delta), idelta)] * innovations.var
  }
  oldClass(z) <- "FDSimulate"
  attr(z, "delta") <- delta
  attr(z, "innov") <- innovations.var
  attr(z, "method") <- "circulant embedding"
  z
}

lmModel <- function (model, variance. = 1, delta = 0.45, alpha = -0.9, 
                     HG = 0.95, HB = 0.95, innovations.var = NULL, Cs = NULL, 
                     bterms = 10, dterms = 10, M = 100) 
{
  fdp.var.to.ivar <- function(var, delta) {
    del.s <- if (delta < 0.5) 
      delta
    else delta - floor(delta + 0.5)
    var * (gamma(1 - del.s))^2/gamma(1 - 2 * del.s)
  }
  fdp.ivar.to.var <- function(ivar, delta) {
    del.s <- if (delta < 0.5) 
      delta
    else delta - floor(delta + 0.5)
    ivar * gamma(1 - 2 * del.s)/(gamma(1 - del.s))^2
  }
  ppl.Cs.to.var <- function(Cs, alpha) {
    ppl.term.m1.to.m3 <- function(n, alpha) {
      tn <- 2 * n
      ((-1)^n) * (pi^tn)/(factorial(tn - 1) * (alpha + 
                                                 1 + tn))
    }
    ppl.term.m3.to.m5 <- function(n) {
      tn <- 2 * n
      ((-1)^(n + 1)) * (pi^(tn + 2)) * (1 - 2^tn)/(factorial(tn - 
                                                               1) * (alpha + 3 + tn))
    }
    if (alpha > -1) 
      Cs/((2^alpha) * (alpha + 1))
    else {
      if (alpha == -1) 
        Cs * 6.59311055481803
      else {
        if (alpha > -3) {
          Cs * (1 + 0.5 * (sum(ppl.term.m1.to.m3(seq(1, 
                                                     15, 2), alpha)) + sum(ppl.term.m1.to.m3(seq(2, 
                                                                                                 16, 2), alpha))))/((2^(alpha - 2)) * (alpha + 
                                                                                                                                         1))
        }
        else {
          if (alpha == -3) 
            Cs * 185.306445416753
          else {
            if (alpha > -5) {
              Cs * (1 - (pi^2/((alpha + 2) * (alpha + 
                                                3))) + (sum(ppl.term.m3.to.m5(seq(1, 
                                                                                  15, 2))) + sum(ppl.term.m3.to.m5(seq(2, 
                                                                                                                       16, 2))))/(2 * (alpha + 2) * (alpha + 
                                                                                                                                                       3)))/((2^(alpha - 4)) * (alpha + 1))
            }
            else {
              if (alpha == -5) 
                Cs * 5508.61024271439
              else stop("sorry - cannot currently handle alpha < -5")
            }
          }
        }
      }
    }
  }
  checkScalarType(model, "character")
  checkScalarType(variance., "numeric")
  if (variance. <= 0) 
    stop("variance must be positive")
  checkScalarType(delta, "numeric")
  checkScalarType(alpha, "numeric")
  checkScalarType(HG, "numeric")
  checkScalarType(HB, "numeric")
  checkScalarType(bterms, "integer")
  checkScalarType(dterms, "integer")
  checkScalarType(M, "integer")
  model <- match.arg(model, c("ppl", "fdp", "fgn", "dfbm"))
  if (model == "fgn" && (HG >= 1 || HG <= 0)) 
    stop("HG must be on the interval (0,1)")
  if (model == "dfbm" && (HB >= 1 || HB <= 0)) 
    stop("HB must be on the interval (0,1)")
  if (model == "ppl") {
    if (!is.null(Cs)) {
      checkScalarType(Cs, "numeric")
      if (Cs <= 0) 
        stop("Cs must be positive")
      variance. <- ppl.Cs.to.var(Cs, alpha)
      innovations.var <- Cs/(2 * exp(1))^alpha
    }
    else {
      if (!is.null(innovations.var)) {
        checkScalarType(, "numeric")
        if (innovations.var <= 0) 
          stop("innovations.var must be positive")
        Cs <- ((2 * exp(1))^alpha) * innovations.var
        variance. <- ppl.Cs.to.var(Cs, alpha)
      }
      else {
        Cs <- variance./ppl.Cs.to.var(1, alpha)
        innovations.var <- Cs/(2 * exp(1))^alpha
      }
    }
  }
  if ((model == "fdp") && !is.null(innovations.var)) {
    checkScalarType(innovations.var, "numeric")
    if (innovations.var <= 0) 
      stop("innovations.var must be positive")
    variance. <- fdp.ivar.to.var(innovations.var, delta)
  }
  else {
    innovations.var <- fdp.var.to.ivar(variance., delta)
  }
  z <- c(list(model = model), switch(model, ppl = list(alpha = alpha, 
                                                       Cs = Cs, bterms = bterms, dterms = dterms, innovations.var = innovations.var, 
                                                       variance = variance.), fdp = list(delta = delta, innovations.var = innovations.var, 
                                                                                         variance = variance.), fgn = list(HG = HG, M = M, variance = variance.), 
                                     dfbm = list(HB = HB, M = M, variance = variance.)))
  oldClass(z) <- "lmModel"
  z
}

library(ifultools)
lmACF <- function (x, lag.max = 32, type = "correlation") 
{
  if (!is(x, "lmModel")) 
    stop("x must be an object of class \"lmModel\"")
  if (x$model == "dfbm") 
    stop("FBM model not supported in ACF estimator")
  checkScalarType(type, "character")
  checkScalarType(lag.max, "integer")
  type <- match.arg(lowerCase(type), c("correlation", "covariance", 
                                       "partial"))
  partial <- is.element(type, "partial")
  if (x$model == "ppl") {
    "Bm" <- function(m, alpha) {
      q <- 4 * m + 1
      (pi^(alpha + q) * ((q + 1) * q * (alpha + q + 2) - 
                           (alpha + q) * pi^2))/((alpha + q) * (alpha + 
                                                                  q + 2) * gamma(q + 2))
    }
    "DeltaTau" <- function(tau, alpha, dterms) {
      x0 <- pi * tau + pi/2
      result <- 0
      an <- alpha * x0^(alpha - 1)
      An <- 2
      for (n in seq(from = 1, by = 2, length = dterms)) {
        result <- result + an * An
        An <- 2 * (n + 2) * (pi/2)^(n + 1) - (n + 2) * 
          (n + 1) * An
        an <- (an * (alpha - n) * (alpha - n - 1))/((n + 
                                                       2) * (n + 1) * x0^2)
      }
      result * (-1)^(tau + 1)
    }
    temp <- 0:lag.max
    temp[1] <- x$variance
    if (lag.max > 0) {
      IandDeltaTau <- 1:lag.max
      IandDeltaTau[1] <- sum(sapply(x$bterms:0, Bm, x$alpha))
      if (lag.max > 1) {
        lags <- 2:lag.max
        IandDeltaTau[lags] <- sapply(lags - 1, DeltaTau, 
                                     x$alpha, x$dterms)
      }
      temp[2:(lag.max + 1)] <- (2 * x$Cs * cumsum(IandDeltaTau))/(2 * 
                                                                    pi * (1:lag.max))^(x$alpha + 1)
    }
    exponent <- x$alpha
  }
  else if (x$model == "fdp") {
    temp <- seq(0, lag.max)
    temp[1] <- x$variance
    if (lag.max) {
      positive.lags <- seq(lag.max)
      temp[2:(lag.max + 1)] <- (positive.lags + x$delta - 
                                  1)/(positive.lags - x$delta)
      temp <- cumprod(temp)
    }
    exponent <- x$delta
  }
  else if (x$model == "fgn") {
    temp <- 0:lag.max
    temp <- (x$variance * (abs(temp + 1)^(2 * x$HG) - 2 * 
                             abs(temp)^(2 * x$HG) + abs(temp - 1)^(2 * x$HG)))/2
    exponent <- x$HG
  }
  if (type == "partial") {
    data <- ACVStoPACS(temp)
    tag <- "PACS"
  }
  else if (type == "correlation") {
    data <- temp/temp[1]
    tag <- "ACF"
  }
  else {
    data <- temp
    tag <- "ACVF"
  }
  create.signalSeries(data, position = list(from = ifelse1(partial, 
                                                           1, 0), by = 1, length = length(data), units = "lag"), 
                      title.data = paste(upperCase(x$model), "(", exponent, 
                                         ") ", tag, sep = ""), documentation = paste(tag, 
                                                                                     " for ", upperCase(x$model), "(", exponent, ") process", 
                                                                                     sep = ""), units = "ACF")
}


## functions needed from ifultools package
checkScalarType <- function (x, isType = "numeric") 
{
  if (!is.character(isType)) 
    stop("isType must be an object of class character")
  if (isType == "integer") {
    if (!is.numeric(x) || length(x) > 1) 
      stop(deparseText(substitute(x)), " must be scalar of class ", 
           isType)
  }
  else {
    if (!eval(parse(text = paste("is.", isType, "(x)", sep = ""))) || 
        length(x) > 1) 
      stop(deparseText(substitute(x)), " must be scalar of class ", 
           isType)
  }
  invisible(NULL)
}

checkVectorType <- function (x, isType = "numeric") 
{
  checkScalarType(isType, "character")
  if (isType == "integer") {
    if (!isVectorAtomic(x) || !is.numeric(x)) 
      stop(deparseText(substitute(x)), " must be a vector of class ", 
           isType)
  }
  else {
    if (!isVectorAtomic(x) || !eval(parse(text = paste("is.", 
                                                       isType, "(x)", sep = "")))) 
      stop(deparseText(substitute(x)), " must be a vector of class ", 
           isType)
  }
  invisible(NULL)
}

