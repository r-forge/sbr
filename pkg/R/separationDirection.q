separationDirection <- function(X, y, method = c("primal", "reduced"),
                                check = FALSE, tolerance = 1e-9)
{
  method <- match.arg(method)

  X <- as.matrix(X)
  dimnames(X) <- NULL
  n <- nrow(X)
  p <- ncol(X)
  y[y != 1] <- -1

  beta <- switch(method,

    primal = {
      mXbar <- -y * X

      obj <- colSums(mXbar)
      lp.control <- list(pivoting = c("firstindex"))
      lb <- rep(-1.0, p)
      ub = rep(1.0, p)
  
      lpSolve(obj, A = mXbar, b = double(n), lb = lb, ub = ub, control = lp.control)$x
    },

    reduced = {
      A <- matrix(0.0, n+1, 2*p)
      A[-(n+1), 1:p] <- -y*X
      A[n+1, 1:p] <- colSums(A[-(n+1), 1:p])
      A[1:p, -(1:p)] <- diag(p)
      btilde <- rep(0.0, n-p)
      ub <- rep(1.0, n-p)
      pivotingRule <- 0

      soln <- .C("reducedLP",
                  A = A,
                  ldA = as.integer(n+1),
                  p = as.integer(p),
                  btilde = as.double(btilde),
                  ub = as.double(ub),
                  beta = double(p+n),
                  pivotingRule = as.integer(pivotingRule),
                  status = as.integer(1),
                  NAOK = TRUE)

      soln$beta[1:p]
    })

  if(check && sqrt(sum(beta^2)) > tolerance) {
    pred <- X %*% matrix(beta, ncol = 1)
    cs <- which(abs(pred) > tolerance)
    if(all.equal(sign(pred[cs]), y[cs]))
      cat("\n    The computed beta quasiseparates the sample points.\n\n")
    else
      stop("the computed beta does not separate the sample points")
  }

  beta
}


