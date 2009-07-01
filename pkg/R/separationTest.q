separationTest <- function(X, y, method = c("dual", "primal", "reduced"))
{
  method <- match.arg(method)

  X <- as.matrix(X)
  dimnames(X) <- NULL
  n <- dim(X)[1]
  p <- dim(X)[2]
  y[y != 1] <- -1

  soln <- switch(method,

    "primal" = {
      mXbar <- -y * X
      obj <- colSums(mXbar)
      lb <- rep(-Inf, p)
      ub = rep(Inf, p)
      lp.control <- list(pivoting = c("firstindex"))
      lpSolve(obj, A = mXbar, double(n), lb = lb, ub = ub, control = lp.control)
    },

    "dual" = {
      mXbarT <- t(-y * X)
      b <- -rowSums(mXbarT)
      obj <- rep(0.0, n)
      lp.control <- list(simplex.type = c("primal", "primal"), pivoting = c("firstindex"))
      lpSolve(obj, A = NULL, b = NULL, Aeq = mXbarT, beq = b, control = lp.control)
    },

    "reduced" = {
      A <- matrix(0.0, n+1, 2*p)
      A[-(n+1), 1:p] <- -y*X
      A[n+1, 1:p] <- colSums(A[-(n+1),1:p])
      A[1:p, -(1:p)] <- diag(p)
      btilde <- rep(0.0, n-p)
      ub <- rep(Inf, n-p)
      pivotingRule <- 0

      .C("reducedLP",
         A = A,
         ldA = as.integer(n+1),
         p = as.integer(p),
         btilde = as.double(btilde),
         ub = as.double(ub),
         beta = double(p+n),
         pivotingRule = as.integer(pivotingRule),
         status = as.integer(1),
         NAOK = TRUE)
    })

  soln$status != 0
}


