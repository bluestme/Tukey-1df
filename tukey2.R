library(fastmatrix)
library(MASS)
library(L1pack)

tukey <- function(Y, S1, S2, loss, gamma.grid)
{
  p <- ncol(S2)
  q <- ncol(S1)
  S <- cbind(S1, S2)
  N <- length(Y)
  tst.stat <- gamma.grid
  grid <- length(tst.stat)
  if (loss == "absolute")
  {
    sigma <- rlm(Y ~ S2, psi = psi.huber, scale.est="proposal 2")$s
    fi <- l1fit(S2, Y)$coefficients
    beta2.est <- fi[2:(p + 1)]
    for (i in 1:grid)
    {
      gamma <- gamma.grid[i]
      w <- 1 + gamma * S2 %*% beta2.est
      res <- Y - S2 %*% beta2.est - fi[1]
      score <- -(t(S1) %*% (sign(res)*w)) / (2 * sigma^2)
      info <- (t(S1) %*% diag(c(w))) %*% t(t(S1) %*% diag(c(w))) / (4 * sigma^4)
      tst.stat[i] <- t(score) %*% solve(info) %*% score
    }
    return (list(TG = gamma.grid[which.max(tst.stat)],
                 Tscore = max(tst.stat),sc = tst.stat))
  }
  if (loss == "logistic")
  {
    fi <- glm(Y ~ S2, family = binomial(link="logit"))$coefficients
    beta2.est <- fi[2:(p+1)]
    prob <- exp(cbind(1, S2) %*% fi) / (1 + exp(cbind(1, S2) %*% fi))
    for (l in 1:grid)
    {
      gamma <- gamma.grid[l]
      res <- (Y - prob)
      w <- (1 + gamma * (S2 %*% beta2.est))
      score <- t(S1) %*% (w * res)
      Z <- cbind(1, S2)
      I11 <- (t(S1) %*% diag(c(w * (1 - prob)))) %*% t(t(S1) %*% diag(c(w * prob)))
      I12 <- (t(S1) %*% diag(c(w*(1 - prob) * prob)))%*%(Z)
      I22 <- (t(Z) %*% diag(c(1 - prob)))%*%t(t(Z)%*%diag(c(prob)))
      info <- solve(I11 - I12 %*% solve(I22) %*% t(I12))
      tst.stat[l] <- t(score) %*% info %*% score
    }
    return (list(TG = gamma.grid[which.max(tst.stat)],
                 Tscore = max(tst.stat), sc = tst.stat))
  }
  if (loss == "squared")
  {
    n <- length(Y)
    fi <- lm(Y ~ S2)$coefficients
    beta2.est <- fi[2:(p+1)]
    sigma <- summary(lm(Y~S2))$sigma
    for (i in 1:grid)
    {
      gamma <- gamma.grid[i]
      w <- 1 + gamma*S2 %*% beta2.est
      res <- Y - S2 %*% beta2.est - fi[1]
      score <- (t(S1) %*% (res*w))/(sigma^2)
      M <- rbind(1, t(S1) %*% diag(c(w)), t(S2), 0)
      I <- M %*% t(M)
      I[nrow(I), ncol(I)] <- n/(2*sigma^2)
      I <- I/sigma^2
      n <- length(Y)
      ind <- 2:(p + 1)
      ind2 <- c(1:(p + q + 2))[-ind]
      I11 <- I[ind, ind]
      I22 <- I[ind2, ind2]
      I12 <- I[ind, ind2]
      I.adj <- I11 - I12 %*% solve(I22) %*% t(I12)
      stat <- t(score) %*% solve(I.adj) %*% score
      tst.stat[i] <- stat
    }
    return (list(TG = gamma.grid[which.max(tst.stat)],
                 Tscore = max(tst.stat),sc = tst.stat))
  }
}
