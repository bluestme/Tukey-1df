library(jmuOutlier)

calc_main0 <- function(N, p, q)
{
  beta1 = c(rep(0, 5)) 
  beta2 = c(rep(0.2, 5)) 
  beta = c(beta1, beta2)
  alpha = 3
  true.gamma = 1
  tmp = rnorm(N * (p + q))
  S = matrix(tmp, nrow = N, ncol = (p + q))
  S1 = S[,1:p]
  S2 = S[,(p + 1):(p + q)]
  main = S %*% beta + alpha + true.gamma * 
    ((S1 %*% beta1) * (S2 %*% beta2))
  return (main)
}

Type1_simu <- function(p = 5, q = 5, 
                       dist, loss, 
                       N, N0, B)
{
  pvalue = rep(0, N0)
  gamma_grid = seq(-5, 5, by = 0.01)
  if (dist == "Cauchy")
  {
    for (i in 1:N0)
    {
      main <- calc_main0(N, p, q)
      res <- rcauchy(N, location = 0, scale = 3)
      Y <- main + res
      pv <- permute(Y, S1, S2, B, loss, gamma_grid)
      pvalue[i] <- mean(pv$null > pv$test)
      print(i)
    }
    return (pvalue)
  }
  if (dist == "Laplace")
  {
    for (i in 1:N0)
    {
      main <- calc_main0(N, p, q)
      res <- rlaplace(N, 0, 3 / sqrt(2))
      Y <- main + res
      pv <- permute(Y, S1, S2, B, loss, gamma_grid)
      pvalue[i] <- mean(pv$null > pv$test)
      print(i)
    }
    return (pvalue)
  }
  if (dist == "Logistic")
  {
    for (i in 1:N0)
    {
      main <- calc_main0(N, p, q)
      res <- rlogis(N, 0, 3 * sqrt(3) / pi)
      Y <- main + res
      pv <- permute(Y, S1, S2, B, loss, gamma_grid)
      pvalue[i] <- mean(pv$null > pv$test)
      print(i)
    }
    return (pvalue)
  }
  if (dist == "Mixed1")
  {
    for (i in 1:N0)
    {
      main <- calc_main0(N, p, q)
      res <- mixed(N, 0.1, 3, 10)
      Y <- main + res
      pv <- permute(Y, S1, S2, B, loss, gamma_grid)
      pvalue[i] <- mean(pv$null > pv$test)
      print(i)
    }
    return (pvalue)
  }
  if (dist == "Mixed2")
  {
    for (i in 1:N0)
    {
      main <- calc_main0(N, p, q)
      res <- mixed(N, 0.2, 3, 20)
      Y <- main + res
      pv <- permute(Y, S1, S2, B, loss, gamma_grid)
      pvalue[i] <- mean(pv$null > pv$test)
      print(i)
    }
  }
  if (dist == "Mixed3")
  {
    for (i in 1:N0)
    {
      main <- calc_main0(N, p, q)
      res <- mixed(N, 0.2, 3, 70)
      Y <- main + res
      pv <- permute(Y, S1, S2, B, loss, gamma_grid)
      pvalue[i] <- mean(pv$null > pv$test)
      print(i)
    }
  }
  if (dist == "Binomial")
  {
    for (i in 1:N0)
    {
      miu <- calc_main0(N, p, q)
      prob <- exp(miu) / (1 + exp(miu))
      Y <- rbinom(N, 1, prob)
      pv <- permute(Y, S1, S2, B, loss, gamma_grid)
      pvalue[i] <- mean(pv$null > pv$test)
    }
  }
}