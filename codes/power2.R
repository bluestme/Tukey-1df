library(jmuOutlier)

calc_main <- function(N, beta1, p, q)
{
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

Power_simu <- function(p = 5, q = 5, 
                       dist, loss, 
                       N, N0, B)
{
  all_beta1 = seq(0.1, 0.8, by = 0.1)
  pvaluetmp = rep(0, N0)
  pvalue = data.frame(matrix(data = NA,
                             nrow = N0, 
                             ncol = 8))
  gamma_grid = seq(-5, 5, by = 0.01)
  if (dist == "Cauchy")
  {
    for (j in 1:length(all_beta1))
    {
      for (i in 1:N0)
      {
        beta1 = rep(all_beta1[j], 5)
        main = calc_main(N, beta1, p, q)
        res = rcauchy(N, location = 0, scale = 3)
        Y = main + res
        p = permute(Y, S1, S2, B, 
                    loss, gamma_grid)
        pvaluetmp[i] = mean(p$null > p$test)
        print(i)
      }
      pvalue[,j] = pvaluetmp
    }
    return (pvalue)
  }
  if (dis == "Laplace")
  {
    for (j in 1:length(beta1))
    {
      for (i in 1:N0)
      {
        beta1 = rep(all_beta1[j], 5)
        main = calc_main(N, beta1, p, q)
        res = rlaplace(N, 0, 3/sqrt(2))
        Y = main + res
        p = permute(Y, S1, S2, B, loss, gamma_grid)
        pvaluetmp[i] = mean(p$null > p$test)
        print(i)
      }
      pvalue[,j] = pvaluetmp
    }
    return (pvalue)
  }
  if (dis == "Logistic")
  {
    for (j in 1:length(beta1))
    {
      for (i in 1:N0)
      {
        beta1 = rep(all_beta1[j], 5)
        main = calc_main(N, beta1, p, q)
        res = rlogis(N, 0, 3 * sqrt(3) / pi)
        Y = main + res
        p = permute(Y, S1, S2, B, loss, gamma_grid)
        pvaluetmp[i] = mean(p$null > p$test)
        print(i)
      }
      pvalue[,j] = pvaluetmp
    }
    return (pvalue)
  }
  if (dis == "Mixed1")
  {
    for (j in 1:length(beta1))
    {
      for (i in 1:300)
      {
        beta1 = rep(all_beta1[j], 5)
        main = calc_main(N, beta1, p, q)
        res = mixed(N, 0.1, 3, 10)
        Y = main + res
        p = permute(Y, S1, S2, B, loss, gamma_grid)
        pvaluetmp[i] = mean(p$null > p$test)
        print(i)
      }
      pvalue[,j] = pvaluetmp
    }
    return (pvalue)
  }
  if (dis=="Mixed2")
  {
    for (j in 1:length(beta1))
    {
      for (i in 1:300)
      {
        beta1 = rep(all_beta1[j], 5)
        main = calc_main(N, beta1, p, q)
        res = mixed(N, 0.2, 3, 20)
        Y = main + res
        p = permute(Y, S1, S2, 200, loss, gamma_grid)
        pvaluetmp[i] = mean(p$null > p$test)
        print(i)
      }
      pvalue[,j] = pvaluetmp
    }
    return (pvalue)
  }
  if (dis == "Mixed3")
  {
    for (j in 1:length(beta1))
    {
      for (i in 1:300)
      {
        beta1 = rep(all_beta1[j], 5)
        main = calc_main(N, beta1, p, q)
        res = mixed(N, 0.2, 3, 70)
        Y = main + res
        p = permute(Y, S1, S2, 200, loss, gamma_grid)
        pvaluetmp[i] = mean(p$null > p$test)
        print(i)
      }
      pvalue[,j] = pvaluetmp
    }
    return (pvalue)
  }
}



