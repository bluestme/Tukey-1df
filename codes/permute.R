tukey_permute <- function(Y, S1, S2, B, loss, gamma)
{
  test.stat <- tukey(Y, S1, S2, 
                     loss, gamma)$Tscore
  Y.star <- Y
  null.stat <- rep(0, B)
  S <- cbind(1, S2)
  fi <- l1fit(S2, Y)$coefficients
  res <- Y - S%*%fi
  for (i in 1:B)
  {
    res <- sample(res, length(res), 
                  replace = FALSE)
    Y.star <- res + S%*%fi
    null.stat[i] <- tukey(Y.star, S1, S2, 
                          loss, gamma)$Tscore
  }
  return(list(null = null.stat, test = test.stat))
}