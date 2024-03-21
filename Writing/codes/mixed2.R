mixed <- function(N, P, var1, var2)
{
  y0 = rnorm(N, 0, var2)
  y1 = rnorm(N, 0, var1)
  index = sample(c(0, 1), N, 
                 replace = T, 
                 prob = c(P, (1 - P)))
  out = y0 * index + y1 * (1 - index)
  return (out)
}
# generate a mixed distribution with P portion of var1 and (1-P) portion of var2