dmixNorm <- nimbleFunction(
  run = function(x = double(), p = double(1), mu = double(1), sigma = double(1), log = logical(0, default=0)) {
    returnType(double())
    K = length(p)
    # out = 0
    # for (j in 1:d) {
    #   out = out + p[j] * dnorm(x, mu[j], sigma[j])
    # }
    tmp = rep(0, K)
    tmp[1:K] = log(p[1:K]) + dnorm(x, mu[1:K], sigma[1:K], log=TRUE)
    a = max(tmp[1:K])
    out = a + log(sum(exp(tmp[1:K] - a)))
    # out = sum(p * dnorm(x, mu, sigma))
    # if (log) {
    #   return(log(out))
    # } else {
    #   return(out)
    # }
    if (log) {
      return(out)
    } else {
      return(exp(out))
    }
  })


# dmixNorm <- nimbleFunction(
#   run = function(x = double(), p = double(1), mu = double(1), sigma = double(1), log = logical(0, default=0)) {
#     returnType(double())
#     d = length(p)
#     tmp = numeric(length=d)
#     for (j in 1:d) {
#       tmp[j] = log(p[j]) + dnorm(x, mu[j], sigma[j], log = TRUE)
#     }
#     A = max(tmp)
#     out <- A + log(sum(exp(tmp[j] - A)))
#     if (log) {
#       return(out)
#     } else {
#       return(exp(out))
#     }
#   })

rmixNorm <- nimbleFunction(
  run = function(n = double(), p = double(1), mu = double(1), sigma = double(1)) {
    returnType(double())
    d = length(p)
    out = rep(0, d)
    idx <- rcat(1, prob=p)
    for (j in 1:d) {
      out[j] = rnorm(1, mu[j], sigma[j])
    }
    # return(sum(p * rnorm(1, mu, sigma)))
    return(out[idx])
  })

registerDistributions(list(
  dmixNorm = list(
    BUGSdist = "dmixNorm(p, mu, sigma)",
    discrete = FALSE,
    range = c(-Inf, Inf),
    types = c('value = double()', 'p = double(1)', 'mu = double(1)', 'sigma = double(1)')
  )))
