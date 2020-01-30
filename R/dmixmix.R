dmixmix <- nimbleFunction(
  run = function(x = double(), phi = double(1), p = double(2), mu = double(2), sigma = double(2), log = logical(0, default=0)) {
    returnType(double())
    d = nimDim(p)[1]
    Kmax = nimDim(p)[2]
    # out = 0
    # for (j in 1:d) {
    #   out = out + phi[j] * dmixNorm(x, p[j, 1:Kmax], mu[j, 1:Kmax], sigma[j, 1:Kmax], log=FALSE)
    # }
    # out = sum(p * dnorm(x, mu, sigma))
    mixN <- nimRep(0, d)
    for (j in 1:d){
      mixN[j] = dmixNorm(x, p[j, ], mu[j, ], sigma[j, ], TRUE)
    }
    tmp = nimRep(0, d)
    tmp[1:d] = log(phi[1:d]) + mixN[1:d]
    a <- max(tmp)
    out <- a + log(sum(exp(tmp[1:d]-a)))
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

rmixmix <- nimbleFunction(
  run = function(n = double(), phi = double(1), p = double(2), mu = double(2), sigma = double(2)) {
    returnType(double())
    d = dim(p)[1]
    Kmax = dim(p)[2]
    idx <- rcat(1, phi)
    out = rmixNorm(1, p[idx, 1:Kmax], mu[idx, 1:Kmax], sigma[idx, 1:Kmax])
    # return(sum(p * rnorm(1, mu, sigma)))
    return(out)
  })

registerDistributions(list(
  dmixmix = list(
    BUGSdist = "dmixmix(phi, p, mu, sigma)",
    discrete = FALSE,
    range = c(-Inf, Inf),
    types = c('value = double()', 'phi = double(1)', 'p = double(2)', 'mu = double(2)', 'sigma = double(2)')
  )))

##
## Shared Kernels
##
dmixmixShared <- nimbleFunction(
  run = function(x = double(), phi = double(1), p = double(2), mu = double(1), sigma = double(1), log = logical(0, default=0)) {
    returnType(double())
    B = length(phi)
    K = nimDim(p)[2]
    # out = 0
    # for (b in 1:B) {
    #   out = out + phi[b] * dmixNorm(x, p[b, 1:Kmax], mu[1:Kmax], sigma[1:Kmax], log=FALSE)
    # }
    mixN <- nimRep(0, B)
    for (b in 1:B){
      mixN[b] = dmixNorm(x, p[b, 1:K], mu[1:K], sigma[1:K], log=TRUE)
    }
    tmp = nimRep(0, B)
    tmp[1:B] = log(phi[1:B]) + mixN[1:B]
    a <- max(tmp)
    out <- a + log(sum(exp(tmp[1:B] - a)))
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

rmixmixShared <- nimbleFunction(
  run = function(n = double(), phi = double(1), p = double(2), mu = double(1), sigma = double(1)) {
    returnType(double())
    B = length(phi)
    Kmax = nimDim(p)[2]
    out = 0
    idx <- rcat(1, phi)
    out = rmixNorm(1, p[idx, 1:Kmax], mu[1:Kmax], sigma[1:Kmax])
    # return(sum(p * rnorm(1, mu, sigma)))
    return(out)
  })

registerDistributions(list(
  dmixmixShared = list(
    BUGSdist = "dmixmixShared(phi, p, mu, sigma)",
    discrete = FALSE,
    range = c(-Inf, Inf),
    types = c('value = double()', 'phi = double(1)', 'p = double(2)', 'mu = double(1)', 'sigma = double(1)')
  )))
