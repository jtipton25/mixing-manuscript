broken_stick <- function (p_rel) {
  n <- length(p_rel)
  p <- rep(0, n+1)
  p[1] <- p_rel[1]
  for (i in 2:n) {
    # p(i) <- p_rel(i) * prod((1.0 - p_rel(span(0, i-1))));
    # equivalent expression, more stable to calculate as a sum rather than a product
    p[i] <- p_rel[i] * (1.0 - sum(p[1:i-1]))
  }
  p[n+1] <- 1.0 - sum(p[1:n])
  return(p);
}


broken_stick_jacobian <- function (p, p_rel, logd = true) {
  ## vector input are vectors p, and p_rel 
  ## (inputs and outputs of broken sticks)
  n <- length(p)
  out <- rep(0, n-1)
  out[1]  <- log(p_rel[1]) + log(1.0 - p_rel[1])
  p_increment <- p[1]
  for (i in 2:(n-1)) {
    out[i] <- log(p_rel[i]) + log(1.0 - p_rel[i]) + log(1.0 - p_increment)
    p_increment <- p_increment + p[i]
  }
  if (logd) {
    return(sum(out));
  } else {
    return(exp(sum(out)));
  }
}
