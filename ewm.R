ewma <- function(x, alpha = NULL, span = NULL) {
  
  if (!is.null(alpha) && !is.null(span)) {
    stop("alpha 와 span 중 하나만 입력해야 합니다.")
  }
  
  if (is.null(alpha) && is.null(span)) {
    stop("alpha 또는 span 중 하나를 반드시 입력해야 합니다.")
  }
  
  if (!is.null(span)) {
    alpha <- 2 / (span + 1)
  }

  n <- length(x)
  res <- numeric(n)
    
  for (t in 1:n) {
    w <- (1 - alpha)^(t:1 - 1) # 가중치 벡터
    res[t] <- sum(w * x[1:t]) / sum(w)
  }
  return(res)
}

ewmsd <- function(x, alpha = NULL, span = NULL, bias = FALSE) {
  
  if (!is.null(alpha) && !is.null(span)) {
    stop("alpha 와 span 중 하나만 입력해야 합니다.")
  }
  if (is.null(alpha) && is.null(span)) {
    stop("alpha 또는 span 중 하나를 반드시 입력해야 합니다.")
  }
  if (!is.null(span)) {
    alpha <- 2 / (span + 1)
  }
  
  n <- length(x)
  ewma_mean <- ewma(x, alpha = alpha)
  res <- numeric(n)

  for (t in 1:n) {
    w <- (1 - alpha)^(t:1 - 1) # 가중치 벡터
    centered_sq <- (x[1:t] - ewma_mean[t])^2
    
    w_sum  <- sum(w)
    w_sum2 <- sum(w^2)
    
    if (bias) {
      denom <- w_sum
    } else {
      denom <- w_sum - w_sum2 / w_sum
    }
    res[t] <- sqrt(sum(w * centered_sq) / denom)
  }
  return(res)
}
