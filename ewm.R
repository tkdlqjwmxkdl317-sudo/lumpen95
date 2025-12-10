ewma_adjust <- function(x, alpha = NULL, span = NULL, adjust = FALSE) {
  
  if (!is.null(alpha) && !is.null(span)) {
    stop("alpha 와 span 중 하나만 입력해야 합니다.")
  }
  
  if (is.null(alpha) && is.null(span)) {
    stop("alpha 또는 span 중 하나를 반드시 입력해야 합니다.")
  }
  
  if (!is.null(span)) {
    alpha <- 2 / (span + 1)
  }
  
  if (!adjust) {
    res <- accumulate(x, function(prev, xi) alpha * xi + (1 - alpha) * prev)
  } else {
    n <- length(x)
    res <- numeric(n)
    
    for (t in 1:n) {
      w <- (1 - alpha)^(t:1 - 1) # 가중치 벡터
      res[t] <- sum(w * x[1:t]) / sum(w)
    }
  }
  return(res)
}

ewmsd_adjust <- function(x, alpha = NULL, span = NULL, adjust = FALSE) {
  
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
  ewma_mean <- ewma_adjust(x, alpha = alpha, adjust = adjust)
  ewma_var <- numeric(n)
  ewma_sd <- numeric(n)
  
  if (!adjust) {
    for (t in 2:n) {
      diff_t <- x[t] - ewma_mean[t - 1]
      ewma_var[t] <- alpha * diff_t^2 + (1 - alpha) * ewma_var[t - 1]
      ewma_sd[t] <- sqrt(ewma_var[t])
    }
  } else {
    for (t in 1:n) {
      w <- (1 - alpha)^(t:1 - 1) # 가중치 벡터
      centered_sq <- (x[1:t] - ewma_mean[t])^2
      ewma_var[t] <- sum(w * centered_sq) / sum(w)
      ewma_sd[t] <- sqrt(ewma_var[t])
    }
  }
  return(ewma_sd)
}
