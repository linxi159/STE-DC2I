
# [0]Main function：dark_causality
dark_causality <- function(input_X, input_Y, embedding_E, delay_tau, horizon_h, is_weighted, 
                           dist_metric = "euclidean", custom_dist_fn = NULL,
                           custom_space_fn = NULL, use_relative = TRUE, show_progress = FALSE)
{
  Zz0Zz(arg1X = input_X, arg1Y = input_Y, arg1E = embedding_E, arg1T = delay_tau, 
        arg1M = dist_metric, arg1H = horizon_h, arg1W = is_weighted, arg1F = custom_dist_fn)# [0]Main function：dark_causality
  Jj9Jj <- Kk8Kk(arg2E = embedding_E, arg2T = delay_tau)#[2]Sub function 2：initialize_components
  Ss7Ss <- Mm6Mm(arg3X = input_X, arg3Y = input_Y, arg3E = embedding_E, arg3T = delay_tau, arg3M = dist_metric, 
                 arg3F = custom_dist_fn, arg3S = custom_space_fn, 
                 arg3R = use_relative, arg3V = show_progress)#[3]Sub function 3：compute_state_spaces
  Cc5Cc <- Vv4Vv(arg4E = embedding_E, arg4T = delay_tau, arg4H = horizon_h, arg4X = input_X, arg4V = show_progress)#[4]Sub function 4：check_causality_points
  if (!Cc5Cc$valid) {
    stop("Insufficient data length for analysis", call. = FALSE)
  }
  Pp3Pp <- Rr2Rr(arg5X = input_X, arg5Y = input_Y, arg5E = embedding_E, arg5P = Cc5Cc$first_point,
                 arg5V = show_progress)#[5]Sub function 5：initialize_matrice
  Aa1Aa <- Bb0Bb(arg6S = Ss7Ss, arg6M = Pp3Pp, arg6C = Jj9Jj,
                 arg6V = Cc5Cc, arg6H = horizon_h, arg6W = is_weighted, arg6P = show_progress)#[6]Sub function 6: analyze_causality_result * 
  Ll9Ll <- Nn8Nn(arg7R = Aa1Aa, arg7W = is_weighted)#[7]Sub function 7: compute_causality_measures * 
  
  final_result <- structure(list(shadow = Ll9Ll$shadow),
                            class = "causality_fit")#Statistical Results
  
  return(final_result)
}


#[1]Sub function 1：validate_inputs
Zz0Zz <- function(arg1X, arg1Y, arg1E, arg1T, arg1M, arg1H, arg1W, arg1F = NULL) {
  if(!is.numeric(arg1X) || !is.numeric(arg1Y)) {
    stop("X and Y must be numeric vectors", call. = FALSE)
  }
  if(length(arg1X) != length(arg1Y)) {
    stop("X and Y must have the same length", call. = FALSE)
  }
  if(any(is.infinite(arg1X)) || any(is.infinite(arg1Y))) {
    stop("X and Y cannot contain infinite values", call. = FALSE)
  }
  
  if(!is.numeric(arg1E) || arg1E <= 1 || arg1E != round(arg1E)) {
    stop("E must be a positive integer greater than 1", call. = FALSE)
  }
  if(!is.numeric(arg1T) || arg1T < 1 || arg1T != round(arg1T)) {
    stop("tau must be a positive integer", call. = FALSE)
  }
  
  if(is.null(arg1F)) {
    if(!is.character(arg1M) || !arg1M %in% c("euclidean", "manhattan", "maximum")) {
      stop("metric must be one of: 'euclidean', 'manhattan', 'maximum'", call. = FALSE)
    }
  }
  
  if(!is.numeric(arg1H) || arg1H < 0 || arg1H != round(arg1H)) {
    stop("h must be a non-negative integer", call. = FALSE)
  }
  if(!is.logical(arg1W)) {
    stop("weighted must be logical", call. = FALSE)
  }
  
  z1z <- (arg1E - 1) * arg1T + arg1H + 1
  if(length(arg1X) < z1z) {
    stop(sprintf("Time series length must be at least %d for given E, tau and h",
                 z1z), call. = FALSE)
  }
  
  if(!is.null(arg1F) && !is.function(arg1F)) {
    stop("distance_fn must be a function", call. = FALSE)
  }
}


#[2]Sub function 2：initialize_components
Kk8Kk <- function(arg2E, arg2T) {
  list(
    a1a = arg2E + 1,
    b2b = (arg2E - 1) * arg2T,
    c3c = Xx7Xx(arg2E = arg2E)$hashes
  )
}

Xx7Xx <- function(arg2E, arg2V = FALSE) {
  if(!is.numeric(arg2E) || arg2E <= 1 || arg2E != round(arg2E)) {
    stop("E must be an integer greater than 1", call. = FALSE)
  }
  
  if(arg2V) {
    cat("Generating patterns for embedding dimension", arg2E, "\n")
  }
  
  d4d <- as.matrix(expand.grid(rep(list(1:3), arg2E - 1)))
  
  if(arg2V) {
    cat("Computing hash values for", nrow(d4d), "patterns\n")
  }
  
  e5e <- apply(d4d, 1, function(f6f) {
    sum(f6f * factorial(seq_along(f6f) + 2))
  })
  
  Qq6Qq(patterns = d4d, hashes = e5e)
}

Qq6Qq <- function(patterns, hashes) {
  structure(
    list(
      patterns = patterns,
      hashes = hashes
    ),
    class = "pattern_object"
  )
}


#[3]Sub function 3：compute_state_spaces
Mm6Mm <- function(arg3X, arg3Y, arg3E, arg3T, arg3M,
                  arg3F = NULL, arg3S = NULL, arg3R = TRUE,
                  arg3V = FALSE) {
  tryCatch({
    if(arg3V) cat("Computing spaces...\n")
    
    g7g <- if(!is.null(arg3S)) arg3S else Tt5Tt
    h8h <- g7g(arg3X, arg3E, arg3T)$matrix
    i9i <- g7g(arg3Y, arg3E, arg3T)$matrix
    if(arg3V) cat("Done\n")
    
    if(arg3V) cat("  - Computing signature spaces... ")
    j0j <- Uu4Uu(h8h, relative = arg3R)
    k1k <- Uu4Uu(i9i, relative = arg3R)
    if(arg3V) cat("Done\n")
    
    if(arg3V) cat("  - Computing pattern spaces... ")
    l2l <- Ww3Ww(j0j)
    m3m <- Ww3Ww(k1k)
    if(arg3V) cat("Done\n")
    
    if(arg3V) cat("  - Computing distance matrices... ")
    n4n <- if(!is.null(arg3F)) {
      function(o5o) if(!is.matrix(o5o)) as.matrix(arg3F(o5o)) else arg3F(o5o)
    } else {
      function(o5o) as.matrix(stats::dist(o5o, method = arg3M))
    }
    p6p <- n4n(h8h)
    q7q <- n4n(i9i)
    if(arg3V) cat("Done\n")
    
    if(arg3E == 2) {
      j0j <- matrix(j0j, ncol = 1)
      k1k <- matrix(k1k, ncol = 1)
      l2l <- matrix(l2l, ncol = 1)
      m3m <- matrix(m3m, ncol = 1)
      if(!is.matrix(h8h)) h8h <- matrix(h8h, ncol = arg3E)
      if(!is.matrix(i9i)) i9i <- matrix(i9i, ncol = arg3E)
      if(!is.matrix(p6p)) p6p <- matrix(p6p, ncol = ncol(h8h))
      if(!is.matrix(q7q)) q7q <- matrix(q7q, ncol = ncol(i9i))
    }
    
    list(
      state_X = h8h, state_Y = i9i,
      signature_X = j0j, signature_Y = k1k,
      pattern_X = l2l, pattern_Y = m3m,
      dist_X = p6p, dist_Y = q7q
    )
  }, error = function(r8r) {
    stop("Error in space computation: ", r8r$message, call. = FALSE)
  })
}

Tt5Tt <- function(s9s, t0t, u1u, v2v = FALSE) {
  if (!is.numeric(s9s)) {
    stop("Time series must be numeric", call. = FALSE)
  }
  if (t0t < 2) {
    stop("Embedding dimension must be greater than 1", call. = FALSE)
  }
  if (u1u < 1) {
    stop("Time delay must be positive", call. = FALSE)
  }
  
  if (v2v) {
    cat("Reconstructing state space...\n")
  }
  
  w3w <- length(s9s) - (t0t - 1) * u1u
  x4x <- matrix(NA_real_, w3w, t0t)
  for (y5y in 1:nrow(x4x)) {
    x4x[y5y,] <- s9s[seq(from = y5y, to = y5y + (t0t - 1) * u1u, by = u1u)]
    if (anyNA(x4x[y5y,])) {
      x4x[y5y,] <- rep(NA_real_, t0t)
    }
  }
  
  result <- structure(
    list(
      matrix = x4x,
      parameters = list(
        E = t0t,
        tau = u1u,
        n_points = w3w
      ),
      original = s9s
    ),
    class = "state_space"
  )
  
  return(result)
}

Uu4Uu <- function(z6z, relative = TRUE) {
  if(!is.matrix(z6z)) {
    stop("Input must be a matrix", call. = FALSE)
  }
  
  a7a <- ncol(z6z)
  if(a7a < 2) {
    stop("State space matrix must have at least 2 columns", call. = FALSE)
  }
  
  if(a7a == 2) {
    b8b <- matrix(apply(z6z, 1, compute_vector_differences, relative = relative), ncol = 1)
  } else {
    b8b <- t(apply(z6z, 1, compute_vector_differences, relative = relative))
  }
  
  b8b[is.na(b8b)] <- NA_real_
  
  return(b8b)
}

compute_vector_differences <- function(c9c, relative = TRUE) {
  if(!is.numeric(c9c)) {
    stop("Input must be a numeric vector", call. = FALSE)
  }
  
  if(relative) {
    d0d <- diff(c9c) / c9c[-length(c9c)]
  } else {
    d0d <- diff(c9c)
  }
  
  return(d0d)
}

Ww3Ww <- function(e1e) {
  if (!is.matrix(e1e) || !is.numeric(e1e)) {
    stop("Input must be a numeric matrix", call. = FALSE)
  }
  
  do.call(rbind, sapply(seq_len(nrow(e1e)), function(f2f) convert_to_pattern(e1e[f2f, ]), simplify = FALSE))
}

convert_to_pattern <- function(g3g) {
  if (anyNA(g3g)) {
    return(rep(NA_real_, length(g3g)))
  }
  
  h4h <- ifelse(g3g > 0, 3,
                ifelse(g3g < 0, 1, 2))
  
  compute_pattern_hash(h4h)
}

compute_pattern_hash <- function(i5i) {
  if (anyNA(i5i)) {
    return(NA_real_)
  }
  sum(i5i * factorial(seq_along(i5i) + 2))
}


#[4] Sub function 4：check_causality_points
Vv4Vv <- function(arg4E, arg4T, arg4H, arg4X, arg4V) {
  if(arg4V) cat("Checking causality points...\n")
  
  j6j <- Yy2Yy(arg4E, arg4T, arg4H, arg4X)
  if(!j6j$valid) {
    return(list(
      valid = FALSE,
      first_point = NA_real_,
      analysis_range = NA_real_
    ))
  }
  
  k7k <- Zz1Zz(arg4E, arg4T, arg4H, arg4X)
  if(is.null(k7k) || is.na(k7k$point)) {
    return(list(
      valid = FALSE,
      first_point = NA_real_,
      analysis_range = NA_real_
    ))
  }
  
  l8l <- length(arg4X) - (arg4E - 1) * arg4T - arg4H
  if(k7k$point > l8l) {
    return(list(
      valid = FALSE,
      first_point = NA_real_,
      analysis_range = NA_real_
    ))
  }
  
  m9m <- k7k$point:l8l
  
  list(
    valid = TRUE,
    first_point = k7k$point,
    analysis_range = m9m
  )
}

Yy2Yy <- function(n0n, o1o, p2p, q3q, r4r = FALSE) {
  if(!is.numeric(n0n) || n0n <= 0 || n0n != round(n0n)) {
    stop("E must be a positive integer", call. = FALSE)
  }
  
  if(!is.numeric(o1o) || o1o <= 0 || o1o != round(o1o)) {
    stop("tau must be a positive integer", call. = FALSE)
  }
  
  if(!is.numeric(p2p) || p2p < 0 || p2p != round(p2p)) {
    stop("h must be a non-negative integer", call. = FALSE)
  }
  
  if(!is.numeric(q3q)) {
    stop("X must be a numeric vector", call. = FALSE)
  }
  
  s5s <- n0n + 1
  t6t <- n0n - 1
  u7u <- p2p
  
  v8v <- s5s + 2 * t6t + u7u
  w9w <- length(q3q)
  
  if(r4r) {
    cat("Checking causality analysis feasibility:\n")
    cat("Required length:", v8v, "\n")
    cat("Available length:", w9w, "\n")
  }
  
  Aa0Aa(
    feasible = (s5s + t6t + u7u < w9w - t6t),
    required_length = v8v,
    available_length = w9w,
    parameters = list(
      E = n0n,
      tau = o1o,
      h = p2p
    )
  )
}

Aa0Aa <- function(feasible, required_length, available_length, parameters) {
  structure(
    list(
      feasible = feasible,
      required_length = required_length,
      available_length = available_length,
      parameters = parameters
    ),
    class = "feasibility_check"
  )
}

Zz1Zz <- function(x1x, y2y, z3z, a4a, b5b = FALSE) {
  if(!is.numeric(x1x) || x1x <= 0 || x1x != round(x1x)) {
    stop("E must be a positive integer", call. = FALSE)
  }
  
  if(!is.numeric(y2y) || y2y <= 0 || y2y != round(y2y)) {
    stop("tau must be a positive integer", call. = FALSE)
  }
  
  if(!is.numeric(z3z) || z3z < 0 || z3z != round(z3z)) {
    stop("h must be a non-negative integer", call. = FALSE)
  }
  
  if(!is.numeric(a4a)) {
    stop("X must be a numeric vector", call. = FALSE)
  }
  
  c6c <- x1x + 1
  d7d <- (x1x - 1) * y2y
  e8e <- z3z
  
  if(b5b) {
    cat("Computing first causality point:\n")
    cat("Nearest neighbor span:", c6c, "\n")
    cat("Common coordinate span:", d7d, "\n")
    cat("Prediction span:", e8e, "\n")
  }
  
  f9f <- 1 + c6c + d7d + e8e
  
  if(c6c + d7d + e8e >= length(a4a) - d7d) {
    stop(
      sprintf(
        paste(
          "Insufficient data for causality analysis.",
          "Required length: %d, Available length: %d\n",
          "Check parameters: E=%d, tau=%d, h=%d"
        ),
        c6c + 2 * d7d + e8e,
        length(a4a),
        x1x, y2y, z3z
      ),
      call. = FALSE
    )
  }
  
  Bb0Bb(
    point = f9f,
    spans = list(
      nn = c6c,
      cc = d7d,
      pred = e8e
    )
  )
}

Bb0Bb <- function(point, spans) {
  structure(
    list(
      point = point,
      spans = spans
    ),
    class = "causality_point"
  )
}


#[5]Sub function 5：initialize_matrice
Rr2Rr <- function(arg5X, arg5Y, arg5E, arg5P, arg5V) {
  if(arg5V) cat("Initializing matrices...\n")
  
  g1g <- length(arg5Y)
  h2h <- 3^(arg5E-1)
  
  i3i <- Cc9Cc(type = "array",
               dimensions = c(h2h, h2h, g1g))
  j4j <- Cc9Cc(type = "array",
               dimensions = c(h2h, h2h, g1g))
  
  k5k <- Cc9Cc(type = "matrix",
               dimensions = c(g1g, arg5E-1))
  l6l <- Cc9Cc(type = "matrix",
               dimensions = c(g1g, arg5E-1))
  m7m <- Cc9Cc(type = "matrix",
               dimensions = c(g1g, arg5E-1))
  
  n8n <- Cc9Cc(type = "vector", dimensions = g1g)
  o9o <- Cc9Cc(type = "vector", dimensions = g1g)
  p0p <- Cc9Cc(type = "vector", dimensions = g1g)
  
  q1q <- structure(
    Cc9Cc(type = "matrix", dimensions = c(g1g, arg5E)),
    dimnames = list(NULL, c("current", rep("predicted", arg5E-1)))
  )
  r2r <- structure(
    Cc9Cc(type = "matrix", dimensions = c(g1g, arg5E)),
    dimnames = list(NULL, c("current", rep("future", arg5E-1)))
  )
  
  list(
    causality_matrices = list(
      predicted = i3i,
      actual = j4j
    ),
    signatures = list(
      predicted = k5k,
      actual = l6l,
      causal = m7m
    ),
    patterns = list(
      predicted = n8n,
      actual = o9o,
      causal = p0p
    ),
    values = list(
      predicted = q1q,
      actual = r2r
    )
  )
}

Cc9Cc <- function(type, dimensions, s3s = FALSE) {
  if(!is.character(type) || length(type) != 1) {
    stop("'type' must be a single character string", call. = FALSE)
  }
  
  if(!is.numeric(dimensions) || any(dimensions <= 0) || any(dimensions != round(dimensions))) {
    stop("'dimensions' must be a vector of positive integers", call. = FALSE)
  }
  
  t4t <- switch(type,
                "array" = {
                  if(s3s) cat("Creating array with dimensions:", paste(dimensions, collapse = "x"), "\n")
                  array(NA_real_, dim = dimensions)
                },
                "vector" = {
                  if(s3s) cat("Creating vector of length:", dimensions[1], "\n")
                  rep(NA_real_, dimensions[1])
                },
                "matrix" = {
                  if(s3s) cat("Creating matrix with dimensions:", paste(dimensions[1:2], collapse = "x"), "\n")
                  matrix(NA_real_, nrow = dimensions[1], ncol = dimensions[2])
                },
                stop("Invalid type specified", call. = FALSE)
  )
  
  return(t4t)
}


#[6]Sub function 6: analyze_causality_result
Bb0Bb <- function(arg6S, arg6M, arg6C, arg6V, arg6H, arg6W, arg6P) {
  
  u5u <- numeric(0)
  
  for(v6v in seq_along(arg6V$analysis_range)) {
    w7w <- arg6V$analysis_range[v6v]
    
    if(!anyNA(c(arg6S$state_X[w7w,], arg6S$state_Y[w7w + arg6H,]))) {
      x8x <- Dd8Dd(
        coord_span = arg6C$b2b,
        neighbor_span = arg6C$a1a,
        state_X = arg6S$state_X,
        dist_X = arg6S$dist_X,
        signature_X = arg6S$signature_X,
        pattern_X = arg6S$pattern_X,
        index = w7w,
        horizon = arg6H
      )
      
      if(!anyNA(x8x$distances) && !anyNA(arg6S$dist_Y[w7w, x8x$indices + arg6H])) {
        u5u <- c(u5u, w7w)
        
        y9y <- Ee7Ee(
          state_Y = arg6S$state_Y,
          dist_Y = arg6S$dist_Y,
          signature_Y = arg6S$signature_Y,
          pattern_Y = arg6S$pattern_Y,
          neighbor_indices = x8x$indices,
          index = w7w,
          horizon = arg6H
        )
        
        arg6M <- Ff6Ff(
          arg6M, arg6S, x8x, y9y,
          w7w, arg6H, arg6W, arg6P, arg6C$c3c
        )
      }
    }
    
    if(arg6P) {
      Gg5Gg(v6v, length(arg6V$analysis_range), "Analyzing causality patterns", arg6P)
    }
  }
  
  if(arg6P) {
    cat("\nComputing final results...\n")
  }
  
  if(length(u5u) > 0) {
    z0z <- Hh4Hh(arg6M$causality_matrices, u5u, arg6C$c3c, arg6S$state_X[,1])
    results <- list(
      no_causality = z0z$predicted$no_causality,
      positive = z0z$predicted$positive,
      negative = z0z$predicted$negative,
      shadow = z0z$predicted$shadow,
      valid_points = u5u
    )
  } else {
    results <- list(
      no_causality = numeric(0),
      positive = numeric(0),
      negative = numeric(0),
      shadow = numeric(0),
      valid_points = u5u
    )
  }
  
  results
}

Dd8Dd <- function(coord_span, neighbor_span, state_X, dist_X, signature_X, pattern_X, index, horizon, a1a = FALSE) {
  if(!is.numeric(coord_span) || coord_span < 0 || coord_span != round(coord_span)) {
    stop("CCSPAN must be a non-negative integer", call. = FALSE)
  }
  
  if(!is.numeric(neighbor_span) || neighbor_span <= 0 || neighbor_span != round(neighbor_span)) {
    stop("NNSPAN must be a positive integer", call. = FALSE)
  }
  
  if(!is.matrix(state_X) || !is.matrix(dist_X) || !is.matrix(signature_X)) {
    stop("State matrices must be matrices", call. = FALSE)
  }
  
  if(index <= coord_span + horizon) {
    stop("Insufficient past data for the given parameters", call. = FALSE)
  }
  
  if(a1a) {
    cat("Finding nearest neighbors for point", index, "\n")
  }
  
  b2b <- dist_X[index, 1:(index - coord_span - horizon)]
  
  c3c <- order(b2b)[1:min(neighbor_span, length(b2b))]
  d4d <- as.numeric(names(b2b[c3c]))
  e5e <- b2b[c3c]
  
  if(a1a) {
    cat("Found", length(d4d), "nearest neighbors\n")
  }
  
  Ii3Ii(
    index = index,
    indices = d4d,
    distances = e5e,
    signatures = signature_X[d4d, ],
    patterns = pattern_X[d4d],
    coordinates = state_X[d4d, ]
  )
}

Ii3Ii <- function(index, indices, distances, signatures, patterns, coordinates, weights = NULL) {
  if(!is.numeric(index) || length(index) != 1) {
    stop("index must be a single numeric value", call. = FALSE)
  }
  if(!is.numeric(indices) || !is.numeric(distances)) {
    stop("indices and distances must be numeric vectors", call. = FALSE)
  }
  
  if(is.null(dim(signatures))) {
    signatures <- matrix(signatures, ncol = 1)
  }
  
  if(is.null(dim(coordinates))) {
    coordinates <- matrix(coordinates, ncol = 2)
  }
  
  if(!is.matrix(signatures) || !is.matrix(coordinates)) {
    stop("signatures and coordinates must be matrices", call. = FALSE)
  }
  if(!is.numeric(patterns)) {
    stop("patterns must be a numeric vector", call. = FALSE)
  }
  if(!is.null(weights) && !is.numeric(weights)) {
    stop("weights must be NULL or a numeric vector", call. = FALSE)
  }
  
  structure(
    list(
      index = index,
      indices = indices,
      distances = distances,
      signatures = signatures,
      patterns = patterns,
      coordinates = coordinates,
      weights = weights
    ),
    class = "neighbors_info"
  )
}

Ee7Ee <- function(state_Y, dist_Y, signature_Y, pattern_Y, neighbor_indices, index, horizon) {
  if(!is.matrix(state_Y) || !is.matrix(dist_Y) || !is.matrix(signature_Y) || !is.matrix(pattern_Y)) {
    stop("All inputs must be matrices", call. = FALSE)
  }
  
  f6f <- Jj2Jj(dist_Y[index, neighbor_indices + horizon])
  
  Ii3Ii(
    index = index,
    indices = neighbor_indices + horizon,
    distances = dist_Y[index, neighbor_indices + horizon],
    weights = as.vector(f6f),
    signatures = signature_Y[neighbor_indices + horizon, ],
    patterns = pattern_Y[neighbor_indices + horizon],
    coordinates = state_Y[neighbor_indices + horizon, ]
  )
}

Jj2Jj <- function(g7g) {
  h8h <- g7g
  i9i <- sum(h8h)
  if (i9i == 0) i9i <- 0.0001
  
  h8h <- h8h / i9i
  exp(-h8h) / sum(exp(-h8h))
}

Ff6Ff <- function(j0j, k1k, l2l, m3m, n4n, o5o, p6p, q7q, r8r) {
  s9s <- Kk1Kk(m3m)
  t0t <- s9s$signature
  u1u <- s9s$pattern[1]
  
  v2v <- k1k$signature_X[n4n,]
  w3w <- k1k$pattern_X[n4n]
  
  x4x <- k1k$signature_Y[n4n + o5o,]
  y5y <- k1k$pattern_Y[n4n + o5o]
  
  j0j$signatures$predicted[n4n,] <- t0t
  j0j$signatures$actual[n4n,] <- x4x
  j0j$signatures$causal[n4n,] <- v2v
  
  j0j$patterns$predicted[n4n] <- u1u
  j0j$patterns$actual[n4n] <- y5y
  j0j$patterns$causal[n4n] <- w3w
  
  z6z <- Ll0Ll(p6p, u1u, y5y,
               t0t, x4x,
               w3w, v2v)
  
  if(!is.null(z6z$predicted)) {
    j0j$causality_matrices$predicted[
      which(r8r == w3w),
      which(r8r == u1u),
      n4n
    ] <- z6z$predicted
  }
  if(!is.null(z6z$actual)) {
    j0j$causality_matrices$actual[
      which(r8r == w3w),
      which(r8r == u1u),
      n4n
    ] <- z6z$actual
  }
  
  j0j
}

Ll0Ll <- function(is_weighted, predicted_pattern, actual_pattern,
                  predicted_signature, actual_signature,
                  causal_pattern, causal_signature,
                  a7a = FALSE) {
  if(!is.logical(is_weighted)) {
    stop("weighted must be TRUE or FALSE", call. = FALSE)
  }
  
  if(!is.numeric(c(predicted_pattern, actual_pattern, causal_pattern))) {
    stop("All patterns must be numeric", call. = FALSE)
  }
  
  if(!is.numeric(c(predicted_signature, actual_signature, causal_signature))) {
    stop("All signatures must be numeric vectors", call. = FALSE)
  }
  
  b8b <- NA_real_
  c9c <- NA_real_
  
  if (!anyNA(c(predicted_pattern, actual_pattern, causal_pattern))) {
    if (length(predicted_pattern) > 0 && length(causal_pattern) > 0) {
      if (a7a) {
        cat("Computing causality strengths:\n")
        cat("Predicted pattern:", predicted_pattern, "\n")
        cat("Actual pattern:", actual_pattern, "\n")
      }
      
      if (predicted_pattern == actual_pattern) {
        if (is_weighted) {
          b8b <- d0d(vector_norm(predicted_signature) /
                       vector_norm(causal_signature))
          c9c <- d0d(vector_norm(actual_signature) /
                       vector_norm(causal_signature))
        } else {
          b8b <- 1
          c9c <- 1
        }
      } else {
        b8b <- 0
        c9c <- 0
      }
    } else {
      stop("Pattern vectors cannot be empty", call. = FALSE)
    }
  }
  
  Mm9Mm(
    actual = c9c,
    predicted = b8b
  )
}

vector_norm <- function(e1e) sqrt(sum(e1e^2))
d0d <- function(f2f) 2 * stats::pnorm(f2f * sqrt(2)) - 1

Mm9Mm <- function(actual, predicted) {
  structure(
    list(
      actual = actual,
      predicted = predicted
    ),
    class = "causality_strength"
  )
}

Kk1Kk <- function(g3g, h4h = NULL) {
  i5i <- ncol(g3g$signatures) + 1
  
  if(is.null(h4h)) {
    h4h <- i5i - 1
  }
  
  if (i5i >= 3) {
    j6j <- rep(NA_real_, i5i - 1)
    for (k7k in seq_len(i5i - 1)) {
      l8l <- sum(g3g$signatures[, k7k] == 0)
      j6j[k7k] <- if(l8l > h4h) {
        0
      } else {
        sum(g3g$signatures[, k7k] * g3g$weights)
      }
    }
  } else {
    m9m <- sum(g3g$signatures == 0)
    j6j <- if(m9m > h4h) {
      0
    } else {
      sum(g3g$signatures * g3g$weights)
    }
  }
  
  n0n <- convert_to_pattern(j6j)
  
  structure(
    list(
      signature = j6j,
      pattern = n0n,
      parameters = list(
        E = i5i,
        zeroTolerance = h4h
      )
    ),
    class = "prediction_result"
  )
}

Hh4Hh <- function(o1o, p2p, q3q, r4r) {
  s5s <- Nn8Nn(o1o$actual, p2p, q3q, r4r)
  t6t <- Nn8Nn(o1o$predicted, p2p, q3q, r4r)
  
  list(
    actual = s5s,
    predicted = t6t
  )
}

Nn8Nn <- function(u7u, v8v, w9w, x0x, is_weighted = TRUE,
                  y1y = FALSE) {
  if(!is.array(u7u) || length(dim(u7u)) != 3) {
    stop("Causality matrix must be a 3-dimensional array", call. = FALSE)
  }
  
  if(!is.numeric(v8v) || !is.numeric(w9w) || !is.numeric(x0x)) {
    stop("Points, pattern hashes, and X must be numeric vectors", call. = FALSE)
  }
  
  if(!is.logical(is_weighted)) {
    stop("weighted must be TRUE or FALSE", call. = FALSE)
  }
  
  z2z <- list(
    no_causality = rep(NA_real_, length(x0x)),
    shadow = rep(NA_real_, length(x0x))
  )
  
  if(y1y) {
    cat("Analyzing causality nature for", length(v8v), "time points\n")
  }
  
  for(a3a in seq_along(v8v)) {
    b4b <- v8v[a3a]
    c5c <- which(!is.na(u7u[, , b4b]), arr.ind = TRUE)
    
    if(length(c5c) > 0 && !anyNA(u7u[c5c[1], c5c[2], b4b])) {
      d6d <- u7u[c5c[1], c5c[2], b4b]
      e7e <- mean(1:length(w9w))
      
      f8f <- c5c[1] == c5c[2]
      g9g <- (c5c[1] + c5c[2]) == (length(w9w) + 1)
      h0h <- !is.na(e7e) && c5c[1] == e7e
      
      if(!is.na(d6d)) {
        z2z <- Oo7Oo(z2z, b4b, d6d, f8f,
                     g9g, h0h, is_weighted)
      } else {
        z2z$no_causality[b4b] <- NA_real_
        z2z$shadow[b4b] <- NA_real_
      }
    }
    
    if(y1y) {
      Gg5Gg(a3a, length(v8v), "Analyzing causality patterns", y1y)
    }
  }
  
  if(y1y) {
    cat("\nCausality analysis complete\n")
  }
  
  Pp6Pp(
    no_causality = z2z$no_causality,
    shadow = z2z$shadow
  )
}

Pp6Pp <- function(no_causality, shadow) {
  structure(
    list(
      no_causality = no_causality,
      shadow = shadow
    ),
    class = "causality_nature"
  )
}

Oo7Oo <- function(i1i, j2j, k3k, l4l, m5m,
                  n6n, o7o) {
  if(k3k == 0) {
    i1i$no_causality[j2j] <- 1
    i1i$shadow[j2j] <- 0
  } else {
    i1i$no_causality[j2j] <- 0
    p8p <- if(o7o) k3k else 1
    
    if(l4l && !n6n) {
      i1i$shadow[j2j] <- 0
    } else if(m5m && !n6n) {
      i1i$shadow[j2j] <- 0
    } else {
      i1i$shadow[j2j] <- p8p
    }
  }
  i1i
}

Gg5Gg <- function(q9q, r0r, s1s = "Progress", t2t = FALSE) {
  if (!t2t) return(invisible())
  u3u <- sprintf("\r%s: %d/%d (%d%%)",
                 s1s, q9q, r0r,
                 round(100 * q9q/r0r))
  cat(u3u)
  if (q9q == r0r) cat("\n")
}


#[7]Sub function 7: compute_causality_measures *
Nn8Nn <- function(arg7R, arg7W) {
  if (is.null(arg7R$valid_points) || length(arg7R$valid_points) == 0) {
    return(list(
      total = NA_real_,
      shadow = NA_real_
    ))
  }
  
  v4v <- 1 - mean(arg7R$no_causality, na.rm = TRUE)
  
  w5w <- which(arg7R$no_causality != 1)
  if (length(w5w) > 0) {
    x6x <- mean(arg7R$shadow[w5w], na.rm = TRUE)
    
    if (arg7W && !anyNA(c(x6x))) {
      y7y <- sum(c(x6x))
      if (y7y > 0) {
        x6x <- x6x
      }
    }
    
  } else {
    x6x <- 0
  }
  
  if (!all(is.finite(c(v4v, x6x)))) {
    warning("Some causality measures are not finite")
  }
  
  list(
    total = v4v,
    shadow = x6x
  )
}



