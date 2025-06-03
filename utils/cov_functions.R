
tridiag <- function(upper, lower, main){
  out <- matrix(0,length(main),length(main))
  diag(out) <- main
  indx <- seq.int(length(upper))
  out[cbind(indx+1,indx)] <- lower
  out[cbind(indx,indx+1)] <- upper
  return(out)
}

mat.ma2 <- function(MAparameters, order){
  if (length(MAparameters) < 2) 
    stop("Must supply two values in MAparameters")
  if (abs(MAparameters[1]) >= 1) 
    stop("MAparameters[1] must be less than 1")
  if (abs(MAparameters[2]) >= 1) 
    stop("MAparameters[2] must be less than 1")
  if ((MAparameters[1] + MAparameters[2] >= 1) | (MAparameters[1] - 
                                                  MAparameters[2] >= 1)) 
    stop("The sum and difference of MAparameters must be less than 1")
  if (order < 3) 
    stop("The order must be 3 or more")
  div <- 1 + MAparameters[1] * MAparameters[1] + MAparameters[2] * 
    MAparameters[2]
  rho1 <- -MAparameters[1] * (1 - MAparameters[2])/div
  rho2 <- -MAparameters[2]/div
  ma2 <- mat.banded(x = c(1, rho1, rho2), nrow = order, ncol = order)
  return(ma2)
}


mat.banded <- function(x, nrow, ncol){
  nband <- length(x)
  if (nband > min(nrow, ncol)) 
    stop("Have supplied values for more than ", min(nrow, 
                                                    ncol), " bands")
  matrix <- matrix(0, nrow = nrow, ncol = ncol)
  for (i in 1:nband) matrix[row(matrix) == col(matrix) + i - 
                              1 | row(matrix) + i - 1 == col(matrix)] <- x[i]
  return(matrix)
}

#------- Covariance structures ####

get_cor_matrix_ar1 <- function (ar, nobs) {
  out <- array(0, dim = c(NROW(ar), nobs, nobs))
  fac <- 1/(1 - ar^2)
  pow_ar <- as.list(rep(1, nobs + 1))
  for (i in seq_len(nobs)) {
    pow_ar[[i + 1]] <- ar^i
    out[, i, i] <- fac
    for (j in seq_len(i - 1)) {
      out[, i, j] <- fac * pow_ar[[i - j + 1]]
      out[, j, i] <- out[, i, j]
    }
  }
  out
}

autocorr.mat <- function(p = 100, rho = 0.9) {
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}

get_sigmaX_ar1 <- function(D, rho){
  cor=get_cor_matrix_ar1(ar=rho, nobs=D)
  temp=matrix(0,nrow=D, ncol=D)
  for(i in 1:D){
    temp[i,]=cor[,,i]
  }
  temp
}


cor2cov <- function(R,S){
  # Given a correlation matrix R and a vector of standard deviations S, find the
  # corresponding covariance matrix
  sweep(sweep(R,1,S,"*"), 2, S, "*") #fast
  #diag(S) %*% R %*% diag(S) #slow
}

gen_covmatrix <- function(params) {
  # Validate required parameters
  if (!is.list(params)) {
    stop("Parameter 'params' must be a list.")
  }
  
  if (is.null(params$p) || !is.numeric(params$p) || params$p <= 0) {
    stop("Parameter 'p' must be a positive integer.")
  }
  
  if (is.null(params$type)) {
    stop("Parameter 'type' must be a string specifying the covariance structure.")
  }
  
  p <- as.integer(params$p)
  type <- tolower(params$type)
  
  switch(type,
         "ar" = {
           if (is.null(params$rho)) stop("Parameter 'rho' must be provided for AR structure.")
           return(generate_ar_covariance(p, params$rho))
         },
         "toe" = {
           if (is.null(params$toe)) stop("Parameter 'toe' must be provided for TOE structure.")
           return(generate_toeplitz_covariance(p, params$toe))
         },
         "ic" = {
           if (is.null(params$rho)) stop("Parameter 'rho' must be provided for IC structure.")
           return(generate_ic_covariance(p, params$rho))
         },
         "ma1" = {
           if (is.null(params$rho)) stop("Parameter 'rho' must be provided for MA1 structure.")
           return(generate_ma1_covariance(p, params$rho))
         },
         "ma2" = {
           if (is.null(params$rho)) stop("Parameter 'rho' must be provided for MA2 structure.")
           return(generate_ma2_covariance(p, params$rho))
         },
         "blockedma1" = {
           if (is.null(params$rho) || is.null(params$block_size)) {
             stop("Parameters 'rho' and 'block_size' must be provided for Blocked MA1 structure.")
           }
           return(generate_blocked_ma1_covariance(p, params$rho))
         },
         "blockedma2" = {
           if (is.null(params$rho) || is.null(params$block_size)) {
             stop("Parameters 'rho' and 'block_size' must be provided for Blocked MA2 structure.")
           }
           return(generate_blocked_ma2_covariance(p, params$rho))
         },
         "blockedar" = {
           if (is.null(params$rho) || is.null(params$block_size)) {
             stop("Parameters 'rho' and 'block_size' must be provided for Blocked AR structure.")
           }
           return(generate_blocked_ar_covariance(p, params$rho))
         },
         "grouped" = {
           if (is.null(params$groups) || !is.list(params$groups)) {
             stop("Parameter 'groups' must be provided as a list for GROUPED structure.")
           }
           return(generate_grouped_covariance(p, params$groups))
         },
         {
           stop("Unknown covariance structure type provided.")
         }
  )
}

generate_ar_covariance <- function(p, rho) {
  if (!is.numeric(rho) || abs(rho) >= 1) {
    stop("Parameter 'rho' must be numeric and between -1 and 1 (exclusive) for AR structure.")
  }
  indices <- 0:(p - 1)
  sigma <- rho ^ abs(outer(indices, indices, "-"))
  return(sigma)
}

generate_toeplitz_covariance <- function(p, toe) {
  if (!is.numeric(toe) || length(toe) != p) {
    stop("Parameter 'toe' must be a numeric vector of length equal to 'p'.")
  }
  sigma <- toeplitz(toe)
  return(sigma)
}

generate_ic_covariance <- function(p, rho) {
  if (!is.numeric(rho) || rho <= -1/(p - 1) || rho >= 1) {
    stop("Parameter 'rho' must be between -1/(p-1) and 1 for IC structure.")
  }
  sigma <- matrix(rho, nrow = p, ncol = p)
  diag(sigma) <- 1
  return(sigma)
}


generate_ma1_covariance <- function(p, rho) {
  if (!is.numeric(rho)) {
    stop("Parameter 'rho' must be numeric for MA1 structure.")
  }
  sigma <- tridiag( upper= rep(-rho, p-1), lower= rep(-rho, p-1), main= rep(1+rho^2, p))
  return(sigma)
}

generate_ma2_covariance <- function(p, rho) {
  if (!is.numeric(rho)) {
    stop("Parameter 'rho' must be numeric for MA2 structure.")
  }
  
  a= rho
  b= (1-rho)*rho
  sigma= mat.ma2(c(a,b), order= p)

  return(sigma)
}

generate_blocked_ma1_covariance <- function(p, rho, block_size = 5, type = "BlockedMA1s") {
  
  if (p %% block_size != 0) {
    stop("Parameter 'p' must be divisible by 'block_size' for Blocked MA1 structure.")
  }
  
  if(type == "BlockedMA1ns"){
    # MA1 blocked at borders and in the middle independent blocks
    blockA= tridiag( upper= rep(-rho, block_size-1), lower= rep(-rho, block_size-1), main= rep(1+rho^2, block_size))
    blockB= diag(1, p-2*block_size)
    blockC= tridiag( upper= rep(-rho, block_size-1), lower= rep(-rho, block_size-1), main= rep(1+rho^2, block_size))
    sigma <- as.matrix(Matrix::bdiag(blockA, blockB, blockC))
    
  } else if(type== "BlockedMA1s"){
    # MA1 blocks 
    num_blocks <- p/ block_size  
    #block <- tridiag( upper= rep(-rho, block_size-1), lower= rep(-rho, block_size-1), main= rep(1+rho^2, block_size)) 
    block <- generate_ma1_covariance(p = block_size, rho = rho)
    sigma <- kronecker(diag(num_blocks), block) #repeat M in block
  }
  
  return(sigma)
}

generate_blocked_ma2_covariance <- function(p, rho, block_size = 5) {
  
  if (p %% block_size != 0) {
    stop("Parameter 'p' must be divisible by 'block_size' for Blocked MA2 structure.")
  }
  
  # MA2 blocks 
  num_blocks <- p/ block_size  
  block <- generate_ma2_covariance(p = block_size, rho =  rho)
  sigma <- kronecker(diag(num_blocks), block) #repeat M in block

  
  return(sigma)
}

generate_blocked_ar_covariance <- function(p, rho, block_size = 5) {
  if (p %% block_size != 0) {
    stop("Parameter 'p' must be divisible by 'block_size' for Blocked AR structure.")
  }
  num_blocks <- p / block_size
  block <- generate_ar_covariance(block_size, rho)
  sigma <- Matrix::bdiag(replicate(num_blocks, block, simplify = FALSE))
  sigma <- as.matrix(sigma)
  return(sigma)
}

generate_grouped_covariance <- function(p, groups) {
  total_size <- sum(sapply(groups, function(g) g$size))
  if (total_size != p) {
    stop("Sum of group sizes must equal 'p' for GROUPED structure.")
  }
  group_matrices <- lapply(groups, function(g) {
    generate_covariance_matrix(g)
  })
  sigma <- bdiag(group_matrices)
  sigma <- as.matrix(sigma)
  return(sigma)
}

