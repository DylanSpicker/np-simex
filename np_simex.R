library(np)

# Helper Function
# Permutations to Generate the Perms of a certain length
permutations <- function(n) {
  if(n == 1) {
    return(matrix(1))
  } else {
    sp <- permutations(n - 1)
    p <- nrow(sp)
    A <- matrix(nrow = n * p, ncol = n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

# Sample from a conditional KDE, with specified bandwidth parameters
sample_conditional_errors <- function(w, lambda, xbw, ybw, W, U) {
    weights <- dnorm(x = w - W, mean = 0, sd = xbw)
    w_probs <- weights/sum(weights)

    rnorm(n = lambda,
          sample(U, 
                 size = lambda,
                 replace = TRUE,
                 prob = w_probs),
          ybw)
}

np_simex <- function(replicates,
                     estimator,
                     U = NULL,
                     val_Xs = NULL,
                     M = 10,
                     B = 50,
                     parallel = TRUE,
                     numCores = parallel::detectCores()/2,
                     est.variance = "none",
                     parPackage = "foreach",
                     smoothed = FALSE,
                     het = FALSE,
                     ...) {
  if(class(estimator) != "function") stop("You must provide your estimator.")
  if(!("X" %in% formalArgs(estimator))) stop("Your estimator must accept the 'X' parameter.")

  if(het && is.null(X)) { stop("To conduct a heteroscedastic inference, 'X' must be provided, alongside 'U'.") }

  if (parPackage == "foreach") {
    library(doMC)
    registerDoMC(numCores)
  }

  if(class(replicates) != "list") {
    ## This method we provide U and replicates is just a numeric vector;
    ## this should work, e.g. with an internal validation sample
    if(is.null(U)) stop("You provided a vector of covariates, but then no set of U. That doesn't work.") 
    Xs <- replicates
    n <- length(Xs)
  } else {
    k <- length(replicates)
    n <- nrow(as.matrix(replicates[[1]]))
    

    if(length(replicates) %% 2 != 0) {
      factor_v <- c(rep(1/(k+1), (k+1)/2), rep(1/(k-1), (k-1)/2))
      Xs <- Reduce("+", lapply(1:k, function(xx){
        replicates[[xx]]*factor_v[xx]
      }))
      contrast_standard <- c(rep(1/(k+1), (k+1)/2), rep(-1/(k-1), (k-1)/2))
      
    } else {
      Xs <- (1/k)*Reduce("+", replicates)
      contrast_standard <- c(rep(1/k, k/2), rep(-1/k, k/2))
    }
    
    
    contrast_matrix <- permutations(k)
    contrast_matrix <- unique(t(apply(contrast_matrix, MARGIN = 1, FUN = function(row) {
      contrast_standard[row]
    })))
    
    
    
    U <- c(apply(contrast_matrix, MARGIN = 1, FUN = function(contrast) {
      us <- 0
      for (ii in 1:k) {
        us <- us + contrast[ii]*replicates[[ii]]
      }
      us 
    }))

  }
  
  ## Compute Standard Estimator
  if (est.variance == "jackknife") {
    e1 <- estimator(X=Xs, ...)
    theta.hat <- list(e1$theta)
    sigma.hat <- list(e1$sigma)
    v.hat <- list(0)
  } else {
    theta.hat <- list(estimator(X=Xs, ...))
  }

  if(smoothed) {
    density_U <- density(U, bw = "SJ")
    h_U <- density_U$bw
  } else if (het) {
    cbw <- npcdensbw(U ~ val_Xs)
    h_U <- cbw$ybw
    h_Xs <- cbw$xbw

    all_pseudo_errors <- t(sapply(Xs, sample_conditional_errors, lambda = M*(M+1)*B/2, xbw = h_Xs, ybw = h_U, W=val_Xs, U=U))
  }

  ## Simulation Stage
  if (! parallel) {
    for (lambda in 1:M) {
      theta.hat[[lambda+1]] <- 0
      if (est.variance == "jackknife") {
        sigma.hat[[lambda+1]] <- 0
      }
      col_start <- B*lambda*(lambda - 1)/2 + 1

      for(b in 1:B) {
        if(smoothed) {
          pseudo_error <- replicate(rnorm(n, sample(U, size = n, replace=TRUE), h_U), n=lambda)
        } else if(het) {
          pseudo_error <- as.matrix(all_pseudo_errors[,col_start:(col_start+lambda-1)])
          col_start <- col_start + lambda
        } else {
          pseudo_error <- replicate(sample(U, size = n, replace=TRUE), n=lambda)
        }
        

        Xs.lambda <- Xs + rowSums(pseudo_error)
        dEst <- estimator(X=Xs.lambda, ...)
        if(est.variance == "jackknife") {
          theta.hat[[lambda+1]] <- theta.hat[[lambda+1]] + dEst$theta
          sigma.hat[[lambda+1]] <- sigma.hat[[lambda+1]] + dEst$sigma
        } else {
          theta.hat[[lambda+1]] <- theta.hat[[lambda+1]] + dEst
        }
      }
      
      if (est.variance == "jackknife") {
        sigma.hat[[lambda+1]] <-(1/B)*sigma.hat[[lambda+1]]
        v.hat[[lambda+1]] <- var(unlist(theta.hat))
      }
      theta.hat[[lambda+1]] <- (1/B)*theta.hat[[lambda+1]]
    }
  } else {
    if (parPackage == "foreach") {
      ### USE FOREACH
      if(est.variance == "jackknife") {
        allObjs <- foreach::foreach(lambda = 1:M) %dopar% {
            col_start <- B*lambda*(lambda - 1)/2 + 1
            objList <- foreach::foreach(b = 1:B) %dopar% {
              
              if(smoothed) {
                pseudo_error <- replicate(rnorm(n, sample(U, size = n, replace=TRUE), h_U), n=lambda)
              } else if(het) { 
                pseudo_error <- as.matrix(all_pseudo_errors[,(col_start+(b-1)*lambda):(col_start+(b)*lambda-1)])
              } else {
                pseudo_error <- replicate(sample(U, size = n, replace=TRUE), n=lambda)
              }

              Xs.lambda <- Xs + rowSums(pseudo_error)
              estimator(X=Xs.lambda, ...)
            }

            thetas <- c()
            sigma.hat.ind <- 0

            for(ob in objList) {
              thetas <- c(thetas, ob$theta)
              sigma.hat.ind <- sigma.hat.ind + ob$sigma
            }

            list(theta = mean(thetas), sigma = sigma.hat.ind / B, v = var(thetas))
          }


        for (idx in 2:(M+1)) {
          theta.hat[[idx]] <- allObjs[[idx-1]]$theta
          sigma.hat[[idx]] <- allObjs[[idx-1]]$sigma
          v.hat[[idx]] <- allObjs[[idx-1]]$v
        }

      } else {
        theta.hat[2:(M+1)] <- foreach::foreach(lambda = 1:M) %dopar% {
            col_start <- B*lambda*(lambda - 1)/2 + 1
            (1/B)*Reduce("+", foreach::foreach(b = 1:B) %dopar% {

              if(smoothed) {
                pseudo_error <- replicate(rnorm(n, sample(U, size = n, replace=TRUE), h_U), n=lambda)
              } else if(het) { 
                  pseudo_error <- as.matrix(all_pseudo_errors[,(col_start+(b-1)*lambda):(col_start+(b)*lambda-1)])
              } else {
                pseudo_error <- replicate(sample(U, size = n, replace=TRUE), n=lambda)
              }

              Xs.lambda <- Xs + rowSums(pseudo_error)
              estimator(X=Xs.lambda, ...)
            })
          }
      }
    } else {
      ### USE MCLAPPLY
      if(est.variance == "jackknife") {
        allObjs <- parallel::mclapply(
          1:M,
          function(lambda){
            col_start <- B*lambda*(lambda - 1)/2 + 1
            objList <- parallel::mclapply(1:B, function(b){

              if(smoothed) {
                pseudo_error <- replicate(rnorm(n, sample(U, size = n, replace=TRUE), h_U), n=lambda)
              } else if(het) { 
                pseudo_error <- as.matrix(all_pseudo_errors[,(col_start+(b-1)*lambda):(col_start+(b)*lambda-1)])
              } else {
                pseudo_error <- replicate(sample(U, size = n, replace=TRUE), n=lambda)
              }

              Xs.lambda <- Xs + rowSums(pseudo_error)
              estimator(X=Xs.lambda, ...)
            }, mc.cores = numCores)

            thetas <- c()
            sigma.hat.ind <- 0

            for(ob in objList) {
              thetas <- c(thetas, ob$theta)
              sigma.hat.ind <- sigma.hat.ind + ob$sigma
            }

            list(theta = mean(thetas), sigma = sigma.hat.ind / B, v = var(thetas))
          },
          mc.cores = numCores
        )


        for (idx in 2:(M+1)) {
          theta.hat[[idx]] <- allObjs[[idx-1]]$theta
          sigma.hat[[idx]] <- allObjs[[idx-1]]$sigma
          v.hat[[idx]] <- allObjs[[idx-1]]$v
        }

      } else {
        theta.hat[2:(M+1)] <- parallel::mclapply(
          1:M,
          function(lambda){
            col_start <- B*lambda*(lambda - 1)/2 + 1
            (1/B)*Reduce("+", parallel::mclapply(1:B, function(b){
              if(smoothed) {
                pseudo_error <- replicate(rnorm(n, sample(U, size = n, replace=TRUE), h_U), n=lambda)
              } else if(het) { 
                pseudo_error <- as.matrix(all_pseudo_errors[,(col_start+(b-1)*lambda):(col_start+(b)*lambda-1)])
              } else {
                pseudo_error <- replicate(sample(U, size = n, replace=TRUE), n=lambda)
              }

              Xs.lambda <- Xs + rowSums(pseudo_error)
              estimator(X=Xs.lambda, ...)
            }, mc.cores = numCores))
          },
          mc.cores = numCores
        )
      }
    }
  }
  
  if(est.variance == "jackknife") {
    df <- data.frame(lambda = 0:M, 
               theta = matrix(unlist(theta.hat), nrow=(M+1), byrow=T), 
               sigma = matrix(unlist(sigma.hat), nrow=(M+1), byrow=T), 
               v = matrix(unlist(v.hat), nrow=(M+1), byrow=T))
    df$variance_extrapolant <- df$sigma - df$v
    df
  } else {
    data.frame(lambda = 0:M, theta = matrix(unlist(theta.hat), nrow=(M+1), byrow=T))
  }
}