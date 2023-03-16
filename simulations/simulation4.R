library(np)
library(fastglm)
library(simex)
library(doMC)

registerDoMC(8)

## Helper Functions
sample_conditional <- function(x, n, xbw, ybw, VX, VY) {
    weights <- dnorm(x = x - VX, mean = 0, sd = xbw)

    if(sum(weights) == 0) {
        w_probs <- rep(1/length(weights), length(weights))
    } else {
        w_probs <- weights/sum(weights)
    }

    rnorm(n = n, sample(VY, size = n, replace = TRUE, prob = w_probs), ybw)
}
expit <- function(w){ 1/(1+exp(-w))}

################################################################################
# Parameters
################################################################################
M <- 10
B <- 500

n <- 1000
prop <- 0.2
muX <- 1
sdX <- 2

relation_factors <- c(0, 0.5, 1, 2)

replicates <- 500

seedoffset <- 31415
################################################################################
# Single Iteration
################################################################################

iteration_replicates <- function(rf) {
    X <- rnorm(n, mean = muX, sd = sdX)
    Y <- rbinom(n, size=1, prob=expit(1-X))
    Xs <- X + rnorm(n, (X-muX)*rf, 1)
    sample_idx <- sample(1:n, n*prop, FALSE)
    val_X <- X[sample_idx]
    val_Xs <- Xs[sample_idx]
    val_U <- val_Xs - val_X

    # Let's do each analysis
    estimator <- function(X, Y) { 
        MX <- model.matrix(~X)
        M <- fastglm(y=Y, x=MX, family='binomial')
        M$coefficients[2]
    }

    ## Find the bandwidth parameters for the three conditional distributions
    f_UX <- npcdensbw(val_U ~ val_X)
    f_UXs <- npcdensbw(val_U ~ val_Xs)
    f_XXs <- npcdensbw(val_X ~ val_Xs)

    # Estimator 1: Drawing from X|X*, and then averaging 'B' times
    method_1_results <- 0
    all_X <- t(sapply(X=Xs, 
                      FUN=sample_conditional, 
                      n = B, 
                      xbw = f_XXs$xbw, 
                      ybw = f_XXs$ybw, 
                      VX = val_Xs, 
                      VY = val_X))

    for(b in 1:B) {
        method_1_results <- method_1_results + estimator(all_X[,b], Y)
    }
    method_1_results <- (1/B)*method_1_results

    # Estimator 2: Drawing from X|X*, and then using this to draw from U|X
    method_2.theta.hat <- list(estimator(Xs, Y))

    for (lambda in 1:M) {
        method_2.theta.hat[[lambda+1]] <- 0
        
        for(b in 1:B) {
            if(lambda == 1) {
                pseudo_error <- as.matrix(sapply(all_X[,b], sample_conditional, n = lambda, xbw = f_UX$xbw, ybw = f_UX$ybw, VX = val_X, VY = val_U))
            } else {
                pseudo_error <- t(sapply(all_X[,b], sample_conditional, n = lambda, xbw = f_UX$xbw, ybw = f_UX$ybw, VX = val_X, VY = val_U))
            }
            Xs.lambda <- Xs + rowSums(pseudo_error)
            method_2.theta.hat[[lambda+1]] <- method_2.theta.hat[[lambda+1]] + estimator(X=Xs.lambda, Y)
        }
        
        method_2.theta.hat[[lambda+1]] <- (1/B)*method_2.theta.hat[[lambda+1]]
    }
    method_2_df <- data.frame(lambda = 0:M, theta = matrix(unlist( method_2.theta.hat), nrow=(M+1), byrow=T))
    fitted_model_nl_2 <- nls(theta ~ a + b/(c + lambda), data = method_2_df, start = list(a = 1, b = 1, c = 1))
    predicted_value_nl_2 <- predict(fitted_model_nl_2, newdata=data.frame(lambda = -1))

    # Estimator 3: Drawing from U|X* directly
    all_pseudo_errors <- t(sapply(Xs, sample_conditional, n = M*(M+1)*B/2, xbw = f_UX$xbw, ybw = f_UX$ybw, VX = val_Xs, VY = val_U))
    method_3.theta.hat <- list(estimator(Xs, Y))

    for (lambda in 1:M) {
        method_3.theta.hat[[lambda+1]] <- 0
        
        col_start <- B*lambda*(lambda - 1)/2 + 1
        
        for(b in 1:B) {
            pseudo_error <- as.matrix(all_pseudo_errors[,col_start:(col_start+lambda-1)])
            Xs.lambda <- Xs + rowSums(pseudo_error)
            method_3.theta.hat[[lambda+1]] <- method_3.theta.hat[[lambda+1]] + estimator(X=Xs.lambda, Y)
        }
        
        method_3.theta.hat[[lambda+1]] <- (1/B)*method_3.theta.hat[[lambda+1]]
    }
    method_3_df <- data.frame(lambda = 0:M, theta = matrix(unlist(method_3.theta.hat), nrow=(M+1), byrow=T))
    fitted_model_nl_3 <- nls(theta ~ a + b/(c + lambda), data = method_3_df, start = list(a = 1, b = 1, c = 1))
    predicted_value_nl_3 <- predict(fitted_model_nl_3, newdata=data.frame(lambda = -1))

    # Estimator 4: Standard SIMEX
    naive_model <- glm(Y ~ Xs, family=binomial)
    sm1 <- simex(naive_model, SIMEXvariable="Xs", measurement.error = sqrt(var(val_U)), lambda = seq(0,2,length.out=M)[2:M], B=B, jackknife.estimation=FALSE, asymptotic=FALSE)

    c(method_1_results, sm1$coefficients[2], predicted_value_nl_2, predicted_value_nl_3)
}

### Run the Simulations
simulation4_results <- foreach::foreach(idx = 1:replicates) %dopar% {
    set.seed(seedoffset + idx)
    try({ lapply(relation_factors, function(rf) { iteration_replicates(rf) }) })
}

save(simulation4_results, file = "simulation4.Rda")
