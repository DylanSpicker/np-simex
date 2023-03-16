library(parallel)
library(foreach)
library(doParallel)
source("np_simex.R")

t1 <- Sys.time()
no_cores <- detectCores() - 1  
registerDoParallel(cores=no_cores)
cl <- makeCluster(no_cores)  
ns <- c(1000, 10000, 100000)
percs <- c(0.05, 0.10, 0.5)
sigma_ratios <- c(0.1, 0.5, 1, 2)
reps <- 1000
B <- 200
M <- 10
results <- foreach (n = ns) %:%
            foreach(perc = percs) %:%
            foreach(sr = sigma_ratios) %:%
            foreach(r = 1:reps) %dopar% {
    set.seed(10*n + 1000*perc + r)
    X <- rgamma(n, shape = 2) # The standard deviation is given by sqrt(2)
    
    # The standard deviation is given by sqrt(2) * ratio (sr)
    U.full <- rgamma(n, shape = 1, scale=sr) - rgamma(n, shape = 1, scale = sr)
    W <- X + U.full
    sample <- sample(1:n, floor(perc*n), FALSE)
    Z <- rnorm(n)
    Y <- rbinom(n, 1, 1/(1+exp(-1*(1 + Z - 1.25*X))))
    U <- U.full[sample]

    sigma.hat <- var(U)
    naive_model <- glm(Y ~ W + Z, family='binomial')
    try({
        sm1 <- simex::simex(naive_model, 
                                SIMEXvariable='W', 
                                measurement.error=sigma.hat, 
                                B=B, 
                                jackknife.estimation=FALSE,
                                asymptotic=FALSE,
                                lambda=seq(0,2,length.out=M)[2:M])
        simex_df <- np_simex(W, function(X){ coef(glm(Y ~ X + Z, family='binomial')) }, U=U, M=M, B=B, parallel = FALSE)
        smoothed_simex_df <- np_simex(W, function(X){ coef(glm(Y ~ X + Z, family='binomial')) }, U=U, M=M, B=B, parallel = FALSE, smoothed = TRUE)
        
        h.0 <- predict(nls(theta.1 ~ a + b/(c + lambda), data = simex_df, start = list(a=1, b=1, c=1), nls.control(maxiter = 5000)), newdata=data.frame(lambda=-1))
        h.1 <- predict(nls(theta.2 ~ a + b/(c + lambda), data = simex_df, start = list(a=1, b=1, c=1), nls.control(maxiter = 5000)), newdata=data.frame(lambda=-1))
        h.2 <- predict(nls(theta.3 ~ a + b/(c + lambda), data = simex_df, start = list(a=1, b=1, c=1), nls.control(maxiter = 5000)), newdata=data.frame(lambda=-1))

        sh.0 <- predict(nls(theta.1 ~ a + b/(c + lambda), data = smoothed_simex_df, start = list(a=1, b=1, c=1), nls.control(maxiter = 5000)), newdata=data.frame(lambda=-1))
        sh.1 <- predict(nls(theta.2 ~ a + b/(c + lambda), data = smoothed_simex_df, start = list(a=1, b=1, c=1), nls.control(maxiter = 5000)), newdata=data.frame(lambda=-1))
        sh.2 <- predict(nls(theta.3 ~ a + b/(c + lambda), data = smoothed_simex_df, start = list(a=1, b=1, c=1), nls.control(maxiter = 5000)), newdata=data.frame(lambda=-1))

        list(
            "Naive" = coef(naive_model),
            "P-SIMEX" = sm1$coefficients,
            "NP-SIMEX" = c(h.0, h.1, h.2),
            "NP-SIMEX-SMOOTH" = c(sh.0, sh.1, sh.2),
            "setting" = list("r" = r, "sr" = sr, "n" = n, "perc" = perc)
        )
    })
}
stopCluster(cl)
t2 <- Sys.time()

results_df <- data.frame(method=NULL, beta0=NULL, beta1=NULL, beta2=NULL, n=NULL, percentage=NULL, sigma_ratio=NULL)

i<-1

for(n in ns) {
    j <- 1
    for(perc in percs) { 
        k <- 1
        for(sr in sigma_ratios) {
            l <- 1
            for (r in 1:reps) { 
                al <- results[[i]][[j]][[k]][[l]]
                if (class(al) != "try-error") {
                    results_df <- rbind(results_df, 
                                    data.frame(method="Naive", beta0=al$Naive[1], beta1=al$Naive[2], beta2=al$Naive[3], n=n, percentage=perc, sigma_ratio=sr),
                                    data.frame(method="P-SIMEX", beta0=al$`P-SIMEX`[1], beta1=al$`P-SIMEX`[2], beta2=al$`P-SIMEX`[3], n=n, percentage=perc, sigma_ratio=sr),
                                    data.frame(method="NP-SIMEX", beta0=al$`NP-SIMEX`[1], beta1=al$`NP-SIMEX`[2], beta2=al$`NP-SIMEX`[3], n=n, percentage=perc, sigma_ratio=sr),
                                    data.frame(method="NP-SIMEX-SMOOTH", beta0=al$`NP-SIMEX-SMOOTH`[1], beta1=al$`NP-SIMEX-SMOOTH`[2], beta2=al$`NP-SIMEX-SMOOTH`[3], n=n, percentage=perc, sigma_ratio=sr)
                                    )
                }
                l <- l + 1
            }
            k <- k + 1
        }
        j <- j + 1
    }
    i <- i + 1
}

rownames(results_df) <- NULL
save(results_df, file="validation_set_results.RDa")