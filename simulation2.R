source("np_simex.R")

estimator <- function(X) { sum(X^4)/length(X) }
est.fun <- function(X, y.filter=NULL) { sum(X^4)/length(X) }

ns <- c(100, 500, 1000, 5000, 10000, 20000, 50000, 100000)
runs <-  1000
B <- 500
M <- 10
models <- list(Y~a + b*X + d*X^2)
starts <- list(list(a = 1, b = 1, d = 1))
simulation_2 <- list()

for(n in ns){
    simulation_2[[paste0("n_", n)]] <- list()
    for (ii in 1:runs) {
        set.seed(314 + ii + 2000 + n)
        cat(paste0("\r Simulation 2 (", n,"):\t", ii, " / ", runs))
        try({

            X <- rnorm(n, mean = 5, sd = 2)

            filter1 <- rbinom(n, size = 1, prob = 1-mix_parm)
            X1 <- X + rt(n, df=5)
            X2 <- X + rt(n, df=5)
            replicates <- list(X1, X2)

            # Regular SIMEX
            sigma.hat <- sqrt(1/(2*n)*sum((X1-X2)^2))
            standard_pred <-simex(replicates, Sigma=sigma.hat, est.fun=est.fun, B=B, M=M, models=models, starts=starts, seed=314)

            simex_df <- np_simex(replicates, estimator)
            fitted_model <- nls(theta ~ a + b*lambda + c*lambda^2, data = simex_df, start = list(a = 1, b = 1, c = 1))
            predicted_value <- predict(fitted_model, newdata=data.frame(lambda = -1))

            
            smoothed_simex_df <- np_simex(replicates, estimator, smoothed = TRUE)
            smoothed_fitted_model <- nls(theta ~ a + b*lambda + c*lambda^2, data = smoothed_simex_df, start = list(a = 1, b = 1, c = 1))
            smoothed_predicted_value <- predict(smoothed_fitted_model, newdata=data.frame(lambda = -1))

            simulation_2[[paste0("mp_", n)]][[ii]] <- list("NP-SIMEX" = list(simex_df = simex_df, pred_theta = predicted_value), 
                                                           "NP-SIMEX-SMOOTHED" = list(simex_df = smoothed_simex_df, pred_theta = smoothed_predicted_value), 
                                                           "P-SIMEX" = list(pred_theta = standard_pred))
        })
    }
    cat("\n")
    save(simulation_2, file = "simulation2.RDa")
}

save(simulation_2, file = "simulation2.RDa")