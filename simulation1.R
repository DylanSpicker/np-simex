library(foreach)
library(doParallel)
source("np_simex.R")

dfs <- c(3, 4, 5, 10, 30)           # Something that contains all of the scenarios to loop through
replicate_start <- 1                # Number of times to repeat the loop
replicate_end <- 200               
bs_replicates <- 500                # Number of times to bootstrap the estimator
n <- 5000                           # Sample Size

### NP SIMEX/SIMEX options
B <- 100
M <- 10

cluster <- parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cluster)

set.seed(3141592) 

replicate_seeds <- sample(1000:100000, replicate_end, FALSE)
bootstrap_seeds <- sample(1000:100000, bs_replicates, FALSE)

estimator <- function(X, Y) { 
    MX <- model.matrix(~X)
    M <- fastglm::fastglm(y=Y, x=MX, family='binomial')
    M$coefficients[2]
}

t1 <- Sys.time()
loop_list <- foreach(df = dfs) %:%
                foreach(ridx = replicate_start:replicate_end) %:%
                foreach(bidx = 1:bs_replicates) %dopar% {
    ## Set the Seed to the Replicate # to Generate the Standard Data
    set.seed(replicate_seeds[ridx]) 
    X <- rnorm(n, mean = 1, sd = 2)
    Y <- rbinom(n, size=1, prob=1/(1+exp(-1+X)))
    X1 <- X + rt(n, df)
    X2 <- X + rt(n, df)

    ## Set the Seed Based on the Bootstrap Replicate
    ## Offset by ridx so that bs_rep 1 for replicate 1 is not the same as 
    ### bs_rep 1 for replicate 2, etc.
    set.seed(replicate_seeds[ridx] + bootstrap_seeds[bidx])

    which_idx <- sample(1:n, n, TRUE)
    Ybs <- Y[which_idx]
    X1bs <- X1[which_idx]
    X2bs <- X2[which_idx]
    replicates_bs <- list(X1bs, X2bs)
    Xsbs <- (X1bs + X2bs)/2

    myl <- tryCatch({
        sigma.hat_bs <- sqrt(1/(2*n)*sum((X1bs-X2bs)^2))
        naive_model_bs <- glm(Ybs ~ Xsbs, family=binomial)
        sm1_bs <- simex::simex(naive_model_bs, 
                               SIMEXvariable="Xsbs", 
                               measurement.error = sigma.hat_bs, 
                               lambda = seq(0,2,length.out=M)[2:M], 
                               B=B, 
                               jackknife.estimation=FALSE, 
                               asymptotic=FALSE)

        #### NP-SIMEX
        simex_df_bs <- np_simex(replicates_bs, estimator, Y=Ybs, parallel=FALSE, M=M, B=B)
        fitted_model_bs <- nls(theta ~ a + b/(c + lambda), data = simex_df_bs, start = list(a = 1, b = 1, c = 1))
        predicted_value_bs <- predict(fitted_model_bs, newdata=data.frame(lambda = -1))
        v1 <- sm1_bs$coefficients[2]
        v2 <- predicted_value_bs
        c(v1, v2)
    },
    error = function(e) {
        return(c(NA, NA))
    })

    myl
}
t2 <- Sys.time()
parallel::stopCluster(cluster)

save(
    loop_list, 
    file=paste0("simulation1.RDa")
)
