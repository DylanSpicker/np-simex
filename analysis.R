library(tidyverse)
library(haven)
library(fastglm)

mydf_main <- read.table(file="pre_generated/data_prim.txt")
mydf_val <- read.table(file="pre_generated/data_valid.txt")

colnames(mydf_main) <- c("Y", "W", "Z")
colnames(mydf_val) <- c("Y", "X", "W", "Z")

source("../np_simex.R")

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

Xs <- c(mydf_main$W, mydf_val$W)
val_Xs <- mydf_val$W
U <- mydf_val$W - mydf_val$X
X <- mydf_val$X
Z <- c(mydf_main$Z, mydf_val$Z)
Y <- c(mydf_main$Y, mydf_val$Y)

estimator <- function(X) { 
    MX <- model.matrix(~X+Z)
    M <- fastglm(y=Y, x=MX, family='binomial')
    M$coefficients
}

################################################################################
# GENERATE PLOT:
#
p <- mydf_val %>% 
        mutate(U = W - X) %>% 
        ggplot(aes(sample = U)) + 
        geom_qq_line() + 
        stat_qq() + 
        stat_qq_line() + 
        theme_minimal() + 
        xlab("Theoretical Quantiles") + 
        ylab("Observed Quantiles") + theme(text = element_text(size = 24))

ggsave(p, filename="qq_klosa.png", width = 12, height = 12, dpi = 600, units = "in", bg="transparent") 

p2 <- mydf_val %>% 
    mutate(U = W - X) %>% 
    ggplot(aes(x = X, y = U)) + 
    geom_smooth() +
    coord_cartesian(ylim = c(-10, 10)) + 
    geom_point() + 
    theme_minimal() + 
    xlab("X") + 
    ylab("U") + theme(text = element_text(size = 24))

ggsave(p2, filename="klosa_het.png", width = 12, height = 12, dpi = 600, units = "in", bg="transparent") 

#### PARAMETERS
M <- 10
B <- 200

##### STANDARD SIMEX
sigma.hat <- var(U)
naive_model <- glm(Y ~ Xs + Z , family=binomial)
sm1 <- simex::simex(naive_model, SIMEXvariable="Xs", measurement.error = sigma.hat, lambda = seq(0,2,length.out=M)[2:M], B=B, jackknife.estimation=FALSE, asymptotic=FALSE)

#### NP-SIMEX
simex_df <- np_simex(Xs, estimator, U=U, B=B, M=M, parallel=FALSE)

intercept_nls <- nls(theta.1 ~ a + b/(c + lambda), data = simex_df, start=list(a = 1, b = 1, c = 1))
beta1_nls <- nls(theta.2 ~ a + b/(c + lambda), data = simex_df, start=list(a = 1, b = 1, c = 1))
beta2_nls <- nls(theta.3 ~ a + b/(c + lambda), data = simex_df, start=list(a = 1, b = 1, c = 1))

intercept.hat <- predict(intercept_nls, data.frame(lambda = -1))
beta1.hat <- predict(beta1_nls, data.frame(lambda = -1))
beta2.hat <- predict(beta2_nls, data.frame(lambda = -1))

data.frame("Naive" = coef(naive_model), "P-SIMEX" = sm1$coefficients, "NP-SIMEX" = c(intercept.hat, beta1.hat, beta2.hat))

################################################################################
# CONDITIONAL ANALYSIS
################################################################################
f_UX <- npcdensbw(U ~ X)
f_UXs <- npcdensbw(U ~ val_Xs)
f_XXs <- npcdensbw(X ~ val_Xs)

# Method 2: 
all_X <- t(sapply(X=Xs, FUN=sample_conditional, n = B, xbw = f_XXs$xbw, ybw = f_XXs$ybw, VX = val_Xs, VY = X))
method_2.theta.hat <- list(estimator(Xs))

for (lambda in 1:M) {
    method_2.theta.hat[[lambda+1]] <- c(0, 0, 0)
    
    for(b in 1:B) {
        cat(paste0("\r Method 2: Lambda=", lambda, " and b=", b, "/", B))
        if(lambda == 1) {
            pseudo_error <- as.matrix(sapply(all_X[,b], sample_conditional, n = lambda, xbw = f_UX$xbw, ybw = f_UX$ybw, VX = X, VY = U))
        } else {
            pseudo_error <- t(sapply(all_X[,b], sample_conditional, n = lambda, xbw = f_UX$xbw, ybw = f_UX$ybw, VX = X, VY = U))
        }
        Xs.lambda <- Xs + rowSums(pseudo_error)
        method_2.theta.hat[[lambda+1]] <- method_2.theta.hat[[lambda+1]] + estimator(X=Xs.lambda)
    }
    
    method_2.theta.hat[[lambda+1]] <- (1/B)*method_2.theta.hat[[lambda+1]]
}

method_2_df <- data.frame(lambda = 0:M, theta = do.call(rbind, method_2.theta.hat))
colnames(method_2_df) <- c("lambda", "beta0", "beta1", "beta2")

plot(beta0 ~ lambda, data = data.frame(method_2_df))
plot(beta1 ~ lambda, data = data.frame(method_2_df))
plot(beta2 ~ lambda, data = data.frame(method_2_df))

beta0_mod_2_nls <- nls(beta0 ~ a + b/(c + lambda), data = method_2_df, start = list(a = 1, b = 1, c = 1))
beta0_mod_2_quad <- nls(beta0 ~ a + b*lambda + c*I(lambda^2), data = method_2_df, start = list(a = -1, b = 1, c = -1))
plot(beta0 ~ lambda, data = data.frame(method_2_df))
lines(predict(beta0_mod_2_nls) ~ method_2_df$lambda, col = 'red')
lines(predict(beta0_mod_2_quad) ~ method_2_df$lambda, col = 'blue')


beta1_mod_2_nls <- nls(beta1 ~ a + b/(c + lambda), data = method_2_df, start = list(a = 1, b = 1, c = 1))
beta1_mod_2_quad <- nls(beta1 ~ a + b*lambda + c*I(lambda^2), data = method_2_df, start = list(a = -1, b = 1, c = -1))
plot(beta1 ~ lambda, data = data.frame(method_2_df))
lines(predict(beta1_mod_2_nls) ~ method_2_df$lambda, col = 'red')
lines(predict(beta1_mod_2_quad) ~ method_2_df$lambda, col = 'blue')


beta2_mod_2_quad <- nls(beta2 ~ a + b*lambda + c*I(lambda^2), data = method_2_df, start = list(a = -1, b = 1, c = -1))
plot(beta2 ~ lambda, data = data.frame(method_2_df))
lines(predict(beta2_mod_2_quad) ~ method_2_df$lambda, col = 'blue')

b0_nl_2 <- predict(beta0_mod_2_nls, newdata=data.frame(lambda = -1))
b0_q_2 <- predict(beta0_mod_2_quad, newdata=data.frame(lambda = -1))
b1_nl_2 <- predict(beta1_mod_2_nls, newdata=data.frame(lambda = -1))
b1_q_2 <- predict(beta1_mod_2_quad, newdata=data.frame(lambda = -1))
b2_q_2 <- predict(beta2_mod_2_quad, newdata=data.frame(lambda = -1))

# Estimator 3: Drawing from U|X* directly
all_pseudo_errors <- t(sapply(Xs, sample_conditional, n = M*(M+1)*B/2, xbw = f_UX$xbw, ybw = f_UX$ybw, VX = val_Xs, VY = U))
method_3.theta.hat <- list(estimator(Xs))

for (lambda in 1:M) {
    method_3.theta.hat[[lambda+1]] <- 0
    
    col_start <- B*lambda*(lambda - 1)/2 + 1
    
    for(b in 1:B) {
        cat(paste0("\r Method 3: Lambda=", lambda, " and b=", b, "/", B))
        pseudo_error <- as.matrix(all_pseudo_errors[,col_start:(col_start+lambda-1)])
        Xs.lambda <- Xs + rowSums(pseudo_error)
        method_3.theta.hat[[lambda+1]] <- method_3.theta.hat[[lambda+1]] + estimator(X=Xs.lambda)
    }
    
    method_3.theta.hat[[lambda+1]] <- (1/B)*method_3.theta.hat[[lambda+1]]
}

method_3_df <- data.frame(lambda = 0:M, theta = matrix(unlist(method_3.theta.hat), nrow=(M+1), byrow=T))
colnames(method_3_df) <- c("lambda", "beta0", "beta1", "beta2")


beta0_mod_3_nls <- nls(beta0 ~ a + b/(c + lambda), data = method_3_df, start = list(a = 1, b = 1, c = 1))
beta0_mod_3_quad <- nls(beta0 ~ a + b*lambda + c*I(lambda^2), data = method_3_df, start = list(a = -1, b = 1, c = -1))
plot(beta0 ~ lambda, data = data.frame(method_3_df))
lines(predict(beta0_mod_3_nls) ~ method_3_df$lambda, col = 'red')
lines(predict(beta0_mod_3_quad) ~ method_3_df$lambda, col = 'blue')


beta1_mod_3_nls <- nls(beta1 ~ a + b/(c + lambda), data = method_3_df, start = list(a = 1, b = 1, c = 1))
beta1_mod_3_quad <- nls(beta1 ~ a + b*lambda + c*I(lambda^2), data = method_3_df, start = list(a = -1, b = 1, c = -1))
plot(beta1 ~ lambda, data = data.frame(method_3_df))
lines(predict(beta1_mod_3_nls) ~ method_3_df$lambda, col = 'red')
lines(predict(beta1_mod_3_quad) ~ method_3_df$lambda, col = 'blue')

beta2_mod_3_quad <- nls(beta2 ~ a + b*lambda + c*I(lambda^2), data = method_3_df, start = list(a = -1, b = 1, c = -1))
plot(beta2 ~ lambda, data = data.frame(method_3_df))
lines(predict(beta2_mod_3_quad) ~ method_3_df$lambda, col = 'blue')

b0_nl_3 <- predict(beta0_mod_3_nls, newdata=data.frame(lambda = -1))
b0_q_3 <- predict(beta0_mod_3_quad, newdata=data.frame(lambda = -1))
b1_nl_3 <- predict(beta1_mod_3_nls, newdata=data.frame(lambda = -1))
b1_q_3 <- predict(beta1_mod_3_quad, newdata=data.frame(lambda = -1))
b2_q_3 <- predict(beta2_mod_3_quad, newdata=data.frame(lambda = -1))

### Bootstrap Analysis
library(foreach)
library(doParallel)

cluster <- parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cluster)
results <- foreach(bs_rep = 1:1000) %dopar% {
    myDf <- tryCatch({
        bs_sample_M <- sample(1:nrow(mydf_main), nrow(mydf_main), TRUE)
        bs_sample_V <- sample(1:nrow(mydf_val), nrow(mydf_val), TRUE)

        bs_main <- mydf_main[bs_sample_M,]
        bs_val <- mydf_val[bs_sample_V,]
                
        Xs <- c(mydf_main$W, mydf_val$W)
        val_Xs <- mydf_val$W
        U <- mydf_val$W - mydf_val$X
        X <- mydf_val$X
        Z <- c(mydf_main$Z, mydf_val$Z)
        Y <- c(mydf_main$Y, mydf_val$Y)

        estimator <- function(X) { 
            MX <- model.matrix(~X+Z)
            M <- fastglm(y=Y, x=MX, family='binomial')
            M$coefficients
        }

        ##### Standard SIMEX
        sigma.hat <- var(U)
        naive_model <- glm(Y ~ Xs + Z , family=binomial)
        sm1 <- simex::simex(naive_model, SIMEXvariable="Xs", measurement.error = sigma.hat, lambda = seq(0,2,length.out=M)[2:M], B=B, jackknife.estimation=FALSE, asymptotic=FALSE)

        #### NP-SIMEX
        simex_df <- np_simex(Xs, estimator, U=U, B=B, M=M, parallel=FALSE)

        intercept_nls <- nls(theta.1 ~ a + b/(c + lambda), data = simex_df, start=list(a = 1, b = 1, c = 1))
        beta1_nls <- nls(theta.2 ~ a + b/(c + lambda), data = simex_df, start=list(a = 1, b = 1, c = 1))
        beta2_nls <- nls(theta.3 ~ a + b/(c + lambda), data = simex_df, start=list(a = 1, b = 1, c = 1))

        intercept.hat <- predict(intercept_nls, data.frame(lambda = -1))
        beta1.hat <- predict(beta1_nls, data.frame(lambda = -1))
        beta2.hat <- predict(beta2_nls, data.frame(lambda = -1))
        data.frame("Naive" = coef(naive_model), "P-SIMEX" = sm1$coefficients, "NP-SIMEX" = c(intercept.hat, beta1.hat, beta2.hat))
    },
    error = function(e) {
        return(NA)
    })
    myDf
}
parallel::stopCluster(cluster)

###############################################################################
# BS Conditional Analysis
stacked_bs_results <- do.call(rbind, lapply(results, function(r1){
    if(length(r1) > 1 || !is.na(r1)) {
        res <- unlist(r1)
        res <- res[c(1,2,4,5,6,7,8,10,11,12)]
        names(res) <- c("b0_nl_2", "b1_nl_2", "b0_q_2", "b1_q_2", "b2_q_2", 
                        "b0_nl_3", "b1_nl_3", "b0_q_3", "b1_q_3", "b2_q_3")
        res
    }
}))

resulting_df <- rbind(
    c(b0_nl_2, b1_nl_2, b0_q_2, b1_q_2, b2_q_2, b0_nl_3, b1_nl_3, b0_q_3, b1_q_3, b2_q_3),
    sqrt(diag(cov(stacked_bs_results)))
)

resulting_df <- data.frame(resulting_df)
resulting_df$estimand <- c("Estimate", "SE")

coef_table <- resulting_df %>% 
    pivot_longer(cols = -estimand, names_sep = "_", names_to = c("Parameter", "Extrapolant", "Method")) %>%
    pivot_wider(id_cols = c(Parameter, Extrapolant, Method), names_from = estimand, values_from = value) %>%
    pivot_wider(id_cols = c(Extrapolant, Method), names_from = Parameter, names_sep = "_", values_from = c(Estimate, SE)) %>%
    select(Method, Extrapolant, Estimate_b0, SE_b0, Estimate_b1, SE_b1, Estimate_b2, SE_b2)