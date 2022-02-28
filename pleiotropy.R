library("dplyr")
library("broom")
library("GWASTools")
source("funs.R")
set.seed(123)

n_sim <- 200
n_obs <- 10000
q <- 0.25

# Taken from https://stats.stackexchange.com/a/477214
rcorrbinom <- function(n, size = 1, prob1, prob2, corr = 0) {
  
  #Check inputs
  if (!is.numeric(n))             { stop('Error: n must be numeric') }
  if (length(n) != 1)             { stop('Error: n must be a single number') }
  if (as.integer(n) != n)         { stop('Error: n must be a positive integer') }
  if (n < 1)                      { stop('Error: n must be a positive integer') }
  if (!is.numeric(size))          { stop('Error: n must be numeric') }
  if (length(size) != 1)          { stop('Error: n must be a single number') }
  if (as.integer(size) != size)   { stop('Error: n must be a positive integer') }
  if (size < 1)                   { stop('Error: n must be a positive integer') }
  if (!is.numeric(prob1))         { stop('Error: prob1 must be numeric') }
  if (length(prob1) != 1)         { stop('Error: prob1 must be a single number') }
  if (prob1 < 0)                  { stop('Error: prob1 must be between 0 and 1') }
  if (prob1 > 1)                  { stop('Error: prob1 must be between 0 and 1') }
  if (!is.numeric(prob2))         { stop('Error: prob2 must be numeric') }
  if (length(prob2) != 1)         { stop('Error: prob2 must be a single number') }
  if (prob2 < 0)                  { stop('Error: prob2 must be between 0 and 1') }
  if (prob2 > 1)                  { stop('Error: prob2 must be between 0 and 1') }
  if (!is.numeric(corr))          { stop('Error: corr must be numeric') }
  if (length(corr) != 1)          { stop('Error: corr must be a single number') }
  if (corr < -1)                  { stop('Error: corr must be between -1 and 1') }
  if (corr > 1)                   { stop('Error: corr must be between -1 and 1') }
  
  #Compute probabilities
  P00   <- (1-prob1)*(1-prob2) + corr*sqrt(prob1*prob2*(1-prob1)*(1-prob2))
  P01   <- 1 - prob1 - P00
  P10   <- 1 - prob2 - P00
  P11   <- P00 + prob1 + prob2 - 1
  PROBS <- c(P00, P01, P10, P11)
  if (min(PROBS) < 0)       { stop('Error: corr is not in the allowable range') }
  
  #Generate the output
  RAND <- array(sample.int(4, size = n*size, replace = TRUE, prob = PROBS),
                dim = c(n, size))
  VALS <- array(0, dim = c(2, n, size))
  OUT  <- array(0, dim = c(2, n))
  for (i in 1:n)    { 
  for (j in 1:size) { 
    VALS[1,i,j] <- (RAND[i,j] %in% c(3, 4))
    VALS[2,i,j] <- (RAND[i,j] %in% c(2, 4)) } 
    OUT[1, i]   <- sum(VALS[1,i,])
    OUT[2, i]   <- sum(VALS[2,i,]) }
  
  #Give output
  OUT 
}

# horizontal pleiotropy Z-X
results <- data.frame()
for (alpha in c(0, 1)){
    for (i in 1:n_sim){
        c <- rnorm(n_obs)
        z <- get_simulated_genotypes(q, n_obs)
        x <- z + c + rnorm(n_obs)
        y <- x + c + z*alpha + rnorm(n_obs)
        e <- lm(x~z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::pull(estimate)
        o <- lm(y~z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::pull(estimate)
        w <- o / e
        obs <- lm(y~x) %>% tidy %>% dplyr::filter(term == "x") %>% dplyr::pull(estimate)
        v <- model(data.frame(z, x), "z", "x")
        results <- rbind(results, data.frame(
            w, obs, alpha, phi_p=v$phi_p
        ))
    }
}

# horizontal pleiotropy Z-X1-X2-Y
results <- data.frame()
for (alpha in c(0, 1)){
    for (i in 1:n_sim){
        c <- rnorm(n_obs)
        z <- get_simulated_genotypes(q, n_obs)
        x1 <- z + c + rnorm(n_obs)
        x2 <- x + c + z*alpha + rnorm(n_obs)
        y <- x2 + c + z*alpha + rnorm(n_obs)
        e <- lm(x2~z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::pull(estimate)
        o <- lm(y~z) %>% tidy %>% dplyr::filter(term == "z") %>% dplyr::pull(estimate)
        w <- o / e
        obs <- lm(y~x2) %>% tidy %>% dplyr::filter(term == "x2") %>% dplyr::pull(estimate)
        v <- model(data.frame(z, x2), "z", "x2")
        results <- rbind(results, data.frame(
            w, obs, alpha, phi_p=v$phi_p
        ))
    }
}

# phantom vQTL
results <- data.frame()
for (i in 1:n_sim){
    # correlation between Z1 and Z2
    rho <- 0.4
    # varaince explained by Z1-X
    z1_r2 <- 0.05
    # simulate SNPs
    q1 <- 0.1; p1 <- 1-q1; v1 <- 2*p1*q1
    q2 <- 0.25; p2 <- 1-q2; v2 <- 2*p2*q2
    z1_b <- sqrt(z1_r2/v1)
    z <- rcorrbinom(n_obs, size = 2, q1, q2, corr = rho)
    # simulate exposure with Z1-X main effect
    x <- z[1,]*z1_b + rnorm(n_obs, sd=sqrt(1-z1_r2))
    # test for Z1-X and Z2-X variance effects
    phi_p1 <- model(data.frame(z=z[1,], x), "z", "x")$phi_p
    phi_p2 <- model(data.frame(z=z[2,], x), "z", "x")$phi_p
    # store result
    results <- rbind(results, data.frame(rho, phi_p1, phi_p2))
}
# QQplot
pdf("qq1.pdf")
GWASTools::qqPlot(results$phi_p1)
dev.off()
pdf("qq2.pdf")
GWASTools::qqPlot(results$phi_p2)
dev.off()

# Z-X horizontal pleiotropy
results <- data.frame()
for (alpha in c(0, 1)){
    for (i in 1:n_sim){
        c <- rnorm(n_obs)
        z <- get_simulated_genotypes(q, n_obs)
        x1 <- z*1 + c + rnorm(n_obs)
        x2 <- z*2 + c + rnorm(n_obs)
        x3 <- z*3 + c + rnorm(n_obs)
        x <- x1*5 + x2*75 + x3*3 + rnorm(n_obs)
        f_exp <- summary(lm(x ~ z))$fstatistic[1] %>% as.numeric
        v <- model(data.frame(x, z), "z", "x")
        results <- rbind(results, data.frame(
            phi_p=v$phi_p, f_exp
        ))
    }
}