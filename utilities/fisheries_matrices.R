fisheries_matrices <- function(states = 0:20,
                               actions = states,
                               observed_states = states,
                               reward_fn = function(x,a) pmin(x,a),
                               f = logistic,
                               sigma_g = 0.1,
                               sigma_m = 0.1,
                               r = 1,
                               K = 1,
                               noise = c("rescaled-lognormal", "lognormal", "uniform", "normal"),
                               C = 0.2){
  
 
  
  noise <- match.arg(noise)
  n_s <- length(states)
  n_a <- length(actions)
  n_z <- length(observed_states)
  transition <- array(0, dim = c(n_s, n_s, n_a))
  reward <- array(0, dim = c(n_s, n_a))
  observation <- array(0, dim = c(n_s, n_z, n_a))
  
  for (k in 1:n_s) {
    for (i in 1:n_a) {
      nextpop <- f(states[k], actions[i],r,K,C)
      if(nextpop <= 0){
        transition[k, , i] <- c(1, rep(0, n_s-1)) # c(.99, 0.01, rep(0, n_s-2))  # permit 1% chance of not going to 0 when f(x) <= 0?
      } else {
        transition[k, , i] <- prob(states, nextpop, sigma_g, noise = noise)
      }
      reward[k, i] <- reward_fn(states[k], actions[i])
    }
  }
  
  for (k in 1:n_a) {
    if(sigma_m <= 0){
      observation[, , k] <- diag(n_s)
    } else {
      for (i in 1:n_s) {
        if(states[i] <= 0)
          observation[i, , k] <- c(1, rep(0, n_z - 1))
        else {
          observation[i, , k] <- prob(observed_states, states[i], sigma_m, noise = noise)
        }
      }
    }
  }
  list(transition = transition, observation = observation, reward = reward)
}

prob <- function(states, mu, sigma, noise = "lognormal"){
  n_s <- length(states)
  if(noise == "rescaled-lognormal"){
    var <- mu * sigma^2 / 3 ## Rescale to the variance of uniform
    meanlog <- log( mu^2 / sqrt(var + mu^2) )
    sdlog <- sqrt( log(1 + var / mu^2) )
    x <- dlnorm(states, meanlog, sdlog)
    N <- plnorm(states[n_s], meanlog, sdlog)
  } else if(noise == "lognormal"){
    meanlog <- log(mu) - sigma ^ 2 / 2
    x <- dlnorm(states, meanlog, sigma)
    N <- plnorm(states[n_s], meanlog, sigma)
  } else if(noise == "normal"){
    x <- dnorm(states, mu, sigma)
    ## Pile negative density onto boundary
    negs <- dnorm(-states[-1], mu, sigma)
    x[1] <- x[1] + sum(negs)
    N <- pnorm(states[n_s], mu, sigma)
  } else if(noise == "uniform"){
    x <- dunif(states, mu * (1 - sigma), mu * (1 + sigma))
    N <- punif(states[n_s], mu * (1 - sigma), mu * (1 + sigma))
  }
  
  ## Handle exceptions (e.g. from mu ~ 0)
  if(sum(x) == 0){  ## nextpop is computationally zero, would create NAs
    x <- c(1, rep(0, n_s - 1))
    
    ## Normalize, pile on boundary
  } else {
    x <- x * N / sum(x)         # normalize densities to  = cdf(boundary)
    x[n_s] <- 1 - N + x[n_s]    # pile remaining probability on boundary
  }
  
  x
}

logistic <- function(x, h, r, K,C){
  S = x - h
  max(r * S * (1 - S / K) + S,0)
}

## Ricker growth function

ricker <- function(x, h, r, K,C){
  S = x - h
  max(S * exp(r * (1 - S / K)),0)
}

## Alen growth function
allen <- function(x, h, r, K,C){
  S = x - h
  #S * exp(r * (1 - S / K) * ((S - C))) 
  max(S * exp(r * (1 - S / K) * ((S - C)/K)),0)
}

## Beverton-Holt growth function
bevertonholt <- function(x, h, r, K,C){
  S = x - h
  max((1 + r) * S / (1 + (r * S) / K),0)
}