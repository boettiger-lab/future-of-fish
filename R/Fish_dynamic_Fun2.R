# Needed functions for Fish_Dynamic_MM

## stupid aliases for some R functions to  ~ their matlab names
zeros <- function(x, y, z=1, fill=0){
  A <- array(fill, dim = c(x, y, z))
  if(z == 1)
    return(A[,,1])
  else
    return(A)
}
ones <- function(x,y,z=1) zeros(x,y,z, fill=1)

size <- function(x, d=NULL){
  if(is.array(x)){
    if(is.null(d))
      return(dim(x))
    else
      return(dim(x)[d])
  } else
    return(length(x))
}


discrete <- function(p,n){
  s = zeros(1,n)
  for (i in 1:n){
    s[i] = sum(cumsum(p)<runif(1))+1
  }
  s
}

normalize <- function(A){
  
  z = sum(A)
  s = z + (z==0)
  A / s
  
}


round_Milad <- function(P){
  
  if(!is.array(P)){
    P <- array(P, dim=c(length(P), 1,1))
  }
  
  num_dig = 4;
  
  for(j in 1:dim(P)[1]){
    for(l in 1:dim(P)[2]){
      P[j,l] = round(P[j,l] * (10 ^ num_dig)) / 10 ^ num_dig
    }
  }
  
  
  for(k in 1:dim(P)[1]){
    if( sum(P[k,]) != 1){
      ind = which.max(P[k,]);
      P[k,ind]=P[k,ind] + (1.0000 - sum(P[k,]))
    }
  }
  
  P
  
}


## Logistic growth function
logistic <- function(x, h, r, K){
  S = x - h
  r * S * (1 - S / K) + S
}

## Ricker growth function

ricker <- function(x, h, r, K){
  S = x - h
  S * exp(r * (1 - S / K))
}

## Alen growth function
allen <- function(x, h, r, C, K){
  S = x - h
  #S * exp(r * (1 - S / K) * ((S - C))) 
  S * exp(r * (1 - S / K) * ((S - C)/K)) 
}

## Beverton-Holt growth function
bevertonholt <- function(x, h, r, K){
  S = x - h
  (1 + r) * S / (1 + (r * S) / K)
}

## Linear growth function
linear <- function(x, h, a, b){
  S = x - h
  a * S + b
}


## Transition Matrix function
create_T <- function(states,actions,sigma_g,f,r,C,K,noise_dist){
  a = 1.064; b = round(1160.686/100) 
  T = zeros(length(states), length(states), length(actions))
  
  if(sigma_g == 0){
    for (i in 1:length(states)) {
      
      for (j in 1:length(actions)){
        
        x = states[i]
        a = actions[j]
        h = min(x,a)
        s = x-h
        if(f == "logistic"){
          G = logistic(x,h,r,K)
        } else if(f == "allen"){
          G = allen(x,h,r,C,K)
        } else if(f == "ricker"){
          G = ricker(x,h,r,K)
        } else if(f == "bevertonholt"){
          G = bevertonholt(x,h,r,K)
        } else if(f == "linear"){
          G = linear(x,h,a,b)
        }
        #G = floor(G)
        
       
        if((G<=length(states)) && (G!=0)){
          cc = which.min(abs(states-G))
          T[i,cc,j] = 1
        } else if(G == 0){
          T[i,1,j] = 1
        } else if(G>length(states)){
          T[i,length(states),j] = 1
        }
      }
    }
  } else if(noise_dist == "uniform"){
    
    range = zeros(length(states), length(actions), 2)
    
    for (i in 1:length(states)) {
      for (j in 1:length(actions)) {
        x = states[i]
        a = actions[j]
        h = min(x, a)
        s = x - h
        if(f == "logistic"){
          G = logistic(x,h,r,K)
        } else if(f == "allen"){
          G = allen(x,h,r,C,K)
        } else if(f == "ricker"){
          G = ricker(x,h,r,K)
        } else if(f == "bevertonholt"){
          G = bevertonholt(x,h,r,K)
        } else if(f == "linear"){
          G = linear(x,h,a,b)
        }
        range[i, j, 1] = ((1-sigma_g) * G)
        range[i, j, 2] = ((1+sigma_g) * G)
        
      }
    }
    
    for (i in 1:size(range, 1)) {
      for (j in 1:size(range, 2)) {
        m = range[i, j, ]
        #m = matrix(m, nrow = size(m, 2), ncol = size(m, 3))
        #m = m + 1
        if(m[1]>states[length(states)]){
          T[i,length(states),j] = 1
        } else{
          nnn = max(states[1], m[1])
          mmm = min(states[length(states)], m[2])
          nnn = to_grid(nnn,states)
          mmm = to_grid(mmm,states)
          mm = 1 / (mmm - nnn + 1)
          T[i, (nnn:mmm), j] = mm
          #T[i, mmm, j] = 1 - ((mmm - nnn + 1) * mm)
        }
        
      }
    }
    
  } else if (noise_dist == "lognormal"){
    z = zeros(length(states), length(states), length(actions))
    # var_g_x = (1+sigma_g-(1-sigma_g))^2 / 12
    # var_g_logn =  log(1+var_g_x)
    # sigma_g_logn = sqrt(var_g_logn)
    # mu_g_logn = log(1) - var_g_logn/2
    
    sigma_g_logn = sigma_g
    for (i in 1:length(states)) {
      for (j in 1:length(actions)) {
        x = states[i]
        a = actions[j]
        h = min(x, a)
        s = x - h
        if(f == "logistic"){
          G = logistic(x,h,r,K)
          #G = floor(G)
        } else if(f == "allen"){
          G = allen(x,h,r,C,K)
          #G = floor(G)
        } else if(f == "ricker"){
          G = ricker(x,h,r,K)
          #G = floor(G)
        } else if(f == "bevertonholt"){
          G = bevertonholt(x,h,r,K)
          #G = floor(G)
        } else if(f == "linear"){
          G = linear(x,h,a,b)
          #G = floor(G)
        }
        
        
        
        if(G <= 0){
          T[i,1,j] = 1
        } else{
          for (ii in 1:length(states)){
            z[i,ii,j] = states[ii] / G
            #T[i,ii,j] = dlnorm(z[i,ii,j],mu_g_logn,sigma_g_logn)
            T[i,ii,j] = dlnorm(states[ii],(log(G) - sigma_g^2/2),sigma_g)
            
          }
        }
        }
      }
   
  } else {
    stop("Noise distribution not recognized")
  }
  # normalizing transition matrix
  for (i in 1:size(T, 1)) {
    for (j in 1:size(T, 3)) {
      T[i, , j] = normalize(T[i, , j])
      
    }
  }
  #for(i in 1:size(T,3)){
  #  T[,,i] = round_Milad(T[,,i])
  #}
  
  output <- T
}


create_T_eve <- function(Num_s,Num_a,sigma_g,f,r,C,K,noise_dist){
  
  T = zeros(Num_s, Num_s, Num_a)
  
  if(sigma_g == 0){
    for (i in 1:Num_s) {
      for (j in 1:Num_a){
        x = round(i/2)
        a = j-1
        h = min(x,a)
        s = x-h
        if(f == "logistic"){
          G = logistic(x,h,r,K)
        } else if(f == "allen"){
          G = allen(x,h,r,C,K)
        } else if(f == "ricker"){
          G = ricker(x,h,r,K)
        } else if(f == "bevertonholt"){
          G = bevertonholt(x,h,r,K)
        }
        G = round(G)
        if((G<=Num_s) && (G!=0)){
          T[i,G,j] = 1
        } else if(G == 0){
          T[i,1,j] = 1
        } else if(G>Num_s){
          T[i,Num_s,j] = 1
        }
      }
    }
  } else if(noise_dist == "uniform"){
    
    range = zeros(Num_s, Num_a, 2)
    
    for (i in 1:Num_s) {
      for (j in 1:Num_a) {
        x = round(i/2)
        a = j - 1
        h = min(x, a)
        s = x - h
        if(f == "logistic"){
          G = logistic(x,h,r,K)
        } else if(f == "allen"){
          G = allen(x,h,r,C,K)
        } else if(f == "ricker"){
          G = ricker(x,h,r,K)
        } else if(f == "bevertonholt"){
          G = bevertonholt(x,h,r,K)
        }
        range[i, j, 1] = floor((1-sigma_g) * G)
        range[i, j, 2] = floor((1+sigma_g) * G)
        
      }
    }
    
    for (i in 1:size(range, 1)) {
      for (j in 1:size(range, 2)) {
        m = range[i, j, ]
        #m = matrix(m, nrow = size(m, 2), ncol = size(m, 3))
        m = m + 1
        if(m[1]>Num_s){
          T[i,Num_s,j] = 1
        } else{
          nnn = max(1, m[1])
          mmm = min(Num_s, m[2])
          mm = 1 / (max(m) - nnn + 1)
          T[i, nnn:(mmm - 1), j] = mm
          T[i, mmm, j] = 1 - ((mmm - nnn + 1) * mm)
        }
        
      }
    }
    
  } else if (noise_dist == "lognormal"){
    z = zeros(Num_s, Num_s, Num_a)
    var_g_x = (1+sigma_g-(1-sigma_g))^2 / 12
    var_g_logn =  log(1+var_g_x)
    sigma_g_logn = sqrt(var_g_logn)
    mu_g_logn = log(1) - var_g_logn/2
    for (i in 1:Num_s) {
      for (j in 1:Num_a) {
        x = round(i/2)
        a = j - 1
        h = min(x, a)
        s = x - h
        if(f == "logistic"){
          G = logistic(x,h,r,K)
          G = floor(G)
        } else if(f == "allen"){
          G = allen(x,h,r,C,K)
          G = floor(G)
        } else if(f == "ricker"){
          G = ricker(x,h,r,K)
          G = floor(G)
        } else if(f == "bevertonholt"){
          G = bevertonholt(x,h,r,K)
          G = floor(G)
        }
        for (ii in 1:Num_s){
          if (G == 0){
            T[i,1,j] = 1
          } else{
            z[i,ii,j] = ii / G
            T[i,ii,j] = dlnorm(z[i,ii,j],mu_g_logn,sigma_g_logn)
          }
        }
      }
    }
  } else {
    stop("Noise distribution not recognized")
  }
  # normalizing transition matrix
  for (i in 1:size(T, 1)) {
    for (j in 1:size(T, 3)) {
      T[i, , j] = normalize(T[i, , j])
      
    }
  }
  for(i in 1:size(T,3)){
    T[,,i] = round_Milad(T[,,i])
  }
  
  output <- T
}


## Emission matrix function
create_O <- function(states,actions,obs,sigma_m,noise_dist){
  
  O = zeros(length(states), length(obs), length(actions))
  
  if(sigma_m == 0){
    for (j in 1:Num_a){
      O[ , ,j] = diag(Num_s)
    }
  } else if(noise_dist == "uniform"){
    
    for (i in 1:length(states)) {
      for (j in 1:length(actions)) {
        m = (c((1-sigma_m),(1+sigma_m)) * states[i])
        nnn = max(states[1], m[1])
        mmm = min(obs[length(obs)], m[2])
        nnn = to_grid(nnn,states)
        mmm = to_grid(mmm,states)
        mm = 1 / (mmm - nnn + 1)
        
        #if (m[2] == states[1]) {
        #  O[i, 1, j] = mm
        #} else {
        O[i, (nnn:mmm), j] = mm
        #O[i, mmm, j] = 1 - ((mmm - nnn + 1) * mm)
        #}
      }
    }
    
  } else if (noise_dist == "lognormal"){
    var_m_x = sigma_m^2 # (1+sigma_m-(1-sigma_m))^2 / 12
    var_m_logn =  log(1+var_m_x)
    sigma_m_logn = sigma_m# sqrt(var_m_logn)
    mu_m_logn = log(1) - var_m_logn/2
    zz = zeros(length(states), length(obs), length(actions))
    for (i in 1:length(states)){
      for (j in 1:length(actions)){
        for (ii in 1:length(obs)){
          if(i == 1){
            O[i,1,j] = 1
          } else{
            zz[i,ii,j] = obs[ii] /states[i]
            O[i,ii,j] = dlnorm(zz[i,ii,j],mu_m_logn,sigma_m_logn)
          }
        }
      }
    }
    
  } else {
    stop("Noise distribution not recognized")
  }
  for (i in 1:size(O, 1)) {
    for (j in 1:size(O, 3)) {
      O[i, , j] = normalize(O[i, , j])
    }
  }
  #for(i in 1:size(O,3)){
  #  O[,,i] = round_Milad(O[,,i])
  #}
  
  output <- O
}

## Reward Function
create_R <- function(states,actions){
  R = zeros(length(states), length(actions),1)
  for (i in 1:length(states)) {
    for (j in 1:length(actions)) {
      x = states[i]
      a = actions[j]
      #if (a <= x)
      R[i, j] = min(x, a)
      #else
      # R[i, j] = -1
    }
  }
  output <- R
}


## Reward Function
create_R2 <- function(states,actions,alpha){
  R = zeros(length(states), length(actions),1)
  for (i in 1:length(states)) {
    for (j in 1:length(states)) {
      x = states[i]
      a = actions[j]
      #if (a <= x)
        R[i, j] = min(x, a) - (alpha * a)
      #else
      #  R[i, j] = -1
    }
  }
  output <- R
}


create_R3 <- function(states,actions,alpha){
  R = zeros(length(states), length(actions),1)
  for (i in 1:length(states)) {
    for (j in 1:length(actions)) {
      x = states[i]
      a = actions[j]
      #if (a <= x)
        R[i, j] = min(x, a) - (alpha * a / x)
      #else
      #  R[i, j] = -1
    }
  }
  output <- R
}

## running POMDP

runPOMDP <- function(Num_s,Num_a,Num_z,initial,T,O,R,GAMMA){
  
  out = vector("numeric", length = Num_z)
  value = vector("numeric", length = Num_z)
  policy = vector("numeric", length = Num_z)
  av = vector("list",length = Num_z)
  aa = vector("list",length = Num_z)
  actions = paste0("a", 1:Num_a)
  XX = paste0("a", 1:800)
  
  for (i in 1:Num_z) {
    print(i)
    belief = initial * O[, i, 1]
    belief = normalize(belief)
    belief = round(belief,4) / sum(round(belief,4))
    if(is.nan(sum(belief))){
      value[i] = 0
      policy[i] = 0
    } else{
      milad_writepomdp(T, O, R, GAMMA, belief, Num_s, Num_a, Num_z, actions, XX)
      command = "/Users/Milad/Documents/appl-0.96/src/pomdpsol input.pomdp --output output.policy --precision 1 --timeout 200";
      run_MOMDP(command)
      out = findpolicy_Milad2(belief)
      value[i] = out[[1]]
      policy[i] = out[[2]]
      av[[i]] = out[[3]]
      aa[[i]] = out[[4]]
    }
    
    
  }
  list(value,policy,av,aa)
  
}



create_T_unc <- function(states,actions,sigma_g,f,r,C,K,noise_dist,w){
  T = zeros(length(states), length(states), length(actions))
  z = zeros(length(states), length(states), length(actions))
  # var_g_x = sigma_g^2# (1+sigma_g-(1-sigma_g))^2 / 12
  # var_g_logn =  log(1+var_g_x)
  # sigma_g_logn = sigma_g#sqrt(var_g_logn)
  # mu_g_logn = log(1) - var_g_logn/2
  # 
  
  sigma_g_logn = sigma_g
 
  
  if(sigma_g == 0){
    for (i in 1:length(states)) {
      for (j in 1:length(actions)) {
        x = states[i]
        a = actions[j]
        h = min(x, a)
        s = x - h
        G = array(0, dim = length(r)*length(K))
        W_g = array(0, dim = length(r)*length(K))
        ll = 1
        for(rr in 1:length(r)){
          for(kk in 1:length(K)){
            if(f == "logistic"){
              G[ll] = logistic(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "allen"){
              G[ll] = allen(x,h,r[rr],C,K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "ricker"){
              G[ll] = ricker(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "bevertonholt"){
              G[ll] = bevertonholt(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "linear"){
              G = linear(x,h,a,b)
              #G = floor(G)
            }
            ll = ll+1
          }
        }
        
        GG = sum(W_g * G)
        if((GG<=states[length(states)]) && (GG!=0)){
          cc = which.min(abs(states-GG))
          T[i,cc,j] = 1
        } else if(GG == 0){
          T[i,1,j] = 1
        } else if(GG>states[length(states)]){
          T[i,length(states),j] = 1
        }
        
        
      }
      
    } 
  } else {
    for (i in 1:length(states)) {
      for (j in 1:length(actions)) {
        x = states[i]
        a = actions[j]
        h = min(x, a)
        s = x - h
        G = array(0, dim = length(r)*length(K))
        W_g = array(0, dim = length(r)*length(K))
        ll = 1
        for(rr in 1:length(r)){
          for(kk in 1:length(K)){
            if(f == "logistic"){
              G[ll] = logistic(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "allen"){
              G[ll] = allen(x,h,r[rr],C,K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "ricker"){
              G[ll] = ricker(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "bevertonholt"){
              G[ll] = bevertonholt(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "linear"){
              G = linear(x,h,a,b)
              #G = floor(G)
            }
            ll = ll+1
          }
        }
        
        GG = sum(W_g * G)
        
        mu_g_logn = log(GG) - sigma_g^2/2
        if(GG == 0){
          T[i,1,j] = 1
        } else{
          for (ii in 1:length(states)){
            #z[i,ii,j] = states[ii] / GG
            #T[i,ii,j] = dlnorm(z[i,ii,j],mu_g_logn,sigma_g_logn)
            T[i,ii,j] = dlnorm(states[ii],mu_g_logn,sigma_g_logn)
          }
          if(sum(T[i,,j]) == 0){
            T[i,1,j] = 1
          }
        }
      }
      
    } 
  }
  
  
  
  
  
  # normalizing transition matrix
  for (i in 1:size(T, 1)) {
    for (j in 1:size(T, 3)) {
      T[i, , j] = normalize(T[i, , j])
      
    }
  }
  #for(i in 1:size(T,3)){
  #  T[,,i] = round_Milad(T[,,i])
  #}
  
  output <- T
}




create_T_unc2 <- function(states,actions,sigma_g,f,r,C,K,noise_dist,w){
  T = zeros(length(states), length(states), length(actions))
  z = zeros(length(states), length(states), length(actions))
  # var_g_x = sigma_g^2# (1+sigma_g-(1-sigma_g))^2 / 12
  # var_g_logn =  log(1+var_g_x)
  # sigma_g_logn = sigma_g#sqrt(var_g_logn)
  # mu_g_logn = log(1) - var_g_logn/2
  # 
  
  sigma_g_logn = sigma_g
  
  
  if(sigma_g == 0){
    for (i in 1:length(states)) {
      for (j in 1:length(actions)) {
        x = states[i]
        a = actions[j]
        h = min(x, a)
        s = x - h
        G = array(0, dim = length(r)*length(K))
        W_g = array(0, dim = length(r)*length(K))
        ll = 1
        for(rr in 1:length(r)){
          for(kk in 1:length(K)){
            if(f == "logistic"){
              G[ll] = logistic(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "allen"){
              G[ll] = allen(x,h,r[rr],C,K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "ricker"){
              G[ll] = ricker(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "bevertonholt"){
              G[ll] = bevertonholt(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "linear"){
              G = linear(x,h,a,b)
              #G = floor(G)
            }
            ll = ll+1
          }
        }
        
        GG = sum(W_g * G)
        if((GG<=states[length(states)]) && (GG!=0)){
          cc = which.min(abs(states-GG))
          T[i,cc,j] = 1
        } else if(GG == 0){
          T[i,1,j] = 1
        } else if(GG>states[length(states)]){
          T[i,length(states),j] = 1
        }
        
        
      }
      
    } 
  } else {
    for (i in 1:length(states)) {
      for (j in 1:length(actions)) {
        x = states[i]
        a = actions[j]
        h = min(x, a)
        s = x - h
        G = array(0, dim = length(r)*length(K))
        W_g = array(0, dim = length(r)*length(K))
        ll = 1
        for(rr in 1:length(r)){
          for(kk in 1:length(K)){
            if(f == "logistic"){
              G[ll] = logistic(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "allen"){
              G[ll] = allen(x,h,r[rr],C,K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "ricker"){
              G[ll] = ricker(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "bevertonholt"){
              G[ll] = bevertonholt(x,h,r[rr],K[kk])
              W_g[ll] = w[rr,kk]
              #G = floor(G)
            } else if(f == "linear"){
              G = linear(x,h,a,b)
              #G = floor(G)
            }
            ll = ll+1
          }
        }
        
        GG = sum(W_g * G)
        
        mu_g_logn = log(GG) - sigma_g^2/2
        
        if(GG <= 0){
          T[i, , j] <- c(.9999, 0.0001, rep(0, length(states)-2))  # permit s chance of not going to 0 when f(x) <= 0
        } else {
          T[i, , j] <- prob(states, GG, sigma_g, noise = noise_dist)
        }
      
        # if(GG == 0){
        #   T[i,1,j] = 1
        # } else{
        #   for (ii in 1:length(states)){
        #     #z[i,ii,j] = states[ii] / GG
        #     #T[i,ii,j] = dlnorm(z[i,ii,j],mu_g_logn,sigma_g_logn)
        #     T[i,ii,j] = dlnorm(states[ii],mu_g_logn,sigma_g_logn)
        #   }
        #   if(sum(T[i,,j]) == 0){
        #     T[i,1,j] = 1
        #   }
        # }
      }
      
    } 
  }
  
  
  # 
  # 
  # 
  # # normalizing transition matrix
  # for (i in 1:size(T, 1)) {
  #   for (j in 1:size(T, 3)) {
  #     T[i, , j] = normalize(T[i, , j])
  #     
  #   }
  # }
  #for(i in 1:size(T,3)){
  #  T[,,i] = round_Milad(T[,,i])
  #}
  
  output <- T
}


prob <- function(states, mu, sigma, noise = "lognormal"){
  n_s <- length(states)
  if(noise == "lognormal"){
    ## Rescaling for apples-to-apples comparison bewtween a lognormal and uniform with parameter sigma
    #  var <- mu * sigma^2 / 3 ## Rescale to the variance of uniform
    #  meanlog <- log( mu^2 / sqrt(var + mu^2) )
    #  sdlog <- sqrt( log(1 + var / mu^2) )
    meanlog <- log(mu) - sigma ^ 2 / 2
    x <- dlnorm(states, meanlog, sigma)
    N <- plnorm(states[n_s], meanlog, sigma)
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


create_O2 <- function(states,actions,obs,sigma_m,noise_dist){
  
  O = zeros(length(states), length(obs), length(actions))
  
  if(sigma_m == 0){
    for (j in 1:Num_a){
      O[ , ,j] = diag(Num_s)
    }
  } else if (noise_dist == "lognormal"){
    
    # var_m_x = sigma_m^2 # (1+sigma_m-(1-sigma_m))^2 / 12
    # var_m_logn =  log(1+var_m_x)
    # sigma_m_logn = sigma_m# sqrt(var_m_logn)
    # mu_m_logn = log(1) - var_m_logn/2
    zz = zeros(length(states), length(obs), length(actions))
    for (i in 1:length(states)){
      for (j in 1:length(actions)){
        for (ii in 1:length(obs)){
          if(i == 1){
            O[i,1,j] = 1
          } else{
            zz[i,ii,j] = obs[ii] /states[i]
            O[i,ii,j] = dlnorm(zz[i,ii,j],(log(1) - sigma_m^2 /2),sigma_m)
          }
        }
      }
    }
    
  } else {
    stop("Noise distribution not recognized")
  }
  for (i in 1:size(O, 1)) {
    for (j in 1:size(O, 3)) {
      O[i, , j] = normalize(O[i, , j])
    }
  }
  #for(i in 1:size(O,3)){
  #  O[,,i] = round_Milad(O[,,i])
  #}
  
  output <- O
}




