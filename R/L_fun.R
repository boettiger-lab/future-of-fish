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

discrete_sample <- function(p,n)
{
  stopifnot(length(p) > 1)
  c = cumsum(p)
  
  icmf = function(pp){
    min(which(pp < c))
  }
  ps = runif(n)
  sapply(ps, icmf)
}

runPOMDP_L <- function(Num_s,Num_a,Num_z,initial,T,O,R,GAMMA){
  
  out = list()
  value = list()
  policy = list()
  alpha_vector = list()
  alpha_action = list()
  actions = paste0("a", 1:Num_a)
  XX = paste0("a", 1:800)
  
  if(is.matrix(initial)){
    siz = size(initial,1)
  } else{
    siz = 1
  }

  for (i in 1:siz) {
    if(siz==1){
      belief = initial
    } else{
      belief = initial[i,]
    }
    belief = normalize(belief)
    belief = round(belief,4) / sum(round(belief,4))
    if(siz==1){
      milad_writepomdp(T,O,R,GAMMA,belief, Num_s, Num_a, Num_z, actions, XX)
    } else{
      milad_writepomdp(T[[i]], O[[i]], R, GAMMA, belief, Num_s, Num_a, Num_z, actions, XX)
    }
    out = run_Milad(belief)
    value[i] = out[[1]]
    policy[i] = out[[2]]
    alpha_vector[[i]] = out[[3]]
    alpha_action[[i]] = out[[4]]
  }
  list(value,policy,alpha_vector,alpha_action)
  
}
Expectation_star <- function(initial,alpha_vector,alpha_action,Num_a){
  V = vector('list',length(alpha_vector))
  for(i in 1:length(V)){
    for(j in 1:Num_a){
      V[[i]][j] = -Inf
    }
  }
  
  for(i in 1:length(alpha_vector)){
      for(k in 1:Num_a){
        if(alpha_action[[i]] == (k)){
          vv = sum(initial*alpha_vector[[i]])
          V[[i]][k] = max(vv,V[[i]][k])
        }
      }
  }
  
  indd = zeros(1,Num_a)
  for(i in 1:length(V)){
    for(j in 1:length(V[[i]])){
      if(V[[i]][j] == -Inf){
        V[[i]][j] = 0
        indd[j] = indd[j] + 1
      }
    }
  }
  
  V_action = zeros(1,Num_a)
  for(i in 1:Num_a){
    summ = 0
    for(j in 1:length(V)){
      summ = summ + (1 / (length(V) - indd[i])) * V[[j]][i]
    }
    V_action[i]  = summ
  }
  for(i in 1:length(V_action)){
    if(is.nan(V_action[i])){
      V_action[i] = - Inf 
    }
  }
  cc = max(V_action)
  ind = which.max(V_action)
  ExpectedV  = cc
  ExpecctedPol = ind
  list(ExpectedV, ExpecctedPol)
}

Expectation_pl <- function(initial,alpha_vector,alpha_action,Num_a,w){
  V = vector('list',size(alpha_vector))
  for(i in 1:size(V)){
    for(j in 1:Num_a){
      V[[i]][j] = -Inf
    }
  }
  
  for(i in 1:size(alpha_vector)){
    for(j in 1:size(alpha_vector[[i]])){
      for(k in 1:Num_a){
        if(alpha_action[[i]][j] == k){
          vv = sum(initial[i,]*alpha_vector[[i]][[j]])
          V[[i]][k] = max(vv,V[[i]][k])
        }
      }
    }
  }
  V_action = zeros(1,Num_a)
  for(i in 1:Num_a){
    summ = 0
    for(j in 1:size(V)){
      if(V[[j]][i] != -Inf){
        summ = summ + V[[j]][i] * w[j]
      }
    }
    V_action[i]  = summ
  }
  for(i in 1:length(V_action)){
    if(is.nan(V_action[i])){
      V_action[i] = -1
    }
  }
  cc = max(V_action)
  ind = which.max(V_action)
  ExpectedV  = cc
  ExpecctedPol = ind
  list(ExpectedV, ExpecctedPol)
}

Expectation_opt <- function(initial,alpha_vector,alpha_action,Num_a){
  
  if(is.list(alpha_action)){
    V = vector('list',size(alpha_vector))
    for(i in 1:size(V)){
      for(j in 1:Num_a){
        V[[i]][j] = -Inf
      }
    }
    for(i in 1:size(alpha_vector)){
      for(j in 1:size(alpha_vector[[i]])){
        for(k in 1:Num_a){
          if(alpha_action[[i]][j] == k){
            if(is.vector(initial)){
              vv = sum(initial*alpha_vector[[i]][[j]])
            } else{
              vv = sum(initial[i,]*alpha_vector[[i]][[j]])
            }
            V[[i]][k] = max(vv,V[[i]][k])
          }
        }
      }
    }
    v = zeros(1,size(V))
    p = zeros(1,size(V))
    for(i in 1:size(V)){
      v[i] = max(V[[i]])
      p[i] = which.max(V[[i]])
    }
  cc = max(v)
  ind = which.max(v)
  ExpectedV  = cc
  ExpecctedPol = p[ind]
  list(ExpectedV, ExpecctedPol)
  
  } else{
    out = Interp_MM(initial,alpha_vector,alpha_action)
    ExpectedV  = out[[1]]
    ExpecctedPol = out[[2]]
    list(ExpectedV, ExpecctedPol)
  }
}

Expectation_pess <- function(initial,alpha_vector,alpha_action,Num_a){
  V = vector('list',size(alpha_vector))
  for(i in 1:size(V)){
    for(j in 1:Num_a){
      V[[i]][j] = -Inf
    }
  }
  
  for(i in 1:size(alpha_vector)){
    for(j in 1:size(alpha_vector[[i]])){
      for(k in 1:Num_a){
        if(alpha_action[[i]][j] == k){
          if(is.vector(initial)){
            vv = sum(initial*alpha_vector[[i]][[j]])
          } else{
            vv = sum(initial[i,]*alpha_vector[[i]][[j]])
          }
          V[[i]][k] = max(vv,V[[i]][k])
        }
      }
    }
  }
  v = zeros(1,size(V))
  p = zeros(1,size(V))
  for(i in 1:size(V)){
    v[i] = max(V[[i]])
    p[i] = which.max(V[[i]])
  }
  cc = min(v)
  ind = which.min(v)
  ExpectedV  = cc
  ExpecctedPol = p[ind]
  list(ExpectedV, ExpecctedPol)
}

fakemodel <- function(T,O,s_0,act,R,Num_s,Num_z){
  t = length(act)
  # transition at each time step
  v_s = zeros(Num_s,1)
  # emission at each time step
  v_z = zeros(Num_z,1)
  
  # outputs
  # states
  v_state = zeros(t,1)
  # observation
  v_obs = zeros(t,1)
  # reward 
  v_rew = zeros(t,1)
  
  s_t1 = s_0
  
  for(k in 1:t){
    a = act[k]
    v_s = T[s_t1, , a]
    
    # new state
    s_t2 = discrete(v_s,1)
    
    v_z = O[s_t2, , a]
    
    # new observation
    z_t2 = discrete(v_z,1)
    
    # recording
    v_state[k] = s_t2
    v_obs[k] = z_t2
    v_rew[k] = R[s_0,a]
    
    s_t1 = s_t2
  }
  list(v_state,v_obs,v_rew)
}


fakemodel_reality <- function(T,O,s_0,act,R,Num_s,Num_z,z){
  t = length(act)
  # transition at each time step
  v_s = zeros(Num_s,1)
  # emission at each time step
  v_z = zeros(Num_z,1)
  
  # outputs
  # states
  v_state = zeros(t,1)
  # observation
  v_obs = zeros(t,1)
  # reward 
  v_rew = zeros(t,1)
  
  s_t1 = s_0
  
  for(k in 1:t){
    a = act[k]
    v_s = T[s_t1, , a]
    
    # new state
    s_t2 = discrete_sample(v_s,1)
    
    #v_z = O[s_t2, , a]
    
    # new observation
    z_t2 = z
    
    # recording
    v_state[k] = s_t2
    v_obs[k] = z_t2
    v_rew[k] = R[s_0,a]
    
    s_t1 = s_t2
  }
  list(v_state,v_obs,v_rew)
}


update_belief <- function(b_0,T,O,a,z){
  b_1 = b_0 %*% T[,,a] * t(O[,z,a])
  b_1 = normalize(b_1)
  b_1
}

ff <- function(initial,T,O,a,z,P_0){
  b = zeros(length(a)+1,size(T,1))
  b[1,] = initial
  for(i in 1:length(a)){
    b[i+1,] = (b[i,] %*% T[,,a[i]]*t(O[,z[i],a[i]]))
    b[i+1,] = normalize(b[i+1,])
  }
  P = (b[length(a),] %*% T[,,a[length(a)]]) %*% O[,z[length(z)],a[length(a)]] %*% P_0
  b1 = b[length(a)+1,]
  
  list(P,b1)
}

OptConf <- function(PP,T,O,b_opt,tt,alpha_vector,alpha_action,cf,m){
  tmp=sort(PP,index.return=TRUE,decreasing = TRUE);
  PP_s=tmp$x; ind=tmp$ix
  cumP = zeros(1,length(PP_s))
  cumP = cumsum(PP_s)
  T_opt = list()
  O_opt = list()
  belief_upd_opt = list()
  av_opt = list()
  aa_opt = list()
  if((sum(PP > cf)) == 1 ){
    indd = which(PP > cf)
    T_opt[[1]] = T[[indd]]
    O_opt[[1]] = O[[indd]]
    belief_upd_opt = b_opt[[indd]][m,tt+1,]
    av_opt[[1]] = alpha_vector[[indd]]
    aa_opt[[1]] = alpha_action[[indd]]
  } else{
    k = 1
    while((1-cumP[k]) >= (1-cf)){
      T_opt[[k]] = T[[ind[k]]]
      O_opt[[k]] = O[[ind[k]]]
      k = k + 1
    }
    for(i in (k-1):(length(PP_s)-1)){
      if(PP_s[i+1] == PP_s[k-1]){
        T_opt[[i+1]] = T[[ind[i+1]]]
        O_opt[[i+1]] = O[[ind[i+1]]]
      }
    }
    belief_upd_opt = matrix(,nrow = size(T_opt), ncol = Num_s)
    for(i in 1:size(T_opt)){
      belief_upd_opt[i,] = b_opt[[ind[i]]][m,tt+1,]
      av_opt[[i]] = alpha_vector[[ind[i]]]
      aa_opt[[i]] = alpha_action[[ind[i]]]
    }
  }
  list(PP_s,cumP,T_opt,O_opt,av_opt,aa_opt,belief_upd_opt,ind)
}


###################################

pl_evaluator <- function(T,R,GAMMA,Pol,s,w,Q){
  V = array(0,dim = c(size(T),size(R,2)))
  for(i in 1:size(T)){
    #TT = T[[i]]; Policy = Pol[i,]
    for(j in 1:size(R,2)){
    #  Policy[s] = j
    #  out = mdp_eval_policy_iterative(TT,R,GAMMA,Policy)
    #  V[i,j] = out[s]
      V[i,j] = Q[i,s,j]
    }
  }
  V_action = array(0,dim = c(size(R,2)))
  w = normalize(w)
  for(i in 1:length(V_action)){
    for(j in 1:size(V,1)){
      V_action[i] = V_action[i] + (V[j,i]*w[j])
    }
   # V_action[i] = V_action[i]/size(V,1)
  }
  EP = which.max(V_action)
  EP
}


###########################

fakemodel_mdp <- function(T,s_0,act,R,Num_s){
  t = length(act)
  # transition at each time step
  v_s = zeros(Num_s,1)

  # outputs
  # states
  v_state = zeros(t,1)
  # reward 
  v_rew = zeros(t,1)
  
  s_t1 = s_0
  
  for(k in 1:t){
    a = act[k]
    v_s = T[s_t1, , a]
    
    # new state
    s_t2 = discrete(v_s,1)
    

    # recording
    v_state[k] = s_t2
    v_rew[k] = R[s_0,a]
    
    s_t1 = s_t2
  }
  list(v_state,v_rew)
}

###########################

ff_mdp <- function(s_0,T,a,s,P_0){
  P = T[s_0,s,a] * P_0
  P
}
###########################

discrete <- function(p,n){
  s = zeros(1,n)
  for (i in 1:n){
    s[i] = sum(cumsum(p)<runif(1))+1
  }
  s
}

###########################

Temp_tran <- function(Num_s,Num_a,sigma_m,noise_dist){
  
  O = zeros(Num_s,Num_s,Num_a)
  
  if(sigma_m == 0){
    for (j in 1:Num_a){
      O[ , ,j] = diag(Num_s)
    }
  } else if(noise_dist == "uniform"){
    
    for (i in 1:Num_s) {
      for (j in 1:Num_a) {
        m = floor(c((1-sigma_m),(1+sigma_m)) * i)
        nnn = max(1, m[1])
        mmm = min(Num_s, m[2])
        mm = 1 / (max(m) - nnn + 1)
        
        if (m[2] == 1) {
          O[i, nnn:mmm, j] = mm
        } else {
          O[i, nnn:(mmm - 1), j] = mm
          O[i, mmm, j] = 1 - ((mmm - nnn + 1) * mm)
        }
      }
    }
    
  } else if (noise_dist == "lognormal"){
    var_m_x = (1+sigma_m-(1-sigma_m)) / 12
    var_m_logn =  log(1+var_m_x)
    sigma_m_logn = sqrt(var_m_logn)
    mu_m_logn = log(1) - var_m_logn/2
    zz = zeros(Num_s, Num_s, Num_a)
    for (i in 1:Num_s){
      for (j in 1:Num_a){
        for (ii in 1:Num_s){
          if(i == 1){
            O[i,1,j] = 1
          } else{
            zz[i,ii,j] = ii / i
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
  for(i in 1:size(O,3)){
    O[,,i] = round_Milad(O[,,i])
  }
  
  output <- O
}

#############################################
runPOMDP_L_new <- function(Num_s,Num_a,Num_z,initial,T,O,R,GAMMA){
  
  out = list()
  value = list()
  policy = list()
  alpha_vector = list()
  alpha_action = list()
  actions = paste0("a", 1:Num_a)
  states = paste0("s", 1:Num_s)
  nnn = 1
  if(is.matrix(initial)){
    siz = size(initial,1)
  } else{
    siz = 1
  }
  
  for (i in 1:siz) {
    if(siz==1){
      belief = initial
    } else{
      belief = initial[i,]
    }
    belief = normalize(belief)
    belief = round(belief,4) / sum(round(belief,4))
    if(siz==1){
      milad_writepomdpx_POMDP(T,O,R,GAMMA,belief,Num_s,Num_a,Num_z,actions,states,nnn)
    } else{
      milad_writepomdpx_POMDP(T[[i]],O[[i]],R,GAMMA,belief,Num_s,Num_a,Num_z,actions,states,nnn=1)
    }
    command = paste0("/Users/Milad/Documents/appl-0.96/src/pomdpsol stationary",
                     nnn,".pomdpx --output statout",nnn,
                     ".policy --precision 1 --timeout 5")
    run_MOMDP(command)
    out = findpolicy_Milad2(belief,file = paste0("statout",nnn,".policy"))
    value[i] = out[[1]]
    policy[i] = out[[2]]
    alpha_vector[[i]] = out[[3]]
    alpha_action[[i]] = out[[4]]
  }
  
  
  list(value,policy,alpha_vector,alpha_action)
  
}

#############################################

discrete <- function(p,n){
  s = zeros(1,n)
  for (i in 1:n){
    s[i] = sum(cumsum(p)<runif(1))+1
  }
  s
}
