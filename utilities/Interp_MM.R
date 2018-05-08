Interp_MM = function(initial, alpha, alpha_action){
  library(purrr)
  
  ## Compute dot product with initial
  a <- alpha %>% map_dbl(function(x) initial %*% matrix(x, ncol=1))
  
  ## Return policy of the vector which has the biggest inner product
  #alpha_action[which.max(a)]
  if(sum(a==0)==length(a)){
    output = list(0,1)
  } else{
    output = list(max(a),alpha_action[which.max(a)])
  }
}
