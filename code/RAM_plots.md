    library(ggplot2)
    library(dplyr)

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    library(tidyr)
    source("../R/Fish_dynamic_Fun2.R")
    source("../R/fisheries_matrices.R")
    source("../R/findpolicy_Milad2.R")
    library(xml2)
    source("../R/Interp_MM.R")
    source("../R/L_fun.R")

Here we are loading the results of forward simulations and creating a
data frame of results for visualization

    load("../data/final_index.RData")
    load("../data/RAM_files3.RData")
    load("../data/RAM_data_us.RData")
    states = seq(0,2,by = 0.02)
    actions <- states
    obs <- states
    ### results of small measurement error

    load("../data/RAM4_small_pomdp_long.RData")
    load("../data/RAM4_small_MSY_real_long.RData")
    load("../data/RAM4_small_BAU_long.RData")
    load("../data/RAM4_small_MDP_long.RData")

    pomdp_proj_st2_s = vector("list", length = length(ind))
    pomdp_proj_rew2_s = vector("list", length = length(ind))
    mdp_proj_st2_s = vector("list", length = length(ind))
    mdp_proj_rew2_s = vector("list", length = length(ind))
    det_8_proj_st2_s = vector("list", length = length(ind))
    det_8_proj_rew2_s = vector("list", length = length(ind))
    pomdp_proj_act2_s = vector("list", length = length(ind))
    mdp_proj_act2_s = vector("list", length = length(ind))
    det_8_proj_act2_s = vector("list", length = length(ind))
    bau_proj_st2_s = vector("list", length = length(ind))
    bau_proj_act2_s = vector("list", length = length(ind))
    bau_proj_rew2_s = vector("list", length = length(ind))
    Ind_new_s = array(0, dim = length(ind))

    k = 1
    for(i in 1:length(ind)){
        pomdp_proj_st2_s[[i]] = pomdp_proj_st[[i]]
        pomdp_proj_rew2_s[[i]] = pomdp_proj_rew[[i]]
        mdp_proj_st2_s[[i]] = mdp_proj_st[[i]]
        mdp_proj_rew2_s[[i]] = mdp_proj_rew[[i]]
        det_8_proj_st2_s[[i]] = msy_proj_st[[i]]
        det_8_proj_rew2_s[[i]] = msy_proj_rew[[i]]
        pomdp_proj_act2_s[[i]] = pomdp_proj_act[[i]]
        mdp_proj_act2_s[[i]] = mdp_proj_act[[i]]
        det_8_proj_act2_s[[i]] = msy_proj_act[[i]]
        bau_proj_st2_s[[i]] = bau_proj_st[[i]]
        bau_proj_act2_s[[i]] = bau_proj_act[[i]]
        bau_proj_rew2_s[[i]] = bau_proj_rew[[i]]
        Ind_new_s[i] = ind[i]
    }

    ### results of large measurement error


    load("../data/RAM4_large_pomdp_long.RData")
    load("../data/RAM4_large_MSY_real_long.RData")
    load("../data/RAM4_large_BAU_long.RData")
    load("../data/RAM4_large_MDP_long.RData")

    pomdp_proj_st2_l = vector("list", length = length(ind))
    pomdp_proj_rew2_l = vector("list", length = length(ind))
    mdp_proj_st2_l = vector("list", length = length(ind))
    mdp_proj_rew2_l = vector("list", length = length(ind))
    det_8_proj_st2_l = vector("list", length = length(ind))
    det_8_proj_rew2_l = vector("list", length = length(ind))
    pomdp_proj_act2_l = vector("list", length = length(ind))
    mdp_proj_act2_l = vector("list", length = length(ind))
    det_8_proj_act2_l = vector("list", length = length(ind))
    bau_proj_st2_l = vector("list", length = length(ind))
    bau_proj_act2_l = vector("list", length = length(ind))
    bau_proj_rew2_l = vector("list", length = length(ind))
    Ind_new_l = array(0, dim = length(ind))

    k = 1
    for(i in 1:length(ind)){
      k = i
        pomdp_proj_st2_l[[k]] = pomdp_proj_st[[i]]
        pomdp_proj_rew2_l[[k]] = pomdp_proj_rew[[i]]
        mdp_proj_st2_l[[k]] = mdp_proj_st[[i]]
        mdp_proj_rew2_l[[k]] = mdp_proj_rew[[i]]
        det_8_proj_st2_l[[k]] = msy_proj_st[[i]]
        det_8_proj_rew2_l[[k]] = msy_proj_rew[[i]]
        pomdp_proj_act2_l[[k]] = pomdp_proj_act[[i]]
        mdp_proj_act2_l[[k]] = mdp_proj_act[[i]]
        det_8_proj_act2_l[[k]] = msy_proj_act[[i]]
        bau_proj_st2_l[[i]] = bau_proj_st[[i]]
        bau_proj_act2_l[[i]] = bau_proj_act[[i]]
        bau_proj_rew2_l[[i]] = bau_proj_rew[[i]]
        Ind_new_l[k] = ind[i]
           k = k + 1
    }


    ### results of medium measurement error


    load("../data/RAM4_medium_pomdp_long.RData")
    load("../data/RAM4_medium_MSY_real_long.RData")
    load("../data/RAM4_medium_BAU_long.RData")
    load("../data/RAM4_medium_MDP_long.RData")

    pomdp_proj_st2_m = vector("list", length = length(ind))
    pomdp_proj_rew2_m = vector("list", length = length(ind))
    mdp_proj_st2_m = vector("list", length = length(ind))
    mdp_proj_rew2_m = vector("list", length = length(ind))
    det_8_proj_st2_m = vector("list", length = length(ind))
    det_8_proj_rew2_m = vector("list", length = length(ind))
    pomdp_proj_act2_m = vector("list", length = length(ind))
    mdp_proj_act2_m = vector("list", length = length(ind))
    det_8_proj_act2_m = vector("list", length = length(ind))
    bau_proj_st2_m = vector("list", length = length(ind))
    bau_proj_act2_m = vector("list", length = length(ind))
    bau_proj_rew2_m = vector("list", length = length(ind))

    Ind_new_m = array(0, dim = length(ind))

    k = 1
    for(i in 1:length(ind)){
      k = i
        pomdp_proj_st2_m[[k]] = pomdp_proj_st[[i]]
        pomdp_proj_rew2_m[[k]] = pomdp_proj_rew[[i]]
        mdp_proj_st2_m[[k]] = mdp_proj_st[[i]]
        mdp_proj_rew2_m[[k]] = mdp_proj_rew[[i]]
        det_8_proj_st2_m[[k]] = msy_proj_st[[i]]
        det_8_proj_rew2_m[[k]] = msy_proj_rew[[i]]
        pomdp_proj_act2_m[[k]] = pomdp_proj_act[[i]]
        mdp_proj_act2_m[[k]] = mdp_proj_act[[i]]
        det_8_proj_act2_m[[k]] = msy_proj_act[[i]]
        bau_proj_st2_m[[i]] = bau_proj_st[[i]]
        bau_proj_act2_m[[i]] = bau_proj_act[[i]]
        bau_proj_rew2_m[[i]] = bau_proj_rew[[i]]
           Ind_new_m[k] = ind[i]
           k = k + 1
    }


    ### results of zero measurement error


    load("../data/RAM4_zero_MSY_real_long.RData")
    load("../data/RAM4_zero_MDP_long.RData")
    load("../data/RAM4_zero_BAU_long.RData")

    pomdp_proj_st2_z = vector("list", length = length(ind))
    pomdp_proj_rew2_z = vector("list", length = length(ind))
    mdp_proj_st2_z = vector("list", length = length(ind))
    mdp_proj_rew2_z = vector("list", length = length(ind))
    det_8_proj_st2_z = vector("list", length = length(ind))
    det_8_proj_rew2_z = vector("list", length = length(ind))
    pomdp_proj_act2_z = vector("list", length = length(ind))
    mdp_proj_act2_z = vector("list", length = length(ind))
    det_8_proj_act2_z = vector("list", length = length(ind))
    bau_proj_st2_z = vector("list", length = length(ind))
    bau_proj_act2_z = vector("list", length = length(ind))
    bau_proj_rew2_z = vector("list", length = length(ind))

    Ind_new_z = array(0, dim = length(ind))

    k = 1
    for(i in 1:length(ind)){
      k = i
        pomdp_proj_st2_z[[k]] = mdp_proj_st[[i]]
        pomdp_proj_rew2_z[[k]] = mdp_proj_rew[[i]]
        mdp_proj_st2_z[[k]] = mdp_proj_st[[i]]
        mdp_proj_rew2_z[[k]] = mdp_proj_rew[[i]]
        det_8_proj_st2_z[[k]] = msy_proj_st[[i]]
        det_8_proj_rew2_z[[k]] = msy_proj_rew[[i]]
        pomdp_proj_act2_z[[k]] = mdp_proj_act[[i]]
        mdp_proj_act2_z[[k]] = mdp_proj_act[[i]]
        det_8_proj_act2_z[[k]] = msy_proj_act[[i]]
        bau_proj_st2_z[[i]] = bau_proj_st[[i]]
        bau_proj_act2_z[[i]] = bau_proj_act[[i]]
        bau_proj_rew2_z[[i]] = bau_proj_rew[[i]]
           Ind_new_z[k] = ind[i]
           k = k + 1
    }




    ### Just some reshaping and calibration of data

    pomdp_proj_st_l = vector("list", length = length(pomdp_proj_st2_l))
    mdp_proj_st_l = vector("list", length = length(pomdp_proj_st2_l))
    det_8_proj_st_l = vector("list", length = length(pomdp_proj_st2_l))
    bau_proj_st_l = vector("list", length = length(pomdp_proj_st2_l))
    pomdp_proj_st_s = vector("list", length = length(pomdp_proj_st2_l))
    mdp_proj_st_s = vector("list", length = length(pomdp_proj_st2_l))
    det_8_proj_st_s = vector("list", length = length(pomdp_proj_st2_l))
    bau_proj_st_s = vector("list", length = length(pomdp_proj_st2_l))
    pomdp_proj_st_m = vector("list", length = length(pomdp_proj_st2_l))
    mdp_proj_st_m = vector("list", length = length(pomdp_proj_st2_l))
    det_8_proj_st_m = vector("list", length = length(pomdp_proj_st2_l))
    bau_proj_st_m = vector("list", length = length(pomdp_proj_st2_l))
    pomdp_proj_st_z = vector("list", length = length(pomdp_proj_st2_l))
    mdp_proj_st_z = vector("list", length = length(pomdp_proj_st2_l))
    det_8_proj_st_z = vector("list", length = length(pomdp_proj_st2_l))
    bau_proj_st_z = vector("list", length = length(pomdp_proj_st2_l))

    for(i in 1:length(pomdp_proj_st2_l)){
      k = i
      pomdp_proj_st_l[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
        mdp_proj_st_l[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
        det_8_proj_st_l[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
        bau_proj_st_l[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         pomdp_proj_st_s[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         mdp_proj_st_s[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         det_8_proj_st_s[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         bau_proj_st_s[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         pomdp_proj_st_m[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         mdp_proj_st_m[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         det_8_proj_st_m[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         bau_proj_st_m[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         pomdp_proj_st_z[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         mdp_proj_st_z[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         det_8_proj_st_z[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))
         bau_proj_st_z[[k]] = array(0, dim = dim(pomdp_proj_st2_l[[i]]))

        for(j in 1:dim(pomdp_proj_st2_l[[i]])[1]){
          pomdp_proj_st_l[[k]][j,] = states[pomdp_proj_st2_l[[i]][j,]]
          mdp_proj_st_l[[k]][j,] = states[mdp_proj_st2_l[[i]][j,]]
          det_8_proj_st_l[[k]][j,] = states[det_8_proj_st2_l[[i]][j,]]
          bau_proj_st_l[[k]][j,] = states[bau_proj_st2_l[[i]][j,]]
          
           pomdp_proj_st_s[[k]][j,] = states[pomdp_proj_st2_s[[i]][j,]]
          mdp_proj_st_s[[k]][j,] = states[mdp_proj_st2_s[[i]][j,]]
           det_8_proj_st_s[[k]][j,] = states[det_8_proj_st2_s[[i]][j,]]
           bau_proj_st_s[[k]][j,] = states[bau_proj_st2_s[[i]][j,]]
          # 
           pomdp_proj_st_m[[k]][j,] = states[pomdp_proj_st2_m[[i]][j,]]
           mdp_proj_st_m[[k]][j,] = states[mdp_proj_st2_m[[i]][j,]]
           det_8_proj_st_m[[k]][j,] = states[det_8_proj_st2_m[[i]][j,]]
           bau_proj_st_m[[k]][j,] = states[bau_proj_st2_m[[i]][j,]]
          # 
           pomdp_proj_st_z[[k]][j,] = states[pomdp_proj_st2_z[[i]][j,]]
           mdp_proj_st_z[[k]][j,] = states[mdp_proj_st2_z[[i]][j,]]
           det_8_proj_st_z[[k]][j,] = states[det_8_proj_st2_z[[i]][j,]]
           bau_proj_st_z[[k]][j,] = states[bau_proj_st2_z[[i]][j,]]
        }
    }

    ## calculating the ratio of biomass /B_MSY

    bb_pomdp_large = vector("list", length = length(pomdp_proj_st2_l))
    bb_pomdp_small = vector("list", length = length(pomdp_proj_st2_l))
    bb_pomdp_m = vector("list", length = length(pomdp_proj_st2_l))
    bb_pomdp_z = vector("list", length = length(pomdp_proj_st2_l))
    bb_mdp_large = vector("list", length = length(pomdp_proj_st2_l))
    bb_mdp_small = vector("list", length = length(pomdp_proj_st2_l))
    bb_mdp_m = vector("list", length = length(pomdp_proj_st2_l))
    bb_mdp_z = vector("list", length = length(pomdp_proj_st2_l))
    bb_msy_large = vector("list", length = length(pomdp_proj_st2_l))
    bb_msy_small = vector("list", length = length(pomdp_proj_st2_l))
    bb_msy_m = vector("list", length = length(pomdp_proj_st2_l))
    bb_msy_z = vector("list", length = length(pomdp_proj_st2_l))
    bb_bau_large = vector("list", length = length(pomdp_proj_st2_l))
    bb_bau_small = vector("list", length = length(pomdp_proj_st2_l))
    bb_bau_m = vector("list", length = length(pomdp_proj_st2_l))
    bb_bau_z = vector("list", length = length(pomdp_proj_st2_l))

    k=1
    for(i in 1:length(pomdp_proj_st2_l)){
      k = i  
      n = Ind_new_l[i]
        S_star = mean_K[n]/2
        B_MSY <- logistic(S_star,0,mean_r[n],mean_K[n])
        bb_pomdp_large[[k]] = pomdp_proj_st_l[[k]] / B_MSY
        bb_pomdp_small[[k]] = pomdp_proj_st_s[[k]] / B_MSY
        bb_pomdp_m[[k]] = pomdp_proj_st_m[[k]] / B_MSY
        bb_pomdp_z[[k]] = pomdp_proj_st_z[[k]] / B_MSY
        bb_mdp_large[[k]] = mdp_proj_st_l[[k]] / B_MSY
        bb_mdp_small[[k]] = mdp_proj_st_s[[k]] / B_MSY
        bb_mdp_m[[k]] = mdp_proj_st_m[[k]] / B_MSY
        bb_mdp_z[[k]] = mdp_proj_st_z[[k]] / B_MSY
        bb_msy_large[[k]] = det_8_proj_st_l[[k]] / B_MSY
        bb_msy_small[[k]] = det_8_proj_st_s[[k]] / B_MSY
        bb_msy_m[[k]] = det_8_proj_st_m[[k]] / B_MSY
        bb_msy_z[[k]] = det_8_proj_st_z[[k]] / B_MSY
        bb_bau_large[[k]] = bau_proj_st_l[[k]] / B_MSY
        bb_bau_small[[k]] = bau_proj_st_s[[k]] / B_MSY
        bb_bau_m[[k]] = bau_proj_st_m[[k]] / B_MSY
        bb_bau_z[[k]] = bau_proj_st_z[[k]] / B_MSY
    }

    ## claculating the percentage of stock above 80% B/B_MSY


    ind_pomdp_large = vector("list", length = length(pomdp_proj_st_l))
    ind_pomdp_small = vector("list", length = length(pomdp_proj_st_l))
    ind_pomdp_m = vector("list", length = length(pomdp_proj_st_l))
    ind_pomdp_z = vector("list", length = length(pomdp_proj_st_l))
    ind_mdp_large = vector("list", length = length(pomdp_proj_st_l))
    ind_mdp_small = vector("list", length = length(pomdp_proj_st_l))
    ind_mdp_m = vector("list", length = length(pomdp_proj_st_l))
    ind_mdp_z = vector("list", length = length(pomdp_proj_st_l))
    ind_msy_large = vector("list", length = length(pomdp_proj_st_l))
    ind_msy_small = vector("list", length = length(pomdp_proj_st_l))
    ind_msy_m = vector("list", length = length(pomdp_proj_st_l))
    ind_msy_z = vector("list", length = length(pomdp_proj_st_l))
    ind_bau_large = vector("list", length = length(pomdp_proj_st_l))
    ind_bau_small = vector("list", length = length(pomdp_proj_st_l))
    ind_bau_m = vector("list", length = length(pomdp_proj_st_l))
    ind_bau_z = vector("list", length = length(pomdp_proj_st_l))

    thr = 0.8
    for(i in 1:length(bb_pomdp_large)){
      ind_pomdp_large[[i]] = bb_pomdp_large[[i]] >=thr
      ind_pomdp_small[[i]] = bb_pomdp_small[[i]] >= thr
      ind_pomdp_m[[i]] = bb_pomdp_m[[i]] >= thr
      ind_pomdp_z[[i]] = bb_pomdp_z[[i]] >= thr
      # 
      ind_mdp_large[[i]] = bb_mdp_large[[i]] >= thr
      ind_mdp_small[[i]] = bb_mdp_small[[i]] >= thr
      ind_mdp_m[[i]] = bb_mdp_m[[i]] >= thr
      ind_mdp_z[[i]] = bb_mdp_z[[i]] >= thr
      
      ind_msy_large[[i]] = bb_msy_large[[i]] >= thr
      ind_msy_small[[i]] = bb_msy_small[[i]] >=thr
      ind_msy_m[[i]] = bb_msy_m[[i]] >= thr
      ind_msy_z[[i]] = bb_msy_z[[i]] >= thr
      
      ind_bau_large[[i]] = bb_bau_large[[i]] >= thr
      ind_bau_small[[i]] = bb_bau_small[[i]] >= thr
      ind_bau_m[[i]] = bb_bau_m[[i]] >= thr
      ind_bau_z[[i]] = bb_bau_z[[i]] >= thr
    }

    Num_sim = 500; Tmax = 500
      
    find_pomdp_large = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_pomdp_small = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_pomdp_m = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_pomdp_z = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_mdp_large = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_mdp_small = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_mdp_m = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_mdp_z = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_msy_large = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_msy_small = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_msy_m = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_msy_z = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_bau_large = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_bau_small = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_bau_m = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))
    find_bau_z = array(0, dim = c(length(pomdp_proj_st_l), Num_sim, Tmax))

    for(i in 1:length(pomdp_proj_st_l)){
      find_pomdp_large[i,,] = ind_pomdp_large[[i]]
      find_pomdp_small[i,,] = ind_pomdp_small[[i]]
      find_pomdp_m[i,,] = ind_pomdp_m[[i]]
      find_pomdp_z[i,,] = ind_pomdp_z[[i]]
      find_mdp_large[i,,] = ind_mdp_large[[i]]
      find_mdp_small[i,,] = ind_mdp_small[[i]]
      find_mdp_m[i,,] = ind_mdp_m[[i]]
      find_mdp_z[i,,] = ind_mdp_z[[i]]
      find_msy_large[i,,] = ind_msy_large[[i]]
      find_msy_small[i,,] = ind_msy_small[[i]]
      find_msy_m[i,,] = ind_msy_m[[i]]
      find_msy_z[i,,] = ind_msy_z[[i]]
      find_bau_large[i,,] = ind_bau_large[[i]]
      find_bau_small[i,,] = ind_bau_small[[i]]
      find_bau_m[i,,] = ind_bau_m[[i]]
      find_bau_z[i,,] = ind_bau_z[[i]]
     }


    ## calculating the average over the stocks

    ff_pomdp_large = apply(find_pomdp_large,c(2,3),mean)
    ff_pomdp_small = apply(find_pomdp_small,c(2,3),mean)
    ff_pomdp_m = apply(find_pomdp_m,c(2,3),mean)
    ff_pomdp_z = apply(find_pomdp_z,c(2,3),mean)
    ff_mdp_large = apply(find_mdp_large,c(2,3),mean)
    ff_mdp_small = apply(find_mdp_small,c(2,3),mean)
    ff_mdp_m = apply(find_mdp_m,c(2,3),mean)
    ff_mdp_z = apply(find_mdp_z,c(2,3),mean)
    ff_msy_large = apply(find_msy_large,c(2,3),mean)
    ff_msy_small = apply(find_msy_small,c(2,3),mean)
    ff_msy_m = apply(find_msy_m,c(2,3),mean)
    ff_msy_z = apply(find_msy_z,c(2,3),mean)
    ff_bau_large = apply(find_bau_large,c(2,3),mean)
    ff_bau_small = apply(find_bau_small,c(2,3),mean)
    ff_bau_m = apply(find_bau_m,c(2,3),mean)
    ff_bau_z = apply(find_bau_z,c(2,3),mean)



    ## Finding the percentage of stock above B/B_MSY for historical data
    Ind = 1:length(data)
    N = array(0, dim = length(Ind))
    y = vector("list", length = length(Ind))
    a = vector("list", length = length(Ind))
    scale = array(0, dim = length(Ind))
    name = character(length = length(Ind))
    for(ii in 1:length(Ind)){
      post <- data[[ii]]
      name[ii] = post$commmonname
      N[ii] <- length(post$year)
      scale[ii] <- max(post$biomass)
      scaled_data <- data.frame(t = 1:N[ii], y = post$biomass / scale[ii], 
                                a = post$catch / scale[ii])
      y[[ii]] <- sapply(scaled_data$y, function(y) which.min(abs(states - y)))
      a[[ii]] <- sapply(scaled_data$a, function(a) which.min(abs(actions - a)))
    }

    find_biomass = array(0, dim = c(length(ind), 28))
    low = 1980; high = 2007
    k = 1
    for(i in 1:length(ind)){
      n = ind[i]
      post = data[[n]]
      year = post$year
      S_star = mean_K[n]/2
      B_MSY <- logistic(S_star,0,mean_r[n],mean_K[n])
      nom = post$year[post$year %in% seq(low,high)] - (low-1)
      find_biomass[k,nom] = (post$biomass[post$year %in% seq(low,high)] / scale[n]) / B_MSY*100#(mean_K[n]/2) * 100
      k = k+1
    }

    cc = apply(find_biomass == 0, 2,sum)
    ind_biomass = find_biomass >= 80
    ff_biomass = array(0, dim = dim(find_biomass)[2])
    for(i in 1:dim(find_biomass)[2]){
      ff_biomass[i] = sum(ind_biomass[,i]) / (dim(find_biomass)[1] - cc[i]) 
    }

here we create a data frame of the results (no need to run, results are
saved in the data folder)

    ## Creating data frame

    low = 1980; high = 2007
    years = seq(low,high)
    Tmax = 500; Num_sim = 500
    t = c(years,seq(years[length(years)]+1,(years[length(years)]+ Tmax-1)))
    model = c("MDP","POMDP", "History","MSY","BAU")
    type = c("0% Error","5% Error","10% Error","20% Error")
    qqq = expand.grid(t,1:Num_sim, model,type)
    traj = array(0, dim = dim(qqq)[1])

    biomass_pomdp_l = array(0,dim = c(Num_sim,length(t)))
    biomass_pomdp_s = array(0,dim = c(Num_sim,length(t)))
    biomass_pomdp_m = array(0,dim = c(Num_sim,length(t)))
    biomass_pomdp_z = array(0,dim = c(Num_sim,length(t)))

    biomass_mdp_l = array(0,dim = c(Num_sim,length(t)))
    biomass_mdp_s = array(0,dim = c(Num_sim,length(t)))
    biomass_mdp_m = array(0,dim = c(Num_sim,length(t)))
    biomass_mdp_z = array(0,dim = c(Num_sim,length(t)))

    biomass_MSY_l = array(0,dim = c(Num_sim,length(t)))
    biomass_MSY_s = array(0,dim = c(Num_sim,length(t)))
    biomass_MSY_m = array(0,dim = c(Num_sim,length(t)))
    biomass_MSY_z = array(0,dim = c(Num_sim,length(t)))

    biomass_BAU_l = array(0,dim = c(Num_sim,length(t)))
    biomass_BAU_s = array(0,dim = c(Num_sim,length(t)))
    biomass_BAU_m = array(0,dim = c(Num_sim,length(t)))
    biomass_BAU_z = array(0,dim = c(Num_sim,length(t)))


    biomass = array(0,dim = c(Num_sim,length(t)))
    for(j in 1:Num_sim){
      biomass[j,1:28] = ff_biomass * 100
    }

        mm = length(low:high)+1; nn = length(low:high)+499 
        biomass_pomdp_l[,(mm:nn)] = ff_pomdp_large[,1:499] * 100#biomass_pomdp_l[,(mm:nn)] + pomdp_proj_st[[i]][,1:99] * scale[n]
        biomass_pomdp_s[,(mm:nn)] = ff_pomdp_small[,1:499] * 100
        biomass_pomdp_m[,(mm:nn)] = ff_pomdp_m[,1:499] * 100
        biomass_pomdp_z[,(mm:nn)] = ff_pomdp_z[,1:499] * 100
      
        biomass_mdp_l[,(mm:nn)] = ff_mdp_large[,1:499] * 100#biomass_pomdp_l[,(mm:nn)] + pomdp_proj_st[[i]][,1:99] * scale[n]
        biomass_mdp_s[,(mm:nn)] = ff_mdp_small[,1:499] * 100
        biomass_mdp_m[,(mm:nn)] = ff_mdp_m[,1:499] * 100
        biomass_mdp_z[,(mm:nn)] = ff_mdp_z[,1:499] * 100
     
        biomass_MSY_l[,(mm:nn)] = ff_msy_large[,1:499] * 100#biomass_MSY[,(mm:nn)] + det_8_proj_st[[i]][,1:99] * scale[n]
        biomass_MSY_s[,(mm:nn)] = ff_msy_small[,1:499] * 100
        biomass_MSY_m[,(mm:nn)] = ff_msy_m[,1:499] * 100
        biomass_MSY_z[,(mm:nn)] = ff_msy_z[,1:499] * 100
        
        biomass_BAU_l[,(mm:nn)] = ff_bau_large[,1:499] * 100
        biomass_BAU_s[,(mm:nn)] = ff_bau_small[,1:499] * 100
        biomass_BAU_m[,(mm:nn)] = ff_bau_m[,1:499] * 100
        biomass_BAU_z[,(mm:nn)] = ff_bau_z[,1:499] * 100
      #}
    #}

    biomass_pomdp_l[,length(low:high)] = biomass[,length(low:high)]
    biomass_mdp_l[,length(low:high)] = biomass[,length(low:high)]
    biomass_pomdp_m[,length(low:high)] = biomass[,length(low:high)]
    biomass_pomdp_z[,length(low:high)] = biomass[,length(low:high)]
    biomass_mdp_m[,length(low:high)] = biomass[,length(low:high)]
    biomass_mdp_z[,length(low:high)] = biomass[,length(low:high)]
    # biomass_mdp[,length(low:high)] = biomass[,length(low:high)]
    biomass_MSY_l[,length(low:high)] = biomass[,length(low:high)]
    biomass_MSY_m[,length(low:high)] = biomass[,length(low:high)]
    biomass_MSY_z[,length(low:high)] = biomass[,length(low:high)]
    biomass_pomdp_s[,length(low:high)] = biomass[,length(low:high)]
    biomass_mdp_s[,length(low:high)] = biomass[,length(low:high)]
    # biomass_mdp[,length(low:high)] = biomass[,length(low:high)]
    biomass_MSY_s[,length(low:high)] = biomass[,length(low:high)]
    biomass_BAU_l[,length(low:high)] = biomass[,length(low:high)]
    biomass_BAU_s[,length(low:high)] = biomass[,length(low:high)]
    biomass_BAU_m[,length(low:high)] = biomass[,length(low:high)]
    biomass_BAU_z[,length(low:high)] = biomass[,length(low:high)]

    for(i in 1:Num_sim){
      traj[(qqq$Var2 == i & qqq$Var3 == "History")] = biomass[i,]

      traj[(qqq$Var2 == i & qqq$Var3 == "POMDP" & qqq$Var4 == "20% Error")] = biomass_pomdp_l[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "POMDP" & qqq$Var4 == "5% Error")] = biomass_pomdp_s[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "POMDP" & qqq$Var4 == "0% Error")] = biomass_pomdp_z[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "POMDP" & qqq$Var4 == "10% Error")] = biomass_pomdp_m[i,]

      traj[(qqq$Var2 == i & qqq$Var3 == "MDP" & qqq$Var4 == "20% Error")] = biomass_mdp_l[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "MDP" & qqq$Var4 == "5% Error")] = biomass_mdp_s[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "MDP" & qqq$Var4 == "0% Error")] = biomass_mdp_z[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "MDP" & qqq$Var4 == "10% Error")] = biomass_mdp_m[i,]

      traj[(qqq$Var2 == i & qqq$Var3 == "MSY" & qqq$Var4 == "20% Error")] = biomass_MSY_l[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "MSY" & qqq$Var4 == "5% Error")] = biomass_MSY_s[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "MSY" & qqq$Var4 == "10% Error")] = biomass_MSY_m[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "MSY" & qqq$Var4 == "0% Error")] = biomass_MSY_z[i,]
      
      
      traj[(qqq$Var2 == i & qqq$Var3 == "BAU" & qqq$Var4 == "20% Error")] = biomass_BAU_l[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "BAU" & qqq$Var4 == "5% Error")] = biomass_BAU_s[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "BAU" & qqq$Var4 == "10% Error")] = biomass_BAU_m[i,]
      traj[(qqq$Var2 == i & qqq$Var3 == "BAU" & qqq$Var4 == "0% Error")] = biomass_BAU_z[i,]
    }

    qqq = data.frame(qqq,traj)

    qqq$traj[qqq$Var3 == "History" &  qqq$Var1 %in% seq((high+1),(high+99))] = NA
    qqq$traj[qqq$Var3 == "POMDP" &  qqq$Var1 %in% seq(low,(high-1))] = NA
    qqq$traj[qqq$Var3 == "MDP" &  qqq$Var1 %in% seq(low,(high-1))] = NA
    qqq$traj[qqq$Var3 == "MSY" &  qqq$Var1 %in% seq(low,(high-1))] = NA
    qqq$traj[qqq$Var3 == "BAU" &  qqq$Var1 %in% seq(low,(high-1))] = NA


    save(qqq, file = "../data/fwd_sim_results.RData")




    GAMMA = 0.99
    gam = GAMMA ^ seq(0,499)

    load("../data/RAM4_large_pomdp_long.RData")
    load("../data/RAM4_large_MDP_long.RData")
    load("../data/RAM4_large_BAU_long.RData")
    load("../data/RAM4_large_MSY_real_long.RData")

    pomdp_large = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    mdp_large = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    msy_large = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    bau_large = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))


    for(i in 1:length(pomdp_proj_st)){
      pomdp_large[i,,] = pomdp_proj_st[[i]]
      mdp_large[i,,] = mdp_proj_st[[i]]
      msy_large[i,,] = msy_proj_st[[i]]
      bau_large[i,,] = bau_proj_st[[i]]
    }



    load("../data/RAM4_medium_pomdp_long.RData")
    load("../data/RAM4_medium_MDP_long.RData")
    load("../data/RAM4_medium_BAU_long.RData")
    load("../data/RAM4_medium_MSY_real_long.RData")

    pomdp_medium = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    mdp_medium = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    msy_medium = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    bau_medium = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))


    for(i in 1:length(pomdp_proj_st)){
      pomdp_medium[i,,] = pomdp_proj_st[[i]]
      mdp_medium[i,,] = mdp_proj_st[[i]]
      msy_medium[i,,] = msy_proj_st[[i]]
      bau_medium[i,,] = bau_proj_st[[i]]
    }





    load("../data/RAM4_small_pomdp_long.RData")
    load("../data/RAM4_small_MDP_long.RData")
    load("../data/RAM4_small_BAU_long.RData")
    load("../data/RAM4_small_MSY_real_long.RData")

    pomdp_small = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    mdp_small = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    msy_small = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    bau_small = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))


    for(i in 1:length(pomdp_proj_st)){
      pomdp_small[i,,] = pomdp_proj_st[[i]]
      mdp_small[i,,] = mdp_proj_st[[i]]
      msy_small[i,,] = msy_proj_st[[i]]
      bau_small[i,,] = bau_proj_st[[i]]
    }





    load("/Users/Milad/Documents/UC Berkeley/Fisheries/R/Ram_fit/RAM4_zero_MDP_long.RData")
    load("/Users/Milad/Documents/UC Berkeley/Fisheries/R/Ram_fit/RAM4_zero_BAU_long.RData")
    load("/Users/Milad/Documents/UC Berkeley/Fisheries/R/Ram_fit/RAM4_zero_MSY_real_long.RData")

    pomdp_zero = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    mdp_zero = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    msy_zero = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))
    bau_zero = array(0, dim = c(length(pomdp_proj_st),Num_sim, Tmax))


    for(i in 1:length(pomdp_proj_st)){
      pomdp_zero[i,,] = mdp_proj_st[[i]]
      mdp_zero[i,,] = mdp_proj_st[[i]]
      msy_zero[i,,] = msy_proj_st[[i]]
      bau_zero[i,,] = bau_proj_st[[i]]
    }

    value_pomdp_large = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_pomdp_medium = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_pomdp_small = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_pomdp_zero = array(0, dim = c(length(pomdp_proj_st),Num_sim))


    value_mdp_large = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_mdp_medium = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_mdp_small = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_mdp_zero = array(0, dim = c(length(pomdp_proj_st),Num_sim))


    value_msy_large = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_msy_medium = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_msy_small = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_msy_zero = array(0, dim = c(length(pomdp_proj_st),Num_sim))


    value_bau_large = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_bau_medium = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_bau_small = array(0, dim = c(length(pomdp_proj_st),Num_sim))
    value_bau_zero = array(0, dim = c(length(pomdp_proj_st),Num_sim))

    scale_2 = scale[ind]


    for(i in 1:dim(pomdp_large)[1]){
      for(j in 1:dim(pomdp_large)[2]){
        value_pomdp_large[i,j] = sum(gam * pomdp_large[i,j,]) * scale[i]
        value_pomdp_medium[i,j] = sum(gam * pomdp_medium[i,j,]) * scale[i]
        value_pomdp_small[i,j] = sum(gam * pomdp_small[i,j,]) * scale[i]
        value_pomdp_zero[i,j] = sum(gam * pomdp_zero[i,j,]) * scale[i]
        
        value_mdp_large[i,j] = sum(gam * mdp_large[i,j,]) * scale[i]
        value_mdp_medium[i,j] = sum(gam * mdp_medium[i,j,]) * scale[i]
        value_mdp_small[i,j] = sum(gam * mdp_small[i,j,]) * scale[i]
        value_mdp_zero[i,j] = sum(gam * mdp_zero[i,j,]) * scale[i]
        
        value_msy_large[i,j] = sum(gam * msy_large[i,j,]) * scale[i]
        value_msy_medium[i,j] = sum(gam * msy_medium[i,j,]) * scale[i]
        value_msy_small[i,j] = sum(gam * msy_small[i,j,]) * scale[i]
        value_msy_zero[i,j] = sum(gam * msy_zero[i,j,]) * scale[i]
        
        value_bau_large[i,j] = sum(gam * bau_large[i,j,]) * scale[i]
        value_bau_medium[i,j] = sum(gam * bau_medium[i,j,]) * scale[i]
        value_bau_small[i,j] = sum(gam * bau_small[i,j,]) * scale[i]
        value_bau_zero[i,j] = sum(gam * bau_zero[i,j,]) * scale[i]

      }
    }


    total_pomdp_large = apply(value_pomdp_large,2,sum)
    total_pomdp_medium = apply(value_pomdp_medium,2,sum)
    total_pomdp_small = apply(value_pomdp_small,2,sum)
    total_pomdp_zero = apply(value_pomdp_small,2,sum)

    total_mdp_large = apply(value_mdp_large,2,sum)
    total_mdp_medium = apply(value_mdp_medium,2,sum)
    total_mdp_small = apply(value_mdp_small,2,sum)
    total_mdp_zero = apply(value_pomdp_small,2,sum)

    total_msy_large = apply(value_msy_large,2,sum)
    total_msy_medium = apply(value_msy_medium,2,sum)
    total_msy_small = apply(value_msy_small,2,sum)
    total_msy_zero = apply(value_msy_zero,2,sum)

    total_bau_large = apply(value_bau_large,2,sum)
    total_bau_medium = apply(value_bau_medium,2,sum)
    total_bau_small = apply(value_bau_small,2,sum)
    total_bau_zero = apply(value_bau_zero,2,sum)

    cc = mean(total_pomdp_large)

    total_pomdp_large = total_pomdp_large / cc
    total_mdp_large = total_mdp_large / cc
    total_msy_large = total_msy_large / cc
    total_bau_large = total_bau_large / cc

    cc = mean(total_pomdp_medium)

    total_pomdp_medium = total_pomdp_medium / cc
    total_mdp_medium = total_mdp_medium / cc
    total_msy_medium = total_msy_medium / cc
    total_bau_medium = total_bau_medium / cc


    cc = mean(total_pomdp_small)

    total_pomdp_small = total_pomdp_small / cc
    total_mdp_small = total_mdp_small /cc 
    total_msy_small = total_msy_small / cc
    total_bau_small = total_bau_small / cc

    cc = mean(total_pomdp_zero)

    total_pomdp_zero = total_pomdp_zero / cc
    total_mdp_zero = total_mdp_zero / cc
    total_msy_zero = total_msy_zero / cc
    total_bau_zero = total_bau_zero /cc


    model = c("POMDP","MDP","MSY","BAU")
    error = c("0% Error","5% Error","10% Error","20% Error")
    rew_d = expand.grid(model,error)

    mean = array(0, dim = dim(rew_d)[1])
    sd = array(0, dim = dim(rew_d)[1])

    mean[rew_d$Var1 == "POMDP" & rew_d$Var2 == "0% Error"] = mean(total_pomdp_small)
    mean[rew_d$Var1 == "MDP" & rew_d$Var2 == "0% Error"] = mean(total_pomdp_small)
    mean[rew_d$Var1 == "MSY" & rew_d$Var2 == "0% Error"] = mean(total_msy_zero)
    mean[rew_d$Var1 == "BAU" & rew_d$Var2 == "0% Error"] = mean(total_bau_zero)

    mean[rew_d$Var1 == "POMDP" & rew_d$Var2 == "5% Error"] = mean(total_pomdp_small)
    mean[rew_d$Var1 == "MDP" & rew_d$Var2 == "5% Error"] = mean(total_mdp_small)
    mean[rew_d$Var1 == "MSY" & rew_d$Var2 == "5% Error"] = mean(total_msy_small)
    mean[rew_d$Var1 == "BAU" & rew_d$Var2 == "5% Error"] = mean(total_bau_small)

    mean[rew_d$Var1 == "POMDP" & rew_d$Var2 == "10% Error"] = mean(total_pomdp_medium)
    mean[rew_d$Var1 == "MDP" & rew_d$Var2 == "10% Error"] = mean(total_mdp_medium)
    mean[rew_d$Var1 == "MSY" & rew_d$Var2 == "10% Error"] = mean(total_msy_medium)
    mean[rew_d$Var1 == "BAU" & rew_d$Var2 == "10% Error"] = mean(total_bau_medium)

    mean[rew_d$Var1 == "POMDP" & rew_d$Var2 == "20% Error"] = mean(total_pomdp_large)
    mean[rew_d$Var1 == "MDP" & rew_d$Var2 == "20% Error"] = mean(total_mdp_large)
    mean[rew_d$Var1 == "MSY" & rew_d$Var2 == "20% Error"] = mean(total_msy_large)
    mean[rew_d$Var1 == "BAU" & rew_d$Var2 == "20% Error"] = mean(total_bau_large)


    sd[rew_d$Var1 == "POMDP" & rew_d$Var2 == "0% Error"] = sd(total_pomdp_zero)
    sd[rew_d$Var1 == "MDP" & rew_d$Var2 == "0% Error"] = sd(total_mdp_zero)
    sd[rew_d$Var1 == "MSY" & rew_d$Var2 == "0% Error"] = sd(total_msy_zero)
    sd[rew_d$Var1 == "BAU" & rew_d$Var2 == "0% Error"] = sd(total_bau_zero)

    sd[rew_d$Var1 == "POMDP" & rew_d$Var2 == "5% Error"] = sd(total_pomdp_small)
    sd[rew_d$Var1 == "MDP" & rew_d$Var2 == "5% Error"] = sd(total_mdp_small)
    sd[rew_d$Var1 == "MSY" & rew_d$Var2 == "5% Error"] = sd(total_msy_small)
    sd[rew_d$Var1 == "BAU" & rew_d$Var2 == "5% Error"] = sd(total_bau_small)

    sd[rew_d$Var1 == "POMDP" & rew_d$Var2 == "10% Error"] = sd(total_pomdp_medium)
    sd[rew_d$Var1 == "MDP" & rew_d$Var2 == "10% Error"] = sd(total_mdp_medium)
    sd[rew_d$Var1 == "MSY" & rew_d$Var2 == "10% Error"] = sd(total_msy_medium)
    sd[rew_d$Var1 == "BAU" & rew_d$Var2 == "10% Error"] = sd(total_bau_medium)

    sd[rew_d$Var1 == "POMDP" & rew_d$Var2 == "20% Error"] = sd(total_pomdp_large)
    sd[rew_d$Var1 == "MDP" & rew_d$Var2 == "20% Error"] = sd(total_mdp_large)
    sd[rew_d$Var1 == "MSY" & rew_d$Var2 == "20% Error"] = sd(total_msy_large)
    sd[rew_d$Var1 == "BAU" & rew_d$Var2 == "20% Error"] = sd(total_bau_large)

    rew_d = data.frame(rew_d,mean,sd)

    save(rew_d, file = "../data/value.RData")

Now let's visualize the results.

    load("../data/fwd_sim_results.RData")
    ### Figure 2: trajectories of the % of biomass above B/B_MSY for each scenario
    font = 22
    qqq %>% filter(Var1 <= 2050) %>% ggplot(aes(Var1, traj)) +
      stat_summary(aes(color = Var3), geom = "line", fun.y = mean, lwd = 1) + 
      stat_summary(aes(fill = Var3), geom = "ribbon", fun.data = mean_sdl, fun.args = list(mult = 2), alpha = .2) +
      theme_bw() + ylab(expression(paste("% of stock above 80% ", B, "/", B[MSY]))) + xlab("Year")  + 
      facet_wrap(~Var4, nrow = 2, strip.position = "top") +
      scale_color_manual(values = c("goldenrod1","#00BFC4","black","#F8766D","blueviolet")) + 
      scale_fill_manual(values = c("goldenrod1","#00BFC4","black","#F8766D","blueviolet")) + 
      theme(legend.title=element_blank(), axis.title.x = element_text(size = font), 
            axis.title.y = element_text(size = font), legend.text = element_text(size = font) ,
            legend.position = "top",axis.text = element_text(size = font), strip.text = element_text(size = font),
            plot.margin=unit(c(1,1,0.5,0.5),"cm")) + coord_fixed(ylim = c(0,100), ratio = 0.56)

    ## Warning: Removed 302000 rows containing non-finite values (stat_summary).

    ## Warning: Removed 302000 rows containing non-finite values (stat_summary).

![](RAM_plots_files/figure-markdown_strict/unnamed-chunk-4-1.png)

    #### Figure 1, plotting by region or species


    N_Atl = c(1  , 2  , 6 , 10  ,11 , 21 , 29 ,30,  31,  43 , 44 , 47 , 74, 100, 101 ,107)
    Indian = c(15 ,88, 96)
    N_Pac = c(3   ,9 , 14 , 22,  23 , 24  ,32  ,33 , 34 , 37  ,54 , 55 , 64 , 72  ,76 , 85 , 94 , 95, 106)
    S_Pac = c(4  ,48 ,108)
    Bering = c(5   ,7 , 12  ,39 , 56 , 58  ,69  ,71,  79,  84,  92,  97 ,102, 103 ,109)
    Medit = c(9 ,156 ,185 ,191 ,194)
    Africa = c(25 ,35, 36, 83)
    Japan =c(20,  27,  28  ,51,  52 , 53 , 78,  82, 105)
    Alaska = c(8  ,38 , 57  ,70 , 75,  80, 104)
    Ross = c(86 ,98)
    NZ = c(13, 16, 17, 18, 19, 42 ,49 ,50, 59 ,60, 61, 62, 63 ,65 ,66 ,67, 68, 77, 87, 89, 90, 91)
    Mexico = c(40)
    North = c(73)
    S_Atl = c(81)
    Biscay = c(45)


    Flounder = c(7  , 8  , 9 ,107)
    Tuna = c(2 , 3 , 4, 10, 11, 14 ,15, 64, 86, 88)
    Mackerel = c(12 ,26 ,27, 28, 35 ,48, 53)
    Rock = c(33 , 34 , 38 , 56 , 57,  79 , 80, 106)
    Cod = c(29 ,30 ,31, 69 ,70)
    Sole = c(37  ,39 , 58 , 76, 109)
    Hake = c(25 ,36 , 46 , 72 ,93 ,100)
    Oreo = c(17, 18, 89 ,90 ,91)
    Anchovy = c(51, 52)
    Plaice = c(5 , 6 ,73)
    Herring = c(47 ,82)
    Pollock = c(74 ,102 ,103 ,104 ,105)
    Haddock = c(43 ,44 ,45)
    Ling = c(59, 60 ,61, 62)
    Whiting = c(87)
    Crab = c(92, 97)
    Other = c(66  ,67 , 83,  20 , 78 , 95 , 96 , 22  ,23 , 24  ,77 , 94 ,101  ,54 ,108 , 32,  84 , 85 ,
              1 ,  6 ,  9 , 38 , 10  ,14  ,16,  53 , 17  ,26,  60,  34 , 40,  44)


    # You can change the areas here by their name
    cc = N_Atl
    ind_f = cc

    pomdp = ff_pomdp_m[ind_f,50]
    mdp = ff_mdp_m[ind_f,50]
    msy = ff_msy_m[ind_f,50]
    bau = ff_bau_m[ind_f,50]
    q <- expand.grid(c("POMDP","MDP","MSY","BAU"))
    q <- data.frame(q,c(mean(pomdp),mean(mdp),mean(msy), mean(bau)))

    names(q)[1] = "Model"
    names(q)[2] = "Perc"

    q %>% ggplot() + 
      geom_hline(yintercept = 0, alpha = 0.9, size = 1.5) +
      geom_hline(yintercept = 25, alpha = 0.9, size = 1.5) +
      geom_hline(yintercept = 50, alpha = 0.9, size = 1.5) +
      geom_hline(yintercept = 75, alpha = 0.9, size = 1.5) +
      geom_hline(yintercept = 100, alpha = 0.9, size = 1.5) +
      geom_col(aes(Model,Perc*100, col = Model, fill = Model), width = .85) +
      scale_color_manual(values = c("#00BFC4","goldenrod1","#F8766D","blueviolet")) + 
      scale_fill_manual(values = c("#00BFC4","goldenrod1","#F8766D","blueviolet")) +
      theme_bw(base_size = 24) +coord_fixed(ylim = c(0,100), ratio = 0.032) +
      theme(legend.title = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(),
            axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
            legend.text = element_blank(),legend.position = "none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            panel.background = element_blank(),panel.border = element_rect(size = 1.5),plot.background = element_rect(fill = "transparent",colour = NA)) 

![](RAM_plots_files/figure-markdown_strict/unnamed-chunk-4-2.png)

    ### Plotting by species

    fish = c("Flounder","Tuna","Mackerel","Rockfish","Cod","Hake","Oreo","Anchovy","Plaice",
             "Herring","Pollock","Crab","Haddock","Whiting","Other")


    model = c("POMDP","MDP","MSY","BAU")

    fish_data = expand.grid(fish,model)

    min_fish = array(0, dim = dim(fish_data)[1])
    max_fish = array(0, dim = dim(fish_data)[1])

    fish2 = c("Flounder","Tuna","Mackerel","Rock","Cod","Hake","Oreo","Anchovy","Plaice",
             "Herring","Pollock","Crab","Haddock","Whiting","Other")

    for(i in 1:length(fish2)){
      
      xx=  eval(parse(text = fish2[i]))
      history = ind_biomass[xx,28]
      if(sum(!history) != 0){
        ind_f = xx[!history]
      } else{
        ind_f = xx[history]
      }
      

      min_fish[fish_data$Var1 == fish[i] & fish_data$Var2 == "MDP"] = mean(ff_mdp_m[ind_f,50]) * 100

      min_fish[fish_data$Var1 == fish[i] & fish_data$Var2 == "POMDP"] = mean(ff_pomdp_m[ind_f,50]) * 100

      min_fish[fish_data$Var1 == fish[i] & fish_data$Var2 == "MSY"] = mean(ff_msy_m[ind_f,50]) * 100
      
      min_fish[fish_data$Var1 == fish[i] & fish_data$Var2 == "BAU"] = mean(ff_bau_m[ind_f,50]) * 100
    }
    fish_data= data.frame(fish_data,min_fish)


    fish_data %>% ggplot() + 
      geom_hline(yintercept = 0, alpha = 0.9, size = 1.5) +
      geom_hline(yintercept = 25, alpha = 0.9, size = 1.5) +
      geom_hline(yintercept = 50, alpha = 0.9, size = 1.5) +
      geom_hline(yintercept = 75, alpha = 0.9, size = 1.5) +
      geom_hline(yintercept = 100, alpha = 0.9, size = 1.5) +
      geom_col(aes(Var1,min_fish, col = Var2, fill = Var2), width = 0.5, position=position_dodge(0.75)) +
      scale_color_manual(values = c("#00BFC4","goldenrod1","#F8766D","blueviolet")) + 
      scale_fill_manual(values = c("#00BFC4","goldenrod1","#F8766D","blueviolet")) +
      theme_bw(base_size = 30) + coord_fixed(ylim = c(0,100), ratio = .07) + 
      ylab(expression(paste("% of stock recovered above 80% ", B, "/", B[MSY]))) +
      theme(legend.title = element_blank(), axis.title.x = element_blank(),
            legend.position = "top",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            panel.background = element_blank(),panel.border = element_rect(size = 1),plot.background = element_rect(fill = "transparent",colour = NA)) +
      geom_vline(xintercept=seq(1.5, length(unique(fish_data$Var1))-0.5, 1), 
                 lwd=.5, colour="grey", linetype = "dashed") + scale_y_discrete(expand = c(0, 0)) 

![](RAM_plots_files/figure-markdown_strict/unnamed-chunk-4-3.png)

    ### Figure 3, the values

    load("../data/value.RData")
    rew_d %>% ggplot(aes(Var1,mean, fill = Var1)) + 
      geom_bar(position=position_dodge(), stat="identity",width = .85) +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.9)) +
      facet_wrap(~Var2, nrow = 2, strip.position = "top") +
      scale_color_manual(values = c("#00BFC4","goldenrod1","#F8766D","bluevioloet")) + 
      scale_fill_manual(values = c("#00BFC4","goldenrod1","#F8766D","blueviolet")) +
      ylab("Long term economic return relative to ORM") + 
      theme_bw(base_size = 24) +coord_fixed(ylim = c(0,1), ratio = 3) +
      theme(legend.title = element_blank(), axis.title.x = element_blank(),
            legend.position = "top",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            panel.background = element_blank(),panel.border = element_rect(size = 1),
            plot.background = element_rect(fill = "transparent",colour = NA),  axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      geom_vline(xintercept=seq(1.5, length(unique(rew_d$Var1))-0.5, 1), 
                 lwd=.5, colour="grey", linetype = "dashed")  + scale_y_continuous(expand = c(0, 0)) 

![](RAM_plots_files/figure-markdown_strict/unnamed-chunk-4-4.png)