#############
#Power curves
#############
power_curve = function(N_val, M, alpha_val, beta_val, path, seq_EPR, sig2E, sig2R, sig2P, theta_val, pi_A,
                       pi_B, pi_P, gamma_A, gamma_B, pi_AC, pi_BC, mu_LA, sig2_L, mu_LP, zeta_0, zeta_1A,
                       zeta_1B, zeta_1P, xi_0, xi_1A, xi_1B, xi_2AC, xi_2AD, xi_2BC, xi_2BD){

  # In the paper, A is same as "a", B is same as "ac", P is same as "p", C is same as "v",
  # and D is same as "m"

  pval13 = rep(NA, M)
  pval13_vec = rep(NA, M)
  pval14 = rep(NA, M)
  pval14_vec = rep(NA, M)
  pval23 = rep(NA, M)
  pval23_vec = rep(NA, M)
  pval24 = rep(NA, M)
  pval24_vec = rep(NA, M)
  pval12 = rep(NA, M)
  pval12_vec = rep(NA, M)
  pval34 = rep(NA, M)
  pval34_vec = rep(NA, M)


  Td1d5d3  = rep(NA, M)
  Td1d5d4  = rep(NA, M)
  Td2d5d3  = rep(NA, M)
  Td2d5d4  = rep(NA, M)
  Td1d5d2  = rep(NA, M)
  Td3d5d4  = rep(NA, M)

  deltadp13 = rep(NA, M)
  deltadp14 = rep(NA, M)
  deltadp23 = rep(NA, M)
  deltadp24 = rep(NA, M)
  deltasp12 = rep(NA, M)
  deltasp34 = rep(NA, M)

  data_sim = vector("list", M)
  for(m in 1:M){
    #Generate the data for fixed N
    data_sim[[m]] = data_smart3arm(N_val, alpha_val, beta_val,  sig2E, sig2R, sig2P, theta_val, pi_A, pi_B,
                                   pi_P, gamma_A, gamma_B, pi_AC, pi_BC, mu_LA, sig2_L, mu_LP, zeta_0,
                                   zeta_1A, zeta_1B, zeta_1P, xi_0, xi_1A, xi_1B, xi_2AC, xi_2AD,
                                   xi_2BC, xi_2BD)


    #Calculate test statistic
    ##Mean corresponding to the generated data
    dbar = rep(NA, 7)
    dbarindex = vector("list", 7)
    for(i in 1: 7){
      if(i==1){
        dbarindex[[1]] = which((data_sim[[m]][[1]][,1] == 1) & (data_sim[[m]][[1]][,2] == 1))
        dbar[1] =  sum(data_sim[[m]][[1]][dbarindex[[1]],3]/(pi_A))
      }else if(i == 2){
        dbarindex[[2]] = which((data_sim[[m]][[1]][,1] == 1) & (data_sim[[m]][[1]][,2] == 3))
        dbar[2] =  sum(data_sim[[m]][[1]][dbarindex[[2]],3]/(pi_A*pi_AC))
      }else if(i == 3){
        dbarindex[[3]] = which((data_sim[[m]][[1]][,1] == 1) & (data_sim[[m]][[1]][,2] == 4))
        dbar[3] =  sum(data_sim[[m]][[1]][dbarindex[[3]],3]/(pi_A*pi_AD))
      }else if(i == 4){
        dbarindex[[4]] = which((data_sim[[m]][[1]][,1] == 2) & (data_sim[[m]][[1]][,2] == 2))
        dbar[4] =  sum(data_sim[[m]][[1]][dbarindex[[4]],3]/(pi_B))
      }else if(i == 5){
        dbarindex[[5]] = which((data_sim[[m]][[1]][,1] == 2) & (data_sim[[m]][[1]][,2] == 3))
        dbar[5] =  sum(data_sim[[m]][[1]][dbarindex[[5]],3]/(pi_B*pi_BC))
      }else if(i == 6){
        dbarindex[[6]] = which((data_sim[[m]][[1]][,1] == 2) & (data_sim[[m]][[1]][,2] == 4))
        dbar[6] =  sum(data_sim[[m]][[1]][dbarindex[[6]],3]/(pi_B*pi_BD))
      }else if(i == 7){
        dbarindex[[7]] = which((data_sim[[m]][[1]][,1] == 3) & (data_sim[[m]][[1]][,2] == 5))
        dbar[7] =  sum(data_sim[[m]][[1]][dbarindex[[7]],3]/(pi_P))
      }
    }

    dbar1AC = sum(dbar[1], dbar[2])/N_val
    dbar2AD = sum(dbar[1], dbar[3])/N_val
    dbar3BC = sum(dbar[4], dbar[5])/N_val
    dbar4BD = sum(dbar[4], dbar[6])/N_val

    ##Test statistics: Distinct path
    Td1d5d3[m] = dbar1AC - (dbar[7]/N_val) - theta_val*(dbar3BC - (dbar[7]/N_val))
    Td1d5d4[m] = dbar1AC - (dbar[7]/N_val) - theta_val*(dbar4BD - (dbar[7]/N_val))
    Td2d5d3[m] = dbar2AD - (dbar[7]/N_val) - theta_val*(dbar3BC - (dbar[7]/N_val))
    Td2d5d4[m] = dbar2AD - (dbar[7]/N_val) - theta_val*(dbar4BD - (dbar[7]/N_val))

    ##Test statistics: Shared path
    Td1d5d2[m] = dbar1AC - (dbar[7]/N_val) - theta_val*(dbar2AD - (dbar[7]/N_val))
    Td3d5d4[m] = dbar3BC - (dbar[7]/N_val) - theta_val*(dbar4BD - (dbar[7]/N_val))


    ##Mean corresponding to the model
    mu_d1 = data_sim[[m]][[5]][1]
    mu_d2 = data_sim[[m]][[5]][2]
    mu_d3 = data_sim[[m]][[5]][3]
    mu_d4 = data_sim[[m]][[5]][4]
    mu_d5 = data_sim[[m]][[5]][5]

    ##Variance corresponding to the model
    mu_AA = data_sim[[m]][[2]][1]
    mu_AC = data_sim[[m]][[2]][2]
    mu_AD = data_sim[[m]][[2]][3]
    mu_BB = data_sim[[m]][[2]][4]
    mu_BC = data_sim[[m]][[2]][5]
    mu_BD = data_sim[[m]][[2]][6]
    mu_PP = data_sim[[m]][[2]][7]

    sig_d1 = (sig2E*(gamma_A*pi_AC + (1 - gamma_A))/(pi_A*pi_AC)) +
      ((mu_AA^2)*gamma_A*(1 - gamma_A*pi_A)/pi_A) +
      ((mu_AC^2)*(1 - gamma_A)*(1 - (1 - gamma_A)*pi_A*pi_AC)/(pi_A*pi_AC)) -
      (2*gamma_A*(1 - gamma_A)*mu_AA*mu_AC)

    sig_d2 = (sig2R*(gamma_A*pi_AD + (1 - gamma_A))/(pi_A*pi_AD)) +
      ((mu_AA^2)*gamma_A*(1 - gamma_A*pi_A)/pi_A) +
      ((mu_AD^2)*(1 - gamma_A)*(1 - (1 - gamma_A)*pi_A*pi_AD)/(pi_A*pi_AD)) -
      (2*gamma_A*(1 - gamma_A)*mu_AA*mu_AD)

    sig_d3 = (sig2R*(gamma_B*pi_BC + (1 - gamma_B))/(pi_B*pi_BC)) +
      ((mu_BB^2)*gamma_B*(1 - gamma_B*pi_B)/pi_B) +
      ((mu_BC^2)*(1 - gamma_B)*(1 - (1 - gamma_B)*pi_B*pi_BC)/(pi_B*pi_BC)) -
      (2*gamma_B*(1 - gamma_B)*mu_BB*mu_BC)

    sig_d4 = (sig2R*(gamma_B*pi_BD + (1 - gamma_B))/(pi_B*pi_BD)) +
      ((mu_BB^2)*gamma_B*(1 - gamma_B*pi_B)/pi_B) +
      ((mu_BD^2)*(1 - gamma_B)*(1 - (1 - gamma_B)*pi_B*pi_BD)/(pi_B*pi_BD)) -
      (2*gamma_B*(1 - gamma_B)*mu_BB*mu_BD)

    sig_d5 = (sig2P/pi_P) + ((mu_PP^2)*(1 - pi_P)/pi_P)

    sig_d1d2 =  (gamma_A*(sig2E + mu_AA^2)/pi_A) -
      ((gamma_A*mu_AA + (1 - gamma_A)*mu_AC)*(gamma_A*mu_AA + (1 - gamma_A)*mu_AD))

    sig_d3d4 =  (gamma_B*(sig2E + mu_BB^2)/pi_B) -
      ((gamma_B*mu_BB + (1 - gamma_B)*mu_BC)*(gamma_B*mu_BB + (1 - gamma_B)*mu_BD))

    #DP: experimental = 1 or 2, reference = 3 or 4
    var_d1d5d3 = (sig_d1 + (theta_val^2)*sig_d3 + ((1 - theta_val)^2)*sig_d5)/N_val
    var_d1d5d4 = (sig_d1 + (theta_val^2)*sig_d4 + ((1 - theta_val)^2)*sig_d5)/N_val
    var_d2d5d3 = (sig_d2 + (theta_val^2)*sig_d3 + ((1 - theta_val)^2)*sig_d5)/N_val
    var_d2d5d4 = (sig_d2 + (theta_val^2)*sig_d4 + ((1 - theta_val)^2)*sig_d5)/N_val

    #SP: experimental = 1 or 23, reference = 2 or 4
    var_d1d5d2 = (sig_d1 + (theta_val^2)*sig_d2 + ((1 - theta_val)^2)*sig_d5 - 2*theta_val*sig_d1d2)/N_val
    var_d3d5d4 = (sig_d3 + (theta_val^2)*sig_d4 + ((1 - theta_val)^2)*sig_d5 - 2*theta_val*sig_d3d4)/N_val

    #delta: effect size corresponding to DP
    deltadp13[m] = (mu_d1 - mu_d5 - theta_val*(mu_d3 - mu_d5))
    deltadp14[m] = (mu_d1 - mu_d5 - theta_val*(mu_d4 - mu_d5))
    deltadp23[m] = (mu_d2 - mu_d5 - theta_val*(mu_d3 - mu_d5))
    deltadp24[m] = (mu_d2 - mu_d5 - theta_val*(mu_d4 - mu_d5))

    #delta: effect size corresponding to SP
    deltasp12[m] = (mu_d1 - mu_d5 - theta_val*(mu_d2 - mu_d5))
    deltasp34[m] = (mu_d3 - mu_d5 - theta_val*(mu_d4 - mu_d5))

    #Standardized test statistic under H_0 corresponding to DP
    test_stat13 = (Td1d5d3[m])/sqrt(var_d1d5d3)
    test_stat14 = (Td1d5d4[m])/sqrt(var_d1d5d4)
    test_stat23 = (Td2d5d3[m])/sqrt(var_d2d5d3)
    test_stat24 = (Td2d5d4[m])/sqrt(var_d2d5d4)

    #Standardized test statistic under H_0 corresponding to SP
    test_stat12 = (Td1d5d2[m])/sqrt(var_d1d5d2)
    test_stat34 = (Td3d5d4[m])/sqrt(var_d3d5d4)

    #P-value
    pval13[m] =  1 - pnorm(test_stat13)
    pval14[m] =  1 - pnorm(test_stat14)
    pval23[m] =  1 - pnorm(test_stat23)
    pval24[m] =  1 - pnorm(test_stat24)

    pval12[m] =  1 - pnorm(test_stat12)
    pval34[m] =  1 - pnorm(test_stat34)


    #If p value < alpha => 1, else 0
    pval13_vec[m] = ifelse(pval13[m] < alpha_val, 1, 0)
    pval14_vec[m] = ifelse(pval14[m] < alpha_val, 1, 0)
    pval23_vec[m] = ifelse(pval23[m] < alpha_val, 1, 0)
    pval24_vec[m] = ifelse(pval24[m] < alpha_val, 1, 0)

    pval12_vec[m] = ifelse(pval12[m] < alpha_val, 1, 0)
    pval34_vec[m] = ifelse(pval34[m] < alpha_val, 1, 0)
  }

  #Sample size and Empirical power
  if(path == "DP" && seq_EPR == "d1d5d3"){
    power_val = mean(pval13_vec)
    stdES = mean(deltadp13)/sqrt((sig_d1 + (theta_val^2)*sig_d3 + ((1 - theta_val)^2)*sig_d5))
  }else if(path == "DP" && seq_EPR == "d1d5d4"){
    power_val = mean(pval14_vec)
    stdES = mean(deltadp14)/sqrt((sig_d1 + (theta_val^2)*sig_d4 + ((1 - theta_val)^2)*sig_d5))
  }else if(path == "DP" && seq_EPR == "d2d5d3"){
    power_val = mean(pval23_vec)
    stdES = mean(deltadp23)/sqrt((sig_d2 + (theta_val^2)*sig_d3 + ((1 - theta_val)^2)*sig_d5))
  }else if(path == "DP" && seq_EPR == "d2d5d4"){
    power_val = mean(pval24_vec)
    stdES = mean(deltadp24)/sqrt((sig_d2 + (theta_val^2)*sig_d4 + ((1 - theta_val)^2)*sig_d5))
  }else if(path == "SP" && seq_EPR == "d1d5d2"){
    power_val = mean(pval12_vec)
    stdES = mean(deltasp12)/sqrt((sig_d1 + (theta_val^2)*sig_d2 + ((1 - theta_val)^2)*sig_d5 - 2*theta_val*sig_d1d2))
  }else if(path == "SP" && seq_EPR == "d3d5d4"){
    power_val = mean(pval34_vec)
    stdES = mean(deltasp34)/sqrt((sig_d3 + (theta_val^2)*sig_d4 + ((1 - theta_val)^2)*sig_d5 - 2*theta_val*sig_d3d4))
  }


  list_pwrcurve = cbind(stdES, power_val)
  return(list_pwrcurve)

}
