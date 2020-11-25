###################
# Data Generation #
###################

####################
#SMART: 3 arm trial#
####################

#Structure
#                 Responder ------> A
#       1. A -->  Nonresponder ---> C --> d1 (AC)
#                 Nonresponder ---> D --> d2 (AD)
# R --> 2. P ---------------------> P --> d5 (PP)
#                 Responder ------> B
#       3. B -->  Nonresponder ---> C --> d3 (BC)
#                 Nonresponder ---> D --> d4 (BD)

################
#Data generation
################
data_smart3arm = function(N_val, alpha_val, beta_val, sig2E, sig2R, sig2P, theta_val, pi_A, pi_B, pi_P,
                          gamma_A, gamma_B, pi_AC, pi_BC, mu_LA, sig2_L, mu_LP, zeta_0, zeta_1A, zeta_1B,
                          zeta_1P, xi_0, xi_1A, xi_1B, xi_2AC, xi_2AD, xi_2BC, xi_2BD){

  # In the paper, A is same as "a", B is same as "ac", P is same as "p", C is same as "v",
  # and D is same as "m"


  #Step 2: Fix N
  N = N_val

  #Step 3: Set latent variable L with mean mu_LA, variance sig2_L, observation N_A
  N_A = ceiling(N/3)
  sample_A = rnorm(N_A, mu_LA, sqrt(sig2_L))

  #Step 4: Define cutoff value nu corresponding  to the responder
  nu_L = mu_LA + sqrt(sig2_L)*qnorm(1 - gamma_A)

  #Step 5: Generate N_B with mean mu_LB, variance sig2_L
  mu_LB = nu_L - sqrt(sig2_L)*qnorm(1 - gamma_B)
  N_B = ceiling(N/3)

  sample_B = rnorm(N_B, mu_LB, sqrt(sig2_L))

  #Step 6: Generate N_P with mean mu_LP, variance sig2_L
  N_P = N - N_A - N_B
  sample_P = rnorm(N_P, mu_LP, sqrt(sig2_L))

  #Step 7: mean of nonresponders
  u_star = (nu_L - mean(sample_A))/sqrt(sig2_L)
  w_star = (nu_L - mean(sample_B))/sqrt(sig2_L)

  mu_LA_NR = mean(sample_A) - sqrt(sig2_L)*dnorm(u_star)/pnorm(u_star)
  mu_LB_NR = mean(sample_B) - sqrt(sig2_L)*dnorm(w_star)/pnorm(w_star)

  #Step 8: Regimen means
  mu_d = rep(NA, 5)
  mu_d[1] = (zeta_0 + zeta_1A*mean(sample_A))*gamma_A + (xi_0 + xi_1A*mean(sample_A) + xi_2AC*mu_LA_NR)*(1 - gamma_A)
  mu_d[2] = (zeta_0 + zeta_1A*mean(sample_A))*gamma_A + (xi_0 + xi_1A*mean(sample_A) + xi_2AD*mu_LA_NR)*(1 - gamma_A)
  mu_d[3] = (zeta_0 + zeta_1B*mean(sample_B))*gamma_B + (xi_0 + xi_1B*mean(sample_B) + xi_2BC*mu_LB_NR)*(1 - gamma_B)
  mu_d[4] = (zeta_0 + zeta_1B*mean(sample_B))*gamma_B + (xi_0 + xi_1B*mean(sample_B) + xi_2BD*mu_LB_NR)*(1 - gamma_B)
  mu_d[5] = (zeta_0 + zeta_1P*mean(sample_P))

  #Step 9: Define effect size
  ES = rep(NA, 6)
  ES[1] = delta_d3d5d1 = mu_d[1] - mu_d[5] - theta_val*(mu_d[3] - mu_d[5])
  ES[2] = delta_d4d5d1 = mu_d[1] - mu_d[5] - theta_val*(mu_d[4] - mu_d[5])
  ES[3] = delta_d3d5d2 = mu_d[2] - mu_d[5] - theta_val*(mu_d[3] - mu_d[5])
  ES[4] = delta_d4d5d2 = mu_d[2] - mu_d[5] - theta_val*(mu_d[4] - mu_d[5])
  ES[5] = delta_d2d5d1 = mu_d[1] - mu_d[5] - theta_val*(mu_d[2] - mu_d[5])
  ES[6] = delta_d4d5d3 = mu_d[3] - mu_d[5] - theta_val*(mu_d[4] - mu_d[5])

  #Step 10: Generate N_T1T2 from Normal(mu_T1T2, sig2)
  N_T1T2 = rep(NA, 7)
  mu_T1T2 = rep(NA, 7)
  sig_T1T2 = rep(NA, 7)

  #T1 = A, T2 = A
  N_T1T2[1] = ceiling(N_A*gamma_A)
  mu_T1T2[1] = zeta_0 + zeta_1A*mean(sample_A)
  sig_T1T2[1] = sig2E

  #T1 = A, T2 = C
  N_T1T2[2] = ceiling(N_A*(1 - gamma_A)*0.5)
  mu_T1T2[2] = xi_0 + xi_1A*mean(sample_A) + xi_2AC*mu_LA_NR
  sig_T1T2[2] = sig2E

  #T1 = A, T2 = D
  N_T1T2[3] = ceiling(N_A*(1 - gamma_A)*0.5)
  mu_T1T2[3] = xi_0 + xi_1A*mean(sample_A) + xi_2AD*mu_LA_NR
  sig_T1T2[3] = sig2R

  #T1 = B, T2 = B
  N_T1T2[4] = ceiling(N_B*gamma_B)
  mu_T1T2[4] = zeta_0 + zeta_1B*mean(sample_B)
  sig_T1T2[4] = sig2R

  #T1 = B, T2 = C
  N_T1T2[5] = ceiling(N_B*(1 - gamma_B)*0.5)
  mu_T1T2[5] = xi_0 + xi_1B*mean(sample_B) + xi_2BC*mu_LB_NR
  sig_T1T2[5] = sig2R

  #T1 = B, T2 = D
  N_T1T2[6] = ceiling(N_B*(1 - gamma_B)*0.5)
  mu_T1T2[6] = xi_0 + xi_1B*mean(sample_B) + xi_2BD*mu_LB_NR
  sig_T1T2[6] = sig2R

  #T1 = P, T2 = P
  N_T1T2[7] = N - sum( N_T1T2[1:6])
  mu_T1T2[7] = zeta_0 + zeta_1P*mean(sample_P)
  sig_T1T2[7] = sig2P


  trt1 = vector("list", 7)
  trt2 = vector("list", 7)
  smaple_T1T2 = vector("list", 7)
  for(i in 1: 7){
    if(i==1){
      trt1[[i]] = rep("A", N_T1T2[i])
      trt2[[i]] = rep("A", N_T1T2[i])
    }else if(i == 2){
      trt1[[i]] = rep("A", N_T1T2[i])
      trt2[[i]] = rep("C", N_T1T2[i])
    }else if(i == 3){
      trt1[[i]] = rep("A", N_T1T2[i])
      trt2[[i]] = rep("D", N_T1T2[i])
    }else if(i == 4){
      trt1[[i]] = rep("B", N_T1T2[i])
      trt2[[i]] = rep("B", N_T1T2[i])
    }else if(i == 5){
      trt1[[i]] = rep("B", N_T1T2[i])
      trt2[[i]] = rep("C", N_T1T2[i])
    }else if(i == 6){
      trt1[[i]] = rep("B", N_T1T2[i])
      trt2[[i]] = rep("D", N_T1T2[i])
    }else if(i == 7){
      trt1[[i]] = rep("P", N_T1T2[i])
      trt2[[i]] = rep("P", N_T1T2[i])
    }

    smaple_T1T2[[i]] = rnorm(N_T1T2[i], mu_T1T2[i], sqrt(sig_T1T2[i]))
  }

  data_all = cbind(as.factor(unlist(trt1)), as.factor(unlist(trt2)), as.numeric(unlist(smaple_T1T2)))
  colnames(data_all) = c("Treatment1", "Treatment2",  "Sample")

  list_all = list(data_all, mu_T1T2, N_T1T2, ES, mu_d)
  names(list_all) = c("Data", "mu_T1T2", "N_T1T2", "Effect_size", "mu_d")
  return(list_all)
}

