# SMART3armNI: R codes for SMART 3 arm NI trial.
R codes for SMART 3 arm NI trial.

The codes are provided the R code to (1) generate the data; (2) plot power curve; (3) calculate sample size and power. 

Background: 

Consider a SMART design where individuals are randomized at the first stage to one of three initial intervention options, namely, A, P (Placebo), B. Based on the individual's progress during the first intervention stage, the participant can be classified as either a responder (R = 1) or a non-responder (R = 0) if they are in the A or B. At the second stage, responders are allowed to continue with the same initial intervention option, whereas non-responders are re-randomized to one of two second-stage options: C or D. We assume that the participants with first stage placebo group continue the same. Five adaptive interventions (AIs) are embedded in this SMART:
d1 = {A, A^R C^{1-R}}, d2 = {A, A^R D^{1-R}}, d3 = {B, B^R C^{1-R}}, d4 = {B, B^R D^{1-R}}, d5 = {P, P}.

For clarity in comparing embedded AIs, we classify any pair of embedded AIs under comparison to be either (1) distinct path (DP): those starting with different initial interventions (e.g., {d1, d5, d3} or {d1, d5, d4} or {d2, d5, d3} or {d2, d5, d4}) or (2) shared path (SP): those starting with the same initial intervention except the placebo (e.g., {d1, d5, d2} or {d3, d5, d4}).

We assume a single continuous primary outcome Y observed at the end of the trial. Let T_1 and T_2 denote the intervention options at stages 1 and 2 and N denote the total number of individuals in the trial. Since the distribution of the primary outcome is indexed by the sequences of intervention options received, we write an individual's outcome as Y_{T_1T_2}. Then we assume E(Y_{T_1T_2}) = mu_{T_1T_2}, Var(Y_{T_1T_2}) = sigma^2, otherwise, we can use unequal variances: sigma_E^2, sigma_R^2, and sigma_P^2 corresponding to experimental, reference, and placebo arms, respectively.

Define the first-stage randomization probability in favor of intervention option T_1 as pi_{T_1} and the second-stage randomization probability for those who started with the first-stage option T_1, in favor of intervention option T_2 as pi_{T_1T_2}. Note that, pi_{P} = 1 - pi_{A} - pi_{B}, pi_{AD} = 1 - pi_{AC}, and pi_{BD} = 1 - pi_{BC}.
See other details in the paper. In the paper, A is same as "a", B is same as "ac", P is same as "p", C is same as "v", and D is same as "m".

Structure of SMART 3 arm trial:


                    Responder ------> A
          1. A -->  Nonresponder ---> C --> d1 (AC)
                    Nonresponder ---> D --> d2 (AD)
    R --> 2. P ---------------------> P --> d5 (PP)
                    Responder ------> B
          3. B -->  Nonresponder ---> C --> d3 (BC)
                    Nonresponder ---> D --> d4 (BD)


We give a brief description of the R files below:

1. data_smart3arm: 

It generates the two stage data for three arm SMART design. It can be used to carry out the sample size and power calculation for two stage three arm non-inferiority SMART trial. 

Input:

    N_val: sample size.

    alpha_val: level of significance.

    beta_val: type II error.

    sig2E: variance corresponding to experimental treatment sequence.

    sig2R: variance corresponding to reference treatment sequence.

    sig2P: variance corresponding to placebo treatment sequence.

    theta_val: non-inferiority margin, [0.5, 1).

    pi_A: probability corresponding to arm A at first stage.

    pi_B: probability corresponding to arm B at first stage.

    pi_P: probability corresponding to arm P at first stage.

    gamma_A: response rate corresponding to A.

    gamma_B: response rate corresponding to B.

    pi_AC: probability corresponding to arm A at 1st stage and C at 2nd stage.

    pi_BC: probability corresponding to arm B at 1st stage and C at 2nd stage.

    mu_LA: mean of the latent variable corresponding to arm A.

    sig2_L: variance of the latent variable.

    mu_LP: mean of the latent variable corresponding to arm A.

    zeta_0, zeta_1A, zeta_1B, zeta_1P: coefficients of same treatment sequences at both stage.

    xi_0, xi_1A, xi_1B, xi_2AC, xi_2AD, xi_2BC, xi_2BD: coefficients of different treatment sequences.

Output: 

data_smart3arm returns a list containing the following components:

    Data: data consisting of three columns: Treatment1, Treatment2, Sample(Y);

    mu_T1T2: means corresponding to seven treatment sequences: AA, AC, AD, BB, BC, BD, PP;

    N_T1T2: sample sizes corresponding to seven treatment sequences: AA, AC, AD, BB, BC, BD, PP;

    Effect_size: effect sizes corresponding to six possible hypotheses containing: {d1, d5, d3}, {d1, d5, d4} , {d2, d5, d3}, {d2, d5, d4}, {d1, d5, d2}, {d3, d5, d4};

    mu_d: regimen means corresponding to five adaptive interventions: d1, d2, d3, d4, d5.

Example:

    ##Specify the values

    N_val = 300

    alpha_val = 0.05 #level alpha

    beta_val = 0.2 #type II error

    #equal variance

    sig2 = 10^2

    sig2E = sig2R = sig2P = sig2

    #NI margin

    theta_val = 0.6

    #Response rates

    gamma_A = 0.4

    gamma_B = 0.3

    #Probabilitites

    pi_A = 1/3 #probability corresponding to arm A at first stage

    pi_B = 1/3 #probability corresponding to arm B at first stage

    pi_P = 1/3 #probability corresponding to arm P at first stage

    pi_AC = 1/2 #probability corresponding to arm A at 1st stage and C at 2nd stage

    pi_AD = 1 - pi_AC #probability corresponding to arm A at 1st stage and D at 2nd stage

    pi_BC = 1/2 #probability corresponding to arm B at 1st stage and C at 2nd stage

    pi_BD = 1 - pi_BC #probability corresponding to arm B at 1st stage and D at 2nd stage

    ##Data generation parameters

    mu_LA = 2 #mean of the latent variable corresponding to arm A

    sig2_L = 0.2^2 #variance of the latent variable

    mu_LP = 2 #mean of the latent variable corresponding to arm A

    #parameters to calculate regimen means

    zeta_0 = 0.1

    zeta_1A = 2

    zeta_1B = 1.5

    zeta_1P = 0.1

    xi_0 = 0.03

    xi_1A = 0.1

    xi_1B = 0.6

    xi_2AD = 0.1

    xi_2BC = 0.2

    xi_2BD = 0.5

    xi_2AC = 3.5

    #Generated 3 arm NI SMART data

    data_smart3arm(N_val, alpha_val, beta_val, sig2E, sig2R, sig2P, theta_val, pi_A, pi_B, pi_P, gamma_A, gamma_B, pi_AC, pi_BC, mu_LA, sig2_L, mu_LP, zeta_0, zeta_1A, zeta_1B, zeta_1P, xi_0, xi_1A, xi_1B, xi_2AC, xi_2AD, xi_2BC, xi_2BD)


2. power_curve:

It draws the power curve of two stage SMART design for three arm non-ineriority trial for predefined sample size.

Input:

    N_val: sample size.

    M: number of simulated datasets.

    alpha_val: level of significance.

    beta_val: type II error.

    path: path options are: "DP" (distinct path), "SP" (shared path).

    seq_EPR: treatment sequnces of experimental (E), placebo (P), and reference (R) arms. If path = "DP", the options are: "d1d5d3", "d1d5d4", "d2d5d3", "d2d5d4"; if path = "SP", the options are: "d1d5d2", "d3d5d4".

    sig2E: variance corresponding to experimental treatment sequence.

    sig2R: variance corresponding to reference treatment sequence.

    sig2P: variance corresponding to placebo treatment sequence.

    theta_val: non-inferiority margin, [0.5, 1).

    pi_A: probability corresponding to arm A at first stage.

    pi_B: probability corresponding to arm B at first stage.

    pi_P: probability corresponding to arm P at first stage.

    gamma_A: response rate corresponding to A.

    gamma_B: response rate corresponding to B.

    pi_AC: probability corresponding to arm A at 1st stage and C at 2nd stage.

    pi_BC: probability corresponding to arm B at 1st stage and C at 2nd stage.

    mu_LA: mean of the latent variable corresponding to arm A.

    sig2_L: variance of the latent variable.

    mu_LP: mean of the latent variable corresponding to arm A.

    zeta_0, zeta_1A, zeta_1B, zeta_1P: coefficients of same treatment sequences at both stage.

    xi_0, xi_1A, xi_1B, xi_2AC, xi_2AD, xi_2BC, xi_2BD: coefficients of different treatment sequences.

Output:

power_curve returns the power curve corresponding to specified values.

Example:

    ##Specify the values

    N_val = 300

    alpha_val = 0.05 #level alpha

    beta_val = 0.2 #type II error

    #Probabilitites

    pi_A = 1/3 #probability corresponding to arm A at first stage

    pi_B = 1/3 #probability corresponding to arm B at first stage

    pi_P = 1/3 #probability corresponding to arm P at first stage

    pi_AC = 1/2 #probability corresponding to arm A at 1st stage and C at 2nd stage

    pi_AD = 1 - pi_AC #probability corresponding to arm A at 1st stage and D at 2nd stage

    pi_BC = 1/2 #probability corresponding to arm B at 1st stage and C at 2nd stage

    pi_BD = 1 - pi_BC #probability corresponding to arm B at 1st stage and D at 2nd stage

    ##Data generation parameters

    mu_LA = 2 #mean of the latent variable corresponding to arm A

    sig2_L = 0.2^2 #variance of the latent variable

    mu_LP = 2 #mean of the latent variable corresponding to arm A

    #M datasets

    M = 1000

    #equal variance

    sig2 = 1^2

    sig2E = sig2

    sig2R = sig2P = sig2

    #NI cutoff

    theta_val = 0.8

    #response rates

    gamma_A = 0.3

    gamma_B = 0.4

    #parameters to calculate regimen means

    zeta_0 = 0.02

    zeta_1A = 0.9

    zeta_1B = 1.2

    zeta_1P = 0.1

    xi_0 = 0.03

    xi_1A = 0.25

    xi_1B = 0.5

    xi_2AC = seq(1.2, 6, length = 20)

    xi_2AD = 1.3

    xi_2BC = 1.3

    xi_2BD = 0.8

    #Power

    powercurve_dp = matrix(NA, length(xi_2AC), 2)

    for(k in 1:length(xi_2AC)){

    powercurve_dp[k, ] = power_curve(N_val, M, alpha_val, beta_val, path, seq_EPR, sig2E, sig2R, sig2P,
    theta_val, pi_A, pi_B, pi_P, gamma_A, gamma_B, pi_AC, pi_BC, mu_LA, sig2_L, mu_LP, zeta_0, zeta_1A,
    zeta_1B, zeta_1P, xi_0, xi_1A, xi_1B, xi_2AC[k], xi_2AD, xi_2BC, xi_2BD)

    }

    #Plot the power curve

    plot(powercurve_dp[,1], powercurve_dp[,2], type = c("l"), col = 1:4, lwd = rep(2, 4), ylab = "Power", xlab = expression(eta^{DP}), main = "Distinct path: d1, d3")

3. sample_power:

sample_power calculates the sample size and estimates the empirical power of two stage SMART design for three arm non-ineriority trial.

Input:

    M: number of simulated datasets.

    alpha_val: level of significance.

    beta_val: type II error.

    path: path options are: "DP" (distinct path), "SP" (shared path).

    seq_EPR: treatment sequnces of experimental (E), placebo (P), and reference (R) arms. If path = "DP", the options are: "d1d5d3", "d1d5d4", "d2d5d3", "d2d5d4"; if path = "SP", the options are: "d1d5d2", "d3d5d4".

    sig2E: variance corresponding to experimental treatment sequence.

    sig2R: variance corresponding to reference treatment sequence.

    sig2P: variance corresponding to placebo treatment sequence.

    theta_val: non-inferiority margin, [0.5, 1).

    pi_A: probability corresponding to arm A at first stage.

    pi_B: probability corresponding to arm B at first stage.

    pi_P: probability corresponding to arm P at first stage.

    gamma_A: response rate corresponding to A.

    gamma_B: response rate corresponding to B.

    pi_AC: probability corresponding to arm A at 1st stage and C at 2nd stage.

    pi_BC: probability corresponding to arm B at 1st stage and C at 2nd stage.

    mu_LA: mean of the latent variable corresponding to arm A.

    sig2_L: variance of the latent variable.

    mu_LP: mean of the latent variable corresponding to arm A.

    zeta_0, zeta_1A, zeta_1B, zeta_1P: coefficients of same treatment sequences at both stage.

    xi_0, xi_1A, xi_1B, xi_2AC, xi_2AD, xi_2BC, xi_2BD: coefficients of different treatment sequences.

Output:

sample_power returns a vector containing the following components:

    SES: standardized effect size;

    N: required sample size;

    power: empirical power.

Example:

    ##Specify the values

    M = 5000 #M datasets

    alpha_val = 0.05 #level alpha

    beta_val = 0.2 #type II error

    #Probabilitites

    pi_A = 1/3 #probability corresponding to arm A at first stage

    pi_B = 1/3 #probability corresponding to arm B at first stage

    pi_P = 1/3 #probability corresponding to arm P at first stage

    pi_AC = 1/2 #probability corresponding to arm A at 1st stage and C at 2nd stage

    pi_AD = 1 - pi_AC #probability corresponding to arm A at 1st stage and D at 2nd stage

    pi_BC = 1/2 #probability corresponding to arm B at 1st stage and C at 2nd stage

    pi_BD = 1 - pi_BC #probability corresponding to arm B at 1st stage and D at 2nd stage

    ##Data generation parameters

    mu_LA = 2 #mean of the latent variable corresponding to arm A

    sig2_L = 0.2^2 #variance of the latent variable

    mu_LP = 2 #mean of the latent variable corresponding to arm A

    #parameters to calculate regimen means

    zeta_0 = 0.1

    zeta_1A = 2

    zeta_1B = 1.5

    zeta_1P = 0.1

    xi_0 = 0.03

    xi_1A = 0.1

    xi_1B = 0.6

    xi_2AC = 3.5

    xi_2AD = 0.1 #For path = "SP", use 0.01

    xi_2BC = 0.2

    xi_2BD = 0.5

    #Path: Options are "DP", "SP"

    path = "DP"

    #Sequences corresponding to DP: "d1d5d3", "d1d5d4", "d2d5d3", "d2d5d4"

    #Sequences corresponding to SP: "d1d5d2", "d3d5d4"

    seq_EPR = "d1d5d3"

    #equal variance

    sig2 = 10^2

    sig2E = sig2

    sig2R = sig2P = sig2

    #NI margin

    theta_val = 0.6

    #Response rates

    gamma_A = 0.4

    gamma_B = 0.3

    #Power and sample size

    sample_power(M, alpha_val, beta_val, path, seq_EPR, sig2E, sig2R, sig2P, theta_val, pi_A, pi_B, pi_P, gamma_A, gamma_B, pi_AC, pi_BC, mu_LA, sig2_L, mu_LP, zeta_0, zeta_1A, zeta_1B, zeta_1P, xi_0, xi_1A, xi_1B, xi_2AC, xi_2AD, xi_2BC, xi_2BD)


