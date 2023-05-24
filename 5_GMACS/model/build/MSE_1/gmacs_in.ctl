# TheHeader
## GMACS Version 2.01.M.01; ** MV **; Previous compilation on:  2023-03-23 10:54:35; Last compilation on:  2023-04-05 17:25:03
# ntheta
52
# Core parameters
## Initial: Initial value for the parameter (must lie between lower and upper)
## Lower & Upper: Range for the parameter
## Phase: Set equal to a negative number not to estimate
## Prior type:
## 0: Uniform   - parameters are the range of the uniform prior
## 1: Normal    - parameters are the mean and sd
## 2: Lognormal - parameters are the mean and sd of the log
## 3: Beta      - parameters are the two beta parameters [see dbeta]
## 4: Gamma     - parameters are the two gamma parameters [see dgamma]
# Initial_value Lower_bound Upper_bound Phase Prior_type Prior_1 Prior_2
0.28 0.15 0.7 -4 1 0.271 0.0154
16.5 -10 20 -2 0 -10 20
16.2343 -10 30 -1 0 10 20
9.38607 -10 30 -1 0 10 20
32.5 7.5 42.5 -4 0 32.5 2.25
1 0.1 10 -4 0 0.1 5
-0.9 -10 0.75 -4 0 -10 0.75
0.75 0.2 1 -2 3 3 2
0.01 0.0001 1 -3 3 1.01 1.01
-0.123787 -20 25 -1 0 10 20
-0.203795 -20 25 -1 0 10 20
0.015563 -20 25 -1 0 10 20
0.469611 -20 25 -1 0 10 20
0.935789 -20 25 -1 0 10 20
1.28392 -20 25 -1 0 10 20
1.32934 -20 25 -1 0 10 20
1.26751 -20 25 -1 0 10 20
1.31349 -20 25 -1 0 10 20
1.35939 -20 25 -1 0 10 20
1.26407 -20 25 -1 0 10 20
1.25457 -20 25 -1 0 10 20
1.33996 -20 25 -1 0 10 20
1.33339 -20 25 -1 0 10 20
0.914101 -20 25 -1 0 10 20
0.448411 -20 25 -1 0 10 20
-0.275944 -20 25 -1 0 10 20
-1.05449 -20 25 -1 0 10 20
-1.38357 -20 25 -1 0 10 20
-1.61397 -20 25 -1 0 10 20
-1.6406 -20 25 -1 0 10 20
0.525981 -20 25 -1 0 10 20
0.755816 -20 25 -1 0 10 20
1.2524 -20 25 -1 0 10 20
2.63572 -20 25 -1 0 10 20
2.88774 -20 25 -1 0 10 20
2.10586 -20 25 -1 0 10 20
1.88484 -20 25 -1 0 10 20
1.63281 -20 25 -1 0 10 20
1.43973 -20 25 -1 0 10 20
1.13356 -20 25 -1 0 10 20
0.646656 -20 25 -1 0 10 20
0.306411 -20 25 -1 0 10 20
0.0701145 -20 25 -1 0 10 20
-0.411507 -20 25 -1 0 10 20
-0.984487 -20 25 -1 0 10 20
-1.36407 -20 25 -1 0 10 20
-1.59499 -20 25 -1 0 10 20
-1.7698 -20 25 -1 0 10 20
-1.88836 -20 25 -1 0 10 20
-1.79584 -20 25 -1 0 10 20
-1.75286 -20 25 -1 0 10 20
-1.73808 -20 25 -1 0 10 20
# lw_type
1
# Using LW_RELATIONSHIP. Parameters are:
# lw_alfa
 0.001
# lw_beta
 3
# maturity
 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1
# legal
 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1
## Options for the growth matrix
## 1: Fixed growth transition matrix (requires molt probability)
## 2: Fixed size transition matrix (molt probability is ignored)
## 3: Growth increment is gamma distributed
## 4: Post-molt size is gamma distributed
## 5: Von Bert.: kappa varies among individuals
## 6: Von Bert.: Linf varies among individuals
## 7: Von Bert.: kappa and Linf varies among individuals
## 8: Growth increment is normally distributed
# bUseCustomGrowthMatrix
2
## Options for the growth increment model matrix
## 1: Linear
## 2: Individual
## 3: Individual (Same as 2)
# bUseGrowthIncrementModel
1
# bUseCustomMoltProbability
1
# nSizeClassRec
 5
# nSizeIncVaries
 1
# Start of the blocks in which molt increment changes (one row for each sex) - the first block starts in 1989
# Note: there is one less year than there are blocks
 # male
# nMoltVaries
 1
# Start of the blocks in which molt probability changes (one row for each sex) - the first block starts in 1989
# Note: there is one less year than there are blocks
# iYrsMoltChanges:
 # male
# BetaParRelative
0
# Growth parameters
# Initial_value Lower_bound Upper_bound Phase Prior_type Prior_1 Prior_2
 2.13986 -5 20 -3 1 2.049 1
 -0.218972 -1 0 -3 1 -0.2258 0.5
 0.25 0.001 5 -3 0 0 999
# Using custom growth matrix
# Using custom size transition matrix
# trans(CustomGrowthMatrix(h,i))
 0.0623394 0.36571 0.449704 0.115913 0.00626254 7.09226e-05 1.68358e-07 8.37717e-11 8.73728e-15 1.91016e-19 8.75343e-25 8.40818e-31 1.69293e-37 7.14485e-45 6.32064e-53 1.17204e-61 4.55556e-71 3.71154e-81 6.33843e-92 2.26895e-103 1.70248e-115 2.67765e-128
 0 0.0343043 0.283344 0.490563 0.178029 0.0135426 0.000215937 7.21716e-07 5.05617e-10 7.42491e-14 2.28547e-18 1.47461e-23 1.9943e-29 5.65353e-36 3.35942e-43 4.1843e-51 1.09244e-59 5.97841e-69 6.85788e-79 1.64895e-89 8.31078e-101 8.77991e-113
 0 0 0.0175288 0.203849 0.496913 0.253902 0.0271937 0.000610499 2.87287e-06 2.83375e-09 5.859e-13 2.53921e-17 2.30669e-22 4.39234e-28 1.75314e-34 1.46673e-41 2.57218e-49 9.45512e-58 7.28529e-67 1.17664e-76 3.98338e-87 2.82668e-98
 0 0 0 0.00831327 0.136119 0.467179 0.336094 0.050682 0.00160199 1.06141e-05 1.47408e-08 4.29114e-12 2.61842e-16 3.34905e-21 8.97881e-27 5.0458e-33 5.9437e-40 1.46757e-47 7.59547e-56 8.23997e-65 1.87375e-74 8.93125e-85
 0 0 0 0 0.00365859 0.0843437 0.407574 0.412834 0.0876516 0.00390084 3.63891e-05 7.11542e-08 2.91638e-11 2.50554e-15 4.51206e-20 1.70319e-25 1.34761e-31 2.23503e-38 7.7699e-46 5.66191e-54 8.64819e-63 2.76887e-72
 0 0 0 0 0 0.00149395 0.0484915 0.329921 0.470512 0.140652 0.00881324 0.000115755 3.18683e-07 1.83905e-10 2.22456e-14 5.64037e-19 2.99769e-24 3.33949e-30 7.79811e-37 3.81692e-44 3.91608e-52 8.4218e-61
 0 0 0 0 0 0 0.000566012 0.0258671 0.24779 0.497548 0.209412 0.0184749 0.000341647 1.32431e-06 1.076e-09 1.83254e-13 6.54198e-18 4.8953e-23 7.67829e-29 2.52444e-35 1.73972e-42 2.5131e-50
 0 0 0 0 0 0 0 0.000198969 0.0128026 0.172673 0.488166 0.289284 0.0359333 0.000935587 5.10606e-06 5.84121e-09 1.40066e-12 7.04013e-17 7.41723e-22 1.63802e-27 7.58244e-34 7.35724e-41
 0 0 0 0 0 0 0 0 6.4895e-05 0.00587915 0.111643 0.444393 0.370779 0.0648453 0.00237715 1.82662e-05 2.9421e-08 9.93299e-12 7.02939e-16 1.04273e-20 3.24218e-26 2.1131e-32
 0 0 0 0 0 0 0 0 0 1.9638e-05 0.0025049 0.0669731 0.37534 0.440925 0.108572 0.00560387 6.06279e-05 1.3749e-07 6.53559e-11 6.51199e-15 1.36006e-19 5.95409e-25
 0 0 0 0 0 0 0 0 0 0 5.51366e-06 0.000990208 0.0372758 0.294132 0.486489 0.168663 0.0122569 0.000186704 5.96135e-07 3.98978e-10 5.59718e-14 1.6459e-18
 0 0 0 0 0 0 0 0 0 0 0 1.43631e-06 0.000363184 0.0192494 0.213858 0.49802 0.243099 0.0248734 0.000533459 2.39818e-06 2.25984e-09 4.46364e-13
 0 0 0 0 0 0 0 0 0 0 0 0 3.47156e-07 0.000123593 0.00922309 0.144269 0.473028 0.325098 0.0468336 0.00141421 8.95132e-06 1.18761e-08
 0 0 0 0 0 0 0 0 0 0 0 0 0 7.78509e-08 3.90232e-05 0.00410012 0.0902997 0.41686 0.403375 0.0818168 0.00347849 3.09995e-05
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.61996e-08 1.14328e-05 0.00169129 0.0524444 0.340874 0.464413 0.132626 0.00793906
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.18135e-09 3.16121e-06 0.00065843 0.0287462 0.263067 0.504624 0.202901
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8.14189e-10 1.13909e-06 0.000334045 0.0205337 0.264573 0.714558
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.41736e-10 1.46108e-06 0.000603271 0.0522115 0.947184
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5.04841e-09 1.40013e-05 0.00813954 0.991846
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.12495e-07 0.00122025 0.998779
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000181854 0.999818
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
#  * Selectivity parameter controls
# 
## Selectivity parameter controls
# ## Selectivity (and retention) types
# ##  <0: Mirror selectivity
# ##   0: Nonparametric selectivity (one parameter per class)
# ##   1: Nonparametric selectivity (one parameter per class, constant from last specified class)
# ##   2: Logistic selectivity (inflection point and slope)
# ##   3: Logistic selectivity (50% and 95% selection)
# ##   4: Double normal selectivity (3 parameters)
# ##   5: Flat equal to zero (1 parameter; phase must be negative)
# ##   6: Flat equal to one (1 parameter; phase must be negative)
# ##   7: Flat-topped double normal selectivity (4 parameters)
# ##   8: Declining logistic selectivity with initial values (50% and 95% selection plus extra)
# ##   9: Cubic-spline (specified with knots and values at knots)
# ##  10: One parameter logistic selectivity (inflection point and slope)
# ## Extra (type 1): number of selectivity parameters to estimated
# #  Pot_Fishery Survey_1
#  # selectivity periods
# slx_nsel_period_in
 1 1
#  # sex specific selectivity (1=Yes, 0=No)
# slx_bsex_in
 0 0
#  # selectivity type (by sex)
# slx_type_in
 2 0
#  # selectivity within another gear
# slx_include_in
 0 0
#  # extra parameters for each pattern
# slx_extra_in
 0 0
#  # retention periods 
# ret_nret_period_in
 1 1
#  # sex specific retention (1=Yes, 0=No)
# ret_bsex_in
 0 0
#  # retention type (by sex)
# ret_type_in
 2 6
#  # retention flag
# slx_nret
 1 0
#  # extra parameters for each pattern
# ret_extra_in
 0 0
#  # Is maximum selectivity at size is forced to equal 1 or not
# slx_max_at_1_in
 1 0
#  
# Selectivity parameters
## Fleet: The index of the fleet  (positive for capture selectivity; negative for retention)
## Index: Parameter count (not used)
## Parameter_no: Parameter count within the current pattern (not used)
## Sex: Sex (not used)
## Initial: Initial value for the parameter (must lie between lower and upper)
## Lower & Upper: Range for the parameter
## Phase: Set equal to a negative number not to estimate
## Prior type:
## 0: Uniform   - parameters are the range of the uniform prior
## 1: Normal    - parameters are the mean and sd
## 2: Lognormal - parameters are the mean and sd of the log
## 3: Beta      - parameters are the two beta parameters [see dbeta]
## 4: Gamma     - parameters are the two gamma parameters [see dgamma]
## Start / End block: years to define the current block structure
## Env_Link: Do environmental impact ? (0/1)
## Env_Link_Var: Which environmental variable to consider for tihs parameter ? (column of Env data)
## Rand_Walk: Do a random walk? (0/1)
## Start_year_RW: Start year of the random walk
## End_year_RW: End year of the random walk
## Sigma_RW: Sigma used for the random walk
# Fleet Index Parameter_no Sex Initial Lower_bound Upper_bound Prior_type Prior_1 Prior_2 Phase Start_block End_block Env_Link Env_Link_Var Rand_Walk Start_year_RW End_year_RW Sigma_RW
 1 1 1 1 105.009 5 186 0 1 999 -4 1989 2019 0 0 0 0 0 0
 1 2 2 1 4.75444 0.01 20 0 1 999 -4 1989 2019 0 0 0 0 0 0
 2 3 1 1 0.0204142 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 4 2 1 0.020414 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 5 3 1 0.0204193 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 6 4 1 0.0205745 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 7 5 1 0.0210718 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 8 6 1 0.02435 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 9 7 1 0.0343391 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 10 8 1 0.0544918 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 11 9 1 0.0836207 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 12 10 1 0.128434 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 13 11 1 0.151843 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 14 12 1 0.183341 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 15 13 1 0.282837 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 16 14 1 0.400038 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 17 15 1 0.347033 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 18 16 1 0.228233 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 19 17 1 0.18586 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 20 18 1 0.173558 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 21 19 1 0.18528 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 22 20 1 0.226014 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 23 21 1 0.252906 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 2 24 22 1 0.272496 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0
 -1 25 1 1 0.862719 1 190 1 96 10 -4 1989 2019 0 0 0 0 0 0
 -1 26 2 1 0.592775 0.001 20 0 1 999 -4 1989 2019 0 0 0 0 0 0
 -2 27 1 0 0.999999 1 999 0 1 999 -3 2015 2019 0 0 0 0 0 0
#Number of asymptotic selectivity parameters
0
# Fleet Sex Year Initial lower_bound upper_bound phase

# Number of environmental parameters
# DevParPhase
-1
#Catchability
# Initial Lower_bound Upper_bound Phase Prior_type Prior_1 Prior_2 Index_lambda Index_lambda
 0.9999 0.01 1 -5 0 0.843136 0.03 0 1 1
# Index CV
# Initial Lower_bound Upper_bound Phase Prior_type Prior_1 Prior_2
 0.0001 1e-05 10 -4 0 1 100
# Additional variance controls
# 0 ignore; >0 use
 0
# Controls on F
# Initial_male_f Initial_female_F Penalty_SD (early phase) Penalty_SD (later Phase) Phase_mean_F_male Phase_mean_F_female Lower_bound_mean_F Upper_bound_mean_F Lower_bound_annual_male_F Upper_bound_annual_male_F Lower_bound_annual_female_F Upper_bound_annual_female_F
 1 0.0505 0.5 45.5 1 1 -12 4 -10 10 -10 10
 0 0 2 20 -1 -1 -12 4 -10 10 -10 10
# Options when fitting size-composition data
## Likelihood types: 
##  1:Multinomial with estimated/fixed sample size
##  2:Robust approximation to multinomial
##  3:logistic normal
##  4:multivariate-t
##  5:Dirichlet

#  Pot_Fishery Pot_Fishery Pot_Fishery Survey_1
#  male male male male
#  retained retained retained retained
#  all_shell all_shell all_shell all_shell
#  immature+mature immature+mature immature+mature immature+mature
 2 2 2 2 # Type of likelihood
 0 0 0 0 # Auto tail compression (pmin)
 1 1 1 1 # Initial value for effective sample size multiplier
 -4 -4 -4 -4 # Phz for estimating effective sample size (if appl.)
 1 1 1 2 # Composition aggregator codes
 1 1 1 2 # Set to 1 for catch-based predictions; 2 for survey or total catch predictions
 1 1 1 1 # Lambda for effective sample size
 1 1 1 1 # Lambda for overall likelihood
# Type of M specification
## 1: Time-invariant M
## 2: Default random walk M
## 3: Cubic spline with time M
## 4: Blocked changes in  M
## 5: Blocked changes in  M (type 2)
## 6: Blocked changes in  M (returns to default)
# m_type
0
# Mdev_phz_def
-4
# m_stdev
0
# m_nNode_sex
0 # male
# Start of the blocks in which M changes (one row for each sex) - the first block starts in 1989
# Note: there is one less year than there are blocks
 # male
# nSizeDevs
0
# Start of the size-class blocks in which M changes (one row for each sex) - the first block start at size-class 1
# Note: there is one less size-class than there are blocks (no input implies M is independent of size
# m_size_nodeyear

# Init_Mdev
0
# # Init_MDev==NO
1 # tag_emphasis
# # maturity specific natural mortality? (yes = 1; no = 0; only for use if nmature > 1)
# m_maturity
1
# # Initial Lower_bound Upper_bound Phase Prior_type Prior_1 Prior_2
 0.133531 -4 4 4 1 0 0.05

# Extra controls
1989 # First year of recruitment estimation
2019 # Last year of recruitment estimation
1 # Consider terminal molting (0 = off, 1 = on). If on, the calc_stock_recruitment_relationship() isn't called in the procedure
1 # Phase for recruitment estimation
2 # Phase for recruitment sex-ratio estimation
0.5 # Initial value for recruitment sex-ratio
-3 # Phase for initial recruitment estimation
1 # VERBOSE FLAG (0 = off, 1 = on, 2 = objective func; 3 diagnostics)
3 # Initial conditions (0 = Unfished, 1 = Steady-state fished, 2 = Free parameters, 3 = Free parameters (revised))
1 # Lambda (proportion of mature male biomass for SPR reference points)
0 # Stock-Recruit-Relationship (0 = none, 1 = Beverton-Holt)
10 # Maximum phase (stop the estimation after this phase)
-1 # Maximum number of function calls
1 # Calculate reference points (0=no)
0 # Use years specified to computed average sex ratio in the calculation of average recruitment for reference points (0 = off -i.e. Rec based on End year, 1 = on)
200 # Years to compute equilibria

# ## Emphasis Factors (Catch: number of catch dataframes) ##
# nCatchDF
3
# catch_emphasis
 1 1 1
# ## Emphasis Factors (Fdev Penalties; number of fleets) ##
# nfleet
2
# Penalty_fdevs
 1 1 0 0
 0 0 0 0
# ## Emphasis Factors (Priors/Penalties: 13 values) ##
10000	#--Log_fdevs
0	#--MeanF
0	#--Mdevs
0	#--Rec_devs
0	#--Initial_devs
0	#--Fst_dif_dev
0	#--Mean_sex_ratio
0	#--Molt_prob
0	#--free selectivity
0	#--Init_n_at_len
0	#--Fvecs
0	#--Fdovss
0	#--Random walk in selectivity
# eof_ctl
9999
