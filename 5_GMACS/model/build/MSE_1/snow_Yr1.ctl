# ============================================================ #
#                  GMACS main control file 
# 
#_*** 
#_  GMACS Version 2.01.M.01 
#_Last GMACS mofification made by:   ** MV ** 
#_Date of writing the control file: 2023-05-22 16:35:14 
#_*** 
# 
#_Stock of interest:  Snow crab 
#_Model name:  MSE_1 
#_Year of assessment:  2019 
# ============================================================ #

# -------------------------------------- #
##_Key parameter controls
# -------------------------------------- #
#_ntheta - Number of leading parameters (guestimated)
52 
#
#_Core parameters
# ************************************** #
#_For each parameter columns are:
#_Init_val: Initial value for the parameter (must lie between lower and upper bounds)
#_Lower_Bd & Upper_Bd: Range for the parameter
#_Phase: Set equal to a negative number not to estimate
#_Available prior types:
#_-> 0 = Uniform   - parameters are the range of the uniform prior
#_-> 1 = Normal    - parameters are the mean and sd
#_-> 2 = Lognormal - parameters are the mean and sd of the log
#_-> 3 = Beta      - parameters are the two beta parameters [see dbeta]
#_-> 4 = Gamma     - parameters are the two gamma parameters [see dgamma]
#_p1; p2: priors
# ************************************** #
# 
#_Init_val_| Lower_Bd_| Upper_Bd_| Phase_| Prior_type_| p1_| p2
0.28 0.15 0.7 -4 1 0.271 0.0154
16.5 -10 20 -2 0 -10 20
16.23426596 -10 30 -1 0 10 20
9.38606864 -10 30 -1 0 10 20
32.5 7.5 42.5 -4 0 32.5 2.25
1 0.1 10 -4 0 0.1 5
-0.9 -10 0.75 -4 0 -10 0.75
0.75 0.2 1 -2 3 3 2
0.01 1e-04 1 -3 3 1.01 1.01
-0.12378693 -20 25 -1 0 10 20
-0.20379527 -20 25 -1 0 10 20
0.01556296 -20 25 -1 0 10 20
0.46961068 -20 25 -1 0 10 20
0.93578927 -20 25 -1 0 10 20
1.283922 -20 25 -1 0 10 20
1.329338 -20 25 -1 0 10 20
1.26751204 -20 25 -1 0 10 20
1.31349274 -20 25 -1 0 10 20
1.35938613 -20 25 -1 0 10 20
1.26406998 -20 25 -1 0 10 20
1.25457069 -20 25 -1 0 10 20
1.33995928 -20 25 -1 0 10 20
1.3333891 -20 25 -1 0 10 20
0.91410052 -20 25 -1 0 10 20
0.44841086 -20 25 -1 0 10 20
-0.27594402 -20 25 -1 0 10 20
-1.05449365 -20 25 -1 0 10 20
-1.38357119 -20 25 -1 0 10 20
-1.61397215 -20 25 -1 0 10 20
-1.64060221 -20 25 -1 0 10 20
0.52598068 -20 25 -1 0 10 20
0.75581555 -20 25 -1 0 10 20
1.25239548 -20 25 -1 0 10 20
2.63572142 -20 25 -1 0 10 20
2.88774023 -20 25 -1 0 10 20
2.10585994 -20 25 -1 0 10 20
1.88483761 -20 25 -1 0 10 20
1.63280933 -20 25 -1 0 10 20
1.43972512 -20 25 -1 0 10 20
1.13355898 -20 25 -1 0 10 20
0.64665609 -20 25 -1 0 10 20
0.30641068 -20 25 -1 0 10 20
0.0701145 -20 25 -1 0 10 20
-0.41150691 -20 25 -1 0 10 20
-0.98448677 -20 25 -1 0 10 20
-1.36407207 -20 25 -1 0 10 20
-1.59498884 -20 25 -1 0 10 20
-1.76980477 -20 25 -1 0 10 20
-1.88835677 -20 25 -1 0 10 20
-1.79584274 -20 25 -1 0 10 20
-1.75286245 -20 25 -1 0 10 20
-1.73807951 -20 25 -1 0 10 20
# -------------------------------------- #

# -------------------------------------- #
##_Allometry
# -------------------------------------- #
#_Length-weight type/method
#_1 = Length-weight relationship (vector of sex specific parameters: w_l = a[s]*l^b[s])
#_2 = Input vector of mean weight-at-size by sex (dim=[1:nclass])
#_3 = Input matrix of mean weight-at-size by sex and year (dim=[nsex*Nyear; nclass])
1 
#_lw_alfa
0.001 
#_lw_beta
3 
# -------------------------------------- #

# -------------------------------------- #
##_Fecundity for MMB/MMA calculation
# -------------------------------------- #
#_Maturity definition: Proportion of mature at size by sex
0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 
#_Legal definition of the proportion of mature at size by sex
0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 
# -------------------------------------- #

# -------------------------------------- #
##_Growth parameter controls
# -------------------------------------- #
 
#_Two lines for each parameter are required if the model considers two sexes, one line is not
 
#_Growth transition matrix definition
# ************************************** #
#_1 = Fixed growth transition matrix (requires molt probability)
#_2 = Fixed size transition matrix (molt probability is ignored)
#_3 = Growth increment is gamma distributed
#_4 = Size after growth is gamma distributed
#_5 = kappa varies among individuals
#_6 = Linf varies among individuals
#_7 = kappa and Ling varies among individuals
#_8 = Growth increment is normally distributed
# ************************************** #
2 
 
#_Growth increment model matrix
# ************************************** #
#_0 = Pre-specified growth increment
#_1 = linear (alpha; beta parameters)
#_2 = Estimated by size-class
#_3 = Pre-specified by size-class (empirical approach)
# ************************************** #
1 
 
#_Molt probability function
# ************************************** #
#_0 = Pre-specified probability of molting
#_1 = Constant probability of molting (flat approach)
#_2 = Logistic function
#_3 = Free estimated parameters
# ************************************** #
#_If the custom growth model option = 1 then the molt probability function must be 1 
1 
#_Maximum of size-classes to which recruitment must occur (males then females)
5 
#_Number of blocks of growth matrix parameters (i.e., number of size-increment period)
1 
#_Year(s) with changes in the growth matrix
#_-> 1 line per sex - blank if no change (i.e., if the number of blocks of growth matrix parameters = 1)
 
#_Number of blocks of molt probability
1 
#_Year(s) with changes in molt probability
#_-> 1 line per sex - blank if no change (i.e., if the number of blocks of growth matrix parameters = 1)
#_Are the beta parameters relative to a base level?
0 

#_Growth increment model controls
# ************************************** #
#_For each parameter columns are:
#_Init_val: Initial value for the parameter (must lie between lower and upper bounds)
#_Lower_Bd & Upper_Bd: Range for the parameter
#_Phase: Set equal to a negative number not to estimate
#_Available prior types:
#_-> 0 = Uniform   - parameters are the range of the uniform prior
#_-> 1 = Normal    - parameters are the mean and sd
#_-> 2 = Lognormal - parameters are the mean and sd of the log
#_-> 3 = Beta      - parameters are the two beta parameters [see dbeta]
#_-> 4 = Gamma     - parameters are the two gamma parameters [see dgamma]
#_p1; p2: priors
# ************************************** #
# 
#_Init_val_| Lower_Bd_| Upper_Bd_| Phase_| Prior_type_| p1_| p2
2.13986432 -5 20 -3 1 2.049 1
-0.21897198 -1 0 -3 1 -0.2258 0.5
0.25 0.001 5 -3 0 0 999

#_Molt probability controls
# ************************************** #
#_For each parameter columns are:
#_Init_val: Initial value for the parameter (must lie between lower and upper bounds)
#_Lower_Bd & Upper_Bd: Range for the parameter
#_Phase: Set equal to a negative number not to estimate
#_Available prior types:
#_-> 0 = Uniform   - parameters are the range of the uniform prior
#_-> 1 = Normal    - parameters are the mean and sd
#_-> 2 = Lognormal - parameters are the mean and sd of the log
#_-> 3 = Beta      - parameters are the two beta parameters [see dbeta]
#_-> 4 = Gamma     - parameters are the two gamma parameters [see dgamma]
#_p1; p2: priors
# ************************************** #
 
#_Init_val_| Lower_Bd_| Upper_Bd_| Phase_| Prior_type_| p1_| p2

#_Custom growth-increment matrix or size-transition matrix (if any)
0.0623393849946927 0.365710289816406 0.449703967011724 0.115912728602741 0.00626253854605257 7.09225866663364e-05 1.68357936215739e-07 8.37717315107914e-11 8.73728018752682e-15 1.91016070064012e-19 8.75343317322169e-25 8.40817626366145e-31 1.6929341522096e-37 7.14485134227217e-45 6.32064210823777e-53 1.17204442939919e-61 4.55556053558838e-71 3.71154252674209e-81 6.33843340980389e-92 2.2689467062e-103 1.70247808512164e-115 2.67764894139527e-128
0 0.0343043022406814 0.283344238227879 0.490563298600434 0.178028934720081 0.013542567374292 0.00021593661475715 7.21716184448137e-07 5.05616772699934e-10 7.42491304863182e-14 2.28547325765748e-18 1.47460529105303e-23 1.99429914834805e-29 5.65352987867979e-36 3.35941727945352e-43 4.18430259474753e-51 1.0924387276187e-59 5.97841369077158e-69 6.8578773213902e-79 1.64895320818506e-89 8.31078072990912e-101 8.77991251003373e-113
0 0 0.0175287710592694 0.203848909982475 0.496912807321698 0.253902431011577 0.0271937061224907 0.000610498797332946 2.87287081617846e-06 2.83375470677977e-09 5.85899710476933e-13 2.53921330092994e-17 2.30669384777136e-22 4.39233647957256e-28 1.75313806443711e-34 1.46673450320219e-41 2.57218244247309e-49 9.45511689674691e-58 7.2852910666456e-67 1.17663505992339e-76 3.98337812422268e-87 2.82667680493393e-98
0 0 0 0.00831327280835128 0.136119418881629 0.467178676694792 0.33609403914467 0.050681969115052 0.00160199450405228 1.06141063731149e-05 1.47407887277596e-08 4.29114193585453e-12 2.61842361854341e-16 3.34905139052156e-21 8.97880595270539e-27 5.04580231768311e-33 5.94369820894068e-40 1.46756768464816e-47 7.59546562079528e-56 8.23996821669851e-65 1.87374955517911e-74 8.93125359618671e-85
0 0 0 0 0.00365858954837132 0.0843437433999934 0.407574427466054 0.412834376959026 0.0876515606199302 0.00390084167700247 3.63891462987967e-05 7.11541577109207e-08 2.91637635930311e-11 2.50554223676688e-15 4.51205889008341e-20 1.7031882896697e-25 1.34761376260374e-31 2.23502859967448e-38 7.76990312390519e-46 5.66191043762749e-54 8.64819071772565e-63 2.76886890304667e-72
0 0 0 0 0 0.00149394661083529 0.0484914744147125 0.329921465928262 0.470511918834521 0.140651885132891 0.0088132351393546 0.000115755072000465 3.18683496553784e-07 1.83905191696965e-10 2.2245556663582e-14 5.64036797404932e-19 2.99768828168514e-24 3.33949174646402e-30 7.79810629498853e-37 3.81691772772471e-44 3.91607787246157e-52 8.42179588913543e-61
0 0 0 0 0 0 0.000566012119331906 0.0258670660082196 0.247789833790785 0.497547638659583 0.209411595909415 0.0184748812528683 0.000341646878270866 1.32430533911869e-06 1.07600347813653e-09 1.83254265564618e-13 6.54198255436381e-18 4.89530260362877e-23 7.67829115532305e-29 2.52443667839869e-35 1.73971961033669e-42 2.51309501272885e-50
0 0 0 0 0 0 0 0.000198969092827061 0.0128026003883006 0.17267347615473 0.488166469827473 0.289284478537305 0.0359333072213628 0.000935586875285143 5.10606010839593e-06 5.84120700379812e-09 1.40066457101881e-12 7.04012810628556e-17 7.41723181091174e-22 1.63801552860768e-27 7.58244116278172e-34 7.35724054345248e-41
0 0 0 0 0 0 0 0 6.48949898683244e-05 0.00587915443006186 0.111643499032767 0.444392665524262 0.370779058861073 0.0648452822625177 0.00237714921944515 1.82662490882194e-05 2.94209829922966e-08 9.93298670837998e-12 7.02938615583718e-16 1.04272508400783e-20 3.24217992653583e-26 2.11309623423737e-32
0 0 0 0 0 0 0 0 0 1.96379513303639e-05 0.00250490354623158 0.066973147211913 0.375340372031194 0.440925000866734 0.108572300682758 0.00560387228740146 6.06278670405095e-05 1.374900348121e-07 6.53559341260866e-11 6.51198966807207e-15 1.36005744245042e-19 5.95409379527723e-25
0 0 0 0 0 0 0 0 0 0 5.5136639358042e-06 0.00099020758595106 0.0372758090596428 0.294132341200861 0.486489367753027 0.168662603000209 0.0122568569090744 0.000186704293613437 5.96134651738823e-07 3.98978322919472e-10 5.59717783963795e-14 1.64590128221878e-18
0 0 0 0 0 0 0 0 0 0 0 1.43631308772299e-06 0.000363183579071934 0.0192494427890518 0.213857749849215 0.498019898493954 0.243099054608909 0.0248733749381637 0.000533458985780199 2.39818247776083e-06 2.25984310403477e-09 4.4636406363647e-13
0 0 0 0 0 0 0 0 0 0 0 0 3.47156315401333e-07 0.000123592832556425 0.00922308844166875 0.144269435563002 0.473028343585017 0.325098426690317 0.046833589318179 0.00141421321419994 8.95132262941658e-06 1.18761152143509e-08
0 0 0 0 0 0 0 0 0 0 0 0 0 7.78508696148867e-08 3.90231939686674e-05 0.0041001243576879 0.0902996605490079 0.416859859989319 0.403374937334308 0.0818168237595407 0.00347849346872974 3.09994965684488e-05
0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.61995655077294e-08 1.14328324786365e-05 0.0016912943014381 0.0524444340619941 0.340874483485497 0.464413092412412 0.132626188354465 0.00793905835214946
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.18134818621336e-09 3.16120606293522e-06 0.000658429500457924 0.028746212229849 0.263067287223623 0.504624041018456 0.202900865640203
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 8.14188961316066e-10 1.13908974416881e-06 0.000334045438306339 0.0205337382306053 0.264572906914485 0.71455816951267
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7.41736400241902e-10 1.46107817355758e-06 0.000603270634875696 0.0522114543593864 0.947183813185828
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5.04841281358734e-09 1.40013480069419e-05 0.00813953719570028 0.99184645640788
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.12494575558112e-07 0.00122024939773035 0.998779438107694
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000181854380362521 0.999818145619637
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1

#_Custom molt probability matrix  (if any)
# 
# -------------------------------------- #

# -------------------------------------- #
##_Vulnerability parameter controls
# 
#_Vulnerability is the combination of selectivity and retention selectivity.
#_Gmacs requires that each gear has a vulnerability.
# 
# -------------------------------------- #
# 
#_For each of the vulnerability component (selectivity and retention), the following need to be specified:
# ************************************** #
#_Component periods: Number of component time periods
#_Sex specific component: 0 = No; 1 = Yes
#_Vulnerability types
#_-> <0 = Mirror vulnerability component
#_-> 0 = Nonparameric component (one parameter per class)
#_-> 1 = Nonparameric component (one parameter per class, constant from last specified class)
#_-> 2 = Logistic component (inflection point and slope)
#_-> 3 = Logistic component (50% and 95% selection)
#_-> 4 = Double normal component (3 parameters)
#_-> 5 = Flat equal to one (1 parameter; phase must be negative)
#_-> 6 = Flat equal to zero (1 parameter; phase must be negative) 
#_-> 7 = Flat-topped double normal component (4 parameters) 
#_-> 8 = Declining logistic component with initial values (50% and 95% selection plus extra) 
#_-> 9 = Cubic-spline (specified with knots and values at knots) 
#_-> 10 = One parameter logistic component (inflection point and slope) 
#_Is the fleet within another? (0 = No; 1 = Yes)
#_Extra parameters for each pattern - 1 line per sex
# 
#_Is the maximum selectivity at size forced to equal 1 or not ?
# ************************************** #
 
#_The number of columns corresponds to the number of fleets (fisheries and surveys)
# Selectivity
#  Gear-1 | Gear-2#  Pot_Fishery | Survey_1 
1 1 #_Number of selectivity time period per fleet
0 0 #_Sex specific selectivity
2 0 #_Selectivity type
0 0 #_Insertion of fleet in another
0 0 #_Extra parameter for each pattern
# 
#_Retention
#_Gear-1 | Gear-2#_Pot_Fishery_| Survey_1
1 1 #_Number of Retention time period per fleet
0 0 #_Sex specific Retention
2 6 #_Selectivity type
1 0 #_retention flag (0 = No, 1 = Yes)
0 0 #_Extra parameter for each pattern
1 0 #_Selectivity for the maximum size class if forced to be 1?
 
# ====================================== #
# ====================================== #

#_Selectivity parameter controls
# ************************************** #
#_For each parameter (for each gear) columns are:
#_Fleet: The index of the fleet (positive for capture selectivity)
#_Index: Parameter count
#_Par_no: Parameter count within the current pattern
#_Sex: 0 = both; 1 = male; 2 = female
#_Init_val: Initial value for the parameter (must lie between lower and upper bounds)
#_Lower_Bd & Upper_Bd: Range for the parameter
#_Available prior types:
#_-> 0 = Uniform   - parameters are the range of the uniform prior
#_-> 1 = Normal    - parameters are the mean and sd
#_-> 2 = Lognormal - parameters are the mean and sd of the log
#_-> 3 = Beta      - parameters are the two beta parameters [see dbeta]
#_-> 4 = Gamma     - parameters are the two gamma parameters [see dgamma]
#_p1; p2: priors
#_Phase: Set equal to a negative number not to estimate
#_Start / End block: years to define the current block structure
#_Env_Link: Is there any environmental link for this parameter (0 = no; 1 = yes)
#_Link_Par: If 'Env_Link'=1; indicate the link to the environmental parameter
#_(i.e., which parameter (column) in the Envdata matrix)
#_Rand_Walk: Is there a random walk (0/1/2)- If so (1/2), which type :
#_1 = First order autoregressive process; 2 = gaussian white noise
#_Start_RdWalk / End_RdWalk: years (start/end) to define the period for random walk deviations
#_Sigma_RdWalk: sigma for the random walk
# ************************************** #
 
#_Fleet_| Index_| Par_no_| Sex_| Init_val_| Lower_Bd_| Upper_Bd_| Prior_type_| p1_| p2_| Phase_| Start_Block_| End_Block_| Env_Link_| Link_Par_| Rand_Walk_| Start_RdWalk_| End_RdWalk_| Sigma_RdWalk
# Pot_Fishery  
1 1 1 1 105.0093 5 186 0 1 999 -4 1989 2019 0 0 0 0 0 0  
1 2 2 1 4.754445 0.01 20 0 1 999 -4 1989 2019 0 0 0 0 0 0  
# Survey_1  
2 3 1 1 0.02041417 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 4 2 1 0.02041405 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 5 3 1 0.02041931 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 6 4 1 0.02057452 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 7 5 1 0.02107175 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 8 6 1 0.02435004 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 9 7 1 0.03433907 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 10 8 1 0.05449176 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 11 9 1 0.0836207 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 12 10 1 0.1284335 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 13 11 1 0.1518434 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 14 12 1 0.1833407 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 15 13 1 0.2828371 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 16 14 1 0.400038 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 17 15 1 0.3470329 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 18 16 1 0.2282325 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 19 17 1 0.1858602 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 20 18 1 0.1735583 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 21 19 1 0.1852795 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 22 20 1 0.2260137 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 23 21 1 0.2529056 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  
2 24 22 1 0.2724963 1e-05 1 0 0 999 -1 2015 2019 0 0 0 0 0 0  

#_Retention parameter controls
# ************************************** #
#_For each parameter (for each gear) columns are:
#_Fleet: The index of the fleet (negative for retention)
#_Index: Parameter count
#_Par_no: Parameter count within the current pattern
#_Sex: 0 = both; 1 = male; 2 = female
#_Init_val: Initial value for the parameter (must lie between lower and upper bounds)
#_Lower_Bd & Upper_Bd: Range for the parameter
#_Available prior types:
#_-> 0 = Uniform   - parameters are the range of the uniform prior
#_-> 1 = Normal    - parameters are the mean and sd
#_-> 2 = Lognormal - parameters are the mean and sd of the log
#_-> 3 = Beta      - parameters are the two beta parameters [see dbeta]
#_-> 4 = Gamma     - parameters are the two gamma parameters [see dgamma]
#_p1; p2: priors
#_Phase: Set equal to a negative number not to estimate
#_Start / End block: years to define the current block structure
#_Env_Link: Is there any environmental link for this parameter (0 = no; 1 = yes)
#_Link_Par: If 'Env_Link'=1; indicate the link to the environmental parameter
#_(i.e., which parameter (column) in the Envdata matrix)
#_Rand_Walk: Is there a random walk (0/1/2)- If so (1/2), which type :
#_1 = First order autoregressive process; 2 = gaussian white noise
#_Start_RdWalk / End_RdWalk: years (start/end) to define the period for random walk deviations
#_Sigma_RdWalk: sigma for the random walk
# ************************************** #
 
#_Fleet_| Index_| Par_no_| Sex_| Init_val_| Lower_Bd_| Upper_Bd_| Prior_type_| p1_| p2_| Phase_| Start_Block_| End_Block_| Env_Link_| Link_Par_| Rand_Walk_| Start_RdWalk_| End_RdWalk_| Sigma_RdWalk
# Pot_Fishery  
-1 25 1 1 0.8627194 1 190 1 96 10 -4 1989 2019 0 0 0 0 0 0  
-1 26 2 1 0.5927754 0.001 20 0 1 999 -4 1989 2019 0 0 0 0 0 0  
# Survey_1  
-2 27 1 0 0.999999 1 999 0 1 999 -3 2015 2019 0 0 0 0 0 0  

#_Number of asymptotic retention parameter
0 
#_Asymptotic parameter controls
# ************************************** #
#_Fleet: The index of the fleet (negative for retention)
#_Sex: 0 = both; 1 = male; 2 = female
#_Year: year of interest 
#_Init_val: Initial value for the parameter (must lie between lower and upper bounds)
#_Lower_Bd & Upper_Bd: Range for the parameter
#_Phase: Set equal to a negative number not to estimate
# ************************************** #
#_Fleet_| Sex_| Year_| Init_val_| Lower_Bd_| Upper_Bd_| Phase 
# -------------------------------------- #


#_Environmental parameters Control
# ************************************** #
#_Init_val: Initial value for the parameter (must lie between lower and upper bounds)
#_Lower_Bd & Upper_Bd: Range for the parameter
#_Phase: Set equal to a negative number not to estimate
#_Init_val_| Lower_Bd_| Upper_Bd_| Phase 
# 
#_One line for each parameter ordered as the parameters are in the
#_control matrices
# ************************************** #

#_Vulnerability impact#_Init_val_| Lower_Bd_| Upper_Bd_| Phase 
# -------------------------------------- #

#_Deviation parameter phase for the random walk in vulnerability parameters
#_Need to be defined
-1 

# -------------------------------------- #
## Priors for catchability
# -------------------------------------- #
 
# ************************************** #
# Init_val: Initial value for the parameter (must lie between lower and upper bounds)
# Lower_Bd & Upper_Bd: Range for the parameter
# Phase: Set equal to a negative number not to estimate
# Available prior types:
# -> 0 = Uniform   - parameters are the range of the uniform prior
# -> 1 = Normal    - parameters are the mean and sd
# -> 2 = Lognormal - parameters are the mean and sd of the log
# -> 3 = Beta      - parameters are the two beta parameters [see dbeta]
# -> 4 = Gamma     - parameters are the two gamma parameters [see dgamma]
# p1; p2: priors
# Q_anal: Do we need to solve analytically Q? (0 = No; 1 = Yes)
# CV_mult: multiplier ofr the input survey CV
# Loglik_mult: weight for the likelihood
# ************************************** #
# Init_val | Lower_Bd | Upper_Bd | Phase | Prior_type | p1 | p2 | Q_anal | CV_mult | Loglik_mult
0.9999 0.01 1 -5 0 0.843136 0.03 0 1 1
# -------------------------------------- #

# -------------------------------------- #
## Additional CV controls
# -------------------------------------- #
 
# ************************************** #
# Init_val: Initial value for the parameter (must lie between lower and upper bounds)
# Lower_Bd & Upper_Bd: Range for the parameter
# Phase: Set equal to a negative number not to estimate
# Available prior types:
# -> 0 = Uniform   - parameters are the range of the uniform prior
# -> 1 = Normal    - parameters are the mean and sd
# -> 2 = Lognormal - parameters are the mean and sd of the log
# -> 3 = Beta      - parameters are the two beta parameters [see dbeta]
# -> 4 = Gamma     - parameters are the two gamma parameters [see dgamma]
# p1; p2: priors
# ************************************** #
# Init_val | Lower_Bd | Upper_Bd | Phase | Prior_type| p1 | p2
1e-04 1e-05 10 -4 0 1 100
 
# Additional variance control for each survey (0 = ignore; >0 = use)
0 
# -------------------------------------- #

# -------------------------------------- #
## Penalties for the average fishing mortality rate
# -------------------------------------- #
 
# ************************************** #
# Fishing mortality controls
# ************************************** #
# Mean_F_male: mean male fishing mortality (base value for the fully-selected F) #
# Female_Offset: Offset between female and male fully-selected F  #
# Pen_std_Ph1 & Pen_std_Ph2: penalties on the fully-selected F during the early and later phase, respectively  #
# Ph_Mean_F_male & Ph_Mean_F_female: Phases to estimate the fishing mortality for males and females, respectively #
# Low_bd_mean_F & Up_bd_mean_F: Range for the mean fishing mortality (lower and upper bounds, respectivly) #
# Low_bd_Y_male_F & Up_bd_Y_male_F: Range for the male fishing mortality (lower and upper bounds, respectivly) #
# Low_bd_Y_female_F & Up_bd_Y_female_F: Range for the female fishing mortality (lower and upper bounds, respectivly)#
# ************************************** #
#  Mean_F_male | Female_Offset | Pen_std_Ph1 | Pen_std_Ph2 | Ph_Mean_F_male | Ph_Mean_F_female | Low_bd_mean_F | Up_bd_mean_F | Low_bd_Y_male_F | Up_bd_Y_male_F | Low_bd_Y_female_F | Up_bd_Y_female_F 
1 0.0505 0.5 45.5 1 1 -12 4 -10 10 -10 10
0 0 2 20 -1 -1 -12 4 -10 10 -10 10
# -------------------------------------- #

# -------------------------------------- #
## Size composition data control
# -------------------------------------- #
 
# ************************************** #
# Available types of likelihood:
# -> 0 = Ignore size-composition data in model fitting
# -> 1 = Multinomial with estimated/fixed sample size
# -> 2 = Robust approximation to multinomial
# -> 5 = Dirichlet
# Auto tail compression (pmin):
# -> pmin is the cumulative proportion used in tail compression
# Type-like prediction (1 = catch-like predictions; 2 = survey-like predictions)
# Lambda: multiplier for the effective sample size
# Emphasis: multiplier for weighting the overall likelihood
# ************************************** #
 
# The number of columns corresponds to the number size-composition data frames
2 2 2 2 # Type of likelihood for the size-composition
0 0 0 0 # Option for the auto tail compression
1 1 1 1 # Initial value for effective sample size multiplier
-4 -4 -4 -4 # Phase for estimating the effective sample size
1 1 1 2 # Composition appender (Should data be aggregated?)
1 1 1 2 # Type-like predictions
1 1 1 1 # Lambda: multiplier for the effective sample size
1 1 1 1 # Emphasis: multiplier for weighting the overall likelihood
# -------------------------------------- #

# -------------------------------------- #
## Time-varying Natural mortality controls
# -------------------------------------- #
 
# ************************************** #
# Available types of M specification:
# -> 0 = Constant natural mortality
# -> 1 = Random walk (deviates constrained by variance in M)
# -> 2 = Cubic Spline (deviates constrained by nodes & node-placement)
# -> 3 = Blocked changes (deviates constrained by variance at specific knots)
# -> 4 = Natural mortality is estimated as an annual deviation
# -> 5 = Deviations in M are estimated for specific periods relatively to the M estimated in the first year of the assessment
# -> 6 = Deviation in M are estimated for specific periods relatively to M during the current year
# ************************************** #
# Type of natural mortality
0 
# Is female M relative to M male?
# 0: No (absolute); 1: Yes (relative) 

# Phase of estimation
-4 
# Standard deviation in M deviations
0 
# Number of nodes for cubic spline or number of step-changes for option 3
# -> One line per sex
0
# Year position of the knots for each sex (vector must be equal to the number of nodes)
# -> One line per sex

# number of breakpoints in M by size
0 
# Size positions of breakpoints in M by size class
 
# Specific initial value for natural mortality deviations
0 
# Natural mortality deviation controls
# ************************************** #
# Init_val: Initial value for the parameter (must lie between lower and upper bounds)
# Lower_Bd & Upper_Bd: Range for the parameter
# Phase: Set equal to a negative number not to estimate
# Size_spec: Are the deviations size-specific ? (integer that specifies which size-class (negative to be considered))
# ************************************** #
# Init_val | Lower_Bd | Upper_Bd | Phase | Size_spec

# -------------------------------------- #

# -------------------------------------- #
## Tagging controls
# -------------------------------------- #
# Emphasis (likelihood weight) on tagging
1 
# -------------------------------------- #

# -------------------------------------- #
##  Immature/mature natural mortality 
# -------------------------------------- #
# maturity specific natural mortality? ( 0 = No; 1 = Yes - only for use if nmature > 1)
1 
# immature/mature natural mortality controls
# ************************************** #
# Init_val: Initial value for the parameter (must lie between lower and upper bounds)
# Lower_Bd & Upper_Bd: Range for the parameter
# Phase: Set equal to a negative number not to estimate
# Available prior types:
# -> 0 = Uniform   - parameters are the range of the uniform prior
# -> 1 = Normal    - parameters are the mean and sd
# -> 2 = Lognormal - parameters are the mean and sd of the log
# -> 3 = Beta      - parameters are the two beta parameters [see dbeta]
# -> 4 = Gamma     - parameters are the two gamma parameters [see dgamma]
# p1; p2: priors
# ************************************** #
# Init_val | Lower_Bd | Upper_Bd | Phase | Prior_type| p1 | p2
0.133531392624523 -4 4 4 1 0 0.05
# -------------------------------------- #

# -------------------------------------- #
## Other (additional) controls
# -------------------------------------- #
# First year of recruitment estimation deviations
1989 
# Last year of recruitment estimation deviations
2019 
# Consider terminal molting? (0 = No; 1 = Yes
1 
# Phase for recruitment estimation
1 
# Phase for recruitment sex-ratio estimation
2 
# Initial value for expected sex-ratio
0.5 
# Phase for initial recruitment estimation
-3 
# Verbose flag (0 = off; 1 = on; 2 = objective function; 3 = diagnostics)
1 
# Initial conditions (1 = unfished, 2 = steady-state, 3 = free params, 4 = free params revised)
3 
# Proportion of mature male biomass for SPR reference points
1 
# Stock-Recruit-Relationship (0 = none, 1 = Beverton-Holt) 
0 
# Maximum phase (stop the estimation after this phase)
10 
# Maximum number of function calls
-1 
# Calculate reference points (0 = No, 1 = Yes)
1 
# Use years specified to computed average sex ratio in the calculation of average recruitment for reference points
# -> 0 = No, i.e. Rec based on End year; 1 = Yes 
0 
# Years to compute equilibrium
200 
# -------------------------------------- #

# -------------------------------------- #
## Emphasis factor (weights for likelihood) controls
# -------------------------------------- #
# Weights on catches for the likelihood component
1 1 1 

# Penalties on deviations
# ************************************** #
#  Fdev_total | Fdov_total | Fdev_year | Fdov_year 
1 1 0 0
0 0 0 0

# Account for priors (penalties)
# ************************************** #
10000 	#_ Log_fdevs 
0 	#_ meanF 
0 	#_ Mdevs 
0 	#_ Rec_devs 
0 	#_ Initial_devs 
0 	#_ Fst_dif_dev 
0 	#_ Mean_sex-Ratio 
0 	#_ Molt_prob 
0 	#_ Free_selectivity 
0 	#_ Init_n_at_len 
0 	#_ Fvecs 
0 	#_ Fdovs 
0 	#_ Vul_devs 

# -------------------------------------- #

# -------------------------------------- #
## End of control file
# -------------------------------------- #
9999
