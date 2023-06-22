## Stock assessment
#------------------
if(print_messages) print("Make Stock assessment")

if(SA == "spatialIPM"){
  
  if(display_print) print("spatial IPM")
  
  source(paste0(project_spatialIPM,"03_spatial_model/run_model_mse.R"))
  
}

if(SA == "nonspatialIPM"){
  
}

if(SA == "GMACS"){
  
  # Update data and control files
  source(file = file.path(dir_GMACS,
                          "functions", "Update_GMACS_file.R", fsep = fsep))
  source(file = file.path(dir_GMACS,
                          "functions", "Run_Gmacs.R", fsep = fsep))
  # Now run GMACS
  runGMACS()
  
}

if(SA == "none"){
  
  n_at_size = c()
  c_at_size = c()
  f_at_size = c()
  for(s in 1:length(sizes)){
    
    # c: catch, m: natural mortality, n: abundance
    # --> c = f / (f + m) * n * (1 - exp(-(f + m)))
    imm_male_M = 0.32 
    mat_male_M = 0.28
    m = 0.32
    
    
    
    ## If temproal resolution is month
    if(temp_resolution == "month"){
      
      c = sum(total_spatial_catch_nb[,,2,s,(t-9):t])
      
      ab_at_t = c()
      catch_at_t = c()
      ab_at_t_mat = matrix(0,nrow = length(lat),ncol = length(lon))
      catch_at_t_mat = matrix(0,nrow = length(lat),ncol = length(lon))
      for(t2 in (t-9):t){
        
        # Abundance
        ab_at_t0 = sum(imm_N_at_Len[,,2,s,t2] + mat_N_at_Len[,,2,s,t2])
        ab_at_t = c(ab_at_t,ab_at_t0)
        
        ab_at_t0_mat = imm_N_at_Len[,,2,s,t2] + mat_N_at_Len[,,2,s,t2]
        ab_at_t_mat = ab_at_t_mat + ab_at_t0_mat
        
        # Catch
        catch_at_t0 = sum(total_spatial_catch_nb[,,2,s,t2])
        catch_at_t = c(catch_at_t,catch_at_t0)
        
        catch_at_t0_mat = total_spatial_catch_nb[,,2,s,t2]
        catch_at_t_mat = catch_at_t_mat + catch_at_t0_mat
        
      }
      
      n = ab_at_t[1]
      n_mat = ab_at_t_mat / 9
      
    }else if(temp_resolution == "year"){
      
      c = sum(total_spatial_catch_nb[,,2,s,t])
      
      ab_at_t = c()
      catch_at_t = c()
      ab_at_t_mat = matrix(0,nrow = length(lat),ncol = length(lon))
      catch_at_t_mat = matrix(0,nrow = length(lat),ncol = length(lon))
      
      # Abundance
      ab_at_t0 = sum(imm_N_at_Len[,,2,s,t] + mat_N_at_Len[,,2,s,t])
      ab_at_t = c(ab_at_t,ab_at_t0)
      
      ab_at_t0_mat = imm_N_at_Len[,,2,s,t] + mat_N_at_Len[,,2,s,t]
      ab_at_t_mat = ab_at_t_mat + ab_at_t0_mat
      
      # Catch
      catch_at_t0 = sum(total_spatial_catch_nb[,,2,s,t])
      catch_at_t = c(catch_at_t,catch_at_t0)
      
      catch_at_t0_mat = total_spatial_catch_nb[,,2,s,t]
      catch_at_t_mat = catch_at_t_mat + catch_at_t0_mat
      
      n = ab_at_t[1]
      n_mat = ab_at_t_mat
      
    }
    
    # Aggregated value
    n_plus_1 = (n - c * exp(m/2)) / exp(m)  # Pope's approximation (see https://www.fao.org/3/x9026e/x9026e06.htm)
    f = log(n / n_plus_1) - m
    
    # # Spatial value
    # n_plus_1_mat = (n_mat - catch_at_t_mat * exp(m/2)) / exp(m)  # Pope's approximation (see https://www.fao.org/3/x9026e/x9026e06.htm)
    # f = log(n_mat / n_plus_1_mat) - m
    
    n_at_size = c(n_at_size,n)
    c_at_size = c(c_at_size,c)
    f_at_size = c(f_at_size,f)
    
  }
  
}
