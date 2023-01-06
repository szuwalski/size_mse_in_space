pars_LHP_setting <-
  function (pref_hab,
            pars_pref_hab,
            x_omega,
            x_epsilon,
            scale_g,
            LHP,
            ad_eff,
            LHP_spatial,
            pars_LHP,
            n_grpar,
            Years_climsc) {
    
    pars_LHP_setting <- list(
      "lon" = lon,
      "lat" = lat,
      #"clim_sc_test" = clim_sc_test,
      "pref_hab" = pref_hab,
      "Years_climsc" = Years_climsc,
      "pars_pref_hab" = pars_pref_hab,
      "x_omega" = x_omega,
      "x_epsilon" = x_epsilon,
      "binsize" = binsize,
      "sizes" = sizes,
      "scale_g" = scale_g,
      "LHP" = LHP,
      "ad_eff" = ad_eff,
      "LHP_spatial" = LHP_spatial,
      "pars_LHP" = pars_LHP,
      "n_grpar" = n_grpar,
      "year_LHP" = NULL,
      "clim_sc" = clim_sc
    )
    
    return(pars_LHP_setting)
  }
