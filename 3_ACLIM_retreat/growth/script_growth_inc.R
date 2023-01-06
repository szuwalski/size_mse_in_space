growth_trans  <- function(sizes,m_last , m_first, L_first, L_last, binsize,beta_scale){

  #L_first=sizes[1]
  #L_last=sizes[n_p]
  #m_first= growth_par[k,1]
  #m_last =growth_par[k,2]
  growmat <- matrix(0, nrow = n_p, ncol = n_p)
  for (i in 1:n_p){
    mean <- sizes[i] + m_first + (m_last-m_first)* (sizes[i]-L_first)/(L_last- L_first)
    alpha = mean/beta_scale
    for (k in i:n_p){
      growmat[i,k] = pgamma((sizes[k]+binsize/2),alpha,scale=beta_scale)-pgamma((sizes[k]-binsize/2),alpha,scale=beta_scale)
    }
  }
  growmat <- growmat/rowSums(growmat, na.rm = T)
  return(growmat)
}


