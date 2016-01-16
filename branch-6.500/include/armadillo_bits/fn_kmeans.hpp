// Copyright (C) 2008-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_kmeans
//! @{


#if defined(ARMA_BAD_COMPILER)


template<typename T1>
inline
bool
kmeans(Mat<typename T1::elem_type>& means, const Base<typename T1::elem_type,T1>&, const uword, const gmm_seed_mode&, const uword, const bool)
  {
  arma_extra_debug_sigprint();
  
  arma_stop("kmeans(): unsupported/inadequate compiler");
  
  means.reset();
  
  return false;
  }


#else


// TODO: restriction on element type
// TODO: allow complex numbers?
template<typename T1>
inline
bool
kmeans(Mat<typename T1::elem_type>& means, const Base<typename T1::elem_type,T1>& data, const uword k, const gmm_seed_mode& seed_mode, const uword n_iter, const bool print_mode)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   U(data.get_ref());
  const Mat<eT>& X = U.M;
  
  uvec valid_indices;
  
  const uword kk = (seed_mode == keep_existing) ? means.n_cols : k;
  
  if(seed_mode == keep_existing)
    {
    arma_debug_check( (means.n_rows != X.n_rows), "kmeans(): dimensionality mismatch between given means and data" );
    
    if( X.is_empty() || (kk == 0) )  { return true; }
    
    gmm_priv::gmm_diag<eT> gmm(means.n_rows, kk);
    
    gmm.set_means(means);
    
    const bool status = gmm.learn(X, kk, eucl_dist, keep_existing, n_iter, 0, eT(0), print_mode);
    
    if(status == false)  { return false; }
    
    means = gmm.means;
    
    valid_indices = find( gmm.hefts > std::numeric_limits<eT>::min() );
    }
  else
    {
    if(kk = 0)  { means.set_size(X.n_rows,0); return true; };
    
    gmm_priv::gmm_diag<eT> gmm;
    
    const bool status = gmm.learn(X, kk, eucl_dist, seed_mode, n_iter, 0, eT(0), print_mode);
    
    if(status == false)  { means.reset(); return false; }
    
    means = gmm.means;
    
    valid_indices = find( gmm.hefts > std::numeric_limits<eT>::min() );
    }
  
  
  // keep only the valid means
  
  if(valid_indices.n_elem == kk)  { return true; }
  if(valid_indices.n_elem == 0 )  { means.set_size(X.n_rows,0); return false; }
  
  means = means.cols(valid_indices);
  
  return true;
  }


#endif




//! @}
