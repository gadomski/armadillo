// Copyright (C) 2009-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup glue_solve
//! @{



template<typename T1, typename T2>
inline
void
glue_solve::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_solve::solve( out, X.A, X.B, X.aux_uword );
  
  if(status == false)
    {
    out.reset();
    arma_bad("solve(): solution not found");
    }
  }



template<typename eT, typename T1, typename T2>
inline
bool
glue_solve::solve(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags)
  {
  arma_extra_debug_sigprint();
  
  const bool nofallback  = (flags & solve_opt::flag_nofallback );
  const bool equilibrate = (flags & solve_opt::flag_equilibrate);
  const bool refine      = (flags & solve_opt::flag_refine     );
  const bool rankdef     = (flags & solve_opt::flag_rankdef    );
  const bool symu        = (flags & solve_opt::flag_symu       );
  const bool syml        = (flags & solve_opt::flag_syml       );
  const bool triu        = (flags & solve_opt::flag_triu       );
  const bool tril        = (flags & solve_opt::flag_tril       );
  
  
  arma_extra_debug_print("enabled flags:");
  
  if(nofallback )  { arma_extra_debug_print("nofallback");  }
  if(equilibrate)  { arma_extra_debug_print("equilibrate"); }
  if(refine     )  { arma_extra_debug_print("refine");      }
  if(rankdef    )  { arma_extra_debug_print("rankdef");     }
  if(symu       )  { arma_extra_debug_print("symu");        }
  if(syml       )  { arma_extra_debug_print("syml");        }
  if(triu       )  { arma_extra_debug_print("triu");        }
  if(tril       )  { arma_extra_debug_print("tril");        }
  
  bool status = false;
  
  if(rankdef)
    {
    const Proxy<T1> PA(A_expr.get_ref());
    
    if(PA.get_n_rows() == PA.get_n_cols())
      {
      arma_extra_debug_print("solve(): detected square system");
      
      // TODO: BUG: matrix A doesn't have the right form if {symu, syml, triu, tril} is given
      status = glue_solve::solve_pinv(out, PA.Q, B_expr);
      }
    else
      {
      arma_extra_debug_print("solve(): detected non-square system");
      
      arma_debug_check( (symu || syml || triu || tril), "solve(): incorrect options for non-square matrix" );
      
      Mat<eT> A = PA.Q;
      status = auxlib::solve_nonsquare_ext(out, A, B_expr);
      }
    
    return status;
    }
  
  if(symu || syml)
    {
    // ensure matrix is square
    // check for equlibrate and refine
    // ...
    // check for nofallback 
    }
  else
  if(triu || tril)
    {
    // ensure matrix is square
    // check for equlibrate and refine
    // ...
    // check for nofallback 
    }
  else
    {
    // general matrix
    // - square matrix
    // - nonsquare matrix
    // check for nofallback 
    }
  
  Mat<eT> A = A_expr.get_ref();
  
  if(A.n_rows == A.n_cols)
    {
    arma_extra_debug_print("solve(): detected square system");
    
    if(equilibrate)
      {
      // TODO: BUG: matrix A doesn't have the right form if {symu, syml, triu, tril} is given
      status = auxlib::solve_square_ext(out, A, B_expr, true);
      }
    else
      {
           if(symu) { status = auxlib::solve_sym(out, A, B_expr, uword(0));    }
      else if(syml) { status = auxlib::solve_sym(out, A, B_expr, uword(1));    }
      else if(triu) { status = auxlib::solve_tri(out, A, B_expr, uword(0));    }  // NOTE: solve_tri() doesn't overwrite A
      else if(tril) { status = auxlib::solve_tri(out, A, B_expr, uword(1));    }
      }
    
    if( (status == false) && (nofallback == false) )
      {
      status = glue_solve::solve_pinv(out, A_expr, B_expr);  // using A_expr as A can be overwritten at this stage
      }
    }
  else
    {
    arma_extra_debug_print("solve(): detected non-square system");
    
    arma_debug_check( (equilibrate), "solve(): option 'equilibrate' not supported for non-square matrix" );
    
    arma_debug_check( (refine), "solve(): option 'refine' not supported for non-square matrix" );
    
    arma_debug_check( (symu || syml || triu || tril), "solve(): options imply a square matrix, but non-square matrix was given" );
    
    status = auxlib::solve_nonsquare(out, A, B_expr);
    
    if( (status == false) && (nofallback == false) )
      {
      status = auxlib::solve_nonsquare_ext(out, A, B_expr);
      }
    }
  
  return status;
  }



template<typename eT, typename T1, typename T2>
inline
bool
glue_solve::solve_pinv(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  // brutal last-resort approximate solver based on pseudo-inverse
  
  Mat<eT> Ai;
  
  bool status = pinv(Ai, A_expr.get_ref());  // use default tolerance
  
  if(status == false)  { return false; }
  
  out = Ai * B_expr.get_ref();
  
  return out.is_finite();
  }



//! @}
