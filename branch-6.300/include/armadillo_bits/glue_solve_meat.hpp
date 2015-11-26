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
  
  const bool sym         = (flags & flag_sym        );
  const bool sympd       = (flags & flag_sympd      );
  const bool tril        = (flags & flag_tril       );
  const bool triu        = (flags & flag_triu       );
  const bool equilibrate = (flags & flag_equilibrate);
  const bool refine      = (flags & flag_refine     );
  const bool rankdef     = (flags & flag_rankdef    );
  const bool fallback    = (flags & flag_fallback   );
  
  
  arma_extra_debug_print("enabled flags:");
  
  if(sym        )  { arma_extra_debug_print("sym");         }
  if(sympd      )  { arma_extra_debug_print("sympd");       }
  if(tril       )  { arma_extra_debug_print("tril");        }
  if(triu       )  { arma_extra_debug_print("triu");        }
  if(equilibrate)  { arma_extra_debug_print("equilibrate"); }
  if(refine     )  { arma_extra_debug_print("refine");      }
  if(rankdef    )  { arma_extra_debug_print("rankdef");     }
  if(fallback   )  { arma_extra_debug_print("fallback");    }
  
  
  bool status = false;
  
  if(rankdef)
    {
    const Proxy<T1> PA(A_expr.get_ref());
    
    if(PA.get_n_rows() == PA.get_n_cols())
      {
      arma_extra_debug_print("solve(): detected square system");
      
      status = glue_solve::solve_pinv(out, PA.Q, B_expr);
      }
    else
      {
      arma_extra_debug_print("solve(): detected non-square system");
      
      arma_debug_check( (sym || sympd || tril || triu), "solve(): incorrect options for non-square matrix" );
      
      Mat<eT> A = PA.Q;
      status = auxlib::solve_nonsquare_ext(out, A, B_expr, equilibrate);
      }
    
    return status;
    }
  
  Mat<eT> A = A_expr.get_ref();
  
  if(A.n_rows == A.n_cols)
    {
    arma_extra_debug_print("solve(): detected square system");
    
         if(equilibrate || refine)  { status = auxlib::solve_square_ext(out, A, B_expr, equilibrate); }
    else if(sympd                )  { status = auxlib::solve_sympd     (out, A, B_expr);              }
    else if(sym                  )  { status = auxlib::solve_sym       (out, A, B_expr);              }
    
    if( (status == false) && (fallback) )
      {
      status = glue_solve::solve_pinv(out, A_expr, B_expr);  // using A_expr as A can be overwritten at this stage
      }
    }
  else
    {
    arma_extra_debug_print("solve(): detected non-square system");
    
    arma_debug_check( (sym || sympd || tril || triu), "solve(): incorrect options for non-square matrix" );
    
    status = auxlib::solve_nonsquare(out, A, B_expr);
    
    if( (status == false) && (fallback) )
      {
      status = auxlib::solve_nonsquare_ext(out, A, B_expr, equilibrate);
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



template<typename T1, typename T2>
inline
void
glue_solve_tr::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_tr>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap_check<T1> A_tmp(X.A, out);
  const unwrap_check<T2> B_tmp(X.B, out);
  
  const Mat<eT>& A = A_tmp.M;
  const Mat<eT>& B = B_tmp.M;
  
  bool  err_state = false;
  char* err_msg   = 0;
  
  arma_debug_set_error( err_state, err_msg, ((&A) == (&B)),           "solve(): A is an alias of B" );
  arma_debug_set_error( err_state, err_msg, (A.n_rows != B.n_rows),   "solve(): number of rows in A and B must be the same" );
  arma_debug_set_error( err_state, err_msg, (A.is_square() == false), "solve(): A is not a square matrix" );
  
  arma_debug_check(err_state, err_msg);
  
  const bool status = auxlib::solve_tr(out, A, B, X.aux_uword);
  
  if(status == false)
    {
    out.reset();
    arma_bad("solve(): solution not found");
    }
  }



//! @}
