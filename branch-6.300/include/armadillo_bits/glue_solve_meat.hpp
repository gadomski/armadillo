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
  
  const bool status = glue_solve::solve( out, X.A, X.B );
  
  if(status == false)
    {
    out.reset();
    arma_bad("solve(): solution not found");
    }
  }



template<typename eT, typename T1, typename T2>
inline
bool
glue_solve::solve(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  bool status = false;
  
  Mat<eT> A = A_expr.get_ref();
  
  if(A.n_rows == A.n_cols)
    {
    status = auxlib::solve_square(out, A, B_expr);
    
    if(status == false)
      {
      arma_debug_warn("solve(): attempting approximate solution via pseudo-inverse");
      
      status = glue_solve::solve_pinv(out, A_expr, B_expr);  // using A_expr as auxlib::solve_square() overwrites A
      }
    }
  else
    {
    arma_extra_debug_print("solve(): detected non-square system");
    status = auxlib::solve_nonsquare(out, A, B_expr);
    
    if(status == false)
      {
      status = auxlib::solve_nonsquare_rankdef(out, A, B_expr);
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



template<typename eT, typename T2>
inline
bool
glue_solve::solve_reinterpreted_inv(Mat<eT>& out, Mat<eT>& A, const Base<eT,T2>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_check( (A.is_square() == false), "inv(): given matrix must be square sized" );
  
  const bool status = auxlib::solve_square(out, A, B_expr);
  
  return status;
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
