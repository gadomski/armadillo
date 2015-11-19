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
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> A = X.A.get_ref();
  
  const unwrap_check<T2> U(X.B.get_ref(), out);
  const Mat<eT>& B =     U.M;
  
  arma_debug_check( (A.n_rows != B.n_rows), "solve(): number of rows in the given matrices must be the same" );
  
  const bool status = glue_solve::solve_robust( out, A, B, (X.aux_uword == 1) );
  
  if(status == false)
    {
    out.reset();
    arma_bad("solve(): solution not found");
    }
  }



template<typename eT>
inline
bool
glue_solve::solve_plain(Mat<eT>& out, Mat<eT>& A, const Mat<eT>& B, const bool slow)
  {
  arma_extra_debug_sigprint();
  
  bool status = false;
  
  if(A.n_rows == A.n_cols)
    {
    status = auxlib::solve(out, A, B, slow);
    }
  else
  if(A.n_rows > A.n_cols)
    {
    arma_extra_debug_print("solve(): detected over-determined system");
    status = auxlib::solve_od(out, A, B);
    }
  else
    {
    arma_extra_debug_print("solve(): detected under-determined system");
    status = auxlib::solve_ud(out, A, B);
    }
  
  return status;
  }



template<typename eT>
inline
bool
glue_solve::solve_robust(Mat<eT>& out, Mat<eT>& A, const Mat<eT>& B, const bool slow)
  {
  arma_extra_debug_sigprint();
  
  bool status = false;
  
  if(A.n_rows == A.n_cols)
    {
    status = auxlib::solve(out, A, B, slow);
    
    if(status == false)
      {
      arma_debug_warn("solve(): matrix appears ill-conditioned; using pseudo-inverse to obtain approximate solution");
      
      status = glue_solve::solve_pinv(out, A, B);
      }
    }
  else
  if(A.n_rows > A.n_cols)
    {
    arma_extra_debug_print("solve(): detected over-determined system");
    status = auxlib::solve_od(out, A, B);
    }
  else
    {
    arma_extra_debug_print("solve(): detected under-determined system");
    status = auxlib::solve_ud(out, A, B);
    }
  
  return status;
  }



template<typename eT>
inline
bool
glue_solve::solve_pinv(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT> Ai;
  
  bool status = pinv(Ai, A);  // use default tolerance
  
  if(status == false)  { return false; }
  
  out = Ai * B;  // TODO: measure quality of the approximate solution?
  
  return true;
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
