// Copyright (C) 2009-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup fn_solve
//! @{



//! Solve a system of linear equations, i.e., A*X = B, where X is unknown.
//! For a square matrix A, this function is conceptually the same as X = inv(A)*B,
//! but is done more efficiently.
//! The number of rows in A and B must be the same.
//! B can be either a column vector or a matrix.
//! This function will also try to provide approximate solutions
//! to under-determined as well as over-determined systems (non-square A matrices).

template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve> >::result
solve
  (
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts&                      settings = solve_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  const uword flags = glue_solve::encode_flags(settings);
  
  return Glue<T1, T2, glue_solve>(A.get_ref(), B.get_ref(), flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve> >::result
solve
  (
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_solve>(A.get_ref(), B.get_ref());
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve> >::result
solve
  (
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_solve>(A.get_ref(), B.get_ref());
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve_tr> >::result
solve
  (
  const Op<T1, op_trimat>&               A,
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts&                      settings = solve_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  // TODO: process settings, with upper_tri or lower_tri set to true
  
  return Glue<T1, T2, glue_solve_tr>(A.m, B.get_ref(), A.aux_uword_a);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve_tr> >::result
solve
  (
  const Op<T1, op_trimat>&               A,
  const Base<typename T1::elem_type,T2>& B,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_solve_tr>(A.m, B.get_ref(), A.aux_uword_a);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve_tr> >::result
solve
  (
  const Op<T1, op_trimat>&               A,
  const Base<typename T1::elem_type,T2>& B,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_solve_tr>(A.m, B.get_ref(), A.aux_uword_a);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts&                      settings = solve_opts()
  )
  {
  arma_extra_debug_sigprint();
  
  const uword flags = glue_solve::encode_flags(settings);
  
  return glue_solve::solve(out, A.get_ref(), B.get_ref(), flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  return glue_solve::solve(out, A.get_ref(), B.get_ref());
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  return glue_solve::solve(out, A.get_ref(), B.get_ref());
  }



//! @}
