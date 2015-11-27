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
  const solve_opts::opts&                opts = solve_opts::none
  )
  {
  arma_extra_debug_sigprint();
  
  return Glue<T1, T2, glue_solve>(A.get_ref(), B.get_ref(), opts.flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve> >::result
solve
  (
  const Op<T1, op_symmat>&               A,  // TODO: complex matrices use cx_symmat
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts::opts&                opts = solve_opts::none
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = opts.flags;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_symu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_syml; }
  
  return Glue<T1, T2, glue_solve>(A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve> >::result
solve
  (
  const Op<T1, op_trimat>&               A,
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts::opts&                opts = solve_opts::none
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = opts.flags;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return Glue<T1, T2, glue_solve>(A.m, B.get_ref(), flags);
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
  
  return Glue<T1, T2, glue_solve>(A.get_ref(), B.get_ref(), solve_opts::flag_none);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve> >::result
solve
  (
  const Op<T1, op_symmat>&               A,
  const Base<typename T1::elem_type,T2>& B,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_symu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_syml; }
  
  return Glue<T1, T2, glue_solve>(A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve> >::result
solve
  (
  const Op<T1, op_trimat>&               A,
  const Base<typename T1::elem_type,T2>& B,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return Glue<T1, T2, glue_solve>(A.m, B.get_ref(), flags);
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
  
  return Glue<T1, T2, glue_solve>(A.get_ref(), B.get_ref(), solve_opts::flag_none);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve> >::result
solve
  (
  const Op<T1, op_symmat>&               A,
  const Base<typename T1::elem_type,T2>& B,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_symu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_syml; }
  
  return Glue<T1, T2, glue_solve>(A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, const Glue<T1, T2, glue_solve> >::result
solve
  (
  const Op<T1, op_trimat>&               A,
  const Base<typename T1::elem_type,T2>& B,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return Glue<T1, T2, glue_solve>(A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const Base<typename T1::elem_type,T1>& A,
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts::opts&                opts = solve_opts::none
  )
  {
  arma_extra_debug_sigprint();
  
  return glue_solve::solve(out, A.get_ref(), B.get_ref(), opts.flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const   Op<T1, op_symmat>&             A,
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts::opts&                opts = solve_opts::none
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = opts.flags;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_symu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_syml; }
  
  return glue_solve::solve(out, A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const   Op<T1, op_trimat>&             A,
  const Base<typename T1::elem_type,T2>& B,
  const solve_opts::opts&                opts = solve_opts::none
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = opts.flags;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return glue_solve::solve(out, A.m, B.get_ref(), flags);
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
  
  return glue_solve::solve(out, A.get_ref(), B.get_ref(), solve_opts::flag_none);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const   Op<T1, op_symmat>&             A,
  const Base<typename T1::elem_type,T2>& B,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_symu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_syml; }
  
  return glue_solve::solve(out, A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const   Op<T1, op_trimat>&             A,
  const Base<typename T1::elem_type,T2>& B,
  const bool   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return glue_solve::solve(out, A.m, B.get_ref(), flags);
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
  
  return glue_solve::solve(out, A.get_ref(), B.get_ref(), solve_opts::flag_none);
  }


template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const   Op<T1, op_symmat>&             A,
  const Base<typename T1::elem_type,T2>& B,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_symu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_syml; }
  
  return glue_solve::solve(out, A.m, B.get_ref(), flags);
  }



template<typename T1, typename T2>
inline
typename enable_if2< is_supported_blas_type<typename T1::elem_type>::value, bool >::result
solve
  (
         Mat<typename T1::elem_type>&    out,
  const   Op<T1, op_trimat>&             A,
  const Base<typename T1::elem_type,T2>& B,
  const char*   // argument kept only for compatibility with old user code
  )
  {
  arma_extra_debug_sigprint();
  
  uword flags = solve_opts::flag_none;
  
  if(A.aux_uword_a == 0)  {  flags |= solve_opts::flag_triu; }
  if(A.aux_uword_a == 1)  {  flags |= solve_opts::flag_tril; }
  
  return glue_solve::solve(out, A.m, B.get_ref(), flags);
  }






//! @}
