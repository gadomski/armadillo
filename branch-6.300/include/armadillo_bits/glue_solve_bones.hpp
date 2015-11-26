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



class glue_solve
  {
  public:
  
  inline static uword encode_flags(const solve_opts& settings);
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve>& X);
  
  template<typename eT, typename T1, typename T2> inline static bool solve(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags);
  
  template<typename eT, typename T1, typename T2> inline static bool solve_pinv(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr);
  
  
  
  private:
  
  // the flags below are internal to Armadillo, and are subject to change without notice.
  // DO NOT USE THESE FLAGS IN USER CODE!
  
  static const uword flag_fallback    = (1u << 0);
  static const uword flag_equilibrate = (1u << 1);
  static const uword flag_refine      = (1u << 2);
  static const uword flag_rankdef     = (1u << 3);
  static const uword flag_sympd       = (1u << 4);
  static const uword flag_symu        = (1u << 5);
  static const uword flag_syml        = (1u << 6);
  static const uword flag_triu        = (1u << 7);
  static const uword flag_tril        = (1u << 8);
  };



class glue_solve_tr
  {
  public:
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_tr>& X);
  };



//! @}
