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
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve>& X);
  
  template<typename eT> inline static bool solve_plain(Mat<eT>& out, Mat<eT>& A, const Mat<eT>& B, const bool slow);
  
  template<typename eT> inline static bool solve_robust(Mat<eT>& out, Mat<eT>& A, const Mat<eT>& B, const bool slow);
  
  template<typename eT> inline static bool solve_pinv(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B);
  };



class glue_solve_tr
  {
  public:
  
  template<typename T1, typename T2> inline static void apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_tr>& X);
  };



//! @}
