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
  
  template<typename eT, typename T1, typename T2> inline static bool solve(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags);
  
  template<typename eT, typename T1, typename T2> inline static bool solve_pinv(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr);
  };




namespace solve_opts
  {
  struct opts
    {
    const uword flags;
    
    inline explicit opts(const uword in_flags);
    
    inline const opts operator|(const opts& rhs) const;
    };
  
  inline
  opts::opts(const uword in_flags)
    : flags(in_flags)
    {}
  
  inline
  const opts
  opts::operator|(const opts& rhs) const
    {
    const opts result( flags | rhs.flags );
    
    return result;
    }
  
  // The values below (eg. 1u << 1) are for only internal Armadillo use.
  // The values can change without notice.
  
  static const uword flag_none        = uword(0      );
  static const uword flag_nofallback  = uword(1u << 0);
  static const uword flag_equilibrate = uword(1u << 1);
  static const uword flag_refine      = uword(1u << 2);
  static const uword flag_rankdef     = uword(1u << 3);
  static const uword flag_symu        = uword(1u << 4);
  static const uword flag_syml        = uword(1u << 5);
  static const uword flag_triu        = uword(1u << 6);
  static const uword flag_tril        = uword(1u << 7);
  
  struct opts_none        : public opts { inline opts_none()        : opts(flag_none       ) {} };
  struct opts_nofallback  : public opts { inline opts_nofallback()  : opts(flag_nofallback ) {} };
  struct opts_equilibrate : public opts { inline opts_equilibrate() : opts(flag_equilibrate) {} };
  struct opts_refine      : public opts { inline opts_refine()      : opts(flag_refine     ) {} };
  struct opts_rankdef     : public opts { inline opts_rankdef()     : opts(flag_rankdef    ) {} };
  struct opts_symu        : public opts { inline opts_symu()        : opts(flag_symu       ) {} };
  struct opts_syml        : public opts { inline opts_syml()        : opts(flag_syml       ) {} };
  struct opts_triu        : public opts { inline opts_triu()        : opts(flag_triu       ) {} };
  struct opts_tril        : public opts { inline opts_tril()        : opts(flag_tril       ) {} };
  
  static const opts_none        none;
  static const opts_nofallback  nofallback;
  static const opts_equilibrate equilibrate;
  static const opts_refine      refine;
  static const opts_rankdef     rankdef;
  static const opts_symu        symu;
  static const opts_syml        syml;
  static const opts_triu        triu;
  static const opts_tril        tril;
  }



//! @}
