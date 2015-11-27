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



//
// glue_solve_gen


template<typename T1, typename T2>
inline
void
glue_solve_gen::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_gen>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_solve_gen::apply( out, X.A, X.B, X.aux_uword );
  
  if(status == false)
    {
    out.reset();
    arma_bad("solve(): solution not found");
    }
  }



template<typename eT, typename T1, typename T2>
inline
bool
glue_solve_gen::apply(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags)
  {
  arma_extra_debug_sigprint();
  
  const bool nofallback  = (flags & solve_opts::flag_nofallback );
  const bool equilibrate = (flags & solve_opts::flag_equilibrate);
  const bool refine      = (flags & solve_opts::flag_refine     );
  const bool rankdef     = (flags & solve_opts::flag_rankdef    );
  
  arma_extra_debug_print("enabled flags:");
  
  if(nofallback )  { arma_extra_debug_print("nofallback");  }
  if(equilibrate)  { arma_extra_debug_print("equilibrate"); }
  if(refine     )  { arma_extra_debug_print("refine");      }
  if(rankdef    )  { arma_extra_debug_print("rankdef");     }
  
  
  bool status = false;
  
  Mat<eT> A = A_expr.get_ref();
  
  if(A.n_rows == A.n_cols)
    {
    arma_extra_debug_print("(square)");
    
    if(equilibrate || refine)
      {
      arma_extra_debug_print("(equilibrate || refine)");
      
      status = auxlib::solve_square_ext(out, A, B_expr, equilibrate);  // A is overwritten
      }
    else
      {
      arma_extra_debug_print("(standard)");
      
      status = auxlib::solve_square(out, A, B_expr.get_ref());  // A is overwritten
      }
    
    
    if( (status == false) && (nofallback == false) )
      {
      arma_extra_debug_print("(fallback)");
      
      status = glue_solve_gen::apply_pinv(out, A_expr.get_ref(), B_expr.get_ref());
      }
    }
  else
    {
    arma_extra_debug_print("(non-square)");
    
    arma_debug_check( (equilibrate || refine), "solve(): options 'equilibrate' and 'refine' not applicable to non-square matrices" );
    
    status = auxlib::solve_nonsquare(out, A, B_expr.get_ref());  // A is overwritten
    
    
    if( (status == false) && (nofallback == false) )
      {
      arma_extra_debug_print("(fallback)");
      
      status = auxlib::solve_nonsquare_ext(out, A, B_expr.get_ref());  // A is overwritten
      }
    }
  
  return status;
  }



template<typename eT, typename T1, typename T2>
inline
bool
glue_solve_gen::apply_pinv(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr)
  {
  arma_extra_debug_sigprint();
  
  // brutal last-resort approximate solver based on pseudo-inverse
  
  Mat<eT> Ai;
  
  bool status = pinv(Ai, A_expr.get_ref());  // use default tolerance
  
  if(status == false)  { return false; }
  
  out = Ai * B_expr.get_ref();
  
  return out.is_finite();
  }



//
// glue_solve_sym


template<typename T1, typename T2>
inline
void
glue_solve_sym::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_sym>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_solve_sym::apply( out, X.A, X.B, X.aux_uword );
  
  if(status == false)
    {
    out.reset();
    arma_bad("solve(): solution not found");
    }
  }



template<typename eT, typename T1, typename T2>
inline
bool
glue_solve_sym::apply(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags)
  {
  arma_extra_debug_sigprint();
  
  const bool nofallback  = (flags & solve_opts::flag_nofallback );
  const bool equilibrate = (flags & solve_opts::flag_equilibrate);
  const bool refine      = (flags & solve_opts::flag_refine     );
  const bool rankdef     = (flags & solve_opts::flag_rankdef    );
  const bool symu        = (flags & solve_opts::flag_symu       );
  const bool syml        = (flags & solve_opts::flag_syml       );
  
  arma_extra_debug_print("enabled flags:");
  
  if(nofallback )  { arma_extra_debug_print("nofallback");  }
  if(equilibrate)  { arma_extra_debug_print("equilibrate"); }
  if(refine     )  { arma_extra_debug_print("refine");      }
  if(rankdef    )  { arma_extra_debug_print("rankdef");     }
  if(symu       )  { arma_extra_debug_print("symu");        }
  if(syml       )  { arma_extra_debug_print("syml");        }
  
  
  bool status = false;
  
  const uword layout = (symu) ? uword(0) : uword(1);
    
  if(equilibrate)
    {
    arma_extra_debug_print("(equilibrate)");
    
    Mat<eT> symA = (symu) ? symmatu( A_expr.get_ref() ) : symmatl( A_expr.get_ref() );
    
    status = auxlib::solve_square_ext(out, symA, B_expr.get_ref(), true);  // A is overwritten
    }
  else
  if(refine)
    {
    arma_extra_debug_print("(refine)");
    
    const unwrap_check<T1> U(A_expr.get_ref(), out);
    
    const Mat<eT>& A = U.M;
    
    arma_debug_check( (A.is_square() == false), "symmatu()/symmatl(): given matrix must be square sized" );
    
    status = auxlib::solve_sym_ext(out, A, B_expr.get_ref(), layout);  // A is not modified
    }
  else
    {
    arma_extra_debug_print("(standard)");
    
    Mat<eT> A = A_expr.get_ref();
    
    arma_debug_check( (A.is_square() == false), "symmatu()/symmatl(): given matrix must be square sized" );
    
    status = auxlib::solve_sym(out, A, B_expr.get_ref(), layout);  // A is overwritten
    }
  
  
  if( (status == false) && (nofallback == false) )
    {
    arma_extra_debug_print("(fallback)");
    
    const Mat<eT> symA = (symu) ? symmatu( A_expr.get_ref() ) : symmatl( A_expr.get_ref() );
    
    status = glue_solve_gen::apply_pinv(out, symA, B_expr.get_ref());
    }
  
  return status;
  }



//
// glue_solve_tri


template<typename T1, typename T2>
inline
void
glue_solve_tri::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_solve_tri>& X)
  {
  arma_extra_debug_sigprint();
  
  const bool status = glue_solve_tri::apply( out, X.A, X.B, X.aux_uword );
  
  if(status == false)
    {
    out.reset();
    arma_bad("solve(): solution not found");
    }
  }



template<typename eT, typename T1, typename T2>
inline
bool
glue_solve_tri::apply(Mat<eT>& out, const Base<eT,T1>& A_expr, const Base<eT,T2>& B_expr, const uword flags)
  {
  arma_extra_debug_sigprint();
  
  const bool nofallback  = (flags & solve_opts::flag_nofallback );
  const bool equilibrate = (flags & solve_opts::flag_equilibrate);
  const bool refine      = (flags & solve_opts::flag_refine     );
  const bool rankdef     = (flags & solve_opts::flag_rankdef    );
  const bool triu        = (flags & solve_opts::flag_triu       );
  const bool tril        = (flags & solve_opts::flag_tril       );
  
  
  arma_extra_debug_print("enabled flags:");
  
  if(nofallback )  { arma_extra_debug_print("nofallback");  }
  if(equilibrate)  { arma_extra_debug_print("equilibrate"); }
  if(refine     )  { arma_extra_debug_print("refine");      }
  if(rankdef    )  { arma_extra_debug_print("rankdef");     }
  if(triu       )  { arma_extra_debug_print("triu");        }
  if(tril       )  { arma_extra_debug_print("tril");        }
  
  
  bool status = false;
  
  if(equilibrate || refine)
    {
    arma_extra_debug_print("(equilibrate || refine)");
    
    Mat<eT> triA = (triu) ? trimatu( A_expr.get_ref() ) : trimatl( A_expr.get_ref() );
    
    status = auxlib::solve_square_ext(out, triA, B_expr.get_ref(), equilibrate);  // A is overwritten
    }
  else
    {
    arma_extra_debug_print("(standard)");
    
    const unwrap_check<T1> U(A_expr.get_ref(), out);
    
    const Mat<eT>& A = U.M;
    
    const uword layout = (triu) ? uword(0) : uword(1);
    
    status = auxlib::solve_tri(out, A, B_expr.get_ref(), layout);  // A is not modified
    }
  
  
  if( (status == false) && (nofallback == false) )
    {
    arma_extra_debug_print("(fallback)");
    
    const Mat<eT> triA = (triu) ? trimatu( A_expr.get_ref() ) : trimatl( A_expr.get_ref() );
    
    status = glue_solve_gen::apply_pinv(out, triA, B_expr.get_ref());
    }
  
  return status;
  }



//! @}
