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
  
  const bool equilibrate = bool(flags & solve_opts::flag_equilibrate);
  const bool refine      = bool(flags & solve_opts::flag_refine     );
  const bool approx      = bool(flags & solve_opts::flag_approx     );
  const bool noapprox    = bool(flags & solve_opts::flag_noapprox   );
  const bool rankdef     = bool(flags & solve_opts::flag_rankdef    );
  
  arma_extra_debug_print("glue_solve_gen::apply(): enabled flags:");
  
  if(equilibrate)  { arma_extra_debug_print("equilibrate"); }
  if(refine     )  { arma_extra_debug_print("refine");      }
  if(approx     )  { arma_extra_debug_print("approx");      }
  if(noapprox   )  { arma_extra_debug_print("noapprox");    }
  if(rankdef    )  { arma_extra_debug_print("rankdef");     }
  
  
  bool status = false;
  
  Mat<eT> A = A_expr.get_ref();
  
  if(A.n_rows == A.n_cols)
    {
    arma_extra_debug_print("glue_solve_gen::apply(): detected square system");
    
    if(equilibrate || refine)
      {
      arma_extra_debug_print("glue_solve_gen::apply(): (equilibrate || refine)");
      
      status = auxlib::solve_square_ext(out, A, B_expr, equilibrate);  // A is overwritten
      }
    else
      {
      arma_extra_debug_print("glue_solve_gen::apply(): (standard)");
      
      status = auxlib::solve_square(out, A, B_expr.get_ref());  // A is overwritten
      }
    
    if( (status == false) && (approx == true) )
      {
      arma_extra_debug_print("glue_solve_gen::apply(): solving rank deficient system");
      
      arma_debug_warn("system appears singular to machine precision; attempting approximate solution");
      
      Mat<eT> AA = A_expr.get_ref();
      status = auxlib::solve_approx(out, AA, B_expr.get_ref());  // A is overwritten
      }
    }
  else
    {
    arma_extra_debug_print("glue_solve_gen::apply(): detected non-square system");
    
    if(equilibrate)  { arma_debug_warn( "solve(): option 'equilibrate' ignored for non-square matrix" ); }
    if(refine     )  { arma_debug_warn( "solve(): option 'refine' ignored for non-square matrix"      ); }
    
    status = auxlib::solve_approx(out, A, B_expr.get_ref());  // A is overwritten
    }
  
  return status;
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
  
  const bool equilibrate = bool(flags & solve_opts::flag_equilibrate);
  const bool refine      = bool(flags & solve_opts::flag_refine     );
  const bool approx      = bool(flags & solve_opts::flag_approx     );
  const bool noapprox    = bool(flags & solve_opts::flag_noapprox   );
  const bool rankdef     = bool(flags & solve_opts::flag_rankdef    );
  const bool symu        = bool(flags & solve_opts::flag_symu       );
  const bool syml        = bool(flags & solve_opts::flag_syml       );
  
  arma_extra_debug_print("glue_solve_sym::apply(): enabled flags:");
  
  if(equilibrate)  { arma_extra_debug_print("equilibrate"); }
  if(refine     )  { arma_extra_debug_print("refine");      }
  if(approx     )  { arma_extra_debug_print("approx");      }
  if(noapprox   )  { arma_extra_debug_print("noapprox");    }
  if(rankdef    )  { arma_extra_debug_print("rankdef");     }
  if(symu       )  { arma_extra_debug_print("symu");        }
  if(syml       )  { arma_extra_debug_print("syml");        }
  
  
  bool status = false;
  
  const uword layout = (symu) ? uword(0) : uword(1);
  
  if(equilibrate)
    {
    arma_extra_debug_print("glue_solve_sym::apply(): (equilibrate)");
    
    Mat<eT> symA = (symu) ? symmatu( A_expr.get_ref() ) : symmatl( A_expr.get_ref() );
    
    status = auxlib::solve_square_ext(out, symA, B_expr.get_ref(), true);  // A is overwritten
    }
  else
  if(refine)
    {
    arma_extra_debug_print("glue_solve_sym::apply(): (refine)");
    
    const unwrap_check<T1> U(A_expr.get_ref(), out);
    
    const Mat<eT>& A = U.M;
    
    arma_debug_check( (A.is_square() == false), "symmatu()/symmatl(): given matrix must be square sized" );
    
    status = auxlib::solve_sym_ext(out, A, B_expr.get_ref(), layout);  // A is not modified
    }
  else
    {
    arma_extra_debug_print("glue_solve_sym::apply(): (standard)");
    
    Mat<eT> A = A_expr.get_ref();
    
    arma_debug_check( (A.is_square() == false), "symmatu()/symmatl(): given matrix must be square sized" );
    
    status = auxlib::solve_sym(out, A, B_expr.get_ref(), layout);  // A is overwritten
    }
  
  
  if( (status == false) && (approx == true) )
    {
    arma_extra_debug_print("glue_solve_sym::apply(): solving rank deficient system");
    
    Mat<eT> symA = (symu) ? symmatu( A_expr.get_ref() ) : symmatl( A_expr.get_ref() );
    
    status = auxlib::solve_approx(out, symA, B_expr.get_ref());  // A is overwritten
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
  
  const bool equilibrate = bool(flags & solve_opts::flag_equilibrate);
  const bool refine      = bool(flags & solve_opts::flag_refine     );
  const bool approx      = bool(flags & solve_opts::flag_approx     );
  const bool noapprox    = bool(flags & solve_opts::flag_noapprox   );
  const bool rankdef     = bool(flags & solve_opts::flag_rankdef    );
  const bool triu        = bool(flags & solve_opts::flag_triu       );
  const bool tril        = bool(flags & solve_opts::flag_tril       );
  
  arma_extra_debug_print("glue_solve_tri::apply(): enabled flags:");
  
  if(equilibrate)  { arma_extra_debug_print("equilibrate"); }
  if(refine     )  { arma_extra_debug_print("refine");      }
  if(approx     )  { arma_extra_debug_print("approx");      }
  if(noapprox   )  { arma_extra_debug_print("noapprox");    }
  if(rankdef    )  { arma_extra_debug_print("rankdef");     }
  if(triu       )  { arma_extra_debug_print("triu");        }
  if(tril       )  { arma_extra_debug_print("tril");        }
  
  
  bool status = false;
  
  if(equilibrate || refine)
    {
    arma_extra_debug_print("glue_solve_tri::apply(): (equilibrate || refine)");
    
    Mat<eT> triA = (triu) ? trimatu( A_expr.get_ref() ) : trimatl( A_expr.get_ref() );
    
    status = auxlib::solve_square_ext(out, triA, B_expr.get_ref(), equilibrate);  // A is overwritten
    }
  else
    {
    arma_extra_debug_print("glue_solve_tri::apply(): (standard)");
    
    const unwrap_check<T1> U(A_expr.get_ref(), out);
    
    const Mat<eT>& A = U.M;
    
    const uword layout = (triu) ? uword(0) : uword(1);
    
    status = auxlib::solve_tri(out, A, B_expr.get_ref(), layout);  // A is not modified
    }
  
  
  if( (status == false) && (approx == true) )
    {
    arma_extra_debug_print("glue_solve_tri::apply(): solving rank deficient system");
    
    Mat<eT> triA = (triu) ? trimatu( A_expr.get_ref() ) : trimatl( A_expr.get_ref() );
    
    status = auxlib::solve_approx(out, triA, B_expr.get_ref());  // A is overwritten
    }
  
  return status;
  }



//! @}
