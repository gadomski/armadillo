// Copyright (C) 2010-2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


//! \addtogroup glue_conv
//! @{


//! rudimentary implementation of 1D convolution operation

template<typename eT>
inline
void
glue_conv::apply_noalias(Mat<eT>& out, const Mat<eT>& A, const Mat<eT>& B, const bool A_is_col)
  {
  arma_extra_debug_sigprint();
  
  const Mat<eT>& h = (A.n_elem <= B.n_elem) ? A : B;
  const Mat<eT>& x = (A.n_elem <= B.n_elem) ? B : A;
  
  const uword   h_n_elem = h.n_elem;
  const uword   x_n_elem = x.n_elem;
  const uword out_n_elem = h_n_elem + x_n_elem - 1;
  
  if( (h_n_elem == 0) || (x_n_elem == 0) )  { out.reset(); return; }
  
  (A_is_col) ? out.set_size(out_n_elem, 1) : out.set_size(1, out_n_elem);
  
  const eT*   h_mem = h.memptr();
  const eT*   x_mem = x.memptr();
        eT* out_mem = out.memptr();
  
  
  for(uword out_i = 0; out_i < (h_n_elem-1); ++out_i)
    {
    eT acc = eT(0);
    
    uword h_i = out_i;
    
    for(uword x_i = 0; x_i <= out_i; ++x_i, --h_i)
      {
      acc += h_mem[h_i] * x_mem[x_i];
      }
    
    out_mem[out_i] = acc;
    }
  
  
  for(uword out_i = h_n_elem-1; out_i < out_n_elem - (h_n_elem-1); ++out_i)
    {
    eT acc = eT(0);
   
    uword h_i = h_n_elem - 1;
    
    for(uword x_i = out_i - h_n_elem + 1; x_i <= out_i; ++x_i, --h_i)
      {
      acc += h_mem[h_i] * x_mem[x_i];
      }
      
    out_mem[out_i] = acc;
    }
  
  
  for(uword out_i = out_n_elem - (h_n_elem-1); out_i < out_n_elem; ++out_i)
    {
    eT acc = eT(0);
    
    uword h_i = h_n_elem - 1;
    
    for(uword x_i = out_i - h_n_elem + 1; x_i < x_n_elem; ++x_i, --h_i)
      {
      acc += h_mem[h_i] * x_mem[x_i];
      }
    
    out_mem[out_i] = acc;
    }
  }



template<typename T1, typename T2>
inline
void
glue_conv::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_conv>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UA(X.A);
  const quasi_unwrap<T2> UB(X.B);
  
  arma_debug_check
    (
    ( ((UA.M.is_vec() == false) && (UA.M.is_empty() == false)) || ((UB.M.is_vec() == false) && (UB.M.is_empty() == false)) ),
    "conv(): given object is not a vector"
    );
  
  const bool A_is_col = ((T1::is_col) || (UA.M.n_cols == 1));
  
  if(UA.is_alias(out) || UB.is_alias(out))
    {
    Mat<eT> tmp;
    
    glue_conv::apply_noalias(tmp, UA.M, UB.M, A_is_col);
    
    out.steal_mem(tmp);
    }
  else
    {
    glue_conv::apply_noalias(out, UA.M, UB.M, A_is_col);
    }
  }



///


//! rudimentary implementation of 2D convolution operation

template<typename T1, typename T2>
inline
void
glue_conv2::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_conv2>& expr)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const quasi_unwrap<T1> UA(expr.A);
  const quasi_unwrap<T2> UB(expr.B);
  
  const Mat<eT>& A = UA.M;
  const Mat<eT>& B = UB.M;
  
  const Mat<eT>& G = (A.n_elem <= B.n_elem) ? A : B;   // unflipped filter coefficients
  const Mat<eT>& W = (A.n_elem <= B.n_elem) ? B : A;   // original 2D image
  
  const uword out_n_rows = ((W.n_rows > 0) || (G.n_rows > 0)) ? (W.n_rows + G.n_rows - 1) : uword(0);
  const uword out_n_cols = ((W.n_cols > 0) || (G.n_cols > 0)) ? (W.n_cols + G.n_cols - 1) : uword(0);
  
  out.zeros( out_n_rows, out_n_cols );
  
  if(G.is_empty() || W.is_empty())  { return; }
  
  Mat<eT> H(G.n_rows, G.n_cols);  // flipped filter coefficients
  
  const uword H_n_rows = H.n_rows;
  const uword H_n_cols = H.n_cols;
  
  const uword H_n_rows_m1 = H_n_rows - 1;
  const uword H_n_cols_m1 = H_n_cols - 1;
  
  for(uword col=0; col < H_n_cols; ++col)
    {
          eT* H_colptr = H.colptr(H_n_cols_m1 - col);
    const eT* G_colptr = G.colptr(col);
    
    for(uword row=0; row < H_n_rows; ++row)
      {
      H_colptr[H_n_rows_m1 - row] = G_colptr[row];
      }
    }
  
  Mat<eT> X( (W.n_rows + 2*(H_n_rows - 1)), (W.n_cols + 2*(H_n_cols - 1)), fill::zeros );
  
  X( H_n_rows-1, H_n_cols-1, size(W) ) = W;  // zero padded version of 2D image
  
  for(uword col=0; col < out_n_cols; ++col)
    {
    eT* out_colptr = out.colptr(col);
    
    for(uword row=0; row < out_n_rows; ++row)
      {
      // out.at(row, col) = accu( H % X(row, col, size(H)) );
      out_colptr[row] = accu( H % X.submat(row, col, (row + H_n_rows_m1), (col + H_n_cols_m1)) );
      }
    }
  }



//! @}
