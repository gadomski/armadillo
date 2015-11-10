// Copyright (C) 2015 National ICT Australia (NICTA)
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// -------------------------------------------------------------------
// 
// Written by Conrad Sanderson - http://conradsanderson.id.au


#include <armadillo>
#include "catch.hpp"

using namespace arma;


TEST_CASE("fn_symmatl_1")
  {
  mat A = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745   0.051408;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153   0.035437;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317  -0.454499;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040   0.373833;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768   0.258704;\
    ";
  
  mat B = symmatu( A(0,0,size(5,5)) );
  mat C = symmatl( A(0,0,size(5,5)) );
  
  mat BB = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745;\
     0.201990   0.058956  -0.149362  -0.045465   0.296153;\
     0.019678  -0.149362   0.314156   0.419733   0.068317;\
    -0.493936  -0.045465   0.419733  -0.393139  -0.135040;\
    -0.126745   0.296153   0.068317  -0.135040  -0.353768;\
    ";
  
  mat CC = 
    "\
     0.061198   0.437242  -0.492474   0.336352   0.239585;\
     0.437242   0.058956  -0.031309   0.411541  -0.428913;\
    -0.492474  -0.031309   0.314156   0.458476  -0.406953;\
     0.336352   0.411541   0.458476  -0.393139  -0.291020;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768;\
    ";
  
  REQUIRE( accu(abs( B - BB )) == Approx(0.0) );
  REQUIRE( accu(abs( C - CC )) == Approx(0.0) );
  
  mat X;
  REQUIRE_THROWS( X = symmatu(A) ); // symmatu() and symmatl() currently handle only square matrices
  }


// TODO: test complex number version, with and without conjugation