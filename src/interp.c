/* interp.c
 * 
 * Copyright (C) 2007, 2021 Michael Carley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/**
 * @file   interp.c
 * @author Michael Carley
 * @date   Wed Nov 14 17:53:17 2007
 * 
 * @brief Interpolation and differentiation using multi-variable
 * orthogonal polynomials
 * 
 * @defgroup interp Interpolation and differentiation
 * @{
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <glib.h>

#include "mop.h"
#include "mopblock.h"

#ifdef _HAVE_CONFIG_H_
#include "config.h"
#endif

/** 
 * Calculate the coefficients of an expansion of a function \a f in
 * the polynomial \a p.
 * 
 * @param p ::mop_polynomial_t for the expansion;
 * @param f value(s) of function(s) at base points of \a p;
 * @param n number of values of \f$f\f$ at each base point;
 * @param c coefficients of expansion, 
 * \f$f\approx\sum_{i}c_{i}P_{i}(\mathbf{x})\f$;
 * @param w a ::mop_polynomial_workspace_t of appropriate size.
 * 
 * @return 0 on success.
 */

gint mop_polynomial_transform(mop_polynomial_t *p, gdouble *f, gint n,
			      gdouble *c, mop_polynomial_workspace_t *w)

{
  gdouble *ws ;
  gint i, j, k, m ;

  ws = w->block ;

  mop_polynomial_eval_base(p, ws) ;

  for ( i = 0 ; i < mop_polynomial_term_number(p)*n ; i ++ ) c[i] = 0.0 ;
  
  for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
    for ( m = 0 ; m < mop_polynomial_term_number(p) ; m ++ ) {
      j = mop_polynomial_index(p,m) ;
      for ( k = 0 ; k < n ; k ++ ) {
	if ( isnan(c[i*n+k] += f[j*n+k]*ws[m*mop_polynomial_term_number(p)+i]*
	  mop_polynomial_weight(p,j)) ) 
	  g_error("%s: NaN error at i=%d, j=%d, k=%d",
		  __FUNCTION__, i, j, k) ;
      }      
    }
  }  

  return 0 ;
}

/** 
 * Evaluate a function expanded in orthogonal polynomials,
 * \f$f\approx\sum_{i}c_{i}P_{i}(\mathbf{x})\f$.
 * 
 * @param p ::mop_polynomial_t for the expansion;
 * @param c coefficients of the expansion, from 
 * ::mop_polynomial_transform;
 * @param n number of function values, as in ::mop_polynomial_transform;
 * @param x coordinates of evaluation point;
 * @param f value(s) of function(s).
 * 
 * @return 0 on success.
 */

gint mop_interpolate(mop_polynomial_t *p, gdouble *c, gint n,
		     gdouble *x, gdouble *f)

{  
  gint i, j ;
  gdouble P ;

  block_powers(x, 1,
	       mop_polynomial_dimension(p),
	       mop_polynomial_order(p),
	       p->xp) ;
  
  for ( j = 0 ; j < n ; j ++ ) f[j] = 0.0 ;

  for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
    for ( (P = 0.0), (j = 0) ; j < mop_polynomial_term_number(p) ; j ++ )
      P += mop_polynomial_coefficient(p,i,j)*
	block_multipower(p->xp, 0, 1, 
			 mop_polynomial_dimension(p),
			 mop_polynomial_monomial(p,j)) ;
    for ( j = 0 ; j < n ; j ++ ) 
      f[j] += P*c[i*n+j] ;
  }

  return 0 ;
}

/** 
 * Compute the weights for direct interpolation at a point, based on
 * values at the base points of a system of orthogonal polynomials, so
 * that \f$f(\mathbf{x})\approx\sum v_{i}f_{i}\f$. Remember to use
 * ::mop_polynomial_index to map base points to real points.
 * 
 * @param p a ::mop_polynomial;
 * @param x interpolation point;
 * @param v interpolation weights;
 * @param w a suitably sized ::mop_polynomial_workspace.
 * 
 * @return 0 on success.
 */

gint mop_interpolation_weights(mop_polynomial_t *p, 
			       gdouble *x, gdouble *v,
			       mop_polynomial_workspace_t *w)

{
  gdouble *P, *Px ;
  gint i, j, n ;

  n = mop_polynomial_term_number(p) ;
  P = w->block ; Px = &(w->block[n*n]) ;

  mop_polynomial_eval(p, x, Px) ;
  mop_polynomial_eval_base(p, P) ;

  for ( j = 0 ; j < mop_polynomial_term_number(p) ; j ++ ) {    
    for ( (v[j] = 0.0), (i = 0) ; i < mop_polynomial_term_number(p) ; i ++ )
      v[j] += Px[i]*P[j*mop_polynomial_term_number(p)+i] ;
    v[j] *= mop_polynomial_weight(p,mop_polynomial_index(p,j)) ;
  }
  
  return 0 ;
}

/** 
 * Evaluate derivatives of a ::mop_polynomial_t at a point \a x. 
 * 
 * @param p ::mop_polynomial_t to evaluate;
 * @param x coordinates of evaluation point;
 * @param d orders of differentiation for each dimension;
 * @param P array containing values of
 * \f$\partial^{d_{1}+d_{2}+\dots}P(\mathbf{x})/\partial
 * x_{1}^{d_{1}}x_{2}^{d_{2}}\dots\f$ at \a x.
 * 
 * @return 0 on success.
 */

gint mop_polynomial_differentiate(mop_polynomial_t *p, gdouble *x,
				  gint *d, gdouble *P)

{
  gint i, j ;

  block_derivatives(x, 1,
		    mop_polynomial_dimension(p),
		    mop_polynomial_order(p),
		    d, p->xp) ;
  
  for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ )
    for ( (P[i] = 0.0), (j = 0) ; j < mop_polynomial_term_number(p) ; j ++ )
      P[i] += mop_polynomial_coefficient(p,i,j)*
	block_multipower(p->xp, 0, 1, 
			 mop_polynomial_dimension(p),
			 mop_polynomial_monomial(p,j)) ;

  return 0 ;
}

/** 
 * Compute the weights for direct differentiation at a point, based on
 * values at the base points of a system of orthogonal polynomials, so
 * that \f$\partial^{d_{1}+d_{2}+\dots}f(\mathbf{x})/\partial
 * x_{1}^{d_{1}}x_{2}^{d_{2}}\dots\approx\sum v_{i}f_{i}\f$. Remember
 * to use ::mop_polynomial_index to map base points to real points.
 * 
 * @param p a ::mop_polynomial;
 * @param x interpolation point;
 * @param d array of derivative orders;
 * @param v differentiation weights;
 * @param w a suitably sized ::mop_polynomial_workspace.
 * 
 * @return 0 on success.
 */


gint mop_differentiation_weights(mop_polynomial_t *p, 
				 gdouble *x, gint *d, gdouble *v,
				 mop_polynomial_workspace_t *w)


{
  gdouble *P, *Px, scp ;
  gint i, j, k, n ;

  n = mop_polynomial_point_number(p) ;
  P = w->block ; Px = &(w->block[n*n]) ;

  mop_polynomial_eval_base(p, P) ;
  mop_polynomial_differentiate(p, x, d, Px) ;

  for ( j = 0 ; j < mop_polynomial_term_number(p) ; j ++ ) {    
    for ( (v[j] = 0.0), (i = 0) ; i < mop_polynomial_term_number(p) ; i ++ )
      v[j] += Px[i]*P[j*mop_polynomial_term_number(p)+i] ;
    k = mop_polynomial_index(p,j) ;
    v[j] *= mop_polynomial_weight(p,k) ;
  }

  for ( (i = 0), (scp = 1.0) ; i < mop_polynomial_dimension(p) ; i ++ ) {
    for ( j = 0 ; j < d[i] ; j ++ ) scp *= p->sc ;
  }

  for ( j = 0 ; j < mop_polynomial_term_number(p) ; j ++ ) v[j] *= scp ;

  return 0 ;
}

/**
 * @}
 * 
 */
