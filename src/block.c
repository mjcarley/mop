/* interp.c
 * 
 * Copyright (C) 2007 Michael Carley
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

/* #include <stdio.h> */
#include <math.h>

#include <glib.h>

#include "mopblock.h"

#include <gsl/gsl_specfunc.h>

gdouble block_multipower(gdouble *xp, gint i, gint n, gint dim,
			 gint *p)

{
  gdouble f ;
  gint j ;

  for ( (f = 1.0), (j = 0) ; j < dim ; j ++ ) f *= xp[p[j]*n*dim+i*dim+j] ;

  return f ;
}

gint block_powers(gdouble *x, gint n, gint dim,
		  gint pmax, gdouble *xp)

{
  gint i, j ;

  i = 0 ; 
  for ( j = 0 ; j < n*dim ; j ++ ) xp[i*n*dim+j] = 1.0 ;

  for ( i = 1 ; i <= pmax ; i ++ )
    for ( j = 0 ; j < n*dim ; j ++ ) 
      if ( isnan(xp[i*n*dim+j] = x[j]*xp[(i-1)*n*dim+j]) )
	g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
	      "%s: NaN error at i=%d, j=%d", __FUNCTION__, i, j) ;

  return 0 ;
}

gint block_derivatives(gdouble *x, gint n,
		       gint dim, gint pmax, gint *d, 
		       gdouble *xp)

{
  gint i, j, k ;
  
  for ( i = 0 ; i < dim ; i ++ ) {
    for ( j = 0 ; j < d[i] ; j ++ ) {
      for ( k = 0 ; k < n ; k ++ ) 
	xp[j*n*dim+k*dim+i] = 0.0 ;
    }
    j = d[i] ;
    for ( k = 0 ; k < n ; k ++ ) xp[j*n*dim+k*dim+i] = 1.0 ; 
    for ( j = d[i]+1 ; j <= pmax ; j ++ ) {
      for ( k = 0 ; k < n ; k ++ )
	if ( isnan(xp[j*n*dim+k*dim+i] = 
		   xp[(j-1)*n*dim+k*dim+i]*x[k*dim+i]) )
	  g_log(G_LOG_DOMAIN, G_LOG_LEVEL_ERROR,
		"%s: NaN error at i=%d, j=%d, k=%d", 
		__FUNCTION__, i, j, k) ;
    }
  }

  for ( i = 0 ; i < dim ; i ++ ) {
    for ( j = d[i] ; j <= pmax ; j ++ ) {
      for ( k = 0 ; k < n ; k ++ )
	xp[j*n*dim+k*dim+i] *= gsl_sf_fact(j)/gsl_sf_fact(j-d[i]) ;
    }
  }

  return 0 ;
}
