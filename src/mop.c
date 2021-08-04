/* mop.c
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
 * @file   mop.c
 * @author Michael Carley
 * @date   Wed Nov 14 17:53:17 2007
 * 
 * @brief  Multi-variable orthogonal polynomials
 * 
 * @defgroup mop Multi-variable orthogonal polynomials
 * @{
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <glib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_specfunc.h>

#include "mop.h"
#include "mopblock.h"

#ifdef _HAVE_CONFIG_H_
#include "config.h"
#endif

#define SDS_TRANSPOSE_S

static void monomial_string(gchar *s, gint *mn, gint dim)

{
  switch ( dim ) {
  default: g_assert_not_reached() ; break ;
  case 1: sprintf(s, "%d", mn[0]) ; break ;
  case 2: sprintf(s, "%d,%d", mn[0], mn[1]) ; break ;
  case 3: sprintf(s, "%d,%d,%d", mn[0], mn[1], mn[2]) ; break ;
  }

  return ;
}

static gint sds_decomp(gsl_matrix *A, gsl_matrix *S, gsl_vector *D,
		       gsl_vector *v)

{
  gdouble tmp ;
  gint i, j, k, m, n ;

  if ( A->size1 != S->size1 || A->size2 != S->size2 ) 
    g_error("%s: matrices A (%lux%lu) and S (%lux%lu) not the same size",
	  __FUNCTION__, A->size1, A->size2, S->size1, S->size2) ;
  if ( A->size1 != A->size2 )
    g_error("%s: matrix A (%lux%lu) must be square",
	    __FUNCTION__, A->size1, A->size2) ;
  
  if ( D->size != A->size1 ) 
    g_error("%s: vector D (%lu) must have as many rows as matrix A (%lux%lu)",
	  __FUNCTION__, D->size, A->size1, A->size2) ;

  n = A->size1 ;

  gsl_vector_set_zero(v) ; gsl_vector_set_zero(D) ;
  gsl_matrix_set_zero(S) ;

  for ( j = 0 ; j < n ; j ++ ) {
    for ( i = 0 ; i <= j-1 ; i ++ ) {
#ifdef SDS_TRANSPOSE_S
      gsl_vector_set(v, i, gsl_matrix_get(S, i, j)*gsl_vector_get(D, i)) ;
#else
      gsl_vector_set(v, i, gsl_matrix_get(S, j, i)*gsl_vector_get(D, i)) ;
#endif /*SDS_TRANSPOSE_S*/
    }
    for ( (tmp = gsl_matrix_get(A, j, j)), (k = 0) ; k <= j-1 ; k ++ ) {
#ifdef SDS_TRANSPOSE_S
      tmp -= gsl_matrix_get(S, k, j)*gsl_vector_get(v,k) ;
#else
      tmp -= gsl_matrix_get(S, j, k)*gsl_vector_get(v,k) ;
#endif /*SDS_TRANSPOSE_S*/
    }
    if ( isnan(tmp) ) 
      g_error("%s: NaN error at j=%d", __FUNCTION__, j) ;
    if ( tmp <= 0.0 ) {
      g_warning("%s: D_j=%lg out of range at j=%d", __FUNCTION__, tmp, j) ;
      return MOP_FAILURE ;
    }
    gsl_vector_set(v, j, tmp) ; gsl_vector_set(D, j, tmp) ;
    
    for ( m = j+1 ; m < n ; m ++ ) {
      for ( (tmp = gsl_matrix_get(A, m, j)), (k = 0) ; k <= j-1 ; k ++ ) {
#ifdef SDS_TRANSPOSE_S
	tmp -= gsl_matrix_get(S, k, m)*gsl_vector_get(v, k) ;
#else
 	tmp -= gsl_matrix_get(S, m, k)*gsl_vector_get(v, k) ;
#endif /*SDS_TRANSPOSE_S*/
	if ( isnan(tmp) ) 
	  g_error("%s: NaN error at i=%d, j=%d, m=%d", __FUNCTION__, i, j, m) ;	
      }

#ifdef SDS_TRANSPOSE_S
      gsl_matrix_set(S, j, m, tmp/gsl_vector_get(v, j)) ;
#else
      gsl_matrix_set(S, m, j, tmp/gsl_vector_get(v, j)) ;
#endif /*SDS_TRANSPOSE_S*/
      if ( isnan(tmp/gsl_vector_get(v, j)) ) 
	g_error("%s: NaN error at j=%d, m=%d, v_j=%lg", 
		__FUNCTION__, j, m, gsl_vector_get(v, j)) ;	
    }
  }

  for ( i = 0 ; i < n ; i ++ ) gsl_matrix_set(S, i, i, 1.0) ;

  return MOP_SUCCESS ;
}

static gint matrix_rank(gsl_matrix *A, gdouble tol, 
			gsl_vector *tau, gsl_vector *norm)

{
  gsl_permutation *p ;
  gint sgn, rank, i, j, m, n ;

  m = A->size1 ; n = A->size2 ;
  p = gsl_permutation_alloc(n) ;

  gsl_linalg_QRPT_decomp(A, tau, p, &sgn, norm) ;

  for ( rank = i = 0 ; i < GSL_MIN(m,n) ; i ++ )
    for ( j = i ; j < GSL_MIN(m,n) ; j ++ )
      if ( fabs(gsl_matrix_get(A,i,j)) > tol ) {
	rank ++ ; break ;
      }

  gsl_permutation_free(p) ;

  return rank ;
}

static gint make_index_list(gint n, gint m, gint *id, gint *ni)

{
  gint i, j, p ;

  if ( m <= 0 || m > 3 ) 
    g_error("%s: m=%d out of range, 0 < m <= 3", __FUNCTION__, m) ;

  switch ( m ) {
  case 1:
    for ( i = 0 ; i <= n ; i ++ ) id[i] = i ;
    *ni = n ; 
    return MOP_SUCCESS ;
    break ;
  case 2:
    for ( i = (*ni) = 0 ; i <= n ; i ++ ) {
      for ( j = 0 ; j <= i ; j ++ ) {
	id[2*(*ni)+0] = i-j ; id[2*(*ni)+1] = j ;
	g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	      "%s: %d %d %d", __FUNCTION__, (*ni), 
	      id[2*(*ni)+0], id[2*(*ni)+1]) ;
	(*ni) ++ ;
      }
    }
    return MOP_SUCCESS ;
    break ;
  case 3:
    for ( p = (*ni) = 0 ; p <= n ; p ++ ) {
      for ( i = 0 ; i <= p ; i ++ ) 
	for ( j = 0 ; j <= p-i ; (j ++), ((*ni) ++) ) {
	  id[(*ni)*3+0] = p-i-j ; id[(*ni)*3+1] = j ; id[(*ni)*3+2] = i ; 
	}
    }
    return MOP_SUCCESS ;
    break ;
  }

  g_assert_not_reached() ;

  return MOP_SUCCESS ;
}

/* start of public functions */

/** 
 * The number of monomials required for a multi-variable polynomial of
 * a given order.
 * 
 * @param dim dimension of system;
 * @param order order of polynomial.
 * 
 * @return number of terms in general polynomial of order \a order in
 * \a dim dimensions.
 */

gint mop_number_of_terms(gint dim, gint order)

{
  switch (dim) {
  case 1: return (order+1) ; break ;
  case 2: return ((order+1)*(order+2)/2) ; break ;
  case 3: return ((order+1)*(order+2)*(order+3)/6) ; break ;
  default: break ;
  }

  g_error("%s: dimension %d is currently unsupported", __FUNCTION__, dim) ;
  
  return MOP_SUCCESS ;
}

/** 
 * Allocate a ::mop_polynomial.
 * 
 * @param np (maximum) number of points in basis;
 * @param dim dimension of system;
 * @param order maximum order of polynomial to be generated.
 * 
 * @return pointer to new ::mop_polynomial.
 */

mop_polynomial_t *mop_polynomial_alloc(gint np, gint dim, gint order)

{
  mop_polynomial_t *p ;

  if ( (dim < 1) || ( dim > 3) ) 
    g_error("%s: dimension %d is currently unsupported (sorry)", 
	  __FUNCTION__, dim) ;

  p = (mop_polynomial_t *)g_malloc(sizeof(mop_polynomial_t)) ;

  mop_polynomial_dimension(p) = dim ;
  mop_polynomial_point_number(p) = 0 ;
  mop_polynomial_point_number_max(p) = np ;
  mop_polynomial_order(p) = -1 ;
  mop_polynomial_order_max(p) = order ;
  mop_polynomial_term_number(p) = 0 ;
  mop_polynomial_term_number_max(p) = mop_number_of_terms(dim, order) ;

  mop_polynomial_powers(p) = (gint *)g_malloc(mop_polynomial_dimension(p)*
					      mop_polynomial_term_number_max(p)*
					      sizeof(gint)) ;
  mop_polynomial_indices(p) = 
    (gint *)g_malloc(mop_polynomial_point_number_max(p)*sizeof(gint)) ;

  mop_polynomial_points(p) = 
    (gdouble *)g_malloc(mop_polynomial_point_number_max(p)*
			mop_polynomial_dimension(p)*sizeof(gdouble)) ;
  mop_polynomial_weights(p) = 
    (gdouble *)g_malloc(mop_polynomial_point_number_max(p)*sizeof(gdouble)) ;
  p->R = (gdouble *)g_malloc((mop_polynomial_point_number_max(p)*
			      mop_polynomial_term_number_max(p))*
			     sizeof(gdouble)) ;
  p->xp = (gdouble *)g_malloc(mop_polynomial_point_number_max(p)*
			      mop_polynomial_term_number_max(p)*
			      mop_polynomial_dimension(p)*
			      sizeof(gdouble)) ;
  p->sc = 1.0 ;

  return p ;
}

/** 
 * Free a ::mop_polynomial_t and associated memory
 * 
 * @param p ::mop_polynomial_t to free.
 * 
 * @return MOP_SUCCESS on success.
 */

gint mop_polynomial_free(mop_polynomial_t *p)

{
  g_free(p->xp) ;
  g_free(p->R) ;
  g_free(mop_polynomial_powers(p)) ;
  g_free(mop_polynomial_indices(p)) ;
  g_free(mop_polynomial_points(p)) ; 
  g_free(mop_polynomial_weights(p)) ; 

/*   g_free(p) ; */

  return MOP_SUCCESS ;
}

/** 
 * Allocate a ::mop_polynomial_workspace_t for use in generating
 * orthogonal polynomials.
 * 
 * @param np (maximum) number of points in basis;
 * @param dim dimension of system;
 * @param order maximum order of polynomial to be generated.
 * 
 * @return pointer to new ::mop_polynomial_workspace.
 */

mop_polynomial_workspace_t *mop_polynomial_workspace_alloc(gint np,
							 gint dim,
							 gint order)

{
  mop_polynomial_workspace_t *w ;
  gint bsize ;

  w = (mop_polynomial_workspace_t *)
    g_malloc(sizeof(mop_polynomial_workspace_t)) ;
  
  w->npmax = np ; w->omax = order ; w->dim = dim ;
  bsize = np*np /*matrix of monomials*/
    + np*dim*(order+1) /*block of powers*/
    + np /*tau vector*/
    + np ; /*norm vector*/
  w->block = (gdouble *)g_malloc(4*bsize*sizeof(gdouble)) ;
  
  w->imax = mop_number_of_terms(dim, order)*dim ;
  w->iblock = (gint *)g_malloc(w->imax*sizeof(gint)) ;

  return w ;
}

/** 
 * Free a ::mop_polynomial_workspace_t and associated memory
 * 
 * @param w ::mop_polynomial_workspace_t to free.
 * 
 * @return MOP_SUCCESS on success.
 */

gint mop_polynomial_workspace_free(mop_polynomial_workspace_t *w)

{
  g_free(w->block) ; g_free(w->iblock) ;

  g_free(w) ;

  return MOP_SUCCESS ;
}

/** 
 * Set points and weights to be used in generating sets of orthogonal
 * polynomials. The data are copied into \a p so the arrays can be
 * reused. If \a w is NULL, unit weights are used.
 * 
 * @param p a ::mop_polynomial_t of appropriate size;
 * @param x array of points of the same dimension as \a p;
 * @param w array of weights, one for each \a x;
 * @param n number of points.
 * 
 * @return MOP_SUCCESS on success.
 */

gint mop_polynomial_set_points(mop_polynomial_t *p, 
			       gdouble *x, gdouble *w, gint n)

{
  gint i ;

  if ( mop_polynomial_point_number_max(p) < n )
    g_error("%s: too many points (%d) for polynomial "
	    "(maximum number of points %d)", __FUNCTION__, 
	    n, mop_polynomial_point_number_max(p)) ;

  g_memmove(mop_polynomial_points(p), x, 
	    n*mop_polynomial_dimension(p)*sizeof(gdouble)) ;
  if ( w != NULL ) {
    for ( i = 0 ; i < n ; i ++ )
      if ( (mop_polynomial_weight(p,i) = w[i]) <= 0.0 )
	g_error("%s: weight %d is not positive, %lg",
		__FUNCTION__, i, w[i]) ;
  } else
    for ( i = 0 ; i < n ; i ++ ) mop_polynomial_weight(p,i) = 1.0 ;
    
  mop_polynomial_point_number(p) = n ;

  for ( i = 0 ; i < n ; i ++ ) mop_polynomial_index(p,i) = i ;

  return MOP_SUCCESS ;
}

/** 
 * Write a ::mop_polynomial_t to a file.
 * 
 * @param p ::mop_polynomial_t to write;
 * @param f file pointer.
 * 
 * @return MOP_SUCCESS on success.
 */

gint mop_polynomial_write(mop_polynomial_t *p, FILE *f)

{
  gint i, j ;

  fprintf(f, "%d %d %d %d\n", 
	  mop_polynomial_dimension(p),
	  mop_polynomial_point_number(p),
	  mop_polynomial_order(p),
	  mop_polynomial_term_number(p)) ;

  for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
    for ( j = 0 ; j < mop_polynomial_dimension(p) ; j ++ ) 
      fprintf(f, " %d", mop_polynomial_monomial_power(p,i,j)) ;
    fprintf(f, "\n") ;
  }

  for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
    for ( j = 0 ; j < mop_polynomial_term_number(p) ; j ++ ) {
      fprintf(f, " %lg", mop_polynomial_coefficient(p,i,j)) ;
    }
    fprintf(f, "\n") ;
  }

  return MOP_SUCCESS ;
}

/** 
 * Generate the basis for a ::mop_polynomial_t by selecting sufficient
 * monomial powers to match the number of points in the
 * ::mop_polynomial, using the method of Xu, Yuan, 2004, `On discrete
 * orthogonal polynomials of several variables', Advances in Applied
 * Mathematics, 33:615--632, doi:10.1016/j.aam.2004.03.002.
 * 
 * @param p a ::mop_polynomial, which should have been initialized by
 * a call to ::mop_polynomial_set_points;
 * @param order maximum order of monomial to be employed;
 * @param tol tolerance to be used in determining rank of basis matrix; 
 * @param w suitably sized ::mop_polynomial_workspace.
 * 
 * @return MOP_SUCCESS on success or MOP_NO_BASIS if no basis can be
 * generated on the points and weights supplied.
 */

gint mop_polynomial_basis_power(mop_polynomial_t *p, gint order,
				gdouble tol, mop_polynomial_workspace_t *w)

{
  gint i, j, k, rank, ni, *mon, offset ;
  gdouble *powers ;
  gsl_matrix_view A ;
  gsl_vector_view tau, norm ;
  gchar mns[32] ;

  if ( mop_polynomial_order_max(p) < order )
    g_error("%s: order too high (%d) for polynomial (maximum order %d)", 
	    __FUNCTION__, order, mop_polynomial_order_max(p)) ;

  if ( w->npmax < mop_polynomial_point_number(p) )
    g_error("%s: workspace is not big enough (%d points) "
	    "for polynomial (%d points)", 
	    __FUNCTION__, w->npmax, mop_polynomial_point_number(p)) ;

  if ( w->dim < mop_polynomial_dimension(p) )
    g_error("%s: workspace dimension (%d) smaller than polynomial's (%d)", 
	    __FUNCTION__, w->dim, mop_polynomial_dimension(p)) ;

  if ( w->omax < order )
    g_error("%s: workspace maximum order (%d) "
	  "too small for specified order (%d)", 
	  __FUNCTION__, w->omax, order) ;

  mon = w->iblock ; 
  mop_polynomial_order(p) = order ;
  make_index_list(mop_polynomial_order(p),
		  mop_polynomial_dimension(p),
		  mon, &ni) ;

  offset = mop_polynomial_point_number(p)*mop_polynomial_point_number(p) ;

  powers = &(w->block[offset+2*mop_polynomial_point_number(p)]) ;

  block_powers(mop_polynomial_points(p),
	       mop_polynomial_point_number(p),
	       mop_polynomial_dimension(p),
	       mop_polynomial_order(p),
	       powers) ;
  
  mop_polynomial_term_number(p) = 1 ;
  for ( i = 0 ; i < mop_polynomial_dimension(p) ; i ++ ) 
    mop_polynomial_monomial_power(p,0,i) = 0 ;

  for ( k = 1 ; k < ni ; k ++ ) {
    A = gsl_matrix_view_array(w->block, 
			      mop_polynomial_term_number(p)+1, 
			      mop_polynomial_point_number(p)) ;
    norm = gsl_vector_view_array(&(w->block[offset]),mop_polynomial_point_number(p)) ;
    tau = gsl_vector_view_array(&(w->block[offset+mop_polynomial_point_number(p)]),
				GSL_MIN(A.matrix.size1,A.matrix.size2)) ;
    for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
      for ( j = 0 ; j < mop_polynomial_point_number(p) ; j ++ ) {
	gsl_matrix_set(&(A.matrix), i, j,
		       block_multipower(powers, j, 
					mop_polynomial_point_number(p),
					mop_polynomial_dimension(p), 
					mop_polynomial_monomial(p,i))) ;
      }
    }
    for ( j = 0 ; j < mop_polynomial_point_number(p) ; j ++ ) {
      gsl_matrix_set(&(A.matrix), mop_polynomial_term_number(p), j, 
		     block_multipower(powers, j,
				      mop_polynomial_point_number(p),
				      mop_polynomial_dimension(p),
				      &(mon[k*mop_polynomial_dimension(p)]))) ;
    }
    rank = matrix_rank(&(A.matrix), tol,
		       &(tau.vector), &(norm.vector)) ;
    monomial_string(mns, &(mon[k*mop_polynomial_dimension(p)]),
		    mop_polynomial_dimension(p)) ;  
		    
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	  "%s: trial monomial %d (%s): rank=%d", __FUNCTION__, k, mns, rank) ;
    if ( rank == A.matrix.size1 ) {
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	    "%s: adding trial monomial %d", __FUNCTION__, k) ;
      for ( j = 0 ; j < mop_polynomial_dimension(p) ; j ++ ) {
	mop_polynomial_monomial_power(p,mop_polynomial_term_number(p),j)
	  = mon[k*mop_polynomial_dimension(p)+j] ;
      }
      mop_polynomial_term_number(p) ++ ;
    } else {
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	    "%s: A rank deficient, monomial %d rejected", 
	    __FUNCTION__, k) ;
    }

  }

  for ( i = 0 ; i < mop_polynomial_point_number(p) ; i ++ ) 
    mop_polynomial_index(p,i) = i ;
  
  return MOP_SUCCESS ;
}

/** 
 * Generate the basis for a ::mop_polynomial_t by selecting sufficient
 * points to match the number of monomial powers supplied. 
 * 
 * @param p a ::mop_polynomial, which should have been initialized by
 * a call to ::mop_polynomial_set_points;
 * @param order maximum order of monomial to employed;
 * @param tol tolerance to be used in determining rank of basis matrix; 
 * @param w suitably sized ::mop_polynomial_workspace.
 * 
 * @return MOP_SUCCESS on success
 */

gint mop_polynomial_basis_points(mop_polynomial_t *p, gint order,
				 gdouble tol, mop_polynomial_workspace_t *w)

{
  gint i, j, k, rank, ni, offset ;
  gdouble *powers ;
  gsl_matrix_view A ;
  gsl_vector_view tau, norm ;

  if ( mop_polynomial_order_max(p) < order )
    g_error("%s: order too high (%d) for polynomial (maximum order %d)", 
	    __FUNCTION__, order, mop_polynomial_order_max(p)) ;

  if ( w->npmax < mop_polynomial_point_number(p) )
    g_error("%s: workspace is not big enough (%d points) "
	    "for polynomial (%d points)", 
	    __FUNCTION__, w->npmax, mop_polynomial_point_number(p)) ;

  if ( w->dim < mop_polynomial_dimension(p) )
    g_error("%s: workspace dimension (%d) smaller than polynomial's (%d)", 
	    __FUNCTION__, w->dim, mop_polynomial_dimension(p)) ;

  if ( w->omax < order )
    g_error("%s: workspace maximum order (%d) "
	    "too small for specified order (%d)", 
	    __FUNCTION__, w->omax, order) ;

  for ( i = 0 ; i < mop_polynomial_point_number(p) ; i ++ )
     mop_polynomial_index(p,i) = -1 ;

  mop_polynomial_order(p) = order ;

  make_index_list(mop_polynomial_order(p),
		  mop_polynomial_dimension(p),
		  mop_polynomial_powers(p),
		  &(mop_polynomial_term_number(p))) ;

  offset = mop_polynomial_point_number(p)*mop_polynomial_point_number(p) ;

  powers = &(w->block[offset+2*mop_polynomial_point_number(p)]) ;

  block_powers(mop_polynomial_points(p),
	       mop_polynomial_point_number(p),
	       mop_polynomial_dimension(p),
	       mop_polynomial_order(p),
	       powers) ;
  
  mop_polynomial_index(p,0) = 0 ; ni = 1 ;

  for ( k = 1 ; 
	(k < mop_polynomial_point_number(p)) &&  (ni < mop_polynomial_term_number(p)) ; 
	k ++ ) {
/*     A = gsl_matrix_view_array(w->block,  */
/* 			      ni+1, */
/* 			      mop_polynomial_term_number(p)) ; */
/*     norm = gsl_vector_view_array(&(w->block[offset]), */
/* 				 mop_polynomial_term_number(p)) ; */
    A = gsl_matrix_view_array(w->block,
			      mop_polynomial_term_number(p), ni+1) ;
    norm = gsl_vector_view_array(&(w->block[offset]), 
				 A.matrix.size2) ;
    tau = gsl_vector_view_array(&(w->block[offset+mop_polynomial_point_number(p)]),
				GSL_MIN(A.matrix.size1,A.matrix.size2)) ;
    for ( i = 0 ; i < ni ; i ++ ) {
      for ( j = 0 ; j < mop_polynomial_term_number(p) ; j ++ ) {
	gsl_matrix_set(&(A.matrix), j, i,
		       block_multipower(powers, mop_polynomial_index(p,i), 
					mop_polynomial_point_number(p),
					mop_polynomial_dimension(p), 
					mop_polynomial_monomial(p,j))) ;
      }
    }
    for ( j = 0 ; j < mop_polynomial_term_number(p) ; j ++ ) {
      gsl_matrix_set(&(A.matrix), j, ni,
		     block_multipower(powers, k,
				      mop_polynomial_point_number(p),
				      mop_polynomial_dimension(p), 
				      mop_polynomial_monomial(p,j))) ;
    }
    rank = matrix_rank(&(A.matrix), tol,
		       &(tau.vector), &(norm.vector)) ;
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	  "%s: trial point %d: rank=%d", __FUNCTION__, k, rank) ;
    if ( rank == A.matrix.size2 ) {
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	    "%s: adding trial point %d", __FUNCTION__, k) ;
      mop_polynomial_index(p,ni) = k ;
      ni ++ ;
    } else {
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	    "%s: A rank deficient, trial point %d rejected", 
	    __FUNCTION__, k) ;
    }
  }

  if ( ni < mop_polynomial_term_number(p) ) return MOP_FAILURE ;
  
  return MOP_SUCCESS ;
}

/** 
 * Generate the discrete orthogonal polynomials associated with \a p,
 * using the method of Xu, Yuan, 2004, `On discrete orthogonal
 * polynomials of several variables', Advances in Applied Mathematics,
 * 33:615--632, doi:10.1016/j.aam.2004.03.002.
 * 
 * @param p a ::mop_polynomial, which should have been initialized by
 * a call to ::mop_polynomial_basis_power or 
 * ::mop_polynomial_basis_points; 
 * @param w a suitably sized ::mop_polynomial_workspace, generated by 
 * ::mop_polynomial_workspace_alloc.
 * 
 * @return MOP_SUCCESS on success.
 */

gint mop_polynomial_make(mop_polynomial_t *p, 
			 mop_polynomial_workspace_t *w)

{
  gsl_matrix_view vM, vS ;
  gsl_vector_view vD, viD, vr ;
  gsl_vector *iD, *D ;
  gsl_matrix *M, *S ;
  gint i, j, k, m, n, *pi, *pj ;
  gsl_permutation *P ;
  gdouble tmp ;

  if ( mop_polynomial_term_number(p) > mop_polynomial_point_number(p) )
    g_error("%s: number of monomial powers (%d) greater "
	    "than number of points (%d)",
	    __FUNCTION__, mop_polynomial_term_number(p),
	    mop_polynomial_point_number(p)) ;

  n = mop_polynomial_term_number(p) ;

  vM = gsl_matrix_view_array(w->block, n, n) ;
  vS = gsl_matrix_view_array(&(w->block[n*n]), n, n) ;
  vD = gsl_vector_view_array(&(w->block[2*n*n]), n) ;
  viD = gsl_vector_view_array(&(w->block[2*n*n+n]), n) ;
  P = gsl_permutation_alloc(n) ;
  
  M = &(vM.matrix) ; S = &(vS.matrix) ;
  D = &(vD.vector) ; iD = &(viD.vector) ; 

  gsl_matrix_set_zero(M) ;

  block_powers(mop_polynomial_points(p),
	       mop_polynomial_point_number(p),
	       mop_polynomial_dimension(p),
	       mop_polynomial_order(p),
	       p->xp) ;

  for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
    pi = mop_polynomial_monomial(p, i) ;
    for ( j = 0 ; j < mop_polynomial_term_number(p) ; j ++ ) {
      pj = mop_polynomial_monomial(p, j) ;
      for ( (tmp = 0.0), (m = 0) ; m < mop_polynomial_term_number(p) ; m ++ ) {
	k = mop_polynomial_index(p,m) ;
	tmp +=
	  block_multipower(p->xp, k, 
			   mop_polynomial_point_number(p), 
			   mop_polynomial_dimension(p),
			   pi)*
	  block_multipower(p->xp, k, 
			   mop_polynomial_point_number(p), 
			   mop_polynomial_dimension(p),
			   pj)*mop_polynomial_weight(p,k) ;

	if ( isnan(tmp) ) 
	  g_error("%s: NaN error at i=%d, j=%d, k=%d", 
		  __FUNCTION__, i, j, k) ;	     
      }
      gsl_matrix_set(M, i, j, tmp) ;
    }
  }

  if ( sds_decomp(M, S, D, iD) != 0 ) 
    return MOP_NO_BASIS ;

  gsl_permutation_init(P) ;
  for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
    vr = gsl_vector_view_array(&(p->R[i*mop_polynomial_term_number(p)]), 
			       mop_polynomial_term_number(p)) ;
    gsl_vector_set_zero(iD) ;
    j = mop_polynomial_index(p,i) ;

    if ( isnan(tmp=1.0/sqrt(mop_polynomial_weight(p,j)*gsl_vector_get(D,i))) )
      g_error("%s: NaN error in iD, i=%d, w=%lg, D_i=%lg", 
	      __FUNCTION__, i, mop_polynomial_weight(p,j),
	      gsl_vector_get(D,i)) ;
    
    gsl_vector_set(iD, i, tmp) ;
		   
    gsl_linalg_LU_solve(S, P, iD, &(vr.vector)) ;
    for ( j = 0 ; j < vr.vector.size ; j ++ ) {
      if ( isnan(gsl_vector_get(&(vr.vector),j)) ) {
	g_error("%s: NaN error in R, i=%d, j=%d",
		__FUNCTION__, i, j) ;
      }
    }
  }

  gsl_permutation_free(P) ;

  return MOP_SUCCESS ;
}

/** 
 * Write a ::mop_polynomial_t as a fragment of LaTeX code.
 * 
 * @param p ::mop_polynomial_t to write;
 * @param f file pointer to write to.
 * 
 * @return MOP_SUCCESS on success.
 */

gint mop_polynomial_write_latex(mop_polynomial_t *p, FILE *f)

{
  gint i, j, k ;
  gboolean plus ;

  for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
    fprintf(f, "P_{%d} &= ", i) ;
    plus = FALSE ;
    for ( j = 0 ; j < mop_polynomial_term_number(p) ; j ++ ) {
      if ( fabs(mop_polynomial_coefficient(p, i, j)) > 1e-9 ) {
	if ( plus && mop_polynomial_coefficient(p, i, j) > 0) 
	  fprintf(f, "+") ;
	fprintf(f, "%.2f", mop_polynomial_coefficient(p, i, j)) ;
	for ( k = 0 ; k < mop_polynomial_dimension(p) ; k ++ ) {
	  if ( mop_polynomial_monomial_power(p,j,k) == 1 )
	    fprintf(f, "x_{%d}", k) ;
	  if ( mop_polynomial_monomial_power(p,j,k) > 1 )
	    fprintf(f, "x_{%d}^{%d}",
		    k, mop_polynomial_monomial_power(p,j,k)) ;
	}
	plus = TRUE ;
      }
    }
    if ( i < mop_polynomial_term_number(p)-1 ) fprintf(f, "\\\\\n") ;
  }

  return MOP_SUCCESS ;
}

/** 
 * Scale the coefficients of a ::mop_polynomial_t to give unit inner
 * product, \f$\sum_{i}P_{j}^{2}(x_{i})w_{i}\equiv 1\f$.
 * 
 * @param p ::mop_polynomial_t to normalize;
 * @param w a ::mop_polynomial_workspace_t of appropriate size.
 * 
 * @return MOP_SUCCESS on success.
 */

gint mop_polynomial_normalize(mop_polynomial_t *p,
			      mop_polynomial_workspace_t *w)

{
  gdouble *ws ;
  gint i, j ;
  gdouble s ;

  ws = w->block ;
  mop_polynomial_eval_base(p, ws) ;

/*   for ( i = 0 ; i < mop_polynomial_point_number(p) ; i ++ ) { */
  for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
    for ( (s = 0.0), (j = 0) ; j < mop_polynomial_term_number(p) ; j ++ )
      s += ws[j*mop_polynomial_term_number(p)+i]*
	ws[j*mop_polynomial_term_number(p)+i]*
	mop_polynomial_weight(p,mop_polynomial_index(p,j)) ;

    s = sqrt(s) ;
    if ( s == 0.0 ) {
      g_error("%s: zero inner product for polynomial %d",
	      __FUNCTION__, i) ;
    }
    for ( j = 0 ; j < mop_polynomial_term_number(p) ; j ++ )
      mop_polynomial_coefficient(p,i,j) /= s ;
  }

  return MOP_SUCCESS ;
}


/** 
 * Evaluate a ::mop_polynomial_t at a point \a x. 
 * 
 * @param p ::mop_polynomial_t to evaluate;
 * @param x coordinates of evaluation point;
 * @param P array containing values of \a p at \a x.
 * 
 * @return MOP_SUCCESS on success.
 */

gint mop_polynomial_eval(mop_polynomial_t *p, gdouble *x,
			 gdouble *P)

{
  gint i, j ;

  block_powers(x, 1, 
	       mop_polynomial_dimension(p),
	       mop_polynomial_order(p),
	       p->xp) ;

  for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
    for ( (P[i] = 0.0), (j = 0) ; j < mop_polynomial_term_number(p) ; j ++ )
      P[i] += mop_polynomial_coefficient(p,i,j)*
	block_multipower(p->xp, 0, 1, 
			 mop_polynomial_dimension(p),
			 mop_polynomial_monomial(p,j)) ;
  }

  return MOP_SUCCESS ;
}

/** 
 * Evaluate a ::mop_polynomial_t at its base points. Note that the base
 * points used are those in the index list of \a p and they are used
 * in the order in that list. To connect a given value of orthogonal
 * polynomial to a particular base point, use ::mop_polynomial_index. 
 * 
 * @param p ::mop_polynomial_t to evaluate;
 * @param P array of \a p evaluated at its base points so that 
 * \f$P_{i}(x_{j})\f$ is P[j*mop_polynomial_term_number(p)+i].
 * 
 * @return MOP_SUCCESS on success.
 */

gint mop_polynomial_eval_base(mop_polynomial_t *p, gdouble *P)

{
  gint i, j, k, n ;

  block_powers(mop_polynomial_points(p),
	       mop_polynomial_point_number(p),
	       mop_polynomial_dimension(p),
	       mop_polynomial_order(p),
	       p->xp) ;
  
  for ( k = 0 ; k < mop_polynomial_term_number(p) ; k ++ )
    for ( i = 0 ; i < mop_polynomial_term_number(p) ; i ++ ) {
      for ( (P[k*mop_polynomial_term_number(p)+i] = 0.0), (n = 0) ; 
	    n < mop_polynomial_term_number(p) ; n ++ ) {
	j = mop_polynomial_index(p,k) ;
	if ( isnan(P[k*mop_polynomial_term_number(p)+i] +=
		   mop_polynomial_coefficient(p,i,n)*
		   block_multipower(p->xp, j,
				    mop_polynomial_point_number(p),
				    mop_polynomial_dimension(p),
				    mop_polynomial_monomial(p,n))) )
	  g_error("%s: NaN error at i=%d, j=%d, k=%d", 
		  __FUNCTION__, i, j, k) ;
      }
    }  

  return MOP_SUCCESS ;
}

gdouble mop_polynomial_point_magnitude2(mop_polynomial_t *p, gint i)

{
  gdouble R2 = 0.0 ;
  gint j ;

  for ( j = 0 ; j < mop_polynomial_dimension(p) ; j ++ ) 
    R2 += p->x[i*mop_polynomial_dimension(p)+j]*
      p->x[i*mop_polynomial_dimension(p)+j] ;

  return R2 ;
}

gint mop_polynomial_points_scale(mop_polynomial_t *p, gint method)

{
  gint i ;
  gdouble R2 ;

/*   for ( (i = 1), (R2 = G_MAXDOUBLE) ; i < mop_polynomial_point_number(p) ; i ++ ) */
/*     R2 = MIN(R2,mop_polynomial_point_magnitude2(p, i)) ; */
  for ( (i = 1), (R2 = 0.0) ; i < mop_polynomial_point_number(p) ; i ++ )
    R2 = MAX(R2,mop_polynomial_point_magnitude2(p, i)) ;
/*     R2 += mop_polynomial_point_magnitude2(p, i) ; */
/*   R2 /= mop_polynomial_point_number(p) ; */

  p->sc = 1.0/sqrt(R2) ;

  g_debug("%s: scale set to %lg", __FUNCTION__, p->sc) ;

  for ( (i = 0) ; 
	i < mop_polynomial_point_number(p)*mop_polynomial_dimension(p) ; i ++ )
    p->x[i] *= p->sc ;

  return MOP_SUCCESS ;
}

/**
 * @}
 * 
 */
