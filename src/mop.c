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

#include <blaswrap.h>

#include "mop.h"
#include "mopblock.h"

#ifdef _HAVE_CONFIG_H_
#include "config.h"
#endif

#define matrix_index_lower(_n,_i,_j) ((_j)*(2*(_n)-(_j)-1)/2+(_i))

static gint ldl_decomp_tp(gdouble *A, gint n, gdouble *v)

/*
 * LDL^T decomposition for symmetric matrix A, using FORTRAN packed
 * ordering. On exit A is overwritten with D in its diagonal entries
 * and the off-diagonal terms filled with the corresponding entries of
 * L.
 *
 * (Golub and van Loan page 138)
 */
  
{
  gint i, j, k ;

  for ( j = 0 ; j < n ; j ++ ) {
    for ( i = 0 ; i < j ; i ++ ) {
      v[i] = A[matrix_index_lower(n,j,i)]*A[matrix_index_lower(n,i,i)] ;
    }

    v[j] = A[matrix_index_lower(n,j,j)] ;
    for ( i = 0 ; i < j ; i ++ ) {
      v[j] -= A[matrix_index_lower(n,j,i)]*v[i] ;
    }

    A[matrix_index_lower(n,j,j)]  = v[j] ;
    for ( i = j+1 ; i < n ; i ++ ) {
      for ( k = 0 ; k < j ; k ++ ) {
	A[matrix_index_lower(n,i,j)] -= A[matrix_index_lower(n,i,k)]*v[k] ;
      }
      A[matrix_index_lower(n,i,j)] /= v[j] ;
    }
  }
  
  return 0 ;
}

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

static gint matrix_rank_svd(gdouble *A, gint m, gint n, gint lda,
			    gdouble tol, gdouble *work, gint lwork)

{
  gdouble *s ;
  gint one = 1, info, rank ;

  /* lwork = MAX(3*MIN(m,n)+MAX(m,n),5*MIN(m,n)-4) ; */
  /* info = 0 ; */
  /* dgesvd_("N", "N", &n, &m, A, &lda, s, NULL, &one, NULL, &one, work, */
  /* 	  &lwork, &info) ; */
  /* lwork = work[0] ; */

  lwork -= MAX(m,n) ;
  s = &(work[lwork]) ;
  info = 1 ;
  
  dgesvd_("N", "N", &m, &n, A, &lda, s, NULL, &one, NULL, &one, work,
	  &lwork, &info) ;

  for ( rank = 1 ; rank < n ; rank ++ ) {
    if ( fabs(s[rank]/s[0]) < tol ) return rank ;
  }
  
  return rank ;
}

static gint monomial_matrix(gdouble *X, gint nr, gint nc,
			    gdouble *powers, gint *p, gint npts,
			    gint d, gint *monomials, gboolean trans)

/*
 * generate matrix X with X_{ij} = x_{i}^p_{j} in FORTRAN ordering (if
 * !trans)
 *
 * powers: powers of individual components of points
 * p:      indices of points to use
 * npts:   number of points in block (not number of rows)
 * d:      dimension of problems
 * monomials: monomial powers
 */
  
{
  gint i, j ;

  if ( !trans ) {
    for ( i = 0 ; i < nr ; i ++ ) {
      for ( j = 0 ; j < nc ; j ++ ) {
	X[j*nr+i] =
	  block_multipower(powers, p[i], npts, d, &(monomials[d*j])) ;
      }
    }    
    
    return 0 ;
  }

  for ( i = 0 ; i < nr ; i ++ ) {
    for ( j = 0 ; j < nc ; j ++ ) {
      X[j*nr+i] =
	block_multipower(powers, p[j], npts, d, &(monomials[d*i])) ;
    }
  }    

    return 0 ;
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
  w->block = (gdouble *)g_malloc0(4*bsize*sizeof(gdouble)) ;
  
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
  gint i, j, k, rank, ni, *mon, offset, np, nt ;
  gdouble *powers, *A ;
  gchar mns[32] ;
  static gint calls = 0 ;

  calls ++ ;
  
  g_assert(mop_polynomial_point_number(p) <=
	   mop_polynomial_term_number_max(p)) ;
  
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

  np = mop_polynomial_point_number(p) ; rank = 0 ;
  for ( k = 1 ; (k < ni) && (mop_polynomial_term_number(p) < np) ; k ++ ) {
    nt = mop_polynomial_term_number(p) ;
    A = w->block ; /*A is size (nt+1)*np and FORTRAN indexed*/

    /*put the trial monomial in the last slot in the list*/
    for ( j = 0 ; j < mop_polynomial_dimension(p) ; j ++ ) {
      mop_polynomial_monomial_power(p,mop_polynomial_term_number(p),j)
	= mon[k*mop_polynomial_dimension(p)+j] ;
    }
    
    monomial_matrix(A, nt+1, np, powers, p->perm, np, 
		    mop_polynomial_dimension(p),
		    mop_polynomial_monomial(p,0), TRUE) ;
    
    gdouble work[1024] ;
    rank = matrix_rank_svd(A, nt+1, np, nt+1, tol, work, 1024) ;
    monomial_string(mns, &(mon[k*mop_polynomial_dimension(p)]),
		    mop_polynomial_dimension(p)) ;  
		    
    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	  "%s: trial monomial %d (%s): rank=%d", __FUNCTION__, k, mns, rank) ;
    if ( rank == nt+1 ) {
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	    "%s: adding trial monomial %d", __FUNCTION__, k) ;
      /*trial monomial is already in the last slot, we just need to
	increment the count to include it*/
      mop_polynomial_term_number(p) ++ ;
    } else {
      g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	    "%s: A rank deficient, monomial %d rejected", 
	    __FUNCTION__, k) ;
    }

  }

  if ( rank < mop_polynomial_point_number(p) ) {
    mop_polynomial_point_number(p) -- ;

    return mop_polynomial_basis_power(p, order, tol, w) ;
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
  gint i, j, k, rank, ni, offset, npts, nterm, nr, nc ;
  gdouble *powers, *A ;

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

  npts = mop_polynomial_point_number(p) ;
  
  for ( i = 0 ; i < npts ; i ++ ) mop_polynomial_index(p,i) = -1 ;

  mop_polynomial_order(p) = order ;

  make_index_list(mop_polynomial_order(p),
		  mop_polynomial_dimension(p),
		  mop_polynomial_powers(p),
		  &(mop_polynomial_term_number(p))) ;

  nterm = mop_polynomial_term_number(p) ;
  offset = npts*npts ;

  powers = &(w->block[offset+2*npts]) ;

  block_powers(mop_polynomial_points(p),
	       npts,
	       mop_polynomial_dimension(p),
	       mop_polynomial_order(p),
	       powers) ;
  
  mop_polynomial_index(p,0) = 0 ; ni = 1 ;

  /*check rank with all the points and powers in place*/
  /* nr = ni+1 ; nc = ni+1 ; */
  
  ni = nr = npts ; nc = nterm = npts ;
  
  for ( i = 0 ; i < ni ; i ++ ) mop_polynomial_index(p,i) = i ;

  A = w->block ; /*A is size nterm*(ni+1) and FORTRAN indexed*/
  monomial_matrix(A, nr, nc, powers, p->perm, npts, 
		  mop_polynomial_dimension(p),
		  mop_polynomial_monomial(p,0), FALSE) ;
    
  gdouble work[1024] ;
  /* rank = matrix_rank_svd(A, ni+1, ni+1, ni+1, tol, work, 1024) ; */
  /* rank = matrix_rank_svd(A, ni+1, nterm, ni+1, tol, work, 1024) ; */
  rank = matrix_rank_svd(A, nr, nc, nr, tol, work, 1024) ;

  if ( rank == nr ) {
    fprintf(stderr, "%s: full rank (%d)\n", __FUNCTION__, rank) ;

    return 0 ;
  } else {
    fprintf(stderr, "%s: rank deficient (%d)\n", __FUNCTION__, rank) ;    
  }
  
  /* g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG, */
  /* 	"%s: trial point %d: rank=%d; ni=%d", __FUNCTION__, k, rank, ni) ; */
  /* if ( rank == ni+1 ) { */
  /*   g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG, */
  /* 	  "%s: adding trial point %d", __FUNCTION__, k) ; */
  /*   mop_polynomial_index(p,ni) = k ; */
  /*   ni ++ ; */
  /* } else { */
  /*   g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG, */
  /* 	  "%s: A rank deficient, trial point %d rejected",  */
  /* 	  __FUNCTION__, k) ; */
  /* } */
  
#if 0  
  for ( k = 1 ; (k < npts) &&  (ni < nterm) ; k ++ ) {
    /* nr = nterm ; nc = ni+1 ; */
    nr = ni+1 ; nc = ni+1 ;
    A = w->block ; /*A is size nterm*(ni+1) and FORTRAN indexed*/
    for ( i = 0 ; i < ni ; i ++ ) {
      for ( j = 0 ; j < nc ; j ++ ) {
      /* for ( j = 0 ; j < ni+1 ; j ++ ) { */
	/* A[i*nterm+j] = */
	/* A[i*(ni+1)+j] = */
	A[j*nr+i] =
	  block_multipower(powers, mop_polynomial_index(p,i),
			   npts,
			   mop_polynomial_dimension(p),
			   mop_polynomial_monomial(p,j)) ;
      }
    }
    for ( j = 0 ; j < nc ; j ++ ) {
    /* for ( j = 0 ; j < ni+1 ; j ++ ) { */
      /* A[ni*nterm+j] = */
      /* A[ni*(ni+1)+j] = */
      A[j*nr+ni] =
	block_multipower(powers, k,
			 npts,
			 mop_polynomial_dimension(p),
			 mop_polynomial_monomial(p,j)) ;
    }

    if ( 0 ) {
      gint ii, jj ;
      fprintf(stdout, "A%d = [\n", k) ;
      for ( ii = 0 ; ii < nr ; ii ++ ) {
    	for ( jj = 0 ; jj < nc ; jj ++ ) {
    	  fprintf(stdout, "%1.16e ", A[jj*nr+ii]) ;
    	}
	fprintf(stdout, "; \n") ;
      }
      fprintf(stdout, "] ; \n") ;
    }
    
    gdouble work[1024] ;
    /* rank = matrix_rank_svd(A, ni+1, ni+1, ni+1, tol, work, 1024) ; */
    /* rank = matrix_rank_svd(A, ni+1, nterm, ni+1, tol, work, 1024) ; */
    rank = matrix_rank_svd(A, nr, nc, nr, tol, work, 1024) ;

    g_log(G_LOG_DOMAIN, G_LOG_LEVEL_DEBUG,
	  "%s: trial point %d: rank=%d; ni=%d", __FUNCTION__, k, rank, ni) ;
    if ( rank == ni+1 ) {
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

  if ( ni < nterm ) {
    /*adjust number of terms in polynomials to match number of points
      available*/
    /* g_assert_not_reached() ; /\*this needs a rethink*\/ */
    mop_polynomial_term_number(p) = ni ;
    /* mop_polynomial_point_number(p) = ni ; */

    return MOP_FAILURE ;
  }
#endif
  
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
  gdouble *M, *v ;
  gint n, i, j, k, m, *pi, *pj, info ;
  static gint calls = 0 ;
  gboolean fail ;
  
  calls ++ ;
  fail = FALSE ;
  
  n = mop_polynomial_term_number(p) ;

  M = &(w->block[0]) ;
  v = &(M[n*(n+1)/2]) ;

  block_powers(mop_polynomial_points(p),
	       mop_polynomial_point_number(p),
	       mop_polynomial_dimension(p),
	       mop_polynomial_order(p),
	       p->xp) ;

  for ( j = 0 ; j < n ; j ++ ) {
    pj = mop_polynomial_monomial(p, j) ;
    for ( i = j ; i < n ; i ++ ) {
      pi = mop_polynomial_monomial(p, i) ;
      M[matrix_index_lower(n,i,j)] = 0.0 ;
      for (m = 0 ; m < mop_polynomial_term_number(p) ; m ++ ) {
	k = mop_polynomial_index(p,m) ;
	if ( k < 0 )
	  g_error("%s: polynomial index (%d) out of bounds", __FUNCTION__, k) ;
	/* fprintf(stderr, "k=%d\n", k) ; */
	M[matrix_index_lower(n,i,j)] +=
	  block_multipower(p->xp, k, 
			   mop_polynomial_point_number(p), 
			   mop_polynomial_dimension(p),
			   pi)*
	  block_multipower(p->xp, k, 
			   mop_polynomial_point_number(p), 
			   mop_polynomial_dimension(p),
			   pj)*mop_polynomial_weight(p,k) ;
      }
    }
  }
  
  memset(p->R, 0, n*n*sizeof(gdouble)) ;

  if ( ldl_decomp_tp(M, n, v) != 0 ) return MOP_NO_BASIS ;

  for ( i = 0 ; i < n ; i ++ ) {
    j = mop_polynomial_index(p,i) ;
    p->R[i*n+i] = 1.0/sqrt(mop_polynomial_weight(p,j)*
			   M[matrix_index_lower(n,i,i)]) ;
  }
  
  dtptrs_("L", "T", "U", &n, &n, M, p->R, &n, &info) ;

  /*NaN check*/
  for ( i = 0 ; i < n ; i ++ ) {
    for ( j = 0 ; j < n ; j ++ ) {
      if ( isnan(p->R[i*n+j]) ) {
	fail = TRUE ;
	g_error("%s: coefficient (%d,%d) NaN", __FUNCTION__, i, j) ;
      }
    }    
  }

  if ( fail ) {
    /* for ( i = 0 ; i < mop_polynomial_point_n */
    g_error("%s: coefficient (%d,%d) NaN", __FUNCTION__, i, j) ;
  }
  
  return 0 ;
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
