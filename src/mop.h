/* mop.h
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

#ifndef _MOP_H_
#define _MOP_H_

#include <glib.h>

enum {
  MOP_SUCCESS = 0,
  MOP_FAILURE = -1,   /*failure for unspecified reasons*/
  MOP_NO_BASIS = 1,   /*unable to find a basis for the specified points*/
/*   MOP_ */
/*   MOP_ */
/*   MOP_ */
} ;

typedef struct _mop_polynomial_t mop_polynomial_t ;

struct _mop_polynomial_t {
  gint d,      /*dimension of system (1, 2, or 3 for now)*/
    np, npmax, /*number of points in the system*/
    o, omax,   /*order (linear, quadratic, etc.)*/
    nt, ntmax, /*number of terms/polynomials*/
    *p,        /*powers in monomials*/
    *perm ;    /*permutation of points*/
  gdouble *x,  /*points*/
    *w,        /*weights*/
    *R,        /*coefficients*/
    sc,        /*length scale*/
    *xp ;      /*powers of x as workspace*/
} ;

typedef struct _mop_polynomial_workspace_t mop_polynomial_workspace_t ;

struct _mop_polynomial_workspace_t {
  gint npmax, imax, *iblock, omax, dim ;  
  gdouble *block ;
} ;

#define mop_polynomial_dimension(p) (p->d)
#define mop_polynomial_point_number(p) (p->np)
#define mop_polynomial_point_number_max(p) (p->npmax)
#define mop_polynomial_order(p) (p->o)
#define mop_polynomial_order_max(p) (p->omax)
#define mop_polynomial_term_number(p) (p->nt)
#define mop_polynomial_term_number_max(p) (p->ntmax)
#define mop_polynomial_powers(p) (p->p)
#define mop_polynomial_monomial(p,_i) (&(p->p[(_i)*(p->d)]))
#define mop_polynomial_monomial_power(p,_i,_j) (p->p[(_i)*(p->d)+(_j)])
#define mop_polynomial_indices(p) (p->perm)
#define mop_polynomial_index(p,_i) (p->perm[(_i)])
#define mop_polynomial_coefficient(p,_i,_j) (p->R[_i*(p->nt)+(_j)])

#define mop_polynomial_points(p) (p->x)
#define mop_polynomial_weights(p) (p->w)
#define mop_polynomial_weight(p,_i) (p->w[(_i)])

mop_polynomial_t *mop_polynomial_alloc(gint np, gint dim, gint order) ;
gint mop_polynomial_free(mop_polynomial_t *p) ;

gint mop_number_of_terms(gint dim, gint order) ;
gint mop_polynomial_set_points(mop_polynomial_t *p, 
			       gdouble *x, gdouble *w, gint n) ;
gint mop_polynomial_points_scale(mop_polynomial_t *p, gint method) ;

gint mop_polynomial_write(mop_polynomial_t *p, FILE *f) ;
gint mop_polynomial_write_latex(mop_polynomial_t *p, FILE *f) ;

gint mop_polynomial_basis_power(mop_polynomial_t *p, gint order,
				gdouble tol, mop_polynomial_workspace_t *w) ;
gint mop_polynomial_basis_points(mop_polynomial_t *p, gint order,
				 gdouble tol, mop_polynomial_workspace_t *w) ;

gint mop_polynomial_make(mop_polynomial_t *p, 
			 mop_polynomial_workspace_t *w) ;
gint mop_polynomial_normalize(mop_polynomial_t *p,
			      mop_polynomial_workspace_t *w) ;
gint mop_polynomial_eval(mop_polynomial_t *p, gdouble *x,
			 gdouble *P) ;
gint mop_polynomial_eval_base(mop_polynomial_t *p, gdouble *P) ;

gint mop_polynomial_transform(mop_polynomial_t *p, gdouble *f, gint n,
			      gdouble *c, mop_polynomial_workspace_t *w) ;
gint mop_interpolate(mop_polynomial_t *p, gdouble *c, gint n,
		     gdouble *x, gdouble *f) ;

gint mop_interpolation_weights(mop_polynomial_t *p, 
			       gdouble *x, gdouble *v,
			       mop_polynomial_workspace_t *w) ;

gint mop_polynomial_differentiate(mop_polynomial_t *p, gdouble *x,
				  gint *d, gdouble *P) ;
gint mop_differentiation_weights(mop_polynomial_t *p, 
				 gdouble *x, gint *d, gdouble *v,
				 mop_polynomial_workspace_t *w) ;
gdouble mop_polynomial_point_magnitude2(mop_polynomial_t *p, gint i) ;


mop_polynomial_workspace_t *mop_polynomial_workspace_alloc(gint np,
							   gint dim,
							   gint order) ;
gint mop_polynomial_workspace_free(mop_polynomial_workspace_t *w) ;

gint mop_logging_init(FILE *f, gchar *p, 
		      GLogLevelFlags log_level,
		      gpointer exit_func) ;
#endif /*_MOP_H_*/
