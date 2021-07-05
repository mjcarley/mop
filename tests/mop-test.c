/* mop-test.c
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <glib.h>

#include <mop.h>

/* static gdouble distance2(gdouble *x, gdouble *y, gint dim) */

/* { */
/*   gdouble R ; */
/*   gint i ; */

/*   for ( (i = 0), (R = 0.0) ; i < dim ; i ++ ) */
/*     R += (x[i]-y[i])*(x[i]-y[i]) ; */

/*   return R ; */
/* } */

/* static gint compare_distance(gconstpointer a, gconstpointer b, gpointer data) */

/* { */
/*   gpointer *f = (gpointer *)data ; */
/*   gint i, j, dim ; */
/*   gdouble *x0, *x, *x1, *x2, R1, R2 ; */

/*   x0 = (gdouble *)f[0] ; */
/*   x = (gdouble *)f[1] ; */
/*   dim = *(gint *)f[2] ; */
/*   i = *(gint *)a ; j = *(gint *)b ; */

/*   x1 = &(x[i*dim]) ; x2 = &(x[j*dim]) ; */
/*   R1 = distance2(x0, x1, dim) ; */
/*   R2 = distance2(x0, x2, dim) ; */
  
/*   if ( R1 < R2 ) return -1 ; */
/*   if ( R1 > R2 ) return 1 ; */

/*   return 0 ; */
/* } */

/* static gint sort_point_list(gint *idx, gdouble *x, gint np, gint nc,  */
/* 		     gdouble *y) */

/* { */
/*   gpointer data[5] ; */

/*   data[0] = y ; data[1] = x ; data[2] = &nc ; */

/*   g_qsort_with_data(idx, np, sizeof(gint), compare_distance, data) ; */

/*   return 0 ; */
/* } */

static gint parse_gradient(gchar *arg, gint *d)

{
  gchar **tokens ;
  gint i ;

  tokens = g_strsplit(arg, ",", 0) ;

  for ( i = 0 ; tokens[i] != NULL ; i ++ ) d[i] = atoi(tokens[i]) ;

  return 0 ;
}

gint main(gint argc, gchar **argv)

{
  gdouble *x, *xs, *wt, *f, *c ;
  gdouble y[32], tf[32], g[32], tol, ranktol ;
  gint pmax ;
  gint np, nc, nw, nf, ny, i, j, k ;
  gint d[32], *idx ;
  gboolean verbose, sort_points, basis_power, rescale_points, normalize_poly ;
  mop_polynomial_t *p ;
  mop_polynomial_workspace_t *w ;
  gchar ch ;
  FILE *input ;
  
  mop_logging_init(stderr, "", G_LOG_LEVEL_MESSAGE, NULL) ;
  /* mop_logging_init(stderr, "", G_LOG_LEVEL_DEBUG, NULL) ; */

  verbose = TRUE ;
  tol = 1e-2 ;
  pmax = 4 ;
  input = stdin ;
  d[0] = d[1] = d[2] = d[3] = d[4] = d[5] = d[6] = d[7] = d[8] = 0 ;
  ranktol = 1e-6 ;
  sort_points = FALSE ;
  basis_power = TRUE ;
  rescale_points = FALSE ;
  normalize_poly = TRUE ;
  
  while ( (ch = getopt(argc, argv, "hd:nt:p:Pr:sS:")) != EOF ) {
    switch (ch) {
    default:
    case'h':
      fprintf(stderr, 
	      "Usage: %s <options> < <input file>\n\n"
	      "Options: \n\n"
	      "  -d #,#,... derivatives in each dimension\n"
	      "  -n do not orthonormalize polynomials\n"
	      "  -p # maximum power/order of polynomials (%d)\n"
	      "  -r # tolerance for rank deficiency test (%lg)\n"
	      "  -s sort points before computing polynomials\n"
	      "  -S rescale points before computing polynomials\n\n"
	      "Input file format:\n\n"
	      "<number of points> <dimension> <number of functions> "
	      "<number of weights>\n"
	      "x_{i} y_{i} ... w_{i} f_{i}\n"
	      "<number of evaluation points>\n"
	      "x_{i} y_{i} ... f_{i}\n",
	      argv[0], pmax, ranktol
	      ) ;
      return 0 ;
      break ;
    case 'd':
      parse_gradient(optarg, d) ;
      break ;
    case 'n': normalize_poly = FALSE ; break ;
    case 't': tol = atof(optarg) ; break ;
    case 'p': pmax = atoi(optarg) ; break ;
    case 'P': basis_power = FALSE ; break ;
    case 'r': ranktol = atof(optarg) ; break ;
    case 's': sort_points = TRUE ; break ;
    case 'S': rescale_points = TRUE ; break ;
    }
  }

  fscanf(input, "%d", &np) ;
  fscanf(input, "%d", &nc) ;
  fscanf(input, "%d", &nf) ;
  fscanf(input, "%d", &nw) ;

  x = (gdouble *)g_malloc(np*nc*sizeof(gdouble)) ;
  xs = (gdouble *)g_malloc(np*nc*sizeof(gdouble)) ;
  idx = (gint *)g_malloc(np*sizeof(gdouble)) ;
  f = (gdouble *)g_malloc(np*nf*sizeof(gdouble)) ;
  c = (gdouble *)g_malloc(np*np*sizeof(gdouble)) ;
  if ( nw == 0 ) wt = NULL ;
  else wt = (gdouble *)g_malloc(np*sizeof(gdouble)) ;

  if ( wt == NULL ) {
    for ( i = 0 ; i < np ; i ++ ) {
      for ( j = 0 ; j < nc ; j ++ ) fscanf(input, "%lg", &x[i*nc+j]) ;
      for ( j = 0 ; j < nf ; j ++ ) fscanf(input, "%lg", &f[i*nf+j]) ;	
    } 
  } else {
    for ( i = 0 ; i < np ; i ++ ) {
      for ( j = 0 ; j < nc ; j ++ ) fscanf(input, "%lg", &x[i*nc+j]) ;
      fscanf(input, "%lg", &wt[i]) ;
      for ( j = 0 ; j < nf ; j ++ ) fscanf(input, "%lg", &f[i*nf+j]) ;	
    }
  }

  p = mop_polynomial_alloc(np, nc, pmax) ;
  w = mop_polynomial_workspace_alloc(np, nc, pmax) ;

  mop_polynomial_set_points(p, x, wt, np) ;
  if ( rescale_points ) mop_polynomial_points_scale(p, 0) ;

  if ( basis_power ) 
    mop_polynomial_basis_power(p, pmax, ranktol, w) ;
  else
    mop_polynomial_basis_points(p, pmax, ranktol, w) ;

  mop_polynomial_make(p, w) ;

  if ( normalize_poly ) mop_polynomial_normalize(p, w) ;

  if ( nf == 0 ) {
    mop_polynomial_write(p, stdout) ;
    
    return 0 ;
  }

  mop_polynomial_transform(p, f, nf, c, w) ;

  fscanf(input, "%d", &ny) ;
  for ( i = 0 ; i < ny ; i ++ ) {
    for ( j = 0 ; j < nc ; j ++ ) fscanf(input, "%lg", &y[j]) ;
    for ( j = 0 ; j < nf ; j ++ ) fscanf(input, "%lg", &tf[j]) ;

    mop_differentiation_weights(p, y, d, c, w) ;
    
    for ( j = 0 ; j < nf ; j ++ ) g[j] = 0.0 ;
    for ( j = 0 ; j < mop_polynomial_nterms(p) ; j ++ ) {
      k = mop_polynomial_index(p,j) ;
      g[0] += c[j]*f[k*nf+0] ;
    }

    for ( j = 0 ; j < nf ; j ++ )
      fprintf(stdout, " %1.16e %1.16e %1.16e", g[j], tf[j], g[j]-tf[j]) ;
    fprintf(stdout, "\n") ;
  }

  mop_polynomial_free(p) ;
  mop_polynomial_workspace_free(w) ;
  
  return 0 ;
}
